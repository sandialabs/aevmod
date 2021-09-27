/* =====================================================================================
aevmod version 1.0
Copyright (2021) NTESS
https://github.com/sandialabs/aevmod

Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS). 
Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
certain rights in this software.

This file is part of aevmod. aevmod is open-source software: you can redistribute it
and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the license is also
provided under the main directory

Questions? Contact Habib Najm at <hnnajm@sandia.gov>

Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <cmath>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <unistd.h>
#include <cstdlib>

#include "aev.h"
#include "util.h"
#include "vio.h"
#include "config.h"

bool verbose = false;
//================================================================================
void aev::get_aev_tags(int& rad_dout, int& ang_dout, std::vector<std::string>& aev_tag){
	rad_dout = rad_par.size() * n_rad;
	ang_dout = ang_par.size() * n_ang;
	aev_tag  = tag;
	return;
}
//================================================================================
aev::aev(const std::vector<std::string>& atom_types, const int& nrho_rad, 
         const int& nrho_ang, const int& nalpha, const std::vector<double>& R_c_in) {

	types = atom_types;
	R_c   = R_c_in;

	// spec for radial SFs
	double drho_rad  = R_c[0]/static_cast<double>(nrho_rad);
	double delta_rad = (2./3.)*drho_rad;
	double eta_rad   = 1.0/pow(delta_rad,2);    // a single eta
	std::vector<double> rhov_rad = linspace(drho_rad/2.0,R_c[0]-drho_rad/2.0,nrho_rad);
	std::vector<double> etav_rad = {eta_rad};
	rad_par = two_d_2col_pack(etav_rad,rhov_rad);

	if (verbose){
		std::cout << "atom_types: " << atom_types << std::endl;
		std::cout << "nrho_rad:   " << nrho_rad << std::endl;
		std::cout << "nrho_ang:   " << nrho_ang << std::endl;
		std::cout << "nalpha:     " << nalpha << std::endl;
		std::cout << "R_c:        " << R_c << std::endl;
		std::cout << "drho_rad:   " << drho_rad << std::endl;
		std::cout << "delta_rad:  " << delta_rad << std::endl;
		std::cout << "eta_rad:    " << eta_rad << std::endl; 
		std::cout << "rahov_rad:  " << rhov_rad << std::endl;
		std::cout << "etav_rad:   " << etav_rad << std::endl;
		std::cout << "rad_par:  \n" << rad_par << std::endl;
	}

	//spec for angular SFs
	double drho_ang  = R_c[1]/static_cast<double>(nrho_ang);
	double delta_ang = (2./3.)*drho_ang;
	double eta_ang   = 1.0/pow(delta_ang,2);   // a single eta

	std::vector<double> rhov_ang = linspace(drho_ang/2.0,R_c[1]-drho_ang/2.0,nrho_ang);
	std::vector<double> etav_ang = {eta_ang};
	std::vector<double> alphav   = linspace(0.0,pi,nalpha);

	double zeta               = 8.0;               // a single zeta
	std::vector<double> zetav = {zeta};

	ang_par = two_d_4col_pack(etav_ang, rhov_ang, zetav, alphav);

	if (verbose){
		std::cout << "drho_ang:   " << drho_ang << std::endl;
		std::cout << "delta_ang:  " << delta_ang << std::endl;
		std::cout << "eta_ang:    " << eta_ang << std::endl;
		std::cout << "rhov_ang:   " << rhov_ang << std::endl;
		std::cout << "etav_ang:   " << etav_ang << std::endl;
		std::cout << "alphav:     " << alphav << std::endl;
		std::cout << "zetav:      " << zetav << std::endl;
		std::cout << "ang_par:  \n" << ang_par << std::endl; 
	}

	int nr  = rad_par.size();
	int na  = ang_par.size();
	int n   = atom_types.size();
	n_rad   = n;
	n_ang   = std::round(n*(n+1)/2);
	dout    = nr * n_rad + na * n_ang;
	rad_typ = atom_types;
	ang_typ = std::vector<std::vector<std::string>> (n_ang,std::vector<std::string>(2));
	for (int i=0, k=0; i<n; ++i)
		for (int j=i; j<n; ++j)
			ang_typ[k++] = std::vector<std::string>({atom_types[i],atom_types[j]});
	tag     = std::vector<std::string>(dout," ");

/*
	n is the number of atom types
	nr is the number of radial SFs (# of pairs of parameters in rad_par)
	na is the number of angular SFs (# of quad-tuples of parameters in ang_par)
	n_rad is the number of radial blocks in the AEV, each block corresponding to an atom type
			   and each block being composed of nr elements
	n_ang is the number of angular blocks in the AEV, each block corresponding to a pair of atom types
		       irrespective of ordering, and each block being composed of na elements
	dout is the full size of the AEV
	rad_typ is the list of atom types in the radial AEV components, a direct copy of atom_types
	ang_typ is the list of atom type pairs in the angular AEV components
	tag is primarily for diagnostic purposes
		* a tag of "H0" indicates that the corresponding element of the AEV is a radial element
		(because only one atom type is included in the tag), that it's the radial element for atom type H
		and that it's the 0-th element in the radial H block. We would have elements H0, H1, ..., H<nr-1>
		* a tag of "HC0" indicates that the corresponding element of the AEV is an angular element
		(because a pair of atom types is included in the tag), that it's the angular element for atom types H and C in any order
		and that it's the 0-th element in the angular block. We would have elemtns HC0, HC1, ..., HC<na-1>
*/

	for (int kk=0; kk<n_rad; ++kk){	      // loop over radial AEV components
		int kknr = kk*nr;
		for (int l=0; l<nr; ++l)
			tag[kknr+l] = atom_types[kk] + std::to_string(l);
	}

	int q_offset = nr * n_rad;
	for (int kk=0; kk<n_ang; ++kk){	// loop over angular AEV components
		int kkna = kk*na;
		for (int l=0; l<na; ++l)
			tag[q_offset+kkna+l] = ang_typ[kk][0] + ang_typ[kk][1] + std::to_string(l);
	}

	if (verbose){
		std::cout << "nr:         " << nr << std::endl;
		std::cout << "na:         " << na << std::endl;
		std::cout << "n:          " << n << std::endl;
		std::cout << "n_rad:      " << n_rad << std::endl;
		std::cout << "n_ang:      " << n_ang << std::endl;
		std::cout << "dout:       " << dout << std::endl;
		std::cout << "rad_typ:    " << rad_typ << std::endl;
		std::cout << "ang_typ:  \n" << ang_typ << std::endl; 
		std::cout << "tag:      \n" << tag << std::endl; 
	}

}

//================================================================================
