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
#include <algorithm>

#include "aev.h"
#include "config.h"
#include "vio.h"

//==============================================================================================================
void aev::bld_index_sets(const std::vector<std::string>& symb, std::vector<std::vector<int>>& rad_ind_set, std::vector<std::vector<std::vector<int>>>& ang_ind_set) {

	// sanity check, check if any atom in the symb list is not among the atom types of this aev
	for (auto key : symb)
		if (std::count(types.begin(), types.end(), key) == 0)
			{std::cerr << "atom type mismatch\n"; exit(EXIT_FAILURE);}

	int N = symb.size();

	// build radial index sets S_\tau
	rad_ind_set = std::vector<std::vector<int>>(n_rad, std::vector<int>());

	for (int i=0; i<N; ++i)               // loop over atoms in conf
		for (int j=0; j<n_rad; ++j)       // loop over elemtns of index set
			if (rad_typ[j] == symb[i])
				rad_ind_set[j].push_back(i);


	// build angular index sets S_{\tau,\kappa}
	ang_ind_set = std::vector<std::vector<std::vector<int>>> (n_ang, std::vector<std::vector<int>>());

	for (int i=0; i<N-1; ++i)     
		for (int j=i+1; j<N; ++j) 
			for (int k=0; k<n_ang; ++k)	{
				std::string tp1 = ang_typ[k][0];
				std::string tp2 = ang_typ[k][1];
				if ( (symb[i] == tp1 && symb[j] == tp2) || (symb[i] == tp2 && symb[j] == tp1)	)
					ang_ind_set[k].push_back(std::vector<int>({i,j}));
			}     

}
//==============================================================================================================
//==============================================================================================================
void aev::build_index_sets(config& conf) {
	bld_index_sets(conf.symb,conf.radial_index_set,conf.angular_index_set);
}
//==============================================================================================================
