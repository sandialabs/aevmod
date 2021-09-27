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
#include <functional>
#include <cmath>
#include <numeric>

#include "aev.h"
#include "config.h"
#include "vio.h"
#include "util.h"

static bool diagnose = false;
static bool verbose  = false;

//==============================================================================================================
// kernel functions for AEV
// cutoff function
double aev::fc(const double& Rij, const double& Rc){
 	double fcv;
	if (Rij <= Rc)
		fcv = 0.5 * cos(pi * Rij/Rc) + 0.5;
	else
		fcv = 0.0;
	return fcv;
}

// kernel for Radial SFs
double aev::rad_kern(const std::vector<double>& par, const double& Rij, const double& Rc){
	double eta = par[0];
	double rho = par[1];
	double rk  = exp( - eta * pow(Rij-rho,2) ) * fc(Rij,Rc);
	return rk;
}

// kernel for Angular SFs
double aev::ang_kern(const std::vector<double>& par, const double& Rij, const double& Rik, const double& theta, const double& Rc){
	double eta   = par[0];
	double rho   = par[1];
	double zeta  = par[2];
	double alpha = par[3];
	double ra    = 0.5 * (Rij+Rik);
	double rk    = exp ( - eta * pow(ra-rho,2) ) * fc(Rij,Rc) * fc(Rik,Rc);
	double ak    = pow(0.5 + 0.5 * cos(theta-alpha), zeta) * rk;
	return ak;
}

void aev::write_aev_to_file(const std::string& fname){
    std::ofstream ofs;
	ofs.open (fname);  
	for (auto yp : val_aev ) {
	  for (auto ypa : yp) {
	    for (auto ypar : ypa)
	      ofs << std::scientific << std::setprecision(14) << std::setw(24) << ypar << " ";
	    ofs << std::endl;
	  }
	  ofs << std::endl;
	}
	ofs.close();
}

void aev::write_jac_to_file(const std::string& fname){
    std::ofstream ofs;
    ofs.open (fname);  
    for (auto Jp : val_jac) {
      for (auto Jpa : Jp){
          for (auto Jpar : Jpa){
              for (auto Jparl : Jpar)
                  ofs << std::scientific << std::setprecision(14) << std::setw(24) << Jparl << " ";
              ofs << std::endl;
          }
          ofs << std::endl;
      }
      ofs << std::endl;
    }
    ofs.close();
}

//==============================================================================================================
// evaluate the aev for npt structures, all being the same config defined by symb in diff geometries
// returns yaev being a 3d vector [ p ] [ i ] [ q ], where 
// p is the structure index 
// i is the atom index within the config
// q is the index pointing to each AEV entry for that atom i

std::vector<std::vector<std::vector<double>>> aev::eval(config& conf){

	const std::vector<std::string>                          symb = conf.symb; 
	const std::vector<std::vector<int>>              rad_ind_set = conf.radial_index_set; 
	const std::vector<std::vector<std::vector<int>>> ang_ind_set = conf.angular_index_set; 
	const std::vector<std::vector<double>>                  xmat = conf.xyz;

	int N   = symb.size();
	int npt = xmat.size();
	//int din = 3*N;
	int n   = types.size();
	std::vector<std::vector<std::vector<double>>> y(npt, std::vector<std::vector<double>>(N, std::vector<double>(dout,0.0)));

	for (int p=0; p<npt; ++p){   // loop over data points
	  if (diagnose)
	    std::cout << "=====================================================\n"
	              << "Data point p:" << p << std::endl;
	  std::vector<std::vector<double>> x = to_2d(xmat[p],3);

	  for (int i=0; i<N; ++i){	 // loop over atoms in configuration
	  	std::string atom = symb[i];
	  	if (diagnose)
	      std::cout << "=====================================================\n"
	                << "centered on atom # i: " << i << ", atom:" << atom << std::endl;
	  	// deal with radial SFs
	  	int nst = rad_ind_set.size();
	  	int nsf = rad_par.size();		// # radial SFs
	  	if (diagnose)
	  	  std::cout << "Evaluating " << nsf << " Radial Symmetry Functions\n" 
	  	            << "Looping over " << nst << " unitary index sets, one for each atom type\n";
 		for (int kk=0; kk<nst; ++kk){   // loop over unitary index sets
 			std::string elmt   = rad_typ[kk];
 			std::vector<int> indset = rad_ind_set[kk];
 			int kknsf  = kk*nsf;
 		  	if (diagnose)
 				std::cout << "Type: " << elmt << ": Index Set: " << indset << "\n"
					  << "Starting index in radial block: " << kknsf << std::endl;
			for (int j : indset){	// loop over indices in index set if any
			  if (j != i){			// ignore case when j=i
			  	if (diagnose) 
			  	  std::cout << "working on atom # j: " << j << std::endl;
			  	double Rij = l2_distance(x[j],x[i]);    // ||x_j-x_i|| distance
			  	for (int l=0; l<nsf; ++l){              // loop over radial SFs
			  		int q = kknsf + l;
			  		std::vector<double> par = rad_par[l];
			  		y[p][i][q] += rad_kern(par,Rij,R_c[0]);
				  	if (diagnose) 
					  std::cout << "p: " << p << ", i: " << i << ", atom: " << atom << ", type: " << elmt
					            << ", j: " << j << ", q: " << q << ", y: " << y[p][i][q] << std::endl;
			  	}
			  }
			}
 		}
 		int q_offset = nsf * n;

 		// deal with angular SFs
		nst = ang_ind_set.size();
		nsf = ang_par.size();		// nsf is # angular SFs
		if (diagnose)
			std::cout << "Evaluating" << nsf << "Angular Symmetry Functions" 
			          << "Looping over" << nst << "pairwise index sets, one for each unordered pair of atom types\n";
		for (int kk=0; kk<nst; ++kk){  		        // loop over pairwise index sets
			std::string elmt1 = ang_typ[kk][0];		// type 1
			std::string elmt2 = ang_typ[kk][1];		// type 2
			std::vector<std::vector<int>> indset = ang_ind_set[kk];		// index set
			int kknsf = kk*nsf;
			if (diagnose)
				std::cout << "Unordered pair of types: " << elmt1 << ", " << elmt2 << ": Index Set: " << indset 
				          << "Starting index in angular block: " << kknsf << std::endl;
			for (auto pair : indset){			// loop over pairwise index set
				int j=pair[0];
				int k=pair[1];
				if (j != i && k != i) {			// ignore when either j=i or k=i
					if (diagnose)
						std::cout << "Working on atom #s j: " << j << ", k: " << k << std::endl;
					std::vector<double> vij = vdiff(x[j],x[i]);	// x_j - x_i
					std::vector<double> vik = vdiff(x[k],x[i]);	// x_k - x_i
					double Rij = l2_length(vij);      // ||x_j-x_i||
					double Rik = l2_length(vik);      // ||x_k-x_i||
					double vdv = std::inner_product(vij.begin(), vij.end(), vik.begin(), 0.0);
					double theta = acos( std::max(-1.0, std::min( vdv/(Rij*Rik), 1.0)) );   // angle theta
					for (int l=0; l<nsf; ++l){	// loop over angular SFs
						int q = q_offset + kknsf + l;
						std::vector<double> par = ang_par[l];
						y[p][i][q] += ang_kern(par,Rij,Rik,theta,R_c[1]);
					  	if (diagnose) 
						  std::cout << "p: " << p << ", i: " << i << ", atom: " << atom << ", pair: " << elmt1 << ", " << elmt2
						            << ", j: " << j << ", k: " << k << ", q: " << q << ", y: " << y[p][i][q] << std::endl;
					}
				}
			}
		}
	  }
	}
	if (verbose)
		std::cout << y;

	val_aev = y;

	return y;
}
//==============================================================================================================
//==============================================================================================================
// returns, for each structure p, and each atom i, the Jacobian of the AEV w.r.t. the xyz coordinates of all atoms in the config
std::vector<std::vector<std::vector<std::vector<double>>>> aev::eval_Jac (config& conf){

	const std::vector<std::string>                          symb = conf.symb; 
	const std::vector<std::vector<int>>              rad_ind_set = conf.radial_index_set; 
	const std::vector<std::vector<std::vector<int>>> ang_ind_set = conf.angular_index_set; 
	const std::vector<std::vector<double>>                  xmat = conf.xyz;

	// npt is the number of structures in xmat
	int npt = xmat.size();

	// N is the number of atoms in the config
	int N   = symb.size();

	// din is the number of xyz dimensions for the config
	int din = 3*N;

	// n is the number of atom types in the system 
	int n   = types.size();

	// N atoms, involving n atom-types, we have                        (:= N in aev document)
	// n_rad radial index sets := rad_ind_set.size()  := n             (:= n in aev document)
	// n_ang angular index sets := ang_ind_set.size() := n(n+1)/2      (:= m in the aev document)
	// The number of parameter combinations in the radial AEV, per radial index set, is rad_par.size()  (:= M_R in the aev document)
	// The number of parameter combinations in the angular AEV, per angular index set, is ang_par.size() (:= M_A in aev document)
	// Length of AEV, dout = M_R * n_rad + M_A * n_ang;     

	// dout is the dimension of the AEV for any atom in the config
	// jac(p,i,r,s): Jac for structure p=0,...,npt , atom i=0,...,N , AEV component r=0,...,dout ; xyz-coordinate s=0,1,2,....,din=3N
	
	std::vector<std::vector<std::vector<std::vector<double>>>> 
		jac (npt, std::vector<std::vector<std::vector<double>>>(N, std::vector<std::vector<double>>(dout,std::vector<double> (din,0.0))));

	for (int p=0; p<npt; ++p){   // loop over data points
	  if (diagnose)
	    std::cout << "=====================================================\n"
	              << "Data point p:" << p << std::endl;
	  std::vector<std::vector<double>> x = to_2d(xmat[p],3);

	  for (int i=0; i<N; ++i){	 // loop over atoms in configuration
	  	std::string atom = symb[i];
	  	if (diagnose)
	      std::cout << "=====================================================\n"
	                << "centered on atom # i: " << i << ", atom:" << atom << std::endl;
	  	// deal with radial SFs
	  	int nst = rad_ind_set.size();   // this is n
	  	int nsf = rad_par.size();		// # radial SFs (this is M_R)
	  	if (diagnose)
	  	  std::cout << "Evaluating " << nsf << " Radial Symmetry Functions\n" 
	  	            << "Looping over " << nst << " unitary index sets, one for each atom type\n";
	  	// ordering follows : 
	  	//   outer loop over the n_rad=n radial index sets, kk=0,...,n-1
	  	//        inner loop (index l below) over M_R radial SFs, i.e. over parameter-combinations, in the index set rad_ind_set[kk], 
	  	// thus the innermost index is q = kk*nsf + l, where nsf = rad_par.size() = # radial SFs = M_R

	  	// the loop over j : indset does the sum over atoms j contributions to the AEV sum for current atom i
	  	// loop over radial (unitary) index sets, i.e. over element types 0,1,...,n-1, being the n high level blocks of the radial aev

 		for (int kk=0; kk<nst; ++kk){   
 		  	std::string elmt   = rad_typ[kk];           // specific element type for this radial index set
 		  	std::vector<int> indset = rad_ind_set[kk];  // specific radial index set
 		  	int kknsf  = kk*nsf;                        // starting index in radial block

 		  	if (diagnose)
 				std::cout << "Type: " << elmt << ": Index Set: " << indset << "\n"
					  << "Starting index in radial block: " << kknsf << std::endl;

			for (int j : indset){	                  // loop over indices in index set if any
			  	if (j != i){			// ignore case when j=i
			  		if (diagnose) 
			  	  		std::cout << "working on atom # j: " << j << std::endl;

				  	double Rij = l2_distance(x[j],x[i]);    // ||x_j-x_i|| distance
					double Rc  = R_c[0];

					if (Rij <= Rc) {                  // If false, all grad contribs are zero below ... skip

						double pioRc    = pi/Rc;
						double oRij     = 1.0/Rij;
						double dfccon   = -0.5*pioRc*sin(Rij*pioRc)*oRij;

						std::vector<double> dRijdi(3);
						for (int kc=0; kc<3; ++kc)
							dRijdi[kc] = (x[i][kc]-x[j][kc]) * oRij;

						std::vector<double> dfcdi(3);
						for (int kc=0; kc<3; ++kc)
							dfcdi[kc] = dfccon * (x[i][kc]-x[j][kc]);

					  	// inner most loop over radial AEV components, q is the index: 0, 1, ..., M_R*n 
					  	for (int l=0; l<nsf; ++l){              // loop over radial SFs

					  		int q = kknsf + l;
					  		std::vector<double> par = rad_par[l];

						  	if (diagnose) 
							  std::cout << "p: " << p << ", i: " << i << ", atom: " << atom << ", type: " << elmt
							            << ", j: " << j << ", q: " << q << std::endl;

							// now loop over 3N xyz coordinates. we do this by looping over the N atoms, then over the 3 coords for each
							// mapping to the aev document for Jac computation (eq. 23)
		                    //       index type                         Here                     AEV Jac document
							//  outer atom index                          i                            i
						    //  elmt type or radial index set index       kk                           p
						    //  rad_par, radial aev param combo index     l                            q
						    //  inner atom index in AEV sum               j                            j 
							//  inner-most AEV index                      q                            r
						    //  inner atom index for Jac computations     t                            l 
						    // 
							// dout is the dimension of the AEV for any atom in the config
							// jac(p,i,q,s): Jac for structure p=0,...,npt , atom i=0,...,N , AEV component q=0,...,dout ; xyz-coordinate s=0,1,2,....,din=3N
						    // s = 0,1,2 , 3,4,5 , ...     being composed of N triplets indexed by t=0,1,...,

							// prepare terms:

							double eta = par[0];
							double rho = par[1];
							double dRr = Rij - rho;
							double rke = exp( - eta * pow(dRr,2) );
							double ctt = -2.0 * eta * dRr * fc(Rij,Rc);

							// only t=i and t=j have non-zero contributions 

							// contributions for t=i
							int s = i*3;
							for (int kc=0; kc<3; ++kc)
								jac[p][i][q][s+kc] += rke * ( ctt * dRijdi[kc] + dfcdi[kc] );

							// contributions for t=j
							// nb. dRijdj[]=-dRijdi[] and dfcdj[]=-dfcdi[]
							s = j*3;
							for (int kc=0; kc<3; ++kc)
								jac[p][i][q][s+kc] += -rke * ( ctt * dRijdi[kc] + dfcdi[kc] );

						}
				  	}

			  	}
			}
 		}

 		int q_offset = nsf * n;

 		// deal with angular SFs
		nst = ang_ind_set.size();
		nsf = ang_par.size();		// nsf is # angular SFs (this is M_A)
		if (diagnose)
			std::cout << "Evaluating" << nsf << "Angular Symmetry Functions" 
			          << "Looping over" << nst << "pairwise index sets, one for each unordered pair of atom types\n";
		for (int kk=0; kk<nst; ++kk){  		        // loop over pairwise index sets
			std::string elmt1 = ang_typ[kk][0];		// type 1
			std::string elmt2 = ang_typ[kk][1];		// type 2
			std::vector<std::vector<int>> indset = ang_ind_set[kk];		// index set
			int kknsf = kk*nsf;
			if (diagnose)
				std::cout << "Unordered pair of types: " << elmt1 << ", " << elmt2 << ": Index Set: " << indset 
				          << "Starting index in angular block: " << kknsf << std::endl;
			for (auto pair : indset){			// loop over pairwise index set
				int j=pair[0];
				int k=pair[1];
				if (j != i && k != i) {			// ignore when either j=i or k=i
					if (diagnose)
						std::cout << "Working on atom #s j: " << j << ", k: " << k << std::endl;
					std::vector<double> vij = vdiff(x[j],x[i]);	// x_j - x_i
					std::vector<double> vik = vdiff(x[k],x[i]);	// x_k - x_i
					double Rij = l2_length(vij);      // ||x_j-x_i||
					double Rik = l2_length(vik);      // ||x_k-x_i||
					double Rc  = R_c[1];

					if ( Rij <= Rc && Rik <= Rc) {  // If false, all grad contribs are zero below ... skip

						double ra       = 0.5 * (Rij+Rik);
						double pioRc    = pi/Rc;

						double oRij     = 1.0/Rij;
						double dfcconj  = -0.5*pioRc*sin(Rij*pioRc)*oRij;

						double oRik     = 1.0/Rik;
						double dfcconk  = -0.5*pioRc*sin(Rik*pioRc)*oRik;

						double vdv   = std::inner_product(vij.begin(), vij.end(), vik.begin(), 0.0);
						double theta = acos( std::max(-1.0, std::min( vdv/(Rij*Rik), 1.0)) );   // angle theta

						double costh = cos(theta);
						double ctden = Rij * Rik * sin(theta);

						double fcij  = fc(Rij,Rc);
						double fcik  = fc(Rik,Rc);
						double fcijk = fcij * fcik;

						// fcij block
						std::vector<double> dRijdi(3);
						for (int kc=0; kc<3; ++kc)
							dRijdi[kc] = (x[i][kc]-x[j][kc]) * oRij;

						std::vector<double> dfcijdi(3);
						for (int kc=0; kc<3; ++kc)
							dfcijdi[kc] = dfcconj * (x[i][kc]-x[j][kc]);

						//dRijdj = -dRijdi and dfcijdj = -dfijcdi
						//dRijdk = 0 and dfcijdk = 0

						// fcik block
						std::vector<double> dRikdi(3);
						for (int kc=0; kc<3; ++kc)
							dRikdi[kc] = (x[i][kc]-x[k][kc]) * oRik;

						std::vector<double> dfcikdi(3);
						for (int kc=0; kc<3; ++kc)
							dfcikdi[kc] = dfcconk * (x[i][kc]-x[k][kc]);

						//dRikdk = -dRikdi and dfcikdk = -dfcikcdi
						//dRikdj = 0 and dfcikdj = 0

						// derivatives of (R_ij * R_ik)
						std::vector<double> dRRdi(3);
						for (int kc=0; kc<3; ++kc)
							dRRdi[kc] = Rik * dRijdi[kc] + Rij * dRikdi[kc];

						std::vector<double> dRRdj(3);
						for (int kc=0; kc<3; ++kc)
							dRRdj[kc] = - Rik * dRijdi[kc];

						std::vector<double> dRRdk(3);
						for (int kc=0; kc<3; ++kc)
							dRRdk[kc] = - Rij * dRikdi[kc];

						// derivatives of x_ij.x_ik
						std::vector<double> dxxdi(3);
						for (int kc=0; kc<3; ++kc)
							dxxdi[kc] = -(x[k][kc]-x[i][kc]+x[j][kc]-x[i][kc]);

						std::vector<double> dxxdj(3);
						for (int kc=0; kc<3; ++kc)
							dxxdj[kc] = x[k][kc]-x[i][kc];

						std::vector<double> dxxdk(3);
						for (int kc=0; kc<3; ++kc)
							dxxdk[kc] = x[j][kc]-x[i][kc];

						// now loop over angular SFs
						for (int l=0; l<nsf; ++l){

							int q = q_offset + kknsf + l;
							std::vector<double> par = ang_par[l];

						  	if (diagnose) 
							  std::cout << "p: " << p << ", i: " << i << ", atom: " << atom << ", pair: " << elmt1 << ", " << elmt2
							            << ", j: " << j << ", k: " << k << ", q: " << q << std::endl;

							double eta   = par[0];
							double rho   = par[1];
							double zeta  = par[2];
							double alpha = par[3];
							double costa = cos(theta-alpha);
							double carga = 0.5 + 0.5 * costa;
							double cta   = pow(carga, zeta);
							double ctda  = -0.5 * zeta * (cta / carga) * sin(theta-alpha);
							double ctb   = exp ( -eta * pow(ra-rho,2) );
							double ctdb  = -eta * ctb * (ra - rho);

							double ctab  = cta  * ctb;
							double ct1   = ctda * ctb  * fcijk;
							double ct2   = cta  * ctdb * fcijk;
							double ct3   = ctab * fcik;
							double ct4   = ctab * fcij;

							// only components 
							// now loop over 3N xyz coordinates. we do this by looping over the N atoms, then over the 3 coords for each
							// mapping to the aev document for Jac computation (eq. 24)
		                    //       index type                         Here                     AEV Jac document
							//  outer atom index                          i                            i
						    //  elmt type or radial index set index       kk                           p
						    //  rad_par, radial aev param combo index     l                            q
						    //  inner atom index in AEV sum               j                            j 
						    //  inner atom index in AEV sum               k                            k 
							//  inner-most AEV index                      q                            r
						    //  inner atom index for Jac computations     t                            l 
						    // 
							// dout is the dimension of the AEV for any atom in the config
							// jac(p,i,q,s): Jac for structure p=0,...,npt , atom i=0,...,N , AEV component q=0,...,dout ; xyz-coordinate s=0,1,2,....,din=3N
						    // s = 0,1,2 , 3,4,5 , ...     being composed of N triplets indexed by t=0,1,...,

							std::vector<double> dthetadi(3);
							for (int kc=0; kc<3; ++kc)
								dthetadi[kc] = ( - dxxdi[kc] + costh * dRRdi[kc] ) / ctden;

							std::vector<double> dthetadj(3);
							for (int kc=0; kc<3; ++kc)
								dthetadj[kc] = ( - dxxdj[kc] + costh * dRRdj[kc] ) / ctden;

							std::vector<double> dthetadk(3);
							for (int kc=0; kc<3; ++kc)
								dthetadk[kc] = ( - dxxdk[kc] + costh * dRRdk[kc] ) / ctden;

							// t = i
							int s = i*3;
							for (int kc=0; kc<3; ++kc)
								jac[p][i][q][s+kc] += ct1 * dthetadi[kc] + ct2 * (dRijdi[kc] + dRikdi[kc]) + 
													  ct3 * dfcijdi[kc]  + ct4 * dfcikdi[kc];

							// t = j
							s = j*3;
							for (int kc=0; kc<3; ++kc)
								jac[p][i][q][s+kc] += ct1 * dthetadj[kc] - ct2 * dRijdi[kc] - ct3 * dfcijdi[kc];

							// t = k
							s = k*3;
							for (int kc=0; kc<3; ++kc)
								jac[p][i][q][s+kc] += ct1 * dthetadk[kc] - ct2 * dRikdi[kc] - ct4 * dfcikdi[kc];


						}
					}
				}
			}
		}

	  }
	}

	val_jac = jac;

	return jac;
}

