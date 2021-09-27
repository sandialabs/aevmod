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
#include "Sacado.hpp"

static bool diagnose = false;
static bool verbose  = false;

template<typename ValueType>
ValueType fc_sac(const ValueType& Rij, const double& Rc);
template<typename ValueType>
ValueType rad_kern_sac(const std::vector<double>& par, const ValueType& Rij, const double& Rc);
template<typename ValueType>
ValueType ang_kern_sac(const std::vector<double>& par, const ValueType& Rij, const ValueType& Rik, const ValueType& theta, const double& Rc);
template<typename ValueType>
ValueType l2_distance_sac(const std::vector<ValueType>& a, const std::vector<ValueType>& b);
template<typename ValueType>
std::vector<ValueType> vdiff_sac(const std::vector<ValueType>& a, const std::vector<ValueType>& b);
template<typename ValueType>
ValueType l2_length_sac(const std::vector<ValueType>& v);
//==============================================================================================================
template<typename ValueType>
std::vector<ValueType> vdiff_sac(const std::vector<ValueType>& a, const std::vector<ValueType>& b){
	std::vector<ValueType> amb(a.size());
	std::transform(a.begin(), a.end(), b.begin(), amb.begin(), std::minus<ValueType>());
	return amb;
}
//==============================================================================================================
template<typename ValueType>
ValueType l2_length_sac(const std::vector<ValueType>& v){
	return std::sqrt(std::inner_product(v.begin(), v.end(), v.begin(), static_cast<ValueType>(0.0)));
}
//==============================================================================================================
template<typename ValueType>
ValueType l2_distance_sac(const std::vector<ValueType>& a, const std::vector<ValueType>& b){
	std::vector<ValueType> d = vdiff_sac(a,b);
	return std::sqrt(std::inner_product(d.begin(), d.end(), d.begin(), static_cast<ValueType>(0.0)));
}
//==============================================================================================================
// incoming with xin being xmat[p], a 3*N long 1d vector of xyz coords of this structure p
// call this separately for each structure

template<typename ValueType>
std::vector<ValueType> aev::evaluateFunction(const std::vector<ValueType> &xin, const int p, const int i, config& conf) {

	const std::vector<std::string>                          symb = conf.symb; 
	const std::vector<std::vector<int>>              rad_ind_set = conf.radial_index_set; 
	const std::vector<std::vector<std::vector<int>>> ang_ind_set = conf.angular_index_set; 

	int N   = symb.size();
	int din = 3*N;
	int n   = types.size();
	std::vector<ValueType> y(dout,static_cast<ValueType>(0.0));

#ifdef VERBOSE
  std::cout << "[===============evaluateFunction====================\n";
  for (auto xv : xin){std::cout << "x:" << xv << std::endl;}
#endif

  std::vector<std::vector<ValueType>> x = to_2d(xin,3);

  	std::string atom = symb[i];

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
			  	ValueType Rij = l2_distance_sac(x[j],x[i]);    // ||x_j-x_i|| distance
			  	for (int l=0; l<nsf && Rij <= R_c[0]; ++l){              // loop over radial SFs if Rij <= R_c[0]
			  		int q = kknsf + l;
			  		std::vector<double> par = rad_par[l];
			  		y[q] += rad_kern_sac(par,Rij,R_c[0]);
				  	if (diagnose) 
					  std::cout << "p: " << p << ", i: " << i << ", atom: " << atom << ", type: " << elmt
					            << ", j: " << j << ", q: " << q << ", y: " << y[q] << std::endl;
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
					std::vector<ValueType> vij = vdiff_sac(x[j],x[i]);	// x_j - x_i
					std::vector<ValueType> vik = vdiff_sac(x[k],x[i]);	// x_k - x_i
					ValueType Rij = l2_length_sac(vij);      // ||x_j-x_i||
					ValueType Rik = l2_length_sac(vik);      // ||x_k-x_i||
					ValueType vdv = std::inner_product(vij.begin(), vij.end(), vik.begin(), static_cast<ValueType>(0.0));
					ValueType theta = acos( std::max(-static_cast<ValueType>(1.0), std::min( vdv/(Rij*Rik), static_cast<ValueType>(1.0))) );   // angle theta
					for (int l=0; l<nsf && Rij <= R_c[1] && Rik <= R_c[1]; ++l){	// loop over angular SFs if Rij and Rik both <= R_c[1]
						int q = q_offset + kknsf + l;
						std::vector<double> par = ang_par[l];
						y[q] += ang_kern_sac(par,Rij,Rik,theta,R_c[1]);
					  	if (diagnose) 
						  std::cout << "p: " << p << ", i: " << i << ", atom: " << atom << ", pair: " << elmt1 << ", " << elmt2
						            << ", j: " << j << ", k: " << k << ", q: " << q << ", y: " << y[q] << std::endl;
					}
				}
			}
		}

	if (verbose)
		std::cout << y;

#ifdef VERBOSE
  for (auto yv : y){std::cout << "y:" << yv << std::endl;}
	std::cout << "================evaluateFunction===================]\n";
#endif

return y;
}

//==============================================================================================================
template<typename ValueType>
std::vector<std::vector<ValueType>> aev::evaluateJacobian(const std::vector<ValueType> &x, const int p, const int ia, config& conf) {

#ifdef VERBOSE
	std::cout << "[===============evaluateJacobian====================\n";
  for (auto xv : x){std::cout << "x:" << xv << std::endl;}
#endif

  using FadType = Sacado::Fad::DFad<ValueType>;
  //using FadType = Sacado::Fad::SLFad<ValueType,100>;   // static version if willing to put hard-wired input space size

  const int m = x.size();  // # input dimensions
  std::vector<FadType> x_fad(m);
  for (int i=0;i<m;++i) 
    x_fad[i] = FadType(m, i, x[i]);

  auto f_fad = evaluateFunction(x_fad,p,ia,conf);
  const int n = f_fad.size();  // # output dimensions

  std::vector<std::vector<ValueType>> g(n,std::vector<ValueType>(m));
  for (int k=0; k<n; ++k)
    for (int i=0; i<m; ++i)
      g[k][i] = f_fad[k].dx(i);

#ifdef VERBOSE
  for (auto gv : g)
    for (auto gvv : gv) {std::cout << "g:" << gvv << std::endl;}
	std::cout << "================evaluateJacobian===================]\n";
#endif
  return g;
}

//==============================================================================================================
template<typename ValueType>
void aev::evaluateHessian(const std::vector<ValueType> &x,
		     std::vector<std::vector<ValueType>> &g,
		     std::vector<std::vector<std::vector<ValueType>>> &h, const int p, const int ia, config& conf) {

#ifdef VERBOSE
	std::cout << "[===============evaluateHessian====================\n";
  for (auto xv : x){std::cout << "x:" << std::scientific << std::setprecision(16) << xv << std::endl;}
#endif

  using FadType = Sacado::Fad::DFad<ValueType>;
  //using FadType = Sacado::Fad::SLFad<ValueType,100>;		// static version if willing to put hard-wired input space size

  const int m = x.size();   // # input dimensions
  std::vector<FadType> x_fad(m);
  for (int i=0;i<m;++i)
    x_fad[i] = FadType(m, i, x[i]);

  auto g_fad = evaluateJacobian(x_fad,p,ia,conf);  // g_fad[k][i]= df_k/dx_i
  const int n = g_fad.size();  // # output dimensions

  g.resize(n,std::vector<ValueType>(m));
  h.resize(n,std::vector<std::vector<ValueType>>(m, std::vector<ValueType>(m)));
  for (int k=0; k<n; ++k){
    for (int i=0; i<m; ++i) {
      g[k][i] = g_fad[k][i].val();
      if (g_fad[k][i].length() > 0)
        for (int j=0; j<m; ++j) 
          h[k][i][j] = g_fad[k][i].fastAccessDx(j);
      else
        for (int j=0; j<m; ++j)
          h[k][i][j] = static_cast<ValueType>(0.0);
    }
  }
#ifdef VERBOSE
  for (auto hv : h)
    for (auto hvv: hv)
      for (auto hvvv: hvv) {std::cout << "h:" << hvvv << std::endl;}
	std::cout << "================evaluateHessian===================]\n";
#endif

}
//==============================================================================================================
// kernel functions for AEV
// cutoff function
template<typename ValueType>
ValueType aev::fc_sac(const ValueType& Rij, const double& Rc){
 	ValueType fcv;
	if (Rij <= Rc)
		fcv = 0.5 * cos(pi * Rij/Rc) + 0.5;
	else
		fcv = 0.0;
	return fcv;
}

// kernel for Radial SFs
template<typename ValueType>
ValueType aev::rad_kern_sac(const std::vector<double>& par, const ValueType& Rij, const double& Rc){
	double eta = par[0];
	double rho = par[1];
	ValueType rk  = exp( - eta * pow(Rij-rho,2) ) * fc_sac(Rij,Rc);
	return rk;
}

// kernel for Angular SFs
template<typename ValueType>
ValueType aev::ang_kern_sac(const std::vector<double>& par, const ValueType& Rij, const ValueType& Rik, const ValueType& theta, const double& Rc){
	double eta   = par[0];
	double rho   = par[1];
	double zeta  = par[2];
	double alpha = par[3];
	ValueType ra    = 0.5 * (Rij+Rik);
	ValueType rk    = exp ( - eta * pow(ra-rho,2) ) * fc_sac(Rij,Rc) * fc_sac(Rik,Rc);
	ValueType ak    = pow(0.5 + 0.5 * cos(theta-alpha), zeta) * rk;
	return ak;
}


//==============================================================================================================
// evaluate the aev for npt structures, all being the same config defined by symb in diff geometries
// returns yaev being a 3d vector [ p ] [ i ] [ q ], where 
// p is the structure index 
// i is the atom index within the config
// q is the index pointing to each AEV entry for that atom i

std::vector<std::vector<std::vector<double>>> aev::eval_sac(config& conf){

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
		  for (int i=0; i<N; ++i){	 // loop over atoms in configuration

	  	std::string atom = symb[i];
  		if (diagnose)
      		std::cout << "=====================================================\n"
                << "centered on atom # i: " << i << ", atom:" << atom << std::endl;

	  	y[p][i] = evaluateFunction(xmat[p],p,i,conf);
	}
	}
	val_aev = y;

#ifdef VERBOSE
  for (auto yv : y) for (auto yvv: yv) {std::cout << "y:" << yvv << std::endl;}
	std::cout << "================eval_sac===================]\n";
#endif

	return y;
}
//==============================================================================================================
//==============================================================================================================
//==============================================================================================================
// returns, for each structure p, and each atom i, the Jacobian of the AEV w.r.t. the xyz coordinates of all atoms in the config
std::vector<std::vector<std::vector<std::vector<double>>>> aev::eval_Jac_sac (config& conf){

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

	bool Jdiagnose=false;
	for (int p=0; p<npt; ++p){   // loop over data points
	  if (Jdiagnose)
	    std::cout << "=====================================================\n"
	              << "Data point p:" << p << std::endl;
	  std::vector<std::vector<double>> x = to_2d(xmat[p],3);

	  for (int i=0; i<N; ++i){	 // loop over atoms in configuration
	  	std::string atom = symb[i];
	  	if (Jdiagnose)
	      std::cout << "=====================================================\n"
	                << "centered on atom # i: " << i << ", atom:" << atom << std::endl;
		jac[p][i] = evaluateJacobian(xmat[p],p,i,conf);

	  }
	}

	val_jac = jac;

	return jac;
}

//==============================================================================================================
std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> aev::eval_Hess_sac (config& conf){

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
	
	std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> 
		Hess(npt, std::vector<std::vector<std::vector<std::vector<double>>>>(N, std::vector<std::vector<std::vector<double>>>(dout,std::vector<std::vector<double>> (din,std::vector<double>(din,0.0)))));

	bool Hdiagnose = false;
	for (int p=0; p<npt; ++p){   // loop over data points
	  if (Hdiagnose)
	    std::cout << "=====================================================\n"
	              << "Data point p:" << p << std::endl;
	  std::vector<std::vector<double>> x = to_2d(xmat[p],3);

	  for (int i=0; i<N; ++i){	 // loop over atoms in configuration
	  	std::string atom = symb[i];
	  	if (Hdiagnose)
	      std::cout << "=====================================================\n"
	                << "centered on atom # i: " << i << ", atom:" << atom << std::endl;
	    std::vector<std::vector<double>> g;
		evaluateHessian(xmat[p],g,Hess[p][i],p,i,conf);
	  }
	}

	val_Hess = Hess;

	return Hess;
}



//==============================================================================================================
