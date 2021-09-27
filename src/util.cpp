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
#include <functional>
#include <algorithm>
#include <cmath>
#include <numeric>

#include "aev.h"
#include "util.h"
#include "vio.h"

//==============================================================================================================
std::vector<double> vdiff(const std::vector<double>& a, const std::vector<double>& b){
	std::vector<double> amb(a.size());
	std::transform(a.begin(), a.end(), b.begin(), amb.begin(), std::minus<double>());  
	return amb;
}
//==============================================================================================================
double l2_length(const std::vector<double>& v){
	return std::sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.0));
}
//==============================================================================================================
double l2_distance(const std::vector<double>& a, const std::vector<double>& b){
	std::vector<double> d = vdiff(a,b); 
	return std::sqrt(std::inner_product(d.begin(), d.end(), d.begin(), 0.0));
}
//================================================================================
// Lay a uniform mesh of <num> points on the interval [a,b] and assign them to a vector
// if num > 1, both points <a> and <b> are included in the mesh, and h=(b-a)/(num-1)
// else, if num = 1, the one point requested is returned at <a>

std::vector<double> linspace(const double& a, const double& b, const int& num) {

	if (num < 1)
		{std::cerr << "ERROR: linspace requires num > 0"; exit(-1);}

	std::vector<double> v(num);

	if (num == 1) { 
		v[0] = a;
	} else {
		double h = (b - a)/static_cast<double>(num-1);
		for (int i = 0; i < num; ++i) 
		  v[i] = a + static_cast<double>(i) * h;
	}
	
	return v;
}
//================================================================================
//====================================================================================================
std::vector<std::vector<double>> two_d_2col_pack(const std::vector<double>& a, const std::vector<double>& b) {
/*
	given two 1d vectors, a and b, where a has na elements and b has nb elements
	return a 2d vector that is (na.nb x 2) being the tensor product of the two 1d arrays
	thus, e.g. say na=2, nb=3, then the returned 6x2 2d vector "pack" has the structure:
	a0 b0
	a0 b1
	a0 b2
	a1 b0
	a1 b1
	a1 b2
*/
	std::vector<std::vector<double>> pack(a.size()*b.size(),std::vector<double>(2));
	int k=0;
	for (auto xa : a)
		for (auto xb : b)
			pack[k++] = std::vector<double>({xa, xb});

	return pack;
}
//====================================================================================================
//====================================================================================================
std::vector<std::vector<double>> two_d_4col_pack(const std::vector<double>& a, const std::vector<double>& b,
												 const std::vector<double>& c, const std::vector<double>& d) {
/*
	given 4 1d vectors, a, b, c, d, being na, nb, nc, nd long
	return a 2d vector that is (na.nb.nc.nd x 2), thus, say na=2, nb=1, nc=3, nd=2, get:
	a0 b0 c0 d0
	a0 b0 c0 d1
	a0 b0 c1 d0
	a0 b0 c0 d1
	a0 b0 c2 d0
	a0 b0 c2 d1
	a1 b0 c0 d0
	a1 b0 c0 d1
	a1 b0 c1 d0
	a1 b0 c0 d1
	a1 b0 c2 d0
	a1 b0 c2 d1
*/
	std::vector<std::vector<double>> pack(a.size()*b.size()*c.size()*d.size(),std::vector<double>(2));
	int k=0;
	for (auto xa : a)
		for (auto xb : b)
			for (auto xc : c)
				for (auto xd : d)
					pack[k++] = std::vector<double>({xa, xb, xc, xd});

	return pack;
}
//====================================================================================================

