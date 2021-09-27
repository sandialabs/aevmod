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

#ifndef UTIL_H
#define UTIL_H
//================================================================================
std::vector<double> linspace(const double& a, const double& b, const int& num);
std::vector<std::vector<double>> two_d_2col_pack(const std::vector<double>& a, const std::vector<double>& b);
std::vector<std::vector<double>> two_d_4col_pack(const std::vector<double>& a, const std::vector<double>& b,
                         const std::vector<double>& c, const std::vector<double>& d);
double l2_distance(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> vdiff(const std::vector<double>& a, const std::vector<double>& b);
double l2_length(const std::vector<double>& v);

//================================================================================
template<typename T>
bool isEqual(const std::vector<T>& v1, const std::vector<T>& v2);
#include "isEqual.h"
//================================================================================
template < typename T >
std::vector< std::vector<T> > to_2d( const std::vector<T>& flat_vec, std::size_t ncols );
#include "to_2d.h"
//================================================================================
#endif