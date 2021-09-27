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
#include <cmath>
#include <numeric>

#include "config.h"
#include "util.h"
#include "vio.h"

//================================================================================
config::config(const std::vector<std::string>& symb_in){
	
	symb = symb_in;
	
	num_atoms = symb.size();
	num_coord = 3*num_atoms;

	spname = std::accumulate(symb.begin(), symb.end(), std::string(), 
    	[](const std::string& a, const std::string& b) -> std::string { 
        	return a + b; 
    	} ); 
}

//================================================================================================
std::vector<std::string> config::get_symb(){
	return symb;
}

//================================================================================================
std::vector<std::vector<double>> config::get_structures(){
	return xyz;
}

//================================================================================================
void config::delete_all_structures(){
	xyz.clear();
}

//================================================================================================
bool config::chk_symb(const std::vector<std::string>& tsymb){
	return isEqual(symb, tsymb);
}

//================================================================================================
size_t config::clear_and_add_structure(const std::vector<double>& xyz_in){

	if (xyz_in.size() != num_coord)
		{std::cerr << "mismatch of xyz_in with config num_atoms"; exit(-1);}

	xyz.clear();
	xyz.push_back(xyz_in);

	return xyz.size();
}

//================================================================================================
size_t config::add_structure(const std::vector<double>& xyz_in){

	if (xyz_in.size() != num_coord)
		{std::cerr << "mismatch of xyz_in with config num_atoms\n"; exit(-1);}

	xyz.push_back(xyz_in);
	return xyz.size();
}

//================================================================================================
size_t config::clear_and_add_structures(const std::vector<std::vector<double>>& xyz_in){

	for (auto v : xyz_in)
		if (v.size() != num_coord)
			{std::cerr << "mismatch of xyz_in with config num_atoms\n"; exit(-1);}

	xyz.clear();

	xyz.insert(std::end(xyz), std::begin(xyz_in), std::end(xyz_in));
	return xyz.size();
}

//================================================================================================
size_t config::add_structures(const std::vector<std::vector<double>>& xyz_in){

	for (auto v : xyz_in)
		if (v.size() != num_coord)
			{std::cerr << "mismatch of xyz_in with config num_atoms\n"; exit(-1);}

	xyz.insert(std::end(xyz), std::begin(xyz_in), std::end(xyz_in));
	return xyz.size();
}

//================================================================================================
void config::get_index_sets(std::vector<std::vector<int>>& rad, std::vector<std::vector<std::vector<int>>>& ang){
	rad = radial_index_set;
	ang = angular_index_set;	
}

//================================================================================================
std::vector<std::vector<int>> config::get_radial_index_set(){
	return radial_index_set;
}

//================================================================================================
std::vector<std::vector<std::vector<int>>> config::get_angular_index_set(){
	return angular_index_set;	
}

//================================================================================