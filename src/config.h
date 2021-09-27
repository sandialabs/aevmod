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

#ifndef CONFIG_H
#define CONFIG_H

//===============================================================
class config {

//--------------------------------------------------------------------------------
private:

	size_t num_atoms;
	size_t num_coord;
	std::vector<std::string> symb;
	std::vector<std::vector<double>> xyz;
	std::vector<std::vector<int>> radial_index_set;
	std::vector<std::vector<std::vector<int>>> angular_index_set;

//--------------------------------------------------------------------------------
public:

	std::string spname;

	config(const std::vector<std::string>& symb_in);

	bool chk_symb(const std::vector<std::string>& tsymb);

	size_t clear_and_add_structure(const std::vector<double>& xyz_in);
	size_t add_structure(const std::vector<double>& xyz_in);
	std::vector<std::string> get_symb();

	std::vector<std::vector<double>> get_structures();

	size_t clear_and_add_structures(const std::vector<std::vector<double>>& xyz_in);
	size_t add_structures(const std::vector<std::vector<double>>& xyz_in);

	void get_index_sets(std::vector<std::vector<int>>& rad, std::vector<std::vector<std::vector<int>>>& ang);

	std::vector<std::vector<int>> get_radial_index_set();

    std::vector<std::vector<std::vector<int>>> get_angular_index_set();

	void delete_all_structures();

  friend class aev;

//--------------------------------------------------------------------------------
};


#endif