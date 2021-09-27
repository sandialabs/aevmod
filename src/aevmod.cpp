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

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
namespace py = pybind11;

PYBIND11_MODULE(aevmod, m) {

    m.doc() = "pybind11 aev module";

    py::class_<aev>(m,"aev")
    .def(py::init<const std::vector<std::string>&, const int&, const int&, const int&, const std::vector<double>&>()
        ,
        py::arg("atom_types")=std::vector<std::string>({"C","H"}),
        py::arg("nrho_rad")=32,
        py::arg("nrho_ang")=8,
        py::arg("nalpha")=8,
        py::arg("R_c_in")=std::vector<double>({4.6,3.1})
        )
    .def_readwrite("dout", &aev::dout)
    .def_readwrite("n_rad", &aev::n_rad)
    .def_readwrite("n_ang", &aev::n_ang)
    .def("build_index_sets",&aev::build_index_sets)
    .def("eval",&aev::eval)
    .def("eval_sac",&aev::eval_sac)
    .def("eval_Jac",&aev::eval_Jac)
    .def("eval_Jac_sac",&aev::eval_Jac_sac)
    .def("eval_Hess_sac",&aev::eval_Hess_sac)
    .def("write_aev_to_file",&aev::write_aev_to_file)
    .def("write_jac_to_file",&aev::write_jac_to_file)
    ;

    py::class_<config>(m,"config")
    .def(py::init<const std::vector<std::string>&>())
    .def("add_structures", &config::add_structures)
    .def("add_structure", &config::add_structure)
    .def("get_structures", &config::get_structures)
    .def("get_radial_index_set", &config::get_radial_index_set)
    .def("get_angular_index_set", &config::get_angular_index_set)
    ;

}

//nb in case of an overloaded function add_structures(), use something like this:
//.def("add_structures", py::overload_cast< const std::vector<std::vector<double>>& >(&config::add_structures))

