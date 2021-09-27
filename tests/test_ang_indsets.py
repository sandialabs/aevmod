"""
=====================================================================================
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
=====================================================================================
"""

import sys
import numpy as np
import math
import aevmod
import pytest
import utils

# ====================================================================================================
def test_ang_indsets():

	types = ['C','H']
	myaev = aevmod.aev(types,8,4,4)

	symb, vxyz = utils.read_xyz("tests/xyz_c2h5.txt")
	cnf  = aevmod.config(symb)
	myaev.build_index_sets(cnf)

	ais = utils.read_ang_indsets("tests/ais_c2h5_8_4_4.txt")
	assert cnf.get_angular_index_set() == ais

# ====================================================================================================
