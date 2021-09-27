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

# nb. source code must be called test_base.py

import sys
import numpy as np
import math
import aevmod
import pytest
import utils

# ====================================================================================================
def test_structures():

	types = ['C','H']
	myaev = aevmod.aev(types,8,4,4)

	symb, vxyz = utils.read_xyz("tests/xyz_c2h5.txt")
	cnf  = aevmod.config(symb)

	npt   = cnf.add_structures(vxyz);
	np.testing.assert_array_equal(vxyz,cnf.get_structures())
# ====================================================================================================
