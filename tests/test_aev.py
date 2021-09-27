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
def test_aev():

	# define atom types
	types = ['C','H']

	# instantiate aev object
	myaev = aevmod.aev(types,8,4,4)

	# read geometry xyz file
	symb, vxyz = utils.read_xyz("tests/xyz_c2h5.txt")

	# instantiate cnf object
	cnf  = aevmod.config(symb)

	# build index sets
	myaev.build_index_sets(cnf)

	# add the structure to the cnf object data
	npt     = cnf.add_structures(vxyz);

	# evaluate aev of each 
	got_aev = myaev.eval(cnf)

	# read true aev data for this system
	nptt, n_atom, dout, tru_aev = utils.read_aev("tests/aev_c2h5_8_4_4.txt");

	# check
	try:
		# confirm that the two files pertain to the same npt value
		assert npt == nptt

		# confirm that the aev we evaluate is close enough to the reference true value
		assert np.allclose(got_aev,tru_aev, rtol=1.e-15, atol=1.e-15)
	except AssertionError:
	    print("got_aev disagrees with tru_aev")
	    errmx = 0.0
	    for p in range(npt):
	        for a in range(n_atom):
	            for r in range(dout):
	                errmx=max(errmx,abs(got_aev[p][a][r]-tru_aev[p][a][r]))
	    print("errmx:",errmx)
	    sys.exit(1)
# ====================================================================================================
