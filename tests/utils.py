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

# ====================================================================================================
def read_hes(fname):
    """ Reads reference Hessian values from file
    """

    try:
        f = open(fname, "r")
    except IOError:
        print("Could not open file:" + fname)
        sys.exit()
    with f:
        dat = f.readlines()

    n_line = len(dat)
    npt    = int(dat[0])
    n_atom = int(dat[1])
    dout   = int(dat[2])
    din    = int(dat[3])

    hes = [ [ [ [ [0]*din for k in range(din)] for r in range(dout)] for a in range(n_atom)] for p in range(npt)]
    line = 4

    for p in range(npt):
    	for a in range(n_atom):
    		for r in range(dout):
	    		for k in range(din):
	    			for l in range(din):
    					hes[p][a][r][k][l]=float(dat[line])
    					line += 1
    return npt, n_atom, dout, din, hes

# ====================================================================================================
def read_jac(fname):
    """ Reads reference Jacobian values from file
    """

    try:
        f = open(fname, "r")
    except IOError:
        print("Could not open file:" + fname)
        sys.exit()
    with f:
        dat = f.readlines()

    n_line = len(dat)
    npt    = int(dat[0])
    n_atom = int(dat[1])
    dout   = int(dat[2])
    din    = int(dat[3])

    jac = [ [ [ [0]*din for r in range(dout)] for a in range(n_atom)] for p in range(npt)]
    line = 4

    for p in range(npt):
    	for a in range(n_atom):
    		for r in range(dout):
	    		for k in range(din):
    				jac[p][a][r][k]=float(dat[line])
    				line += 1
    return npt, n_atom, dout, din, jac

# ====================================================================================================
def read_aev(fname):
    """ Reads aev values from file
    """

    try:
        f = open(fname, "r")
    except IOError:
        print("Could not open file:" + fname)
        sys.exit()
    with f:
        aevd = f.readlines()

    n_line = len(aevd)
    npt    = int(aevd[0])
    n_atom = int(aevd[1])
    dout   = int(aevd[2])

    aev = [ [ [0]*dout for a in range(n_atom)] for p in range(npt)]
    line = 3
    for p in range(npt):
    	for a in range(n_atom):
    		for i in range(dout):
    			aev[p][a][i]=float(aevd[line])
    			line += 1
    return npt, n_atom, dout, aev

# ====================================================================================================
def read_rad_indsets(fname):
    """ Reads radial index sets from file
    """

    try:
        f = open(fname, "r")
    except IOError:
        print("Could not open file:" + fname)
        sys.exit()
    with f:
        risd = f.readlines()

    n_line = len(risd)

    ris = [ [] for l in range(n_line) ]
    for line in range(n_line):
        ris[line] = (" ".join(risd[line].split())).split(" ")
        ris[line] = [ int(r) for r in ris[line] ]

    return ris

# ====================================================================================================
# ====================================================================================================
def read_ang_indsets(fname):
    """ Reads angular index sets from file
    """

    try:
        f = open(fname, "r")
    except IOError:
        print("Could not open file:" + fname)
        sys.exit()
    with f:
        aisd = f.readlines()

    n_line = len(aisd)
    ais = [ [ [] ] for l in range(n_line) ]
    for line in range(n_line):
        ais[line] = (" ".join(aisd[line].split())).split(" ")
        ais[line] = [ [int(ais[line][r]),int(ais[line][r+1])] for r in range(0,len(ais[line]),2) ]
    return ais

# ====================================================================================================

# ====================================================================================================
def read_indsets(rfname,afname):
    """ Reads index sets from two files
    """

    try:
        f = open(rfname, "r")
    except IOError:
        print("Could not open file:" + rfname)
        sys.exit()
    with f:
        risd = f.readlines()

    try:
        f = open(afname, "r")
    except IOError:
        print("Could not open file:" + afname)
        sys.exit()
    with f:
        aisd = f.readlines()

    n_line = len(risd)

    ris = [ [] for l in range(n_line) ]
    for line in range(n_line):
        ris[line] = (" ".join(risd[line].split())).split(" ")
        ris[line] = [ int(r) for r in ris[line] ]

    n_line = len(aisd)
    ais = [ [ [] ] for l in range(n_line) ]
    for line in range(n_line):
        ais[line] = (" ".join(aisd[line].split())).split(" ")
        ais[line] = [ [int(ais[line][r]),int(ais[line][r+1])] for r in range(0,len(ais[line]),2) ]
    return ris, ais

# ====================================================================================================
# ====================================================================================================
def parse_xyz(name):
    """ Reads xyz file and returns a list containing two items:
		1) a list of strings, being chemical symbols of each atom in the system
		2) a 2d numpy array with n_atom rows and 3 columns where
		   each row is the (x,y,z) coordinates of an atom in the system
        NOTE: this assumes the file contains a single molecule/configuration geometry spec
	"""
    try:
        f = open(name, "r")
    except IOError:
        print("Could not open file:" + name)
        sys.exit()
    with f:
        xyz = f.readlines()

    n_line = len(xyz)
    n_atom = int(xyz[0])
    assert (n_line == n_atom + 2)

    symb = [" "] * n_atom
    x = np.zeros((n_atom, 3))

    for line in range(2, n_line):
        lst = (" ".join(xyz[line].split())).split(" ")
        symb[line - 2] = lst[0]
        x[line - 2] = np.array(lst[1:4])

    return symb, x

# ====================================================================================================
# ====================================================================================================
def read_xyz(name):
    """
	 Reads an xyz file using parse_xyz() and repackages the 2d array of
	 atom positions returned from it, converting it from a n_atom x 3 2d array
	 to a 1 x 3n_atom 2d array
	"""
    symb, x = parse_xyz(name)
    x = np.array([x.flatten()])

    return symb, x
# ====================================================================================================
