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

# A simple test code for computing the AEV and its Jacobian & Hessian
import numpy as np
import aevmod

# define types of atom in system
types = ['C','H']

# define AEV structure object
nrho_rad = 32  # number of radial aev radial shells
nrho_ang = 8   # number of angular aev radial shells
nalpha   = 8   # number of angular aev angular sectors dividing [0,pi]
R_c_rad  = 4.6 # radial aev cutoff radius (Angstroms)
R_c_ang  = 3.1 # angular aev cutoff radius (Angstroms)
myaev    = aevmod.aev(types, nrho_rad, nrho_ang, nalpha, [R_c_rad,R_c_ang])
print("built aev object, AEV size:",myaev.dout)

# define CH2 molecule symbol list
symb  = ['C','H','H']
print("configuration:",symb)

# define CH2 molecule object
cnf   = aevmod.config(symb)

# build index sets
myaev.build_index_sets(cnf)
print("radial  index sets:",cnf.get_radial_index_set())
print("angular index sets:",cnf.get_angular_index_set())

# define numpy array composed of two CH2 xyz structures
vxyz  = np.array([[0., 0., 0., 0., 0., 1.10771, 1.08378, 0., -0.22899],
                  [0., 0., 0., 0., 0., 1.07672, 0.77665, 0., -0.74575]
                 ])
print("structures:")
with np.printoptions(precision=4, suppress=True):
    for v in vxyz:
        print(v)

# add vxyz structures to cnf
npt   = cnf.add_structures(vxyz)
print("number of structures:",npt)

# evaluate AEVs for the array of structures
# got_aev[j] is the AEV for structure j
# got_aev[j][k] is the AEV for atom k in structure j
got_aev = myaev.eval(cnf)

# write out AEVs
fname = "aev.out"
print("printing AEVs to",fname)
np.savetxt(fname,["# AEV:"],fmt='%s')
with open(fname, "a") as f:
    for j in range(npt):
        for k in range(len(symb)):
            tag="# structure:"+str(j)+", atom:"+symb[k]+"\n"
            f.write(tag)
            np.savetxt(f,got_aev[j][k])

# evaluate Jacobians of AEVs
got_jac = myaev.eval_Jac(cnf)

# write out Jacobians
fname = "jac.out"
print("printing Jacobians to",fname)
np.savetxt(fname,["# Jac:"],fmt='%s')
with open(fname, "a") as f:
    for j in range(npt):
        for k in range(len(symb)):
            tag="# structure:"+str(j)+", atom:"+symb[k]+"\n"
            f.write(tag)
            np.savetxt(f,got_jac[j][k])

# evaluate Hessians of AEVs
got_hes = myaev.eval_Hess_sac(cnf)

# write out Hessians
fname = "hes.out"
print("printing Hessians to",fname)
np.savetxt(fname,["# Hes:"],fmt='%s')
with open(fname, "a") as f:
    for j in range(npt):
        for k in range(len(symb)):
            for l in range(myaev.dout):
                tag="# structure:"+str(j)+", atom:"+symb[k]+", AEV element:"+str(l)+"\n"
                f.write(tag)
                np.savetxt(f,got_hes[j][k][l])

