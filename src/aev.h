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

#ifndef AEV_H
#define AEV_H

#include "config.h"

//===============================================================
class aev {

//--------------------------------------------------------------------------------
public:

  aev(const std::vector<std::string>& atom_types={"C","H"}, const int& nrho_rad=32, 
      const int& nrho_ang=8, const int& nalpha=8, const std::vector<double>& R_c={4.6,3.1});

  void build_index_sets(config& conf);

  void bld_index_sets(
    const std::vector<std::string>& symb, 
    std::vector<std::vector<int>>& rad_ind_set, 
    std::vector<std::vector<std::vector<int>>>& ang_ind_set);

  std::vector<std::vector<std::vector<double>>> eval(config& conf);
  std::vector<std::vector<std::vector<double>>> eval_sac(config& conf);

  std::vector<std::vector<std::vector<std::vector<double>>>> eval_Jac (config& conf);
  std::vector<std::vector<std::vector<std::vector<double>>>> eval_Jac_sac (config& conf);
  std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> eval_Hess_sac (config& conf);

  // cutoff function
  double fc(const double& Rij, const double& Rc);

  // kernel for Radial SFs
  double rad_kern(const std::vector<double>& par, const double& Rij, const double& Rc);

  // kernel for Angular SFs
  double ang_kern(const std::vector<double>& par, const double& Rij, const double& Rik, const double& theta, const double& Rc);

  template<typename ValueType>
  std::vector<ValueType> evaluateFunction(const std::vector<ValueType> &xin, const int p, const int i, config& conf);
  template<typename ValueType>
  std::vector<std::vector<ValueType>> evaluateJacobian(const std::vector<ValueType> &x, const int p, const int ia, config& conf);
  template<typename ValueType>
  void evaluateHessian(const std::vector<ValueType> &x,
         std::vector<std::vector<ValueType>> &g,
         std::vector<std::vector<std::vector<ValueType>>> &h, const int p, const int ia, config& conf);

  template<typename ValueType>
  ValueType fc_sac(const ValueType& Rij, const double& Rc);
  template<typename ValueType>
  ValueType rad_kern_sac(const std::vector<double>& par, const ValueType& Rij, const double& Rc);
  template<typename ValueType>
  ValueType ang_kern_sac(const std::vector<double>& par, const ValueType& Rij, const ValueType& Rik, const ValueType& theta, const double& Rc);

  void write_aev_to_file(const std::string& fname);
  void write_jac_to_file(const std::string& fname);
  void get_aev_tags(int& rad_dout, int& ang_dout, std::vector<std::string>& aev_tag);

  int n_rad;
  int n_ang;
  int dout;

  std::vector<std::vector<std::vector<double>>> val_aev;
  std::vector<std::vector<std::vector<std::vector<double>>>> val_jac;
  std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> val_Hess;

//--------------------------------------------------------------------------------
private:
  const double pi = 4.0*atan(1.0);
  std::vector<std::vector<double>> rad_par;
  std::vector<std::vector<double>> ang_par;
  std::vector<std::string> types;
  std::vector<double> R_c;

  std::vector<std::string> rad_typ;
  std::vector<std::vector<std::string>> ang_typ;
  std::vector<std::string> tag;

//--------------------------------------------------------------------------------
};

//================================================================================
#endif
