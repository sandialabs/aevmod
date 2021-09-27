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

#ifndef ISEQUAL_H
#define ISEQUAL_H
//===============================================================================================
// utility function to check if two vectors of any type are equivalent
// it checks that the two vectors are of the same size, and that their entries are equal
/* inputs:
	vectors: 		v1, v2
   outputs: 
   	(returns):		bool ... true if equivalent else false
*/

template<typename T>
bool isEqual(const std::vector<T>& v1, const std::vector<T>& v2)
{
    return (v1.size() == v2.size() &&
            std::equal(v1.begin(), v1.end(), v2.begin()));
}
//===============================================================================================
#endif