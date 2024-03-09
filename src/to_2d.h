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

#ifndef TO_2D_H
#define TO_2D_H
#include <stdexcept>

//===============================================================================================
// utility function to build a 2D vector-vector from a 1D vector of any type
/* inputs: 
	flat_vec: 1D vector
	ncols   	: number of columns desired in the resulting 2D vector
			  	the length of flat_vec has to be divisible by ncols, otherwise will throw an error
   outputs: 
	(returns) 	: an nrows x ncols 2D vector-vector 
 */

template < typename T >
std::vector< std::vector<T> > to_2d( const std::vector<T>& flat_vec, std::size_t ncols )
{
    // sanity check
    if( ncols == 0 || flat_vec.size()%ncols != 0 ) throw std::domain_error( "bad #cols" ) ;

    const auto nrows = flat_vec.size() / ncols ;

    std::vector< std::vector<T> > mtx ;
    const auto begin = std::begin(flat_vec) ;

    for( std::size_t row = 0 ; row < nrows ; ++row ) 
    	mtx.push_back( { begin + row*ncols, begin + (row+1)*ncols } ) ;

    return mtx ;
}
//===============================================================================================

#endif
