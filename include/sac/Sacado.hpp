// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef SACADO_HPP
#define SACADO_HPP

// Kokkos::View specialization for Sacado AD classes
// #include "Kokkos_View_Fad.hpp"

// Version string
#include "Sacado_Version.hpp"

// Declarations of all overloaded math functions
#include "Sacado_MathFunctions.hpp"

// Traits for all of the Sacado classes -- Include these first so they are all
// defined before any nesting of AD classes
#ifdef SACADO_ENABLE_NEW_DESIGN
#include "Sacado_Fad_Exp_ExpressionTraits.hpp"
#include "Sacado_Fad_Exp_GeneralFadTraits.hpp"
#endif
#include "Sacado_Fad_ExpressionTraits.hpp"
#include "Sacado_Fad_DFadTraits.hpp"
#include "Sacado_Fad_SFadTraits.hpp"
#include "Sacado_Fad_SLFadTraits.hpp"

// Standard forward AD classes
#ifdef SACADO_ENABLE_NEW_DESIGN
#include "Sacado_Fad_Exp_DFad.hpp"
#include "Sacado_Fad_Exp_SFad.hpp"
#include "Sacado_Fad_Exp_SLFad.hpp"
#include "Sacado_Fad_Exp_ViewFad.hpp"
#endif
#include "Sacado_Fad_DFad.hpp"
#include "Sacado_Fad_SFad.hpp"
#include "Sacado_Fad_SLFad.hpp"

#endif // SACADO_KOKKOS_HPP
