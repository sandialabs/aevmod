/* Define if want to build teuchos-complex */
/* #undef HAVE_SACADO_COMPLEX */

/* Define if want to build with uninit */
#define HAVE_SACADO_UNINIT

/* Define if ADOL-C is enabled */
/* #undef HAVE_ADOLC */

/* Define if ADIC is enabled */
/* #undef HAVE_ADIC */

/* define if we want to use MPI */
/* #undef HAVE_MPI */

/* Define if want to build sacado-examples */
/* #undef HAVE_SACADO_EXAMPLES */

/* Define if want to build sacado-tests */
/* #undef HAVE_SACADO_TESTS */

/* Define if want to build with teuchos-core enabled */
/* #define HAVE_SACADO_TEUCHOSCORE */

/* Define if want to build with teuchos-comm enabled */
/* #define HAVE_SACADO_TEUCHOSCOMM */

/* Define if want to build with teuchos-numerics enabled */
/* #define HAVE_SACADO_TEUCHOSNUMERICS */

/* Define if want to build with kokkos-core enabled */
/* #define HAVE_SACADO_KOKKOSCORE */
#if !defined(HAVE_SACADO_KOKKOSCORE)
#define KOKKOS_INLINE_FUNCTION inline
#endif

/* Define if want to build with kokkos-containers enabled */
/* #define HAVE_SACADO_KOKKOSCONTAINERS */

/* Define if want to build with TeuchosKokkosComm subpackage enabled */
/* #define HAVE_SACADO_TEUCHOSKOKKOSCOMM */
#ifdef HAVE_SACADO_TEUCHOSKOKKOSCOMM
// For backwards compatibility
#  define HAVE_SACADO_KOKKOSMPICOMM
#endif // HAVE_SACADO_TEUCHOSKOKKOSCOMM

/* Define if want to enable Kokkos view specializations for Sacado */
#define HAVE_SACADO_VIEW_SPEC

/* define if the compiler is confused by std::sin, etc., within namespace
   Sacado::Rad */
/* #undef RAD_NO_USING_STDCC */

/* Define to enable extra debugging checks */
/* #undef SACADO_DEBUG */

/* Define if compiler supports c99 tr1 cmath functions */
#define HAS_C99_TR1_CMATH

/* Define to enable C++11 support*/
#define HAVE_SACADO_CXX11

/* Define to enable hierarchical parallelism support */
#if !defined(SACADO_VIEW_CUDA_HIERARCHICAL) && !defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
/* #undef SACADO_VIEW_CUDA_HIERARCHICAL */
#endif

/* Define to enable hierarchical parallelism support specific to DFad */
#if !defined(SACADO_VIEW_CUDA_HIERARCHICAL) && !defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
/* #undef SACADO_VIEW_CUDA_HIERARCHICAL_DFAD */
#endif

/* Define to enable Kokkos memory pool in Sacado */
#ifndef SACADO_KOKKOS_USE_MEMORY_POOL
/* #undef SACADO_KOKKOS_USE_MEMORY_POOL */
#endif

/* Define if want to make the new Fad design the default, replacing the old one */
#define SACADO_NEW_FAD_DESIGN_IS_DEFAULT

/* Force SFad (in the new design) to initialize value and derivative components in the default constructor (adds additional runtime cost, but protects against uninitialized use). */
#define SACADO_SFAD_INIT_DEFAULT_CONSTRUCTOR
