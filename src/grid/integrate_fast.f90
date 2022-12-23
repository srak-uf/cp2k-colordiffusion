! ideally users test for the optional __GRID_CORE macro definition
! otherwise, we provide some reasonable defaults
! the included fortran files represent the same kernel, with a fortran representation that optimizes best for a given compiler

#ifdef __HAS_LIBGRID
! Nothing here, the libgrid.a is present
#else

#if !defined(__GRID_CORE)
#define __GRID_CORE 4
#endif

#if __GRID_CORE == 1

#include "integrate_fast_1.f90"

#elif __GRID_CORE == 2

#include "integrate_fast_2.f90"

#elif __GRID_CORE == 3

#include "integrate_fast_3.f90"

#elif __GRID_CORE == 4

#include "integrate_fast_4.f90"

#elif __GRID_CORE == 5

#include "integrate_fast_5.f90"

#elif __GRID_CORE == 6

#include "integrate_fast_6.f90"

#else

This is an error, and unknown definition of GRID_CORE(__GRID_CORE) has been used

#endif

#endif
