# Author: Alfio Lazzaro, alazzaro@cray.com (2013)

#
# target compiler... these are the options used for building the library.
# They should be aggessive enough to e.g. perform vectorization for the specific CPU (e.g. -ftree-vectorize -march=native),
# and allow some flexibility in reordering floating point expressions (-ffast-math).
# Higher level optimisation (in particular loop nest optimization) should not be used.
#
target_compile="gfortran -O2 -funroll-loops -ffast-math -ftree-vectorize -march=native -cpp -finline-functions -fopenmp"

#
# target dgemm link options... these are the options needed to link blas (e.g. -lblas)
# blas is used as a fall back option for sizes not included in the library or in those cases where it is faster
# the same blas library should thus also be used when libsmm is linked.
#
blas_linking="-L${INTEL_PATH}/mkl/lib/intel64 -Wl,--start-group -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group"

#
# host compiler... this is used only to compile a few tools needed to build
# the library. The library itself is not compiled this way.
# This compiler needs to be able to deal with some Fortran2003 constructs.
#
host_compile="gfortran -O2"
