#!/bin/csh
set echo
module load intel impi
module load netcdf netcdf-fortran
mkdir obj
cp intel.mk Makefile.fms_spectral_solo Makefile input.nml obj
cd obj
make
cd ../../src/mppnccombine
./mppnccombine_compile