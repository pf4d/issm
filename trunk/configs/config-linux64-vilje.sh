#!/bin/sh
#Vilje configuration script
#this uses system stuff, you need to load the following:
#module load intelcomp/15.0.1 mpt/2.10 petsc/3.5.3d
#module load parmetis/4.0.3 mumps/5.0.0
#module load automake

./configure \
    --prefix=$ISSM_DIR \
    --without-kriging \
    --with-wrappers=no \
    --with-petsc-dir="/sw/sdev/Modules/petsc/petsc-3.5.3" \
    --with-triangle-dir="$ISSM_DIR/externalpackages/triangle/install" \
    --with-mpi-include="/sw/sgi/mpt/mpt-2.10/include/"  \
    --with-mpi-libflags="-L/sw/sgi/mpt/mpt-2.10/lib -lmpi -lmpi++" \
    --with-metis-dir="/sw/sdev/Modules/parmetis/parmetis-4.0.3/include" \
    --with-mumps-dir="/sw/sdev/Modules/mumps/mumps-5.0.0/lib" \
    --with-cxxoptflags="-O2 -xAVX" \
    --with-m1qn3-dir="$ISSM_DIR/externalpackages/m1qn3/install" \
    --enable-development \
    CC=/sw/sdev/Modules/intelcomp/psxe_2015/bin/icc \
    CXX=/sw/sdev/Modules/intelcomp/psxe_2015/bin/icpc \
    F77=/sw/sdev/Modules/intelcomp/psxe_2015/bin/ifort
