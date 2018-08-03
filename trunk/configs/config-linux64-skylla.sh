#!/bin/sh

#External packages versions:
#petsc 3.4
#mpich 3.0

./configure \
 --prefix=$ISSM_DIR \
 --with-matlab-dir=$ISSM_DIR/externalpackages/matlab/install \
 --with-triangle-dir=$ISSM_DIR/externalpackages/triangle/install \
 --with-metis-dir=$ISSM_DIR/externalpackages/petsc/install \
 --with-petsc-dir=$ISSM_DIR/externalpackages/petsc/install \
 --with-mpi-include=$ISSM_DIR/externalpackages/mpich/install/include  \
 --with-mpi-libflags="-L$ISSM_DIR/externalpackages/mpich/install/lib/ -lmpich -lmpl " \
 --with-blas-lapack-dir=$ISSM_DIR/externalpackages/petsc/install \
 --with-tao-dir=$ISSM_DIR/externalpackages/tao/install/ \
 --with-scalapack-dir=$ISSM_DIR/externalpackages/petsc/install/ \
 --with-mumps-dir=$ISSM_DIR/externalpackages/petsc/install/ \
 --with-fortran-lib="-L/usr/lib/gcc/x86_64-redhat-linux/4.1.2/ -lgfortran" \
 --with-graphics-lib="/usr/lib64/libX11.so.6" \
 --with-cxxoptflags="-g -O2" \
 --with-numthreads=8
