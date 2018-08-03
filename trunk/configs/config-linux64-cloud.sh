#!/bin/sh

#External packages versions:
#petsc 3.1 or 3.2
#mpich 1.4

./configure \
 --prefix=$ISSM_DIR \
 --without-wrappers \
 --with-triangle-dir=$ISSM_DIR/externalpackages/triangle/install \
 --with-metis-dir=$ISSM_DIR/externalpackages/petsc/src/externalpackages/metis-5.0.2-p3 \
 --with-petsc-dir=$ISSM_DIR/externalpackages/petsc/install \
 --with-tao-dir=$ISSM_DIR/externalpackages/tao/install \
 --with-mpi-include=$ISSM_DIR/externalpackages/mpich/install/include  \
 --with-mpi-libflags="-L$ISSM_DIR/externalpackages/mpich/install/lib/ -lmpich -lmpl " \
 --with-scalapack-dir=$ISSM_DIR/externalpackages/petsc/install/ \
 --with-mumps-dir=$ISSM_DIR/externalpackages/petsc/install/ \
 --with-cxxoptflags=" -O2 -fpermissive" 
