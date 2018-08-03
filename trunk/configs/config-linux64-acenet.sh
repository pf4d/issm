#!/bin/bash

./configure \
 --prefix=/home/klemorza/issm/trunk-jpl \
 --with-metis-dir=$ISSM_DIR/externalpackages/metis/install \
 --with-metis-dir=$ISSM_DIR/externalpackages/petsc/install
 --with-petsc-dir=$ISSM_DIR/externalpackages/petsc/install \
 --with-mpi-include=$OPENMPI/include  \
 --with-mpi-libflags="-L$OPENMPI/lib -lmpi" \
 --with-blas-lapack-dir=$ISSM_DIR/externalpackages/petsc/install \
 --with-scalapack-dir=$ISSM_DIR/externalpackages/petsc/install/ \
 --with-mumps-dir=$ISSM_DIR/externalpackages/petsc/install/ \
 --with-cxxoptflags="-O2" \
 --with-numthreads=32 \
 --with-serial=no \
 --with-modules=no \
 --enable-debugging CC=mpicc CXX=mpiCC F77=mpif77
