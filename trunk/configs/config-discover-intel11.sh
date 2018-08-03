#!/bin/csh

#PETSc 3.2
#MPI /usr/local/intel/mpi/4.0.3.008/lib64/

./configure \
 --prefix=$ISSM_DIR \
 --with-matlab-dir="$ISSM_DIR/externalpackages/matlab/install" \
 --with-chaco-dir="$ISSM_DIR/externalpackages/chaco/install" \
 --with-triangle-dir=$ISSM_DIR/externalpackages/triangle/install \
 --with-metis-dir=$ISSM_DIR/externalpackages/metis/install \
 --with-petsc-dir=$ISSM_DIR/externalpackages/petsc/install \
 --with-mkl-dir="/usr/local/intel/mkl/10.1.2.024/lib/64/" \
 --with-mpi-include="/usr/local/intel/mpi/4.0.3.008/include64/" \
 --with-mpi-libflags="-L/usr/local/intel/mpi/4.0.3.008/lib64/ -lmpi -lmpiif" \
 --with-dakota-dir=$ISSM_DIR/externalpackages/dakota/install \
 --with-scalapack-dir="/usr/local/intel/mkl/10.1.2.024/lib/64/" \
 --with-blas-lapack-dir="/usr/local/intel/mkl/10.1.2.024/lib/64/" \
 --with-mumps-dir=$ISSM_DIR/externalpackages/petsc/install \
 --with-tao-dir=$ISSM_DIR/externalpackages/tao/install \
 --with-graphics-lib=/usr/lib64/libX11.so \
 --with-cxxoptflags="-O3 -xS -DMPICH_IGNORE_CXX_SEEK" \
 --with-vendor=intel-discover

