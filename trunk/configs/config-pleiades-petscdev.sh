#!/bin/csh

./configure \
 --prefix=$ISSM_DIR \
 --with-wrappers=no \
 --with-triangle-dir=$ISSM_DIR/externalpackages/triangle/install \
 --with-metis-dir=$ISSM_DIR/externalpackages/petsc/src/externalpackages/metis-5.0.2-p3 \
 --with-petsc-dir=$ISSM_DIR/externalpackages/petsc/install \
 --with-tao-dir=$ISSM_DIR/externalpackages/tao/install \
 --with-mpi-include=/nasa/sgi/mpt/2.04/include \
 --with-mpi-libflags="-L/nasa/sgi/mpt/2.04/lib/ -lmpi -lpthread -lgfortran" \
 --with-dakota-dir=$ISSM_DIR/externalpackages/dakota/install \
 --with-mkl-dir=/nasa/intel/mkl/10.0.011/lib/64/ \
 --with-mumps-dir=$ISSM_DIR/externalpackages/petsc/install/ \
 --with-scalapack-dir=$ISSM_DIR/externalpackages/petsc/install/ \
 --with-spai-dir=$ISSM_DIR/externalpackages/petsc/install/ \
 --with-hypre-dir=$ISSM_DIR/externalpackages/petsc/install/ \
 --with-superlu-dir=$ISSM_DIR/externalpackages/petsc/install/ \
 --with-spooles-dir=$ISSM_DIR/externalpackages/petsc/install/ \
 --with-pastix-dir=$ISSM_DIR/externalpackages/petsc/install/ \
 --with-graphics-lib=/usr/lib64/libX11.so \
 --with-cxxoptflags="-O3" \
 --with-vendor=intel-pleiades
