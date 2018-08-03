#!/bin/sh

#petsc 3.6
#mpich 3.0 yosemite
#boost 1.55 yosemite
#dakota 5.3.1 yosemite
export F77='/usr/local/gfortran/bin/gfortran'
export CC=/usr/local/gfortran/bin/gcc
export CXX=/usr/local/gfortran/bin/g++

./configure \
   --prefix=$ISSM_DIR \
	--with-matlab-dir=$MATLAB_DIR \
	--with-triangle-dir=$ISSM_DIR/externalpackages/triangle/install \
	--with-metis-dir=$ISSM_DIR/externalpackages/petsc/install \
	--with-petsc-dir=$ISSM_DIR/externalpackages/petsc/install  \
	--with-boost-dir=$ISSM_DIR/externalpackages/boost/install \
	--with-dakota-dir=$ISSM_DIR/externalpackages/dakota/install \
	--with-mpi-include=$ISSM_DIR/externalpackages/mpich/install/include  \
	--with-mpi-libflags=" $ISSM_DIR/externalpackages/mpich/install/lib/libpmpich.a $ISSM_DIR/externalpackages/mpich/install/lib/libmpich.a $ISSM_DIR/externalpackages/mpich/install/lib/libmpl.a " \
	--with-m1qn3-dir=$ISSM_DIR/externalpackages/m1qn3/install \
	--with-scalapack-dir=$ISSM_DIR/externalpackages/petsc/install \
	--with-mumps-dir=$ISSM_DIR/externalpackages/petsc/install \
	--with-chaco-dir="$ISSM_DIR/externalpackages/chaco/install" \
	--with-numthreads=8

