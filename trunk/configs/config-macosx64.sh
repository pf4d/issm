#!/bin/sh

#petsc 3.5
#mpich 3.0.4

./configure \
	--prefix=$ISSM_DIR \
	--with-matlab-dir="$ISSM_DIR/externalpackages/matlab/install" \
	--with-triangle-dir=$ISSM_DIR/externalpackages/triangle/install \
	--with-metis-dir=$ISSM_DIR/externalpackages/petsc/install \
	--with-petsc-dir=$ISSM_DIR/externalpackages/petsc/install  \
	--with-mpi-include=$ISSM_DIR/externalpackages/mpich/install/include  \
	--with-mpi-libflags=" $ISSM_DIR/externalpackages/mpich/install/lib/libpmpich.a $ISSM_DIR/externalpackages/mpich/install/lib/libmpich.a $ISSM_DIR/externalpackages/mpich/install/lib/libmpl.a " \
	--with-scalapack-dir=$ISSM_DIR/externalpackages/petsc/install \
	--with-mumps-dir=$ISSM_DIR/externalpackages/petsc/install \
	--with-fortran-lib="-L/usr/local/gfortran/lib/ -lgcc_s.1"\
	--with-graphics-lib="/usr/X11/lib/libX11.dylib" \
	--with-numthreads=4
