#!/bin/csh
./configure \
	--prefix=$ISSM_DIR \
	--without-kriging \
	--with-wrappers=no \
	--with-triangle-dir="$ISSM_DIR/externalpackages/triangle/install" \
	--with-mpi-include=$ISSM_DIR/externalpackages/mpich/install/include  \
	--with-mpi-libflags="-L$ISSM_DIR/externalpackages/mpich/install/lib/ -lmpich -lmpl " \
	--with-adolc-dir="$ISSM_DIR/externalpackages/adolc/install"\
	--with-ampi-dir="$ISSM_DIR/externalpackages/adjoinablempi/install"\
	--with-metis-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-blas-lapack-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-gsl-dir="$ISSM_DIR/externalpackages/gsl/install" \
	--with-mumps-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-scalapack-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-hypre-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-graphics-lib="/usr/lib64/libX11.so" \
	--with-cxxoptflags=" -O3 -march=corei7-avx" \
	--with-vendor="intel-pleiades-gcc"  \
	CFLAGS="-O3 -march=corei7-avx" \
	CXXFLAGS="-O3 -march=corei7-avx" \
	FFLAGS="-O3 -march=corei7-avx" 
