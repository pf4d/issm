#!/bin/csh
./configure \
	--prefix=$ISSM_DIR \
	--without-kriging \
	--with-matlab-dir="$ISSM_DIR/externalpackages/matlab/install" \
	--with-petsc-dir=$ISSM_DIR/externalpackages/petsc/install \
	--with-triangle-dir="$ISSM_DIR/externalpackages/triangle/install" \
	--with-mpi-include="/nasa/mvapich2/1.6.sles11/gcc/include" \
	--with-mpi-libflags="-L/nasa/mvapich2/1.6.sles11/gcc/lib -lmpich -lopa -lpthread -lrdmacm -libverbs -libumad -ldl -lrt" \
	--with-blas-lapack-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-metis-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-mumps-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-scalapack-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-hypre-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-graphics-lib="/usr/lib64/libX11.so" \
	--with-cxxoptflags="-g -O2" \
	--with-vendor="intel-pleiades-gcc"
