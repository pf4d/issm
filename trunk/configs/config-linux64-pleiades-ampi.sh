#!/bin/csh
./configure \
	--prefix=$ISSM_DIR \
	--without-kriging \
	--with-wrappers=no \
	--with-triangle-dir=$ISSM_DIR/externalpackages/triangle/install \
	--with-mpi-include="/nasa/sgi/mpt/2.06rp16/include" \
	--with-mpi-libflags="-L/nasa/sgi/mpt/2.06rp16/ -lmpi -lmpi++" \
	--with-adolc-dir="$ISSM_DIR/externalpackages/adolc/install"\
	--with-ampi-dir="$ISSM_DIR/externalpackages/adjoinablempi/install"\
	--with-mkl-dir="/nasa/intel/Compiler/2013.1.117/composer_xe_2013.1.117/mkl/lib/intel64" \
	--with-gsl-dir="$ISSM_DIR/externalpackages/gsl/install" \
	--with-metis-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-mumps-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-scalapack-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-hypre-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-graphics-lib="/usr/lib64/libX11.so" \
	--with-cxxoptflags="-O3" \
	--with-vendor="intel-pleiades"
