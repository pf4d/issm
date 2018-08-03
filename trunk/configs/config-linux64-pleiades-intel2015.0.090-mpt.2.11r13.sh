#!/bin/csh
./configure \
	--prefix=$ISSM_DIR \
	--without-kriging \
	--with-wrappers=no \
	--with-petsc-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-triangle-dir="$ISSM_DIR/externalpackages/triangle/install" \
	--with-mpi-include="/nasa/sgi/mpt/2.11r13/include" \
	--with-mpi-libflags="-L/nasa/sgi/mpt/2.11r13/ -lmpi" \
	--with-mkl-dir="/nasa/intel/Compiler/2015.0.090/composer_xe_2015.0.090/mkl/" \
	--with-metis-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-mumps-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-scalapack-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-cxxoptflags="-O3 -axAVX" \
	--with-vendor="intel-pleiades"
