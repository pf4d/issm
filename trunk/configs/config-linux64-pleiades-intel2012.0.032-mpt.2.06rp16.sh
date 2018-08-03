#!/bin/csh
./configure \
	--prefix=$ISSM_DIR \
	--without-kriging \
	--with-matlab-dir="$ISSM_DIR/externalpackages/matlab/install" \
	--with-petsc-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-triangle-dir="$ISSM_DIR/externalpackages/triangle/install" \
	--with-mpi-include="/nasa/sgi/mpt/2.06rp16/include" \
	--with-mpi-libflags="-L/nasa/sgi/mpt/2.06rp16/ -lmpi" \
	--with-mkl-dir="/nasa/intel/Compiler/2012.0.032/composer_xe_2011_sp1/mkl/lib/intel64"\
	--with-metis-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-mumps-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-scalapack-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-hypre-dir="$ISSM_DIR/externalpackages/petsc/install" \
	--with-graphics-lib="/usr/lib64/libX11.so" \
	--with-cxxoptflags="-g -O2" \
	--with-vendor="intel-pleiades"
