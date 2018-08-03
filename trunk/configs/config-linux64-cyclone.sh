#!/bin/sh

#Default configuration script
#automake 1.14
#petsc 3.6
#mpich 3.0

./configure \
		--prefix=$ISSM_DIR \
		--with-matlab-dir="$ISSM_DIR/externalpackages/matlab/install" \
		--with-triangle-dir="$ISSM_DIR/externalpackages/triangle/install" \
		--with-m1qn3-dir="$ISSM_DIR/externalpackages/m1qn3/install" \
		--with-mpi-include="$ISSM_DIR/externalpackages/mpich/install/include"  \
		--with-mpi-libflags="-L$ISSM_DIR/externalpackages/mpich/install/lib/ -lmpich -lmpl" \
		--with-petsc-dir="$ISSM_DIR/externalpackages/petsc/install" \
		--with-metis-dir="$ISSM_DIR/externalpackages/petsc/install" \
		--with-scalapack-dir="$ISSM_DIR/externalpackages/petsc/install/" \
		--with-mumps-dir="$ISSM_DIR/externalpackages/petsc/install/" \
		--with-fortran-lib="-L/opt/intel/intelcompiler-12.04/composerxe-2011.4.191/compiler/lib/intel64 -lifport -lifcore -limf -lsvml -lintlc "\
		--enable-development \
		CC=mpicc \
		CXX=mpic++