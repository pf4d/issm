#!/bin/sh

#need to swap to GNU environment
#need to change the petsc libflag to lcraypetsc...
./configure \
		--prefix=/work/bfl022/issm_install \
		--with-m1qn3-dir=$ISSM_DIR/externalpackages/m1qn3/install \
		--with-petsc-dir=$PETSC_DIR \
		--with-metis-dir=$CRAY_TPSL_PREFIX_DIR \
		--with-mumps-dir=$CRAY_TPSL_PREFIX_DIR \
		--with-mpi-include="$CRAY_MPICH2_DIR/include" \
		--with-mpi-libflags="-L$CRAY_MPICH2_DIR/lib -lmpich -lmpl -lfmpich -lmpichcxx -lmpichf90" \
		--with-gsl-dir="/work/apps/gsl/1.16-gnu/" \
		--enable-development
