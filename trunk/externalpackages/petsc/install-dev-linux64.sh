#!/bin/bash

STEP=1

if [ $STEP -eq 1 ]; then
	# Used Mercurial to get code
	rm -rf src
	hg clone http://petsc.cs.iit.edu/petsc/petsc-dev src
	cd src
	hg clone http://petsc.cs.iit.edu/petsc/BuildSystem config/BuildSystem
fi

# To update (via Mercurial):
#      cd petsc-dev
#      hg pull -u
#      cd config/BuildSystem
#      hg pull -u

# configure script
# Note: 
#  Metis: -using metis from externalpackages did not work...
#         -for now download new metis
#         -rename metis in externalpackages to metis2
#
# SuperLU: -If download-..-=yes does not work try downloading from
#    --download-superlu=http://crd.lbl.gov/~xiaoye/SuperLU/superlu_4.3.tar.gz \


if [ $STEP -eq 2 ]; then
	rm -rf install
	cd src
	./config/configure.py \
	--prefix="$ISSM_DIR/externalpackages/petsc/install" \
	--with-mpi-dir="$ISSM_DIR/externalpackages/mpich/install" \
	--with-clanguage=C++ \
	--PETSC_ARCH=linux-gnu-amd64 \
	--PETSC_DIR="$ISSM_DIR/externalpackages/petsc/src" \
	--with-debugging=0 \
	--with-shared-libraries=1 \
	--download-mumps=yes \
	--download-scalapack=yes \
	--download-blacs=yes  \
	--download-blas=yes \
	--download-f-blas-lapack=yes \
	--download-parmetis=yes \
	--download-metis=yes \
	--download-trilinos=yes \
	--download-euclid=yes \
	--download-spooles=yes \
	--download-spai=yes \
	--download-superlu=http://crd.lbl.gov/~xiaoye/SuperLU/superlu_4.3.tar.gz \
	--download-hypre=yes \
	--download-c2html=yes
#	--with-pic=1

	#Compile petsc and install it
	make
	make install
fi
