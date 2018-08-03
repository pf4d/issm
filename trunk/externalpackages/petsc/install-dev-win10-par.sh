#!/bin/bash

STEP=1

rm -rf src

if [ $STEP -eq 1 ]; then
	git clone -b maint https://bitbucket.org/petsc/petsc src
fi

export PETSC_DIR=`cygpath -u "$ISSM_DIR/externalpackages/petsc/src"`
export PREFIX_DIR=`cygpath -u "$ISSM_DIR/externalpackages/petsc/install"`

#configure
cd src
./config/configure.py  \
	--prefix=$PREFIX_DIR \
	--PETSC_ARCH=cygwin-intel \
	--PETSC_DIR=$PETSC_DIR \
	--with-mpi-dir="/cygdrive/c/Programs/MPICH2/"\
	--with-debugging=0 \
	--with-valgrind=0 \
	--with-x=0 \
	--with-ssl=0 \
	--download-f2cblaslapack=yes \
	--with-cc='win32fe cl' \
	--with-fc=0 \
	--with-cxx='win32fe cl' \
	--with-clanguage=cxx 

#Compile petsc and install it
make
make install

patch ../install/include/petscfix.h ../configs/3.1/win7/petscfix.h.patch
