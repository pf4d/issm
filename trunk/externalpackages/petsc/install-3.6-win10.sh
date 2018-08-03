#!/bin/bash
set -eu

#Some cleanup
rm -rf install petsc-3.6.2 src
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/petsc-lite-3.6.2.tar.gz' 'petsc-3.6.2.tar.gz'

#Untar and move petsc to install directory
tar -zxvf  petsc-3.6.2.tar.gz
mv petsc-3.6.2/* src/
rm -rf petsc-3.6.2

export PETSC_DIR=`cygpath -u "$ISSM_DIR/externalpackages/petsc/src"`
export PREFIX_DIR=`cygpath -u "$ISSM_DIR/externalpackages/petsc/install"`

#configure
cd src
./config/configure.py  \
	--with-parallel-no \
	--prefix=$PREFIX_DIR \
	--PETSC_ARCH=cygwin-intel \
	--PETSC_DIR=$PETSC_DIR \
	--with-mpi=0 \
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
