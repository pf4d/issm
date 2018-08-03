#!/bin/bash
set -eu

#Some cleanup
rm -rf install petsc-3.6.3 src
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/petsc-lite-3.6.3.tar.gz' 'petsc-3.6.3.tar.gz'

#Untar and move petsc to install directory
tar -zxvf  petsc-3.6.3.tar.gz
mv petsc-3.6.3/* src/
rm -rf petsc-3.6.3

#patch
cat src/config/BuildSystem/config/packages/metis.py | sed -e "s/metis-5.1.0-p1/metis-5.1.0-p3/g" > TEMP
mv TEMP src/config/BuildSystem/config/packages/metis.py

#configure
cd src
./config/configure.py \
	--prefix="$ISSM_DIR/externalpackages/petsc/install" \
	--with-mpi-dir="$ISSM_DIR/externalpackages/mpich/install" \
	--PETSC_DIR="$ISSM_DIR/externalpackages/petsc/src" \
	--with-debugging=0 \
	--with-valgrind=0 \
	--with-x=0 \
	--with-ssl=0 \
	--with-shared-libraries=1 \
	--download-metis=1 \
	--download-parmetis=1 \
	--download-mumps=1 \
	--download-scalapack=1

#Compile and intall
make
make install
