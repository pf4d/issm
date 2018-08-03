#!/bin/bash
set -eu

#Some cleanup
rm -rf install petsc-3.4.3 src
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/petsc-lite-3.4.3.tar.gz' 'petsc-3.4.3.tar.gz'

#Untar and move petsc to install directory
tar -zxvf  petsc-3.4.3.tar.gz
mv petsc-3.4.3/* src/
rm -rf petsc-3.4.3

#configure
cd src
./config/configure.py \
 --prefix="$ISSM_DIR/externalpackages/petsc/install" \
 --PETSC_ARCH="linux-gnu-amd64" \
 --PETSC_DIR="$ISSM_DIR/externalpackages/petsc/src" \
 --with-debugging=0 \
 --with-shared-libraries=1 \
 --known-mpi-shared-libraries=1 \
 --with-mpi-dir=$OPENMPI \
 --with-blas-lapack-dir=/usr/local/intel/composer_xe_2011_sp1.13.367/mkl/ \
 --download-scalapack=yes \
 --download-mumps=yes \
 --download-metis=yes \
 --download-parmetis=yes 

make
make install
