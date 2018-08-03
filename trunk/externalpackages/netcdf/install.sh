#!/bin/bash
set -eu
#you need hdf5 compiled

#Some cleanup
rm -rf src install netcdf-4.3.2
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py "http://issm.jpl.nasa.gov/files/externalpackages/netcdf-4.3.2.tar.gz" "netcdf-4.3.2.tar.gz"

#Untar 
tar -zxvf  netcdf-4.3.2.tar.gz

#Move netcdf to install directory
rm -rf src/*
mv netcdf-4.3.2/* src/
rm -rf netcdf-4.3.2

#Configure and compile
cd src
./configure \
 --prefix="$ISSM_DIR/externalpackages/netcdf/install"  \
 --disable-doxygen
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
