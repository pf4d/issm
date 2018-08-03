#!/bin/bash
set -eu

#Some cleanup
rm -rf src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/emsdk-portable.tar.gz' 'emsdk-portable.tar.gz'

#Untar and move petsc to install directory
tar -zxvf  emsdk-portable.tar.gz
mv emsdk_portable src

cd src

export CXX=g++
export CC=gcc

./emsdk update
./emsdk install latest
./emsdk activate latest
