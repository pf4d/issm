#!/bin/bash
set -eu

#Some cleanup
rm -rf src install hdf5-1.8.9
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/hdf5-1.8.9.tar.gz' 'hdf5-1.8.9.tar.gz'

#Untar 
tar -zxvf  hdf5-1.8.9.tar.gz

#Move hdf5 to src directory
rm -rf src/*
mv hdf5-1.8.9/* src/
rm -rf hdf5-1.8.9

# This project uses C code with C++ style comment default C standard used by 
# GNU's C compiler's default C standard  does not support C++ style comments.
# As such, we need to specify a standard that does.
export CFLAGS='-std=c99'

#Configure and compile
cd src
./configure  --prefix="$ISSM_DIR/externalpackages/hdf5/install" 
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
