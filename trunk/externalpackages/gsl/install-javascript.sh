#!/bin/bash
set -eu

#Some cleanup
rm -rf src-javascript install-javascript gsl-1.15
mkdir src-javascript install-javascript

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/gsl-1.15.tar.gz' 'gsl-1.15.tar.gz'

#Untar 
tar -zxvf  gsl-1.15.tar.gz

#Move gsl into src directory
mv gsl-1.15/* src-javascript
rm -rf gsl-1.15

#Configure gsl
cd src-javascript

export CC=emcc
export CXX=em++

# Issue with variadic function signatures.
export CFLAGS=-DSTDC_HEADERS

./configure --prefix="$ISSM_DIR/externalpackages/gsl/install-javascript" 

#Compile gsl
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install 
