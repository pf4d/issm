#!/bin/bash
set -eu 
VER="2.8.11.2"

#Some cleanup
rm -rf install cmake-$VER
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py "http://issm.jpl.nasa.gov/files/externalpackages/cmake-$VER.tar.gz" "cmake-$VER.tar.gz"

#Untar 
tar -zxvf  cmake-$VER.tar.gz

#Move cmake into install directory
mv cmake-$VER/* install
rm -rf cmake-$VER

LDFLAGS="-L/usr/lib/ -lstdc++ -L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin10/4.6.2/ -lgfortran"
F77="/usr/local/gfortran/bin/x86_64-apple-darwin10-gfortran"

#Compile cmake
cd install 
./bootstrap --prefix=$ISSM_DIR/externalpackages/cmake/install
if [ $# -eq 0 ]; then
	make
else 
	make -j $1; 
fi
make install
