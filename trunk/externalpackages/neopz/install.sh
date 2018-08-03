#!/bin/bash
set -eu

#Erase install
rm -rf install
mkdir install

#Download from GitHub server
#GIT must be installed first. See $ISSM_DIR/externalpackages/git
git clone https://github.com/labmec/neopz.git

#Untar and set src directory
mv neopz/ install/

#Set neopz CMake variables
#CMake must be installed first. See $ISSM_DIR/externalpackages/git
export PROJECT_SOURCE_DIR=$ISSM_DIR/externalpackages/neopz/install/neopz
export PROJECT_BINARY_DIR=$ISSM_DIR/externalpackages/neopz/install/

#Configure neopz using cmake
cd $PROJECT_SOURCE_DIR
cmake -DCMAKE_INSTALL_PREFIX:PATH=$PROJECT_BINARY_DIR

cd $PROJECT_SOURCE_DIR
#Compile and install 
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
cd $PROJECT_BINARY_DIR/pzlib
mv lib ../
mv include ../
cd ..
rm -rf pzlib
cd ..
