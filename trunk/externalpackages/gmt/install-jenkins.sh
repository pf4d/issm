#!/bin/bash
set -eu

#Erase install
rm -rf install  src gmt

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/gmt-5.1.1.tar.gz' 'gmt-5.1.1.tar.gz'

#install directory
mkdir src
tar -zxvf gmt-5.1.1.tar.gz 
mv gmt-5.1.1/* src
rm -rf gmt-5.1.1

#configure: 
cp configs/ConfigUser.cmake-jenkins ./src/cmake/ConfigUser.cmake

cd src
mkdir build
cd build
cmake ../

#compile
if [ $# -eq 0 ]; then
	make install
else
	make -j $1 install
fi

#come back: 
cd ../../
