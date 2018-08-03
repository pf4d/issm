#!/bin/bash
set -eu

#Some cleanup
rm -rf src install curl-7.39.0
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/curl-7.39.0.tar.gz' 'curl-7.39.0.tar.gz'

#Untar 
tar -zxvf  curl-7.39.0.tar.gz

#Move curl into src directory
mv curl-7.39.0/* src
rm -rf curl-7.39.0

#Configure curl
cd src

export CFLAGS=" -arch x86_64"

./configure \
	--prefix="$ISSM_DIR/externalpackages/curl/install" 

#Compile curl
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install 
