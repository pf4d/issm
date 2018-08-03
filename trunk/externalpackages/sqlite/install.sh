#!/bin/bash
set -eu

#Some cleanup
rm -rf src
rm -rf install
rm -rf sqlite-autoconf-3071300
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/sqlite-autoconf-3071300.tar.gz' 'sqlite-autoconf-3071300.tar.gz'

#Untar 
tar -zxvf  sqlite-autoconf-3071300.tar.gz

#Move sqlite-autoconf into src directory
mv sqlite-autoconf-3071300/* src
rm -rf sqlite-autoconf-3071300

#Configure sqlite-autoconf
cd src
./configure  --prefix="$ISSM_DIR/externalpackages/sqlite/install" 

#Compile and install sqlite-autoconf
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
