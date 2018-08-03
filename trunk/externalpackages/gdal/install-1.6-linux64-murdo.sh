#!/bin/bash
set -eu

#Some cleanup
rm -rf src
rm -rf install
rm -rf gdal-1.6.0
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/gdal-1.6.0.tar.gz' 'gdal-1.6.0.tar.gz'

#Untar 
tar -zxvf  gdal-1.6.0.tar.gz

#Move gdal into src directory
mv gdal-1.6.0/* src
rm -rf gdal-1.6.0

#Configure gdal
cd src
./configure --prefix="$ISSM_DIR/externalpackages/gdal/install" \
	--without-python \
	--without-png \
	--with-netcdf=no \
	--with-jasper=no \
	--without-ld-shared \
	--with-unix-stdio-64=no 

#Patch GDALmake.opt
patch GDALmake.opt ../configs/GDALmake.opt.patch

#Compile and install gdal
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
