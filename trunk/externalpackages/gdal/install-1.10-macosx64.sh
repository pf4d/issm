#!/bin/bash
set -eu

#Some cleanup
rm -rf src
rm -rf install
rm -rf gdal-1.10.0
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/gdal-1.10.0.tar.gz' 'gdal-1.10.0.tar.gz'

#Untar 
tar -zxvf  gdal-1.10.0.tar.gz

#Move gdal into src directory
mv gdal-1.10.0/* src
rm -rf gdal-1.10.0

#Configure gdal
cd src
./configure \
	--prefix="$ISSM_DIR/externalpackages/gdal/install" \
	--with-python \
	--with-netcdf=no \
	--with-libiconv-prefix="/usr/lib/" \
	--without-hdf5 \
	--with-expat=$ISSM_DIR/externalpackages/expat/install
#Note:
# look for libiconv in /usr/lib/ because there seemed to be a conflict with port's library

#Compile and install gdal
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
