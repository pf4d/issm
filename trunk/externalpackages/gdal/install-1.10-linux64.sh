#!/bin/bash
set -eu

#WARNING: you need to have python installed in externalpackages

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

export CFLAGS=-D_HAVE_STRNDUP
export CXXFLAGS=-D_HAVE_STRNDUP

#Configure gdal
cd src
./configure \
	--prefix="$ISSM_DIR/externalpackages/gdal/install" \
	--with-python="$ISSM_DIR/externalpackages/python/install/bin/python" \
	--with-netcdf=no \
	--with-jasper=no \
	--without-hdf5

#Compile and install gdal
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install

#For some reasons, on thwaites, one needs to do the following to get the python bindings:
#cd src/swig/python/
# python setup.py build
# python setup.py install
