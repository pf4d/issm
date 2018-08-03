#!/bin/bash
set -eu

#Some cleanup
rm -rf install
rm -rf swig-2.0.4
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/swig-2.0.4.tar.gz' 'swig-2.0.4.tar.gz'

#Untar and move python into install directory
tar -zxvf  swig-2.0.4.tar.gz
mv swig-2.0.4/* install
rm -rf swig-2.0.4

#Copy pcre prototype in include directory
cd install 
#cp $ISSM_DIR/externalpackages/pcre/install/pcre.h  ./Source/Include/
#cp $ISSM_DIR/externalpackages/pcre/install/.libs/*  ./Source/Include/
export CFLAGS="-I$ISSM_DIR/externalpackages/pcre/install"
export LDFLAGS="-L$ISSM_DIR/externalpackages/pcre/install/.libs/"
export LIBS="-lpcre"
#Configure swig
./configure \
 --prefix="$ISSM_DIR/externalpackages/swig/install" \
 --with-pcre-prefix="$ISSM_DIR/externalpackages/pcre/install" \
 --with-pcre-exec-prefix="$ISSM_DIR/externalpackages/pcre/install" \
 --with-python="$ISSM_DIR/externalpackages/python/install/bin/python"

#Compile and install gdal
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
