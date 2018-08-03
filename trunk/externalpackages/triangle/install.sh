#!/bin/bash
set -eu

#Some cleanup 
rm -rf install src triangle
mkdir install src ./src/m4

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/triangle.zip' 'triangle.zip'

#Untar 
cd src
mkdir src
cp ../triangle.zip ./
unzip triangle.zip

rm ./makefile
mv ./*.c ./src
mv ./*.h ./src

cp ./../configs/libtool/configure.ac ./
cp ./../configs/libtool/Makefile.am ./
cp ./../configs/libtool/src/Makefile.am ./src/

autoreconf -ivf
./configure --prefix="${HOME}/externalpackages/triangle/install" --disable-executables

make 
make install

#Patch triangle.h
patch ${HOME}/externalpackages/triangle/install/include/triangle.h ${HOME}/externalpackages/triangle/triangle.h.patch
