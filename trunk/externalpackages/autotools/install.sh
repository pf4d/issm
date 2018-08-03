#!/bin/bash
set -eu

rm -rf install
mkdir install

export PATH="$ISSM_DIR/externalpackages/autotools/install/bin:$PATH"

#install m4
echo " === INSTALLING M4 =="
rm -rf src
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/m4-1.4.17.tar.gz' 'm4-1.4.17.tar.gz'
tar -zxvf m4-1.4.17.tar.gz
mv m4-1.4.17 src
cd src
./configure --prefix="$ISSM_DIR/externalpackages/autotools/install"
make
make install
cd ..

#install autoconf
echo " === INSTALLING AUTOCONF =="
rm -rf src
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/autoconf-2.69.tar.gz' 'autoconf-2.69.tar.gz'
tar -zxvf autoconf-2.69.tar.gz
mv autoconf-2.69 src
cd src 
./configure --prefix="$ISSM_DIR/externalpackages/autotools/install" 
make  
make install
cd ..

#install automake
echo " === INSTALLING AUTOMAKE =="
rm -rf src
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/automake-1.14.tar.gz' 'automake-1.14.tar.gz'
tar -zxvf  automake-1.14.tar.gz
mv automake-1.14 src
cd src 
./configure --prefix="$ISSM_DIR/externalpackages/autotools/install" 
make  
make install
cd ..

#install libtool
echo " === INSTALLING LIBTOOL =="
rm -rf src
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/libtool-2.4.2.tar.gz' 'libtool-2.4.2.tar.gz'
tar -zxvf  libtool-2.4.2.tar.gz
rm libtool-2.4.2.tar.gz
mv libtool-2.4.2 src
cd src 
./configure --prefix="$ISSM_DIR/externalpackages/autotools/install" 
make  
make install
cd ..
