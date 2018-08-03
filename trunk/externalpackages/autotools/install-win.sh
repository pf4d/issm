#!/bin/bash
set -eu

rm -rf install
mkdir install

export PATH="$ISSM_DIR/externalpackages/autotools/install/bin:$PATH"

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

# This patch takes care of removing options passed to the linker that casuse 
# the build to fail, as well as changing some flags to match up to Microsoft
# compilers.
patch ./install/share/aclocal/libtool.m4 < ./patches/libtool.m4.patch

# Small change to Automake in order to get the right flags for Microsoft's
# compiler.
patch ./install/bin/automake < ./patches/automake.patch

# This patch is for ar-lib, and removes carriage return characters that cause
# commands to overwrite themselves and be misinterpreted during linking on 
# Cygwin Windows.
patch ./install/share/automake-1.14/ar-lib < ./patches/ar-lib.patch
