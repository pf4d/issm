#!/bin/bash
set -eu

export ESMF_DIR=$ISSM_DIR/externalpackages/esmf/esmf
export ESMF_INSTALL_PREFIX=$ISSM_DIR/externalpackages/esmf/install
export ESMF_COMPILER=gfortran
export ESMF_COMM=mpich2

#Some cleanup
rm -rf esmf_6_3_0rp1
rm -rf esmf

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/esmf_6_3_0rp1_src.tar.gz' 'esmf_6_3_0rp1_src.tar.gz'
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/docs/ESMF_refdoc.pdf' 'ESMF_refdoc.pdf'
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/docs/ESMF_usrdoc.pdf' 'ESMF_usrdoc.pdf'

#Untar 
tar -zxvf  esmf_6_3_0rp1_src.tar.gz

#Configure esmf
cd esmf

#Compile and install esmf
if [ $# -eq 0 ]; then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
