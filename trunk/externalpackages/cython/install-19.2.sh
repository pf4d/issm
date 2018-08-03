#!/bin/bash
set -eu

#clean up
rm -rf Cython-0.19.2

#download numpy first
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/Cython-0.19.2.tar.gz' 'cython.tar.gz'

#install numpy
tar -zxvf cython.tar.gz
cd Cython-0.19.2
python setup.py install
