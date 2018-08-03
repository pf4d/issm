#!/bin/bash
set -eu

#Erase install
rm -rf install  findbugs-1.3.9
mkdir install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/findbugs-1.3.9.tar.gz' 'findbugs-1.3.9.tar.gz'
tar -zxvf findbugs-1.3.9.tar.gz 

mv findbugs-1.3.9/* install
rm -rf findbugs-1.3.9
