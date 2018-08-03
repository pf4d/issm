#!/bin/bash
set -eu

rm -rf gssh-gmt-2.3.4.tar.gz  src install

#get gssh database from noaa's website:  http://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html
#curl http://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-gmt-2.3.4.tar.gz > gshhg-gmt-2.3.4.tar.gz

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/gshhg-gmt-2.3.4.tar.gz' 'gshhg-gmt-2.3.4.tar.gz'

#untar: 
tar -zxvf gshhg-gmt-2.3.4.tar.gz 

#move: 
mv gshhg-gmt-2.3.4 install
