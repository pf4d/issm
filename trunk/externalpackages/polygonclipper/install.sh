#!/bin/bash
set -eu

#Some cleanup
rm -rf install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/PolygonClipper.zip' 'PolygonClipper.zip'

#install
mkdir install
cd install
cp ../PolygonClipper.zip .

#uncompress
unzip PolygonClipper.zip

#Make
mex gpc.c gpc_mexfile.c -O -output PolygonClip
