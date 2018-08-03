#!/bin/bash
set -eu

#Some cleanup 
rm -rf install-javascript triangle
mkdir install-javascript

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/triangle.zip' 'triangle.zip'

#Untar 
cd install-javascript
cp ../triangle.zip ./
unzip triangle.zip

#copy new makefile
cp ../configs/javascript/configure.make ./
cp ../configs/javascript/makefile ./

#Patch triangle.h
patch triangle.h ../triangle.h.patch.js

#Compile triangle
make


