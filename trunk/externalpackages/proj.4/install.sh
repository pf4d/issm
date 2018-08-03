#!/bin/bash
set -eu

#Some cleanup
rm -rf install

#Download trunk
svn checkout http://svn.osgeo.org/metacrs/proj/trunk/proj install

#compile
cd install
./configure --prefix="$ISSM_DIR/externalpackages/proj.4/install"
make 
make install
cd ..
