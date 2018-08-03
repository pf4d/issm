#!/bin/bash
set -eu

#Some cleanup
rm -rf install src
mkdir install

#Download latest version
svn co --username gmsh --password gmsh https://geuz.org/svn/gmsh/trunk src

#Configure
cd install
cmake ../src -DCMAKE_INSTALL_PREFIX="$ISSM_DIR/externalpackages/gmsh/install"

#Compile and install
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install
