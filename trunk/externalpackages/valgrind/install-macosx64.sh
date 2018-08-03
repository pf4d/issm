#!/bin/bash
set -eu

#Some cleanup
rm -rf install

#Download development version, the current release never supports the latest OS X releases
svn co svn://svn.valgrind.org/valgrind/trunk install

#configure
cd install
./autogen.sh
./configure  --prefix="$ISSM_DIR/externalpackages/valgrind/install" --enable-only64bit

#Compile valgrind

make  -j 8
make install
