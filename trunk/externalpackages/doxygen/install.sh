#!/bin/bash
set -eu

#Some cleanup
rm -rf install src
mkdir install

#Download latest version
svn co https://svn.code.sf.net/p/doxygen/code/trunk src

#Configure doxygen
cd src && ./configure --prefix "$ISSM_DIR/externalpackages/doxygen/install"
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi

#Install doxygen
make install
make install_docs
