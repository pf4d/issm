#!/bin/bash
set -eu

# Some cleanup
rm -rf Chaco-2.2
rm -rf src 
rm -rf install 
mkdir src install 

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/Chaco-2.2.tar.gz' 'Chaco-2.2.tar.gz'
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/docs/chaco_guide.pdf' 'chaco_guide.pdf'

# Untar 
tar -xvzf Chaco-2.2.tar.gz

# Move chaco to src directory
mv Chaco-2.2/* src
rm -rf Chaco-2.2

# Apply patches (all at once)
# (written by diff -rc src ~/Libs/Chaco-2.2 > chaco.patch)
patch -R -p0 < chaco.patch

# Patch src/code/Makefile
patch ./src/code/Makefile ./patches/Makefile.patch


# Build chaco
cd src/code
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make chacominusblas.a

# Clean up objects (but not library or executable)
make clean
cd ../..

# Populate install directory
cp -p src/exec/README install
cp -p src/exec/User_Params install
cp -p src/exec/*.coords install
cp -p src/exec/*.graph install
mkdir install/include
cp -p src/code/main/defs.h install/include/defs.h
cp -p src/code/main/params.h install/include/params.h
cp -p chaco.h install/include/chaco.h
mkdir install/lib
mv src/code/chaco.a install/lib/libchaco.a
mv src/code/chacominusblas.a install/lib/libchacominusblas.a
mkdir install/exec
mv src/exec/chaco install/exec
