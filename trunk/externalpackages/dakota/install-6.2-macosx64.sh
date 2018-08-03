#!/bin/bash
set -eu

#Some cleanup
rm -rf Dakota
rm -rf src 
rm -rf build 
rm -rf install 
mkdir src build install 

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/dakota-6.2-public.src.tar.gz' 'dakota-6.2-public-src.tar.gz'

#Untar 
tar -zxvf dakota-6.2-public-src.tar.gz

#Move Dakota to src directory
mv dakota-6.2.0.src/* src
rm -rf dakota-6.2.0.src

#Set up Dakota cmake variables and config
export DAK_SRC=$ISSM_DIR/externalpackages/dakota/src
export DAK_BUILD=$ISSM_DIR/externalpackages/dakota/build
export MPIHOME=$ISSM_DIR/externalpackages/mpich/install
cp $DAK_SRC/cmake/BuildDakotaTemplate.cmake $DAK_SRC/cmake/BuildDakotaCustom.cmake
patch $DAK_SRC/cmake/BuildDakotaCustom.cmake configs/6.2/BuildDakotaCustom.cmake.mac.patch
patch $DAK_SRC/cmake/DakotaDev.cmake configs/6.2/DakotaDev.cmake.patch
patch $DAK_SRC/CMakeLists.txt configs/6.2/CMakeLists.txt.patch

#Apply patches
patch src/src/NonDSampling.cpp configs/6.2/NonDSampling.cpp.patch
patch src/src/NonDLocalReliability.cpp configs/6.2/NonDLocalReliability.cpp.patch
patch src/packages/pecos/src/pecos_global_defs.hpp configs/6.2/pecos_global_defs.hpp.patch

export BOOST_ROOT=$ISSM_DIR/externalpackages/boost/install

#Configure dakota
# Set your local gcc compiler here
cd $DAK_BUILD
cmake -DBoost_NO_BOOST_CMAKE=TRUE \
	-DBoost_NO_SYSTEM_PATHS=TRUE \
	-DBOOST_ROOT:PATHNAME=$BOOST_ROOT \
	-DBoost_LIBRARY_DIRS:FILEPATH=${BOOST_ROOT}/lib \
	-D CMAKE_C_COMPILER=$ISSM_DIR/externalpackages/mpich/install/bin/mpicc \
	-D CMAKE_CXX_COMPILER=$ISSM_DIR/externalpackages/mpich/install/bin/mpicxx \
	-D CMAKE_Fortran_COMPILER=$ISSM_DIR/externalpackages/mpich/install/bin/mpif77 \
	-D CMAKE_CXX_FLAGS=-fdelayed-template-parsing \
	-DHAVE_ACRO=off \
	-DHAVE_JEGA=off \
	-C $DAK_SRC/cmake/BuildDakotaCustom.cmake \
	-C $DAK_SRC/cmake/DakotaDev.cmake \
	$DAK_SRC
cd ..

# Snowleopard: Mpi should be made with these compilers
#-DCMAKE_CXX_COMPILER=/usr/bin/g++ -DCMAKE_CC_COMPILER=/usr/bin/gcc \
#-DCMAKE_Fortran_COMPILER=/usr/local/gfortran/bin/x86_64-apple-darwin10-gfortran \

#Compile and install dakota
cd $DAK_BUILD
if [ $# -eq 0 ];
then
	make
	make install
else
	make -j $1
	make -j $1 install
fi
cd ..
