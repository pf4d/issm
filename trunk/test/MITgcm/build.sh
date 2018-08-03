#!/bin/bash
#This script compiles and links MITgcm

#recover hostname and model path:
hostname="$1"
modelpath="$2"

if [ -e ~/.bashrc ]; then
	source ~/.bashrc
fi

# Clean up build directory
cd $modelpath
if [ ! -d "build" ]; then mkdir build; fi
\rm -f build/*

# Get MITgcm code, if needed
if [ ! -d "$SLR_DIR/components/mitgcm/install" ]; then
   cd $modelpath/../MITgcm
   source install.sh
   cd $modelpath
fi

#create MITgcm makefile and code links for this run
cd build
if [ ! -f Makefile ]; then

	if [ $hostname == "pleiades" ]; then 
		$SLR_DIR/components/mitgcm/install/tools/genmake2 -of $SLR_DIR/models/ice-ocean/configs/linux_amd64_gfortran+mpi_ice_nas -mo ../code -rd $SLR_DIR/components/mitgcm/install
	else
		$modelpath/../MITgcm/install/tools/genmake2 -mpi -mo $modelpath/../MITgcm/code -rd $modelpath/../MITgcm/install
	fi
    make depend
fi
make -j 8
