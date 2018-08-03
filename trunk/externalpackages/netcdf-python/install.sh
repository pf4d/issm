#!/bin/bash
set -eu

#Some cleanup
rm -rf install netCDF4-1.0
mkdir install 

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py "http://issm.jpl.nasa.gov/files/externalpackages/netCDF4-1.0.tar.gz" "netCDF4-1.0.tar.gz"

#Untar 
tar -zxvf  netCDF4-1.0.tar.gz

#for later: 
rm -rf ISSMDIR
echo $ISSM_DIR | sed 's/\//\\\//g' > ISSMDIR
ISSMDIR=`cat ISSMDIR` && rm -rf ISSMDIR

#Move netCDF4 to install directory
rm -rf install/*
mv netCDF4-1.0/* install/
rm -rf netCDF4-1.0

#Configur and compile
cd install
#edit setup.cfg to point to hdf5 and netcdf
cat setup.cfg.template  | sed "s/\#netCDF4_dir = \/usr\/local/netCDF4_dir = $ISSMDIR\/externalpackages\/netcdf\/install/g"  | sed "s/\#HDF5_dir = \/usr\/local/HDF5_DIR = $ISSMDIR\/externalpackages\/hdf5\/install/g" > setup.cfg

python setup.py build 
python setup.py install
