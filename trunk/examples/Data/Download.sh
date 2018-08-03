#!/bin/bash

echo "Downloading Square shelf dataset"
$ISSM_DIR/scripts/DownloadExternalPackage2.py 'http://issm.jpl.nasa.gov/files/workshop2014/SquareShelf.nc' 'SquareShelf.nc'

echo "Downloading SeaRISE dataset - Antarctica"
$ISSM_DIR/scripts/DownloadExternalPackage2.py 'http://websrv.cs.umt.edu/isis/images/c/cc/Antarctica_5km_withshelves_v0.75.nc' 'Antarctica_5km_withshelves_v0.75.nc'

echo "Downloading InSAR Antarctic velocities"
$ISSM_DIR/scripts/DownloadExternalPackage2.py 'ftp://n5eil01u.ecs.nsidc.org/SAN/MEASURES/NSIDC-0484.001/1996.01.01/antarctica_ice_velocity_900m.nc' 'Antarctica_ice_velocity.nc' 

echo "Downloading PIG errors"
$ISSM_DIR/scripts/DownloadExternalPackage2.py 'http://issm.jpl.nasa.gov/files/workshop2014/CrossOvers2009.mat' 'CrossOvers2009.mat'

echo "Downloading SeaRISE dataset - Greenland"
$ISSM_DIR/scripts/DownloadExternalPackage2.py 'http://websrv.cs.umt.edu/isis/images/e/e9/Greenland_5km_dev1.2.nc' 'Greenland_5km_dev1.2.nc'

echo "Downloading Jason Box's SMB"
$ISSM_DIR/scripts/DownloadExternalPackage2.py 'http://issm.jpl.nasa.gov/files/examples/Box_Greenland_SMB_monthly_1840-2012_5km_cal_ver20141007.nc' 'Box_Greenland_SMB_monthly_1840-2012_5km_cal_ver20141007.nc'

echo "Downloading IceBridge Jakobshavn bedrock"
$ISSM_DIR/scripts/DownloadExternalPackage2.py 'https://data.cresis.ku.edu/data/grids/old_versions/Jakobshavn_2008_2011_Composite.zip' 'Jakobshavn_2008_2011_Composite.zip'
unzip Jakobshavn_2008_2011_Composite.zip
mv Jakobshavn_2008_2011_Composite/grids/Jakobshavn_2008_2011_Composite_XYZGrid.txt .
rm -rf Jakobshavn_2008_2011_Composite Jakobshavn_2008_2011_Composite.zip
