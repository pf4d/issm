#!/bin/bash
set -eu

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/kml_shapefile.zip' 'kml_shapefile.zip'
