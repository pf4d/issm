#!/bin/bash
set -eu

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/win7.sdk7.1.exe' 'win7.sdk7.1.exe'
