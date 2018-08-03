#!/bin/bash
set -eu

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://python-distribute.org/distribute_setup.py' 'distribute_setup.py'
python distribute_setup.py
