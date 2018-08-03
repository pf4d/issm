#!/bin/sh
#  this shell script will automatically reconfigure the entire ISSM 
#  archive.

# As it turns out, the autoreconf script provided by Autotools
# encompasses the functionality of this script. As such, the 
# following two lines will replace the remainder of the script.
# If all goes well, then the script will be shortened in the future.

cd $ISSM_DIR
autoreconf -iv
