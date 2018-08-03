#!/bin/bash

if [ -e ~/.bashrc ]; then 
	source ~/.bashrc
fi

# Download fresh copy of MITgcm
cd ../MITgcm/
source install.sh
