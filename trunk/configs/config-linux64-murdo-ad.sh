#!/bin/sh

./configure \
	 --prefix=$ISSM_DIR \
	 --without-kriging \
	 --without-kml \
	--with-gsl-dir=$ISSM_DIR/externalpackages/gsl/install \
	--with-adolc-dir=$ISSM_DIR/externalpackages/adolc/install \
	--with-matlab-dir="$ISSM_DIR/externalpackages/matlab/install" \
	--with-triangle-dir=$ISSM_DIR/externalpackages/triangle/install
