#!/bin/sh

./configure \
	--prefix=$ISSM_DIR \
	--without-wrappers\
	--with-gsl-dir=$ISSM_DIR/externalpackages/gsl/install\
	--without-thermal \
	--without-control \
	--without-hydrology \
	--without-stressbalance \
	--without-balanced \
	--without-responses \
	--without-slope \
	--without-rifts \
	--without-steadystate \
	--without-transient \
	--without-3d \
	--without-groundingline
