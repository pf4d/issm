#!/bin/sh

./configure \
   --prefix=$ISSM_DIR\
   --libdir=$ISSM_DIR/../mobile/android/ISSM_APP/libs/armeabi\
   --build="i386-apple-darwin10.8.0"\
   --host="arm-linux-androideabi"\
   --enable-shared\
   --with-android=jni\
   --with-android-ndk=$ISSM_DIR/externalpackages/android/android-ndk/install\
   --without-fortran\
   --without-wrappers\
   --without-kriging\
   --with-gsl-dir=$ISSM_DIR/externalpackages/gsl/install
