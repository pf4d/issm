#!/bin/bash

#comprehensive check, except unusedFunction
#cppcheck -j 32 --include=$ISSM_DIR/config.h -DHAVE_CONFIG_H -D_HAVE_ADOLC_ -D_HAVE_DAKOTA_ --enable=all $ISSM_DIR/src/c 2> CPPCHECK.err

#unused function only (slow)
cppcheck --include=$ISSM_DIR/config.h -DHAVE_CONFIG_H -D_HAVE_ADOLC_ -D_HAVE_DAKOTA_ --enable=unusedFunction $ISSM_DIR/src 2> CPPCHECK.err
