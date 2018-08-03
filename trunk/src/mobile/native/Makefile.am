#Generate JNI library
lib_LTLIBRARIES = libIssmJni.la

libIssmJni_la_SOURCES = Main.cpp
libIssmJni_la_LIBADD = $(ISSM_DIR)/src/c/libISSMCore.a $(ISSM_DIR)/externalpackages/gsl/install/lib/libgsl.a $(ISSM_DIR)/externalpackages/gsl/install/lib/libgslcblas.a

AM_LDFLAGS = -llog
