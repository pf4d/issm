ALLCXXFLAGS= -fPIC $(CXXFLAGS) $(CXXOPTFLAGS) 

AM_CPPFLAGS =  -DNO_TIMER ${CPP_LINUX} -DTRILIBRARY

lib_LTLIBRARIES = libtriangle.la

if EXECUTABLES
bin_PROGRAMS = tricall showme
endif 

include_HEADERS= ./triangle.h

libtriangle_la_SOURCES = ./triangle.c
libtriangle_la_LIBADD = -lm

if EXECUTABLES
tricall_SOURCES = ./tricall.c
tricall_LDADD = ./libtriangle.la

showme_SOURCES = ./showme.c
showme_LDADD = -lX11
endif
