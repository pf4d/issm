AM_CPPFLAGS = @DAKOTAINCL@ @PETSCINCL@ @MPIINCL@ @SPOOLESINCL@ @METISINCL@ @TRIANGLEINCL@ @CHACOINCL@ @SCOTCHINCL@ @SHAPELIBINCL@ @AMPIINCL@
AUTOMAKE_OPTIONS = subdir-objects

EXEEXT=$(JAVASCRIPTWRAPPEREXT)

#define prefix (from http://www.gnu.org/software/autoconf/manual/autoconf-2.67/html_node/Defining-Directories.html)
AM_CPPFLAGS+=  -DISSM_PREFIX='"$(prefix)"'

js_scripts = ${ISSM_DIR}/src/wrappers/TriMesh/TriMesh.js  \
			 ${ISSM_DIR}/src/wrappers/NodeConnectivity/NodeConnectivity.js\
			 ${ISSM_DIR}/src/wrappers/ContourToMesh/ContourToMesh.js\
			 ${ISSM_DIR}/src/wrappers/ElementConnectivity/ElementConnectivity.js\
			 ${ISSM_DIR}/src/wrappers/InterpFromMeshToMesh2d/InterpFromMeshToMesh2d.js\
			 ${ISSM_DIR}/src/wrappers/IssmConfig/IssmConfig.js\
			 ${ISSM_DIR}/src/wrappers/Issm/issm.js

bin_SCRIPTS =  issm-prebin.js

issm-prebin.js: ${js_scripts}
	cat ${js_scripts}  > issm-prebin.js
	
#javascript io{{{
if !WINDOWS
lib_LTLIBRARIES = libISSMJavascript.la
else
noinst_LTLIBRARIES = libISSMJavascript.la
lib_LTLIBRARIES = 
endif

io_sources=   ./io/WriteJavascriptData.cpp\
				./io/FetchJavascriptData.cpp

ALLCXXFLAGS= -fPIC -D_WRAPPERS_ $(CXXFLAGS) $(CXXOPTFLAGS) 

libISSMJavascript_la_SOURCES = $(io_sources)
libISSMJavascript_la_CXXFLAGS= $(ALLCXXFLAGS)
#}}}
#api io{{{
if !WINDOWS
lib_LTLIBRARIES += libISSMApi.la
else
noinst_LTLIBRARIES += libISSMApi.la
endif

api_sources= ./io/ApiPrintf.cpp

libISSMApi_la_SOURCES = $(api_sources)
libISSMApi_la_CXXFLAGS= $(ALLCXXFLAGS)
#}}}
#Wrappers {{{
bin_PROGRAMS = 		 IssmModule
#}}}

# Dependencies {{{

#Triangle library
AM_CXXFLAGS =  -DTRILIBRARY -DANSI_DECLARATORS -DNO_TIMER -D_WRAPPERS_
AM_CXXFLAGS +=  -D_HAVE_JAVASCRIPT_MODULES_ -fPIC

deps = ./libISSMJavascript.la ../../c/libISSMModules.la ../../c/libISSMCore.la ./libISSMApi.la

#Optimization flags:
AM_CXXFLAGS += $(CXXOPTFLAGS) 
#}}}
# Module sources and dependencies {{{
if !WINDOWS
libISSMJavascript_la_LIBADD = ./../../c/libISSMCore.la ./../../c/libISSMModules.la $(MPILIB) $(PETSCLIB) $(MKLLIB) $(GSLLIB) $(PROJ4LIB) $(MATHLIB) 
endif

if STANDALONE_LIBRARIES
libISSMJavascript_la_LDFLAGS = -static 
deps += $(PETSCLIB) $(TAOLIB) $(M1QN3LIB) $(PLAPACKLIB) $(MUMPSLIB) $(SUPERLULIB) $(SPOOLESLIB) $(TRIANGLELIB) $(SCALAPACKLIB) $(BLACSLIB) $(HYPRELIB) $(SPAILIB) $(PROMETHEUSLIB) $(PASTIXLIB) $(MLLIB) $(DAKOTALIB) $(METISLIB) $(CHACOLIB) $(SCOTCHLIB) $(BLASLAPACKLIB) $(MKLLIB) $(MPILIB) $(MATHLIB) $(GRAPHICSLIB) $(MULTITHREADINGLIB) $(OSLIBS) $(GSLLIB)   $(ADOLCLIB) $(AMPILIB) $(METEOIOLIB) $(SNOWPACKLIB)
endif

if !WINDOWS
libISSMApi_la_LIBADD = $(MPILIB) $(PETSCLIB) $(GSLLIB) $(PROJ4LIB) $(MATHLIB) $(MEXLIB)
endif

if STANDALONE_LIBRARIES
libISSMApi_la_LDFLAGS = -static 
endif

IssmModule_SOURCES = ../TriMesh/TriMesh.cpp \
					 ../NodeConnectivity/NodeConnectivity.cpp\
					 ../ContourToMesh/ContourToMesh.cpp\
					 ../ElementConnectivity/ElementConnectivity.cpp\
					 ../InterpFromMeshToMesh2d/InterpFromMeshToMesh2d.cpp\
					 ../IssmConfig/IssmConfig.cpp\
					 ../Issm/issm.cpp

IssmModule_CXXFLAGS= -fPIC -D_DO_NOT_LOAD_GLOBALS_  --memory-init-file 0 $(AM_CXXFLAGS) $(CXXFLAGS) $(CXXOPTFLAGS) $(COPTFLAGS) -s EXPORTED_FUNCTIONS="['_TriMeshModule','_NodeConnectivityModule','_ContourToMeshModule','_ElementConnectivityModule','_InterpFromMeshToMesh2dModule','_IssmConfigModule','_IssmModule']"  -s DISABLE_EXCEPTION_CATCHING=0 -s ALLOW_MEMORY_GROWTH=1 -s INVOKE_RUN=0
IssmModule_LDADD = ${deps} $(TRIANGLELIB)  $(GSLLIB)
#}}}
