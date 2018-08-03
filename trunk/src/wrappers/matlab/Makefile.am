AM_CPPFLAGS = @NEOPZINCL@ @DAKOTAINCL@ @MATLABINCL@ @PETSCINCL@ @MPIINCL@ @SPOOLESINCL@ @METISINCL@ @TRIANGLEINCL@ @CHACOINCL@ @SCOTCHINCL@ @SHAPELIBINCL@ @AMPIINCL@
AUTOMAKE_OPTIONS = subdir-objects

EXEEXT=$(MATLABWRAPPEREXT)

#define prefix (from http://www.gnu.org/software/autoconf/manual/autoconf-2.67/html_node/Defining-Directories.html)
AM_CPPFLAGS+=  -DISSM_PREFIX='"$(prefix)"'

#matlab io{{{
if !WINDOWS
lib_LTLIBRARIES = libISSMMatlab.la
else
noinst_LTLIBRARIES = libISSMMatlab.la
lib_LTLIBRARIES = 
endif

io_sources= ./io/CheckNumMatlabArguments.cpp\
				./io/WriteMatlabData.cpp\
				./io/FetchMatlabData.cpp

ALLCXXFLAGS= -fPIC -D_WRAPPERS_ $(CXXFLAGS) $(CXXOPTFLAGS) 

libISSMMatlab_la_SOURCES = $(io_sources)
libISSMMatlab_la_CXXFLAGS= $(ALLCXXFLAGS)
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
lib_LTLIBRARIES += BamgMesher_matlab.la\
						 BamgConvertMesh_matlab.la\
						 BamgTriangulate_matlab.la\
						 ContourToMesh_matlab.la\
						 ContourToNodes_matlab.la\
						 DistanceToMaskBoundary_matlab.la\
						 ElementConnectivity_matlab.la\
						 ExpSimplify_matlab.la\
						 ExpToLevelSet_matlab.la\
						 InterpFromGridToMesh_matlab.la\
						 InterpFromMeshToMesh2d_matlab.la\
						 InterpFromMeshToMesh3d_matlab.la\
						 InterpFromMeshToGrid_matlab.la\
						 InterpFromMesh2d_matlab.la\
						 IssmConfig_matlab.la\
						 NodeConnectivity_matlab.la\
						 M1qn3_matlab.la\
						 MeshPartition_matlab.la\
						 MeshProfileIntersection_matlab.la\
						 PointCloudFindNeighbors_matlab.la\
						 PropagateFlagsFromConnectivity_matlab.la\
						 TriMesh_matlab.la\
						 TriMeshProcessRifts_matlab.la\
						 Scotch_matlab.la

if CHACO
lib_LTLIBRARIES += Chaco_matlab.la
endif
if KRIGING
lib_LTLIBRARIES +=  Kriging_matlab.la
endif
if KML
lib_LTLIBRARIES += ShpRead_matlab.la
endif
#}}}

# Dependencies {{{
deps =  $(MATHLIB)

#Triangle library
AM_CXXFLAGS =  -DTRILIBRARY -DANSI_DECLARATORS -DNO_TIMER -D_WRAPPERS_

#Matlab part
AM_LDFLAGS = -module $(MEXLINK) -shrext ${EXEEXT} --export-dynamic -rdynamic -no-undefined
if VERSION
AM_LDFLAGS +=
else
AM_LDFLAGS += -avoid-version
endif

if STANDALONE_MODULES
AM_LDFLAGS += -Wl,-static
endif

AM_CXXFLAGS +=  -D_HAVE_MATLAB_MODULES_ -fPIC

deps += ./libISSMMatlab.la ../../c/libISSMModules.la ../../c/libISSMCore.la ./libISSMApi.la

if WINDOWS
deps += $(METISLIB)
endif

if ADOLC
deps += $(ADOLCLIB)
endif

if FORTRAN
deps += $(FLIBS) $(FORTRANLIB)
endif

#MEXLIB needs to be the last (undefined references on larsen)
deps += ${MEXLIB}

#Optimization flags:
AM_CXXFLAGS += $(CXXOPTFLAGS) 
#}}}
# Module sources and dependencies {{{
if !WINDOWS
libISSMMatlab_la_LIBADD = ./../../c/libISSMCore.la ./../../c/libISSMModules.la $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(MKLLIB) $(GSLLIB) $(PROJ4LIB) $(MATHLIB) $(MEXLIB) 
endif

if STANDALONE_LIBRARIES
libISSMMatlab_la_LDFLAGS = -static 
deps += $(PETSCLIB) $(TAOLIB) $(NEOPZLIB) $(M1QN3LIB) $(PLAPACKLIB) $(MUMPSLIB) $(SUPERLULIB) $(SPOOLESLIB) $(TRIANGLELIB) $(SCALAPACKLIB) $(BLACSLIB) $(HYPRELIB) $(SPAILIB) $(PROMETHEUSLIB) $(PASTIXLIB) $(MLLIB) $(DAKOTALIB) $(METISLIB) $(CHACOLIB) $(SCOTCHLIB) $(BLASLAPACKLIB) $(MKLLIB) $(MPILIB) $(MATHLIB) $(GRAPHICSLIB) $(MULTITHREADINGLIB) $(OSLIBS) $(GSLLIB)   $(ADOLCLIB) $(AMPILIB) $(METEOIOLIB) $(SNOWPACKLIB)
endif

if !WINDOWS
libISSMApi_la_LIBADD = $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(GSLLIB) $(PROJ4LIB) $(MATHLIB) $(MEXLIB)
endif

if STANDALONE_LIBRARIES
libISSMApi_la_LDFLAGS = -static 
endif

BamgMesher_matlab_la_SOURCES = ../BamgMesher/BamgMesher.cpp
BamgMesher_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
BamgMesher_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(GSLLIB) $(PROJ4LIB)

BamgConvertMesh_matlab_la_SOURCES = ../BamgConvertMesh/BamgConvertMesh.cpp
BamgConvertMesh_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
BamgConvertMesh_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(GSLLIB) $(PROJ4LIB)

BamgTriangulate_matlab_la_SOURCES = ../BamgTriangulate/BamgTriangulate.cpp
BamgTriangulate_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
BamgTriangulate_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(GSLLIB) $(PROJ4LIB)

if CHACO
Chaco_matlab_la_SOURCES = ../Chaco/Chaco.cpp
Chaco_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
Chaco_matlab_la_LIBADD = ${deps} $(MPILIB) $(CHACOLIB) $(GSLLIB) $(PROJ4LIB) $(PETSCLIB) $(NEOPZLIB)
endif

ContourToMesh_matlab_la_SOURCES = ../ContourToMesh/ContourToMesh.cpp
ContourToMesh_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
ContourToMesh_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(MULTITHREADINGLIB) $(GSLLIB) $(PROJ4LIB)

ExpToLevelSet_matlab_la_SOURCES = ../ExpToLevelSet/ExpToLevelSet.cpp
ExpToLevelSet_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
ExpToLevelSet_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(MULTITHREADINGLIB) $(GSLLIB) $(PROJ4LIB)

ContourToNodes_matlab_la_SOURCES = ../ContourToNodes/ContourToNodes.cpp
ContourToNodes_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
ContourToNodes_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(GSLLIB) $(PROJ4LIB)

DistanceToMaskBoundary_matlab_la_SOURCES = ../DistanceToMaskBoundary/DistanceToMaskBoundary.cpp
DistanceToMaskBoundary_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
DistanceToMaskBoundary_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(GSLLIB) $(PROJ4LIB)

ElementConnectivity_matlab_la_SOURCES = ../ElementConnectivity/ElementConnectivity.cpp
ElementConnectivity_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
ElementConnectivity_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(GSLLIB) $(PROJ4LIB)

InterpFromGridToMesh_matlab_la_SOURCES = ../InterpFromGridToMesh/InterpFromGridToMesh.cpp
InterpFromGridToMesh_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
InterpFromGridToMesh_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(MULTITHREADINGLIB) $(GSLLIB) $(PROJ4LIB)

InterpFromMeshToMesh2d_matlab_la_SOURCES = ../InterpFromMeshToMesh2d/InterpFromMeshToMesh2d.cpp
InterpFromMeshToMesh2d_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
InterpFromMeshToMesh2d_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(MULTITHREADINGLIB) $(GSLLIB) $(PROJ4LIB)

InterpFromMeshToMesh3d_matlab_la_SOURCES = ../InterpFromMeshToMesh3d/InterpFromMeshToMesh3d.cpp
InterpFromMeshToMesh3d_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
InterpFromMeshToMesh3d_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(GSLLIB) $(PROJ4LIB)

InterpFromMeshToGrid_matlab_la_SOURCES = ../InterpFromMeshToGrid/InterpFromMeshToGrid.cpp
InterpFromMeshToGrid_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
InterpFromMeshToGrid_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(GSLLIB) $(PROJ4LIB)

InterpFromMesh2d_matlab_la_SOURCES = ../InterpFromMesh2d/InterpFromMesh2d.cpp
InterpFromMesh2d_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
InterpFromMesh2d_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(MULTITHREADINGLIB) $(GSLLIB) $(PROJ4LIB)

IssmConfig_matlab_la_SOURCES = ../IssmConfig/IssmConfig.cpp
IssmConfig_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
IssmConfig_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB)

ExpSimplify_matlab_la_SOURCES = ../ExpSimplify/ExpSimplify.cpp
ExpSimplify_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
ExpSimplify_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(GSLLIB) $(PROJ4LIB)

Kriging_matlab_la_SOURCES = ../Kriging/Kriging.cpp
Kriging_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
Kriging_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(GSLLIB) $(PROJ4LIB) $(MULTITHREADINGLIB) $(GSLLIB) $(PROJ4LIB)

MeshPartition_matlab_la_SOURCES = ../MeshPartition/MeshPartition.cpp
MeshPartition_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
MeshPartition_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(METISLIB) $(GSLLIB) $(PROJ4LIB)

M1qn3_matlab_la_SOURCES = ../M1qn3/M1qn3.cpp
M1qn3_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
M1qn3_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(METISLIB) $(M1QN3LIB) $(GSLLIB) $(PROJ4LIB)

MeshProfileIntersection_matlab_la_SOURCES = ../MeshProfileIntersection/MeshProfileIntersection.cpp
MeshProfileIntersection_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
MeshProfileIntersection_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(GSLLIB) $(PROJ4LIB)

NodeConnectivity_matlab_la_SOURCES = ../NodeConnectivity/NodeConnectivity.cpp
NodeConnectivity_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
NodeConnectivity_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(GSLLIB) $(PROJ4LIB)

PointCloudFindNeighbors_matlab_la_SOURCES = ../PointCloudFindNeighbors/PointCloudFindNeighbors.cpp
PointCloudFindNeighbors_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
PointCloudFindNeighbors_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(MULTITHREADINGLIB) $(GSLLIB) $(PROJ4LIB)

PropagateFlagsFromConnectivity_matlab_la_SOURCES = ../PropagateFlagsFromConnectivity/PropagateFlagsFromConnectivity.cpp
PropagateFlagsFromConnectivity_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
PropagateFlagsFromConnectivity_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(GSLLIB) $(PROJ4LIB)

Scotch_matlab_la_SOURCES = ../Scotch/Scotch.cpp
Scotch_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
Scotch_matlab_la_LIBADD = ${deps} $(SCOTCHLIB) $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(BLASLAPACKLIB)

ShpRead_matlab_la_SOURCES = ../ShpRead/ShpRead.cpp
ShpRead_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
ShpRead_matlab_la_LIBADD = ${deps} $(SHAPELIBLIB) $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(GSLLIB) $(PROJ4LIB)

TriMesh_matlab_la_SOURCES = ../TriMesh/TriMesh.cpp
TriMesh_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
TriMesh_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(TRIANGLELIB) $(GSLLIB) $(PROJ4LIB)

TriMeshProcessRifts_matlab_la_SOURCES = ../TriMeshProcessRifts/TriMeshProcessRifts.cpp
TriMeshProcessRifts_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
TriMeshProcessRifts_matlab_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(NEOPZLIB) $(GSLLIB) $(PROJ4LIB)
#}}}
