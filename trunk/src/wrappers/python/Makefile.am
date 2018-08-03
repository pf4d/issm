AM_CPPFLAGS = @DAKOTAINCL@ @PETSCINCL@ @MPIINCL@ @SPOOLESINCL@ @METISINCL@ @TRIANGLEINCL@ @CHACOINCL@ @SCOTCHINCL@ @SHAPELIBINCL@ @PYTHONINCL@ @PYTHON_NUMPYINCL@ @AMPIINCL@
AUTOMAKE_OPTIONS = subdir-objects

EXEEXT=$(PYTHONWRAPPEREXT)

#define prefix (from http://www.gnu.org/software/autoconf/manual/autoconf-2.67/html_node/Defining-Directories.html)
AM_CPPFLAGS+=  -DISSM_PREFIX='"$(prefix)"'

#python io{{{
lib_LTLIBRARIES = libISSMPython.la

io_sources= ./io/WritePythonData.cpp\
				./io/CheckNumPythonArguments.cpp\
				./io/FetchPythonData.cpp

ALLCXXFLAGS= -fPIC -D_WRAPPERS_ $(CXXFLAGS) $(CXXOPTFLAGS) 
libISSMPython_la_SOURCES = $(io_sources)
libISSMPython_la_CXXFLAGS= $(ALLCXXFLAGS)
#}}}
#api io{{{
lib_LTLIBRARIES += libISSMApi.la

api_sources= ./io/ApiPrintf.cpp

libISSMApi_la_SOURCES = $(api_sources)
libISSMApi_la_CXXFLAGS= $(ALLCXXFLAGS)
#}}}
#Wrappers {{{
if WRAPPERS
lib_LTLIBRARIES += BamgConvertMesh_python.la\
						BamgMesher_python.la\
						BamgTriangulate_python.la\
						ContourToMesh_python.la\
						ContourToNodes_python.la\
						ElementConnectivity_python.la\
						InterpFromMeshToMesh2d_python.la\
						InterpFromMeshToMesh3d_python.la\
						InterpFromGridToMesh_python.la\
						InterpFromMeshToGrid_python.la\
						IssmConfig_python.la\
						MeshProfileIntersection_python.la\
						NodeConnectivity_python.la\
						TriMesh_python.la\
						TriMeshProcessRifts_python.la
endif 
#}}}
#Flags and libraries {{{
deps =  $(MATHLIB) ${PYTHONLIB}

#Triangle library
AM_CXXFLAGS =  -DTRILIBRARY -DANSI_DECLARATORS -DNO_TIMER -D_WRAPPERS_

#Python part
AM_LDFLAGS   = $(PYTHONLINK) -shrext ${EXEEXT} -module
if VERSION
AM_LDFLAGS +=
else
AM_LDFLAGS += -avoid-version
endif

if STANDALONE_MODULES
AM_LDFLAGS += -Wl,-static
endif

deps += ./libISSMPython.la 
AM_LDFLAGS += --no-warnings 

AM_CXXFLAGS +=  -D_HAVE_PYTHON_MODULES_   -fPIC
if PYTHON3
AM_CXXFLAGS +=  -DNPY_NO_DEPRECATED_API 
endif

deps += ../../c/libISSMModules.la ../../c/libISSMCore.la
if ADOLC
deps += $(ADOLCLIB)
endif

deps += ./libISSMApi.la 

if STANDALONE_LIBRARIES
libISSMPython_la_LDFLAGS = -static 
endif
if STANDALONE_LIBRARIES
libISSMApi_la_LDFLAGS = -static 
endif

#Optimization flags:
AM_CXXFLAGS += $(CXXOPTFLAGS) 
#}}}
#Bin sources {{{
BamgConvertMesh_python_la_SOURCES = ../BamgConvertMesh/BamgConvertMesh.cpp
BamgConvertMesh_python_la_CXXFLAGS = ${AM_CXXFLAGS}
BamgConvertMesh_python_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(GSLLIB) $(PROJ4LIB)

BamgMesher_python_la_SOURCES = ../BamgMesher/BamgMesher.cpp
BamgMesher_python_la_CXXFLAGS = ${AM_CXXFLAGS}
BamgMesher_python_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(GSLLIB) $(PROJ4LIB)

BamgTriangulate_python_la_SOURCES = ../BamgTriangulate/BamgTriangulate.cpp
BamgTriangulate_python_la_CXXFLAGS = ${AM_CXXFLAGS}
BamgTriangulate_python_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(GSLLIB) $(PROJ4LIB)

ContourToMesh_python_la_SOURCES = ../ContourToMesh/ContourToMesh.cpp
ContourToMesh_python_la_CXXFLAGS = ${AM_CXXFLAGS}
ContourToMesh_python_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(GSLLIB) $(PROJ4LIB)

ContourToNodes_python_la_SOURCES = ../ContourToNodes/ContourToNodes.cpp
ContourToNodes_python_la_CXXFLAGS = ${AM_CXXFLAGS}
ContourToNodes_python_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(GSLLIB) $(PROJ4LIB)

ElementConnectivity_python_la_SOURCES = ../ElementConnectivity/ElementConnectivity.cpp
ElementConnectivity_python_la_CXXFLAGS = ${AM_CXXFLAGS}
ElementConnectivity_python_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(GSLLIB) $(PROJ4LIB)

InterpFromMeshToMesh2d_python_la_SOURCES = ../InterpFromMeshToMesh2d/InterpFromMeshToMesh2d.cpp
InterpFromMeshToMesh2d_python_la_CXXFLAGS = ${AM_CXXFLAGS}
InterpFromMeshToMesh2d_python_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(MULTITHREADINGLIB) $(GSLLIB) $(PROJ4LIB)

InterpFromMeshToMesh3d_python_la_SOURCES = ../InterpFromMeshToMesh3d/InterpFromMeshToMesh3d.cpp
InterpFromMeshToMesh3d_python_la_CXXFLAGS = ${AM_CXXFLAGS}
InterpFromMeshToMesh3d_python_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(MULTITHREADINGLIB) $(GSLLIB) $(PROJ4LIB)

InterpFromGridToMesh_python_la_SOURCES = ../InterpFromGridToMesh/InterpFromGridToMesh.cpp
InterpFromGridToMesh_python_la_CXXFLAGS = ${AM_CXXFLAGS}
InterpFromGridToMesh_python_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(MULTITHREADINGLIB) $(GSLLIB) $(PROJ4LIB)

InterpFromMeshToGrid_python_la_SOURCES = ../InterpFromMeshToGrid/InterpFromMeshToGrid.cpp
InterpFromMeshToGrid_python_la_CXXFLAGS = ${AM_CXXFLAGS}
InterpFromMeshToGrid_python_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(MULTITHREADINGLIB) $(GSLLIB) $(PROJ4LIB)

IssmConfig_python_la_SOURCES = ../IssmConfig/IssmConfig.cpp
IssmConfig_python_la_CXXFLAGS = ${AM_CXXFLAGS}
IssmConfig_python_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB)

MeshProfileIntersection_python_la_SOURCES = ../MeshProfileIntersection/MeshProfileIntersection.cpp
MeshProfileIntersection_python_la_CXXFLAGS = ${AM_CXXFLAGS}
MeshProfileIntersection_python_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(GSLLIB) $(PROJ4LIB)

NodeConnectivity_python_la_SOURCES = ../NodeConnectivity/NodeConnectivity.cpp
NodeConnectivity_python_la_CXXFLAGS = ${AM_CXXFLAGS}
NodeConnectivity_python_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(GSLLIB) $(PROJ4LIB)

TriMesh_python_la_SOURCES = ../TriMesh/TriMesh.cpp
TriMesh_python_la_CXXFLAGS = ${AM_CXXFLAGS}
TriMesh_python_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(TRIANGLELIB) $(GSLLIB) $(PROJ4LIB)

TriMeshProcessRifts_python_la_SOURCES = ../TriMeshProcessRifts/TriMeshProcessRifts.cpp
TriMeshProcessRifts_python_la_CXXFLAGS = ${AM_CXXFLAGS}
TriMeshProcessRifts_python_la_LIBADD = ${deps} $(MPILIB) $(PETSCLIB) $(GSLLIB) $(PROJ4LIB)
#}}}
