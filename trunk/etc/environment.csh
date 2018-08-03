#ISSM_DIR should have been defined already in your shell settings file (.bashrc, .cshrc, etc ...)

#Load ISSM scripts
setenv PATH {$PATH}:{$ISSM_DIR}/scripts

#MPI
set MPI_DIR="$ISSM_DIR/externalpackages/mpich/install"
if (-d $MPI_DIR) then
	setenv MPI_DIR {$MPI_DIR}
	setenv PATH {$MPI_DIR}/bin:{$PATH}
	setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:$MPI_DIR/lib
	setenv MANPATH {$MANPATH}:$MPI_DIR/man
endif

#PETSC
set PETSC_DIR="$ISSM_DIR/externalpackages/petsc/install"
if (-d $PETSC_DIR) then
	setenv PETSC_DIR {$PETSC_DIR}
	setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:$PETSC_DIR/lib
endif

#SLEPC
set SLEPC_DIR="$ISSM_DIR/externalpackages/slepc/install"
if (-d $SLEPC_DIR) then
	setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:$SLEPC_DIR/lib
endif

#PETSC
set TAO_DIR="$ISSM_DIR/externalpackages/tao/install"
if (-d $TAO_DIR) then
	setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:$TAO_DIR/lib
endif

#Dakota
set DAKOTA_DIR="$ISSM_DIR/externalpackages/dakota/install"
if (-d $DAKOTA_DIR) then
	setenv PATH {$DAKOTA_DIR}/bin:{$PATH}
	setenv MANPATH {$MANPATH}:{$MPI_DIR}/man:{$DAKOTA_DIR}/docs/man:{$DAKOTA_DIR}/docs/man-ref
endif

#Boost
set BOOST_DIR="$ISSM_DIR/externalpackages/boost/install"
set BOOSTROOT="$ISSM_DIR/externalpackages/boost/install"
if (-d $BOOST_DIR) then
   setenv PATH {$BOOST_DIR}/bin:{$PATH}
endif

#Doxygen
set DOXYGEN_DIR="$ISSM_DIR/externalpackages/doxygen/install"
if (-d $DOXYGEN_DIR) then
	setenv MANPATH {$MANPATH}:{$DOXYGEN_DIR}/man
	setenv PATH {$PATH}:{$DOXYGEN_DIR}/bin
endif

#AUTOTOOLS
set AUTOTOOLS_DIR="$ISSM_DIR/externalpackages/autotools/install"
if (-d $AUTOTOOLS_DIR) then
	setenv PATH {$AUTOTOOLS_DIR}/bin:{$PATH}
endif

#SSH
set SSH_DIR="$ISSM_DIR/externalpackages/ssh"
if (-d $SSH_DIR) then
	setenv PATH {$PATH}:{$SSH_DIR}
endif

#VALGRIND
set VALGRIND_DIR="$ISSM_DIR/externalpackages/valgrind/install/bin"
if (-d $VALGRIND_DIR) then
	setenv PATH {$PATH}:{$VALGRIND_DIR}
endif

#MERCURIAL
set MERCURIAL_DIR="$ISSM_DIR/externalpackages/mercurial/install"
if (-d $MERCURIAL_DIR) then
	setenv PYTHONPATH {$MERCURIAL_DIR}/mercurial/pure/
	setenv PATH {$PATH}:{$MERCURIAL_DIR}
endif

#GSL
set GSL_DIR="$ISSM_DIR/externalpackages/gsl/install"
if (-d $GSL_DIR) then
	setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:{$GSL_DIR}/lib
endif

#CMAKE
set CMAKE_DIR="$ISSM_DIR/externalpackages/cmake/install"
if (-d $CMAKE_DIR) then
	setenv PATH {$CMAKE_DIR}/bin:{$PATH}
endif

#YAMS
set YAMS_DIR="$ISSM_DIR/externalpackages/yams/install"
if (-d $YAMS_DIR) then
	setenv PATH {$PATH}:{$YAMS_DIR}/bin
endif

#SHELL2JUNIT
set SHELL2JUNIT_DIR="$ISSM_DIR/externalpackages/shell2junit"
if (-d $SHELL2JUNIT_DIR) then
	setenv PATH {$SHELL2JUNIT_DIR}/install:{$PATH}
endif

#GMT
set GMT_DIR="$ISSM_DIR/externalpackages/gmt"
if (-d $GMT_DIR) then
	setenv PATH {$GMT_DIR}/install/bin/:{$PATH}
endif
