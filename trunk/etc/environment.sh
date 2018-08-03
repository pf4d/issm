#ISSM_DIR and ISSM_ARCH should have been defined already in your shell settings file (.bashrc, .cshrc, etc ...)

pathprepend(){ #{{{
	if [ -d "$1" ] && [[ ":$PATH:" != *":$1:"* ]]; then
		name=$1
		if [[ "$ISSM_ARCH" == "cygwin-intel" ]]; then
			#export path using the cygwin convention
			name=`cygpath -u $1`
		fi
		export PATH="$name:$PATH"
	fi
} #}}}
pathappend(){ #{{{
	if [ -d "$1" ] && [[ ":$PATH:" != *":$1:"* ]]; then
		name=$1
		if [[ "$ISSM_ARCH" == "cygwin-intel" ]]; then
			#export path in cygwin convention
			name=`cygpath -u $1`
		fi
		export PATH="$PATH:$name"
	fi
} #}}}
libpathprepend(){ #{{{
	if [ -d "$1" ] && [[ ":$LD_LIBRARY_PATH:" != *":$1:"* ]]; then
		export LD_LIBRARY_PATH="$1:$LD_LIBRARY_PATH"
	fi
	if [ -d "$1" ] && [[ ":$LD_RUN_PATH:" != *":$1:"* ]]; then
		export LD_RUN_PATH="$1:$LD_RUN_PATH"
	fi
} #}}}
libpathappend(){ #{{{
	if [ -d "$1" ] && [[ ":$LD_LIBRARY_PATH:" != *":$1:"* ]]; then
		export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$1"
	fi
	if [ -d "$1" ] && [[ ":$LD_RUN_PATH:" != *":$1:"* ]]; then
		export LD_RUN_PATH="$LD_RUN_PATH:$1"
	fi
} #}}}
dylibpathprepend(){ #{{{
	if [ -d "$1" ] && [[ ":$DYLD_LIBRARY_PATH:" != *":$1:"* ]]; then
		export DYLD_LIBRARY_PATH="$1:$DYLD_LIBRARY_PATH"
	fi
	if [ -d "$1" ] && [[ ":$LD_RUN_PATH:" != *":$1:"* ]]; then
		export LD_RUN_PATH="$1:$LD_RUN_PATH"
	fi
} #}}}
dylibpathappend(){ #{{{
	if [ -d "$1" ] && [[ ":$DYLD_LIBRARY_PATH:" != *":$1:"* ]]; then
		export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:$1"
	fi
	if [ -d "$1" ] && [[ ":$LD_RUN_PATH:" != *":$1:"* ]]; then
		export LD_RUN_PATH="$LD_RUN_PATH:$1"
	fi
} #}}}

#FIXME: during installation packages are installed one by one but environment.sh was sourced
#before so new packages are NOT in the path.
#may resource environment.sh with:
#if [ -z $(echo "$PATH" | grep "$MATLAB_DIR") ]; then export $PATH...; fi

#Windows compilers: 
if [[ "$ISSM_ARCH" == "cygwin-intel" ]]; then
	source $ISSM_DIR/externalpackages/windows/windows_environment.sh
fi

#Load ISSM scripts
pathappend "$ISSM_DIR/scripts"

#GMT 
GMT_DIR="$ISSM_DIR/externalpackages/gmt/install"
if [ -d "$GMT_DIR" ]; then
	export GMT_DIR
	pathprepend   "$GMT_DIR/bin" 
	libpathappend "$GMT_DIR/lib"
fi

#legacy mpich2 (To be removed)
MPI_DIR="$ISSM_DIR/externalpackages/mpich2/install"
if [ -d "$MPI_DIR" ]; then
	export MPI_DIR
	pathprepend   "$MPI_DIR/bin"
	libpathappend "$MPI_DIR/lib"
fi

MPI_DIR="$ISSM_DIR/externalpackages/mpich/install"
if [ -d "$MPI_DIR" ]; then
	export MPI_DIR
	export MPI_INC_DIR="$MPI_DIR/include"
	pathprepend   "$MPI_DIR/bin"
	libpathappend "$MPI_DIR/lib"
fi

PETSC_DIR="$ISSM_DIR/externalpackages/petsc/install"
if [ -d "$PETSC_DIR" ]; then
	export PETSC_DIR
	libpathappend "$PETSC_DIR/lib"
fi

SLEPC_DIR="$ISSM_DIR/externalpackages/slepc/install"
libpathappend "$SLEPC_DIR/lib"

TAO_DIR="$ISSM_DIR/externalpackages/tao/install"
libpathappend "$TAO_DIR/lib"

DAKOTA_DIR="$ISSM_DIR/externalpackages/dakota/install"
pathappend "$DAKOTA_DIR/bin"
libpathappend "$DAKOTA_DIR/lib"
dylibpathprepend "$DAKOTA_DIR/lib"

DOXYGEN_DIR="$ISSM_DIR/externalpackages/doxygen/install"
pathprepend "$DOXYGEN_DIR/bin"

AUTOTOOLS_DIR="$ISSM_DIR/externalpackages/autotools/install"
pathprepend "$AUTOTOOLS_DIR/bin"

SDK_DIR="C:/MicrosoftVisualStudio 9.0/Microsoft Visual C++ 2008 Express Edition with SP1 - ENU"
pathappend "$SDK_DIR"

SSH_DIR="$ISSM_DIR/externalpackages/ssh"
pathappend "$SSH_DIR"

VALGRIND_DIR="$ISSM_DIR/externalpackages/valgrind/install"
pathprepend "$VALGRIND_DIR/bin"

CPPCHECK_DIR="$ISSM_DIR/externalpackages/cppcheck/install"
pathappend "$CPPCHECK_DIR/bin"

GDAL_DIR="$ISSM_DIR/externalpackages/gdal/install"
pathprepend "$GDAL_DIR/bin"
libpathappend "$GDAL_DIR/lib"

PROJ4_DIR="$ISSM_DIR/externalpackages/proj.4/install"
dylibpathprepend "$PROJ4_DIR/lib"
libpathprepend "$PROJ4_DIR/lib"

MERCURIAL_DIR="$ISSM_DIR/externalpackages/mercurial/install"
if [ -d "$MERCURIAL_DIR" ]; then
	export PYTHONPATH="$PYTHONPATH:$MERCURIAL_DIR/mercurial/pure/"
	pathappend "$MERCURIAL_DIR"
fi

BOOST_DIR="$ISSM_DIR/externalpackages/boost/install"
BOOSTROOT="$ISSM_DIR/externalpackages/boost/install"
if [ -d "$BOOST_DIR" ]; then
	export BOOSTROOT
	export BOOST_DIR
	libpathprepend   "$BOOST_DIR/lib"
	dylibpathprepend "$BOOST_DIR/lib"
	pathprepend      "$BOOST_DIR/bin"
fi

XERCESROOT="$ISSM_DIR/externalpackages/xerces/install"
XERCESCROOT="$ISSM_DIR/externalpackages/xerces/src"
if [ -d "$XERCESROOT" ]; then
	export XERCESROOT 
	export XERCESCROOT
fi

XAIF_DIR="$ISSM_DIR/externalpackages/xaifbooster/xaifBooster"
XAIFBOOSTERROOT="$ISSM_DIR/externalpackages/xaifbooster/"
XAIFBOOSTER_HOME="$ISSM_DIR/externalpackages/xaifbooster/xaifBooster"
PLATFORM="x86-Linux"
if [ -d "$XAIF_DIR" ]; then
	export XAIFBOOSTERROOT
	export XAIFBOOSTER_HOME
	export XAIF_DIR
	export PLATFORM
fi

ANGELROOT="$ISSM_DIR/externalpackages/angel/angel"
if [ -d "$ANGELROOT" ]; then
	export ANGELROOT
fi

OPENANALYSISROOT="$ISSM_DIR/externalpackages/openanalysis/install"
if [ -d "$OPENANALYSISROOT" ]; then
	export OPENANALYSISROOT
	libpathappend "$OPENANALYSISROOT/lib"
fi

JVM_DIR="/usr/local/gcc/4.3.2/lib64/gcj-4.3.2-9/"
libpathappend "$JVM_DIR"

BBFTP_DIR="$ISSM_DIR/externalpackages/bbftp/install"
pathappend "$BBFTP_DIR/bin"

ADIC_DIR="$ISSM_DIR/externalpackages/adic/install"
pathappend "$ADIC_DIR/bin"
libpathappend "$ADIC_DIR/lib"

COLPACK_DIR="$ISSM_DIR/externalpackages/colpack/install"
libpathappend "$COLPACK_DIR/lib"

ECLIPSE_DIR="$ISSM_DIR/externalpackages/eclipse/install"
pathappend "$ECLIPSE_DIR"

APPSCAN_DIR="$ISSM_DIR/externalpackages/appscan/install"
pathappend "$APPSCAN_DIR/bin"

RATS_DIR="$ISSM_DIR/externalpackages/rats/install"
pathappend "$RATS_DIR/bin"

DYSON_DIR="$ISSM_DIR/externalpackages/dyson/"
pathappend "$DYSON_DIR"

CMAKE_DIR="$ISSM_DIR/externalpackages/cmake/install"
pathprepend "$CMAKE_DIR/bin"

SHAPELIB_DIR="$ISSM_DIR/externalpackages/shapelib/install"
pathappend "$SHAPELIB_DIR/exec"

CCCL_DIR="$ISSM_DIR/externalpackages/cccl/install"
pathappend "$CCCL_DIR/bin"

PACKAGEMAKER_DIR="$ISSM_DIR/externalpackages/packagemaker/install"
pathappend "$PACKAGEMAKER_DIR"

#android-dev-dir
export ANDROID_DIR="$ISSM_DIR/externalpackages/android"

export ANDROID_NDK_DIR="$ANDROID_DIR/android-ndk/install"
pathappend "$ANDROID_NDK_DIR/arm-linux-android-install/bin"

export ANDROID_SDK_DIR="$ANDROID_DIR/android-sdk/install"
pathappend "$ANDROID_SDK_DIR/"

GSL_DIR="$ISSM_DIR/externalpackages/gsl/install"
libpathappend "$GSL_DIR/lib"

GMAKE_DIR="$ISSM_DIR/externalpackages/gmake/install"
pathprepend "$GMAKE_DIR/bin"

PYTHON_DIR="$ISSM_DIR/externalpackages/python/install"
if [ -d "$PYTHON_DIR" ]; then
	export PYTHONPATH="$PYTHONPATH:$ISSM_DIR/lib"
	pathprepend    "$PYTHON_DIR/bin"
	libpathprepend "$PYTHON_DIR/lib"
	libpathprepend "$ISSM_DIR/lib"
fi

MODELE_DIR="$ISSM_DIR/externalpackages/modelE/install"
pathappend "$MODELE_DIR/src/exec"

GIT_DIR="$ISSM_DIR/externalpackages/git/install"
pathprepend "$GIT_DIR/bin"

NCVIEW_DIR="$ISSM_DIR/externalpackages/ncview/install"
pathappend "$NCVIEW_DIR"

TCLX_DIR="$ISSM_DIR/externalpackages/tclx/install/lib/tclx8.4"
libpathappend "$TCLX_DIR"

ASPELL_DIR="$ISSM_DIR/externalpackages/aspell/install"
pathappend "$ASPELL_DIR/bin"

HDF5_DIR="$ISSM_DIR/externalpackages/hdf5/install"
dylibpathappend "$HDF5_DIR/lib"
libpathappend "$HDF5_DIR/lib"
if [ -d "$HDF5_DIR" ]; then
	export LIBRARY_PATH="$LIBRARY_PATH:$HDF5_DIR/lib"
	export C_INCLUDE_PATH="$C_INCLUDE_PATH:$HDF5_DIR/include"
fi

NETCDF_DIR="$ISSM_DIR/externalpackages/netcdf/install"
pathappend "$NETCDF_DIR/bin"
dylibpathappend "$NETCDF_DIR/lib"
libpathappend "$NETCDF_DIR/lib"
if [ -d "$NETCDF_DIR" ]; then
	export LIBRARY_PATH="$LIBRARY_PATH:$NETCDF_DIR/lib"
	dylibpathappend "$NETCDF_DIR/lib"
	libpathappend "$NETCDF_DIR/lib"
	export C_INCLUDE_PATH="$C_INCLUDE_PATH:$NETCDF_DIR/include"
fi

NETCDF_CXX_DIR="$ISSM_DIR/externalpackages/netcdf-cxx/install"
libpathappend "$NETCDF_CXX_DIR/lib"

SVN_DIR="$ISSM_DIR/externalpackages/svn/install"
pathprepend   "$SVN_DIR/bin"
libpathappend "$SVN_DIR/lib"

CVS_DIR="$ISSM_DIR/externalpackages/cvs/install"
pathprepend   "$CVS_DIR/bin"

APR_DIR="$ISSM_DIR/externalpackages/apr/install"
pathappend    "$APR_DIR/bin"
libpathappend "$APR_DIR/lib"

APR_UTIL_DIR="$ISSM_DIR/externalpackages/apr-util/install"
pathappend   "$APR_UTIL_DIR/bin:$PATH"
libpathappend "$APR_UTIL_DIR/lib"

SQLITE_DIR="$ISSM_DIR/externalpackages/sqlite/install"
pathappend   "$SQLITE_DIR/bin"
libpathappend "$SQLITE_DIR/lib"

YAMS_DIR="$ISSM_DIR/externalpackages/yams/install"
pathappend   "$YAMS_DIR"

SWIG_DIR="$ISSM_DIR/externalpackages/swig/install"
pathappend   "$SWIG_DIR"

#AUX-CONFIG
pathappend   "$ISSM_DIR/aux-config"

#INISHELL
pathappend   "$ISSM_DIR/externalpackages/inishell/install"

#SHELL2JUNIT
pathappend   "$ISSM_DIR/externalpackages/shell2junit/install"

#EXPAT
libpathprepend   "$ISSM_DIR/externalpackages/expat/install"
dylibpathprepend   "$ISSM_DIR/externalpackages/expat/install"

#GMSH
pathappend   "$ISSM_DIR/externalpackages/gmsh/install"

#CURL
libpathprepend   "$ISSM_DIR/externalpackages/curl/install/lib"
dylibpathprepend   "$ISSM_DIR/externalpackages/curl/install/lib"
pathprepend "$ISSM_DIR/externalpackages/curl/install/bin"

#GMT
pathprepend "$ISSM_DIR/externalpackages/gmt/install/bin"

#NEOPZ
NEOPZ_DIR="$ISSM_DIR/externalpackages/neopz/install"
if [ -d "$NEOPZ_DIR" ]; then
	export REFPATTERNDIR="$NEOPZ_DIR/include/refpatterns"
fi
