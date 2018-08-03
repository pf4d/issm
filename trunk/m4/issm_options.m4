dnl ISSM Options

AC_DEFUN([ISSM_OPTIONS],[

	AC_MSG_NOTICE(============================================================================)
	AC_MSG_NOTICE(=                      Checking ISSM specific options                      =)
	AC_MSG_NOTICE(============================================================================)

	dnl ISSM's internal options
	dnl Build info{{{
	
	dnl build date
	AC_PATH_PROGS(DATE, date)
	AC_MSG_CHECKING(for build date)
	if test "$DATE" ; then
		PACKAGE_DATE=`date`
	else
		PACKAGE_DATE="unknown"
	fi
	AC_DEFINE_UNQUOTED(PACKAGE_BUILD_DATE,"$PACKAGE_DATE", Build date)
	AC_MSG_RESULT($PACKAGE_DATE)

	dnl user name
	AC_MSG_CHECKING([user name])
	if test -n "$USER"
	then
		user_name="$USER"
	else
		if test -n "$LOGNAME"
		then
			user_name="$LOGNAME"
		else
		   user_name=`(whoami) 2>/dev/null` || user_name=unknown
		fi
	fi
	AC_DEFINE_UNQUOTED(USER_NAME, "$user_name", Build user name)
	AC_MSG_RESULT($user_name)

	AC_MSG_CHECKING([host full OS name and version])
	dnl normalize some host OS names
	case ${host_os} in
	  dnl linux is linux is linux, regardless of RMS.
	  linux-gnu* | lignux* )	host_os=linux ;;
	esac
	AC_DEFINE_UNQUOTED(HOST_OS, "$host_os", Host full OS name and version)
	AC_MSG_RESULT($host_os)

  AC_MSG_CHECKING([host cpu])
  AC_DEFINE_UNQUOTED(HOST_CPU, "$host_cpu",Host cpu)
  AC_MSG_RESULT($host_cpu)

  AC_MSG_CHECKING([vendor])
  AC_DEFINE_UNQUOTED(HOST_VENDOR, "$host_vendor",Host vendor)
  AC_MSG_RESULT($host_vendor)

  AC_MSG_CHECKING([host OS name])
  host_os_name=`echo $host_os | sed 's/\..*//g'`
  dnl normalize some OS names
  case ${host_os_name} in
	dnl linux is linux is linux, regardless of RMS.
	linux-gnu* | lignux* )	host_os_name=linux ;;
  esac
  AC_DEFINE_UNQUOTED(HOST_OS_NAME, "$host_os_name", Host OS name)
  AC_MSG_RESULT($host_os_name)

	dnl parse out the OS version of the host
  AC_MSG_CHECKING([host OS version])
  host_os_version=`echo $host_os | sed 's/^[[^0-9]]*//g'`
  if test -z "$host_os_version"
  then
	host_os_version=`(uname -r) 2>/dev/null` || host_os_version=unknown
  fi
  AC_DEFINE_UNQUOTED(HOST_OS_VERSION, "$host_os_version", Host OS version)
  AC_MSG_RESULT($host_os_version)


	dnl figure out host architecture (different than CPU)
  AC_MSG_CHECKING([host OS architecture])
  host_arch=`(uname -m) 2>/dev/null` || host_arch=unknown
  dnl normalize some names
  case ${host_arch} in
	sun4* )	host_arch=sun4 ;;
	sun3x )	host_arch=sun3 ;;
	sun )	host_arch=`(arch) 2>/dev/null` || host_arch=unknown ;;
	i?86 )	host_arch=i386 ;; # all x86 should show up as i386
  esac
  AC_DEFINE_UNQUOTED(HOST_ARCH, "$host_arch",Host Archictecture)
  AC_MSG_RESULT($host_arch)

	dnl }}}
	dnl Debugging {{{
	AC_ARG_ENABLE([debugging],                                        dnl feature
		AS_HELP_STRING([--enable-debugging],[turn debug support on]),  dnl help string
		[enable_debugging=$enableval],                                 dnl action if given
		[enable_debugging=no])                                         dnl action if not given

	AC_MSG_CHECKING(for debugging support)
	if test "x$enable_debugging" = xyes; then
		AC_DEFINE([_ISSM_DEBUG_],[1],[Macro to enable debugging in ISSM])
	fi
	AC_MSG_RESULT($enable_debugging)
	dnl }}}
	dnl Development{{{
	AC_ARG_ENABLE([development],                                      dnl feature
		AS_HELP_STRING([--enable-development],[turn development on]),  dnl help string
		[enable_development=$enableval],                                 dnl action if given
		[enable_development=no])                                      dnl action if not given

	AC_MSG_CHECKING(for development support)
	if test "x$enable_development" = xyes; then
		AC_DEFINE([_DEVELOPMENT_],[1],[Macro to enable development version in ISSM])
	fi
	AM_CONDITIONAL([DEVELOPMENT], [test x$enable_development = xyes])
	AC_MSG_RESULT($enable_development)
	 dnl }}}
    dnl Standalone Options {{{
    AC_ARG_ENABLE([standalone-modules],                                                      dnl feature
        AS_HELP_STRING([--enable-standalone-modules], [produce standalone modules]),         dnl help string
        [enable_standalone_modules=$enableval],                                              dnl action if given
        [enable_standalone_modules=no])                                                      dnl action if not given
	AC_MSG_CHECKING(for standalone modules build)
    AM_CONDITIONAL([STANDALONE_MODULES], [test x$enable_standalone_modules = xyes])
	AC_MSG_RESULT($enable_standalone_modules)

    AC_ARG_ENABLE([standalone-executables],                                                  dnl feature
        AS_HELP_STRING([--enable-standalone-executables], [produce standalone executables]), dnl help string
        [enable_standalone_executables=$enableval],                                          dnl action if given
        [enable_standalone_executables=no])                                                  dnl action if not given
	AC_MSG_CHECKING(for standalone executables build)
    AM_CONDITIONAL([STANDALONE_EXECUTABLES], [test x$enable_standalone_executables = xyes])
	AC_MSG_RESULT($enable_standalone_executables)

    AC_ARG_ENABLE([standalone-libraries],                                                    dnl feature
        AS_HELP_STRING([--enable-standalone-libraries], [produce standalone libraries]),     dnl help string
        [enable_standalone_libraries=$enableval],                                            dnl action if given
        [enable_standalone_libraries=no])                                                    dnl action if not given
	AC_MSG_CHECKING(for standalone libraries build)
    AM_CONDITIONAL([STANDALONE_LIBRARIES], [test x$enable_standalone_libraries = xyes])
	AC_MSG_RESULT($enable_standalone_libraries)
    dnl }}}
    dnl Version{{{
    AC_ARG_ENABLE([version],                                   dnl feature
    AS_HELP_STRING([--enable-version],[produce libISSM.so.0]), dnl help string
    [enable_version=$enableval],                               dnl action if given
    [enable_version=no])                                       dnl action if not given
    AM_CONDITIONAL([VERSION], [test x$enable_VERSION = xyes])
    dnl }}}
	dnl Wrappers build {{{
	AC_ARG_WITH([wrappers],                                           dnl feature
	AS_HELP_STRING([--with-wrappers = value],[wrappers compilation]), dnl help string
	[WRAPPERS_VALUE=$withval],                                        dnl action if given
	[WRAPPERS_VALUE="yes"])                                           dnl action if not given
	AC_MSG_CHECKING(for wrappers compilation)
	AM_CONDITIONAL([WRAPPERS], [test x$WRAPPERS_VALUE = xyes])
	AC_MSG_RESULT($WRAPPERS_VALUE) 
	dnl }}}
	dnl Extensions{{{
	ISSMEXT=".exe"
	AC_SUBST([ISSMEXT])
	dnl }}}

	dnl ISSM's externalpackages
	dnl vendor{{{
	AC_ARG_WITH([vendor],                                              dnl feature
	AS_HELP_STRING([--with-vendor = VENDOR],[vendor name, ex: intel]), dnl help string
	[VENDOR=$withval],                                                 dnl action if given
	[VENDOR=""])                                                       dnl action if not given
	AC_MSG_CHECKING(for vendor compilers)
	if test -n "$VENDOR"; then
		if  test $VENDOR = intel-win32; then
			export CC=icl
			export CXX=icl
			export CFLAGS="-DWIN32 -D_INTEL_WIN_"
			export CXXFLAGS="-DWIN32 -D_INTEL_WIN_"
			IS_WINDOWS=yes
		elif  test $VENDOR = intel-win7-32; then
			export CC=cl
			export CXX=cl
			export CXXFLAGS="-DWIN32 -D_INTEL_WIN_ -EHsc"
			export CFLAGS="-DWIN32 -D_INTEL_WIN_ -EHsc"
			export AR="ar-lib lib"
			export OS_LDFLAG="-Wl,"
			export RANLIB=true
			IS_WINDOWS=yes
			OSLIBS="-Wl,kernel32.lib -Wl,user32.lib -Wl,gdi32.lib -Wl,winspool.lib -Wl,comdlg32.lib -Wl,advapi32.lib -Wl,shell32.lib -Wl,ole32.lib -Wl,oleaut32.lib -Wl,uuid.lib -Wl,odbc32.lib -Wl,odbccp32.lib"
		elif  test $VENDOR = intel-win7-64; then
			export CC=cl
			export CXX=cl
			export CXXFLAGS="-DWIN32 -D_INTEL_WIN_ -EHsc"
			export CFLAGS="-DWIN32 -D_INTEL_WIN_ -EHsc"
			export AR="ar-lib lib"
			export OS_LDFLAG="-Wl,"
			export RANLIB=true
			IS_WINDOWS=yes
			OSLIBS="-Wl,kernel32.lib -Wl,user32.lib -Wl,gdi32.lib -Wl,winspool.lib -Wl,comdlg32.lib -Wl,advapi32.lib -Wl,shell32.lib -Wl,ole32.lib -Wl,oleaut32.lib -Wl,uuid.lib -Wl,odbc32.lib -Wl,odbccp32.lib"
		elif  test $VENDOR = MSVC-Win64; then
			export CC=cl
			export CXX=cl
			export CXXFLAGS="-DWIN32 -D_INTEL_WIN_ -D_HAVE_PETSC_MPI_ -EHsc"
			export CFLAGS="-DWIN32 -D_INTEL_WIN_ -D_HAVE_PETSC_MPI_ -EHsc"
			export AR="ar-lib lib"
			export OS_LDFLAG="-Wl,"
			export RANLIB=true
			IS_WINDOWS=yes
			OSLIBS="-Wl,kernel32.lib -Wl,user32.lib -Wl,gdi32.lib -Wl,winspool.lib -Wl,comdlg32.lib -Wl,advapi32.lib -Wl,shell32.lib -Wl,ole32.lib -Wl,oleaut32.lib -Wl,uuid.lib -Wl,odbc32.lib -Wl,odbccp32.lib"
		elif  test $VENDOR = MSVC-Win64-par; then
			export CC=cl
			export CXX=cl
			export CXXFLAGS="-DWIN32 -D_INTEL_WIN_ -EHsc"
			export CFLAGS="-DWIN32 -D_INTEL_WIN_ -EHsc"
			export AR="ar-lib lib"
			export OS_LDFLAG="-Wl,"
			export RANLIB=true
			IS_WINDOWS=yes
			OSLIBS="-Wl,kernel32.lib -Wl,user32.lib -Wl,gdi32.lib -Wl,winspool.lib -Wl,comdlg32.lib -Wl,advapi32.lib -Wl,shell32.lib -Wl,ole32.lib -Wl,oleaut32.lib -Wl,uuid.lib -Wl,odbc32.lib -Wl,odbccp32.lib"
		elif test $VENDOR = intel-linux; then
			export CC=icc
			export CXX=icpc
			export CFLAGS=" -D_INTEL_LINUX_"
			export CXXFLAGS=" -D_INTEL_LINUX_"
		elif test $VENDOR = intel-gp; then
			export CC=icc
			export CXX=icpc
			export CFLAGS=" -D_INTEL_LINUX_"
			export CXXFLAGS=" -D_INTEL_LINUX_"
		elif test $VENDOR = intel-lonestar; then
			export CC=icc
			export CXX=icpc
		elif test $VENDOR = intel-discover; then
			export CC=icc
			export CXX=icpc
			export CXXFLAGS=" -O3 -D_INTEL_LINUX_ -DMPICH_IGNORE_CXX_SEEK"
			export CFLAGS=" -O3 -D_INTEL_LINUX_ -DMPICH_IGNORE_CXX_SEEK"
		elif test $VENDOR = intel-pleiades; then
			export CC=icc
			export CXX=icpc
			export CXXFLAGS=" -O3 -D_INTEL_LINUX_ "
			export CFLAGS=" -O3 -D_INTEL_LINUX_ "
		elif test $VENDOR = intel-acenet; then
			export CC=icc
			export CXX=icpc
			export CXXFLAGS=" -D_INTEL_LINUX_ "
			export CFLAGS=" -D_INTEL_LINUX_ "
		elif test $VENDOR = intel-pleiades-gcc; then
			export CC=gcc
			export CXX=g++
			export CXXFLAGS="-O3 -march=corei7-avx"
			export CFLAGS="-O3 -march=corei7-avx"
        else
		AC_MSG_ERROR([unknown compiler vendor!])
		fi
	fi
	AC_SUBST([OSLIBS]) 
	AC_MSG_RESULT(done)
	dnl }}}
	dnl matlab{{{

	dnl 1. See if matlab has been provided
	AC_ARG_WITH([matlab-dir],                                         dnl feature
	AS_HELP_STRING([--with-matlab-dir=DIR],[matlab root directory.]), dnl help string
	[MATLAB_ROOT=$withval],                                           dnl action if given
	[MATLAB_ROOT="no"])                                               dnl action if not given

	AC_MSG_CHECKING([whether matlab is enabled])
	if test "x$MATLAB_ROOT" = "xno" ; then
		 HAVE_MATLAB=no
	else
		HAVE_MATLAB=yes
		if ! test -d "$MATLAB_ROOT"; then
		  AC_MSG_ERROR([matlab directory provided ($MATLAB_ROOT) does not exist]);
		fi
		if ! test -f "$MATLAB_ROOT/extern/include/mex.h"; then
			AC_MSG_ERROR([Couldn't find mex.h... check your installation of matlab])
	   fi
	fi
	AC_MSG_RESULT($HAVE_MATLAB)
	AM_CONDITIONAL([MATLAB], [test x$HAVE_MATLAB = xyes])

	dnl 2. Get Matlab libraries
	if test "x$HAVE_MATLAB" = "xyes"; then

		AC_DEFINE([_HAVE_MATLAB_],[1],[with matlab in ISSM src])
		
		dnl 4. get MEXLIB MEXLINK and MEXEXT (experimental) except for windows
		AC_MSG_CHECKING([matlab's mex compilation flags])
  		case "${host_os}" in
  			*cygwin*) 
  				if  test $VENDOR = intel-win7-32; then
  					MEXLIB="-Wl,libmx.lib -Wl,libmex.lib -Wl,libmat.lib ${OSLIBS} -Wl,libf2cblas.lib -Wl,libf2clapack.lib" 
               MEXLINK="-Wl,/LIBPATH:`cygpath -m ${MATLAB_ROOT}/extern/lib/win32/microsoft` -Wl,/link -Wl,/EXPORT:mexFunction -Wl,/DLL"
					MEXEXT=`$MATLAB_ROOT/bin/mexext.bat`
					MEXEXT=".$MEXEXT"
  				elif test $VENDOR = intel-win7-64; then
  					MEXLIB="-Wl,libmx.lib -Wl,libmex.lib -Wl,libmat.lib ${OSLIBS} -Wl,libf2cblas.lib -Wl,libf2clapack.lib" 
               MEXLINK="-Wl,/LIBPATH:`cygpath -m ${MATLAB_ROOT}/extern/lib/win64/microsoft` -Wl,/link -Wl,/EXPORT:mexFunction -Wl,/DLL" 
					MEXEXT=".mexw64"
  				elif test $VENDOR = MSVC-Win64 || test $VENDOR = MSVC-Win64-par; then
  					MEXLIB="-Wl,libmx.lib -Wl,libmex.lib -Wl,libmat.lib ${OSLIBS} -Wl,libf2cblas.lib -Wl,libf2clapack.lib" 
               MEXLINK="-Wl,/link -Wl,/LIBPATH:`cygpath -m ${MATLAB_ROOT}/extern/lib/win64/microsoft` -Wl,/link -Wl,/EXPORT:mexFunction -Wl,/DLL" 
  					MATLABINCL="-I`cygpath -m $MATLAB_ROOT/extern/include/`"
					MEXEXT=".mexw64"
  				fi
  			;;
		   *)
           MATLABINCL="-I$MATLAB_ROOT/extern/include/"
           MEXLINK=$($MATLAB_ROOT/bin/mex -v 2>&1 < /dev/null | grep LDFLAGS     | sed -e "s/         LDFLAGS            = //g")
			  MEXLIB=$( $MATLAB_ROOT/bin/mex -v 2>&1 < /dev/null | grep CXXLIBS     | sed -e "s/         CXXLIBS            = //g")
		     MEXEXT=$( $MATLAB_ROOT/bin/mex -v 2>&1 < /dev/null | grep LDEXTENSION | sed -e "s/         LDEXTENSION        = //g")
				dnl version 2014 and up
				if test "x$MEXEXT" = "x" ; then
					 echo "#include <mex.h>" > conftest.cpp
					 echo "void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){}" >> conftest.cpp
					 $MATLAB_ROOT/bin/mex -v -lmex conftest.cpp > conftest.tmp 2>&1 
					 rm -f conftest.cpp
					 MEXLINK=$(cat conftest.tmp | grep LDFLAGS  | sed -e "s/LDFLAGS ://g")
					 MEXLIB=$( cat conftest.tmp | grep LINKLIBS | sed -e "s/LINKLIBS ://g")
					 MEXEXT=$( cat conftest.tmp | grep LDEXT    | sed -e "s/LDEXT ://g" | awk '{print $[1]}')
					 rm -f conftest.tmp
				fi

				dnl Make sure mexFunction.map is not in MEXLIB to avoid problems with global variables
				dnl MEXLINK=$(echo $MEXLINK | sed -e "s/,-expo.*mexFunction\\.map\"//g" | sed -e "s/-[[^ ]]*mexFunction\\.map//g")
				MEXLINK="" dnl We actually don't need MEXLINK????

  			;;
      esac
		AC_MSG_RESULT(done)
	   if test "x$MEXEXT" = "x" ; then
			AC_MSG_ERROR([Couldn't find mex... check your installation of matlab])
	   fi

		AC_MSG_CHECKING([whether matlab is enabled])

		AC_SUBST([MATLABINCL])
		MATLABWRAPPEREXT=$MEXEXT
		AC_SUBST([MATLABWRAPPEREXT])
	   AC_SUBST([MEXLIB]) 
		AC_SUBST([MEXLINK])
	fi
	dnl }}}
	dnl windows {{{
	AC_MSG_CHECKING([Checking if this is a Win build... ])
	AM_CONDITIONAL([WINDOWS], [test x$IS_WINDOWS = xyes])
	AC_MSG_RESULT(done)
	dnl }}}
	dnl javascript{{{
	AC_ARG_WITH([javascript],
	  AS_HELP_STRING([--with-javascript], [compile javascript wrappers? default is no.]),
	  [JAVASCRIPT=$withval],[JAVASCRIPT="no"]) 

	dnl Check whether javascript wrappers are desired
	AC_MSG_CHECKING([for javascript])
	if test "x$JAVASCRIPT" = "xno" ; then
		HAVE_JAVASCRIPT=no
	else
		HAVE_JAVASCRIPT=yes
		AC_DEFINE([_HAVE_JAVASCRIPT_],[1],[with javascript])
	fi
	AC_MSG_RESULT($HAVE_JAVASCRIPT)
	AM_CONDITIONAL([JAVASCRIPT],[test x$HAVE_JAVASCRIPT = xyes])
	JAVASCRIPTWRAPPEREXT=.js
	AC_SUBST([JAVASCRIPTWRAPPEREXT])

	dnl }}}
	dnl triangle {{{
	AC_ARG_WITH([triangle-dir],
			  AS_HELP_STRING([--with-triangle-dir=DIR], [triangle root directory.]),
			 [TRIANGLE_ROOT=$withval],[TRIANGLE_ROOT="no"]) 

  dnl Check whether triangle is enabled
	AC_MSG_CHECKING([for triangle])
	if test "x$TRIANGLE_ROOT" = "xno" ; then
		HAVE_TRIANGLE=no
	else
		HAVE_TRIANGLE=yes
		if ! test -d "$TRIANGLE_ROOT"; then
			AC_MSG_ERROR([triangle directory provided ($TRIANGLE_ROOT) does not exist]);
		fi
		if ! test -f "$TRIANGLE_ROOT/triangle.h" ; then
			AC_MSG_ERROR([Couldn't find triangle.h... check your installation of triangle])
		fi
	fi
	AC_MSG_RESULT($HAVE_TRIANGLE)
	AM_CONDITIONAL([TRIANGLE],[test x$HAVE_TRIANGLE = xyes])

	dnl library and header files
	if test "x$HAVE_TRIANGLE" = "xyes"; then
		TRIANGLEINCL=-I$TRIANGLE_ROOT/
		case "${host_os}" in
				*cygwin*)
				TRIANGLEINCL="/I`cygpath -m $TRIANGLE_ROOT/`"
				TRIANGLELIB="-Wl,`cygpath -m $TRIANGLE_ROOT/`triangle.lib"
				;;
				*linux*)
				TRIANGLELIB=$TRIANGLE_ROOT/triangle.a
				if test "x$HAVE_JAVASCRIPT" = "xyes"; then
					dnl go to the bit code, not the library.
					TRIANGLELIB=$TRIANGLE_ROOT/triangle.o 
				else
					TRIANGLELIB=$TRIANGLE_ROOT/triangle.a
				fi
				;;
				*darwin*)
				if test "x$HAVE_JAVASCRIPT" = "xyes"; then
					dnl go to the bit code, not the library.
					TRIANGLELIB=$TRIANGLE_ROOT/triangle.o 
				else
					TRIANGLELIB=$TRIANGLE_ROOT/triangle.a
				fi
				;;
			esac
		AC_DEFINE([_HAVE_TRIANGLE_],[1],[with Triangle in ISSM src])
		AC_SUBST([TRIANGLEINCL])
		AC_SUBST([TRIANGLELIB])
	fi
	dnl }}}
	dnl boost{{{
	AC_ARG_WITH([boost-dir],
	  AS_HELP_STRING([--with-boost-dir=DIR], [boost root directory.]),
	  [BOOST_ROOT=$withval],[BOOST_ROOT="no"]) 

	dnl Check whether boost is enabled
	AC_MSG_CHECKING([for boost])
	if test "x$BOOST_ROOT" = "xno" ; then
		HAVE_BOOST=no
	else
		HAVE_BOOST=yes
		if ! test -d "$BOOST_ROOT"; then
			AC_MSG_ERROR([boost directory provided ($BOOST_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_BOOST)
	AM_CONDITIONAL([BOOST],[test x$HAVE_BOOST = xyes])

	dnl library and header files
	if test "x$HAVE_BOOST" = "xyes"; then
		BOOSTINCL=-I$BOOST_ROOT/include
		BOOSTLIB="-L$BOOST_ROOT/lib -lboost_python"
		AC_DEFINE([_HAVE_BOOST_],[1],[with Boost in ISSM src])
		AC_SUBST([BOOSTINCL])
		AC_SUBST([BOOSTLIB])
	fi
	dnl }}}
	dnl dakota{{{
	AC_ARG_WITH([dakota-dir],
	  AS_HELP_STRING([--with-dakota-dir=DIR], [dakota root directory.]),
	  [DAKOTA_ROOT=$withval],[DAKOTA_ROOT="no"]) 
	
	dnl Check whether dakota is enabled
	AC_MSG_CHECKING([for dakota])
	if test "x$DAKOTA_ROOT" = "xno" ; then
		HAVE_DAKOTA=no
	else
		HAVE_DAKOTA=yes
		if ! test -d "$DAKOTA_ROOT"; then
			AC_MSG_ERROR([dakota directory provided ($DAKOTA_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_DAKOTA)
	AM_CONDITIONAL([DAKOTA],[test x$HAVE_DAKOTA = xyes])

	dnl library and header files
	if test "x$HAVE_DAKOTA" = "xyes"; then
		DAKOTAINCL=-I$DAKOTA_ROOT/include

		AC_MSG_CHECKING(for dakota version)
		if test -f "$DAKOTA_ROOT/VERSION"; then
		 DAKOTA_VERSION=`cat $DAKOTA_ROOT/VERSION | grep 'DAKOTA Version' | sed 's/.*DAKOTA Version //' | sed 's/ .*//' `
		else if test -f "$DAKOTA_ROOT/../src/src/CommandLineHandler.C"; then
		 DAKOTA_VERSION=`cat $DAKOTA_ROOT/../src/src/CommandLineHandler.C | grep 'DAKOTA version' | grep 'release' | grep -v // | sed 's/.*DAKOTA version //' | sed 's/ .*//' `
		else if test -f "$DAKOTA_ROOT/../src/src/CommandLineHandler.cpp"; then
		 DAKOTA_VERSION=`cat $DAKOTA_ROOT/../src/src/CommandLineHandler.cpp | grep 'DAKOTA version' | grep 'release' | grep -v // | sed 's/.*DAKOTA version //' | sed 's/ .*//' `
		else
		 AC_MSG_ERROR([Dakota CommandLineHandler.C or CommandLineHandler.cpp file not found to determine DAKOTA_VERSION!]);
		fi
		fi
		fi
		AC_MSG_RESULT($DAKOTA_VERSION)
		AC_DEFINE_UNQUOTED([_DAKOTA_VERSION_],"$DAKOTA_VERSION",[Dakota version number])

		DAKOTAFLAGS=""
		case "${host_os}" in
			*cygwin*)
				if test x$DAKOTA_VERSION = x5.1 || test x$DAKOTA_VERSION = x5.2; then
					DAKOTALIB="-L$DAKOTA_ROOT/lib -L$BOOST_ROOT/lib -ldakota -lteuchos -lpecos -llhs -lsparsegrid -lsurfpack -lconmin -lddace -lfsudace -ljega -lcport -loptpp -lpsuade -lncsuopt -lcolin -linterfaces -lmomh -lscolib -lpebbl -ltinyxml -lutilib -l3po -lhopspack -lnidr -lamplsolver -lboost_signals -lboost_regex -lboost_filesystem"
				else if test x$DAKOTA_VERSION = x6.1 || test x$DAKOTA_VERSION = x6.2; then
				   DAKOTAFLAGS="-DDISABLE_DAKOTA_CONFIG_H -DBOOST_MULTI_INDEX_DISABLE_SERIALIZATION -DDAKOTA_PLUGIN -DBOOST_DISABLE_ASSERTS -DDAKOTA_HAVE_BOOST_FS -DHAVE_UNISTD_H -DHAVE_SYSTEM -DHAVE_WORKING_FORK -DHAVE_WORKING_VFORK -DHAVE_SYS_WAIT_H -DHAVE_USLEEP -DDAKOTA_F90 -DDAKOTA_HAVE_MPI -DHAVE_PECOS -DHAVE_SURFPACK -DDAKOTA_COLINY -DDAKOTA_UTILIB -DHAVE_ADAPTIVE_SAMPLING -DHAVE_CONMIN -DDAKOTA_DDACE -DHAVE_FSUDACE -DDAKOTA_HOPS -DHAVE_JEGA -DHAVE_NCSU -DHAVE_NL2SOL -DHAVE_OPTPP -DDAKOTA_OPTPP -DHAVE_PSUADE -DHAVE_AMPL"
					DAKOTALIB="-L$DAKOTA_ROOT/lib -L$BOOST_ROOT/lib -ldakota_src -ldream -lfsudace -lddace -lnomad -lpecos_src -lscolib -ljega_fe -llhs -lpebbl -lcolin -linterfaces -llhs_mods -lmoga -loptpp -lsoga -lsurfpack -lutilib -lconmin -ldakota_src_fortran -llhs_mod -lncsuopt -lsurfpack_fortran -lteuchos -l3po -lamplsolver -lcport -ldfftpack -leutils -lfsudace -lhopspack -ljega -lnidr -lpecos -lpsuade -ltinyxml -lutilities -lsparsegrid -lboost_serialization -lboost_signals -lboost_regex -lboost_filesystem -lboost_system"
					AC_DEFINE([DISABLE_DAKOTA_CONFIG_H],[1],[disabling DAKOTA_CONFIG_H])
					AC_DEFINE([DAKOTA_HAVE_MPI],[1],[enabling parallel MPI])
				else
					AC_MSG_ERROR([Dakota version not found or version ($DAKOTA_VERSION) not supported!]);
				fi
				fi
			;;
			*linux*)
				if test x$DAKOTA_VERSION = x5.1 || test x$DAKOTA_VERSION = x5.2; then
					DAKOTALIB="-L$DAKOTA_ROOT/lib -ldakota -lteuchos -lpecos -llhs -lsparsegrid -lsurfpack -lconmin -lddace -lfsudace -ljega -lcport -loptpp -lpsuade -lncsuopt -lcolin -linterfaces -lmomh -lscolib -lpebbl -ltinyxml -lutilib -l3po -lhopspack -lnidr -lamplsolver -lboost_signals -lboost_regex -lboost_filesystem -lboost_system -ldl"
				else if test x$DAKOTA_VERSION = x5.3 || test x$DAKOTA_VERSION = x5.3.1; then
					DAKOTAFLAGS="-DDISABLE_DAKOTA_CONFIG_H -DBOOST_MULTI_INDEX_DISABLE_SERIALIZATION -DDAKOTA_PLUGIN -DBOOST_DISABLE_ASSERTS -DDAKOTA_HAVE_BOOST_FS -DHAVE_UNISTD_H -DHAVE_SYSTEM -DHAVE_WORKING_FORK -DHAVE_WORKING_VFORK -DHAVE_SYS_WAIT_H -DHAVE_USLEEP -DDAKOTA_F90 -DDAKOTA_HAVE_MPI -DHAVE_PECOS -DHAVE_SURFPACK -DDAKOTA_COLINY -DDAKOTA_UTILIB -DHAVE_ADAPTIVE_SAMPLING -DHAVE_CONMIN -DDAKOTA_DDACE -DHAVE_FSUDACE -DDAKOTA_HOPS -DHAVE_JEGA -DHAVE_NCSU -DHAVE_NL2SOL -DHAVE_OPTPP -DDAKOTA_OPTPP -DHAVE_PSUADE -DHAVE_AMPL"
					DAKOTALIB="-L$DAKOTA_ROOT/lib -L$BOOST_ROOT/lib -ldakota_src -lpecos_src -lscolib -ljega_fe -llhs -lpebbl -lcolin -linterfaces -lmods -lmoga -loptpp -lsampling -lsoga -lsurfpack -lutilib -lconmin -ldakota_src_fortran -lmod -lncsuopt -lsurfpack_fortran -lteuchos -l3po -lamplsolver -lanalyzer -lbose -lcport -ldace -ldfftpack -leutils -lfsudace -lhopspack -ljega -lnidr -lpecos -lpsuade -lrandom -ltinyxml -lutilities -lsparsegrid -lboost_signals -lboost_regex -lboost_filesystem -lboost_system"
					AC_DEFINE([DISABLE_DAKOTA_CONFIG_H],[1],[disabling DAKOTA_CONFIG_H])
					AC_DEFINE([DAKOTA_HAVE_MPI],[1],[enabling parallel MPI])
				else if test x$DAKOTA_VERSION = x6.1 || test x$DAKOTA_VERSION = x6.2; then
				   DAKOTAFLAGS="-DDISABLE_DAKOTA_CONFIG_H -DBOOST_MULTI_INDEX_DISABLE_SERIALIZATION -DDAKOTA_PLUGIN -DBOOST_DISABLE_ASSERTS -DDAKOTA_HAVE_BOOST_FS -DHAVE_UNISTD_H -DHAVE_SYSTEM -DHAVE_WORKING_FORK -DHAVE_WORKING_VFORK -DHAVE_SYS_WAIT_H -DHAVE_USLEEP -DDAKOTA_F90 -DDAKOTA_HAVE_MPI -DHAVE_PECOS -DHAVE_SURFPACK -DDAKOTA_UTILIB -DHAVE_ADAPTIVE_SAMPLING -DHAVE_CONMIN -DDAKOTA_DDACE -DHAVE_FSUDACE -DDAKOTA_HOPS -DHAVE_NCSU -DHAVE_NL2SOL -DHAVE_OPTPP -DDAKOTA_OPTPP -DHAVE_PSUADE -DHAVE_AMPL"
					DAKOTALIB="-L$DAKOTA_ROOT/lib -L$BOOST_ROOT/lib -ldakota_src -ldream -lfsudace -lddace -lnomad -lpecos_src -llhs -llhs_mods -loptpp -lsurfpack -lconmin -ldakota_src_fortran -llhs_mod -lncsuopt -lsurfpack_fortran -lteuchos -lamplsolver -lcport -ldfftpack -lfsudace -lhopspack -lnidr -lpecos -lpsuade -lsparsegrid -lboost_serialization -lboost_signals -lboost_regex -lboost_filesystem -lboost_system"
					AC_DEFINE([DISABLE_DAKOTA_CONFIG_H],[1],[disabling DAKOTA_CONFIG_H])
					AC_DEFINE([DAKOTA_HAVE_MPI],[1],[enabling parallel MPI])
				else
					AC_MSG_ERROR([Dakota version not found or version ($DAKOTA_VERSION) not supported!]);
				fi
				fi
				fi
			;;
			*darwin*)
				if test x$DAKOTA_VERSION = x5.1 || test x$DAKOTA_VERSION = x5.2; then
					DAKOTALIB="-L$DAKOTA_ROOT/lib -ldakota -lteuchos -lpecos -llhs -lsparsegrid -lsurfpack -lconmin -lddace -lfsudace -ljega -lcport -loptpp -lpsuade -lncsuopt -lcolin -linterfaces -lmomh -lscolib -lpebbl -ltinyxml -lutilib -l3po -lhopspack -lnidr -lamplsolver -lboost_signals -lboost_regex -lboost_filesystem -lboost_system"
					dnl DAKOTALIB+= "-lgslcblas -L/usr/lib -lblas -llapack"
				else if test x$DAKOTA_VERSION = x5.3 || test x$DAKOTA_VERSION = x5.3.1; then
					DAKOTAFLAGS="-DDISABLE_DAKOTA_CONFIG_H -DBOOST_MULTI_INDEX_DISABLE_SERIALIZATION -DDAKOTA_PLUGIN -DBOOST_DISABLE_ASSERTS -DDAKOTA_HAVE_BOOST_FS -DHAVE_UNISTD_H -DHAVE_SYSTEM -DHAVE_WORKING_FORK -DHAVE_WORKING_VFORK -DHAVE_SYS_WAIT_H -DHAVE_USLEEP -DDAKOTA_F90 -DDAKOTA_HAVE_MPI -DHAVE_PECOS -DHAVE_SURFPACK -DDAKOTA_COLINY -DDAKOTA_UTILIB -DHAVE_ADAPTIVE_SAMPLING -DHAVE_CONMIN -DDAKOTA_DDACE -DHAVE_FSUDACE -DDAKOTA_HOPS -DHAVE_JEGA -DHAVE_NCSU -DHAVE_NL2SOL -DHAVE_OPTPP -DDAKOTA_OPTPP -DHAVE_PSUADE -DHAVE_AMPL"
					DAKOTALIB="-L$DAKOTA_ROOT/lib -L$BOOST_ROOT/lib -ldakota_src -lpecos_src -lscolib -ljega_fe -llhs -lpebbl -lcolin -linterfaces -lmods -lmoga -loptpp -lsampling -lsoga -lsurfpack -lutilib -lconmin -ldakota_src_fortran -lmod -lncsuopt -lsurfpack_fortran -lteuchos -l3po -lamplsolver -lanalyzer -lbose -lcport -ldace -ldfftpack -leutils -lfsudace -lhopspack -ljega -lnidr -lpecos -lpsuade -lrandom -ltinyxml -lutilities -lsparsegrid -lboost_signals -lboost_regex -lboost_filesystem -lboost_system"
					AC_DEFINE([DISABLE_DAKOTA_CONFIG_H],[1],[disabling DAKOTA_CONFIG_H])
					AC_DEFINE([DAKOTA_HAVE_MPI],[1],[enabling parallel MPI])
				else if test x$DAKOTA_VERSION = x6.1 || test x$DAKOTA_VERSION = x6.2; then
				   DAKOTAFLAGS="-DDISABLE_DAKOTA_CONFIG_H -DBOOST_MULTI_INDEX_DISABLE_SERIALIZATION -DDAKOTA_PLUGIN -DBOOST_DISABLE_ASSERTS -DDAKOTA_HAVE_BOOST_FS -DHAVE_UNISTD_H -DHAVE_SYSTEM -DHAVE_WORKING_FORK -DHAVE_WORKING_VFORK -DHAVE_SYS_WAIT_H -DHAVE_USLEEP -DDAKOTA_F90 -DDAKOTA_HAVE_MPI -DHAVE_PECOS -DHAVE_SURFPACK -DDAKOTA_UTILIB -DHAVE_ADAPTIVE_SAMPLING -DHAVE_CONMIN -DDAKOTA_DDACE -DHAVE_FSUDACE -DDAKOTA_HOPS -DHAVE_NCSU -DHAVE_NL2SOL -DHAVE_OPTPP -DDAKOTA_OPTPP -DHAVE_PSUADE -DHAVE_AMPL"
					DAKOTALIB="-L$DAKOTA_ROOT/lib -L$BOOST_ROOT/lib -ldakota_src -ldream -lfsudace -lddace -lnomad -lpecos_src -llhs -llhs_mods -loptpp -lsurfpack -lconmin -ldakota_src_fortran -llhs_mod -lncsuopt -lsurfpack_fortran -lteuchos -lamplsolver -lcport -ldfftpack -lfsudace -lhopspack -lnidr -lpecos -lpsuade -lsparsegrid -lboost_serialization -lboost_signals -lboost_regex -lboost_filesystem -lboost_system"
					AC_DEFINE([DISABLE_DAKOTA_CONFIG_H],[1],[disabling DAKOTA_CONFIG_H])
					AC_DEFINE([DAKOTA_HAVE_MPI],[1],[enabling parallel MPI])
				else
					AC_MSG_ERROR([Dakota version not found or version ($DAKOTA_VERSION) not supported!]);
				fi
				fi
				fi
			;;
		esac

		case $DAKOTA_VERSION in
			@<:@1-9@:>@.@<:@0-9@:>@.@<:@0-9@:>@)
				DAKOTA_MAJOR=`echo $DAKOTA_VERSION | sed -e 's/^\(@<:@0-9@:>@*\)\..*/\1/'`
				DAKOTA_MINOR=`echo $DAKOTA_VERSION | sed -e 's/^@<:@0-9@:>@*\.\(@<:@0-9@:>@*\)\..*/\1/'`
				DAKOTA_BUILD=`echo $DAKOTA_VERSION | sed -e 's/^@<:@0-9@:>@*\.@<:@0-9@:>@*\.\(@<:@0-9@:>@*\).*/\1/'`
			;;
			@<:@1-9@:>@.@<:@0-9@:>@)
				DAKOTA_MAJOR=`echo $DAKOTA_VERSION | sed -e 's/^\(@<:@0-9@:>@*\)\..*/\1/'`
				DAKOTA_MINOR=`echo $DAKOTA_VERSION | sed -e 's/^@<:@0-9@:>@*\.\(@<:@0-9@:>@*\).*/\1/'`
				DAKOTA_BUILD=0
			;;
			@<:@1-9@:>@.@<:@0-9@:>@+)
				DAKOTA_MAJOR=`echo $DAKOTA_VERSION | sed -e 's/^\(@<:@0-9@:>@*\)\..*/\1/'`
				DAKOTA_MINOR=`echo $DAKOTA_VERSION | sed -e 's/^@<:@0-9@:>@*\.\(@<:@0-9@:>@*\).*/\1/'`
				DAKOTA_BUILD=0
			;;
			*)
				AC_MSG_ERROR([Dakota version ($DAKOTA_VERSION) not supported!]);
		   ;;
		esac
		AC_MSG_CHECKING(for dakota major version)
		AC_MSG_RESULT($DAKOTA_MAJOR)
		AC_DEFINE_UNQUOTED([_DAKOTA_MAJOR_],$DAKOTA_MAJOR,[Dakota major version number])
		AC_MSG_CHECKING(for dakota minor version)
		AC_MSG_RESULT($DAKOTA_MINOR)
		AC_DEFINE_UNQUOTED([_DAKOTA_MINOR_],$DAKOTA_MINOR,[Dakota minor version number])
		AC_MSG_CHECKING(for dakota build version)
		AC_MSG_RESULT($DAKOTA_BUILD)
		AC_DEFINE_UNQUOTED([_DAKOTA_BUILD_],$DAKOTA_BUILD,[Dakota build version number])

		AC_DEFINE([_HAVE_DAKOTA_],[1],[with Dakota in ISSM src])
		AC_SUBST([DAKOTAINCL])
		AC_SUBST([DAKOTAFLAGS])
		AC_SUBST([DAKOTALIB])
	fi
	AM_CONDITIONAL([ISSM_DAKOTA],[test x$DAKOTA_MAJOR = x6])
	dnl }}}
	dnl python{{{
	AC_ARG_WITH([python-dir],
	  AS_HELP_STRING([--with-python-dir=DIR], [python root directory.]),
	  [PYTHON_ROOT=$withval],[PYTHON_ROOT="no"]) 

	dnl Check whether python is enabled
	AC_MSG_CHECKING([for python])
	if test "x$PYTHON_ROOT" = "xno" ; then
		HAVE_PYTHON=no
		HAVE_PYTHON3=no
	else
		HAVE_PYTHON=yes
		if ! test -d "$PYTHON_ROOT"; then
			AC_MSG_ERROR([python directory provided ($PYTHON_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_PYTHON)
	AM_CONDITIONAL([PYTHON],[test x$HAVE_PYTHON = xyes])

	dnl python specifics
	if test "x$HAVE_PYTHON" = "xyes"; then
		AC_MSG_CHECKING([for python version])
		dnl Query Python for its version number.  Getting [:3] seems to be the
		dnl best way to do this; it's what "site.py" does in the standard library.
		PYTHON_VERSION=$($PYTHON_ROOT/bin/python -c "import sys; print sys.version[[:3]]")
		AC_MSG_RESULT($PYTHON_VERSION)

		dnl recover major 
		PYTHON_MAJOR=${PYTHON_VERSION%.*}
		AC_DEFINE_UNQUOTED([_PYTHON_MAJOR_],$PYTHON_MAJOR,[python version major])
		if test "x$PYTHON_MAJOR" = "x3"; then
			HAVE_PYTHON3="yes"
		else
			HAVE_PYTHON3="no"
		fi

		AC_MSG_CHECKING([for python header file Python.h])
		dnl Python.h mighty be in different locations:
		if test -f "$PYTHON_ROOT/include/Python.h"; then
			PYTHONINCL=-I$PYTHON_ROOT/include
		else if test -f "$PYTHON_ROOT/include/python$PYTHON_VERSION/Python.h"; then
			PYTHONINCL=-I$PYTHON_ROOT/include/python$PYTHON_VERSION
		else
			AC_MSG_ERROR([Python.h not found, locate this file and contact ISSM developers]);
		fi
		fi
		AC_MSG_RESULT(found)

		PYTHONLIB="-L$PYTHON_ROOT/lib -lpython$PYTHON_VERSION"
		PYTHONEXT=.so
		case "${host_os}" in
			*cygwin*)
			PYTHONLINK="-shared"
			;;
			*linux*)
			PYTHONLINK="-shared"
			;;
			*darwin*)
			PYTHONLINK="-dynamiclib"
			;;
		esac
		AC_DEFINE([_HAVE_PYTHON_],[1],[with python in ISSM src])
		AC_SUBST([PYTHONINCL])
		AC_SUBST([PYTHONLIB])
		PYTHONWRAPPEREXT=$PYTHONEXT
		AC_SUBST([PYTHONWRAPPEREXT])
		AC_SUBST([PYTHONLINK])
	fi
	AM_CONDITIONAL([PYTHON3], [test x$HAVE_PYTHON3 = xyes])
	dnl }}}
	dnl python-numpy{{{
	AC_ARG_WITH([python-numpy-dir],
	  AS_HELP_STRING([--with-python-numpy-dir=DIR], [python-numpy root directory.]),
	  [PYTHON_NUMPY_ROOT=$withval],[PYTHON_NUMPY_ROOT="no"]) 
	
	dnl Check whether numpy is enabled
	AC_MSG_CHECKING(for python-numpy)
	if test "x$PYTHON_NUMPY_ROOT" = "xno" ; then
		HAVE_PYTHON_NUMPY=no
	else
		HAVE_PYTHON_NUMPY=yes
		if ! test -d "$PYTHON_NUMPY_ROOT"; then
			AC_MSG_ERROR([numpy directory provided ($PYTHON_NUMPY_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_PYTHON_NUMPY)

	dnl numpy lib
	if test "x$HAVE_PYTHON_NUMPY" = "xyes"; then
		PYTHON_NUMPYINCL="-I$PYTHON_NUMPY_ROOT -I$PYTHON_NUMPY_ROOT/core/include/numpy"
		AC_DEFINE([_HAVE_PYTHON_NUMPY_],[1],[with Python-Numpy in ISSM src])
		AC_SUBST([PYTHON_NUMPYINCL])
	fi
	dnl }}}
	dnl chaco{{{
	AC_ARG_WITH([chaco-dir],
	  AS_HELP_STRING([--with-chaco-dir=DIR], [chaco root directory.]),
	  [CHACO_ROOT=$withval],[CHACO_ROOT="no"]) 
	
	dnl Check whether chaco is enabled
	AC_MSG_CHECKING([for chaco])
	if test "x$CHACO_ROOT" = "xno" ; then
		HAVE_CHACO=no
	else
		HAVE_CHACO=yes
		if ! test -d "$CHACO_ROOT"; then
			AC_MSG_ERROR([chaco directory provided ($CHACO_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_CHACO)
	AM_CONDITIONAL([CHACO],[test x$HAVE_CHACO = xyes])

	dnl library and header files
	if test "x$HAVE_CHACO" = "xyes"; then
		CHACOINCL=-I$CHACO_ROOT/include
		CHACOLIB="-L$CHACO_ROOT/lib -lchacominusblas"
		AC_DEFINE([_HAVE_CHACO_],[1],[with Chaco in ISSM src])
		AC_SUBST([CHACOINCL])
		AC_SUBST([CHACOLIB])
	fi
	dnl }}}
	dnl scotch{{{
	AC_ARG_WITH([scotch-dir],
	  AS_HELP_STRING([--with-scotch-dir=DIR], [scotch root directory.]),
	  [SCOTCH_ROOT=$withval],[SCOTCH_ROOT="no"]) 

  dnl Check whether scotch is enabled
	AC_MSG_CHECKING([for scotch])
	if test "x$SCOTCH_ROOT" = "xno" ; then
		HAVE_SCOTCH=no
	else
		HAVE_SCOTCH=yes
		if ! test -d "$SCOTCH_ROOT"; then
			AC_MSG_ERROR([scotch directory provided ($SCOTCH_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_SCOTCH)
	AM_CONDITIONAL([SCOTCH],[test x$HAVE_SCOTCH = xyes])
	
	dnl scotch libraries
	if test "x$HAVE_SCOTCH" = "xyes"; then
		SCOTCHINCL="-DNOFILEIO -I$SCOTCH_ROOT/include -DSCOTCH_VERSION=\\\"UNKNOWN\\\""
		SCOTCHLIB="-L$SCOTCH_ROOT/lib -lnfioscotch -lnfioscotcherr -lnfioscotcherrexit -lscotchmetis"
		AC_DEFINE([_HAVE_SCOTCH_],[1],[with Scotch in ISSM src])
		AC_SUBST([SCOTCHINCL])
		AC_SUBST([SCOTCHLIB])
	fi
	dnl }}}
	dnl esmf{{{
	AC_ARG_WITH([esmf-dir],
		AS_HELP_STRING([--with-esmf-dir=DIR], [esmf root directory.]),
		[ESMF_ROOT=$withval],[ESMF_ROOT="no"]) 

	dnl Check whether esmf is enabled
	AC_MSG_CHECKING([for esmf])
	if test "x$ESMF_ROOT" = "xno" ; then
		HAVE_ESMF=no
	else
		HAVE_ESMF=yes
		if ! test -d "$ESMF_ROOT"; then
			AC_MSG_ERROR([esmf directory provided ($ESMF_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_ESMF)
	
	dnl esmf headers and libraries
	if test "x$HAVE_ESMF" == "xyes"; then
		ESMFINCL="-I$ESMF_ROOT/include"
		ESMFLIB="-L$ESMF_ROOT/lib -lesmf"
		AC_DEFINE([_HAVE_ESMF_],[1],[with esmf in ISSM src])
		AC_SUBST([ESMFINCL])
		AC_SUBST([ESMFLIB])
	fi
	AM_CONDITIONAL([ESMF], [test x$HAVE_ESMF = xyes])
	dnl }}}
	dnl adolc{{{
	AC_ARG_WITH([adolc-dir],
		AS_HELP_STRING([--with-adolc-dir=DIR], [adolc root directory.]),
		[ADOLC_ROOT=$withval],[ADOLC_ROOT="no"]) 

	dnl Check whether adolc is enabled
	AC_MSG_CHECKING([for adolc])
	if test "x$ADOLC_ROOT" = "xno" ; then
		HAVE_ADOLC=no
	else
		HAVE_ADOLC=yes
		if ! test -d "$ADOLC_ROOT"; then
			AC_MSG_ERROR([adolc directory provided ($ADOLC_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_ADOLC)
	
	dnl adolc headers and libraries
	if test "x$HAVE_ADOLC" == "xyes"; then
		ADOLCINCL="-I$ADOLC_ROOT/include"
		dnl ADOLCLIB="-L$ADOLC_ROOT/lib64 -ladolc" used to be the path
		ADOLCLIB="-L$ADOLC_ROOT/lib -ladolc"
		AC_DEFINE([_HAVE_ADOLC_],[1],[with adolc in ISSM src])
		AC_DEFINE([_HAVE_AD_],[1],[with AD in ISSM src])
		AC_SUBST([ADOLCINCL])
		AC_SUBST([ADOLCLIB])
	fi
	AM_CONDITIONAL([ADOLC], [test x$HAVE_ADOLC = xyes])
        AM_COND_IF(ADOLC,[CXXFLAGS+=" -std=c++11"])
	dnl }}}
	dnl adolc-version{{{
	AC_ARG_WITH([adolc-version],
		AS_HELP_STRING([--with-adolc-version=number], [adolc version.]),
		[ADOLC_VERSION=$withval],[ADOLC_VERSION=2]) 
	AC_MSG_CHECKING(for adolc-version) 

	AC_DEFINE_UNQUOTED([_ADOLC_VERSION_],$ADOLC_VERSION,[ADOLC version])
	AC_MSG_RESULT($ADOLC_VERSION)
	dnl }}}
	dnl adic2{{{
	AC_ARG_WITH([adic2-dir],
	  AS_HELP_STRING([--with-adic2-dir=DIR], [adic2 root directory.]),
	  [ADIC2_ROOT=$withval],[ADIC2_ROOT="no"]) 

	dnl Check whether adic2 is enabled
	AC_MSG_CHECKING([for adic2])
	if test "x$ADIC2_ROOT" = "xno" ; then
		HAVE_ADIC2=no
	else
		HAVE_ADIC2=yes
		if ! test -d "$ADIC2_ROOT"; then
			AC_MSG_ERROR([adic2 directory provided ($ADIC2_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_ADIC2)

	dnl adic2 headers and libraries
	if test "x$HAVE_ADIC2" == "xyes"; then
		ADIC2INCL="-DADIC2_DENSE -I$ADIC2_ROOT/include -I$ADIC2_ROOT/share/runtime_dense/"
		ADIC2LIB=""
		AC_DEFINE([_HAVE_ADIC2_],[1],[with adic2 in ISSM src])
		AC_SUBST([ADIC2INCL])
		AC_SUBST([ADIC2LIB])
	fi
	AM_CONDITIONAL([ADIC2], [test x$HAVE_ADIC2 = xyes])
	dnl }}}
	dnl atlas{{{
	AC_ARG_WITH([atlas-dir],
	  AS_HELP_STRING([--with-atlas-dir=DIR],[atlas root directory]),
	  [ATLAS_ROOT=$withval],[ATLAS_ROOT="no"])
			  
	dnl Check whether atlas is enabled
	AC_MSG_CHECKING(for atlas and cblas libraries)
	if test "x$ATLAS_ROOT" = "xno" ; then
		HAVE_ATLAS=no
	else
		HAVE_ATLAS=yes
		if ! test -d "$ATLAS_ROOT"; then
			AC_MSG_ERROR([atlas directory provided ($ATLAS_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_ATLAS)

	dnl atlas headers and libraries
	if test "x$HAVE_ATLAS" == "xyes"; then
		dnl: branch on whether we are running on windows or linux.
		case "${host_os}" in
			*cygwin*)
			ATLASLIB="-L`cygpath -m $ATLAS_ROOT` -Wl,libatlas.lib  -Wl,libcblas.lib"
			;;
			*linux*)
			ATLASLIB=-L"$ATLAS_ROOT/lib -lcblas -latlas -lm " 
			;;
			*darwin*)
			ATLASLIB=-L"$ATLAS_ROOT/lib -lcblas -latlas -lm" 
			;;
		esac
		AC_DEFINE([_HAVE_ATLAS_],[1],[with ATLAS in ISSM src])
		AC_SUBST([ATLASLIB])
	fi
	dnl }}}
	dnl gsl{{{
	AC_ARG_WITH([gsl-dir],
	  AS_HELP_STRING([--with-gsl-dir=DIR], [gsl root directory.]),
	  [GSL_ROOT=$withval],[GSL_ROOT="no"]) 

	dnl Check whether gsl is enabled
	AC_MSG_CHECKING([for gsl])
	if test "x$GSL_ROOT" = "xno" ; then
		HAVE_GSL=no
	else
		HAVE_GSL=yes
		if ! test -d "$GSL_ROOT"; then
			AC_MSG_ERROR([gsl directory provided ($GSL_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_GSL)
	
	dnl gsl headers and libraries
	if test "x$HAVE_GSL" == "xyes"; then
		GSLINCL="-I$GSL_ROOT/include"
		if test "x$HAVE_ATLAS" = "xyes" ; then
			GSLLIB="-dy -L$GSL_ROOT/lib -lgsl -L$ATLAS_ROOT/lib -lcblas -latlas -lm"
		else
			GSLLIB="-L$GSL_ROOT/lib -lgsl -lgslcblas -lm"
		fi
		AC_DEFINE([_HAVE_GSL_],[1],[with gsl in ISSM src])
		AC_SUBST([GSLINCL])
		AC_SUBST([GSLLIB])
	fi
	AM_CONDITIONAL([GSL], [test x$HAVE_GSL = xyes])
	dnl }}}
	dnl adjoinable-mpi{{{
	AC_ARG_WITH([ampi-dir],
	  AS_HELP_STRING([--with-ampi-dir=DIR], [adjoinable mpi root directory.]),
	  [AMPI_ROOT=$withval],[AMPI_ROOT="no"]) 

	dnl Check whether ampi is enabled
	AC_MSG_CHECKING([for ampi])
	if test "x$AMPI_ROOT" = "xno" ; then
		HAVE_AMPI=no
	else
		HAVE_AMPI=yes
		if ! test -d "$AMPI_ROOT"; then
			AC_MSG_ERROR([ampi directory provided ($AMPI_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_AMPI)
	
	dnl ampi headers and libraries
	if test "x$HAVE_AMPI" == "xyes"; then
		AMPIINCL="-I$AMPI_ROOT/include"
		if test "x$ADOLC_ROOT" == "xno"; then
			AC_MSG_ERROR([cannot run adjoinable mpi without adolc]);
		fi
		dnl AMPILIB="-dy -L$AMPI_ROOT/lib -lampiCommon -L$ADOLC_ROOT/lib -ladolc -L$AMPI_ROOT/lib -lampiCommon -lampiBookkeeping -lampiTape"
		dnl AMPILIB="-dy -L$AMPI_ROOT/lib  -L$ADOLC_ROOT/lib -Wl,--start-group,-lampiCommon,-ladolc,-lampiCommon,-lampiBookkeeping,-lampiTape,-lampiPlainC,-lampiADtoolStubsST,--end-group"
		dnl AMPILIB="-L$AMPI_ROOT/lib  -L$ADOLC_ROOT/lib -Wl,--start-group -lampiCommon -ladolc -lampiCommon -lampiBookkeeping -lampiTape -lampiPlainC -lampiADtoolStubsST -Wl,--end-group"
		dnl AMPILIB="$AMPI_ROOT/lib/libampiCommon.so $ADOLC_ROOT/lib/libadolc.so  $AMPI_ROOT/lib/libampiCommon.so $AMPI_ROOT/lib/libampiBookkeeping.so $AMPI_ROOT/lib/libampiTape.so $AMPI_ROOT/lib/libampiPlainC.so  $AMPI_ROOT/lib/libampiADtoolStubsST.so"
		dnl AMPILIB="-dy -L$AMPI_ROOT/lib  -L$ADOLC_ROOT/lib -lampiCommon -ladolc -lampiCommon -lampiBookkeeping -lampiTape -lampiPlainC -lampiADtoolStubsST"
		AMPILIB="-dy -L$AMPI_ROOT/lib  -lampiCommon -lampiBookkeeping -lampiTape"
		AC_DEFINE([_HAVE_AMPI_],[1],[with adjoinable mpi in ISSM src])
		AC_SUBST([AMPIINCL])
		AC_SUBST([AMPILIB])
	fi
	AM_CONDITIONAL([AMPI], [test x$HAVE_AMPI = xyes])
	dnl }}}
	dnl rose{{{
	AC_ARG_WITH([rose-dir],
	  AS_HELP_STRING([--with-rose-dir=DIR], [rose root directory.]),
	  [ROSE_ROOT=$withval],[ROSE_ROOT="no"]) 

	dnl Check whether rose is enabled
	AC_MSG_CHECKING([for rose])
	if test "x$ROSE_ROOT" = "xno" ; then
		HAVE_ROSE=no
	else
		HAVE_ROSE=yes
		if ! test -d "$ROSE_ROOT"; then
			AC_MSG_ERROR([rose directory provided ($ROSE_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_ROSE)
	AM_CONDITIONAL([ROSE],[test x$HAVE_ROSE = xyes])

	dnl library and header files
	if test "x$HAVE_ROSE" = "xyes"; then
		ROSEINCL="-I$ROSE_ROOT/include"
		ROSELIB=""
		AC_DEFINE([_HAVE_ROSE_],[1],[with rose in ISSM src])
		AC_SUBST([ROSEINCL])
		AC_SUBST([ROSELIB])
	fi
	dnl }}}
	dnl mpi{{{
	AC_MSG_CHECKING(for mpi)
	
	AC_ARG_WITH([mpi-include],
	  AS_HELP_STRING([--with-mpi-include=DIR],[mpi include directory, necessary for parallel build]),
	  [MPI_INCLUDE=$withval],[MPI_INCLUDE=""])

	AC_ARG_WITH([mpi-libdir],
	  AS_HELP_STRING([--with-mpi-libdir=DIR],[mpi lib directory, necessary for parallel build]),
	  [MPI_LIBDIR=$withval],[MPI_LIBDIR=""])

	AC_ARG_WITH([mpi-libflags],
	  AS_HELP_STRING([--with-mpi-libflags=LIBS],[mpi libraries to be used, necessary for parallel build]),
	  [MPI_LIBFLAGS=$withval],[MPI_LIBFLAGS=""])

	
	if test -z "$MPI_INCLUDE" ; then
		HAVE_MPI=no
	else
		HAVE_MPI=yes

		dnl Processing for windows
		if  test x$VENDOR = xintel-win7-32; then
			MPI_LIBDIR=`cygpath -m $MPI_LIBDIR`
			MPI_INCLUDE=`cygpath -m $MPI_INCLUDE`
		elif test x$VENDOR = xintel-win7-64; then
			MPI_LIBDIR="/I`cygpath -m $MPI_LIBDIR`"
			MPI_INCLUDE=`cygpath -m $MPI_INCLUDE`
		elif test x$VENDOR = xMSVC-Win64 || test x$VENDOR = xMSVC-Win64-par; then
			MPI_LIBDIR=`cygpath -m $MPI_LIBDIR`
			MPI_INCLUDE=`cygpath -m $MPI_INCLUDE`
		fi

		if test -z "$MPI_LIBDIR"; then
			MPILIB="$MPI_LIBFLAGS"
		else
			MPILIB="-L$MPI_LIBDIR $MPI_LIBFLAGS"
		fi

		if  test x$IS_WINDOWS = xyes; then
			MPIINCL=/I"$MPI_INCLUDE"
		else 
			MPIINCL=-I"$MPI_INCLUDE"
		fi

		AC_DEFINE([_HAVE_MPI_],[1],[with Mpi in ISSM src])
		AC_DEFINE([HAVE_MPI],[1],[Mpi Flag for Dakota (DO NOT REMOVE)])
		AC_SUBST([MPIINCL])
		AC_SUBST([MPILIB])
	fi
	AM_CONDITIONAL([MPI], [test x$HAVE_MPI = xyes])
	AC_MSG_RESULT($HAVE_MPI)
	dnl }}}
	dnl petsc{{{
	AC_ARG_WITH([petsc-dir],
	  AS_HELP_STRING([--with-petsc-dir=DIR],[PETSc root directory, necessary for parallel build]),
	  [PETSC_ROOT=$withval],[PETSC_ROOT="no"])
		
	dnl Check whether petsc is enabled
	AC_MSG_CHECKING([for petsc])
	if test "x$PETSC_ROOT" = "xno" ; then
		HAVE_PETSC=no
	else
		HAVE_PETSC=yes
		if ! test -d "$PETSC_ROOT"; then
			AC_MSG_ERROR([petsc directory provided ($PETSC_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_PETSC)
	AM_CONDITIONAL([PETSC],[test x$HAVE_PETSC = xyes])

	dnl library and header files
	if test "x$HAVE_PETSC" = "xyes"; then
		AC_MSG_CHECKING(for petsc version)
	   if ! test -f "$PETSC_ROOT/include/petscversion.h"; then
			AC_MSG_ERROR([PETSc not instaled corretly: file ($PETSC_ROOT/include/petscversion.h) does not exist]);
		fi
		PETSC_MAJOR=`cat $PETSC_ROOT/include/petscversion.h | grep "#define PETSC_VERSION_MAJOR" | sed 's/#define PETSC_VERSION_MAJOR//' | sed 's/ //g'`
		PETSC_MINOR=`cat $PETSC_ROOT/include/petscversion.h | grep "#define PETSC_VERSION_MINOR" | sed 's/#define PETSC_VERSION_MINOR//' | sed 's/ //g'`
		AC_DEFINE_UNQUOTED([_PETSC_MAJOR_],$PETSC_MAJOR,[PETSc version major])
		AC_DEFINE_UNQUOTED([_PETSC_MINOR_],$PETSC_MINOR,[PETSc version minor])
		AC_MSG_RESULT($PETSC_MAJOR.$PETSC_MINOR)

		PETSC_VERSION_DATE_HG=`cat $PETSC_ROOT/include/petscversion.h | grep "#define PETSC_VERSION_DATE_HG" | sed 's/#define PETSC_VERSION_DATE_HG//' | sed 's/ //g' | sed -e 's/\"//g' `
		PETSC_RELEASE=`cat $PETSC_ROOT/include/petscversion.h | grep "#define PETSC_VERSION_RELEASE" | sed 's/#define PETSC_VERSION_RELEASE//' | sed 's/ //g'`

		AC_MSG_CHECKING(whether petsc is the development version)
		dnl if test x$PETSC_VERSION_DATE_HG = xunknown; then
		if test "$PETSC_RELEASE" = "0"; then
		   AC_DEFINE([_HAVE_PETSCDEV_],[1],[with PETSc-dev])
			AC_MSG_RESULT(yes)
		else
			AC_MSG_RESULT(no)
		fi
	
		AC_ARG_WITH([petsc-arch],
		  AS_HELP_STRING([--with-petsc-arch=DIR],[PETSc arch, necessary for PETSc < 3.0]),
		  [PETSC_ARCH=$withval],[PETSC_ARCH=""])

		AC_MSG_CHECKING(for petsc headers and libraries in $PETSC_ROOT)
		dnl To ge PETSc's libraries:
		dnl cd externalpackages/petsc/src
		dnl make getlinklibs
		PETSCINCL=" -I$PETSC_ROOT/include"
		dnl Add other location (not needed anymore since at least PETSc 3.0)
		if test "x$PETSC_ARCH" != "x" && test -d "$PETSC_ROOT/$PETSC_ARCH/include"; then
		 PETSCINCL+=" $PETSC_ROOT/$PETSC_ARCH/include"
		fi
		if test "x$PETSC_ARCH" != "x" && test -d "$PETSC_ROOT/include/$PETSC_ARCH"; then
		 PETSCINCL+=" $PETSC_ROOT/include/$PETSC_ARCH"
		fi
		
		case "${host_os}" in
				*cygwin*)
				if test $PETSC_MAJOR -lt 3 ; then
					PETSCLIB=-Wl,/LIBPATH:`cygpath -w $PETSC_ROOT/lib`  -Wl,libpetscksp.lib  -Wl,libpetscdm.lib  -Wl,libpetscmat.lib  -Wl,libpetscvec.lib  -Wl,libpetscsnes.lib  -Wl,libpetscts.lib  -Wl,libmpiuni.lib  -Wl,libpetsc.lib
				else
					PETSCLIB="/link -Wl,/LIBPATH:`cygpath -m $PETSC_ROOT/lib`  -Wl,libpetsc.lib"
					PETSCINCL="/I`cygpath -m $PETSC_ROOT/include`"
				fi
				;;
				*linux*)
				if test $PETSC_MAJOR -lt 3 ; then
					PETSCLIB="-L$PETSC_ROOT/lib -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetscsnes -lpetscts -lmpiuni -lpetsc"
				else
					PETSCLIB="-L$PETSC_ROOT/lib -lpetsc -ldl"
					if test $PETSC_MAJOR -gt 3 || test $PETSC_MINOR -ge 3; then PETSCLIB+=" -lparmetis -lmetis"; fi
				fi
				if test "x$host_os_version" = "x3.0.101-0.31.1_1.0502.8394-cray_gem_s" ; then
					PETSCLIB="-L$PETSC_ROOT/lib -lcraypetsc_gnu_real -lmetis"
				fi 
				;;
				*darwin*)
				if test $PETSC_MAJOR -lt 3 ; then
					PETSCLIB="-L$PETSC_ROOT/lib -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetscsnes -lpetscts -lpetsc"
				else
					PETSCLIB="-L$PETSC_ROOT/lib -lpetsc"
					if test $PETSC_MAJOR -gt 3 || test $PETSC_MINOR -ge 3; then PETSCLIB+=" -lmetis"; fi
				fi
				;;
		esac
		AC_MSG_RESULT(done)
		AC_DEFINE([_HAVE_PETSC_],[1],[with PETSc in ISSM src])
		AC_SUBST([PETSCINCL])
		AC_SUBST([PETSCLIB])
	fi
	dnl }}}
	dnl metis{{{
	if test "$HAVE_PETSC" = "yes" && test "x$PETSC_MAJOR" = "x3" && test $PETSC_MINOR -ge 3 && test "x$VENDOR" != "xMSVC-Win64" && test "x$VENDOR" != "xMSVC-Win64-par"; then
		dnl in petsc >=3.3, metis is provided
		HAVE_METIS="yes"
		AC_DEFINE([_METIS_VERSION_],[5],[ Metis version number])
		AC_DEFINE([_HAVE_METIS_],[1],[with Metis in ISSM src])
	else
		AC_ARG_WITH([metis-dir],
		  AS_HELP_STRING([--with-metis-dir=DIR],[metis root directory. necessary for serial build]),
		  [METIS_ROOT=$withval],[METIS_ROOT="no"])

		dnl Check whether metis is enabled
		AC_MSG_CHECKING([for metis])
		if test "x$METIS_ROOT" = "xno" ; then
			HAVE_METIS=no
		else
			HAVE_METIS=yes
			if ! test -d "$METIS_ROOT"; then
				AC_MSG_ERROR([metis directory provided ($METIS_ROOT) does not exist]);
			fi
		fi
		AC_MSG_RESULT($HAVE_METIS)
		AM_CONDITIONAL([METIS],[test x$HAVE_METIS = xyes])

		dnl library and header files
		if test "x$HAVE_METIS" = "xyes"; then

			AC_MSG_CHECKING(for metis headers and libraries in $METIS_ROOT)
			dnl first figure out version of metis: does the VERSION file exist?
			if test -e "$METIS_ROOT/VERSION"; then
				METIS_VERSION=4
			else
				METIS_VERSION=5
			fi

			if test "$METIS_VERSION" = "4" ; then
				METISINCL=-I"$METIS_ROOT/Lib" 
				case "${host_os}" in
					*cygwin*)
					METISINCL="/I`cygpath -m $METIS_ROOT/Lib`" 
					METISLIB="-Wl,/link -Wl,/LIBPATH:`cygpath -m $METIS_ROOT` -Wl,libmetis.lib"
					;;
					*linux*)
					METISLIB=-L"$METIS_ROOT/ -lmetis"
					;;
					*darwin*)
					METISLIB=-L"$METIS_ROOT/ -lmetis"
					;;
				esac
				AC_DEFINE([_METIS_VERSION_],[4],[ Metis version number])
			fi

			if test "$METIS_VERSION" = "5" ; then
				case "${host_os}" in
					*cygwin*)
					METISLIB="-L$METIS_ROOT libmetis.lib"
					;;
					*linux*)
					METISLIB=-L"$METIS_ROOT/lib -lmetis"
					;;
					*darwin*)
					METISLIB=-L"$METIS_ROOT/lib -lmetis"
					;;
				esac
				METISINCL=-I"$METIS_ROOT/include" 
				AC_DEFINE([_METIS_VERSION_],[5],[ Metis version number])
			fi

			AC_DEFINE([_HAVE_METIS_],[1],[with Metis in ISSM src])
			AC_SUBST([METISINCL])
			AC_SUBST([METISLIB])
		fi
	fi
	AM_CONDITIONAL([METIS],[test x$HAVE_METIS = xyes])
	dnl }}}
	dnl tao{{{
	AC_ARG_WITH([tao-dir],
		AS_HELP_STRING([--with-tao-dir=DIR], [tao root directory.]),
		[TAO_ROOT=$withval],[TAO_ROOT="no"]) 

	dnl Check whether tao is enabled
	AC_MSG_CHECKING([for tao])

	if test "$HAVE_PETSC" = "yes" && test "x$PETSC_MAJOR" = "x3" && test $PETSC_MINOR -ge 5; then
		dnl in petsc >=3.5, tao is provided
		HAVE_TAO="yes"
		AC_DEFINE([_HAVE_TAO_],[1],[with Tao in ISSM src])
		AC_MSG_RESULT($HAVE_TAO)
	else

		if test "x$TAO_ROOT" = "xno" ; then
			HAVE_TAO=no
		else
			HAVE_TAO=yes
			if ! test -d "$TAO_ROOT"; then
				AC_MSG_ERROR([tao directory provided ($TAO_ROOT) does not exist]);
			fi
		fi
		AC_MSG_RESULT($HAVE_TAO)

		dnl tao headers and libraries
		if test "x$HAVE_TAO" == "xyes"; then
		  TAOINCL="-I$TAO_ROOT/ -I$TAO_ROOT/include -I$TAO_ROOT/bmake/ "
		  TAOLIB="-L$TAO_ROOT/lib -ltao -lpetsc"
		  AC_DEFINE([_HAVE_TAO_],[1],[with Tao in ISSM src])
		  AC_SUBST([TAOINCL])
		  AC_SUBST([TAOLIB])
		fi
	fi
	dnl }}}
	dnl m1qn3{{{
	AC_ARG_WITH([m1qn3-dir],
		AS_HELP_STRING([--with-m1qn3-dir=DIR], [m1qn3 root directory.]),
		[M1QN3_ROOT=$withval],[M1QN3_ROOT="no"]) 

	dnl Check whether m1qn3 is enabled
	AC_MSG_CHECKING([for m1qn3])
	if test "x$M1QN3_ROOT" = "xno" ; then
		HAVE_M1QN3=no
	else
		HAVE_M1QN3=yes
		if ! test -d "$M1QN3_ROOT"; then
			AC_MSG_ERROR([m1qn3 directory provided ($M1QN3_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_M1QN3)
	
	dnl m1qn3 headers and libraries
	if test "x$HAVE_M1QN3" == "xyes"; then
	  M1QN3LIB="$M1QN3_ROOT/libm1qn3.a $M1QN3_ROOT//libddot.a"
	  AC_DEFINE([_HAVE_M1QN3_],[1],[with M1QN3 in ISSM src])
	  AC_SUBST([M1QN3LIB])
	fi
	dnl }}}
	dnl proj.4{{{
	AC_ARG_WITH([proj4-dir],
		AS_HELP_STRING([--with-proj4-dir=DIR], [proj4 root directory.]),
		[PROJ4_ROOT=$withval],[PROJ4_ROOT="no"]) 

	dnl Check whether proj4 is enabled
	AC_MSG_CHECKING([for proj.4])
	if test "x$PROJ4_ROOT" = "xno" ; then
		HAVE_PROJ4=no
	else
		HAVE_PROJ4=yes
		if ! test -d "$PROJ4_ROOT"; then
			AC_MSG_ERROR([proj4 directory provided ($PROJ4_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_PROJ4)
	
	dnl proj4 headers and libraries
	if test "x$HAVE_PROJ4" == "xyes"; then
	  PROJ4INCL="-I$PROJ4_ROOT/include "
	  PROJ4LIB="-L$PROJ4_ROOT/lib -lproj"
	  AC_DEFINE([_HAVE_PROJ4_],[1],[with PROJ4 in ISSM src])
	  AC_SUBST([PROJ4INCL])
	  AC_SUBST([PROJ4LIB])
	fi
	AM_CONDITIONAL([PROJ4],[test x$HAVE_PROJ4 = xyes])
	dnl }}}
	dnl slepc{{{
	AC_ARG_WITH([slepc-dir],
	  AS_HELP_STRING([--with-slepc-dir=DIR],[slepc root directory]),
	  [SLEPC_ROOT=$withval],[SLEPC_ROOT="no"])
			  
	dnl Check whether slepc is enabled
	AC_MSG_CHECKING([for slepc])
	if test "x$SLEPC_ROOT" = "xno" ; then
		HAVE_SLEPC=no
	else
		HAVE_SLEPC=yes
		if ! test -d "$SLEPC_ROOT"; then
			AC_MSG_ERROR([slepc directory provided ($SLEPC_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_SLEPC)
	
	dnl slepc headers and libraries
	if test "x$HAVE_SLEPC" == "xyes"; then
		SLEPCINCL=-I"$SLEPC_ROOT/include"
		SLEPCLIB=-L"$SLEPC_ROOT/lib/ -lslepc"
		AC_DEFINE([_HAVE_SLEPC_],[1],[with Slepc in ISSM src])
		AC_SUBST([SLEPCINCL])
		AC_SUBST([SLEPCLIB])
	fi
	dnl }}}
	dnl shapelib{{{
	AC_ARG_WITH([shapelib-dir],
	  AS_HELP_STRING([--with-shapelib-dir=DIR], [shapelib root directory]),
	  [SHAPELIB_ROOT=$withval],[SHAPELIB_ROOT="no"])
			  
	dnl Check whether shapelib is enabled
	AC_MSG_CHECKING([for shapelib])
	if test "x$SHAPELIB_ROOT" = "xno" ; then
		HAVE_SHAPELIB=no
	else
		HAVE_SHAPELIB=yes
		if ! test -d "$SHAPELIB_ROOT"; then
			AC_MSG_ERROR([shapelib directory provided ($SHAPELIB_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_SHAPELIB)
	
	dnl shapelib headers and libraries
	if test "x$HAVE_SHAPELIB" == "xyes"; then
		SHAPELIBINCL=-I"$SHAPELIB_ROOT/include"
		SHAPELIBLIB=-L"$SHAPELIB_ROOT/lib/ -lshape"
		AC_DEFINE([_HAVE_SHAPELIB_],[1],[with Shapelib in ISSM src])
		AC_SUBST([SHAPELIBINCL])
		AC_SUBST([SHAPELIBLIB])
	fi
	dnl }}}
	dnl scalapack{{{
	AC_ARG_WITH([scalapack-dir],
	  AS_HELP_STRING([--with-scalapack-dir=DIR],[scalapack root directory]),
	  [SCALAPACK_ROOT=$withval],[SCALAPACK_ROOT="no"])
			  
	dnl Check whether scalapack is enabled
	AC_MSG_CHECKING([for scalapack])
	if test "x$SCALAPACK_ROOT" = "xno" ; then
		HAVE_SCALAPACK=no
	else
		HAVE_SCALAPACK=yes
		if ! test -d "$SCALAPACK_ROOT"; then
			AC_MSG_ERROR([scalapack directory provided ($SCALAPACK_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_SCALAPACK)
	
	dnl scalapack headers and libraries
	if test "x$HAVE_SCALAPACK" == "xyes"; then
		if test x$VENDOR = xintel-discover; then
		 SCALAPACKLIB=-L"$SCALAPACK_ROOT/ -lmkl_scalapack_lp64"
		else
		 SCALAPACKLIB=-L"$SCALAPACK_ROOT/ -lscalapack"
		fi
		AC_DEFINE([_HAVE_SCALAPACK_],[1],[with Scalapack in ISSM src])
		AC_SUBST([SCALAPACKLIB])
	fi
	dnl }}}
	dnl blas-lapack{{{
	AC_ARG_WITH([blas-lapack-dir],
	  AS_HELP_STRING([--with-blas-lapack-dir=DIR],[blas-lapack root directory]),
	  [BLASLAPACK_ROOT=$withval],[BLASLAPACK_ROOT="no"])
			  
	dnl Check whether blas-lapack is enabled
	AC_MSG_CHECKING([for blas-lapack])
	if test "x$BLASLAPACK_ROOT" = "xno" ; then
		HAVE_BLASLAPACK=no
	else
		HAVE_BLASLAPACK=yes
		if ! test -d "$BLASLAPACK_ROOT"; then
			AC_MSG_ERROR([blas-lapack directory provided ($BLASLAPACK_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_BLASLAPACK)
	
	dnl blas-lapack headers and libraries
	if test "x$HAVE_BLASLAPACK" == "xyes"; then
		BLASLAPACKINCL=""
		if test x$VENDOR = xintel-discover; then
		   dnl works for intel 11
			dnl BLASLAPACKLIB=-L"$BLASLAPACK_ROOT -lmkl_lapack -lmkl -lguide -lpthread "
			dnl dnl works for intel 13
			BLASLAPACKLIB=-L"$BLASLAPACK_ROOT -lmkl_lapack95_lp64 -lmkl_rt " 
		else
			dnl: branch on whether we are running on windows or linux.
			case "${host_os}" in
				*cygwin*)
				BLASLAPACKLIB="-L`cygpath -m $BLASLAPACK_ROOT` -Wl,libf2cblas.lib  -Wl,libf2clapack.lib"
				;;
				*linux*)
				BLASLAPACKLIB=-L"$BLASLAPACK_ROOT/lib -lflapack -lfblas " 
				;;
				*darwin*)
				BLASLAPACKLIB=-L"$BLASLAPACK_ROOT/lib -lflapack -lfblas " 
				;;
			esac
		fi
		AC_DEFINE([_HAVE_BLASLAPACK_],[1],[with blas lapack in ISSM src])
		AC_SUBST([BLASLAPACKLIB])
		AC_SUBST([BLASLAPACKINCL])
	fi
	dnl }}}
	dnl mkl{{{
		AC_ARG_WITH([mkl-libflags],
					AS_HELP_STRING([--with-mkl-libflags=LIBS],[mkl libraries to be used]),
					[MKL_LIBFLAGS=$withval],[MKL_LIBFLAGS="no"])

		  dnl Check whether mkl is enabled
		  AC_MSG_CHECKING([for mkl])
		  if test "x$MKL_LIBFLAGS" = "xno" ; then
				HAVE_MKL=no
		  else
				HAVE_MKL=yes
				MKLLIB="$MKL_LIBFLAGS"
				AC_DEFINE([_HAVE_MKL_],[1],[with mkl in ISSM src])
				AC_SUBST([MKLLIB])
				AC_SUBST([MKLINCL])
			fi
			AC_MSG_RESULT($HAVE_MKL)
	dnl }}}
	dnl plapack{{{
	AC_MSG_CHECKING(for plapack)
	
	AC_ARG_WITH([plapack-lib],
	  AS_HELP_STRING([--with-plapack-lib = lib],[plapack library]),
	  [PLAPACK_LIB=$withval],[PLAPACK_LIB=""])
	
	AC_ARG_WITH([plapack-include],
			  AS_HELP_STRING([--with-plapack-include = include],
							 [plapack include ]),
			  [PLAPACK_INCLUDE=$withval],[PLAPACK_INCLUDE=""])
	  
	if test -n "$PLAPACK_LIB"; then
		if test -n "$PLAPACK_INCLUDE"; then
		
			HAVE_PLAPACK=yes
			PLAPACKINCL="$PLAPACK_INCLUDE"
			PLAPACKLIB="$PLAPACK_LIB"
			AC_DEFINE([_HAVE_PLAPACK_],[1],[with Plapack in ISSM src])
			AC_SUBST([PLAPACKINCL])
			AC_SUBST([PLAPACKLIB])
		else
			HAVE_PLAPACK=no
		fi
	else
		HAVE_PLAPACK=no
	fi
	AC_MSG_RESULT($HAVE_PLAPACK)
	dnl }}}
	dnl mumps{{{
	AC_ARG_WITH([mumps-dir],
	  AS_HELP_STRING([--with-mumps-dir=DIR],[mumps root directory]),
	  [MUMPS_ROOT=$withval],[MUMPS_ROOT="no"])
			  
	dnl Check whether mumps is enabled
	AC_MSG_CHECKING([for mumps])
	if test "x$MUMPS_ROOT" = "xno" ; then
		HAVE_MUMPS=no
	else
		HAVE_MUMPS=yes
		if ! test -d "$MUMPS_ROOT"; then
			AC_MSG_ERROR([mumps directory provided ($MUMPS_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_MUMPS)
	
	dnl mumps headers and libraries
	if test "x$HAVE_MUMPS" == "xyes"; then
		MUMPSINCL=-I"$MUMPS_ROOT/include"
		if test "$PETSC_MAJOR" = "2" ; then
			MUMPSLIB=-L"$MUMPS_ROOT/lib "
		else
			MUMPSLIB=-L"$MUMPS_ROOT/lib -ldmumps -lcmumps  -lmumps_common -lpord -lparmetis -lzmumps -lmetis"
			dnl MUMPSLIB=-L"$MUMPS_ROOT/lib "
		fi
		AC_DEFINE([_HAVE_MUMPS_],[1],[with Mumps in ISSM src])
		AC_SUBST([MUMPSINCL])
		AC_SUBST([MUMPSLIB])
	fi
	AM_CONDITIONAL([MUMPS], [test x$HAVE_MUMPS = xyes])
	dnl }}}
	dnl blacs{{{
	AC_ARG_WITH([blacs-dir],
		AS_HELP_STRING([--with-blacs-dir=DIR],[blacs root directory]),
			  [BLACS_ROOT=$withval],[BLACS_ROOT="no"])
			  
	dnl Check whether blacs is enabled
	AC_MSG_CHECKING([for blacs])
	if test "x$BLACS_ROOT" = "xno" ; then
		HAVE_BLACS=no
	else
		HAVE_BLACS=yes
		if ! test -d "$BLACS_ROOT"; then
			AC_MSG_ERROR([blacs directory provided ($BLACS_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_BLACS)
	
	dnl blacs headers and libraries
	if test "x$HAVE_BLACS" == "xyes"; then
		BLACSINCL=""
		if test x$VENDOR = xintel-discover; then
		 BLACSLIB=-L"$BLACS_ROOT/ -lmkl_blacs_intelmpi_lp64"
		else
		 BLACSLIB=-L"$BLACS_ROOT/ -lblacs"
		fi
		AC_DEFINE([_HAVE_BLACS_],[1],[with Blacs in ISSM src])
		AC_SUBST([BLACSINCL])
		AC_SUBST([BLACSLIB])
	fi
	dnl }}}
	dnl hypre{{{
	AC_ARG_WITH([hypre-dir],
	  AS_HELP_STRING([--with-hypre-dir=DIR],[hypre root directory]),
			  [HYPRE_ROOT=$withval],[HYPRE_ROOT="no"])
			  
	dnl Check whether hypre is enabled
	AC_MSG_CHECKING([for hypre])
	if test "x$HYPRE_ROOT" = "xno" ; then
		HAVE_HYPRE=no
	else
		HAVE_HYPRE=yes
		if ! test -d "$HYPRE_ROOT"; then
			AC_MSG_ERROR([hypre directory provided ($HYPRE_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_HYPRE)

	dnl hypre headers and libraries
	if test "x$HAVE_HYPRE" == "xyes"; then
		HYPREINCL=""
		HYPRELIB=-L"$HYPRE_ROOT/lib -lHYPRE"
		AC_DEFINE([_HAVE_HYPRE_],[1],[with Hypre in ISSM src])
		AC_SUBST([HYPREINCL])
		AC_SUBST([HYPRELIB])
	fi
	dnl }}}
	dnl prometheus{{{
	AC_ARG_WITH([prometheus-dir],
				AS_HELP_STRING([--with-prometheus-dir=DIR],[prometheus root directory]),
				[PROMETHEUS_ROOT=$withval],[PROMETHEUS_ROOT="no"])

		dnl Check whether prometheus is enabled
		AC_MSG_CHECKING([for prometheus])
		if test "x$PROMETHEUS_ROOT" = "xno" ; then
			HAVE_PROMETHEUS=no
		else
			HAVE_PROMETHEUS=yes
			if ! test -d "$PROMETHEUS_ROOT"; then
				AC_MSG_ERROR([prometheus directory provided ($PROMETHEUS_ROOT) does not exist]);
			fi
		fi
		AC_MSG_RESULT($HAVE_PROMETHEUS)

		dnl prometheus headers and libraries
		if test "x$HAVE_PROMETHEUS" == "xyes"; then
			 PROMETHEUSINCL=-I"$PROMETHEUS_ROOT/include"
			 PROMETHEUSLIB=-L"$PROMETHEUS_ROOT/lib -lpromfei -lprometheus -lparmetis"
			 AC_DEFINE([_HAVE_PROMETHEUS_],[1],[with Prometheus in ISSM src])
			 AC_SUBST([PROMETHEUSINCL])
			 AC_SUBST([PROMETHEUSLIB])
	   fi
		dnl }}}
dnl spai{{{
	AC_ARG_WITH([spai-dir],
				AS_HELP_STRING([--with-spai-dir=DIR],[spai root directory]),
				[SPAI_ROOT=$withval],[SPAI_ROOT="no"])

		dnl Check whether spai is enabled
		AC_MSG_CHECKING([for spai])
		if test "x$SPAI_ROOT" = "xno" ; then
			HAVE_SPAI=no
		else
			HAVE_SPAI=yes
			if ! test -d "$SPAI_ROOT"; then
				AC_MSG_ERROR([spai directory provided ($SPAI_ROOT) does not exist]);
			fi
		fi
		AC_MSG_RESULT($HAVE_SPAI)

		dnl spai headers and libraries
		if test "x$HAVE_SPAI" == "xyes"; then
			SPAIINCL=-I"$SPAI_ROOT/include"
			SPAILIB=-L"$SPAI_ROOT/lib -lspai"
			AC_DEFINE([_HAVE_SPAI_],[1],[with Spai in ISSM src])
			AC_SUBST([SPAIINCL])
			AC_SUBST([SPAILIB])
		fi
	  dnl }}}
dnl superlu{{{ 
	AC_ARG_WITH([superlu-dir],
				AS_HELP_STRING([--with-superlu-dir=DIR],[superlu root directory]),
				[SUPERLU_ROOT=$withval],[SUPERLU_ROOT="no"])

	dnl Check whether superlu is enabled
	AC_MSG_CHECKING([for superlu])
	if test "x$SUPERLU_ROOT" = "xno" ; then
		HAVE_SUPERLU=no
	else
		HAVE_SUPERLU=yes
		if ! test -d "$SUPERLU_ROOT"; then
			AC_MSG_ERROR([superlu directory provided ($SUPERLU_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_SUPERLU)
	
	dnl superlu headers and libraries
	if test "x$HAVE_SUPERLU" == "xyes"; then
		  SUPERLUINCL=-I"$SUPERLU_ROOT/include"
		  SUPERLULIB=-L"$SUPERLU_ROOT/lib -lsuperlu_dist"
		  AC_DEFINE([_HAVE_SUPERLU_],[1],[with Superlu in ISSM src])
		  AC_SUBST([SUPERLUINCL])
		  AC_SUBST([SUPERLULIB])
	 fi
	 dnl }}}
dnl spooles{{{ 
	AC_ARG_WITH([spooles-dir],
				AS_HELP_STRING([--with-spooles-dir=DIR],[spooles root directory]),
				[SPOOLES_ROOT=$withval],[SPOOLES_ROOT="no"])

	dnl Check whether spooles is enabled
	AC_MSG_CHECKING([for spooles])
	if test "x$SPOOLES_ROOT" = "xno" ; then
		HAVE_SPOOLES=no
	else
		HAVE_SPOOLES=yes
		if ! test -d "$SPOOLES_ROOT"; then
			AC_MSG_ERROR([spooles directory provided ($SPOOLES_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_SPOOLES)
	
	dnl spooles headers and libraries
	if test "x$HAVE_SPOOLES" == "xyes"; then
		  SPOOLESINCL=-I"$SPOOLES_ROOT/include"
		  SPOOLESLIB=-L"$SPOOLES_ROOT/lib -lspooles"
		  AC_DEFINE([_HAVE_SPOOLES_],[1],[with Spooles in ISSM src])
		  AC_SUBST([SPOOLESINCL])
		  AC_SUBST([SPOOLESLIB])
	 fi
	 dnl }}}
dnl pastix{{{ 
	AC_ARG_WITH([pastix-dir],
				AS_HELP_STRING([--with-pastix-dir=DIR],[pastix root directory]),
				[PASTIX_ROOT=$withval],[PASTIX_ROOT="no"])

	dnl Check whether pastix is enabled
	AC_MSG_CHECKING([for pastix])
	if test "x$PASTIX_ROOT" = "xno" ; then
		HAVE_PASTIX=no
	else
		HAVE_PASTIX=yes
		if ! test -d "$PASTIX_ROOT"; then
			AC_MSG_ERROR([pastix directory provided ($PASTIX_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_PASTIX)
	
	dnl pastix headers and libraries
	if test "x$HAVE_PASTIX" == "xyes"; then
		  PASTIXINCL=-I"$PASTIX_ROOT/include"
		  PASTIXLIB=-L"$PASTIX_ROOT/lib -lpastix_XXbit_mpi_smp_nobubble_int32_simple_real_scotch_i686_pc_linux -lptscotch -lptscotcherr -lpastix"
		  AC_DEFINE([_HAVE_PASTIX_],[1],[with Pastix in ISSM src])
		  AC_SUBST([PASTIXINCL])
		  AC_SUBST([PASTIXLIB])
  fi
  dnl }}}
	dnl ml{{{
	AC_ARG_WITH([ml-dir],
	  AS_HELP_STRING([--with-ml-dir=DIR],[ml root directory]),
			  [ML_ROOT=$withval],[ML_ROOT="no"])
			  
	dnl Check whether ml is enabled
	AC_MSG_CHECKING([for ml])
	if test "x$ML_ROOT" = "xno" ; then
		HAVE_ML=no
	else
		HAVE_ML=yes
		if ! test -d "$ML_ROOT"; then
			AC_MSG_ERROR([ml directory provided ($ML_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_ML)
	
	dnl ml headers and libraries
	if test "x$HAVE_ML" == "xyes"; then
		MLINCL=-I"$ML_ROOT/include"
		MLLIB=-L"$ML_ROOT/lib -lml"
		AC_DEFINE([_HAVE_ML_],[1],[with Blacs in ISSM src])
		AC_SUBST([MLINCL])
		AC_SUBST([MLLIB])
	fi
	dnl }}}
	dnl umfpack{{{
		AC_ARG_WITH([umfpack-dir],
		  AS_HELP_STRING([--with-umfpack-dir=DIR],[UMFPACK root directory]),
					[UMFPACK_ROOT=$withval],[UMFPACK_ROOT="no"])
			  
	dnl Check whether umfpack is enabled
	AC_MSG_CHECKING([for umfpack])
	if test "x$UMFPACK_ROOT" = "xno" ; then
		HAVE_UMFPACK=no
	else
		HAVE_UMFPACK=yes
		if ! test -d "$UMFPACK_ROOT"; then
			AC_MSG_ERROR([umfpack directory provided ($UMFPACK_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_UMFPACK)
	
	dnl umfpack headers and libraries
	if test "x$HAVE_UMFPACK" == "xyes"; then
		UMFPACKINCL=""
		UMFPACKLIB=-L"$UMFPACK_ROOT/lib -lumfpack -lumfpack.5.5.1"
		AC_DEFINE([_HAVE_UMFPACK_],[1],[with UMFPACK in ISSM src])
		AC_SUBST([UMFPACKINCL])
		AC_SUBST([UMFPACKLIB])
	fi
	dnl }}}
dnl math{{{
	AC_MSG_CHECKING(for math library)
	AC_ARG_WITH([math-lib],
	  AS_HELP_STRING([--with-math-lib = otions],[math options, for ex: "/usr/lib/libm.a]),
	  [MATH_LIB=$withval],[MATH_LIB=""])

	dnl check that --with-math-lib may have been provided
	if test -n "$MATH_LIB" ; then
		HAVE_MATH=yes
		MATHLIB="$MATH_LIB"
		AC_DEFINE([_HAVE_MATH_],[1],[with MATH in ISSM src])
		AC_SUBST([MATHLIB])
	fi
	AC_MSG_RESULT(done)
	dnl }}}
	dnl math77{{{
		AC_ARG_WITH([math77-dir],
					AS_HELP_STRING([--with-math77-dir=DIR], [math77 root directory.]),
					[MATH77_ROOT=$withval],[MATH77_ROOT="no"]) 
		  
	dnl Check whether math77 is enabled
	AC_MSG_CHECKING([for math77])
	if test "x$MATH77_ROOT" = "xno" ; then
		HAVE_MATH77=no
	else
		HAVE_MATH77=yes
		if ! test -d "$MATH77_ROOT"; then
			AC_MSG_ERROR([math77 directory provided ($MATH77_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_MATH77)
	
	dnl math77 headers and libraries
	if test "x$HAVE_MATH77" == "xyes"; then
		MATH77LIB="-L$MATH77_ROOT/ -lmath77"
		AC_DEFINE([_HAVE_MATH77_],[1],[with math77 in ISSM src])
		AC_SUBST([MATH77LIB])
   fi
	dnl }}}
	dnl fortran{{{
	AC_ARG_WITH([fortran],
		AS_HELP_STRING([--with-fortran = YES], [do we compile fortran code (default is yes)]),
		[FORTRAN=$withval],[FORTRAN=yes]) 
	AC_MSG_CHECKING(for fortran compilation)
	if test "x$FORTRAN" = "xyes"; then
		dnl defaults
		HAVE_FORTRAN=yes
		AC_DEFINE([_HAVE_FORTRAN_],[1],[with fortran capability])
	else
		HAVE_FORTRAN=no
	fi
	AM_CONDITIONAL([FORTRAN], [test x$FORTRAN = xyes])
	AC_MSG_RESULT($FORTRAN)

	if test "x$FORTRAN" = "xyes"; then
		dnl fortran library  option
		AC_MSG_CHECKING(for fortran library)
		AC_ARG_WITH([fortran-lib],
		  AS_HELP_STRING([--with-fortran-lib = options],[fortran options, for ex: "/usr/lib/gfortran.a]),
			[FORTRAN_LIB=$withval],[FORTRAN_LIB=""])

		dnl check that --with-fortran-lib may have been provided
		if test -n "$FORTRAN_LIB" ; then
			dnl check that library provided EXISTS!
		   FORTRAN_DIR=$(echo $FORTRAN_LIB | sed -e "s/-L//g" | awk '{print $[1]}')
			if test -d "$FORTRAN_DIR" || test -f "$FORTRAN_DIR"; then
				FORTRANLIB="$FORTRAN_LIB"
				AC_DEFINE([_HAVE_FORTRAN_],[1],[with FORTRAN in ISSM src])
				AC_SUBST([FORTRANLIB])
			else
			 if test "x$HAVE_MPI" = "xyes"; then
				FORTRANLIB=$(mpif77 -print-file-name="libgfortran.a")
				if test -f "$FORTRANLIB"; then
					 AC_MSG_ERROR([fortran library provided ($FORTRAN_LIB) does not exist, MPI suggests the following library: $FORTRANLIB]);
				fi
			 fi
				AC_MSG_ERROR([fortran library provided ($FORTRAN_LIB) does not exist!]);
			fi
		fi
		AC_MSG_RESULT(done)
	fi
	dnl }}}
	dnl graphics{{{
	AC_MSG_CHECKING(for graphics library)
	AC_ARG_WITH([graphics-lib],
	  AS_HELP_STRING([--with-graphics-lib = options],[graphics options, for ex: "/usr/X11/lib/libX11.a]),
	  [GRAPHICS_LIB=$withval],[GRAPHICS_LIB=""])

	dnl check that --with-graphics-lib may have been provided
	if test -n "$GRAPHICS_LIB" ; then
		dnl check that library provided EXISTS!
		GRAPHICS_DIR=$(echo $GRAPHICS_LIB | sed -e "s/-L//g" | awk '{print $[1]}')
		if test -d "$GRAPHICS_DIR" || test -f "$GRAPHICS_DIR"; then
			HAVE_GRAPHICS=yes
			GRAPHICSLIB="$GRAPHICS_LIB"
			AC_DEFINE([_HAVE_GRAPHICS_],[1],[with GRAPHICS in ISSM src])
			AC_SUBST([GRAPHICSLIB])
		else
			if test -f "$PETSC_ROOT/conf/petscvariables"; then
				GRAPHICSLIB=$(cat $PETSC_ROOT/conf/petscvariables | grep X_LIB)
				AC_MSG_ERROR([graphics library provided ($GRAPHICS_LIB) does not exist, PETSc suggests the following library: $GRAPHICSLIB]);
			fi
			AC_MSG_ERROR([graphics library provided ($GRAPHICS_LIB$) does not exist!]);
		fi
	fi
	AC_MSG_RESULT(done)
	dnl }}}
	dnl meteoio{{{
	AC_ARG_WITH([meteoio-dir],
	  AS_HELP_STRING([--with-meteoio-dir=DIR], [use meteoio in conjunction with snowpack model.]),
	  [METEOIO_ROOT=$withval],[METEOIO_ROOT="no"]) 

	dnl Check whether meteoio is enabled
	AC_MSG_CHECKING([for meteoio])
	if test "x$METEOIO_ROOT" = "xno" ; then
		HAVE_METEOIO=no
	else
		HAVE_METEOIO=yes
		if ! test -d "$METEOIO_ROOT"; then
			AC_MSG_ERROR([meteoio directory provided ($METEOIO_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_METEOIO)
	
	dnl meteoio headers and libraries
	if test "x$HAVE_METEOIO" == "xyes"; then
		METEOIOINCL="-I$METEOIO_ROOT/include"
		METEOIOLIB="-dy -L$METEOIO_ROOT/lib  -lmeteoio "

		AC_DEFINE([_HAVE_METEOIO_],[1],[with meteoio])
		AC_SUBST([METEOIOINCL])
		AC_SUBST([METEOIOLIB])
	fi
	AM_CONDITIONAL([METEOIO], [test x$HAVE_METEOIO = xyes])
	dnl }}}
	dnl snowpack{{{
	AC_ARG_WITH([snowpack-dir],
	  AS_HELP_STRING([--with-snowpack-dir=DIR], [use snowpack for surface mass balance model.]),
	  [SNOWPACK_ROOT=$withval],[SNOWPACK_ROOT="no"]) 

	dnl Check whether snowpack is enabled
	AC_MSG_CHECKING([for snowpack])
	if test "x$SNOWPACK_ROOT" = "xno" ; then
		HAVE_SNOWPACK=no
	else
		HAVE_SNOWPACK=yes
		if ! test -d "$SNOWPACK_ROOT"; then
			AC_MSG_ERROR([snowpack directory provided ($SNOWPACK_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_SNOWPACK)
	
	dnl snowpack headers and libraries
	if test "x$HAVE_SNOWPACK" == "xyes"; then
		SNOWPACKINCL="-I$SNOWPACK_ROOT/include"
		SNOWPACKLIB="-dy -L$SNOWPACK_ROOT/lib  -lsnowpack "

		AC_DEFINE([_HAVE_SNOWPACK_],[1],[with snowpack for surface mass balance model])
		AC_SUBST([SNOWPACKINCL])
		AC_SUBST([SNOWPACKLIB])
	fi
	AM_CONDITIONAL([SNOWPACK], [test x$HAVE_SNOWPACK = xyes])
	dnl }}}
	dnl neopz{{{
	AC_ARG_WITH([neopz-dir],
		AS_HELP_STRING([--with-neopz-dir=DIR], [neopz root directory.]),
		[NEOPZ_ROOT=$withval],[NEOPZ_ROOT="no"]) 

	dnl Check whether neopz is enabled
	AC_MSG_CHECKING([for neopz])
	if test "x$NEOPZ_ROOT" = "xno" ; then
		HAVE_NEOPZ=no
	else
		HAVE_NEOPZ=yes
		if ! test -d "$NEOPZ_ROOT"; then
			AC_MSG_ERROR([neopz directory provided ($NEOPZ_ROOT) does not exist]);
		fi
	fi
	AC_MSG_RESULT($HAVE_NEOPZ)
	
	dnl neopz headers and libraries
	if test "x$HAVE_NEOPZ" == "xyes"; then
	  NEOPZLIB="$NEOPZ_ROOT/lib/libpz.a"
     NEOPZINCL=-I"$NEOPZ_ROOT/include"
	  AC_DEFINE([_HAVE_NEOPZ_],[1],[with NEOPZ in ISSM src])
	  AC_SUBST([NEOPZINCL])
	  AC_SUBST([NEOPZLIB])
	fi
	AM_CONDITIONAL([NEOPZ], [test x$HAVE_NEOPZ = xyes])
	dnl }}}

	dnl Capabilities
	dnl with-bamg{{{
	AC_ARG_WITH([bamg],
		AS_HELP_STRING([--with-bamg = YES],[compile with bamg capabilities (default is yes)]),
		[BAMG=$withval],[BAMG=yes]) 
	AC_MSG_CHECKING(for bamg capability compilation)

	HAVE_BAMG=no
	if test "x$BAMG" = "xyes"; then
		HAVE_BAMG=yes
		AC_DEFINE([_HAVE_BAMG_],[1],[with bamg meshing capability])
	fi
	AM_CONDITIONAL([BAMG], [test x$HAVE_BAMG = xyes])
	AC_MSG_RESULT($HAVE_BAMG)
	dnl }}}
	dnl with-ocean{{{
	AC_ARG_WITH([ocean],
		AS_HELP_STRING([--with-ocean = YES],[compile with ice/ocean coupling (default is no)]),
		[OCEAN=$withval],[OCEAN=no]) 
	AC_MSG_CHECKING(for ice/ocean capability compilation)

	HAVE_OCEAN=no
	if test "x$OCEAN" = "xyes"; then
		HAVE_OCEAN=yes
		AC_DEFINE([_HAVE_OCEAN_],[1],[with ice/ocean coupling capability])
	fi
	AM_CONDITIONAL([OCEAN], [test x$HAVE_OCEAN = xyes])
	AC_MSG_RESULT($HAVE_OCEAN)
	dnl }}}
	dnl with-kml{{{
	AC_ARG_WITH([kml],
		AS_HELP_STRING([--with-kml = YES],[compile with kml capabilities (default is yes)]),
		[KML=$withval],[KML=yes]) 
	AC_MSG_CHECKING(for kml capability compilation)

	HAVE_KML=no
	if test "x$KML" = "xyes"; then
		HAVE_KML=yes
		AC_DEFINE([_HAVE_KML_],[1],[with kml capability])
	fi
	AM_CONDITIONAL([KML], [test x$HAVE_KML = xyes])
	AC_MSG_RESULT($HAVE_KML)
	dnl }}}
	dnl with-kriging{{{
	AC_ARG_WITH([kriging],
		AS_HELP_STRING([--with-kriging = YES],[compile with kriging capabilities (default is yes)]),
		[KRIGING=$withval],[KRIGING=yes]) 
	AC_MSG_CHECKING(for kriging capability compilation)

	HAVE_KRIGING=no
	if test "x$KRIGING" = "xyes"; then
		HAVE_KRIGING=yes
		AC_DEFINE([_HAVE_KRIGING_],[1],[with kriging capability])
	fi
	AM_CONDITIONAL([KRIGING], [test x$HAVE_KRIGING = xyes])
	AC_MSG_RESULT($HAVE_KRIGING)
	dnl }}}
	AX_ANALYSES_SELECTION

	dnl Platform specifics
	dnl with-ios{{{
	AC_ARG_WITH([ios],
		AS_HELP_STRING([--with-ios = YES], [compile with iOS capabilities (default is no, alternatives are yes)]),
		[IOS=$withval],[IOS=no]) 
	AC_MSG_CHECKING(for iOS compilation)

	if test "x$IOS" = "xyes"; then
		dnl defaults
		HAVE_IOS=yes

		AC_DEFINE([_HAVE_IOS_],[1],[with android capability])
	elif test "x$IOS" = "xno"; then
		HAVE_IOS=no
	else
	  AC_MSG_ERROR([--with-ios should be either no or yes])
	fi
	AM_CONDITIONAL([IOS], [test x$HAVE_IOS != xno])
	AC_MSG_RESULT($HAVE_IOS)
	dnl }}}
	dnl with-android{{{
	AC_ARG_WITH([android],
		AS_HELP_STRING([--with-android = EXE], [compile with android capabilities (default is no, alternatives are exe and jni)]),
		[ANDROID=$withval],[ANDROID=no]) 
	AC_MSG_CHECKING(for android capability compilation)

	if test "x$ANDROID" = "xjni"; then

		dnl defaults
		HAVE_ANDROID=jni
		AC_DEFINE([_HAVE_ANDROID_],[1],[with android capability])
		AC_DEFINE([_HAVE_ANDROID_JNI_],[1],[with android jni])
	elif test "x$ANDROID" = "xexe"; then
		dnl defaults
		HAVE_ANDROID=exe

		AC_DEFINE([_HAVE_ANDROID_],[1],[with android capability])
	elif test "x$ANDROID" = "xno"; then
		HAVE_ANDROID=no
	else
	  AC_MSG_ERROR([--with-android should be either no, exe or jni])
	fi
	AM_CONDITIONAL([ANDROID], [test x$HAVE_ANDROID != xno])
	AM_CONDITIONAL([ANDROIDJNI], [test x$HAVE_ANDROID = xjni])
	AM_CONDITIONAL([ANDROIDEXE], [test x$HAVE_ANDROID = xexe])
	AC_MSG_RESULT($HAVE_ANDROID)
	dnl }}}
	dnl with-android-ndk{{{
	AC_ARG_WITH([android-ndk],
	  AS_HELP_STRING([--with-android-ndk=DIR], [android-ndk root directory.]),
	  [ANDROID_NDK_ROOT=$withval],[ANDROID_NDK_ROOT=""]) 
	AC_MSG_CHECKING(with android ndk)
	
	if test -d "$ANDROID_NDK_ROOT"; then
		dnl defaults
		HAVE_ANDROID_NDK=yes
		ANDROID_NDKINCL="-I$ANDROID_NDK_ROOT/arm-linux-android-install/sysroot/usr/include"

		AC_DEFINE([_HAVE_ANDROID_NDK_],[1],[with android ndk in ISSM src])
		AC_SUBST([ANDROID_NDKINCL])
	else
		HAVE_ANDROID_NDK=no
	fi
	AC_MSG_RESULT($HAVE_ANDROID_NDK)
	dnl }}}

	dnl other options
	dnl optimization{{{
	dnl bypass standard optimization -g -O2 ? 
	AC_ARG_WITH([cxxoptflags],
	  AS_HELP_STRING([--with-cxxoptflags = CXXOPTFLAGS], [optimization using CXX flags, ex: --with-cxxoptflags=-march=opteron -O3]),
	  [CXXOPTFLAGS=$withval],[CXXOPTFLAGS="-g -O2 -fPIC"]) 
	AC_MSG_CHECKING(for c++ optimization flags)
	AC_SUBST([CXXOPTFLAGS])
	AC_MSG_RESULT(done)

	dnl }}}
	dnl multithreading{{{
	AC_ARG_WITH([numthreads],
	  AS_HELP_STRING([--with-numthreads = NUMTHREADS_VALUE],[numthreads, default is 1. ]),
	  [NUMTHREADS_VALUE=$withval],[NUMTHREADS_VALUE=1])
	AC_MSG_CHECKING(for number of threads)
	dnl defaults
	MULTITHREADING=no
	MULTITHREADINLIB=""
	if test "$NUMTHREADS_VALUE" != "1"; then
		
		MULTITHREADINGLIB="-lpthread -lrt"
		case "${host_os}" in
		*cygwin*)
		MULTITHREADINGLIB="-lpthread -lrt"
		;;
		*linux*)
		MULTITHREADINGLIB="-lpthread -lrt"
		;;
		*darwin*)
		MULTITHREADINGLIB="-lpthread"
		;;
		esac
		AC_DEFINE([_MULTITHREADING_],[1],[with numthreads enabled])
	fi
	AC_DEFINE_UNQUOTED([_NUMTHREADS_],[$NUMTHREADS_VALUE],[number of threads])
	AC_SUBST([MULTITHREADINGLIB])
	AC_MSG_RESULT($NUMTHREADS_VALUE) 
	dnl }}}
	dnl 64bit {{{
	AC_ARG_WITH([64bit-indices],
	  AS_HELP_STRING([--with-64bit-indices = bool], [use 64 bit integers, default 0, ex: --with-64bit-indices=1]),
	  [USE_64BIT_INDICES=$withval],[USE_64BIT_INDICES=0]) 
	AC_MSG_CHECKING(for 64 bit indices)

	if test "$USE_64BIT_INDICES" == "1"; then
	AC_DEFINE([ISSM_USE_64BIT_INDICES],[1],[with 64 bits indices])
	else
	AC_DEFINE([ISSM_USE_64BIT_INDICES],[0],[with 64 bits indices])
	fi
	AC_MSG_RESULT($USE_64BIT_INDICES)
	dnl }}}

	dnl Checks
	dnl checks{{{
		AC_MSG_CHECKING(consistency between all libraries)

		  dnl check that if petsc is requested , mpi should be specified
		  if test "$HAVE_PETSC" = "yes" ; then
				if test "$HAVE_MPI" = "NO";  then
					 AC_MSG_ERROR([petsc requires mpi!]);
				fi
		  fi

		  dnl check that we have either python or matlab support if we compile the modules
		  if test "$MODULES_VALUE" = "yes"  && test "$HAVE_MATLAB" = "no" && test "$HAVE_PYTHON" = "no"; then
				AC_MSG_ERROR([need at least python or matlab support to compile modules (or use --with-modules=no)]);
		  fi

		  dnl check that if we have MPI, we have metis
		  if test "$HAVE_METIS" = "yes"  && test "$HAVE_MPI" = "no" ; then
				AC_MSG_ERROR([need mpi if using the metis partitioner!]);
		  fi
		
		dnl check that if we run adolc, we don't compile krigging.exe
		if test "$HAVE_ADOLC" = "yes"  && test "$HAVE_KRIGING" = "yes" ; then
			AC_MSG_ERROR([cannot compile kriging.exe under adolc conditions!]);
		fi
		dnl check that if we run adolc, we don't use PETSc for now
		if test "$HAVE_ADOLC" = "yes"  && test "$HAVE_PETSC" = "yes" ; then
			AC_MSG_ERROR([cannot compile ISSM with both PETSc and adolc]);
		fi
		dnl check that if we run meteoio, we have snowpack also
		if test "$HAVE_METEOIO" = "yes"  && test "$HAVE_SNOWPACK" = "no" ; then
			AC_MSG_ERROR([cannot compile MeteoIO package without Snowpack!]);
		fi
		dnl check that if we run snowpack, we have meteoio also
		if test "$HAVE_METEOIO" = "no"  && test "$HAVE_SNOWPACK" = "yes" ; then
			AC_MSG_ERROR([cannot compile Snowpack package without MeteoIO!]);
		fi

		AC_MSG_RESULT(done)
		dnl }}}
])
