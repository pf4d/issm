# Process this file with autoconf to produce a configure script.

#AUTOCONF
AC_INIT([Triangle],[1.0]) #Initializing configure
AC_CONFIG_AUX_DIR([./aux-config])         #Put config files in aux-config
AC_CONFIG_MACRO_DIR([m4])                 #m4 macros are located in m4
AC_CANONICAL_TARGET                       #Determine the system type
AC_CANONICAL_HOST                      #Determine the system type

case ${host_os} in
	*linux* )
		CPP_LINUX='-DLINUX'
		;;
	*)
		CPP_LINUX=''
esac

AC_SUBST(CPP_LINUX)

AC_ARG_ENABLE([executables],                                                 dnl feature
   AS_HELP_STRING([--enable-executables], [build executables]),              dnl help string
   [enable_executables=$enableval],                                          dnl action if given
   [enable_executables=yes])                                                 dnl action if not given
AC_MSG_CHECKING(for executables)
AM_CONDITIONAL([EXECUTABLES], [test x$enable_executables = xyes])
AC_MSG_RESULT($enable_executables)

#Compilers
AC_PROG_CC([icc cccl cl icl gcc])
AC_PROG_CPP
AC_PROG_CXX([icpc cccl cl icl g++])

#Libraries and linking
AC_PATH_XTRA      #figure out X library and include paths

#AUTOMAKE
#Initialize automake and declare foreign so that we don't need a ChangeLog, INSTALL, etc
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_HEADERS([./config.h])   #Config file must be config.h
AM_SILENT_RULES([yes])           #Do not show compilation command by default
AM_PROG_AR

#Libtool
LT_INIT

#List all Makefiles
AC_CONFIG_FILES([
					  Makefile
					  src/Makefile
					  ])

#End of configure.ac
AC_OUTPUT
