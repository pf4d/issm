#!/bin/bash

# ISSM_DIR and MATLAB_DIR must be set correctly.
# {{{
if [ "${ISSM_DIR}x" == "x" ]; then
   echo "ISSM_DIR is not set!" >&2
   exit 1;
elif [ -d ${ISSM_DIR} ]; then
   echo "ISSM_DIR: ${ISSM_DIR}"
else
   echo "ISSM_DIR: ${ISSM_DIR} does not exist!" >&2
   exit 1;
fi

if [ "${MATLAB_DIR}x" == "x" ]; then
   echo "MATLAB_DIR is not set!"
   exit 1;
elif [ -d ${MATLAB_DIR} ]; then
   echo "MATLAB_DIR: ${MATLAB_DIR}"
else
   echo "MATLAB_DIR: ${MATLAB_DIR} does not exist!" >&2
   exit 1;
fi
# }}}

#List of external pakages to be installed and their installation scripts
EXTERNALPACKAGES="autotools install.sh                
						mpich     install-3.0-linux64.sh
                  cmake     install.sh                
                  petsc     install-3.5-linux64-static.sh
						m1qn3     install.sh    
                  triangle  install-linux64.sh "

# Install Externalpackages
# {{{

#Files source environment to make sure installed packages are in PATH
source $ISSM_DIR/etc/environment.sh

#number of packages: 
NUMPACKAGES=$(($(echo $EXTERNALPACKAGES | wc -w )/2))

for ((i=1;i<=$NUMPACKAGES;i++))
do
	NUM1=$((2*$i-1))
	NUM2=$((2*$i))
	PACKAGENAME=$(echo $EXTERNALPACKAGES | cut -d " " -f $NUM1-$NUM1)
	PACKAGEINST=$(echo $EXTERNALPACKAGES | cut -d " " -f $NUM2-$NUM2)

	cd $ISSM_DIR/externalpackages/$PACKAGENAME

	if [[ $PACKAGEINST -nt  install ]]; then 
		#go ahead and reinstall. 
		echo "Triggering new install of $PACKAGENAME"
		install_test=1
	else
		#ok, we want to skip, unless the package is not installed: 
		if [ -d install ]; then 
			#could be empty, signaling a failed previous install: 
			if [ "$(ls -A install)" ];then
				echo "and install directory exists, so skipping install of $PACKAGENAME"
				install_test=0;
			else
				echo "and install directory exists, however, it is empty, so triggering install of $PACKAGENAME"
				install_test=1;
			fi
		else
			echo "However, install directory does not exist, so triggering install of $PACKAGENAME"
			install_test=1;
		fi
	fi

	if [[ $install_test == 1 ]]; then 
		echo "======================================================";
		echo "       Installing $PACKAGENAME                        ";
		echo "======================================================";
		./$PACKAGEINST |  tee compil.log
		if [ $? -ne 0 ]; then
			echo "======================================================";
			echo "    ERROR: installation of $PACKAGENAME failed        ";
			echo "======================================================";
			#erase install directory, so that next time, we still try and compile this!
			rm -rf install
		fi
		source $ISSM_DIR/etc/environment.sh
	else
		echo "======================================================";
		echo "       Skipping install of $PACKAGENAME                        ";
		echo "======================================================";
	fi
	cd ..
done
# }}}

# Compile ISSM
#{{{
cd $ISSM_DIR
echo "Aureconf..."
autoreconf -if
echo "Configuring..."
./configure \
	--prefix=$ISSM_DIR \
	--disable-static \
	--enable-standalone-executables \
	--enable-standalone-libraries \
	--with-matlab-dir="/usr/local/matlab80/" \
	--with-triangle-dir=$ISSM_DIR/externalpackages/triangle/install \
	--with-metis-dir=$ISSM_DIR/externalpackages/petsc/install \
	--with-petsc-dir=$ISSM_DIR/externalpackages/petsc/install  \
	--with-scalapack-dir=$ISSM_DIR/externalpackages/petsc/install \
	--with-mumps-dir=$ISSM_DIR/externalpackages/petsc/install \
	--with-blas-lapack-dir=$ISSM_DIR/externalpackages/petsc/install \
	--with-mpi-include=$ISSM_DIR/externalpackages/mpich/install/include  \
	--with-mpi-libflags=" -L$ISSM_DIR/externalpackages/mpich/install/lib -lmpich -lmpl " \
	--with-m1qn3-dir=$ISSM_DIR/externalpackages/m1qn3/install \
	--with-fortran-lib="-L/usr/lib/gcc/x86_64-linux-gnu/4.6/ -lgfortran" \
	--enable-debugging \
	--with-numthreads=4


if [ $? -ne 0 ]; then echo "FAILED TO CONFIGURE!" && exit 1; fi

echo "Building..."
make -j 8
if [ $? -ne 0 ]; then echo "FAILED TO BUILD!" && exit 1; fi

echo "Installing..."
make install 
if [ $? -ne 0 ]; then echo "FAILED TO INSTALL!" && exit 1; fi
#}}}
