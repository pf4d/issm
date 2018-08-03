#!/bin/bash
set -eu

#Some cleanup
rm -rf install petsc-3.6.3 src
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/petsc-lite-3.6.3.tar.gz' 'petsc-3.6.3.tar.gz'

#Untar and move petsc to install directory
tar -zxvf  petsc-3.6.3.tar.gz
mv petsc-3.6.3/* src/
rm -rf petsc-3.6.3

#--with-cc=icc --with-cxx=icpc --with-fc=ifort --with-f77=ifort \

#configure
cd src
./config/configure.py \
	--prefix="$ISSM_DIR/externalpackages/petsc/install" \
	--PETSC_DIR="$ISSM_DIR/externalpackages/petsc/src" \
	--with-blas-lapack-dir="/usr/local/intel/Composer/composer_xe_2015.0.090/mkl/" \
	--with-mpi-lib="/usr/local/intel/mpi/4.0.3.008/lib64/libmpi.so" \
	--with-mpi-include="/usr/local/intel/mpi/4.0.3.008/intel64/include/" \
	--known-mpi-shared-libraries=1 \
	--with-debugging=0 \
	--with-valgrind=0 \
	--with-x=0 \
	--with-ssl=0 \
	--with-batch=1  \
	--with-shared-libraries=1 \
	--download-metis=1 \
	--download-parmetis=1 \
	--download-mumps=1 \
	--download-scalapack=1 

#prepare script to reconfigure petsc
cat > script.queue << EOF
#!/bin/bash
#SBATCH -J petscinstall # Job Name
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A s1010
#SBATCH -t 00:01:00 # Run time (hh:mm:ss) - 1.5 hours
#SBATCH --qos=debug
#SBATCH -o petscinstall.outlog
#SBATCH -e petscinstall.errlog

. /usr/share/modules/init/bash
module load comp/intel-15.0.0.090
module load mpi/impi-4.0.3.008

export PATH="$PATH:."
export MPI_GROUP_MAX=64
mpirun -np 1 ./conftest-arch-linux2-c-opt
EOF

#print instructions
echo "== Now: cd src/ "
echo "== sbatch script.queue "
echo "== Then run reconfigure script generated by PETSc and follow instructions"
