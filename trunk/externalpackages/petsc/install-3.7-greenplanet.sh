#!/bin/bash
set -eu

#WARNING: make sure you have the right mpi

#Some cleanup
rm -rf install petsc-3.7.6 src
mkdir install src

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/petsc-lite-3.7.6.tar.gz' 'petsc-3.7.6.tar.gz'

#Untar and move petsc to install directory
tar -zxvf  petsc-3.7.6.tar.gz
mv petsc-3.7.6/* src/
rm -rf petsc-3.7.6

#configure
cd src
./config/configure.py \
 --prefix="$ISSM_DIR/externalpackages/petsc/install" \
 --PETSC_DIR="$ISSM_DIR/externalpackages/petsc/src" \
 --with-blas-lapack-dir="/sopt/Intel/compilers_and_libraries_2016.2.181/linux/mkl/" \
 --with-mpi-dir="/sopt/mpi/openmpi-1.10.2/intel_16.0.2/" \
 --known-mpi-shared-libraries=1 \
 --with-debugging=0 \
 --with-valgrind=0 \
 --with-x=0 \
 --with-ssl=0 \
 --with-batch=1  \
 --with-shared-libraries=1 \
 --download-metis=1 \
 --download-parmetis=1 \
 --download-scalapack=1 \
 --download-mumps=1 

#prepare script to reconfigure petsc
cat > script.queue << EOF
#!/bin/bash
#SBATCH -p c6145
#SBATCH -N 1 -n 1
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=10
#SBATCH --job-name=test

module load compiler/intel/16.0.2
module load mpi/openmpi/1.10.2/intel_16.0.2

cd $(echo $ISSM_DIR)/externalpackages/petsc/src/
mpiexec -np 1 ./conftest-arch-linux2-c-opt
EOF

#print instructions
echo "== Now: cd src/ "
echo "== sbatch script.queue "
echo "== Then run reconfigure script generated by PETSc and follow instructions"
