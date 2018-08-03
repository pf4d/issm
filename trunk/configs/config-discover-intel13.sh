#!/bin/csh

#PETSc 3.2
#MPI /usr/local/other/SLES11.1/mvapich2/1.8.1/intel-13.1.2.183/lib

./configure \
 --prefix=$ISSM_DIR \
 --with-matlab-dir="$ISSM_DIR/externalpackages/matlab/install" \
 --with-triangle-dir=$ISSM_DIR/externalpackages/triangle/install \
 --with-metis-dir=$ISSM_DIR/externalpackages/metis/install \
 --with-petsc-dir=$ISSM_DIR/externalpackages/petsc/install \
 --with-mkl-dir="/usr/local/intel/Composer/composer_xe_2013.3.163/mkl/" \
 --with-mpi-include="/usr/local/other/SLES11.1/mvapich2/1.8.1/intel-13.1.2.183/include" \
 --with-mpi-libflags="-L/usr/local/other/SLES11.1/mvapich2/1.8.1/intel-13.1.2.183/lib -lmpich -lopa -lmpl -lfmpich -lmpichcxx -lmpichf90 -lpthread -lrdmacm -libverbs -libumad -lrt -lnuma" \
 --with-petsc-arch=$ISSM_ARCH \
 --with-blacs-dir="/usr/local/intel/Composer/composer_xe_2013.3.163/mkl/" \
 --with-blas-lapack-dir="/usr/local/intel/Composer/composer_xe_2013.3.163/mkl/" \
 --with-mumps-dir=$ISSM_DIR/externalpackages/petsc/install \
 --with-cxxoptflags="-O3 -xS -DMPICH_IGNORE_CXX_SEEK" \
 --with-vendor=intel-discover

 # --with-esmf-dir=$ESMF_INSTALL_DIRECTORY \

