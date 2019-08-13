#!/bin/bash
set -e  # Exit immediately if any command fails

#module load gcc/7.3.0
#module load cuda/10.1.168
#module swap craype-{haswell,mic-knl}

echo "> module list"
module list

LOCAL_BUILD_DIR=builddir_hsw
LOCAL_INSTALL_DIR=instdir_hsw
CC=$(which icc)
CXX=$(which icpc)
CXXHOST=$(which icpc)
ORIG_DIR=$(pwd)

SUNDIALS_DIST=$(pwd)/sundials
SUNDIALS_WWW=https://github.com/LLNL/sundials

SUITESPARSE_WWW=http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-5.4.0.tar.gz
SUITESPARSE_DIR=$(pwd)/SuiteSparse
if [ ! -d ${SWUITESPARSE_DIR} ]
then
  echo "Building SuiteSparse in ${SUITESPARSE_DIR}"
  wget ${SUITESPARSE_WWW}
  tar zxf SuiteSparse-5.4.0.tar.gz
  cd ${SUITESPARSE_DIR}
  make
  make install
  cd ${ORIG_DIR}
fi
echo "SuiteSparse installed in ${SUITESPARSE_DIR}/{lib,include}"

BUILD_DIR=${SUNDIALS_DIST}/${LOCAL_BUILD_DIR}
echo "Building sundials in: " ${BUILD_DIR}

INSTALL_PREFIX=${SUNDIALS_DIST}/${LOCAL_INSTALL_DIR}
echo "Installing sundials in: " $INSTALL_PREFIX

if [ ! -d ${SUNDIALS_DIST} ]
then
  git clone ${SUNDIALS_WWW} ${SUNDIALS_DIST}
  cd ${SUNDIALS_DIST}
  git checkout develop
  git apply ../sundials-5.0.0-dev.1+cuda_nvector_allocators.patch
  git apply ../sundials-5.0.0-dev.1+cuda_nvector_allocators_part2.patch
  cd ${ORIG_DIR}
fi

cd ${SUNDIALS_DIST}

if [ ! -d ${INSTALL_PREFIX} ]; then
  mkdir ${INSTALL_PREFIX}
fi

if [ ! -d ${BUILD_DIR} ]; then
  mkdir ${BUILD_DIR}
fi

cd ${BUILD_DIR}

cmake     \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}     \
    -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON     \
    -DCMAKE_C_COMPILER=${CC}     \
    -DCMAKE_CXX_COMPILER=${CXX}     \
    -DCUDA_HOST_COMPILER=${CXXHOST}    \
    -DEXAMPLES_INSTALL_PATH=${INSTALL_PREFIX}/examples     \
    -DCMAKE_BUILD_TYPE=Release     \
    -DCMAKE_C_FLAGS_RELEASE="-O3 -DNDEBUG"     \
    -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG"     \
    -DCUDA_ENABLE=ON     \
    -DMPI_ENABLE=OFF     \
    -DKLU_ENABLE=ON \
    -DKLU_INCLUDE_DIR=${SUITESPARSE_DIR}/include \
    -DKLU_LIBRARY_DIR=${SUITESPARSE_DIR}/lib \
    -DOPENMP_ENABLE=ON \
    -DEXAMPLES_ENABLE_C=OFF ../

#build with -j fails on:  Building NVCC (Device) object examples/nvector/mpicuda/CMakeFiles/test_nvector_mpicuda.dir/test_nvector_mpicuda_generated_test_nvector_mpicuda.cu.o
make
make install

echo ""
echo "SUNDIALS build successful.  To use with AMReX codes, set CVODE_LIB_DIR=${INSTALL_PREFIX}/lib64"
module list 2> ${INSTALL_PREFIX}/install_modules.txt

cd ${ORIG_DIR}

