#!/bin/bash
set -e
set -o pipefail

# Check if MKLROOT is set
if [ -z "$MKLROOT" ]; then
  echo "MKLROOT is not set. Please set MKLROOT environment variable."
  echo "Try source <intel-compiler-dir>/setvars.sh"
  exit 1
fi

# Defining compilers
CC=icx
CXX=icpx
FC=ifx
MPIFC="mpiifort -fc=ifx"
MPICC="mpiicc -cc=icx"
MPICXX="mpiicpc -cxx=icpx"

# Defining compiler optimisations
FOPTFLAGS="-O2 -march=native -mtune=native"
COPTFLAGS="-O2 -march=native -mtune=native"
CXXOPTFLAGS="-O2 -march=native -mtune=native"

# Defining directories
E4D_DIR=$(pwd)
VENDOR_DIR="${E4D_DIR}/vendor"
PETSC_DIR="${VENDOR_DIR}/petsc-3.19.5"
TRI_DIR="${VENDOR_DIR}/triangle-1.6"
TET_DIR="${VENDOR_DIR}/tetgen-1.6.0"

# Create required directories if they do not exist
mkdir -p "${E4D_DIR}/bin"
mkdir -p "$VENDOR_DIR"

# Remove directories in vendor to get rid of previous configurations
rm -rf "${PETSC_DIR}" "${TRI_DIR}" "${TET_DIR}"

# Unzip vendor packages
cd "$VENDOR_DIR" || exit 1
tar -xzf petsc-3.19.5.tar.gz
tar -xzf triangle-1.6.tar.gz
tar -xzf tetgen-1.6.0.tar.gz

# Configure petsc
cd "${PETSC_DIR}" || exit 1
./configure  --with-debugging=0 --with-fc="${MPIFC}" --with-cc="${MPICC}" --with-cxx="${MPICXX}" --with-blaslapack-dir="${MKLROOT}" FOPTFLAGS="${FOPTFLAGS}" COPTFLAGS="${COPTFLAGS}" CXXOPTFLAGS="${CXXOPTFLAGS}"
make PETSC_DIR="${PETSC_DIR}" PETSC_ARCH=arch-linux-c-opt all
make PETSC_DIR="${PETSC_DIR}" PETSC_ARCH=arch-linux-c-opt check

# Compile E4D. Remove .mod and .o files after compiling E4D
cd "${E4D_DIR}/src" || exit 1
make FC="${MPIFC}" FFLAGS="${FOPTFLAGS} -r8 -heap-arrays" PETSC_DIR="${PETSC_DIR}"
mv e4d "${E4D_DIR}/bin"
find . -name '*.mod' -delete
find . -name '*.o' -delete

# Compile triangle
cd "${TRI_DIR}" || exit 1
make CC="${CC}"
mv triangle "${E4D_DIR}/bin"
find . -name '*.o' -delete

# Compile tetgen
cd "${TET_DIR}" || exit 1
make CXX="${CXX}"
mv tetgen "${E4D_DIR}/bin"
find . -name '*.o' -delete

# Set permissions
chmod 755 -R "${E4D_DIR}/bin"

echo "Installation completed successfully!"
