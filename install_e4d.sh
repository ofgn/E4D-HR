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
FOPTFLAGS="-O2 -xHost"
COPTFLAGS="-O2 -xHost"
CXXOPTFLAGS="-O2 -xHost"

# Defining directories
E4D_DIR=$(pwd)
VENDOR_DIR="${E4D_DIR}/vendor"
PETSC_DIR="${VENDOR_DIR}/petsc-3.21.4"
TRI_DIR="${VENDOR_DIR}/triangle-1.6"
TET_DIR="${VENDOR_DIR}/tetgen-1.6.0"

# Create required directories if they do not exist
mkdir -p "${E4D_DIR}/bin"
mkdir -p "$VENDOR_DIR"


# Compile E4D. Remove .mod and .o files after compiling E4D
cd "${E4D_DIR}/src" || exit 1
make FC="${MPIFC}" FFLAGS="${FOPTFLAGS}" PETSC_DIR="${PETSC_DIR}"
mv e4d "${E4D_DIR}/bin"
find . -name '*.mod' -delete
find . -name '*.o' -delete

# Set permissions
chmod 755 -R "${E4D_DIR}/bin"

echo "Installation completed successfully!"