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

# Defining compiler optimisations
FOPTFLAGS="-O2 -march=native -mtune=native"
COPTFLAGS="-O2 -march=native -mtune=native"
CXXOPTFLAGS="-O2 -march=native -mtune=native"

# Defining directories
E4D_DIR=$(pwd)
BX_DIR="${E4D_DIR}/utilities/bx"
VENDOR_DIR="${E4D_DIR}/vendor"
NETCDF_DIR="${VENDOR_DIR}/netcdf-c-4.9.0"
EXO_DIR="${VENDOR_DIR}/exodus-5.09/exodus"


# Create required directories if they do not exist
mkdir -p "${E4D_DIR}/bin"
mkdir -p "${VENDOR_DIR}"

# # Remove directories in vendor to get rid of previous configurations
cd "$VENDOR_DIR" || exit 1
sudo rm -rf "${EXO_DIR}" "${NETCDF_DIR}"

# Unzip vendor packages
tar -xzf netcdf-c-4.9.0.tar.gz
tar -xzf exodus-5.09.tar.gz

# Configure netcdf-c
cd "${NETCDF_DIR}/cmake" || exit 1
cmake -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DENABLE_NETCDF_4=OFF -DENABLE_DAP=OFF -DCMAKE_C_FLAGS="-O0" -DCMAKE_CXX_FLAGS="-O0" ..
make 
sudo make install

# Configure exodus
cd "${EXO_DIR}" || exit 1
make -f Makefile.standalone USING_NETCDF4='"NO"' CC="${CC}" FC="${FC}" NETCDF="${NETCDF_DIR}"

# Compile bx. Remove .mod and .o files after compiling BX
cd "${BX_DIR}" || exit 1
make f95="${FC}" FFLAGS="${FOPTFLAGS} -r8 -heap-arrays" EXODUS="${EXO_DIR}" NETCDF="${NETCDF_DIR}"
mv bx "${E4D_DIR}/bin"
find . -name '*.mod' -delete
find . -name '*.o' -delete

# Set permissions
chmod 755 -R "${E4D_DIR}/bin"

echo "Installation completed successfully!"
