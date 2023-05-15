#!/bin/bash

CC=icc
CXX=icpc
FC=ifort
MPIFC=mpiifort
MPICC=mpiicc
MPICXX=mpiicpc

dir=$(pwd)

cd ${dir}/
mkdir bin
mkdir third_party

# Download third party dependencies
cd ${dir}/third_party/
deps=(
    petsc 
    netcdf 
    exodus 
    triangle
    tetgen)

repos=(
    https://gitlab.com/petsc/petsc.git
    "--depth 1 --branch v4.9.0 https://github.com/Unidata/netcdf-c"
    https://github.com/certik/exodus.git
    https://github.com/jdumas/triangle
    "--depth 1 --branch upstream/1.5.0 https://salsa.debian.org/science-team/tetgen.git"
)

for i in "${!repos[@]}"; do
    if [ ! -d "${dir}/third_party/${deps[i]}" ] 
    then
        echo "Downloading latest stable version of ${deps[i]}"
        git clone --depth=1 ${repos[i]} ${deps[i]}
    else
        echo "Latest stable version of ${deps[i]} already downloaded"
    fi 
done
    
#Configure petsc
cd ${dir}/third_party/petsc
./configure --with-cc=mpiicc --with-fc=mpiifort --with-cxx=mpiicpc --with-blaslapack-dir=$MKLROOT --with-debugging=0 >${dir}/setup.log
make PETSC_DIR=${dir}/third_party/petsc PETSC_ARCH=arch-linux-c-opt all
make PETSC_DIR=${dir}/third_party/petsc PETSC_ARCH=arch-linux-c-opt check
cp ${dir}/third_party/petsc/arch-linux-c-opt/bin/mpirun ${dir}/bin

#Compile E4D
cd ${dir}/src
petscdir=$dir/third_party/petsc/
make PETSC_DIR=$petscdir
cp ${dir}/src/e4d ${dir}/bin 

# Compile netcdf
cd ${dir}/third_party/netcdf
./configure --disable-netcdf-4 --disable-dap >${dir}/setup.log
make
make install

# Compile exodus
cd ${dir}/third_party/exodus/exodus
netcdfdir=${dir}/third_party/netcdf
make -f Makefile.standalone USING_NETCDF4='"NO"' NETCDF=$netcdfdir

# Compile bx
cd ${dir}/utilities/bx
exodusdir=${dir}/third_party/exodus/exodus
make f95=mpif90 NETCDF=$netcdfdir EXODUS=$exodusdir
cp bx ${dir}/bin

# Compile triangle
cd ${dir}/third_party/triangle
make CC=mpiicc
cp triangle ${dir}/bin

## Compile tetgen
cd ${dir}/third_party/tetgen
make cxx=mpiicpc
cp tetgen ${dir}/bin

chmod 755 -R ${dir}/bin
