#!/bin/bash

version=2018.7
tgz_file=gromacs-$version.tar.gz
wget ftp://ftp.gromacs.org/pub/gromacs/$tgz_file
tar xfz $tgz_file

dirname=~/opt/gromacs/$version/
suffix=_$version

echo $dirname

cd gromacs-$version
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX=$dirname/single -DGMX_DEFAULT_SUFFIX=OFF -DGMX_DOUBLE=OFF -DGMX_BINARY_SUFFIX=$suffix -DGMX_LIBS_SUFFIX=$suffix
make -j4
make install
make clean

cmake .. -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX=$dirname/double -DGMX_DEFAULT_SUFFIX=OFF -DGMX_DOUBLE=ON -DGMX_BINARY_SUFFIX=$suffix_d -DGMX_LIBS_SUFFIX=$suffix_d
make -j4
make install
make clean

cd ..
