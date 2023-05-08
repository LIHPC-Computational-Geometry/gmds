#!/bin/bash

brew install autoconf automake libtool gcc

git clone --depth=1 -b releases/latest https://github.com/spack/spack.git

#replace the default configuration file of spack by a simpler one without hash, compiler versions, tags and so on 
cp /Users/runner/work/gmds/gmds/.github/workflows/misc/config.yaml /Users/runner/work/gmds/gmds/spack/etc/spack/defaults/
. ./spack/share/spack/setup-env.sh

git clone --branch gmds_temp --depth=1 https://github.com/LIHPC-Computational-Geometry/spack_recipes_meshing.git

spack repo add ./spack_recipes_meshing/meshing_repo
spack repo add ./spack_recipes_meshing/supersede_repo

spack external find cmake
spack install lcov
spack install pybind11
spack install glpk
spack install googletest
spack install --only dependencies gmds+blocking ^cgns~mpi
