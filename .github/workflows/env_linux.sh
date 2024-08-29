#!/bin/bash

git clone --depth=1 -b v0.20.3 https://github.com/spack/spack.git
sed -i 's#"{architecture}/{compiler.name}-{compiler.version}/{name}-{version}-{hash}"#"{name}"#g' spack/etc/spack/defaults/config.yaml
source ./spack/share/spack/setup-env.sh

git clone --depth=1 https://github.com/LIHPC-Computational-Geometry/spack_recipes.git

spack repo add ./spack_recipes/meshing
spack repo add ./spack_recipes/meshing_supersede

spack external find cmake
spack install lcov
spack install py-pybind11
spack install glpk
spack install googletest
spack install --only dependencies gmds+kmds+blocking+cgns ^kokkos+openmp ^cgns~mpi
