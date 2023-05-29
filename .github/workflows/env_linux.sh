#!/bin/bash

git clone --depth=1 -b releases/latest https://github.com/spack/spack.git
cp /Users/runner/work/gmds/gmds/.github/workflows/misc/config.yaml /Users/runner/work/gmds/gmds/spack/etc/spack/defaults/
#sed -i 's#"${ARCHITECTURE}/${COMPILERNAME}-${COMPILERVER}/${PACKAGE}-${VERSION}-${HASH}"#"${PACKAGE}"#g' spack/etc/spack/defaults/config.yaml
. ./spack/share/spack/setup-env.sh

git clone --branch gmds_temp --depth=1 https://github.com/LIHPC-Computational-Geometry/spack_recipes_meshing.git

spack repo add ./spack_recipes_meshing/meshing_repo
spack repo add ./spack_recipes_meshing/supersede_repo

spack external find cmake
spack install lcov
spack install py-pybind11
spack install glpk
spack install googletest
spack install cgal
spack install --only dependencies gmds+kmds+blocking ^kokkos+openmp ^cgns~mpi
