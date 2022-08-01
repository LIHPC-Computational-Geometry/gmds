#!/bin/bash

git clone --depth=1 -b releases/latest https://github.com/spack/spack.git
sed -i 's#"${ARCHITECTURE}/${COMPILERNAME}-${COMPILERVER}/${PACKAGE}-${VERSION}-${HASH}"#"${PACKAGE}"#g' spack/etc/spack/defaults/config.yaml
. ./spack/share/spack/setup-env.sh

git clone --depth=1 https://github.com/LIHPC-Computational-Geometry/spack_recipes_meshing.git

spack repo add ./spack_recipes_meshing/meshing_repo
spack repo add ./spack_recipes_meshing/supersede_repo

spack external find cmake
spack install lcov
spack install --only dependencies gmds+kmds+blocking ^kokkos+openmp
