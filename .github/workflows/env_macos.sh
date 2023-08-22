#!/bin/bash


git clone --depth=1 -b v0.19.2 https://github.com/spack/spack.git

#replace the default configuration file of spack by a simpler one without hash, compiler versions, tags and so on 
cp /Users/runner/work/gmds/gmds/.github/workflows/misc/config.yaml /Users/runner/work/gmds/gmds/spack/etc/spack/defaults/
. ./spack/share/spack/setup-env.sh

git clone --branch gmds_temp --depth=1 https://github.com/LIHPC-Computational-Geometry/spack_recipes.git

spack repo add ./spack_recipes/meshing_repo
spack repo add ./spack_recipes/supersede_repo

spack external find cmake
spack install --only dependencies gmds+blocking~cgns
