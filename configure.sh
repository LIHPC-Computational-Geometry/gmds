#!/bin/bash

#==========================================#
mkdir build
cd build
#==========================================#
#First we clone spack and change his way for installing
git clone --depth=1  https://github.com/spack/spack.git
cp ../config-spack.yaml spack/etc/spack/defaults/config.yaml
. ./spack/share/spack/setup-env.sh
#==========================================
#mandatory if you have already used spack on your computer
spack clean -m
#==========================================
git clone --branch gmds_temp --depth=1 https://github.com/LIHPC-Computational-Geometry/spack_recipes_meshing.git
spack repo add ./spack_recipes_meshing/meshing_repo
spack repo add ./spack_recipes_meshing/supersede_repo
spack compiler find
spack external find cmake
#==========================================
spack install --no-checksum gmds@main+python+blocking~cgns
#==========================================
# testing the install
export PYTHONPATH=spack/opt/spack/gmds/lib:$PYTHONPATH
spack load py-pytest
pytest gmds/pygmds/tst gmds/test_samples
#==========================================