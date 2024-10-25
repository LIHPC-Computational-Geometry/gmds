#!/bin/bash

git clone --depth=1 -b v0.22.2 https://github.com/spack/spack.git

#replace the default configuration file of spack by a simpler one without hash, compiler versions, tags and so on 
#cp /Users/runner/work/gmds/gmds/.github/workflows/misc/config-0.21.1.yaml /Users/runner/work/gmds/gmds/spack/etc/spack/defaults/config.yaml
source ./spack/share/spack/setup-env.sh

git clone --depth=1 https://github.com/LIHPC-Computational-Geometry/spack_recipes.git

spack repo add ./spack_recipes/meshing
spack repo add ./spack_recipes/meshing_supersede

spack external find cmake
spack install --only dependencies gmds+blocking~cgns+python
