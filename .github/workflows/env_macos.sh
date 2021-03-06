#!/bin/bash

brew install autoconf automake libtool gcc

git clone --depth=1 -b releases/latest https://github.com/spack/spack.git
. ./spack/share/spack/setup-env.sh


#replace the default configuration file of spack by a simpler one without hash, compiler versions, tags and so on 
cp /Users/runner/work/gmds/gmds/.github/workflows/misc/config.yaml /Users/runner/work/gmds/gmds/spack/etc/spack/defaults/

. ./spack/share/spack/setup-env.sh
spack external find cmake
spack install lcov
spack install glpk
spack install cgal@5.0.3
ls /Users/runner/work/gmds/gmds/spack/opt/spack/glpk
ls -l /Users/runner/work/gmds/gmds/spack/opt/spack/glpk/lib
