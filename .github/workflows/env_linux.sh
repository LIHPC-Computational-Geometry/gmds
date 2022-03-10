#!/bin/bash

git clone --depth=1 -b releases/latest https://github.com/spack/spack.git
sed -i 's#"${ARCHITECTURE}/${COMPILERNAME}-${COMPILERVER}/${PACKAGE}-${VERSION}-${HASH}"#"${PACKAGE}"#g' spack/etc/spack/defaults/config.yaml
. ./spack/share/spack/setup-env.sh
spack external find cmake
spack install lcov
spack install kokkos+openmp
spack install glpk
