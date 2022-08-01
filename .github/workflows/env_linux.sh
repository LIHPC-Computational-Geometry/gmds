#!/bin/bash

# will return to the latest of spack when a more recent version of CGAL is included
# (my commit is in branch develop at the moment)
#git clone --depth=1 -b releases/latest https://github.com/spack/spack.git
git clone --branch develop https://github.com/spack/spack.git
sed -i 's#"${ARCHITECTURE}/${COMPILERNAME}-${COMPILERVER}/${PACKAGE}-${VERSION}-${HASH}"#"${PACKAGE}"#g' spack/etc/spack/defaults/config.yaml
. ./spack/share/spack/setup-env.sh

spack external find cmake
spack install lcov
spack install kokkos+openmp
spack install glpk
spack install cgal
