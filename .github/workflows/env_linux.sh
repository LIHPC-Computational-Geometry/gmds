#!/bin/bash

apt-get update
apt-get install -y gcc g++ cmake lcov autoconf libtool git python3 curl
git clone --depth=1 https://github.com/spack/spack.git
sed -i 's#"${ARCHITECTURE}/${COMPILERNAME}-${COMPILERVER}/${PACKAGE}-${VERSION}-${HASH}"#"${PACKAGE}"#g' spack/etc/spack/defaults/config.yaml
. ./spack/share/spack/setup-env.sh
spack external find cmake
spack install kokkos+openmp
