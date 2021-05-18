#!/bin/bash

git clone --depth=1 https://github.com/spack/spack.git
echo"${ARCHITECTURE}/${COMPILERNAME}-${COMPILERVER}/${PACKAGE}-${VERSION}-${HASH}"
echo ${ARCHITECTURE}/${COMPILERNAME}-${COMPILERVER}/${PACKAGE}-${VERSION}-${HASH}
sed 's/"${ARCHITECTURE}/${COMPILERNAME}-${COMPILERVER}/${PACKAGE}-${VERSION}-${HASH}"/"${PACKAGE}"/g' spack/etc/spack/defaults/config.yaml

. ./spack/share/spack/setup-env.sh
spack external find cmake
spack install lcov
