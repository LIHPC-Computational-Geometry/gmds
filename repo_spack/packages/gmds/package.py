# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

# ----------------------------------------------------------------------------
# If you submit this package back to Spack as a pull request,
# please first remove this boilerplate and all FIXME comments.
#
# This is a template package file for Spack.  We've put "FIXME"
# next to all the things you'll want to change. Once you've handled
# them, you can save this file and test your package like this:
#
#     spack install gmds
#
# You can edit this file again by typing:
#
#     spack edit gmds
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack import *


class Gmds(CMakePackage):
    """GMDS: Generic Mesh Data and Services."""

    homepage = "https://www.example.com"
    url      = "https://www.example.com/gmds-1.0.0.tar"

    version('1.0.0')

    variant('kmds', default=False, description='Build with Kokkos')
    variant('elg3d', default=False, description='Build Elg3D')

    depends_on('glpk')
    # necessary to build the internal glpk
    depends_on('libtool', type='build')
    depends_on('eigen')

    depends_on('kokkos', when='+kmds')
    depends_on('gts', when='+elg3d')
    # necessary to find gts
    depends_on('pkg-config', type='build', when='+elg3d')
    depends_on('exodusii', when='+elg3d')

    conflicts('+elg3d', when='~kmds',
             msg='elg3d cannot be built without kmds.')

    def cmake_args(self):
        args = []
        args.append(self.define_from_variant('ENABLE_KMDS', 'kmds'))
        args.append(self.define_from_variant('ENABLE_ELG3D', 'elg3d'))
        return args
