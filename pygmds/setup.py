import os
import re
import subprocess
import platform
import sys

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


def log(msg):
    sys.stderr.write(f"pygmds setup.py: {str(msg)}\n")

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(f"CMake must be installed to build the following extensions: {','.join(e.name for e in self.extensions)}")

        self.cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
        if self.cmake_version < '3.12.1':
            raise RuntimeError("CMake >= 3.12.1 is required")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext: CMakeExtension) -> None:
        env = os.environ.copy()

        ext_dir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        build_dir = os.path.abspath(self.build_temp)
        install_dir = ext_dir

        if not ext_dir.endswith(os.path.sep):
            ext_dir += os.path.sep

        rpaths = [
            '$ORIGIN/../lib', '$ORIGIN/../lib64',
            '$ORIGIN/lib', '$ORIGIN/lib64'
        ]

        cmake_args = [
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_INSTALL_PREFIX={install_dir}",
            f"-DCMAKE_INSTALL_RPATH={':'.join(rpaths)}",
            f"-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON",
        ]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = [ '--config', cfg ]

        if platform.system() == "Windows":
            if sys.maxsize > 2**32:
                cmake_args += [ '-A', 'x64' ]
            build_args += [ '--', '/m' ]
        else:
            cmake_args += [ f"-DCMAKE_BUILD_TYPE={cfg}" ]
            import multiprocessing
            build_args += [ '--', '-j', str(multiprocessing.cpu_count()) ]

        env['CXXFLAGS'] = f'{env.get("CXXFLAGS", "")} -DVERSION_INFO="{self.distribution.get_version()}"'

        if not os.path.exists(build_dir):
            os.makedirs(build_dir)

        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=build_dir)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=build_dir)

        if self.cmake_version >= '3.15.0':
            subprocess.check_call(['cmake', '--install', '.'], cwd=build_dir)
        else:
            if platform.system() == 'Windows':
                raise RuntimeError('"make install" Windows equivalent not implemented! please use cmake >= 3.15')
            else:
                subprocess.check_call(['make', 'install'], cwd=build_dir)

try:
    pygmds_version = open('version.txt', 'r').read().strip()
except:
    raise RuntimeError("Could not find pygmds's 'version.txt' file. Please make sure you are in pygmds's root directory.")

from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="pygmds",
    version=pygmds_version,
    author="Franck Ledoux",
    author_email="franck.ledoux@cea.fr",
    description="Python bindings for GMDS meshing library",
    long_description=long_description,
    ext_modules=[CMakeExtension("pygmds", sourcedir="src")],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
    extras_require={"test": ["pytest>=6.0"]},
    python_requires=">=3.7",
)
