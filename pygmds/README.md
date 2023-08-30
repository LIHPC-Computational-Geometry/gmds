# Python binding for gmds

Python binding is a 2023 feature of gmds to help writing high-level meshing algorithms. We mainly focus
on machine learning and blocking subject right now.

## Structure
The python binding provide a python module **gmds** and several submodules relative to gmds components:
- [blocking](src/binding_blocking.cpp) - provides blocking operations.
- [geometry](src/binding_geometry.cpp) - provides access to the geometrical kernel to queries faceted geometric model to write meshing algorithms.
- [math](src/binding_math.cpp) - gives access to some basic mathematical concepts handled by gmds.
- [mesh](src/binding_mesh.cpp) - provides meshing operations including vtk IO.


## How to use the gmds Python API?

The python module packaging is not up-to-date right now. We advise you to apply the following
procedure on Linux platform.

```bash
#==========================================
#First we clone spack and change his way for installing
git clone --depth=1 -b v0.19.2  https://github.com/spack/spack.git
sed -i 's#"${ARCHITECTURE}/${COMPILERNAME}-${COMPILERVER}/${PACKAGE}-${VERSION}-${HASH}"#"${PACKAGE}"#g' spack/etc/spack/defaults/config.yaml
. ./spack/share/spack/setup-env.sh
#==========================================
#mandatory if you have already used back on your computer
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
export PYTHONPATH=spack/opt/spack/gmds/lib:$PYTHONPATH
spack load gmds
```

Our procedure relies on [spack](https://spack.io/), a package manager for supercomputers, Linux and macOS. With the previous
set of commands, you are going to build a complete gmds environment that you can use in your Python code.
