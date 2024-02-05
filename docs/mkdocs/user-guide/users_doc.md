GMDS is a set of C++ libraries that can be used to develop meshing algorithms. It requires features of C++14.

The build system is based on CMAKE (https://cmake.org) and unit tests are performed using google tests and google benchmarks.

In order to perform some linear algebra operations, we rely onto the Eigen library (http://eigen.tuxfamily.org/index.php?title=Main_Page) that provides us matrices, vectors, numerical solvers, and related algorithms.

# How to compile, install and link against GMDS
GMDS build system is CMake. In order to build GMDS, the direct way consists in the following stages

## Cloning the git repository
- Extracting GMDS from the gitub repository by cloning it (best option if you want to contribute afterward) or copying it. Let us assume that it is extracted in */home/gmds* by cloning the github repository.
```Shell
cd /home/
git clone https://github.com/LIHPC-Computational-Geometry/gmds.git
```
## Build directory for CMake purpose
In the */home/$ directory, we create a directory where the project will be built. Let us name it *build*. We will then get into this subdirectory to prepare the project to be built using CMake.
```Shell
cd /home/
mkdir build
cd build
cmake ../gmds/ 
```

## Compilation toolchain 
Depending on the options you provide to the cmake command line, you will get different way of building it. If you chose to create a Makefile, you will then have to build and (optionally) to install it
```Shell
make -j4
make install
```

## Example of toy code that uses gmds


# GMDS Project structure

GMDS project is split into a set of "*modules*", which are connected each to others. Each module is structured as follows:

- One subdirectory **inc**, which contains all the headers files of the module;
- One subdirectory **src**, which contains all the source files of the module;
- One subdirectory **tst**, which contains a set of unit tests written using google test. The aim of this **tst** subdirectory is both to (1) dynamically ensure the right and expected behaviour of the module and (2) provide some usage examples for other developers.
- A file named *CMakeLists.txt* used by CMAKE to build the module.

GMDS modules are of two kinds : *core* ones and *optional* ones

## Core modules
Core modules gather generic data structures, math classes, utils, basic algorithms and some Read and Write interfaces to be implemented. We have:

- [cad](cad.md): contains all the data structure and main query operations required for writing meshing algorithms.
- [math](math.md): mathematical simple concepts that are implemented and reused in many algorithms.
- [io](io.md): read and write classes with an interface definition to be implemented for new mesh data structure you would like to read and write.
- [ig](ig.md): serial generic mesh data structure based on an incidence-graph representation.

## Optional modules
Optional modules can be activated using the cmake option ENABLE_XXX where XXX is the module name. For instance for building the FRAME module, you have to add the option -DENABLE_FRAME=true:
```Shell
cmake ../gmds/ -DENABLE_FRAME=true
```
Available modules are :
- [frame](../../frame/README.md)
- [frame3d](../../frame3d/README.md)
- [kmds](kmds.md): concurrent thread-based mesh data structure. "k" means "Kokkos" the underlying framework we used for managing concurrency.
- [Elg3D](../../Elg3D/README.md)
- [hybridMeshAdapt](../../hybridMeshAdapt/README.md)
- [singGraphBuild](../../singGraphBuild/README.md)
- ...

The lists of all the available optional modules can be seen in the root CMakeLists.txt file or using the *ccmake* command instead of the *cmake* command.

