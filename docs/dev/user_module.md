### How to use gmds

**gmds** is a set of libraries that you can build and use separatly to your purposes. 

It requires features of C++14. The build system is based on CMAKE (https://cmake.org) and unit tests are performed using google tests and google benchmarks.

In order to perform some linear algebra operations, we rely onto the Eigen library (http://eigen.tuxfamily.org/index.php?title=Main_Page) that provides us matrices, vectors, numerical solvers, and related algorithms.

# How to compile and install GMDS
In order to build GMDS, the simple way consists in:
1. Extracting **gmds** from the github repository by cloning it (best option if you want to contribute afterward) or copying it. Let us assume that it is extracted in */home/gmds* by cloning the gitlab repository.
```Shell
cd /home/
git clone https://github.com/LIHPC-Computational-Geometry/gmds.git
```
2. In the */home/$ directory, we create a directory where the project will be built. Let us name it *build*. We will then get into this subdirectory to prepare the project to be built using CMake.
```Shell
cd /home/
mkdir build
cd build
cmake ../gmds/ 
```
3. Depending on the options you provide to the cmake command line, you will get different way of building it. If you chose to create a Makefile, you will then have to build and (optionally) to install it
```Shell
make -j4
make install
```
# Gmds project structure

**gmds** project is split into a set of "*modules*", which are connected each to others. Each module is structured as follows:
- One subdirectory **inc**, which contains all the headers files of the module;
- One subdirectory **src**, which contains all the source files of the module;
- One subdirectory **tst**, which contains a set of unit tests written using google test. The aim of this **tst** subdirectory is both to (1) dynamically ensure the right and expected behaviour of the module and (2) provide some usage examples for other developers.
- A file named *CMakeLists.txt* used by CMAKE to build the module.

**gmds** modules are of two kinds : *core* ones and *optional* ones

## Core modules
Core modules gather generic data structures, math classes, utils, basic algorithms and some Read and Write interfaces to be implemented. We have:
- **cad**: contains all the data structure and main query operations required for writing meshing algorithms.
- **math**: mathematical simple concepts that are implemented and reused in many algorithms.
- **io**: read and write classes with an interface definition to be implemented for new mesh data structure you would like to read and write.
- **ig**: serial generic mesh data structure based on an incidence-graph representation.
- **kmds**: concurrent thread-based mesh data structure. "k" means "Kokkos" the underlying framework we used for managing concurrency.

## Optional modules
Optional modules can be activated using the cmake option ENABLE_XXX where XXX is the module name. For instance for building the FRAME module, you have to add the option -DENABLE_FRAME=true:
```Shell
cmake ../gmds/ -DENABLE_FRAME=true
```
Available modules are :
- FRAME
- FRAME_3D
- DUAL_BLOCKING
- MEDUSA
- ...

The lists of all the available optional modules can be seen in the root CMakeLists.txt file or using the *ccmake* command instead of the *cmake* command.



# The **ig** module
The **ig** module provides basic data structure for handling meshes based on an "*Incident Graph*" representation (which explains why it is called *ig*). Traditionally a graph G=(V,A) is made of vertices and arcs that connect vertices. For giving a model to represent mesh data structures, we use such a graph model, where vertices are the dimensions of cells explicitly stored in the data structure. In other words, vertices are denoted:
1. **R** for regions or 3-dimensional cells (for short 3-cells),
2. **F** for faces or 2-cells,
3. **E** for edges or 1-cells,
4. **N** for nodes or 0-cells.

Connections, ie. graph arcs, are cell connections and are written **X2Y** with **X** and **Y** that are **R**, **F**, **E** or **N**. For instance **F2N** means that we store the connection from faces to nodes, while **N2R** indicates that we have to store the connection from nodes to regions. A *MeshModel* object can be built on this graph structure and given to a *Mesh* object to initialise it.
```cpp
MeshModel model(DIM3|F|E|N|F2N|E2N|N2F);
Mesh m(model);
```

# Utilitary modules **utils** and **math**