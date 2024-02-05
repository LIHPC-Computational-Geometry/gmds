# GMDS - Generic Mesh Data and Services
**GMDS**, for **G**eneric **M**esh **D**ata & **S**ervices,  is a set of C++ library written to provide mesh data 
structures and algorithms to developers that intend to design meshing algorithms and build pipelines of algorithms.

The development of this library started a few years ago to provide a generic way of designing data structures 
representing unstructured 2D and 3D meshes. Such meshes are defined as collections of cells that are topologically 
connected. Cells can be:

- Nodes, or 0-dimensional cells (0-cells),
- Edges, or 1-dimensional cells (1-cells),
- Faces, or 2-dimensional cells (2-cells),
- Regions, or 3-dimensional cells (3-cells).

Depending on the meshing algorithm a developer has to write, he must decide which cells are mandatory and which 
topological connections. Indeed, for an algorithm you may need to store edges and/or faces with the relation from nodes 
to edges and to faces and vice-versa, while for another you may require regions and the topological relation from faces 
to regions. GMDS provide flexible mechanisms to handle a huge variety of models and type of cells (triangles, 
quadrilaterals, tetrahedra, ...).

As we are mainly concerned about structured meshes in our team, most of proposed algorithms are dedicated to 
quadrilateral and hexahedral meshes.

## A generic mesh data structure
The data structures we provide are based onto the definition of a mesh model, which describes the available cells and 
connections in the mesh. For instance, in the next example code (line 2), a model is defined with flags `DIM3|F|E|N|F2N|E2N|N2F`. 
It means that the mesh is a 3D one (`DIM3`) made of faces (`F`), edges (`E`) and nodes (`N`), and connections from 
faces to nodes (`F2N`), edges to nodes (`E2N`) and nodes to faces (`N2F`) are explicitly stored.

Then we define a mesh with this model (line 3). We read nodes and faces of a vtk file in the block of instructions going
from lines 5 to 9 (see the [io module user guide](user-guide/io.md) for more details). Line 11 to 14 show how we can 
iterate on the mesh nodes. 

```cpp linenums="1"
using namespace gmds;
MeshModel model(DIM3|F|E|N|F2N|E2N|N2F); 
Mesh m(model);

IGMeshIOService ioService(&m);
VTKReader vtkReader(&ioService); 
vtkReader.setCellOptions(N|F);
vtkReader.setDataOptions(N|F);
vtkReader.read("B21.vtk");

for(auto node_id: m.nodes()){
    Node n = m.get<Node>(node_id);
    Point p = n.point();
}
```
## Algorithms as services

*gmds* provides data structures, io services and math functionalities for developing meshing algorithms. 

### Object-oriented approach for meshing algorithms
algorithms are gathered in modules, and we apply
an object-oriented approach where each *main algorithm* is handled by a single class instance. For instance, in the 
next code, an instance of the `CrossFieldGeneration2D` class is created. It acts on a mesh `m`. We run the algorithm
with the `execute` method.

```cpp 
CrossFieldGeneration2D field_generator(&m);
field_generator.execute(CrossFieldGeneration2D::laplace_solve);
```

### Main algorithms 
Our current interest is about structured quadrilateral and hexahedral meshing. To generate such meshes, we focus on 
several *technologies* including the following ones:

- **Frame Fields**. 3 modules are currenly dedicated to the usage of frame fields fo meshing
    - The [frame](user-guide/frame.md) module provides algorithms for 2D meshing. It relies on the notion of cross 
fields (see the [math](math/README.md) component for cross definitions). Output of this module are 2D cross fields 
defined on an input simplex mesh.
    - The [singGraphBuild](user-guide/singgraphbuild.md) module provides algorithm to extract the **base complex** structure of a 2D frame field
    - The [frame3d](user-guide/frame3d.md) module provides algorithms for 3D frame field generation. Unlike the 2D case, we are not able to generate a full block structure but such fields are used to drive hybrid mesh generation and point generation algorithms.
- **[Overlay grids](user-guide/elg3d.md)** algorithms, where an object *O* is embedded into a regular grid that is progressively adapted to capture
the boundary features of *O*.
- [Sheet operations](user-guide/sheet.md). This module provides sheet operations for quad and hex meshes.


## How to use gmds
Two options are possible:
## Users and developers documentation

Documentation is under construction. we just start to write it. It is split between:
- [Users documentation](user-guide/users_doc.md), which is dedicated to people who want to use **gmds** as a set of libraries but do not expect to contribute to it.
- [Developers documentation](dev-guide/developers_doc.md), which is dedicated to developers who would like to create a new *gmds* module for instance. In particular, we explain the [git workflow](docs/mkd/git_workflow.md) that we adopted.
- Gitub pages are under construction and available [here](https://lihpc-computational-geometry.github.io/gmds).

## Installation
gmds depends on many external components. On linux systems, and macos, we rely on  [spack](https://spack.io/) 
for installing gmds depencies but also gmds. In a nutshell, spack allows you to install a set of libraries in a
specific directory. You can see it as an equivalent of *Python environment*. 


## Example Project using gmds
We provide a [blank project](https://github.com/LIHPC-Computational-Geometry/gmds/tree/main/docs/example) example 
showing how to use **gmds** as a library with [CMake](https://cmake.org). You can run it with the following command 
lines

```shell
 cd docs/example
 cmake -DCMAKE_PREFIX_PATH=<path_to_gmds_install_dir> .
 make
 ./examplegmds
```

Feel free to copy or fork this project as a way of starting a new personal project using **gmds**.

# Coding Guidelines and Tips
GMDS follows strict coding guidelines, please take a look here before submitting your pull requests. We also have a set of general coding tips on how to code a geometry processing research project.
