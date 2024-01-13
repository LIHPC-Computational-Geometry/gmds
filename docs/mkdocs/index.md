# GMDS Documentation  

# Last changes
- We currently change our installation procedure. In particular, we externalize some of gmds depencies. See the [developer documentation](docs/mkd/developers_doc.md) for more details.
- Important work is done on the Python API, see the the [pygmds](pygmds/README.md) module.
- Our blocking structure is evolving to be more robust and efficient for the blocking procedure. See the [blocking](blocking/README.md) module.

## GMDS in a nuthsell
**GMDS**, for **G**eneric **M**esh **D**ata & **S**ervices,  is a C++ library written to provide mesh data structures and algorithms to developers that intend to design meshing algorithms and build pipelines of those algorithms.

The development of this library started a few years ago to provide a generic way of designing data structures representing unstructured 2D and 3D meshes. Such meshes are defined as collections of cells that are topologically connected. Cells can be:
- Nodes, or 0-dimensional cells (0-cells)
- Edges, or 1-dimensional cells (1-cells)
- Faces, or 2-dimensional cells (2-cells)
- Regions, or 3-dimensional cells (3-cells)

Depending on the meshing algorithm a developer has to write, he must decide which cells are mandatory and which topological connections. Indeed, for an algorithm you may need to store edges and/or faces with the relation from nodes to edges and to faces and vice-versa, while for another you may require regions and the topological relation from faces to regions. GMDS provide flexible mechanisms to handle a huge variety of models and type of cells (triangles, quadrilaterals, tetrahedra, ...).

As we are mainly concerned about structured meshes in our team, most of proposed algorithms are dedicated to quadrilateral and hexahedral meshes.
### A Generic mesh data structure

The data structures we provide are based onto the definition of a mesh model, which describes the available cells and connections in the mesh. For instance, in the next example code, a model is defined with flags DIM3|F|E|N|F2N|E2N|N2F. It means that the mesh is a 3D one (DIM3) made of faces (F), edges (E) and nodes (N) and connections from faces to nodes (F2N), edges to nodes (E2N) and nodes to faces (N2F) are explicitly stored.
```cpp
MeshModel model(DIM3|F|E|N|F2N|E2N|N2F);
Mesh m(model);
// load a file into m

for(auto node_id: m->nodes()){
    Node n = m.get<Node>(node_id);
    Point p = n.point();
}
```
*gmds* provides a framework for developing new algorithms and our current interest is about structured quadrilateral and hexahedral meshing. To generate such meshes, we focus on the following *technologies*:
- **Frame Fields**. 3 modules are currenly dedicated to the usage of frame fields fo meshing
    - the [frame](frame.md) module provides algorithms for 2D meshing. It relies on the notion of cross fields (see the [math](math/README.md) component for cross definitions). Output of this module are 2D cross fields defined on an input simplex mesh.
    - the [singGraphBuild](singGraphBuild/README.md) module provides algorithm to extract the **base complex** structure of a 2D frame field
    - the [frame3d](frame/README.md) module provides algorithms for 3D frame field generation. Unlike the 2D case, we are not able to generate a full block structure but such fields are used to drive hybrid mesh generation and point generation algorithms.
- [Overlay grids algorithms](Elg3D/README.md).
- [Sheet operations](sheet/README.md). This module provides sheet operations for quad and hex meshes.
### A service-based approach
In order to build and prototype secure pipeline algorithms, we propose a **service** module to assemble our algorithms into a verified and dynamically-secured pipeline. We strongly believe that a main drawback of research but also production codes is that they're are written by researchers in mathematics, physics or computer science who focuses on the application "business" without taking care of "software engineering". This is quite usual and understandable but such a behaviour has 2 main consequences:
1. Codes are not
### Algorithms for quad and hex meshing
GMDS is our ma robust enough because the specification of their input and output are often fuzzy, only known by the main developer;
2. It is difficult to reuse an algorithm written by someone else.

The **service** module is an answer to this issue. Input and output of each service must be totally specified using a constraint system that is dynamically checked at execution time.

## Users and developers documentation

Documentation is under construction. we just start to write it. It is split between:
- [Users documentation](docs/mkd/users_doc.md), which is dedicated to people who want to use **gmds** as a set of libraries but do not expect to contribute to it.
- [Developers documentation](docs/mkd/developers_doc.md), which is dedicated to developers who would like to create a new *gmds* module for instance. In particular, we explain the [git workflow](docs/mkd/git_workflow.md) that we adopted.
- Gitub pages are under construction and available [here](https://lihpc-computational-geometry.github.io/gmds).
- The asso


# Coding Guidelines and Tips
GMDS follows strict coding guidelines, please take a look here before submitting your pull requests. We also have a set of general coding tips on how to code a geometry processing research project.
# Built with MkDocs
This documentation uses MkDocs. For full documentation visit [mkdocs.org](https://www.mkdocs.org).
