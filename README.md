# GMDS: a C++ library for writing meshing algorithms

![spack-ci](https://github.com/LIHPC-Computational-Geometry/qqualif/actions/workflows/spack-ci.yml/badge.svg)

![CI Ubuntu](https://github.com//LIHPC-Computational-Geometry/gmds/actions/workflows/continuous-ubuntu.yml/badge.svg)
![CI Macoc](https://github.com//LIHPC-Computational-Geometry/gmds/actions/workflows/continuous-macos.yml/badge.svg)
![CI Windows](https://github.com//LIHPC-Computational-Geometry/gmds/actions/workflows/continuous-windows.yml/badge.svg)
![CI Doxygen](https://github.com//LIHPC-Computational-Geometry/gmds/actions/workflows/code_docs.yml/badge.svg)
![CI_GH](https://github.com//LIHPC-Computational-Geometry/gmds/actions/workflows/deploy-gh.yml/badge.svg)

[![GitHub issues](https://img.shields.io/github/issues/LIHPC-Computational-Geometry/gmds)](https://github.com/LIHPC-Computational-Geometry/gmds/issues)
[![GitHub license](https://img.shields.io/github/license/LIHPC-Computational-Geometry/gmds)](https://github.com/LIHPC-Computational-Geometry/gmds/blob/main/LICENSE)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/bf6ad23f6a3c452ab6f7f3f63d9fdb89)](https://www.codacy.com/gh/LIHPC-Computational-Geometry/gmds/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=LIHPC-Computational-Geometry/gmds&amp;utm_campaign=Badge_Grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/bf6ad23f6a3c452ab6f7f3f63d9fdb89)](https://www.codacy.com/gh/LIHPC-Computational-Geometry/gmds/dashboard?utm_source=github.com&utm_medium=referral&utm_content=LIHPC-Computational-Geometry/gmds&utm_campaign=Badge_Coverage)
[![codecov](https://codecov.io/gh/LIHPC-Computational-Geometry/gmds/branch/main/graph/badge.svg?token=QA3AS0MLDN)](https://codecov.io/gh/LIHPC-Computational-Geometry/gmds)

This project is part of the [magix3d](https://github.com/LIHPC-Computational-Geometry/magix3d) ecosystem and conforms to its [CI policy](https://github.com/LIHPC-Computational-Geometry/spack_recipes#development-in-magix3d-ecosystem-projects).

## Last changes
- We changed our installation procedure. In particular, we externalized some of gmds depencies. See the [developer documentation](docs/mkdocs/dev-guide/developers_doc.md) for more details.
- Important work is done on the Python API, see the [pygmds](pygmds/README.md) module.
- Our blocking structure is evolving to be more robust and efficient for the blocking procedure. See the [blocking](core/blocking/README.md) module.

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
    - the [frame](docs/mkdocs/user-guide/frame.md) module provides algorithms for 2D meshing. It relies on the notion of cross fields (see the [math](core/math/README.md) component for cross definitions). Output of this module are 2D cross fields defined on an input simplex mesh.
    - the [singGraphBuild](docs/mkdocs/user-guide/singgraphbuild.md) module provides algorithm to extract the **base complex** structure of a 2D frame field
    - the [frame3d](docs/mkdocs/user-guide/frame3d.md) module provides algorithms for 3D frame field generation. Unlike the 2D case, we are not able to generate a full block structure but such fields are used to drive hybrid mesh generation and point generation algorithms. 
- [Overlay grids algorithms](docs/mkdocs/user-guide/elg3d.md).
- [Sheet operations](docs/mkdocs/user-guide/sheet.md). This module provides sheet operations for quad and hex meshes.

## Users and developers documentation

Documentation is under construction, you can see it [here](https://lihpc-computational-geometry.github.io/gmds). The associated doygen documenation is [here](https://lihpc-computational-geometry.github.io/gmds/doxygen/index.html).
