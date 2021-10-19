# A C++ library for writing meshing algorithms

[![GitHub issues](https://img.shields.io/github/issues/LIHPC-Computational-Geometry/gmds)](https://github.com/LIHPC-Computational-Geometry/gmds/issues)
[![GitHub license](https://img.shields.io/github/license/LIHPC-Computational-Geometry/gmds)](https://github.com/LIHPC-Computational-Geometry/gmds/blob/main/LICENSE)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/bf6ad23f6a3c452ab6f7f3f63d9fdb89)](https://www.codacy.com/gh/LIHPC-Computational-Geometry/gmds/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=LIHPC-Computational-Geometry/gmds&amp;utm_campaign=Badge_Grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/bf6ad23f6a3c452ab6f7f3f63d9fdb89)](https://www.codacy.com/gh/LIHPC-Computational-Geometry/gmds/dashboard?utm_source=github.com&utm_medium=referral&utm_content=LIHPC-Computational-Geometry/gmds&utm_campaign=Badge_Coverage)
[![codecov](https://codecov.io/gh/LIHPC-Computational-Geometry/gmds/branch/main/graph/badge.svg?token=QA3AS0MLDN)](https://codecov.io/gh/LIHPC-Computational-Geometry/gmds)

**GMDS**, for **G**eneric **M**esh **D**ata & **S**ervices,  is a C++ library written to provide mesh data structures and algorithms to developers that intend to design meshing algorithms and build pipelines of those algorithms.

The development of this library started a few years ago to provide a generic way of designing data structures representing unstructured 2D and 3D meshes. Such meshes are defined as collections of cells that are topologically connected. Cells can be:
- Nodes, or 0-dimensional cells (0-cells)
- Edges, or 1-dimensional cells (1-cells)
- Faces, or 2-dimensional cells (2-cells)
- Regions, or 3-dimensional cells (3-cells)

Depending on the meshing algorithm a developer has to write, he must decide which cells are mandatory and which topological connections. Indeed, for an algorithm you may need to store edges and/or faces with the relation from nodes to edges and to faces and vice-versa, while for another you may require regions and the topological relation from faces to regions. GMDS provide flexible mechanisms to handle a huge variety of models and type of cells (triangles, quadrilaterals, tetrahedra, ...).

As we are mainly concerned about structured meshes in our team, most of proposed algorithms are dedicated to quadrilateral and hexahedral meshes. 

## A Generic mesh data structure

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
## Algorithms for quad and hex meshing
GMDS is our main framework for developing new algorithms and our current interest is about structured quadrilateral and hexahedral meshing. To generate such meshes, we focus on the following *technologies*:
- [Frame fields](https://gitlab.com/franck.ledoux/gmds/wikis/wiki_doc/Frame_Module). In particular, the **frame** module provides algorithms based on frame fields.
- [Overlay grids algorithms](https://gitlab.com/franck.ledoux/gmds/wikis/wiki_doc/OGrid_Module).
- [Sheet operations](sheet/README.md). This module provides sheet operations for quad and hex meshes.
## A service-based approach
In order to build and prototype secure pipeline algorithms, we propose a **service** module to assemble our algorithms into a verified and dynamically-secured pipeline. We strongly believe that a main drawback of research but also production codes is that they're are written by researchers in mathematics, physics or computer science who focuses on the application "business" without taking care of "software engineering". This is quite usual and understandable but such a behaviour has 2 main consequences:
1. Codes are not robust enough because the specification of their input and output are often fuzzy, only known by the main developer;
2. It is difficult to reuse an algorithm written by someone else.

The **service** module is an answer to this issue. Input and output of each service must be totally specified using a constraint system that is dynamically checked at execution time.

## User and developer documentation

All the documentation is currently available in the [GMDS Wiki](https://gitlab.com/franck.ledoux/gmds/wikis/GMDS-Wiki) that is modified asap.

### Git usage
In order to add a complete repository in the external repo, which contains external libraries use the *git add -f* option
