# A C++ library for writing meshing algorithms

[![GitHub issues](https://img.shields.io/github/issues/LIHPC-Computational-Geometry/gmds)](https://github.com/LIHPC-Computational-Geometry/gmds/issues)
[![GitHub license](https://img.shields.io/github/license/LIHPC-Computational-Geometry/gmds)](https://github.com/LIHPC-Computational-Geometry/gmds/blob/main/LICENSE)

**GMDS**, for **G**eneric **M**esh **D**ata & **S**ervices,  is a C++ library written to provide mesh data structures and algorithms to developers that intend to design meshing algorithms and build pipelines of those algorithms.

The development of this library started a few years ago to provide a generic way of designing data structures representing unstructured 2D and 3D meshes. Such meshes are defined as collections of cells that are topologically connected. Cells can be:
- Nodes, or 0-dimensional cells (0-cells)
- Edges, or 1-dimensional cells (1-cells)
- Faces, or 2-dimensional cells (2-cells)
- Regions, or 3-dimensional cells (3-cells)

Depending on the meshing algorithm a developer has to write, he must decide which cells are mandatory and which topological connections. Indeed, for an algorithm you may need to store edges and/or faces with the relation from nodes to edges and to faces and vice-versa, while for another you may require regions and the topological relation from faces to regions. GMDS provide flexible mechanisms to handle a huge variety of models and type of cells (triangles, quadrilaterals, tetrahedra, ...).

As we are mainly concerned about structured meshes in our team, most of proposed algorithms are dedicated to quadrilateral and hexahedral meshes. 

For more information, get to
the associated [pages]( https://lihpc-computational-geometry.github.io/gmds/). 