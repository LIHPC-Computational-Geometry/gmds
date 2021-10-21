#IG module

**IG**, for incidence graph, is the main and historic mesh data structure provided in *gmds*.

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