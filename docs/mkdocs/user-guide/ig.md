# IG module

**IG**, for incidence graph, is the main and legacy mesh data structure provided in *gmds*.

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

## Block-structures meshes
In order to handle block structures meshes, we complete the **IG** representation with a "quasi-structured" or "block-structured" view
built on the traditional **IG** model (see the class [Mesh](inc/gmds/Mesh.h))

The blocking2D structure is a full unstructured mesh made of hex blocks, quad faces, edges and nodes. Block, faces and edges only exist at the macro level, or block level. Nodes are block corners 
but also inner-block nodes. Each node knows it is a corner-block (block level) or an inner-node
using the node on-node variable called "embedding". An embedding is a Cell::Data that corresponds to the block entity the node is
embedded in. It is made of a dim and an id. 
It can be on a corner-block(0), face-block(2), edge-block (1) and inner-block(3). Each block as a grid structure. It must be so possible to use bracket notation [i+,j-1] to acces to a neighbor nodes. This 
type of traversal should be possible almost everywhere and allows the user to traverse several blocks in a row. The only issue is when you meet singular corners (ie with valence not equal to 4).

