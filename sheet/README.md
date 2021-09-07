# Sheet module description

In this module, we provide the usual operations that can 
be safely performed on quad or hex mesh. By safely, we
means that the mesh remains a full quad, respectively hex
mesh after performing a sheet operation.

We provide the 3 basic operations, which are:
- **Sheet selection**. Giving 2 nodes, or an edge, it returns 
the set of cells that belongs to the corresponding sheet.
- **Sheet collapse**. Giving 2 nodes, or an edge, it removes
 the corresponding sheet from the mesh.
- **Sheet pillow**. Giving a path of edges or a convex set of
cells, it inserts a sheet inside the mesh.


## Requirements
Sheet collapse and pillow are operations that modify the
mesh. They are quite easy to implement if you only consider
the mesh structure but taking geometric classification 
makes it more tricky to handle.

Geometric classification requires to handle a few special 
cases where faces (in 3D) and edges (in 2D) are required
along the mesh boundary. It means that we restrict our
implementation to mesh where:
- in 2D, boundary edges must be defined in the mesh and connection E2N at least. 
We do not use the N2E connectivity because it enforces all
the mesh nodes to store an indirection to a potential list of adjacent edges.
- in 23D, boundary faces and edges must be defined in the mesh and connections F2N and E2N at least