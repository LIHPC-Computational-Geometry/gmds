# IGAlgo module

This module provides some basic algorithms for the incidence graph (*ig*) data structure:

- Most important algorithms are relative to retrieving and building some boundary informations: boundaryOperator, boundaryExtractor.
- The *GridBuilder* class provides algoriths to build 2D and 3D structured grids for unit testing purposes mainly.
- The *THexBuilder* create a hexahedral mesh by splitting each tet of an input tetrahedral mesh.
- The *MeshQualityCompute* go through all the mesh cells and assign variables depending on a set of selected quality criteria
- The *SurfaceReorient* class reorient surface mesh (in 2D and 3D) in order to have all faces oriented in the same direction. If
the mesh is split in different connex components, the orientation between different components can be different.
In practice, the orientation is consistent in 2D.