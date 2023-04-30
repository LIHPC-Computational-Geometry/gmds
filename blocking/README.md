# Blocking module

This module provides a data structure for a blocking representation. It relies on
the *n*-G-map model and we specifically use the 3-G-Map implementation 
provided by CGAL library (see https://doc.cgal.org/latest/Generalized_map/index.html).

The development of this module is in progress. The final purpose is to provide a robust and reliable curved block 
structure for hex blocking mesh representation. 

## Compilation and installation

This module is built using the global `CMakeLists.txt`of **gmds**. An extra option (`WITH_CGNS`) is provided to allow
CGNS export. It requires then to install the adequate CGNS library. 

*This procedure must be described in details and the spack recipes must be updated.*


## Main concepts

The blocking structure is described in  the class [CurvedBlocking](inc/CurvedBlocking.h). It relies on the 3-G-map model
and provide several blocking queries and gives access to the underlying 3-G-map.

The [CurvedBlockingClassifier](inc/CurvedBlockingClassifier.h)provides the first prototype to classify a block structure
onto a geometrical model. Up to now, only the node classification is available. It classify nodes of blocks on 
geometrical points, curves and surfaces using two parameters: 
- a **maximal distance** *MD* for projecting a node. If the distance to any geometrical cell is greater than *MD*, then the node is not moved and classified.
- A **snapping distance** *SD* that is used during a correction stage. After projecting a node on a curve or a surface, we check if this node could not be snapped onto a geometrical point with a distance lower to *SD*.

**Remark:** This two parameters are use-case specific and should be automatically computed in the future. Is it possible