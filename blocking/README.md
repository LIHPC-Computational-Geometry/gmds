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

