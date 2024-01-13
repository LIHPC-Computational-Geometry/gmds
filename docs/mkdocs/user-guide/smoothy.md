# Smoothy module

The aim of this module is to gather some smoothing algorithms. Mesh smoothing algorithms consist in moving mesh nodes without changing mesh topology. We classify smoothing algorithms in two categories: 
1. Those that are aware of the geometric classification. They are the less academic ones, but those used in practice. Most of the time they provide functions to smooth nodes on curves, surfaces and inside the volumes.
2. The other ones that smooth all the nodes excepted a bunch of ones that are locked to their positions. They derive from academic papers that deal with such low constraints.

## Elliptic smoothing
The elliptic smoother provides an implementation of the ACM TOG paper **Foldover-free maps in 50 lines of code** (https://dl.acm.org/doi/abs/10.1145/3450626.3459847).
This algoriths allows to smooth a 2D mesh with some nodes that are locked. Example of usages are given in the unit tests.