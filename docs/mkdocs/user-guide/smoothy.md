# Smoothy module

The aim of this module is to gather some smoothing algorithms. Mesh smoothing algorithms consist in moving mesh nodes without changing mesh topology. We classify smoothing algorithms in two categories: 
1. Those that are aware of the geometric classification. They are the less academic ones, but those used in practice. Most of the time they provide functions to smooth nodes on curves, surfaces and inside the volumes.
2. The other ones that smooth all the nodes excepted a bunch of ones that are locked to their positions. They derive from academic papers that deal with such low constraints.

## How to use a smoothing algorithm? 
Like other meshing algorithms in **gmds**, smoothing algorithms provided an `isValid` method to check if the mesh we work
on is adapted to the algorithm and a `smooth` method to perform the algorithm. As an example, let us consider the simple 2D 
laplacian smoother without geometric classification (see the function `grid_2D_smooth_UC` in the test file [LaplacianSmootherTestSuite.h](../../../smoothy/tst/LaplacianSmootherTestSuite.h)). To call this 
smoothing algorithm, we have a next series of instruction to write:
```c++
smoothy::LaplacianSmoother2UC smoother(&m);
ASSERT_TRUE(smoother.isValid());
smoother.setNbIterations(1);
...
smoother.setNodes(to_move);
...
smoother.smooth();
```
First, we initialize the *smoother* object as an instance of `LaplacianSmoother2UC`. It works on mesh `m`. The `2UC` 
postfix means that the algorithm is 2D (`2`) and doesn't consider the classification (`UC`stands for *unclassified*).
We check if the mesh `m` is valid for the algorithm with `smoother.isValid()`and then we indicate the number of iterations of the algorithm and the nodes to smooth. The `setNodes` method takes as parameters
a vector of node ids.

For an algorithm that use the geometric classification, we also need to provide the linker object as usual. See 
for instance the test method `tet_in_cube` in the file [LaplacianSmootherTestSuite.h](../../../smoothy/tst/LaplacianSmootherTestSuite.h).
## Elliptic smoothing
The elliptic smoother provides an implementation of the ACM TOG paper **Foldover-free maps in 50 lines of code** (https://dl.acm.org/doi/abs/10.1145/3450626.3459847).
This algoriths allows to smooth a 2D mesh with some nodes that are locked. Example of usages are given in the unit tests.