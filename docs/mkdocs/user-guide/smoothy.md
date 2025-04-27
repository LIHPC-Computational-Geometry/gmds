# Smoothy module

The aim of this module is to gather some smoothing algorithms. Mesh smoothing algorithms consist in moving mesh nodes without changing mesh topology. We classify smoothing algorithms in two categories: 
1. Those that are aware of the geometric classification. They are the less academic ones, but those used in practice. Most of the time they provide functions to smooth nodes on curves, surfaces and inside the volumes. 
2. The other ones that smooth all the nodes excepted a bunch of ones that are locked to their positions. They derive from academic papers that deal with such low constraints.

Every algorithm is defined in a gmds class whose the naming indicates the type of mesh it acts on. More specifically, a class name that ends with **2C** or **3C** respectively means that the algorithm acts on a 2D classified mesh (**2C**) or a 3D classifed mesh (**3C**). If the postfix is **2UC** or **3UC**, it means that the algorithms acts on an unclassified mesh (**UC**).

## How to use a smoothing algorithm? 

Like other meshing algorithms in **gmds**, smoothing algorithms provided an `isValid` method to check if the mesh we work
on is adapted to the algorithm and a `smooth` method to perform the algorithm. As an example, let us consider the simple 2D laplacian smoother without geometric classification (see the function `grid_2D_smooth_UC` in the test file [LaplacianSmootherTestSuite.h](../../../core/smoothy/tst/LaplacianSmootherTestSuite.h)). To call this 
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

First, we initialize the *smoother* object as an instance of `LaplacianSmoother2UC`.  
We check if the mesh `m` is valid for the algorithm with `smoother.isValid()`and then we indicate the number of iterations (`setNBIterations`) of the algorithm and the nodes to smooth (`setNodes`). The `setNodes` method takes as parameters a vector of node ids.

For an algorithm that use the geometric classification, we also need to provide the linker object as usual. See 
for instance the test method `tet_in_cube` in the file [LaplacianSmootherTestSuite.h](../../../core/smoothy/tst/LaplacianSmootherTestSuite.h).

## Elliptic smoothing

The elliptic smoother provides an implementation of the ACM TOG paper **Foldover-free maps in 50 lines of code** (https://dl.acm.org/doi/abs/10.1145/3450626.3459847).
This algoriths allows to smooth a 2D mesh with some nodes that are locked. Example of usages are given in the unit tests.
