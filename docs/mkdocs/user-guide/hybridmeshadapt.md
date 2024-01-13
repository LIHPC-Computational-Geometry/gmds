# HybridMeshAdapt module

This module proposes an algorithm to adapt a tetrahedral mesh into an hybrid one using frame fiels and metrics.

To achieve a sufficient level of performances, a mesh data structure dedicated to simplicial meshes is used in this
component instead of those proposed by the [ig](../ig/README.md) and [kmds](../kmds/README.md) modules.

## Stage 1
The first executable to use is main_METRIC_FF_PointGeneration : The frontal process is subdividide in 3 step
    The node generation on the curve, the surface and the volume. The function nodesSpreading is used to spread the node
    on the surface and then on the volume. The result of this algorithm is a mesh containing the node and the related hexaedron.

    * metricXYZ_functors is a vector of function containing the metric informations
    * frameXYZ_functors is a vector of function containing the frame informations (if the function return a nul size vector, the frame field of the mesh is used.)

    Once The two vector have been specified the constructor MetricFFPointgeneration p(&simplexMesh, name_final_mesh, interpolation_factor) is called,
    interpolation_factor is the value of the interpolation between the frame and the metric orientation. The frontal process will generate node in the direction of the
    frame field when interpolation_factor tend to 1. If interpolation_factor tend to 0 the node will follow the directions of the metric directions, (see the function MetricFFPointgeneration::metricInterpolationWithDistorsion(const Eigen::Matrix3d  & metric, const Eigen::Matrix3d  & frameField)).


## Stage 2
Once the previous process is done, we call main_HEX_GENERATION with the initial mesh and the result mesh of the main_METRIC_FF_PointGeneration.
    * First, we perform a Delaunay node insertion, then the node that has not been inserted are inserted using only the simplex containing the current node (l:131-247)
    * Second, we remove the initial node of the mesh, this step contribute to build the shell of the hexaedron. (l:267-329)
    * Third, we build the remaining edge between the inserted node (in order to build the shell of the future hexaedron). (l:332-368)
    * To finish we force the faces of the futur hexaedrons. If the 6 faces exist the hexaedron is built, and the tetraedron forming the hex are marked in order to not be built
      in the final mesh (l:378-454)


## Stage 3
The last step will build the hexaedron from the remaining tetradron, the executable main_TET2HEX is called to perform this.
