#ifndef MESH_TRANSFORMATION_H_
#define MESH_TRANSFORMATION_H_

/*----------------------------------------------------------------------------*/
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
#include <gmds/math/Matrix.h>
#include <gmds/utils/BitVector.h>
#include <gmds/utils/VariableManager.h>
#include <gmds/utils/Variable.h>
#include <gmds/utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
#include <Eigen/Sparse>

namespace gmds
{
  namespace hybrid
  {
    class SimplexMesh;

    class MeshTransformation
    {
      public:
        MeshTransformation(SimplexMesh* simplexMesh);

        void transformation();

        Eigen::VectorXd resolveSystem(const Eigen::VectorXd& constraints, const TInt wichDimension);

        void resolveSystemByIterativeGradient();

        Eigen::Matrix3d interpolateMetric(const math::Point& currentNode, const math::Point& nodeStart, const math::Point& nodeEnd, const Eigen::Matrix3d& Mi, const Eigen::Matrix3d& Mj);

      private:
        SimplexMesh* m_simplexMesh;
    };
  }
}

#endif //MESH_TRANSFORMATION_H_
