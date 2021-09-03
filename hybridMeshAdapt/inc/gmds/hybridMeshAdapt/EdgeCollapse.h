#ifndef EDGE_COLLAPSE_H_
#define EDGE_COLLAPSE_H_
/******************************************************************************/
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
#include <gmds/hybridMeshAdapt/SimplicesNode.h>
#include <gmds/hybridMeshAdapt/SimplicesCell.h>
/******************************************************************************/

namespace gmds
{
  namespace hybrid
  {
    namespace operators
    {
      class CriterionRAIS;

      class EdgeCollapse
      {
      public:

        EdgeCollapse() = delete;

        EdgeCollapse(SimplexMesh* simplexMesh, const simplicesNode::SimplicesNode& simpliceNodeA, const simplicesNode::SimplicesNode& simpliceNodeB, const CriterionRAIS& criterion);

        ~EdgeCollapse(){};

      };
    }
  }
}

#endif //EDGE_COLLAPSE_H_
