#ifndef EDGE_INSERTION_H_
#define EDGE_INSERTION_H_
/******************************************************************************/
#include "gmds/hybridMeshAdapt/SimplexMesh.h"
#include "gmds/hybridMeshAdapt/SimplicesNode.h"
#include "gmds/hybridMeshAdapt/SimplicesCell.h"
#include <gmds/hybridMeshAdapt/ICriterion.h>
#include <gmds/hybridMeshAdapt/CavityOperator.h>
#include <gmds/hybridMeshAdapt/PointInsertion.h>
/******************************************************************************/
namespace gmds
{
  namespace hybrid
  {
    namespace operators
    {
      class EdgeInsertion
      {
      public:
        EdgeInsertion() = delete;

        EdgeInsertion(SimplexMesh* simplexMesh, const simplicesNode::SimplicesNode& simpliceNodeA, const simplicesNode::SimplicesNode& simpliceNodeB, const CriterionRAIS& criterionRAIS, std::vector<TInt> markedNodes = std::vector<TInt>{});

        bool insertionEdge(SimplexMesh* simplexMesh, const simplicesNode::SimplicesNode& simpliceNodeA, const simplicesNode::SimplicesNode& simpliceNodeB, bool& cavityRebuild);

        void fillNodeAdj(SimplexMesh* simplexMesh, const simplicesNode::SimplicesNode& simpliceNodeA, const simplicesNode::SimplicesNode& simpliceNodeB, std::vector<TInt>& nodes);

        ~EdgeInsertion(){};
      };
    }
  }
}

#endif //EDGE_INSERTION_H_
