#ifndef POINT_INSERTION_H_
#define POINT_INSERTION_H_
/******************************************************************/
#include "gmds/hybridMeshAdapt/CavityOperator.h"
#include "gmds/hybridMeshAdapt/SimplicesNode.h"
#include "gmds/hybridMeshAdapt/SimplicesCell.h"
/******************************************************************/
namespace gmds
{
  namespace hybrid
  {
    class SimplexMesh;

    namespace operators
    {
          class CriterionRAIS;

          class PointInsertion
          {
            public :

              PointInsertion() = delete;

              PointInsertion(SimplexMesh* simplexMesh, const simplicesNode::SimplicesNode& simpliceNode, const CriterionRAIS& criterion, std::vector<TSimplexID> initialCavity = std::vector<TSimplexID>{}, std::vector<TSimplexID> markedSimplex = std::vector<TSimplexID>{});

              ~PointInsertion(){};
          };
    }
  }
}

#endif // POINT_INSERTION_H_
