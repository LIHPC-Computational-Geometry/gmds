#ifndef DELAUNAY_POINT_INSERTION_H
#define DELAUNAY_POINT_INSERTION_H
/******************************************************************/
#include "gmds/hybridMeshAdapt/CavityOperator.h"
#include "gmds/hybridMeshAdapt/SimplicesNode.h"
#include "gmds/hybridMeshAdapt/SimplicesCell.h"
#include "gmds/hybridMeshAdapt/PointInsertion.h"
#include "gmds/hybridMeshAdapt/Metric.h"
/******************************************************************/
#include<gmds/math/Vector.h>
#include<gmds/math/Point.h>
#include<gmds/utils/CommonTypes.h>
/******************************************************************/
namespace gmds
{
  namespace hybrid
  {
    class SimplexMesh;

    namespace operators
    {
      class CriterionRAIS;

      class DelaunayPointInsertion
      {
      public:
        DelaunayPointInsertion() = delete;

        DelaunayPointInsertion(SimplexMesh* simplexMesh, const simplicesNode::SimplicesNode& simpliceNode, const CriterionRAIS& criterion, std::vector<TSimplexID>& initialCavity, bool& status, const gmds::BitVector& markedNodes /*= gmds::BitVector{}*/, std::vector<TSimplexID> markedSimplex = std::vector<TSimplexID>{});

        ~DelaunayPointInsertion(){};

        bool isNodeInCircumSphere(SimplexMesh* simplexMesh, const simplicesNode::SimplicesNode& node, const TSimplexID& simplex);
      };
    }
  }
}
#endif
