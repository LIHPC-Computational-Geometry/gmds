#ifndef POINT_INSERTION_H_
#define POINT_INSERTION_H_
/******************************************************************/
#include "gmds/hybridMeshAdapt/SimplicesNode.h"
#include "gmds/hybridMeshAdapt/SimplicesCell.h"
#include <gmds/hybridMeshAdapt/ISimplexMeshIOService.h>
/******************************************************************/
#include <gmds/io/VTKWriter.h>
/******************************************************************/
#include <map>
/******************************************************************/
namespace gmds
{
  namespace hybrid
  {
    class SimplexMesh;

    namespace operators
    {
          class CriterionRAIS;

          struct CavityReduction
          {
            std::vector<TSimplexID> cavInitial{};
          };

          class PointInsertion
          {
            public :
              PointInsertion();

              PointInsertion(SimplexMesh* simplexMesh, const simplicesNode::SimplicesNode& simpliceNode, const CriterionRAIS& criterion, bool& status,const  std::vector<TSimplexID>& initialCavity,const gmds::BitVector& markedNodes,
                             std::vector<TInt>& deletedNodes, const std::multimap<TInt, std::pair<TInt, TInt>>& facesAlreadyBuilt, std::vector<TSimplexID> markedSimplex = {});

              ~PointInsertion(){};
          };
    }
  }
}

#endif // POINT_INSERTION_H_
