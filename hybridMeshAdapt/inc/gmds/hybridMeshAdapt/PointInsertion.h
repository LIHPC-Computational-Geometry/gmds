#ifndef POINT_INSERTION_H_
#define POINT_INSERTION_H_
/******************************************************************/
#include <gmds/hybridMeshAdapt/CavityOperator.h>
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
                         std::vector<TInt>& deletedNodes, const std::multimap<TInt, TInt>& facesAlreadyBuilt, std::vector<TSimplexID>& createdCells, std::vector<TSimplexID> markedSimplex = {});

          ~PointInsertion(){};

          void facesIdxPatch(SimplexMesh* simplexMesh, CavityOperator::CavityIO& cavityIO, std::vector<std::vector<TInt>>& facesIdx, const simplicesNode::SimplicesNode & node) const;

          void normalsPatch(SimplexMesh* simplexMesh, CavityOperator::CavityIO& cavityIO, std::vector<math::Vector3d>& normals, const simplicesNode::SimplicesNode & node) const;

      };
    }
  }
}

#endif // POINT_INSERTION_H_
