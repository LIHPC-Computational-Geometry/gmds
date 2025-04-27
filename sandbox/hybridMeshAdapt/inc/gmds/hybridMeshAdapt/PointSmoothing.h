#ifndef POINTSMOOTHING_H_
#define POINTSMOOTHING_H_
/******************************************************************/
#include "gmds/hybridMeshAdapt/SimplexMesh.h"
#include "gmds/hybridMeshAdapt/SimplicesNode.h"
#include "gmds/hybridMeshAdapt/SimplicesCell.h"
/******************************************************************/
namespace gmds
{

    namespace hybrid
    {

      namespace operators
      {
        class CriterionRAIS;

        class PointSmoothing
        {

        public:
          PointSmoothing() = delete;

          PointSmoothing(hybrid::SimplexMesh* simplexMesh, const simplicesNode::SimplicesNode& previousNode, const simplicesNode::SimplicesNode& nextNode, const CriterionRAIS& criterion);

          ~PointSmoothing(){};
        };
      }
    }
}


#endif //POINTSMOOTHING_H_
