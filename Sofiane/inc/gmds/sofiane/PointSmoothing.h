#ifndef POINTSMOOTHING_H_
#define POINTSMOOTHING_H_
/******************************************************************/
#include "gmds/sofiane/SimplexMesh.h"
#include "gmds/sofiane/SimplicesNode.h"
#include "gmds/sofiane/SimplicesCell.h"
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
