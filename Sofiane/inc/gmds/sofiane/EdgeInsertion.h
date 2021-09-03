#ifndef EDGE_INSERTION_H_
#define EDGE_INSERTION_H_
/******************************************************************************/
#include "gmds/sofiane/SimplexMesh.h"
#include "gmds/sofiane/SimplicesNode.h"
#include "gmds/sofiane/SimplicesCell.h"
/******************************************************************************/
namespace gmds
{
  namespace hybrid
  {
    namespace operators
    {
      template<class Criterion>
      class EdgeInsertion
      {
      public:
        EdgeInsertion() = delete;

        EdgeInsertion(SimplexMesh* simplexMesh, const simplicesNode::SimplicesNode& simpliceNodeA, const simplicesNode::SimplicesNode& simpliceNodeB);

        ~EdgeInsertion(){};
      };
    }
  }
}


#include "../../../src/EdgeInsertion.cpp"

#endif //EDGE_INSERTION_H_
