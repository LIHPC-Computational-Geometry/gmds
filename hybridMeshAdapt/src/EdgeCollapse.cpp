/******************************************************************************/
#include <gmds/hybridMeshAdapt/EdgeCollapse.h>
#include <gmds/hybridMeshAdapt/CavityOperator.h>
#include <gmds/hybridMeshAdapt/ICriterion.h>
/******************************************************************************/
using namespace gmds;
using namespace hybrid;
using namespace operators;
using namespace simplicesNode;
/******************************************************************************/
EdgeCollapse::EdgeCollapse(SimplexMesh* simplexMesh, const simplicesNode::SimplicesNode& simpliceNodeA, const simplicesNode::SimplicesNode& simpliceNodeB, const CriterionRAIS& criterion)
{
  if(simplexMesh != nullptr)
  {
    int boundary = std::numeric_limits<int>::min();
    std::vector<TSimplexID>&& initCavity = simpliceNodeB.ballOf();

    if(simpliceNodeB.shell(simpliceNodeA).size()  > 0)
    {
      if(std::find(initCavity.begin(), initCavity.end(), boundary) == initCavity.end())
      {
        CavityOperator cavOp(simplexMesh);

        //cavity first !!! and after point to connect
        CavityOperator::CavityIO&&  cavity = cavOp.cavity(initCavity, simpliceNodeA, criterion);
        std::vector<TInt>&               nodesInCav       = cavity.nodeInCavity();
        std::vector<std::vector<TInt>>&& pointsToConnect  = cavity.pointToConnect();

        const std::vector<TInt> & cavv = cavity.simplexInCavityUnmodifiable();
        for(auto const & simplexInCavity : cavity.simplexInCavityUnmodifiable())
        {
            simplexMesh->deleteTetra(simplexInCavity);

        }
        //delete the node in the cavity
        for(auto const & nodeInCavity : cavity.getNodeInCavity())
        {
          simplexMesh->deleteNode(nodeInCavity);
        }



        for(auto const & points2Connect : pointsToConnect)
        {
          if(points2Connect[0] != simpliceNodeA.getGlobalNode() &&
             points2Connect[1] != simpliceNodeA.getGlobalNode() &&
             points2Connect[2] != simpliceNodeA.getGlobalNode())
             {
               simplexMesh->addTetraedre(points2Connect[0],
                                         points2Connect[1],
                                         points2Connect[2],
                                         simpliceNodeA.getGlobalNode());
             }
        }
      }

    }
    else
    {
      /*TODO*/
      /*mettre une exeption ... A et B ne sont pas conjoint..*/
      std::cout << "SIMPLICE NODE : " << simpliceNodeA << "IS NOT ADJACENT TO " << simpliceNodeB << std::endl;
    }
  }
}
