/******************************************************************/
#include <gmds/hybridMeshAdapt/PointSmoothing.h>
#include <gmds/hybridMeshAdapt/CavityOperator.h>
#include <gmds/hybridMeshAdapt/ICriterion.h>
/******************************************************************/
using namespace gmds;
using namespace hybrid;
using namespace operators;
using namespace simplicesCell;
using namespace simplicesNode;
/******************************************************************/
PointSmoothing::PointSmoothing(hybrid::SimplexMesh* simplexMesh, const simplicesNode::SimplicesNode& previousNode, const simplicesNode::SimplicesNode& nextNode, const CriterionRAIS& criterion)
{

  std::vector<TSimplexID>&& initCavity = previousNode.ballOf();
  CavityOperator cavOp(simplexMesh);
  //cavity first !!! and after point to connect
  CavityOperator::CavityIO&&       cavity          = cavOp.cavity(initCavity, previousNode, nextNode, criterion);
  std::vector<TInt>&              nodesInCav       = cavity.nodeInCavity();


  //On ajoute le node previous si il est pas contenant dans la cavité (pour des cas comme un node a deplacé alors qu'il est en bords de mesh)
  if(std::find(cavity.getNodeInCavity().begin(), cavity.getNodeInCavity().end(), previousNode.getGlobalNode()) == cavity.getNodeInCavity().end())
  {
    cavity.addNodeInCavity(previousNode.getGlobalNode());
  }
  std::vector<std::vector<TInt>>&& pointsToConnect  = cavity.pointToConnect();

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

    TSimplexID tetraCreated= simplexMesh->addTetraedre(points2Connect[0],
                                 points2Connect[1],
                                 points2Connect[2],
                                 nextNode.getGlobalNode());
  }
}
/******************************************************************************/
