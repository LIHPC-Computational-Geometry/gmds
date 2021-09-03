/******************************************************************************/
#include <gmds/hybridMeshAdapt/PointInsertion.h>
#include <gmds/hybridMeshAdapt/CavityOperator.h>
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
#include <gmds/hybridMeshAdapt/ICriterion.h>
/******************************************************************************/
using namespace gmds;
using namespace hybrid;
using namespace operators;
using namespace simplicesCell;
using namespace simplicesNode;
/******************************************************************************/
PointInsertion::PointInsertion(SimplexMesh* simplexMesh, const SimplicesNode& simpliceNode, const CriterionRAIS& criterion, std::vector<TSimplexID> initialCavity, std::vector<TSimplexID> markedSimplex)
{
    TSimplexID simplexContainingNode = 0;
    if(simplexMesh != nullptr)
    {
      /*Si simpliceNode n'est pas a linterrieur de simplexMeshon ne fait rien*/
      if(simplexMesh->simplexContaining(simpliceNode, simplexContainingNode))
      {
        std::vector<TSimplexID> initCavity;
        if(initialCavity.size() == 0)
        {
          initCavity.push_back(simplexContainingNode);
        }
        else
        {
          initCavity = initialCavity;
        }

        initCavity.erase(std::remove_if(initCavity.begin(), initCavity.end(),
        [&](TSimplexID simplex)
        {
            bool flag = false;
            if(std::find(markedSimplex.begin(), markedSimplex.end(), simplex) != markedSimplex.end())
            {
              flag = true;
            }
            return flag;
        }), initCavity.end());
        CavityOperator cavOp(simplexMesh);

        //cavity first !!! and after point to connect
        CavityOperator::CavityIO&&  cavity = cavOp.cavity(initCavity, simpliceNode, criterion);

        std::vector<TInt>&               nodesInCav       = cavity.nodeInCavity();

        std::vector<std::vector<TInt>>&& pointsToConnect  = cavity.pointToConnect();

        //delete the tetra in the cavity
        for(auto const & simplexInCavity : cavity.simplexInCavityUnmodifiable())
        {
          if(std::find(markedSimplex.begin(), markedSimplex.end(), simplexInCavity) != markedSimplex.end())
          {
            std::cout << "simplexInCavity STOP " << std::endl;
          }
          simplexMesh->deleteTetra(simplexInCavity);
        }
        //delete the node in the cavity
        for(auto const & nodeInCavity : cavity.getNodeInCavity())
        {
            simplexMesh->deleteNode(nodeInCavity);
        }

        for(auto const & points2Connect : pointsToConnect)
        {
          simplexMesh->addTetraedre(points2Connect[0],
                                    points2Connect[1],
                                    points2Connect[2],
                                    simpliceNode.getGlobalNode());
        }
      }
      else
      {
        std::cout << "SIMPLICE NODE : " << simpliceNode << "IS NOT ON ANY TETRAEDRE" << std::endl;
      }
    }
}
/******************************************************************************/
