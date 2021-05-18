/******************************************************************************/
#include <chrono>
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
using namespace simplicesTriangle;
using namespace simplicesNode;
/******************************************************************************/
PointInsertion::PointInsertion()
{

}
/******************************************************************************/
PointInsertion::PointInsertion(SimplexMesh* simplexMesh, const SimplicesNode& simpliceNode, const CriterionRAIS& criterion, bool& status, const std::vector<TSimplexID>& initialCavity, const gmds::BitVector& markedNodes, std::vector<TSimplexID> markedSimplex, CavityReduction cavReduction)
{
    if(simplexMesh != nullptr)
    {
      int border = std::numeric_limits<int>::min();
      /*Si simpliceNode n'est pas a linterrieur de simplexMeshon ne fait rien*/
      std::vector<TSimplexID> initCavity;
      if(initialCavity.size() == 0)
      {
        SimplicesNode::nodeNeighborInfo nodeInfo;
        if(!simplexMesh->checkSimplicesContenaing(simpliceNode.getCoords(), initCavity))
        {
          std::cout << "SIMPLICE NODE : " << simpliceNode << "IS NOT ON ANY TETRAEDRE" << std::endl;
          return;
        }
      }
      else
      {
        initCavity = initialCavity;
      }

      bool flag = true;
      if(markedSimplex.size() != 0)
      {
        for(auto const & simplexInit : initCavity)
        {
          if(std::find(markedSimplex.begin(), markedSimplex.end(), simplexInit) != markedSimplex.end() )
          {
            flag = false;
          }
        }
      }

      if(flag && initCavity.size() != 0)
      {

        CavityOperator cavOp(simplexMesh);
        std::vector<TSimplexID> initialCavityCell{};
        std::vector<TSimplexID> initialCavityTriangle{};
        std::copy_if(initialCavity.begin(), initialCavity.end(), std::back_inserter(initialCavityCell), [&](const TSimplexID simplex){
          bool flag = false;
          if(simplex >= 0)
          {
            flag = true;
          }
          else
          {
            initialCavityTriangle.push_back(simplex);
          }
          return flag;
        });

        CavityOperator::CavityIO cavityIO(simplexMesh);
        if(cavReduction.cavInitial.size() == 0)
        {
          cavOp.cavityEnlargement(cavityIO, initialCavityCell, initialCavityTriangle, simpliceNode, criterion, markedSimplex);
        }
        else
        {
          cavOp.cavityReduction(cavityIO, initCavity, simpliceNode, criterion, cavReduction, markedSimplex);
        }
        ////////////////////////////////////////////////////////////////////////////////
        ///////////////////////finding the node inside the cavity///////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        const std::vector<TInt>& nodesInsideCavity = cavityIO.nodeInCavity();
        cavityIO.nodesReconnection();

        //check if node in nodsIncavity are marked
        for(auto const & nodeInCavity : cavityIO.getNodeInCavity())
        {
          if(markedNodes[nodeInCavity] == 1)
          {
            status = false;
            return;
          }
        }

        //deleteThe simplexin the cavity
        for(auto const & simplexInCavity : cavityIO.cellInCavity())
        {
          if(std::find(markedSimplex.begin(), markedSimplex.end(), simplexInCavity) != markedSimplex.end())
          {
            return;
          }
          simplexMesh->deleteTetra(simplexInCavity);
        }

        //delete the surface triangle connected to simpliceNode
        for(auto const & triangleConnectedToP : cavityIO.getTrianglesConnectedToPInCavity())
        {
          simplexMesh->deleteTriangle(triangleConnectedToP);

        }

        //delete the node inside the cavity
        for(auto const & nodeInCavity : cavityIO.getNodeInCavity())
        {
          //simplexMesh->deleteNode(nodeInCavity);
        }

        simplexMesh->rebuildCavity(cavityIO, simpliceNode.getGlobalNode());
        status = true;
      }
    }
}
/******************************************************************************/
