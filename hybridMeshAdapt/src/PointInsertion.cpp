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
PointInsertion::PointInsertion(SimplexMesh* simplexMesh, const SimplicesNode& simpliceNode, const CriterionRAIS& criterion, bool& status, const std::vector<TSimplexID>& initialCavity, const gmds::BitVector& markedNodes, std::vector<TSimplexID>& deletedSimplex, std::vector<TSimplexID> markedSimplex)
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
        cavOp.cavityEnlargement(cavityIO, initialCavityCell, initialCavityTriangle, simpliceNode, criterion, markedSimplex);
        //cavOp.cavityReduction(cavityIO, initCavity, simpliceNode, criterion, cavReduction, markedSimplex);

        //test sur les triangles non connecté a P pour ne pas créer de retournement topologique
        for(auto const triNotCo : cavityIO.getTrianglesNotConnectedToPInCavity())
        {
          const SimplicesTriangle triangle = SimplicesTriangle(simplexMesh, triNotCo);
          math::Orientation::Sign orientation =  triangle.orientation(simpliceNode.getCoords());
          if(orientation < 0)
          {
            status = false;
            return;
          }
        }
        ////////////////////////////////////////////////////////////////////////////////
        ///////////////////////finding the node inside the cavity///////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        const std::vector<TInt>& nodesInsideCavity = cavityIO.nodeInCavity();
        cavityIO.nodesReconnection();

        //check if node in nodsIncavity are marked
        for(auto const & nodeInCavity : cavityIO.getNodeInCavity())
        {
          if(markedNodes[nodeInCavity] == 1 && nodeInCavity != simpliceNode.getGlobalNode())
          {
            status = false;
            return;
          }
        }

        for(auto const & surfaceNodeInCavity : cavityIO.getSurfaceNodeInCavity())
        {
          if(markedNodes[surfaceNodeInCavity] == 1 && surfaceNodeInCavity != simpliceNode.getGlobalNode())
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
          deletedSimplex.push_back(simplexInCavity);
        }

        //delete the surface triangle connected to simpliceNode
        for(auto const & triangleConnectedToP : cavityIO.getTrianglesConnectedToPInCavity())
        {
          simplexMesh->deleteTriangle(triangleConnectedToP);
          //deletedSimplex.push_back(triangleConnectedToP);
        }

        simplexMesh->rebuildCavity(cavityIO, simpliceNode.getGlobalNode());
        status = true;
      }
    }
}
/******************************************************************************/
