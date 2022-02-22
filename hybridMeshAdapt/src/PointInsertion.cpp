/******************************************************************************/
#include <chrono>
/******************************************************************************/
#include <gmds/hybridMeshAdapt/PointInsertion.h>
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
PointInsertion::PointInsertion(SimplexMesh* simplexMesh, const SimplicesNode& simpliceNode, const CriterionRAIS& criterion, bool& status, const std::vector<TSimplexID>& initialCavity,
                               const gmds::BitVector& markedNodes, std::vector<TInt>& deletedNodes, const std::multimap<TInt, TInt>& facesAlreadyBuilt, std::vector<TSimplexID> markedSimplex)
{
    if(simplexMesh != nullptr)
    {
      Variable<int>* BND_SURFACE_COLOR = simplexMesh->getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
      Variable<int>* BND_TRIANGLES     = simplexMesh->getVariable<int,SimplicesTriangle>("BND_TRIANGLES");

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
        std::copy_if(/*initialCavity*/initCavity.begin(), /*initialCavity*/initCavity.end(), std::back_inserter(initialCavityCell), [&](const TSimplexID simplex){
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

        //std::cout << "CavityEnlargment START " <<std::endl;
        if(cavOp.cavityEnlargement(cavityIO, initialCavityCell, initialCavityTriangle, simpliceNode, criterion, facesAlreadyBuilt, markedSimplex))
        {

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
          //double duration2;
          //start = std::clock();
          //std::cout << "nodeIncavity START " <<std::endl;
          if(!cavityIO.nodeInCavity(simpliceNode.getGlobalNode())){
            status = false;
            return;
          }
          const std::vector<TInt>& nodesInsideCavity = cavityIO.getNodeInCavity();
          //duration2 = (std::clock()-start)/(double)CLOCKS_PER_SEC;
          //std::cout << "nodeInCavity duration --> " << duration2 << std::endl;

          //start = std::clock();
          //std::cout << "nodeReconnection START " <<std::endl;
          cavityIO.nodesReconnection();
          //duration2 = (std::clock()-start)/(double)CLOCKS_PER_SEC;
          //std::cout << "nodesReconnection duration --> " << duration2 << std::endl;
          ////////////////////////////ADRIEN IDEA///////////////////////////////////////
          //this section is here in order to optimize the futurs normals of the created
          //surface for the edgeRemove algorithm
          /*if((*BND_SURFACE_COLOR)[simpliceNode.getGlobalNode()] != 0)
          {
            std::vector<TSimplexID> ball = simpliceNode.ballOf();
            if(ball.size() != 0)
            {
              std::vector<std::vector<TInt>> facesPatchIdx{};
              std::vector<TSimplexID> surfaceCell{};
              std::vector<std::vector<TInt>> surfaceIdx{};

              facesIdxPatch(simplexMesh, cavityIO, facesPatchIdx, simpliceNode);
              std::copy_if(ball.begin(), ball.end(), std::back_inserter(surfaceCell), [&](TSimplexID cellIdx){
                return (cellIdx < 0);
              });

              for(auto const idx : surfaceCell)
              {
                const SimplicesTriangle triange = SimplicesTriangle(simplexMesh, idx);
                std::vector<TInt> nodes         = triange.getNodes();
                surfaceIdx.push_back(nodes);
              }

              double eps0 = simplexMesh->subSurfaceFactor(surfaceIdx);
              double eps1 = simplexMesh->subSurfaceFactor(facesPatchIdx);

              std::cout << "eps0 -> " << eps0 << std::endl;
              std::cout << "eps1 -> " << eps1 << std::endl;
              if(eps1 < eps0 * 0.99)
              {
                return;
              }

            }
          }*/
          ////////////////////////////////////////////////////////////////////////////////
          ////////////////////////////////////////////////////////////////////////////////
          ////////////////////////////////////////////////////////////////////////////////

          ///////////////////////////FRANK IDEA///////////////////////////////////////
          //this section is here in order to optimize the futurs normals of the created
          //surface for the edgeRemove algorithm
          if((*BND_SURFACE_COLOR)[simpliceNode.getGlobalNode()] != 0)
          {
            std::vector<TSimplexID> ball = simpliceNode.ballOf();
            if(ball.size() != 0)
            {
              math::Vector3d interpolationNormal{};
              std::vector<math::Vector3d> normals{};
              normalsPatch(simplexMesh, cavityIO, normals, simpliceNode);

              //on teste la forme de la futur surface qui va être reconstruire et on la compare a l'ancienne
              //pour éviter une mauvaise reconstruction surfacique...
              if(normals.size() != 0)
              {
                unsigned int triangleNbr = 0;
                for(auto const simplex : ball)
                {
                  if(simplex < 0 && simplex != border)
                  {
                    triangleNbr++;
                    math::Vector3d n = SimplicesTriangle(simplexMesh, -simplex).getNormal();
                    n.normalize();
                    interpolationNormal = n + interpolationNormal;
                    interpolationNormal.normalize();
                  }
                }
              }

              double minDot = 0.5;
              for(auto const n : normals)
              {
                if(n.dot(interpolationNormal) < minDot)
                {
                  return;
                }
              }
            }
          }
          ////////////////////////////////////////////////////////////////////////////////
          ////////////////////////////////////////////////////////////////////////////////
          ////////////////////////////////////////////////////////////////////////////////
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
          ////deletedNode insertion
          for(auto const & nodeInCavity : cavityIO.getNodeInCavity())
          {
            if(nodeInCavity != simpliceNode.getGlobalNode())
            {
              deletedNodes.push_back(nodeInCavity);
            }
          }
          for(auto const & surfaceNodeInCavity : cavityIO.getSurfaceNodeInCavity())
          {
            if(surfaceNodeInCavity != simpliceNode.getGlobalNode())
            {
              deletedNodes.push_back(surfaceNodeInCavity);
            }
          }
          //

          //deleteThe simplexin the cavity
          for(auto const & simplexInCavity : cavityIO.cellInCavity())
          {
            if(std::find(markedSimplex.begin(), markedSimplex.end(), simplexInCavity) != markedSimplex.end())
            {
              return;
            }
            simplexMesh->deleteTetra(simplexInCavity);
            //deletedSimplex.push_back(simplexInCavity);
          }

          //delete the surface triangle connected to simpliceNode
          for(auto const & triangleConnectedToP : cavityIO.getTrianglesConnectedToPInCavity())
          {
            //std::cout << "triangleConnectedToP to delete --> " << SimplicesTriangle(simplexMesh, -triangleConnectedToP) << std::endl;
            simplexMesh->deleteTriangle(triangleConnectedToP);
          }
          //double duration3;
          //start = std::clock();
          //std::cout << "rebuildCavity START " << std::endl;
          simplexMesh->rebuildCavity(cavityIO, simpliceNode.getGlobalNode());
          //duration2 = (std::clock()-start)/(double)CLOCKS_PER_SEC;
          //std::cout << "rebuildCavity duration --> " << duration2 << std::endl;
          status = true;
        }
      }
    }
}
/******************************************************************************/
void PointInsertion::facesIdxPatch(SimplexMesh* simplexMesh, CavityOperator::CavityIO& cavityIO, std::vector<std::vector<TInt>>& facesIdx, const SimplicesNode & node) const
{
  facesIdx.clear();
  const std::vector<std::vector<TInt>>& borderEdges = cavityIO.getBorderEdges();
  const std::vector<TInt>&     triangleConnectedToP = cavityIO.getTrianglesConnectedToPInCavity();
  if(triangleConnectedToP.size() != 0)
  {
    for(auto const edge : borderEdges)
    {
      TInt nodeA = edge.front();
      TInt nodeB = edge.back();

      std::vector<TInt> faceIdx{node.getGlobalNode(), nodeA, nodeB};
      facesIdx.push_back(faceIdx);
    }
  }
}
/******************************************************************************/
void PointInsertion::normalsPatch(SimplexMesh* simplexMesh, CavityOperator::CavityIO& cavityIO, std::vector<math::Vector3d>& normals, const SimplicesNode & node) const
{
  normals.clear();
  const std::vector<std::vector<TInt>>& borderEdges = cavityIO.getBorderEdges();
  const std::vector<TInt>&     triangleConnectedToP = cavityIO.getTrianglesConnectedToPInCavity();
  if(triangleConnectedToP.size() != 0)
  {
    for(auto const edge : borderEdges)
    {
      TInt nodeA = edge.front();
      TInt nodeB = edge.back();

      math::Point nodeACoord = SimplicesNode(simplexMesh, nodeA).getCoords();
      math::Point nodeBCoord = SimplicesNode(simplexMesh, nodeB).getCoords();
      math::Point nodeCoord  = node.getCoords();

      math::Vector3d vec0 = nodeACoord - nodeCoord;
      math::Vector3d vec1 = nodeBCoord - nodeACoord;

      math::Vector3d n = vec1.cross(vec0);
      n.normalize();
      normals.push_back(n);
    }
  }
}
/******************************************************************************/
////////////////////////////////////////////////////////////////////////////////
//this section is here in order to optimize the futurs normals of the created
//surface for the edgeRemove algorithm
/*if((*BND_SURFACE_COLOR)[simpliceNode.getGlobalNode()] != 0)
{
  std::vector<TSimplexID> ball = simpliceNode.ballOf();
  if(ball.size() != 0)
  {
    math::Vector3d interpolationNormal{};
    std::vector<math::Vector3d> normals{};
    normalsPatch(simplexMesh, cavityIO, normals, simpliceNode);

    //on teste la forme de la futur surface qui va être reconstruire et on la compare a l'ancienne
    //pour éviter une mauvaise reconstruction surfacique...
    if(normals.size() != 0)
    {
      unsigned int triangleNbr = 0;
      for(auto const simplex : ball)
      {
        if(simplex < 0 && simplex != border)
        {
          triangleNbr++;
          math::Vector3d n = SimplicesTriangle(simplexMesh, -simplex).getNormal();
          n.normalize();
          interpolationNormal += n;
        }
      }
      interpolationNormal.normalize();
      if(triangleNbr != 0)
      {
        interpolationNormal = interpolationNormal / triangleNbr;
      }
      else
      {
        std::cout << "triangleNbr == 0" << std::endl;
        return;
      }
    }

    double minDot = 0.15;
    for(auto const n : normals)
    {
      if(n.dot(interpolationNormal) < minDot)
      {
        std::cout << "n.dot(interpolationNormal) < minDot --> " << n.dot(interpolationNormal) << std::endl;
        return;
      }
    }
  }
}*/
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
