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
                               const gmds::BitVector& markedNodes, std::vector<TInt>& deletedNodes, const std::multimap<TInt, TInt>& facesAlreadyBuilt, std::vector<TSimplexID>& createdCells, std::vector<TSimplexID> markedSimplex)
{
    if(simplexMesh != nullptr)
    {
      auto my_make = [=](const TInt a, const TInt b){
        return (a < b)? std::make_pair(a, b) : std::make_pair(b, a);
      };

        createdCells.clear();
        Variable<int>* BND_VERTEX_COLOR   = nullptr;
        Variable<int>* BND_CURVE_COLOR   = nullptr;
        Variable<int>* BND_SURFACE_COLOR = nullptr;
        Variable<int>* BND_TRIANGLES     = nullptr;
      try{
        BND_CURVE_COLOR   = simplexMesh->getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
        BND_SURFACE_COLOR = simplexMesh->getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
        BND_VERTEX_COLOR  = simplexMesh->getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
        BND_TRIANGLES     = simplexMesh->getVariable<int,SimplicesTriangle>("BND_TRIANGLES");

      } catch(gmds::GMDSException e)
      {
        throw gmds::GMDSException(e);
      }

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
        const std::vector<TSimplexID>& base = simplexMesh->getBase();

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
          if(!cavityIO.nodeInCavity(simpliceNode.getGlobalNode())){
            status = false;
            return;
          }

          cavityIO.nodesReconnection(simpliceNode.getGlobalNode());

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
                  status = false;
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
          std::vector<TInt> curveNodeToDel{};
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
              curveNodeToDel.push_back(surfaceNodeInCavity);
            }
          }

          //deleteThe simplexin the cavity
          const gmds::BitVector& bitVectorTet = simplexMesh->getBitVectorTet();
          std::vector<TSimplexID> ball = simpliceNode.ballOf();
          std::vector<std::vector<TInt>> deleted_Tet{};

          for(auto const & simplexInCavity : cavityIO.cellInCavity())
          {
            if(std::find(markedSimplex.begin(), markedSimplex.end(), simplexInCavity) != markedSimplex.end())
            {
              return;
            }
            const std::vector<TInt> nodes = SimplicesCell(simplexMesh, simplexInCavity).getNodes();
            deleted_Tet.push_back(nodes);
            simplexMesh->deleteTetra(simplexInCavity);
            //deletedSimplex.push_back(simplexInCavity);
          }



          //delete the surface triangle connected to simpliceNode
          //use delete_tri in order to rebuild the edge if it was destroy
          //If the node is on ridge update edgeStructure
          unsigned int labelCurve = (*BND_CURVE_COLOR)[simpliceNode.getGlobalNode()];
          unsigned int labelCorner = (*BND_VERTEX_COLOR)[simpliceNode.getGlobalNode()];
          const std::map<unsigned int, std::vector<unsigned int>>& cornerEdgeConnexion = simplexMesh->getCornerEdgeConnexion();
          unsigned int sizeFace = 3;
          if(labelCorner != 0)
          {
            std::set<std::pair<TInt, TInt>> mapNodesOnRidge{};
            std::multimap<TInt, std::pair<TInt,TInt>>& edgeStructure = simplexMesh->getEdgeStructure();
            std::vector<unsigned int> edgesIndices = cornerEdgeConnexion.at(labelCorner);
            //check the edge that will be destroy
            if(edgesIndices.size() == 0)
            {
              throw gmds::GMDSException("edgesIndices.size() == 0");
            }
            for(auto const edgeIdx : edgesIndices)
            {
              for(auto const tri : cavityIO.getTrianglesConnectedToPInCavity())
              {
                const std::vector<TInt> nodes = SimplicesTriangle(simplexMesh, tri).getNodes();
                std::vector<std::vector<TInt>> nodesOnRidge{};

                for(unsigned int lN = 0 ; lN < sizeFace ; lN++)
                {
                  TInt n0 = nodes[lN];
                  TInt n1 = nodes[(lN + 1) % sizeFace];
                  if(labelCorner == (*BND_VERTEX_COLOR)[n0] || (*BND_CURVE_COLOR)[n0] == edgeIdx)
                  {
                    if(labelCorner == (*BND_VERTEX_COLOR)[n1] || (*BND_CURVE_COLOR)[n1] == edgeIdx)
                    {
                      std::pair<TInt, TInt> p = my_make(n0, n1);
                      mapNodesOnRidge.insert(p);
                    }
                  }
                }
              }

              if(mapNodesOnRidge.size() != 0)
              {
                auto it = edgeStructure.equal_range(edgeIdx);
                std::vector<TInt> nodesEdge{};
                std::vector<TInt> borderNodesEdge{};
                for(auto const edge : mapNodesOnRidge)
                {
                  TInt nodeA = edge.first;
                  TInt nodeB = edge.second;
                  for(auto itr = it.first ; itr != it.second; ++itr)
                  {
                    //destroy the node on the edge STRUCTURE
                    if(nodeA == itr->second.first && nodeB == itr->second.second)
                    {
                      //delete an element of the map during the iteration over this map
                      itr = edgeStructure.erase(itr);
                      nodesEdge.push_back(nodeA);
                      nodesEdge.push_back(nodeB);
                      it = edgeStructure.equal_range(edgeIdx);
                      break;
                    }
                  }
                }

                std::sort(nodesEdge.begin(), nodesEdge.end());
                for(unsigned int i = 0 ; i < nodesEdge.size();)
                {
                  TInt node_i = nodesEdge[i];
                  TInt node_j = nodesEdge[i+1];
                  if(node_i != node_j)
                  {
                    borderNodesEdge.push_back(node_i);
                    if(i == nodesEdge.size() -2)
                    {
                      borderNodesEdge.push_back(node_j);
                      break;
                    }
                    ++i;
                  }
                  else
                  {
                    if(i == nodesEdge.size() - 2)
                    {
                        break;
                    }
                    else
                    {
                        i = i+2;
                    }
                  }
                }
                //rebuild the edge structure with the node to connect and the border edge STRUCTURE
                if(borderNodesEdge.size() == 2)
                {
                  if(borderNodesEdge.front() == simpliceNode.getGlobalNode() || borderNodesEdge.back() == simpliceNode.getGlobalNode())
                  {
                    std::pair<TInt, TInt> p = my_make(borderNodesEdge.front(), borderNodesEdge.back());
                    edgeStructure.insert(std::pair<TInt, std::pair<TInt, TInt>>(edgeIdx, p));
                  }
                  else
                  {
                    std::pair<TInt, TInt> p0 = my_make(simpliceNode.getGlobalNode(), borderNodesEdge.back());
                    std::pair<TInt, TInt> p1 = my_make(simpliceNode.getGlobalNode(), borderNodesEdge.front());
                    edgeStructure.insert(std::pair<TInt, std::pair<TInt, TInt>>(edgeIdx, p0));
                    edgeStructure.insert(std::pair<TInt, std::pair<TInt, TInt>>(edgeIdx, p1));
                  }
                }
              }
            }
          }
          else if(labelCurve != 0)
          {
            std::set<std::pair<TInt, TInt>> mapNodesOnRidge{};
            std::multimap<TInt, std::pair<TInt,TInt>>& edgeStructure = simplexMesh->getEdgeStructure();
            /*for(auto const p : edgeStructure)
            {
              std::cout << "edge idx -> " << p.first << " | nodes -> " << p.second.first << " : " << p.second.second << std::endl;
            }*/
            //check the edge that will be destroy
            for(auto const tri : cavityIO.getTrianglesConnectedToPInCavity())
            {
              const std::vector<TInt> nodes = SimplicesTriangle(simplexMesh, tri).getNodes();
              std::vector<TInt> nodesOnRidge{};
              for(unsigned int lN = 0 ; lN < sizeFace ; lN++)
              {
                TInt n0 = nodes[lN];
                TInt n1 = nodes[(lN + 1) % sizeFace];
                std::vector<unsigned int> edgesIndices;

                if((*BND_VERTEX_COLOR)[n0] != 0 && (*BND_VERTEX_COLOR)[n1] != 0)
                {
                  std::vector<unsigned int> edgesIndices0 = cornerEdgeConnexion.at((*BND_VERTEX_COLOR)[n0]);
                  std::vector<unsigned int> edgesIndices1 = cornerEdgeConnexion.at((*BND_VERTEX_COLOR)[n1]);
                  if(std::find(edgesIndices0.begin(), edgesIndices0.end(), labelCurve) == edgesIndices0.end() ||
                      std::find(edgesIndices1.begin(), edgesIndices1.end(), labelCurve) == edgesIndices1.end())
                  {
                    continue;
                  }
                }
                else if((*BND_VERTEX_COLOR)[n0] != 0)
                {
                  edgesIndices = cornerEdgeConnexion.at((*BND_VERTEX_COLOR)[n0]);
                }
                else if((*BND_VERTEX_COLOR)[n1] != 0)
                {
                  edgesIndices = cornerEdgeConnexion.at((*BND_VERTEX_COLOR)[n1]);
                }

                if(edgesIndices.size() == 0)
                {
                  if(labelCurve == (*BND_CURVE_COLOR)[n0] && labelCurve == (*BND_CURVE_COLOR)[n1])
                  {
                    std::pair<TInt, TInt> p = my_make(n0, n1);
                    mapNodesOnRidge.insert(p);
                  }
                }
                else
                {
                  if(labelCurve == (*BND_CURVE_COLOR)[n0] || std::find(edgesIndices.begin(), edgesIndices.end(), labelCurve) != edgesIndices.end())
                  {
                    if(labelCurve == (*BND_CURVE_COLOR)[n1] || std::find(edgesIndices.begin(), edgesIndices.end(), labelCurve) != edgesIndices.end())
                    {
                      std::pair<TInt, TInt> p = my_make(n0, n1);
                      mapNodesOnRidge.insert(p);
                    }
                  }
                }
              }
            }
            for(auto const edge : mapNodesOnRidge)
            {
              TInt nodeA = edge.first;
              TInt nodeB = edge.second;
            }
            if(mapNodesOnRidge.size() != 0)
            {
              auto it = edgeStructure.equal_range(labelCurve);
              std::vector<TInt> nodesEdge{};
              std::vector<TInt> borderNodesEdge{};
              for(auto const edge : mapNodesOnRidge)
              {
                TInt nodeA = edge.first;
                TInt nodeB = edge.second;

                for(auto itr = it.first ; itr != it.second ; ++itr)
                {
                  //destroy the edge on the edge STRUCTURE
                  if(nodeA == itr->second.first && nodeB == itr->second.second)
                  {
                    //delete an element of the map during the iteration over this map
                    edgeStructure.erase(itr);
                    nodesEdge.push_back(nodeA);
                    nodesEdge.push_back(nodeB);
                    it = edgeStructure.equal_range(labelCurve);
                    break;
                  }
                }
              }
              std::sort(nodesEdge.begin(), nodesEdge.end());
              for(unsigned int i = 0 ; i < nodesEdge.size();)
              {
                TInt node_i = nodesEdge[i];
                TInt node_j = nodesEdge[i+1];
                if(node_i != node_j)
                {
                  borderNodesEdge.push_back(node_i);
                  if(i == nodesEdge.size() -2)
                  {
                    borderNodesEdge.push_back(node_j);
                    break;
                  }
                  ++i;
                }
                else
                {
                  if(i == nodesEdge.size() - 2)
                  {
                      break;
                  }
                  else
                  {
                      i = i+2;
                  }
                }
              }
              //for(auto const node : borderNodesEdge){std::cout << " borderNodesEdge **** " << node << std::endl;}
              //rebuild the edge structure with the node to connect and the border edge STRUCTURE
              if(borderNodesEdge.size() != 2)
              {
                gmds::ISimplexMeshIOService ioServiceMesh(simplexMesh);
                gmds::VTKWriter vtkWriterCS(&ioServiceMesh);
                vtkWriterCS.setCellOptions(gmds::N|gmds::R|gmds::F);
                vtkWriterCS.setDataOptions(gmds::N|gmds::R|gmds::F);
                vtkWriterCS.write("no_borderNodesEdge.vtk");
                std::cout << "node being inserted -> " << simpliceNode.getGlobalNode() << std::endl;
                std::cout << "surface -> " << (*BND_SURFACE_COLOR)[simpliceNode.getGlobalNode()] << std::endl;
                std::cout << "curve -> " << (*BND_CURVE_COLOR)[simpliceNode.getGlobalNode()] << std::endl;
                std::cout << "borderNodesEdge.size() -> " << borderNodesEdge.size() << std::endl;
                throw gmds::GMDSException("borderNodesEdge.size() != 2");
              }

              if(borderNodesEdge.front() == simpliceNode.getGlobalNode() || borderNodesEdge.back() == simpliceNode.getGlobalNode())
              {
                std::pair<TInt, TInt> p = my_make(borderNodesEdge.front(), borderNodesEdge.back());
                edgeStructure.insert(std::pair<TInt, std::pair<TInt, TInt>>(labelCurve, p));
              }
              else
              {
                std::pair<TInt, TInt> p0 = my_make(simpliceNode.getGlobalNode(), borderNodesEdge.back());
                std::pair<TInt, TInt> p1 = my_make(simpliceNode.getGlobalNode(), borderNodesEdge.front());
                edgeStructure.insert(std::pair<TInt, std::pair<TInt, TInt>>(labelCurve, p0));
                edgeStructure.insert(std::pair<TInt, std::pair<TInt, TInt>>(labelCurve, p1));
              }
            }
          }
          std::multimap<TInt, std::pair<TInt,TInt>>& edgeStructure = simplexMesh->getEdgeStructure();

          std::vector<std::vector<TInt>> deleted_Tri{};
          for(auto const & triangleConnectedToP : cavityIO.getTrianglesConnectedToPInCavity())
          {
            const std::vector<TInt> nodes = SimplicesTriangle(simplexMesh, triangleConnectedToP).getNodes();
            deleted_Tri.push_back(nodes);
            simplexMesh->deleteTriangle(triangleConnectedToP);
          }

          //In case if the node don't have any information of the base[node] for the reinsertion
          //check pair.second.size() == 1 in rebuildCavity
          if(base[simpliceNode.getGlobalNode()] == border)
          {
            std::vector<TSimplexID> ball = simpliceNode.ballOf();
            for(auto const simplex : ball)
            {
              if(simplex >= 0)
              {
                if(bitVectorTet[simplex] != 0)
                {
                  simplexMesh->setBase(simpliceNode.getGlobalNode(), simplex);
                }
              }
            }
          }
          simplexMesh->rebuildCavity(cavityIO, deleted_Tet, deleted_Tri, simpliceNode.getGlobalNode(), createdCells);

          //simplexMesh->rebuildCav(cavityIO, deleted_Tet, deleted_Tri, simpliceNode.getGlobalNode(), createdCells);
          simplexMesh->checkMeshCavity(createdCells);
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

      //if this condition happen that's mean an reinsertion is occuring
      if(nodeA == node.getGlobalNode() || nodeB == node.getGlobalNode())
      {
        continue;
      }
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
