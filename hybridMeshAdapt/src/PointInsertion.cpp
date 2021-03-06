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
        Variable<int>* BND_CURVE_COLOR   = nullptr;
        Variable<int>* BND_SURFACE_COLOR = nullptr;
        Variable<int>* BND_TRIANGLES     = nullptr;
      try{
        BND_CURVE_COLOR   = simplexMesh->getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
        BND_SURFACE_COLOR = simplexMesh->getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
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

          //test sur les triangles non connect?? a P pour ne pas cr??er de retournement topologique
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
          const std::vector<TInt>& nodesInsideCavity = cavityIO.getNodeInCavity();

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

              //on teste la forme de la futur surface qui va ??tre reconstruire et on la compare a l'ancienne
              //pour ??viter une mauvaise reconstruction surfacique...
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


          /*if(simpliceNode.getGlobalNode() == 271729){
            //simplexMesh->deleteAllSimplicesBut(cavityIO.cellInCavity());
            for(auto const & triangleConnectedToP : cavityIO.getTrianglesConnectedToPInCavity())
            {
              std::cout << "triangleConnectedToP to delete --> " << SimplicesTriangle(simplexMesh, -triangleConnectedToP) << std::endl;
            }
            for(auto const & triangleNotConnectedToP : cavityIO.getTrianglesNotConnectedToPInCavity())
            {
              std::cout << "triangleNotConnectedToP --> " << SimplicesTriangle(simplexMesh, -triangleNotConnectedToP) << std::endl;
            }
            gmds::ISimplexMeshIOService ioService0(simplexMesh);
            gmds::VTKWriter vtkWriter0(&ioService0);
            vtkWriter0.setCellOptions(gmds::N|gmds::R|gmds::F);
            vtkWriter0.setDataOptions(gmds::N|gmds::R|gmds::F);
            vtkWriter0.write("MESH_BUG_" + std::to_string(simpliceNode.getGlobalNode()) + ".vtk");
            throw gmds::GMDSException("OKOK");
          }*/
          //SimplexMesh beforeBugMesh = *simplexMesh;
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


          //If the node is on ridge update edgeStructure
          unsigned int label = (*BND_CURVE_COLOR)[simpliceNode.getGlobalNode()];
          if(label != 0)
          {
            std::multimap<TInt, std::pair<TInt,TInt>>& edgeStructure = simplexMesh->getEdgeStructure();
            std::multimap<TInt, std::pair<TInt,TInt>> edgeStructureCopy = edgeStructure;

            auto it = edgeStructure.equal_range(label);
            if(curveNodeToDel.size() == 0)
            {
              const std::pair<TInt, TInt> edge = cavityIO.getEdgeContainingNode();
              //deletion of the edge in edgestructure
              std::pair<TInt, TInt> newEdge0{};
              std::pair<TInt, TInt> newEdge1{};
              for(auto itr = it.first ; itr != it.second ; itr++)
              {
                if(itr->second.first == edge.first && itr->second.second == edge.second)
                {
                  edgeStructure.erase(itr);
                  newEdge0 = std::make_pair(std::min(simpliceNode.getGlobalNode(), edge.first), std::max(simpliceNode.getGlobalNode(), edge.first));
                  newEdge1 = std::make_pair(std::min(simpliceNode.getGlobalNode(), edge.second), std::max(simpliceNode.getGlobalNode(), edge.second));
                  if(newEdge0.first != newEdge0.second){edgeStructure.insert(std::make_pair(label, newEdge0));}
                  if(newEdge1.first != newEdge1.second){edgeStructure.insert(std::make_pair(label, newEdge1));}
                  break;
                }
              }
              for(auto const data : edgeStructure)
              {
                if(data.second.first == data.second.second)
                {
                  //if DEBUG
                  std::cout << "EDGE STRUCTURE BEFORE MODIFICATION FAILED" << std::endl;
                  for(auto const dataBis : edgeStructureCopy)
                  {
                    std::cout << "data --> " << dataBis.first << " | [" << dataBis.second.first << " : " << dataBis.second.second << "]" << std::endl;
                  }
                  std::cout << std::endl;
                  std::cout << "edge being deleted -> " << edge.first << " | " << edge.second << std::endl;
                  std::cout << "NEW EDGE CREATED -> " << newEdge0.first << " | " << newEdge0.second << std::endl;
                  std::cout << "NEW EDGE CREATED -> " << newEdge1.first << " | " << newEdge1.second << std::endl;
                  std::cout << "DATA TROUBLE | label --> " << data.first << " | [" << data.second.first << " : " << data.second.second << "]" << std::endl;

                    gmds::ISimplexMeshIOService ioService(simplexMesh);
                    gmds::VTKWriter vtkWriter(&ioService);
                    vtkWriter.setCellOptions(gmds::N|gmds::R|gmds::F);
                    vtkWriter.setDataOptions(gmds::N|gmds::R|gmds::F);
                    vtkWriter.write("EDGE_STRUCTURE_MODIFICATION_BUG_" + std::to_string(simpliceNode.getGlobalNode()) + ".vtk");

                  throw gmds::GMDSException("data.second.first == data.second.second");
                }
              }
            }
            else
            {
              std::vector<std::multimap<TInt, std::pair<TInt,TInt>>::iterator> iteratorsToErase{};
              std::vector<std::multimap<TInt, std::pair<TInt,TInt>>::iterator> iteratorsToReconnect{};

              const gmds::BitVector & nodes_ids = simplexMesh->getBitVectorNodes();
              for(auto const node : curveNodeToDel)
              {
                if(nodes_ids[node] != 0)
                {
                  if(label == (*BND_CURVE_COLOR)[node])
                  {
                    for(auto itr = it.first ; itr != it.second ; itr++)
                    {
                      if(itr->second.first == node)
                      {
                        itr->second.first = border;
                      }
                      if(itr->second.second == node)
                      {
                        itr->second.second = border;
                      }

                      TInt minNode = std::min(itr->second.first, itr->second.second);
                      TInt maxNode = std::max(itr->second.first, itr->second.second);

                      itr->second.first  = minNode;
                      itr->second.second = maxNode;

                      if(itr->second.first == border && itr->second.second == border)
                      {
                        iteratorsToErase.push_back(itr);
                      }
                      else if(itr->second.first == border)
                      {
                        iteratorsToReconnect.push_back(itr);
                      }
                    }
                  }
                }
              }

              if(iteratorsToReconnect.size() == 1 || iteratorsToReconnect.size() > 2)
              {
                //Problem de structure invalide mais impossible aprevoir (voir node 286 mesh S23.vtk) 
                //if DEBUG
                /*std::cout << "EDGE STRUCTURE BEFORE MODIFICATION FAILED" << std::endl;
                for(auto const dataBis : edgeStructureCopy)
                {
                  std::cout << "data --> " << dataBis.first << " | [" << dataBis.second.first << " : " << dataBis.second.second << "]" << std::endl;
                }
                for(auto const itr : iteratorsToReconnect)
                {
                  std::cout << "iteratorsToReconnect --> " << itr->second.first << " | " << itr->second.second << std::endl;
                }
                std::cout << std::endl;

                SimplexMesh nodeMesh = SimplexMesh();
                nodeMesh.addNode(simpliceNode.getCoords()[0], simpliceNode.getCoords()[1], simpliceNode.getCoords()[2]);
                nodeMesh.addTetraedre(0, 0, 0, 0);
                gmds::ISimplexMeshIOService ioService(&nodeMesh);
                gmds::VTKWriter vtkWriter(&ioService);
                vtkWriter.setCellOptions(gmds::N|gmds::R|gmds::F);
                vtkWriter.setDataOptions(gmds::N|gmds::R|gmds::F);
                vtkWriter.write("NODE_" + std::to_string(simpliceNode.getGlobalNode()) + ".vtk");

                std::cout << "iteratorsToReconnect.size() -- > " << iteratorsToReconnect.size() << std::endl;
                throw gmds::GMDSException("iteratorsToReconnect.size() PROBLEM");*/
              }
              else
              {
                for(auto itr : iteratorsToReconnect)
                {
                    TInt node = itr->second.second;
                    if(node == simpliceNode.getGlobalNode()){iteratorsToErase.push_back(itr); continue;}
                    TInt minNode = std::min(simpliceNode.getGlobalNode(), itr->second.second);
                    TInt maxNode = std::max(simpliceNode.getGlobalNode(), itr->second.second);
                    itr->second.first  = minNode;
                    itr->second.second = maxNode;
                    if(minNode == maxNode){iteratorsToErase.push_back(itr);}
                }
              }
              for(auto itr : iteratorsToErase){edgeStructure.erase(itr);}


              for(auto const data : edgeStructure)
              {
                if(data.second.first == data.second.second)
                {
                  std::cout << "data --> " << data.first << " | [" << data.second.first << " : " << data.second.second << "]" << std::endl;
                  throw gmds::GMDSException("data.second.first == data.second.second");
                }
              }
            }
          }

          simplexMesh->rebuildCavity(cavityIO, simpliceNode.getGlobalNode());
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

    //on teste la forme de la futur surface qui va ??tre reconstruire et on la compare a l'ancienne
    //pour ??viter une mauvaise reconstruction surfacique...
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
