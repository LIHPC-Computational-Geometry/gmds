/******************************************************************/
#include <gmds/hybridMeshAdapt/CavityOperator.h>
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
#include <gmds/hybridMeshAdapt/ICriterion.h>
#include "gmds/hybridMeshAdapt/PointInsertion.h"
/******************************************************************/
#include <chrono>
/******************************************************************/
using namespace gmds;
using namespace hybrid;
using namespace operators;
using namespace simplicesCell;
using namespace simplicesNode;
using namespace simplicesTriangle;
/******************************************************************/
CavityOperator::CavityOperator(SimplexMesh* simplexMesh):
m_simplex_mesh(simplexMesh)
{
  m_tetSelectedIds = BitVector(m_simplex_mesh->tetCapacity());
}
/******************************************************************/
CavityOperator::CavityIO::CavityIO(SimplexMesh* simplexMesh, const std::vector<TSimplexID>& cavityIn, const std::vector<TSimplexID>& cavityTriangleConnectedToP, const std::vector<TSimplexID>& cavityTriangleNotConnectedToP):
m_simplex_mesh(simplexMesh)
{
  m_cavityCellIn = cavityIn;
  m_cavityTriangleConnectedToP = cavityTriangleConnectedToP;
  m_cavityTriangleNotConnectedToP = cavityTriangleNotConnectedToP;
  //cavityIn.clear();
}
/******************************************************************/
CavityOperator::CavityIO::CavityIO(CavityIO&& cavityIO)
{
  m_simplex_mesh = cavityIO.m_simplex_mesh;
  m_cavityCellIn  = cavityIO.m_cavityCellIn;
  m_cavityTriangleConnectedToP = cavityIO.m_cavityTriangleConnectedToP;
  m_cavityTriangleNotConnectedToP = cavityIO.m_cavityTriangleNotConnectedToP;

  cavityIO.m_cavityCellIn.clear();
  cavityIO.m_cavityTriangleConnectedToP.clear();
  cavityIO.m_cavityTriangleNotConnectedToP.clear();
  cavityIO.m_simplex_mesh = nullptr;
}
/******************************************************************/
CavityOperator::CavityIO& CavityOperator::CavityIO::CavityIO::operator=(CavityIO&& cavityIO)
{
  if(this != &cavityIO)
  {
    m_cavityCellIn.clear();
    /*A VOIR SI JAMAIS IL YA UN PROBLEME LORS D'UNE COPIE*/
    m_simplex_mesh = nullptr;

    m_cavityCellIn  = cavityIO.m_cavityCellIn;
    m_cavityTriangleConnectedToP = cavityIO.m_cavityTriangleConnectedToP;
    m_cavityTriangleNotConnectedToP = cavityIO.m_cavityTriangleNotConnectedToP;
    m_simplex_mesh = cavityIO.m_simplex_mesh;

    cavityIO.m_cavityCellIn.clear();
    cavityIO.m_cavityTriangleConnectedToP.clear();
    cavityIO.m_cavityTriangleNotConnectedToP.clear();
    cavityIO.m_simplex_mesh = nullptr;
  }
  return *this;
}
/******************************************************************/
bool CavityOperator::CavityIO::CavityIO::nodeInCavity(const TInt node)
{
  m_nodeInCavity.clear();
  m_nodesToReconnect.clear();
  m_oppositeCell.clear();

  std::vector<TInt> nodes;
  TInt border = std::numeric_limits<int>::min();
  unsigned int sizeCell = 4;
  unsigned int sizeTriangle = 3;

  gmds::BitVector cellBitvector(m_simplex_mesh->tetCapacity());
  gmds::BitVector cellBitvectorBis(m_simplex_mesh->tetCapacity());
  gmds::BitVector nodesBitvector(m_simplex_mesh->nodesCapacity());
  gmds::BitVector borderNodesBitvector(m_simplex_mesh->nodesCapacity());


  std::vector<TInt> allNodes{};
  std::vector<std::vector<unsigned int>> localsBorderNode{};
  std::vector<std::vector<TSimplexID>> oppisitesSimplex{};

  //////////////////////////////////////////////////////////////////////////////
  ////////////////////this section found the interior nodes with cells//////////
  //////////////////////////////////////////////////////////////////////////////
  for(auto const & tet : cellInCavity())
  {
    if(tet != border)
    {
      cellBitvector.assign(tet);
      cellBitvectorBis.assign(tet);
      SimplicesCell cell = SimplicesCell(m_simplex_mesh, tet);
      const std::vector<TInt>&& nodes = cell.getNodes();
      for(auto const node : nodes)
      {
        if(nodesBitvector[node] == 0)
        {
          allNodes.push_back(node);
          nodesBitvector.assign(node);
        }
      }
    }
  }

  cellBitvectorBis.clear();
  for(auto const & tet : cellInCavity())
  {
    cellBitvectorBis.assign(tet);
    SimplicesCell cell = SimplicesCell(m_simplex_mesh, tet);
    for(unsigned int localIndex = 0 ; localIndex < sizeCell ; localIndex++)
    {
      //look if the oppositeCell is in the cavity with the help of tetBitVector
      TSimplexID oppositeCell = cell.oppositeTetraIdx(localIndex);
      const std::vector<TInt>&& nodes =  cell.getOrderedFace(localIndex);

      if(oppositeCell >= 0)
      {
        if(cellBitvector[oppositeCell] == 0)
        {
          m_nodesToReconnect.push_back(nodes);
          m_oppositeCell.push_back(oppositeCell);
        }
      }
      else
      {
        m_nodesToReconnect.push_back(nodes);
        m_oppositeCell.push_back(oppositeCell);
      }
    }
  }

  for(auto const nodes : m_nodesToReconnect)
  {
    for(auto const node : nodes)
    {
      if(borderNodesBitvector[node] == 0)
      {
        borderNodesBitvector.assign(node);
      }
    }
  }

  //now we extract the interior nodes of the cavity with the help of allNodes & borderNodesBitvector
  std::copy_if(allNodes.begin(), allNodes.end(), std::back_inserter(m_nodeInCavity), [&](const TInt node)
  {
    return (borderNodesBitvector[node] == 0);
  });

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  ///////this section found the interior nodes with help of triangleCoToP///////
  //////////////////////////////////////////////////////////////////////////////
  borderNodesBitvector.clear();
  nodesBitvector.clear();
  borderNodesBitvector.resize(m_simplex_mesh->nodesCapacity());
  nodesBitvector.resize(m_simplex_mesh->nodesCapacity());

  std::vector<std::vector<unsigned int>> localsBorderNodeSurface{};

  std::vector<TInt> allNodeInSurface{};
  gmds::BitVector triBitvector(m_simplex_mesh->triCapacity());
  for(auto const & tri : getTrianglesConnectedToPInCavity())
  {
    if(tri != border)
    {
      triBitvector.assign(-tri);
      SimplicesTriangle triangle = SimplicesTriangle(m_simplex_mesh, tri);
      const std::vector<TInt>&& nodes = triangle.getNodes();
      for(auto const node : nodes)
      {
        if(nodesBitvector[node] == 0)
        {
          allNodeInSurface.push_back(node);
          nodesBitvector.assign(node);
        }
      }
    }
  }

  Variable<int>* BND_TRIANGLES      = m_simplex_mesh->getVariable<int,SimplicesTriangle>("BND_TRIANGLES");
  Variable<int>* BND_VERTEX_COLOR   = m_simplex_mesh->getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  Variable<int>* BND_CURVE_COLOR    = m_simplex_mesh->getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR  = m_simplex_mesh->getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  for(auto const & tri : getTrianglesConnectedToPInCavity())
  {
    SimplicesTriangle triangle = SimplicesTriangle(m_simplex_mesh, -tri);

    for(unsigned int localIndex = 0 ; localIndex < sizeTriangle ; localIndex++)
    {
      //look if the oppositeCell is in the cavity with the help of tetBitVector
      TSimplexID oppositeTriangle = triangle.neighborTriangle(localIndex);
      if(oppositeTriangle != border)
      {
        if(triBitvector[oppositeTriangle] == 0)
        {
          //fill m_borderSurfaceNode
          std::vector<TInt> borderEdge = triangle.getOppositeEdge(localIndex);

          //condition to look if an indesirable triangle will be built during the cavity's rebuilding
          if((*BND_CURVE_COLOR)[node] != 0 && ((*BND_CURVE_COLOR)[borderEdge.front()] == (*BND_CURVE_COLOR)[node]) && ((*BND_CURVE_COLOR)[borderEdge.front()] == (*BND_CURVE_COLOR)[borderEdge.back()]))
          {
            return false;
          }
          for(auto const node : borderEdge)
          {
            if(borderNodesBitvector[node] == 0)
            {
              borderNodesBitvector.assign(node);
            }
          }
          m_borderSurfaceNode.push_back(borderEdge);
        }
      }
    }
  }

  //now we extract the interior nodes of the cavity with the help of allNodes & borderNodesBitvector
  std::copy_if(allNodeInSurface.begin(), allNodeInSurface.end(), std::back_inserter(m_surfaceNodeInCavity), [&](const TInt node)
  {
    return (borderNodesBitvector[node] == 0);
  });
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  return true;
}
/******************************************************************************/
void CavityOperator::CavityIO::CavityIO::nodesReconnection()
{
  Variable<int>* BND_TRIANGLES = m_simplex_mesh->getVariable<int,SimplicesTriangle>("BND_TRIANGLES");

  m_localsNodeForReconnectionWithTriangle.clear();
  m_oppositeTriangle.clear();
  m_triangleIndices.clear();

  unsigned int sizeFace = 3;
  int border = std::numeric_limits<int>::min();
  gmds::BitVector triBitvector(m_simplex_mesh->triCapacity());
  for(auto const & tri : getTrianglesConnectedToPInCavity())
  {
    if(tri != border)
    {
      if(triBitvector[-tri] == 0)
      {
        triBitvector.assign(-tri);
      }
    }
  }

  //looking for the node border surface node contain in the m_cavityTriangleConnectedToP vector
  std::vector<TInt> initVec{border, border, border};
  m_localsNodeForReconnectionWithTriangle.resize(m_nodesToReconnect.size(), initVec);
  m_oppositeTriangle.resize(m_nodesToReconnect.size(), initVec);
  m_triangleIndices.resize(m_nodesToReconnect.size(), initVec);

  unsigned int faceCpt = 0;
  for(auto const& nodes : m_nodesToReconnect)
  {
    TInt node0 = nodes[0];
    TInt nodeA = nodes[1];
    TInt nodeB = nodes[2];

    for(int localNode = 0 ; localNode < sizeFace ; localNode++)
    {
      TInt node0 = nodes[localNode];
      TInt nodeA = nodes[(localNode + 1) % sizeFace];
      TInt nodeB = nodes[(localNode + 2) % sizeFace];

      for(auto const borderSurface : m_borderSurfaceNode)
      {
        TInt tmpNodeA = borderSurface.front();
        TInt tmpNodeB = borderSurface.back();

        if(tmpNodeA == nodeB && tmpNodeB == nodeA)
        {
          m_localsNodeForReconnectionWithTriangle[faceCpt][localNode] = localNode;
          const std::vector<TSimplexID>&& shell = SimplicesNode(m_simplex_mesh, nodeA).shell(SimplicesNode(m_simplex_mesh, nodeB));
          for(auto const simplex : shell)
          {
            if(simplex < 0)
            {
              if(triBitvector[-simplex] == 0)
              {
                m_oppositeTriangle[faceCpt][localNode] = simplex;
              }
              if(triBitvector[-simplex] == 1)
              {
                m_triangleIndices[faceCpt][localNode] = (*BND_TRIANGLES)[-simplex];
              }
            }
          }
        }
      }
    }
    faceCpt++;
  }
}
/******************************************************************************/
std::vector<std::vector<TInt>> CavityOperator::CavityIO::CavityIO::pointToConnect(std::vector<TSimplexID>& orderedAdjTet, std::set<TSimplexID>& complementarySimplex,const std::vector<TInt>& nodesInsideCavity, std::vector<TInt>& cellOfPointToConnect, std::vector<TInt>& nodeNotConnected, const TInt nodeToInsert) const
{
  unsigned int sizeNodeInSimplex = 4;
  std::vector<std::vector<TInt>> v{};
  std::vector<std::vector<TInt>> faceToConnect{};

  TInt errorId = std::numeric_limits<int>::min();
  std::set<TInt> allNodes{};
  std::set<TInt> nodesToConnect{};
  std::set<TInt> nodesLimitToRemovefromBuildSimplexTMP{};
  std::set<TInt> nodesLimitToRemovefromBuildSimplex{};

  for(const auto & simplexInCavity : cellInCavity())
  {
    std::vector<TInt> pointToConnect =
    {
      SimplicesCell(m_simplex_mesh,simplexInCavity).getNode(0).getGlobalNode(),
      SimplicesCell(m_simplex_mesh,simplexInCavity).getNode(1).getGlobalNode(),
      SimplicesCell(m_simplex_mesh,simplexInCavity).getNode(2).getGlobalNode(),
      SimplicesCell(m_simplex_mesh,simplexInCavity).getNode(3).getGlobalNode()
    };

    allNodes.insert(SimplicesCell(m_simplex_mesh,simplexInCavity).getNode(0).getGlobalNode());
    allNodes.insert(SimplicesCell(m_simplex_mesh,simplexInCavity).getNode(1).getGlobalNode());
    allNodes.insert(SimplicesCell(m_simplex_mesh,simplexInCavity).getNode(2).getGlobalNode());
    allNodes.insert(SimplicesCell(m_simplex_mesh,simplexInCavity).getNode(3).getGlobalNode());

    for(TInt node = 0; node < sizeNodeInSimplex; node++)
    {
        TInt node1 = pointToConnect[(node + 1) % 4];
        TInt node2 = pointToConnect[(node + 2) % 4];
        TInt node3 = pointToConnect[(node + 3) % 4];

        TInt node4 = pointToConnect[(node + 4) % 4];


        ///////
        ///////
        if(std::find_if(nodesInsideCavity.begin(), nodesInsideCavity.end(), [=](const TInt nodeInCav)
        {
          bool nodeInsideCavity = false;
            if(nodeInCav == node1 || nodeInCav == node2 || nodeInCav == node3)
            {
              nodeInsideCavity = true;
            }
            return nodeInsideCavity;

        }) != nodesInsideCavity.end())
        {
          continue;
        }
        ///////
        ///////
        TSimplexID simplexOpp = m_simplex_mesh->getOppositeCell(pointToConnect[node], simplexInCavity);
        if(simplexOpp != errorId)
        {
          if(simplexOpp < 0) // triangle was founded
          {
            std::vector<TSimplexID>&& simplexAdjtoTri = SimplicesTriangle(m_simplex_mesh, -simplexOpp).neighborSimplex();
            simplexAdjtoTri.erase(std::find(simplexAdjtoTri.begin(), simplexAdjtoTri.end(), simplexInCavity),simplexAdjtoTri.end());
            simplexOpp = simplexAdjtoTri[0];
            orderedAdjTet.push_back(simplexOpp);
            cellOfPointToConnect.push_back(simplexInCavity);
          }
          else
          {
            if(node1 != nodeToInsert &&
                node2 != nodeToInsert &&
                  node3 != nodeToInsert )
            {
              if(std::find(cellInCavity().begin(), cellInCavity().end(), simplexOpp) == cellInCavity().end())
              {
                v.push_back(std::vector<TInt>{node1, node2, node3});
                nodesToConnect.insert(node1);
                nodesToConnect.insert(node2);
                nodesToConnect.insert(node3);
                orderedAdjTet.push_back(simplexOpp);
                cellOfPointToConnect.push_back(simplexInCavity);
              }
            }
            if(node1 == nodeToInsert ||
                node2 == nodeToInsert ||
                  node3 == nodeToInsert)
            {
              if(std::find(cellInCavity().begin(), cellInCavity().end(), simplexOpp) == cellInCavity().end())
              {
                complementarySimplex.insert(simplexOpp);
              }
            }
        }
      }
    }
  }

  for(auto const & node : nodesToConnect)
  {
    allNodes.erase(node);
  }
  std::move(allNodes.begin(), allNodes.end(), std::back_inserter(nodeNotConnected));

  return std::move(v);
}
/******************************************************************************/
bool CavityOperator::CavityIO::simplexInborderOfCavity(const TSimplexID simplex, std::vector<TInt>& nodesLocal)
{
  int border =  std::numeric_limits<int>::min();

  const unsigned int cellNodeSize = 4;
  SimplicesCell cell = SimplicesCell(m_simplex_mesh, simplex);
  for(TInt nodeLocal = 0 ; nodeLocal < cellNodeSize ; nodeLocal++)
  {

    TSimplexID adjSimplex = cell.oppositeTetraIdx(nodeLocal);
    if(adjSimplex == border)
    {
      nodesLocal.push_back(nodeLocal);
    }
    else
    {
      if(std::find(cellInCavity().begin(), cellInCavity().end(), adjSimplex) == cellInCavity().end())
      {
        nodesLocal.push_back(nodeLocal);
      }
    }
  }
  return (nodesLocal.size() > 0)?true:false;
}
/*******************************************************************************/
void CavityOperator::CavityIO::findExtSimplex(std::vector<TSimplexID>& extSimplex)
{
  for(auto const & simplex : cellInCavity())
  {
    std::vector<TInt> nodesLocal{};
    if(simplexInborderOfCavity(simplex,  nodesLocal))
    {
      for(auto const & nodeLocal : nodesLocal)
      {
        extSimplex.push_back(SimplicesCell(m_simplex_mesh, simplex).oppositeTetraIdx(nodeLocal));
      }
    }
  }
}
/**********************************************************************/
bool CavityOperator::CavityIO::isTetragonalizableFrom(const TInt nodeToInsert)
{
  bool flag = true;
  std::vector<TSimplexID> extSimplexBorder{};
  std::set<TSimplexID> complementarySimplex{};
  std::vector<TSimplexID> cellOfPointToConnect{};
  std::vector<TSimplexID> nodeNotConnected{};
  const std::vector<TInt>/*&*/ nodesInsideCavity;// = nodeInCavity(nodeToInsert);
  std::vector<std::vector<TInt>>&& pointsToConnect  = pointToConnect(extSimplexBorder, complementarySimplex, nodesInsideCavity, cellOfPointToConnect, nodeNotConnected, nodeToInsert);

  //if(pointsToConnect.size() > 0)
  {
    unsigned int cptCell = 0;
    for(auto const & face : pointsToConnect)
    {
      SimplicesNode node(m_simplex_mesh, nodeToInsert);
      //get The normal vector;
      SimplicesCell cellOfFace = SimplicesCell(m_simplex_mesh, cellOfPointToConnect[cptCell]);
      math::Vector3d normal = cellOfFace.normalOfFace(face);
      math::Vector3d vecFacePt = node.getCoords() - SimplicesNode(m_simplex_mesh, face[0]).getCoords();
      flag = node.isFaceVisible(normal, vecFacePt);

      if(flag == false)
      {
        break;
      }
      ++cptCell;
    }
  }
  //else
  {
    //flag = false;
  }

  return flag;
}
/******************************************************************************/
bool CavityOperator::cavityEnlargement(CavityIO& cavityIO, std::vector<TSimplexID>& initCavityCell, std::vector<TSimplexID>& initCavityTriangle,
                                        const simplicesNode::SimplicesNode& node, const CriterionRAIS& criterion, const std::multimap<TInt, std::pair<TInt, TInt>>& facesAlreadyBuilt,
                                         const std::vector<TSimplexID> markedSimplex)
{
  std::vector<TSimplexID> trianglesConnectedToP{};
  std::vector<TSimplexID> trianglesNotConnectedToP{};

  int errorId =  std::numeric_limits<int>::min();
  gmds::math::Point nextPt = node.getCoords();
  std::vector<TSimplexID> bords;
  std::vector<TSimplexID> simplexToRemove{};
  size_t sizeLocalNodeCell = 4;
  size_t sizeLocalNodeTriangle = 3;

  gmds::BitVector bitMarkedCell(m_simplex_mesh->tetCapacity());
  gmds::BitVector bitCellInCavity(m_simplex_mesh->tetCapacity());
  for(auto const cell : initCavityCell){bitCellInCavity.assign(cell);}
  for(auto const cell : markedSimplex){bitMarkedCell.assign(cell);}
  std::vector<TSimplexID> cavityCell{};
  std::set<int> s( initCavityCell.begin(), initCavityCell.end() );
  cavityCell.assign( s.begin(), s.end() );
  const gmds::BitVector& tetBitVector =  m_simplex_mesh->getBitVectorTet();

  gmds::BitVector triBitVector(m_simplex_mesh->triCapacity());
  for(auto const tri : initCavityTriangle){triBitVector.assign(-tri);}

  Variable<int>* BND_TRIANGLES      = m_simplex_mesh->getVariable<int,SimplicesTriangle>("BND_TRIANGLES");

  Variable<int>* BND_VERTEX_COLOR     = m_simplex_mesh->getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  Variable<int>* BND_CURVE_COLOR      = m_simplex_mesh->getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR    = m_simplex_mesh->getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
  unsigned int dimNode = ((*BND_VERTEX_COLOR)[node.getGlobalNode()] != 0)? SimplexMesh::topo::CORNER : ((*BND_CURVE_COLOR)[node.getGlobalNode()] != 0)? SimplexMesh::topo::RIDGE : ((*BND_SURFACE_COLOR)[node.getGlobalNode()] != 0)? SimplexMesh::topo::SURFACE : SimplexMesh::topo::VOLUME;
  unsigned int indexNode;

  if     (dimNode == SimplexMesh::topo::CORNER ){  indexNode =  (*BND_VERTEX_COLOR )[node.getGlobalNode()] ; }
  else if(dimNode == SimplexMesh::topo::RIDGE  ){  indexNode =  (*BND_CURVE_COLOR  )[node.getGlobalNode()] ; }
  else if(dimNode == SimplexMesh::topo::SURFACE){  indexNode =  (*BND_SURFACE_COLOR)[node.getGlobalNode()] ; }
  else                                          {  indexNode =  0;}

  bool flag = true;
  //CELL EXPANSION
  //for(auto const simplex : cavityCell){std::cout << "initCavCell -> " << simplex << std::endl;}
  while(flag)
  {
    flag = false;
    for(TInt idx = 0 ; idx < cavityCell.size() ; idx++)
    {
      TSimplexID simplexId = cavityCell[idx];
      {
        if(tetBitVector[simplexId] != 0)
        {
          SimplicesCell cell = SimplicesCell(m_simplex_mesh, simplexId);
          //if(!cell.containNode(node))
          {
            for(TInt nodeIndexLocal = 0; nodeIndexLocal < sizeLocalNodeCell; nodeIndexLocal++)
            {
              TSimplexID nextSimplexToAdd = cell.oppositeTetraIdx(nodeIndexLocal);
              if(nextSimplexToAdd != errorId)
              {
                std::vector<TInt> face = cell.getOrderedFace(nodeIndexLocal);
                std::sort(face.begin(), face.end());

                auto it = facesAlreadyBuilt.find(face.front());
                if(it != facesAlreadyBuilt.end())
                {
                  if(it->second.first == face[1] && it->second.second == face[2]) // face can not be pass throught
                  {
                    std::cout << "face  -> " << it->first << " | " << it->second.first << " | " << it->second.second << " already built" << std::endl;
                    return false;
                  }
                }

                if(nextSimplexToAdd >= 0)
                {
                  if(criterion.execute(m_simplex_mesh, simplexId, nodeIndexLocal, nextPt))
                  {
                    if(bitCellInCavity[nextSimplexToAdd] == 0 && bitMarkedCell[nextSimplexToAdd] == 0)
                    {
                      cavityCell.resize(cavityCell.size() + 1);
                      cavityCell[cavityCell.size() - 1] = nextSimplexToAdd;
                      bitCellInCavity.assign(nextSimplexToAdd);
                    }
                  }
                }
                else
                {
                  if(triBitVector[-nextSimplexToAdd] == 0)
                  {
                    initCavityTriangle.push_back(nextSimplexToAdd);
                    triBitVector.assign(-nextSimplexToAdd);
                  }
                }
              }
            }
          }
        }
      }
    }

    gmds::BitVector indexedTriangle(m_simplex_mesh->getBitVectorTri().capacity());
    if(dimNode == SimplexMesh::topo::RIDGE)
    {
      /*std::cout << "SimplexMesh::topo::RIDGE | labelNode--> " << indexNode <<std::endl;
      for(auto const tri : initCavityTriangle)
      {
        std::cout << "triangle --> " << tri <<std::endl;
      }
      for(auto const tet : cavityCell)
      {
        std::cout << "tetra --> " << tet <<std::endl;
      }*/
      struct Edge{
        TInt node0;
        TInt node1;
      };

      const std::map<unsigned int, std::pair<unsigned int, unsigned int>>& mapEdgeTriangleIndices =  m_simplex_mesh->getEdgeTianglesIndices();
      const std::pair<unsigned int, unsigned int>& trianglesIndices = mapEdgeTriangleIndices.at(indexNode);
      unsigned int index0 = trianglesIndices.first;
      unsigned int index1 = trianglesIndices.second;

      BitVector trianglesindexedByIndex0(m_simplex_mesh->getBitVectorTri().capacity());
      BitVector trianglesindexedByIndex1(m_simplex_mesh->getBitVectorTri().capacity());
      for(auto const tri : initCavityTriangle)
      {
        int indexTri = (*BND_TRIANGLES)[-tri];
        if(indexTri == index0){trianglesindexedByIndex0.assign(-tri);}
        else if(indexTri == index1){trianglesindexedByIndex1.assign(-tri);}
      }

      std::vector<Edge> edges{};
      for(unsigned int triangleIndex0ID = 0 ; triangleIndex0ID < trianglesindexedByIndex0.capacity() ; triangleIndex0ID++)
      {
        if(trianglesindexedByIndex0[triangleIndex0ID] == 1)
        {
          SimplicesTriangle triangle(m_simplex_mesh, triangleIndex0ID);
          for(unsigned int edgeLocal = 0 ; edgeLocal < 3 ; edgeLocal++)
          {
            TSimplexID adjTriangle = triangle.neighborTriangle(edgeLocal);
            if(trianglesindexedByIndex1[adjTriangle] == 1)
            {
              std::vector<TInt> edgeNode = triangle.getOppositeEdge(edgeLocal);
              Edge edge;
              edge.node0 = edgeNode.front();
              edge.node1 = edgeNode.back();
              edges.push_back(edge);
            }
          }
        }
      }
      //std::cout << "edges.size() --> " << edges.size() << std::endl;
      double epsilon = 1E-3;
      math::Point pC = node.getCoords();
      Edge e;

      if(edges.size() == 0)
      {
        throw GMDSException("edges.size() == 0, node edges was found for the ridge node being inserted...");
      }
      else if(edges.size() == 1)
      {
        e.node0 = edges.front().node0;
        e.node1 = edges.front().node1;
      }
      else
      {
        for(auto const & edge : edges)
        {
          math::Point pA = SimplicesNode(m_simplex_mesh, edge.node0).getCoords();
          math::Point pB = SimplicesNode(m_simplex_mesh, edge.node1).getCoords();

          math::Vector3d vecAB = pB - pA ;
          math::Vector3d vecAC = pC - pA ;

          if((vecAB.cross(vecAC)).norm() < epsilon)
          {
            double K_AC = vecAB.dot(vecAC);
            double K_AB = vecAB.dot(vecAB);

            if(K_AC >= 0 && K_AC <= K_AB)
            {
              e.node0 = edge.node0;
              e.node1 = edge.node1;
            }
          }
        }
      }

      const std::vector<TSimplexID> shell = SimplicesNode(m_simplex_mesh, e.node0).shell(SimplicesNode(m_simplex_mesh, e.node1));
      std::vector<TSimplexID> firstTriangles{};
      for(auto const simplex : shell){if(simplex < 0){firstTriangles.push_back(-simplex);}}

      if(firstTriangles.size() != 2){
        throw GMDSException("firstTriangles.size() != 2, the surface mesh is broken ...");
      }

      for(auto const & triangle : firstTriangles){
        selectConnexTriangle(triangle, triBitVector, indexedTriangle);
      }
    }

    //set type for the triangle see the paper : Robust Boundary Layer Mesh Generation
    for(auto const tri : initCavityTriangle)
    {
      if(dimNode > SimplexMesh::topo::SURFACE)
      {
        trianglesNotConnectedToP.push_back(tri);
      }
      else
      {
        int indexTri = (*BND_TRIANGLES)[-tri];
        if(dimNode == SimplexMesh::topo::SURFACE)
        {
          if(indexTri == indexNode)
          {
            trianglesConnectedToP.push_back(tri);
          }
          else
          {
            trianglesNotConnectedToP.push_back(tri);
          }
        }
        else if(dimNode == SimplexMesh::topo::RIDGE)
        {
          if(indexedTriangle[-tri] == 1){
            trianglesConnectedToP.push_back(tri);
          }
          else
          {
            trianglesNotConnectedToP.push_back(tri);
          }
        }
      }
    }

    /*for(auto const triNotCo : trianglesNotConnectedToP)
    {
      std::cout << "triNotCo Before --> " << triNotCo << std::endl;
    }
    for(auto const triNotCo : trianglesConnectedToP)
    {
      std::cout << "tri Co Before --> " << triNotCo << std::endl;
    }*/

    std::vector<TSimplexID> triangleToRemoveFromVec{};
    for(auto const triNotCo : trianglesNotConnectedToP)
    {
      const std::vector<TSimplexID> neighborTriNotCo = SimplicesTriangle(m_simplex_mesh, triNotCo).adjacentTriangle();
      TSimplexID neighborTriA = neighborTriNotCo[0];
      TSimplexID neighborTriB = neighborTriNotCo[1];
      TSimplexID neighborTriC = neighborTriNotCo[2];
      TInt indexA = (*BND_TRIANGLES)[neighborTriA];
      TInt indexB = (*BND_TRIANGLES)[neighborTriB];
      TInt indexC = (*BND_TRIANGLES)[neighborTriC];
      TInt index = (*BND_TRIANGLES)[-triNotCo];
      for(auto const triCo : trianglesConnectedToP)
      {
        if((neighborTriA == -triCo && indexA == index) || (neighborTriB == -triCo && indexB == index)|| (neighborTriC == -triCo && indexC == index))
        {
          triangleToRemoveFromVec.push_back(triNotCo);
        }
      }
    }

    for(auto const tri : triangleToRemoveFromVec)
    {
      trianglesNotConnectedToP.erase(std::remove(trianglesNotConnectedToP.begin(), trianglesNotConnectedToP.end(), tri), trianglesNotConnectedToP.end());
      trianglesConnectedToP.push_back(tri);
    }
    //TODO check if all face of trianglesNotConnectedToP (i) produce a valide tetraedron
    //TRIANGLE EXPANSION
    for(TInt idx = 0 ; idx < trianglesConnectedToP.size() ; idx++)
    {
      TSimplexID simplexId = -trianglesConnectedToP[idx];
      SimplicesTriangle triangle = SimplicesTriangle(m_simplex_mesh, simplexId);
      std::vector<TInt> nodes = triangle.getNodes();

      for(TInt nodeIndexLocal = 0; nodeIndexLocal < sizeLocalNodeTriangle; nodeIndexLocal++)
      {
        const TInt nodeA = nodes[nodeIndexLocal];
        const TInt nodeB = nodes[(nodeIndexLocal + 1) % sizeLocalNodeTriangle];
        if(nodeA == node.getGlobalNode() || nodeB == node.getGlobalNode())
        {
          continue;
        }
        else
        {
          const TInt nodeC = nodes[(nodeIndexLocal + 2) % sizeLocalNodeTriangle];
          TSimplexID triangleAdj = triangle.neighborTriangle(triangle.getLocalNode(nodeC));
          if(triBitVector[triangleAdj] == 0) // AB is a boundary Edge
          {
            if(0/*triangle is noty valide*/)
            {
              //add Shell AB in the cavity initCavityCell
              flag = true;
            }
          }
        }
      }
    }

    //we check if the cavity do not pocess any curve discontinuity
    if((*BND_CURVE_COLOR)[node.getGlobalNode()] != 0)
    {/*
      gmds::BitVector triConnectedToPBitVector(m_simplex_mesh->triCapacity());
      for(auto const triangle : trianglesConnectedToP)
      {
        triConnectedToPBitVector.assign(-triangle);
      }

      for(auto const triangle : trianglesConnectedToP)
      {
        unsigned int triangleIndex = (*BND_TRIANGLES)[-triangle];
        for(unsigned int triangleEdgeLocal = 0 ; triangleEdgeLocal < sizeLocalNodeTriangle ; triangleEdgeLocal++)
        {
          TSimplexID adjTriangle =  SimplicesTriangle(m_simplex_mesh, -triangle).neighborTriangle(triangleEdgeLocal);
          if(triangleIndex != (*BND_TRIANGLES)[adjTriangle])
          {
            if(triConnectedToPBitVector[adjTriangle] == 0)
            {
              std::cout << "false" << std::endl;
              return false;
            }
          }
        }
      }*/
    }
  }
  cavityIO.setSimplexCavity(cavityCell, trianglesConnectedToP, trianglesNotConnectedToP);
  return true;
}
/**********************************************************************/
/*void CavityOperator::cavityReduction(CavityIO& cavityIO, std::vector<TSimplexID>& initCavity, const simplicesNode::SimplicesNode& node, const CriterionRAIS& criterion, const CavityReduction& cavityReduction, const std::vector<TSimplexID> v)
{

}*/
/**********************************************************************/
void CavityOperator::selectConnexTriangle(const TSimplexID& firstTriangle, const gmds::BitVector& triangleInCav, gmds::BitVector& connexTriangle)
{
  if(connexTriangle[firstTriangle] == 0)
  {
    connexTriangle.assign(firstTriangle);
    Variable<int>* BND_TRIANGLES = m_simplex_mesh->getVariable<int,SimplicesTriangle>("BND_TRIANGLES");
    unsigned int indexTriangle = (*BND_TRIANGLES)[firstTriangle];
    //std::cout << "  firstTriangle -> " << firstTriangle << std::endl;
    const std::vector<TSimplexID>&& adjSimplex = SimplicesTriangle(m_simplex_mesh, firstTriangle).adjacentTriangle();
    for(auto const & triangle : adjSimplex){
      if(indexTriangle == (*BND_TRIANGLES)[triangle] && triangleInCav[triangle] == 1)
      {
        selectConnexTriangle(triangle, triangleInCav, connexTriangle);
      }
    }
  }
}
/**********************************************************************/
void CavityOperator::fillSelectedIds(const std::vector<TSimplexID>& cavity)
{
  for(auto const & tet : cavity)
  {
    m_tetSelectedIds.assign(tet);
  }
}
/**********************************************************************/
void CavityOperator::clearSelectedIds()
{
  m_tetSelectedIds.clear();
}
/**********************************************************************/
