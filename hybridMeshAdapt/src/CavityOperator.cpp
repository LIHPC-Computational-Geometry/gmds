/******************************************************************/
#include <gmds/hybridMeshAdapt/CavityOperator.h>
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
#include <gmds/hybridMeshAdapt/ICriterion.h>
#include "gmds/hybridMeshAdapt/PointInsertion.h"
/******************************************************************/
#include <gmds/io/VTKWriter.h>
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
  TInt border = std::numeric_limits<TInt>::min();
  m_cavityCellIn = cavityIn;
  m_cavityTriangleConnectedToP = cavityTriangleConnectedToP;
  m_cavityTriangleNotConnectedToP = cavityTriangleNotConnectedToP;
  m_edgeContainingNode = std::make_pair(border, border);
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
  gmds::BitVector triangleNotCoToP_BitVector(m_simplex_mesh->triCapacity());
  gmds::BitVector nodesBitvector(m_simplex_mesh->nodesCapacity());
  gmds::BitVector borderNodesBitvector(m_simplex_mesh->nodesCapacity());


  std::vector<TInt> allNodes{};
  std::vector<std::vector<unsigned int>> localsBorderNode{};
  std::vector<std::vector<TSimplexID>> oppisitesSimplex{};

  //////////////////////////////////////////////////////////////////////////////
  ////////////////////this section found the interior nodes with cells//////////
  //////////////////////////////////////////////////////////////////////////////

  std::vector<TSimplexID> triangleNotCoToP = getTrianglesNotConnectedToPInCavity();
  for(auto const triangle : triangleNotCoToP)
  {
    if(triangle < 0 && triangle != border)
    {
      triangleNotCoToP_BitVector.assign(-triangle);
    }
  }

  for(auto const & tet : cellInCavity())
  {
    if(tet != border)
    {
      cellBitvector.assign(tet);
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
    else{
      throw gmds::GMDSException("tet == border");
    }
  }
  for(auto const & tet : cellInCavity())
  {
    SimplicesCell cell = SimplicesCell(m_simplex_mesh, tet);
    for(unsigned int localIndex = 0 ; localIndex < sizeCell ; localIndex++)
    {
      //look if the oppositeCell is in the cavity with the help of tetBitVector
      TSimplexID oppositeCell = cell.oppositeTetraIdx(localIndex);
      const std::vector<TInt>&& nodes =  cell.getOrderedFace(localIndex);
      if(!(nodes[0] == node || nodes[1] == node || nodes[2] == node))
      {
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
          if(triangleNotCoToP_BitVector[-oppositeCell] != 0)
          {
            m_nodesToReconnect.push_back(nodes);
            m_oppositeCell.push_back(oppositeCell);
          }
        }
      }
    }
  }
   for(auto const nodes : m_nodesToReconnect)
  {
    for(auto const node : nodes)
    {
      //if(borderNodesBitvector[node] == 0)
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
    else
    {
      throw gmds::GMDSException("tri == border");
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
          std::vector<TSimplexID> borderEdge = triangle.getOppositeEdge(localIndex);

          //condition to look if an indesirable triangle will be built during the cavity's rebuilding (triangle with it's 3 vertex on the same ridge)
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
      else
      {
        throw gmds::GMDSException("oppositeTriangle == border");
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
void CavityOperator::CavityIO::CavityIO::nodesReconnection(const TInt node)
{
  const std::vector<TSimplexID>& base = m_simplex_mesh->getBase();
  Variable<int>* BND_TRIANGLES = m_simplex_mesh->getVariable<int,SimplicesTriangle>("BND_TRIANGLES");

  m_localsNodeForReconnectionWithTriangle.clear();
  m_oppositeTriangle.clear();
  m_triangleIndices.clear();

  unsigned int sizeFace = 3;
  int border = std::numeric_limits<int>::min();
  gmds::BitVector triBitVector(m_simplex_mesh->triCapacity());
  for(auto const & tri : getTrianglesConnectedToPInCavity())
  {
    if(tri != border)
    {
      if(triBitVector[-tri] == 0)
      {
        triBitVector.assign(-tri);
      }
    }
  }

  //looking for the border surface node contain in the m_cavityTriangleConnectedToP vector
  std::vector<TSimplexID> initVec{border, border, border};
  m_localsNodeForReconnectionWithTriangle.resize(m_nodesToReconnect.size(), initVec);
  m_oppositeTriangle.resize(m_nodesToReconnect.size(), initVec);
  m_triangleIndices.resize(m_nodesToReconnect.size(), initVec);

  for(unsigned int triangle = 0 ; triangle < triBitVector.capacity() ; triangle++)
  {
    if(triBitVector[triangle] != 0)
    {
      SimplicesTriangle t     = SimplicesTriangle(m_simplex_mesh, triangle);
      std::vector<TInt> nodes = t.getNodes();
      for(int localNode = 0 ; localNode < sizeFace ; localNode++)
      {
        TInt nodeA = nodes[localNode];
        TInt nodeB = nodes[(localNode + 1) % sizeFace];
        TInt nodeC = nodes[(localNode + 2) % sizeFace];
        if(nodeA != node && nodeB != node)
        {
          TInt lNode = t.getLocalNode(nodeC);
          TSimplexID neighborTri = t.neighborTriangle(lNode);
          if(neighborTri != border)
          {
            if(triBitVector[neighborTri] == 0)
            {
              m_nodesToReconnect_Tri.push_back(std::vector<TInt>{nodeA, nodeB});
              m_oppositeTri.push_back(neighborTri);
              m_indices.push_back((*BND_TRIANGLES)[triangle]);
            }
          }
          else
          {
            throw gmds::GMDSException("neighborTri == border");
          }
        }
      }
    }
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
                                        const simplicesNode::SimplicesNode& node, const CriterionRAIS& criterion, const std::multimap<TInt, TInt>& facesAlreadyBuilt,
                                         const std::vector<TSimplexID> markedSimplex)
{
  auto my_make = [=](TInt nodeA, TInt nodeB){
    return (nodeA > nodeB)? std::make_pair(nodeB, nodeA):std::make_pair(nodeA, nodeB);
  };

  std::map<std::pair<TInt, TInt>, TSimplexID> mapForReinsertion{};
  std::map<std::pair<TInt, TInt>, std::pair<unsigned int, TSimplexID>> mapForTriangleColor{};
  std::vector<TSimplexID> trianglesConnectedToP{};
  std::vector<TSimplexID> trianglesNotConnectedToP{};

  int errorId =  std::numeric_limits<int>::min();
  gmds::math::Point nextPt = node.getCoords();
  std::vector<TSimplexID> bords;
  std::vector<TSimplexID> simplexToRemove{};
  size_t sizeLocalNodeCell = 4;
  size_t sizeLocalNodeTriangle = 3;
  unsigned int sizeEdge = 3;

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
  //CELL EXPANSION START
  while(flag)
  {
    flag = false;
    //std::cout << std::endl;

    for(TInt idx = 0 ; idx < cavityCell.size() ; idx++)
    {
      TSimplexID simplexId = cavityCell[idx];
      if(tetBitVector[simplexId] != 0)
      {
        SimplicesCell cell = SimplicesCell(m_simplex_mesh, simplexId);
        for(TInt nodeIndexLocal = 0; nodeIndexLocal < sizeLocalNodeCell; nodeIndexLocal++)
        {
          TSimplexID nextSimplexToAdd = cell.oppositeTetraIdx(nodeIndexLocal);
          if(nextSimplexToAdd != errorId)
          {
            if(nextSimplexToAdd >= 0)
            {
              //if simplexId is a sliver we add the adjacent tet to the vector
              const SimplicesCell cell = SimplicesCell(m_simplex_mesh, simplexId);
              if(cell.isSliver())
              {
                std::vector<TSimplexID> adjTet = cell.adjacentTetra();
                //std::vector<TSimplexID> adjTet = cell.directConnectedSimplex();
                //std::cout << "sliver -> " << simplexId << std::endl;
                for(auto const simplex : adjTet)
                {
                  //std::cout << "simplex -> " << simplex << std::endl;
                  if(simplex < 0)
                  {
                    if(triBitVector[-simplex] == 0)
                    {
                      initCavityTriangle.push_back(simplex);
                      triBitVector.assign(-simplex);
                    }
                  }
                  else
                  {
                    if(bitCellInCavity[simplex] == 0 && bitMarkedCell[simplex] == 0)
                    {
                      cavityCell.resize(cavityCell.size() + 1);
                      cavityCell[cavityCell.size() - 1] = simplex;
                      bitCellInCavity.assign(simplex);
                    }
                  }
                }
                continue;
              }

              std::vector<TInt> orderedFace = cell.getOrderedFace(nodeIndexLocal);
              if(orderedFace[0] == node.getGlobalNode() || orderedFace[1] == node.getGlobalNode() || orderedFace[2] == node.getGlobalNode())
              {
                //fill adj tets
                std::vector<TInt> adjNodes{};
                std::copy_if(orderedFace.begin(), orderedFace.end(), std::back_inserter(adjNodes), [&] (TInt adjNode){
                  return (adjNode != node.getGlobalNode());
                });
                //reinsertion on the volume
                std::pair<TInt, TInt> p = my_make(adjNodes.front(), adjNodes.back());
                mapForReinsertion.insert(std::pair<std::pair<TInt, TInt>, TSimplexID>(p, nextSimplexToAdd));
                continue;
              }



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
          else
          {
            throw gmds::GMDSException("nextSimplexToAdd == errorId");
          }
        }
      }
      else
      {
        std::cout << "tetBitVector[ " << simplexId << "] == 0" << std::endl;
        throw gmds::GMDSException("unvalid simplex in cavity");
      }
    }

    /*for(unsigned int tri = 0 ; tri < triBitVector.size() ; tri++)
    {
      const TSimplexID oppositeTriangle = SimplicesTriangle(m_simplexMesh, triBitVector[tri]).adjacentTriangle()
      for(auto const oppoTri : oppositeTriangle)
      {
        const SimplicesTriangle triangle(m_simplexMesh, oppoTri);
        if(triangle.isEdge())
        {
          triBitVector.resize(triBitVector.size() + 1);
          triBitVector[triBitVector.size() - 1] = -oppoTri;
        }
      }
    }*/
  }


  //CELL EXPANSION END
  gmds::BitVector indexedTriangle(m_simplex_mesh->getBitVectorTri().capacity());
  if(dimNode == SimplexMesh::topo::CORNER)
  {
    std::vector<TSimplexID> firstTriangles{};
    const std::vector<TSimplexID> ball = node.ballOf();
    for(auto const simplex : ball){if(simplex < 0){firstTriangles.push_back(-simplex);}}

    for(auto const & triangle : firstTriangles){
      selectConnexTriangle(triangle, triBitVector, indexedTriangle, node.getGlobalNode());
    }
  }
  else if(dimNode == SimplexMesh::topo::RIDGE)
  {
    bool alreadyBelongingToAnEdge = false;
    struct Edge{
      TInt node0;
      TInt node1;
    };
    Edge e;
    e.node0 = errorId;
    e.node1 = errorId;

    const std::multimap<TInt, std::pair<TInt,TInt>>& edgeStructure =  m_simplex_mesh->getConstEdgeStructure();
    auto it = edgeStructure.equal_range(indexNode);

    /*test pour voir si le node appartient pas deja a ridge ex : durant une phase de reinsertion*/
    for(auto itr = it.first ; itr != it.second ; itr++)
    {
      TInt nodeEdge0 = itr->second.first;
      TInt nodeEdge1 = itr->second.second;
      if(nodeEdge0 == node.getGlobalNode() || nodeEdge1 == node.getGlobalNode())
      {
        alreadyBelongingToAnEdge = true;
        cavityIO.setNodeInfoEdge(alreadyBelongingToAnEdge);
        break;
      }
    }
    std::vector<TSimplexID> firstTriangles{};
    if(!alreadyBelongingToAnEdge)
    {
      double epsilon = 1E-3;
      math::Point pC = node.getCoords();

      for(auto itr = it.first ; itr != it.second ; itr++)
      {
        math::Point pA = SimplicesNode(m_simplex_mesh, itr->second.first).getCoords();
        math::Point pB = SimplicesNode(m_simplex_mesh, itr->second.second).getCoords();

        math::Vector3d vecAB = pB - pA ;
        math::Vector3d vecAC = pC - pA ;

        if((vecAB.cross(vecAC)).norm() < epsilon)
        {
          double K_AC = vecAB.dot(vecAC);
          double K_AB = vecAB.dot(vecAB);

          if(K_AC >= 0 && K_AC <= K_AB)
          {
            e.node0 = std::min(itr->second.first, itr->second.second);
            e.node1 = std::max(itr->second.first, itr->second.second);
            cavityIO.setEdgeContainingNode(e.node0, e.node1);
            break;
          }
        }
      }

      if(e.node0 == errorId || e.node1 == errorId)
      {
        //the node is not on a curve
        return false;
      }

      const std::vector<TSimplexID> shell = SimplicesNode(m_simplex_mesh, e.node0).shell(SimplicesNode(m_simplex_mesh, e.node1));
      for(auto const simplex : shell){if(simplex < 0){firstTriangles.push_back(-simplex);}}
      /////////////////////////

      ///////////////////////////
      if(firstTriangles.size() != 2){

        for(auto const data : edgeStructure)
        {
          std::cout << "data --> " << data.first << " | [" << data.second.first << " : " << data.second.second << "]" << std::endl;
        }
        std::cout << "e.node0 -> " << e.node0 << std::endl;
        std::cout << "e.node1 -> " << e.node1 << std::endl;
        std::cout << "NODE BEING INSERTED -> " << node.getGlobalNode() << std::endl;
        std::cout << "LABEL -> " << indexNode << std::endl;
        std::cout << "firstTriangles.size() -> " << firstTriangles.size() << std::endl;

        //IF DEBUG MODE
          gmds::ISimplexMeshIOService ioService(m_simplex_mesh);
          gmds::VTKWriter vtkWriterDI(&ioService);
          vtkWriterDI.setCellOptions(gmds::N|gmds::R|gmds::F);
          vtkWriterDI.setDataOptions(gmds::N|gmds::R|gmds::F);
          vtkWriterDI.write("EDGE_BUG_" + std::to_string(node.getGlobalNode()) + ".vtk");

          SimplexMesh nodeMesh = SimplexMesh();
          std::cout << "node.getCoords() -> " << node.getCoords() << std::endl;
          nodeMesh.addNode(node.getCoords()[0], node.getCoords()[1], node.getCoords()[2]);
          nodeMesh.addTetraedre(0, 0, 0, 0);
          gmds::ISimplexMeshIOService ioServiceNode(&nodeMesh);
          gmds::VTKWriter vtkWriterNode(&ioServiceNode);
          vtkWriterNode.setCellOptions(gmds::N|gmds::R);
          vtkWriterNode.setDataOptions(gmds::N|gmds::R);
          vtkWriterNode.write("NODE_" + std::to_string(node.getGlobalNode()) + ".vtk");

        throw GMDSException("firstTriangles.size() != 2, the surface mesh is broken ...");
      }
    }
    else
    {
      const std::vector<TSimplexID> ball = node.ballOf();
      for(auto const simplex : ball){if(simplex < 0){firstTriangles.push_back(-simplex);}}
    }

    for(auto const & triangle : firstTriangles){
      selectConnexTriangle(triangle, triBitVector, indexedTriangle, node.getGlobalNode());
    }
  }
  else if(dimNode == SimplexMesh::topo::SURFACE)
  {
    double distance = std::numeric_limits<TInt>::max();
    std::vector<TSimplexID> firstTriangles{};
    TSimplexID selectedTriangle = std::numeric_limits<TSimplexID>::max();
    math::Point projectedPoint;

    for(auto const tri : initCavityTriangle)
    {
      if(indexNode == (*BND_TRIANGLES)[-tri])
      {
        const SimplicesTriangle triangle(m_simplex_mesh, -tri);
        const std::vector<TInt> nodes = triangle.getNodes();
        const math::Point coord0 = SimplicesNode(m_simplex_mesh, nodes[0]).getCoords();
        const math::Point coord1 = SimplicesNode(m_simplex_mesh, nodes[1]).getCoords();
        const math::Point coord2 = SimplicesNode(m_simplex_mesh, nodes[2]).getCoords();
        double distanceBis = std::numeric_limits<TInt>::max();

        if(m_simplex_mesh->pointInTriangle(node.getCoords(), coord0, coord1, coord2, distanceBis, projectedPoint))
        {
          //std::cout << "  test for triangle -> " << tri << std::endl;
          //std::cout << "  distance -> " << distanceBis << std::endl;

          if(distanceBis < distance)
          {
            selectedTriangle = -tri;
            distance = distanceBis;
          }
        }
      }
    }

    //std::cout << "selected triangle -> " << -selectedTriangle << std::endl;
    if(selectedTriangle != std::numeric_limits<TSimplexID>::max())
    {
      firstTriangles.push_back(selectedTriangle);
    }


    if(firstTriangles.size() == 0){
      return false; //the node was not find on any triangle surface (due to epsilon)
    }

    for(auto const & triangle : firstTriangles){
      selectConnexTriangle(triangle, triBitVector, indexedTriangle, node.getGlobalNode());
    }
  }
  //std::cout << "SETTING THE TRIANGLE " << std::endl;
  //set type for the triangle see the paper : Robust Boundary Layer Mesh Generation
  if(dimNode > SimplexMesh::topo::SURFACE) //if node is on volume
  {
    for(auto const tri : initCavityTriangle)
    {
      trianglesNotConnectedToP.push_back(tri);
    }
  }
  else //node one curve or surface
  {
    for(auto const tri : initCavityTriangle)
    {
      if(indexedTriangle[-tri] == 1){
        //teste sur les node pour voir si le triangle appartient bien a la zone de triangle a detruire
        trianglesConnectedToP.push_back(tri);
      }
      else
      {
        trianglesNotConnectedToP.push_back(tri);
      }
    }
  }
  //////////////////////////////
  //////////////////////////////
  //////////////////////////////

  std::map<TInt, TSimplexID> mapForReinsertionSurface{};
  gmds::BitVector triCoToP_BitVector(m_simplex_mesh->getBitVectorTri().capacity());
  //std::cout << std::endl;
  //std::cout << std::endl;
  for(auto const t : trianglesConnectedToP)
  {
    triCoToP_BitVector.assign(-t);
    //std::cout << "triagnle to delete -> " << t << std::endl;
  }
  for(auto const t : trianglesConnectedToP)
  {
    const SimplicesTriangle triangle = SimplicesTriangle(m_simplex_mesh, t);
    //std::cout << "triangle -> " << t << std::endl;
    for(unsigned int edge = 0; edge < sizeEdge ; edge++)
    {
      TSimplexID oppositeTriangle = triangle.neighborTriangle(edge);
      if(triCoToP_BitVector[oppositeTriangle] != 1)
      {
        TInt nodeA = triangle.getNodes()[(edge + 1) % sizeEdge];
        TInt nodeB = triangle.getNodes()[(edge + 2) % sizeEdge];
        TInt nodeC = triangle.getNodes()[(edge + 3) % sizeEdge];

        /*std::cout << "node " << node.getGlobalNode() << std::endl;
        std::cout << "nodeA -> " << nodeA << std::endl;
        std::cout << "nodeB -> " << nodeB << std::endl;
        std::cout << "NODE_C -> " << nodeC << std::endl;
        std::cout << "oppositeTriangle -> " << oppositeTriangle << std::endl;

        if(nodeA == 6 && nodeB == 50)
        {
          std::cout << SimplicesTriangle(m_simplex_mesh, oppositeTriangle) << std::endl;
          std::cout << SimplicesCell(m_simplex_mesh, SimplicesTriangle(m_simplex_mesh, oppositeTriangle).neighborTetra()[0]) << std::endl;

          gmds::ISimplexMeshIOService ioServiceMesh(m_simplex_mesh);
          gmds::VTKWriter vtkWriterCS(&ioServiceMesh);
          vtkWriterCS.setCellOptions(gmds::N|gmds::R|gmds::F);
          vtkWriterCS.setDataOptions(gmds::N|gmds::R|gmds::F);
          vtkWriterCS.write("ComputeSlicing_EdgeRemove_TEST.vtk");
        }*/
        if(nodeA == node.getGlobalNode() || nodeB == node.getGlobalNode())
        {
          //reinsertion on the surface
          TInt n = (nodeA == node.getGlobalNode())? nodeB: nodeA;
          mapForReinsertionSurface.insert(std::pair<TInt, TSimplexID>(n, oppositeTriangle));
          continue;
        }
        unsigned int indice = (*BND_TRIANGLES)[-t];
        std::pair<unsigned int, TSimplexID> r(indice, oppositeTriangle);
        std::pair<TInt, TInt> p = my_make(nodeA, nodeB);
        mapForTriangleColor.insert(std::pair<std::pair<TInt, TInt>, std::pair<unsigned int, TSimplexID>>(p, r));
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //Check to see if the surface is a variety or not
  std::vector<unsigned int> cpt_VecNodes(m_simplex_mesh->getBitVectorNodes().capacity(),0);
  gmds::BitVector allTriangleInSurface(m_simplex_mesh->getBitVectorTri().capacity());
  for(unsigned int t = 0 ; t < allTriangleInSurface.capacity() ; t++)
  {
    allTriangleInSurface.unselect(t);
  }
  for(auto const t : trianglesConnectedToP)
  {
    allTriangleInSurface.assign(-t);
  }
  for(auto const t : trianglesNotConnectedToP)
  {
    allTriangleInSurface.assign(-t);
  }
  for(unsigned int t = 1 ; t < allTriangleInSurface.capacity() ; t++)
  {
    if(allTriangleInSurface[t] == 1)
    {
      const SimplicesTriangle triangle = SimplicesTriangle(m_simplex_mesh, t);
      for(unsigned int edge = 0; edge < sizeEdge ; edge++)
      {
        TSimplexID oppositeTriangle = triangle.neighborTriangle(edge);
        if(allTriangleInSurface[oppositeTriangle] == 0)
        {
          TInt nodeA = triangle.getNodes()[(edge + 1) % sizeEdge];
          TInt nodeB = triangle.getNodes()[(edge + 2) % sizeEdge];

          cpt_VecNodes[nodeA] += 1;
          cpt_VecNodes[nodeB] += 1;
        }
      }
    }
  }
  //check if the border surface is a variety by checking the cpt
  for(auto const cpt : cpt_VecNodes)
  {
    if(cpt > 2)
    {
      //std::cout << "CPT > 2 --> " << cpt << std::endl;
      return false;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  cavityIO.setTrianglesColor(mapForTriangleColor);
  cavityIO.setReinsertionSurfaceData(mapForReinsertionSurface);
  cavityIO.setReinsertionData(mapForReinsertion);
  cavityIO.setSimplexCavity(cavityCell, trianglesConnectedToP, trianglesNotConnectedToP);

  /*
  for(auto const tet : cavityCell){std::cout << "tet -> " << tet << std::endl;}
  for(auto const triCo : trianglesConnectedToP){std::cout << "triCo -> " << triCo << " | " << (*BND_TRIANGLES)[-triCo] << std::endl;}
  for(auto const triNotCo : trianglesNotConnectedToP){std::cout << "triNotCo -> " << triNotCo << " | " << (*BND_TRIANGLES)[-triNotCo] << std::endl;}
  */

  /*std::cout << "cavityCell.size() -> " << cavityCell.size() << std::endl;
  std::cout << "trianglesConnectedToP.size() -> " << trianglesConnectedToP.size() << std::endl;
  std::cout << "trianglesNotConnectedToP.size() -> " << trianglesNotConnectedToP.size() << std::endl;*/

  return true;
}
/**********************************************************************/
/*void CavityOperator::cavityReduction(CavityIO& cavityIO, std::vector<TSimplexID>& initCavity, const simplicesNode::SimplicesNode& node, const CriterionRAIS& criterion, const CavityReduction& cavityReduction, const std::vector<TSimplexID> v)
{

}*/
/**********************************************************************/
void CavityOperator::selectConnexTriangle(const TSimplexID& firstTriangle, const gmds::BitVector& triangleInCav, gmds::BitVector& connexTriangle, const TInt nodeToConnect)
{
  if(connexTriangle[firstTriangle] == 0)
  {
    const TInt border = std::numeric_limits<TInt>::min();
    connexTriangle.assign(firstTriangle);
    Variable<int>* BND_TRIANGLES = nullptr;
    Variable<int>* BND_CURVE_COLOR = nullptr;
    Variable<int>* BND_VERTEX_COLOR = nullptr;
    Variable<int>* BND_SURFACE_COLOR = nullptr;
    try{
        BND_TRIANGLES = m_simplex_mesh->getVariable<int,SimplicesTriangle>("BND_TRIANGLES");
        BND_CURVE_COLOR = m_simplex_mesh->getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
        BND_VERTEX_COLOR = m_simplex_mesh->getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
        BND_SURFACE_COLOR = m_simplex_mesh->getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
    }catch(gmds::GMDSException e)
    {
      throw gmds::GMDSException(e);
    }

    unsigned int sizeTriangle = 3;
    unsigned int indexTriangle = (*BND_TRIANGLES)[firstTriangle];
    unsigned int dimNode    = ((*BND_VERTEX_COLOR)[nodeToConnect] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[nodeToConnect] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[nodeToConnect] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
    unsigned int labelNode  = (dimNode == SimplexMesh::topo::CORNER)?(*BND_VERTEX_COLOR)[nodeToConnect]:(dimNode == SimplexMesh::topo::RIDGE)?(*BND_CURVE_COLOR)[nodeToConnect]:(dimNode == SimplexMesh::topo::SURFACE)?(*BND_SURFACE_COLOR)[nodeToConnect]:0;

    const SimplicesTriangle currentTri(m_simplex_mesh, firstTriangle);
    //adjSimplex reprensent the adjacent local's ordered triangles 0, 1, 2
    //adjEdges represent the adjacent local's ordered edges 0, 1, 2
    const std::vector<TSimplexID>&& adjSimplex = currentTri.adjacentTriangle();

    std::vector<std::vector<TInt>> adjEdges{};
    adjEdges.push_back(currentTri.getOppositeEdge(0));
    adjEdges.push_back(currentTri.getOppositeEdge(1));
    adjEdges.push_back(currentTri.getOppositeEdge(2));

    for(unsigned int idx = 0 ; idx < sizeTriangle ; idx++)
    {
      const std::vector<TInt> edge = adjEdges[idx];
      if(edge.size() == 2)
      {
        if(edge.front() != border && edge.back() != border)
        {
          //if(((*BND_SURFACE_COLOR)[edge.front()] != 0 || (*BND_SURFACE_COLOR)[edge.back()] != 0) ||
          //  ((*BND_CURVE_COLOR)[edge.front()] != (*BND_CURVE_COLOR)[edge.back()] && (*BND_CURVE_COLOR)[edge.front()] != 0 && (*BND_CURVE_COLOR)[edge.back()] != 0))
          {
            const TSimplexID triangle = adjSimplex[idx];
            if(indexTriangle == (*BND_TRIANGLES)[triangle] && triangleInCav[triangle] == 1)
            {
              selectConnexTriangle(triangle, triangleInCav, connexTriangle, nodeToConnect);
            }
          }

          /////
          /*const TSimplexID triangle = adjSimplex[idx];
          if(dimNode == SimplexMesh::topo::SURFACE)
          {
            if(indexTriangle == (*BND_TRIANGLES)[triangle] && triangleInCav[triangle] == 1)
            {
              selectConnexTriangle(triangle, triangleInCav, connexTriangle, nodeToConnect);
            }
          }
          else if(dimNode == SimplexMesh::topo::RIDGE)
          {
            const std::map<unsigned int, std::pair<unsigned int, unsigned int>>& C2S = m_simplex_mesh->getEdgeTianglesIndices();

            std::pair<unsigned int, unsigned int> p{};
            try{
              p = C2S.at(labelNode);
            }catch(gmds::GMDSException e)
            {
              throw gmds::GMDSException(e);
            }

            if(((*BND_TRIANGLES)[triangle] == p.first || (*BND_TRIANGLES)[triangle] == p.second) && triangleInCav[triangle] == 1)
            {
              selectConnexTriangle(triangle, triangleInCav, connexTriangle, nodeToConnect);
            }
          }*/
          /////
        }
        else
        {
          throw gmds::GMDSException("edge.first() != border && edge.back() != border");
        }
      }
      else
      {
        throw gmds::GMDSException("edge.size() != 2, bad connex triangle ");
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
