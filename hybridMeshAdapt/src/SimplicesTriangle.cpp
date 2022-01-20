#include "gmds/hybridMeshAdapt/SimplicesTriangle.h"
#include "gmds/hybridMeshAdapt/SimplexMesh.h"
/*---------------------------------------------------------------------------*/
using namespace gmds;
using namespace hybrid;
using namespace math;
using namespace simplicesCell;
using namespace simplicesTriangle;
using namespace simplicesNode;
/*---------------------------------------------------------------------------*/
TreeTrianglePath::TreeTrianglePath()
{

}
/*---------------------------------------------------------------------------*/
TreeTrianglePath::TreeTrianglePath(SimplexMesh* simplexMesh, BitVector& trianglesAlreadyInTree, TSimplexID currentTriangle, TSimplexID triangleToStop):m_currentTriangle(currentTriangle), m_simplex_mesh(simplexMesh)
{
  const std::vector<TSimplexID> adjSimplices = SimplicesTriangle(m_simplex_mesh, currentTriangle).adjacentTriangle();
  for(auto const adjSimplex : adjSimplices)
  {
    if(trianglesAlreadyInTree[triangleToStop] != 1)
    {
      if(trianglesAlreadyInTree[adjSimplex] == 0)
      {
        m_currentTriangle = adjSimplex;
        trianglesAlreadyInTree.assign(adjSimplex);
        TreeTrianglePath* treePath = new TreeTrianglePath(simplexMesh, trianglesAlreadyInTree, adjSimplex, triangleToStop);
        treePath->m_root = this;
        m_leafs.push_back(treePath);
      }
    }
  }
}
/*---------------------------------------------------------------------------*/
TreeTrianglePath::~TreeTrianglePath()
{
  for(auto const tree : m_leafs)
  {
    if(tree != nullptr)
    {
      delete tree;
    }
  }
}
/*---------------------------------------------------------------------------*/
TreeTrianglePath* TreeTrianglePath::findLeafContaining(const TSimplexID triangleToStop) const
{
  for(auto const & leaf : m_leafs)
  {
    if(leaf != nullptr)
    {
      if(leaf->m_currentTriangle == triangleToStop)
      {
        return leaf;
      }
      return leaf->findLeafContaining(triangleToStop);
    }
  }
}
/*---------------------------------------------------------------------------*/
std::vector<TSimplexID> TreeTrianglePath::pathToTriangle(const TSimplexID triangleToStop) const
{
  std::vector<TSimplexID> path{};
  TreeTrianglePath* treePath = findLeafContaining(triangleToStop);

  if(treePath != nullptr)
  {
    for(;;)
    {
      if(treePath != nullptr)
      {
        path.push_back(treePath->m_currentTriangle);
        treePath = treePath->m_root;
        if(treePath == nullptr)
        {
          break;
        }
      }
    }
  }
  else
  {
    GMDSException("treePath != nullptr, the triangleToStop was not found...");
  }

  return std::move(path);
}
/*---------------------------------------------------------------------------*/
bool TreeTrianglePath::isValidPath(const int indexCurve, const TSimplexID triangleToStop) const
{
  const std::vector<TSimplexID>&& path = pathToTriangle(triangleToStop);
  unsigned int sizeEdges = 3;

  Variable<int>* BND_CURVE_COLOR      = m_simplex_mesh->getVariable<int,SimplicesNode>("BND_CURVE_COLOR");

  for(unsigned int iter = path.size() -1 ; iter == 1 ; iter--)
  {
    TSimplexID currentTriangle  = path[iter];
    TSimplexID previousTriangle = path[iter - 1];

    const SimplicesTriangle tri = SimplicesTriangle(m_simplex_mesh, currentTriangle);
    const std::vector<TSimplexID>&& adjTriangleId         = tri.adjacentTriangle();
    TInt localNode = (adjTriangleId[0] == previousTriangle)?0 : (adjTriangleId[1] == previousTriangle)? 1 : (adjTriangleId[2] == previousTriangle)?2:-1;
    if(localNode == -1)
    {
      GMDSException("local node == -1, during the valid Path algorithm...");
    }

    TInt E0 = (localNode + 1) % sizeEdges;
    TInt E1 = (localNode + 1) % sizeEdges;
    if((*BND_CURVE_COLOR)[E0] != 0 && (*BND_CURVE_COLOR)[E0] != indexCurve){return false;}
    if((*BND_CURVE_COLOR)[E1] != 0 && (*BND_CURVE_COLOR)[E1] != indexCurve){return false;}
  }
  return true;
}
/*---------------------------------------------------------------------------*/
SimplicesTriangle::SimplicesTriangle()
{

}
/*---------------------------------------------------------------------------*/
SimplicesTriangle::~SimplicesTriangle()
{

}
/*---------------------------------------------------------------------------*/
SimplicesTriangle::SimplicesTriangle(SimplexMesh* simplexMesh, const TSimplexID indexTriangle)
{
  if(simplexMesh->m_tri_ids[std::abs(indexTriangle)] != 0)
  {
    m_simplex_mesh = simplexMesh;
    m_simplexId = std::abs(indexTriangle);
  }
  else
  {
    /*TODO exception*/
    std::cout << "Creer le triangle " << indexTriangle <<  " avant de l'utiliser !!" << std::endl;
  }
}

double SimplicesTriangle::getAreaOfCell()
{
  if(m_simplex_mesh != nullptr)
  {
    const Point pt0 = m_simplex_mesh->m_coords[m_simplexId][0];//FAUX je pense
    const Point pt1 = m_simplex_mesh->m_coords[m_simplexId][1];
    const Point  pt2 = m_simplex_mesh->m_coords[m_simplexId][2];

    Vector3d Vec0;
    Vector3d Vec1;

    Vec0.setXYZ(pt0[0] - pt1[0], pt0[1] - pt1[1], pt0[2] - pt1[2]);
    Vec1.setXYZ(pt0[0] - pt2[0], pt0[1] - pt2[1], pt0[2] - pt2[2]);

    const Vector3d OrthVec = Vec0.cross(Vec1);

    return 0.5 * OrthVec.norm();
  }
}

TSimplexID SimplicesTriangle::neighborTriangle (const TInt indexNodeLocal) const
{
  TInt border = std::numeric_limits<int>::min();
  TSimplexID triangle = border;

  if(m_simplex_mesh != nullptr)
  {
    if(!(indexNodeLocal < 0 && indexNodeLocal > 2))
    {
      triangle = m_simplex_mesh->m_tri_adj[m_simplexId][indexNodeLocal];
    }
  }

  return triangle;
}

bool SimplicesTriangle::containNode(const simplicesNode::SimplicesNode& simplicesNode)
{
  if(m_simplex_mesh != nullptr)
  {
    return false;
  }
}

void SimplicesTriangle::reorientTriangle()
{
  if(m_simplex_mesh != nullptr)
  {
      /*TODO*/
  }
}

/*---------------------------------------------------------------------------*/
std::vector<unsigned int> SimplicesTriangle::nodes() const
{
  std::vector<unsigned int> v{};
  v.resize(3);
  v[0] = (unsigned int)m_simplex_mesh->m_tri_nodes[m_simplexId][0];
  v[1] = (unsigned int)m_simplex_mesh->m_tri_nodes[m_simplexId][1];
  v[2] = (unsigned int)m_simplex_mesh->m_tri_nodes[m_simplexId][2];

  return std::move(v);
}
/*---------------------------------------------------------------------------*/
std::vector<TInt> SimplicesTriangle::getNodes() const
{
  std::vector<TInt> v{};
  v.resize(3);
  v[0] = m_simplex_mesh->m_tri_nodes[m_simplexId][0];
  v[1] = m_simplex_mesh->m_tri_nodes[m_simplexId][1];
  v[2] = m_simplex_mesh->m_tri_nodes[m_simplexId][2];

  return std::move(v);
}
/*---------------------------------------------------------------------------*/
math::Vector3d SimplicesTriangle::getNormal()  const
{
  const TInt Node0 = m_simplex_mesh->m_tri_nodes[m_simplexId][0];
  const TInt Node1 = m_simplex_mesh->m_tri_nodes[m_simplexId][1];
  const TInt Node2 = m_simplex_mesh->m_tri_nodes[m_simplexId][2];

  const math::Point pt0 = m_simplex_mesh->m_coords[Node0];
  const math::Point pt1 = m_simplex_mesh->m_coords[Node1];
  const math::Point pt2 = m_simplex_mesh->m_coords[Node2];

  const math::Vector3d vec0 = pt1 - pt0;
  const math::Vector3d vec1 = pt2 - pt1;

  return vec1.cross(vec0);
}
/*---------------------------------------------------------------------------*/
std::vector<TSimplexID> SimplicesTriangle::neighborSimplex() const
{
  std::vector<TSimplexID> v{};
  v.resize(2);
  v[0] = m_simplex_mesh->m_tri_nodes[m_simplexId][3];
  v[1] = m_simplex_mesh->m_tri_adj[m_simplexId][3];
  return std::move(v);
}
/*---------------------------------------------------------------------------*/
std::vector<TSimplexID> SimplicesTriangle::adjacentTriangle() const
{
  std::vector<TSimplexID> v{};
  v.push_back(m_simplex_mesh->m_tri_adj[m_simplexId][0]);
  v.push_back(m_simplex_mesh->m_tri_adj[m_simplexId][1]);
  v.push_back(m_simplex_mesh->m_tri_adj[m_simplexId][2]);

  return std::move(v);
}
/*---------------------------------------------------------------------------*/
std::vector<TSimplexID> SimplicesTriangle::buildclockWiseTrianglesbyShell(const TInt nodeA, const TInt nodeB, const TInt nodeC) const
{
  TInt triIdx = -m_simplexId;
  std::vector<TSimplexID> clockWiseOrdererdTriangles{};
  TInt border = std::numeric_limits<int>::min();

  /*fill the cloWiseVector*/
  TInt firstNodeGlobal = nodeC;
  std::vector<TInt> nodeFace{nodeA, nodeB, firstNodeGlobal};
  TSimplexID currentSimplex = m_simplex_mesh->m_tri_nodes[-triIdx][3];
  TSimplexID previousSimplex = border;
  for(;;)
  {
    if(currentSimplex == border)
    {
      break;
    }
    else if(currentSimplex >= 0)
    {
      TInt localNode  = SimplicesCell(m_simplex_mesh, currentSimplex).getLocalNode(firstNodeGlobal);
      TSimplexID nextSimplex = SimplicesCell(m_simplex_mesh, currentSimplex).oppositeTetraIdx(localNode);

      previousSimplex = currentSimplex;
      std::vector<TInt> interNodeFace{};
      if(nextSimplex != border)
      {
        if(nextSimplex >= 0)
        {
          interNodeFace = SimplicesCell(m_simplex_mesh, currentSimplex).intersectionNodes(SimplicesCell(m_simplex_mesh, nextSimplex));
        }
        else
        {
          interNodeFace = SimplicesTriangle(m_simplex_mesh, -nextSimplex).getNodes();
        }

        if(interNodeFace.size() == 3)
        {
          interNodeFace.erase(std::remove(interNodeFace.begin(), interNodeFace.end(), nodeA));
          interNodeFace.erase(std::remove(interNodeFace.begin(), interNodeFace.end(), nodeB));
          if(interNodeFace.size() == 1)
          {
            firstNodeGlobal = interNodeFace.front();
          }
          else
          {
            std::cout << "interNodeFace.size() == 1" << std::endl;
          }
        }
        else
        {
          std::cout << "interNodeFace.size() == 3" << std::endl;
        }
      }
      currentSimplex  = nextSimplex;
    }
    else if(currentSimplex < 0)
    {
      if(currentSimplex == triIdx)
      {
        break;
      }
      else
      {
        clockWiseOrdererdTriangles.push_back(currentSimplex);
        TSimplexID simplexA = m_simplex_mesh->m_tri_nodes[-currentSimplex][3];
        TSimplexID simplexB = m_simplex_mesh->m_tri_adj[-currentSimplex][3];
        currentSimplex = (simplexA == previousSimplex)? simplexB : simplexA ;
      }
    }
  }
  return std::move(clockWiseOrdererdTriangles);
}
/*---------------------------------------------------------------------------*/
std::vector<TSimplexID> SimplicesTriangle::findclockWiseTrianglesbyShell(const TInt nodeA, const TInt nodeB, const TInt nodeC) const
{
  std::vector<TSimplexID> v{};
  const TInt border = std::numeric_limits<TInt>::min();
  TInt firstNode = nodeC;
  TSimplexID currentSimplex = m_simplexId;

  for(;;)
  {
    SimplicesTriangle cell (m_simplex_mesh, currentSimplex);
    TInt localNode = cell.getLocalNode(firstNode);
    TSimplexID nextSimplex = m_simplex_mesh->m_tri_adj[currentSimplex][localNode];
    if(nextSimplex == m_simplexId || nextSimplex == border)
    {
      break;
    }
    v.push_back(nextSimplex);
    std::vector<TInt> otherNodes = SimplicesTriangle(m_simplex_mesh, nextSimplex).otherNodesInTriangle(cell);

    if(otherNodes.size() == 1)
    {
      firstNode = otherNodes.front();
    }
    else
    {
      //TODO exception
      std::cout << "otherNodes.size() != 1" << std::endl;
    }
    currentSimplex = nextSimplex;
  }

  return std::move(v);
}
/*---------------------------------------------------------------------------*/
std::vector<TInt> SimplicesTriangle::getOtherNodeInSimplex(const std::vector<TInt>& v) const
{
  std::vector<TInt> res{};
  if(m_simplex_mesh != nullptr)
  {
    std::vector<TInt> nodeInTri{m_simplex_mesh->m_tri_nodes[m_simplexId][0],
                                m_simplex_mesh->m_tri_nodes[m_simplexId][1],
                                m_simplex_mesh->m_tri_nodes[m_simplexId][2]};


    for(auto const & nodeIndx : nodeInTri)
    {
      if(std::find(v.begin(), v.end(), nodeIndx) == v.end())
      {
        res.push_back(nodeIndx);
      }
    }
  }

    return std::move(res);
}

  TInt SimplicesTriangle::getLocalNode(const TInt generalIndex) const
  {
    int errorId = std::numeric_limits<int>::min();
    TInt localNode = errorId;
    size_t sizeSimplexNodes = 3;
    for(TInt iter = 0; iter < sizeSimplexNodes ; iter++)
    {
      if(m_simplex_mesh->m_tri_nodes[m_simplexId][iter] == generalIndex)
      {
        localNode = iter;
        break;
      }
    }

    return localNode;
  }

  std::vector<TSimplexID> SimplicesTriangle::neighborTetra() const
  {
    std::vector<TSimplexID> v{};
    TInt border = std::numeric_limits<int>::min();

    if(m_simplex_mesh->m_tri_nodes[m_simplexId][3] != border){v.push_back(m_simplex_mesh->m_tri_nodes[m_simplexId][3]);}
    if(m_simplex_mesh->m_tri_adj[m_simplexId][3] != border){v.push_back(m_simplex_mesh->m_tri_adj[m_simplexId][3]);}

    return std::move(v);
  }

  std::vector<TInt> SimplicesTriangle::getOppositeEdge(const unsigned int localNode) const
  {
    std::vector<TInt> v{};
    unsigned int sizeTriangle = 3;

    if(m_simplex_mesh != nullptr)
    {
      if(!(localNode < 0 || localNode > 2))
      {
        v.reserve(2);
        v.push_back(m_simplex_mesh->m_tri_nodes[m_simplexId][(localNode+ 1) % sizeTriangle]);
        v.push_back(m_simplex_mesh->m_tri_nodes[m_simplexId][(localNode+ 2) % sizeTriangle]);
      }
      else
      {
        //todo exception
      }
    }

    return std::move(v);
  }

  std::vector<TInt> SimplicesTriangle::otherNodesInTriangle(const SimplicesTriangle& simpliceTriangle) const
  {
    std::vector<TInt> v{};
    std::vector<TInt> nodesSimplexA = getNodes();
    std::vector<TInt> nodesSimplexB = simpliceTriangle.getNodes();

    for(auto const &node : nodesSimplexA)
    {
      if(node != nodesSimplexB[0] && node != nodesSimplexB[1] && node != nodesSimplexB[2])
      {
        v.push_back(node);
      }
    }
    return v;
  }
  /*---------------------------------------------------------------------------*/
  math::Orientation::Sign SimplicesTriangle::orientation(const gmds::math::Point& pt, bool inverseOrientation) const
  {
    math::Orientation::Sign sign;
    const TInt Node0 = m_simplex_mesh->m_tri_nodes[m_simplexId][0];
    const TInt Node1 = m_simplex_mesh->m_tri_nodes[m_simplexId][1];
    const TInt Node2 = m_simplex_mesh->m_tri_nodes[m_simplexId][2];

    const math::Point pt0 = m_simplex_mesh->m_coords[Node0];
    const math::Point pt1 = m_simplex_mesh->m_coords[Node1];
    const math::Point pt2 = m_simplex_mesh->m_coords[Node2];
    if(inverseOrientation)
    {
      sign = math::Orientation::orient3d(pt0, pt2, pt1, pt);
    }
    else
    {
      sign = math::Orientation::orient3d(pt0, pt1, pt2, pt);
    }

    return sign;
  }
