#include <chrono>
/*****************************************************************************/
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
#include <gmds/hybridMeshAdapt/ICriterion.h>
#include <gmds/hybridMeshAdapt/PointInsertion.h>
#include <gmds/hybridMeshAdapt/EdgeInsertion.h>
#include <map>
/*****************************************************************************/
using namespace gmds;
using namespace gmds::hybrid;
using namespace math;
using namespace simplicesCell;
using namespace simplicesNode;
using namespace simplicesTriangle;
using namespace operators;
//#define EPSILON 0.000001
#define EPSILON 10E-5
/*---------------------------------------------------------------------------*/
SimplexMesh::SimplexMesh()
{
  const size_t initSize = 0;
  const TInt border = std::numeric_limits<TInt>::min();
  std::vector<TInt> borders{border, border, border, border};

  m_node_ids = gmds::BitVector(initSize);
  m_tet_ids = gmds::BitVector(initSize);
  m_tri_ids = gmds::BitVector(initSize);

  //init the tet tri and node size with the same size of ids ...
  m_coords.resize(initSize);
  m_tet_nodes.resize(initSize);
  m_tri_nodes.resize(initSize);

  m_node_variable_manager = new VariableManager();
  m_tet_variable_manager = new VariableManager();
  m_tri_variable_manager = new VariableManager();

  if(m_tri_ids.top() + 1 >=  m_tri_ids.capacity())
  {
    //resize of the m_tri_node vector..
    const TInt previousCapacityBitVector = m_tri_ids.capacity() + 1;
    m_tri_nodes.resize(previousCapacityBitVector*2, borders);
    m_tri_adj.resize(previousCapacityBitVector*2, borders);
  }
  TInt idx0 = m_tri_ids.getFreeIndex();
  m_tri_ids.assign(idx0);
  m_tri_variable_manager->addEntry(idx0);

}
/*---------------------------------------------------------------------------*/
SimplexMesh::SimplexMesh(const SimplexMesh& simplexMesh)
{
  nodeIndx     = simplexMesh.nodeIndx;
  m_node_ids   = simplexMesh.m_node_ids;
  m_coords     = simplexMesh.m_coords;
  m_base       = simplexMesh.m_base;

  tetIndx      = simplexMesh.tetIndx;
  m_tet_ids    = simplexMesh.m_tet_ids;
  m_tet_nodes  = simplexMesh.m_tet_nodes;
  m_tet_adj    = simplexMesh.m_tet_adj;

  triIndx      = simplexMesh.triIndx;
  m_tri_ids    = simplexMesh.m_tri_ids;
  m_tri_nodes  = simplexMesh.m_tri_nodes;
  m_tri_adj    = simplexMesh.m_tri_adj;

  /*creer un operateur d'affectation dans variable manager?*/
  m_node_variable_manager     = simplexMesh.m_node_variable_manager;
  m_tet_variable_manager      = simplexMesh.m_tet_variable_manager;
  m_tri_variable_manager = simplexMesh.m_tri_variable_manager;
}
/*---------------------------------------------------------------------------*/
SimplexMesh& SimplexMesh::operator=(const SimplexMesh& simplexMesh)
{
  if(this != &simplexMesh)
  {
    nodeIndx     = simplexMesh.nodeIndx;
    m_node_ids   = simplexMesh.m_node_ids;
    m_coords     = simplexMesh.m_coords;
    m_base       = simplexMesh.m_base;

    tetIndx      = simplexMesh.tetIndx;
    m_tet_ids    = simplexMesh.m_tet_ids;
    m_tet_nodes  = simplexMesh.m_tet_nodes;
    m_tet_adj    = simplexMesh.m_tet_adj;

    triIndx      = simplexMesh.triIndx;
    m_tri_ids    = simplexMesh.m_tri_ids;
    m_tri_nodes  = simplexMesh.m_tri_nodes;

    m_tri_adj    = simplexMesh.m_tri_adj;

    /*creer un operateur d'affectation dans variable manager?*/
    m_node_variable_manager     = simplexMesh.m_node_variable_manager;
    m_tet_variable_manager      = simplexMesh.m_tet_variable_manager;
    m_tri_variable_manager = simplexMesh.m_tri_variable_manager;
  }
  return *this;
}
/*---------------------------------------------------------------------------*/
SimplexMesh::SimplexMesh(SimplexMesh&& simplexMesh)
{
  nodeIndx     = simplexMesh.nodeIndx;
  m_node_ids   = simplexMesh.m_node_ids;
  m_coords     = simplexMesh.m_coords;
  m_base       = simplexMesh.m_base;

  tetIndx      = simplexMesh.tetIndx;
  m_tet_ids    = simplexMesh.m_tet_ids;
  m_tet_nodes  = simplexMesh.m_tet_nodes;
  m_tet_adj    = simplexMesh.m_tet_adj;

  triIndx      = simplexMesh.triIndx;
  m_tri_ids    = simplexMesh.m_tri_ids;
  m_tri_nodes  = simplexMesh.m_tri_nodes;
  m_tri_adj    = simplexMesh.m_tri_adj;

  /*creer un operateur d'affectation dans variable manager?*/
  m_node_variable_manager     = simplexMesh.m_node_variable_manager;
  m_tet_variable_manager      = simplexMesh.m_tet_variable_manager;
  m_tri_variable_manager = simplexMesh.m_tri_variable_manager;

  /*Clear all the simplexMesh instance data*/
  simplexMesh.clear();


  /*creer une fonction qui clear le vector dans variable manager ?*/
  /*simplexMesh.m_node_variable_manager;
  simplexMesh.m_tet_variable_manager;
  simplexMesh.m_triangle_variable_manager;
  */
}
/*---------------------------------------------------------------------------*/
void SimplexMesh::clear()
{
  nodeIndx = 0;
  m_node_ids.clear();
  m_coords.clear();
  m_base.clear();

  tetIndx = 0;
  m_tet_ids.clear();
  m_tet_nodes.clear();
  m_tet_adj.clear();

  triIndx = 0;
  m_tri_ids.clear();
  m_tri_nodes.clear();
  m_tri_adj.clear();

  const size_t initSize = 0;
  m_node_ids = gmds::BitVector(initSize);
  m_tet_ids = gmds::BitVector(initSize);
  m_tri_ids = gmds::BitVector(initSize);

  //init the tet tri and node size with the same size of ids ...
  m_coords.resize(initSize);
  m_tet_nodes.resize(initSize);
  m_tri_nodes.resize(initSize);
}
/*****************************************************************************/
SimplexMesh::~SimplexMesh()
{
  /*if(m_node_variable_manager != nullptr)
  {
      delete m_node_variable_manager;
  }

  if(m_tet_variable_manager != nullptr)
  {
      delete m_tet_variable_manager;
  }

  if(m_triangle_variable_manager != nullptr)
  {
      delete m_triangle_variable_manager;
  }*/
}
/*****************************************************************************/
void SimplexMesh::deleteAllSimplicesBut(const std::vector<TSimplexID> & simplices)
{
  for(unsigned int tet = 0 ; tet < m_tet_ids.capacity(); tet++)
  {
    if(std::find(simplices.begin(), simplices.end(), tet) == simplices.end())
    {
      m_tet_ids.unselect(tet);
    }
  }

  /*for(unsigned int tri = 1 ; tri < m_tri_ids.capacity(); tri++)
  {
    if(m_tri_ids[tri] != 0)
    {
      if(std::find(simplices.begin(), simplices.end(), -tri) == simplices.end())
      {
        m_tri_ids.unselect(tri);
      }
    }
  }*/
}
/*****************************************************************************/
void SimplexMesh::deleteAllTrianglesBut(const std::vector<TSimplexID> & triangles)
{
  for(unsigned int tri = 0 ; tri < m_tri_ids.capacity(); tri++)
  {
    if(std::find(triangles.begin(), triangles.end(), tri) == triangles.end())
    {
      m_tri_ids.unselect(tri);
    }
  }
}
/*****************************************************************************/
void SimplexMesh::deleteAllTriangle()
{
  for(unsigned int tri = 1 ; tri < m_tri_ids.capacity(); tri++)
  {
    if(m_tri_ids[tri] != 0)
    {
      m_tri_ids.unselect(tri);
    }
  }
}
/*****************************************************************************/
template <size_t N, size_t M>
bool SimplexMesh::containNodes(const TSimplexID* array1, const TSimplexID* array2)
{
  size_t compt = 0;
  for (size_t nodeArrayLocal1 = 0; nodeArrayLocal1 < N; nodeArrayLocal1++)
  {
    for (size_t nodeArrayLocal2 = 0; nodeArrayLocal2 < M; nodeArrayLocal2++)
    {
      if(array1[nodeArrayLocal1 % N] == array2[nodeArrayLocal2 % M])
      {
        ++compt;
      }
    }
  }
  return (compt == N)?true:false;
}
/*****************************************************************************/
template<typename T>
void SimplexMesh::setComparator(std::vector<T>& set1, const std::vector<T>& set2)
{
  for(auto const & valuesConst : set2)
  {
    if(set1.size() != 0)
    {
      set1.erase(std::remove(set1.begin(), set1.end(), valuesConst), set1.end());
    }
    else
    {
      break;
    }
  }
}
/*****************************************************************************/
TInt SimplexMesh::getNbSimplex()
{
  return getNbTriangle() + getNbTetra();
}
/*****************************************************************************/
TInt SimplexMesh::getNbNodes()
{
  TInt cpt = 0;
  for(TInt idx = 0; idx < m_node_ids.capacity(); idx++)
  {
    if(m_node_ids[idx] == 1)
    {
      cpt++;
    }
  }
  return cpt;
}
/*****************************************************************************/
TInt SimplexMesh::getNbTriangle()
{
  TInt cpt = 0;
  for(TInt idx = 0; idx < m_tri_ids.capacity(); idx++)
  {
    if(m_tri_ids[idx] == 1)
    {
      cpt++;
    }
  }
  return cpt;
}
/*---------------------------------------------------------------------------*/
TInt SimplexMesh::getNbTetra()
{
  TInt cpt = 0;
  for(TInt idx = 0; idx < m_tet_ids.capacity(); idx++)
  {
    if(m_tet_ids[idx] == 1)
    {
      cpt++;
    }
  }
  return cpt;
}
/*---------------------------------------------------------------------------*/
TSimplexID SimplexMesh::getOppositeCell(const TInt ANodeID, const TSimplexID ATetID)
{
  const int border  = std::numeric_limits<int>::min();
  size_t localIndex = SimplicesCell(this, ATetID).getLocalNode(ANodeID);
  //si la valeur retournée est < 0 --> c'est la valeur negative de l'index d'un triangle dans m_tri_nodes..
  if(localIndex == border)
  {
    //si il n'existe ni triangle ni tetra opposé cherché
    return border;
  }
  else
  {
    //if(m_tet_adj[ATetID] != nullptr)
    TSimplexID adjTet = m_tet_adj[ATetID][localIndex];

    if(adjTet != border)
    {
      if(adjTet >= 0)
      {
        if(m_tet_ids[adjTet] != 0)
        {
            return adjTet;
        }
      }
      else
      {
        if(m_tri_ids[-adjTet] != 0)
        {
          return adjTet;
        }
      }
    }
    else
    {
      return border;
    }
  }
}
/*---------------------------------------------------------------------------*/
TSimplexID SimplexMesh::getSimplexFromBase(const TInt ANodeID)
{
  return m_base[ANodeID];
}
/*---------------------------------------------------------------------------*/
TInt SimplexMesh::addNodeAndcheck(const math::Point& pt, std::vector<TSimplexID>& tetraContainingPt, bool& alreadyAdd, TSimplexID simplexToCheckFirst)
{
  bool flag = false;
  //check if the point pt is on any other node already in the mesh;
  if(tetraContainingPt.size() == 0)
  {
    flag = checkSimplicesContenaing(pt, tetraContainingPt, simplexToCheckFirst);
  }

  if(flag)
  {
    double epsilon = 10E-4;
    for(auto const & tet : tetraContainingPt)
    {
      const std::vector<TInt>& nodesOfTet = SimplicesCell(this, tet).getNodes();
      for(auto const nodeOfTet : nodesOfTet)
      {
        math::Vector3d VecBetweenPt = (SimplicesNode(this, nodeOfTet).getCoords() - pt);
        double lenght = VecBetweenPt.norm();
        if(lenght < epsilon)//pt is on nodeOfTet
        {
          alreadyAdd = true;
          return nodeOfTet;
        }
      }
    }
  }
  std::cout << "near Node ? ended" << std::endl;

  if(m_node_ids.top() + 1 >= m_node_ids.capacity() )
  {
    //resize of the m_node vector..
    const TInt previousCapacityBitVector = m_node_ids.capacity() + 1;
    m_coords.resize(previousCapacityBitVector*2);
    m_base.resize(previousCapacityBitVector*2, std::numeric_limits<int>::min()); //ajout
  }
  const TInt idx = m_node_ids.getFreeIndex();
  m_node_ids.assign(idx);
  m_coords[idx] = pt;
  m_node_variable_manager->addEntry(idx);

  return idx;
}
/*---------------------------------------------------------------------------*/
TInt SimplexMesh::addNode(const math::Point& pt)
{
  //check if the point pt is on any other node already in the mesh;
  TSimplexID tetraContainingPt;
  bool flag = simplexContaining(pt, tetraContainingPt);
  if(flag)
  {
      double epsilon = 10E-4;
      const std::vector<TInt>& nodesOfTet = SimplicesCell(this, tetraContainingPt).getNodes();
      for(auto const nodeOfTet : nodesOfTet)
      {
        math::Vector3d VecBetweenPt = (SimplicesNode(this, nodeOfTet).getCoords() - pt);
        double lenght = VecBetweenPt.norm();
        if(lenght <  epsilon)//pt is on nodeOfTet
        {
            return nodeOfTet;
        }
      }
  }

  if(m_node_ids.top() + 1 >= m_node_ids.capacity() )
  {
    //resize of the m_node vector..
    const TInt previousCapacityBitVector = m_node_ids.capacity() + 1;
    m_coords.resize(previousCapacityBitVector*2);
    m_base.resize(previousCapacityBitVector*2, std::numeric_limits<int>::min()); //ajout
  }
  const TInt idx = m_node_ids.getFreeIndex();

  m_node_ids.assign(idx);
  m_coords[idx] = pt;
  //m_node_variable_manager->addEntry(idx);
  return idx;
}
/*****************************************************************************/
TInt SimplexMesh::addNode(const Point&& pt)
{

  //check if the point pt is on any other node already in the mesh;
  TSimplexID tetraContainingPt;
  bool flag = simplexContaining(pt, tetraContainingPt);
  if(flag)
  {
      double epsilon = 10E-4;
      const std::vector<TInt>& nodesOfTet = SimplicesCell(this, tetraContainingPt).getNodes();
      for(auto const nodeOfTet : nodesOfTet)
      {
        math::Vector3d VecBetweenPt = (SimplicesNode(this, nodeOfTet).getCoords() - pt);
        double lenght = VecBetweenPt.norm();
        if(lenght <  epsilon)//pt is on nodeOfTet
        {
            return nodeOfTet;
        }
      }
  }

  if(m_node_ids.top() + 1 >= m_node_ids.capacity() )
  {

    //resize of the m_node vector..
    const TInt previousCapacityBitVector = m_node_ids.capacity() + 1;
    m_coords.resize(previousCapacityBitVector*2);
    m_base.resize(previousCapacityBitVector*2, std::numeric_limits<int>::min()); //ajout
  }
  const TInt idx = m_node_ids.getFreeIndex();
  m_node_ids.assign(idx);
  m_coords[idx] = pt;
  //TODO//
  //m_node_variable_manager->addEntry(idx);

  return idx;
}
/*---------------------------------------------------------------------------*/
TInt SimplexMesh::addNode(const TCoord X, const TCoord Y, const TCoord Z)
{
  return addNode(Point(X, Y, Z));
}
/*---------------------------------------------------------------------------*/
bool SimplexMesh::deleteNode(const SimplicesNode& simpliceNode, bool eraseNode)
{
  return deleteNode(simpliceNode.getGlobalNode(), eraseNode);
}
/*****************************************************************************/
bool SimplexMesh::deleteNode(const TInt indexNode, bool eraseNode)
{
  bool flag = false;
  TInt border = std::numeric_limits<int>::min();
  if(indexNode > m_node_ids.capacity() || indexNode < 0)
  {
    /*TODO mettre une exeption a la place d'un bool*/
    /*Node index  <0 || > m_node_ids.capacity()*/
    flag = false;
  }
  else
  {
    if(m_node_ids[indexNode] != 0)
    {
      const std::vector<TSimplexID>&& ballOfNode = SimplicesNode(this, indexNode).ballOf();
      /*Dans un premier temps on delete les tetra liée au node indexNode*/
      for(const auto& simplex : ballOfNode)
      {
        if(simplex != border)
        {
          if(simplex >= 0)
          {
            if(m_tet_ids[simplex] != 0)
            {
                deleteTetra(simplex);
            }
          }
          else
          {
            if(m_tri_ids[-simplex] != 0)
            {
                deleteTetra(-simplex);
            }
          }
        }
      }

      if(eraseNode)
      {
        m_node_ids.unselect(indexNode);
        m_node_variable_manager->removeEntry(indexNode);
      }
      flag = true;
    }
  }
  return flag;
}
/*---------------------------------------------------------------------------*/
TSimplexID SimplexMesh::addTriangle(const simplicesNode::SimplicesNode& ANode0,
                       const simplicesNode::SimplicesNode& ANode1,
                       const simplicesNode::SimplicesNode& ANode2,
                       bool flag)
{
  int border = std::numeric_limits<int>::min();
  TSimplexID idx;
  if(m_tri_ids.top() + 1 >= m_tri_ids.capacity())
  {
    //resize of the m_tri_node vector..
    const TInt previousCapacityBitVector = m_tri_ids.capacity() + 1;
    m_tri_nodes.resize(previousCapacityBitVector*2);
    m_tri_adj.resize(previousCapacityBitVector*2);
  }

  idx = m_tri_ids.getFreeIndex();
  m_tri_ids.assign(idx);

  std::vector<TInt> triangle{ANode0.getGlobalNode(), ANode1.getGlobalNode(), ANode2.getGlobalNode(), border};
  std::vector<TInt> triangleAdj{border, border, border, border};


  m_tri_nodes[idx] = triangle;
  m_tri_adj[idx]   = triangleAdj;
  m_tri_variable_manager->addEntry(idx);

  //TODO//
  m_tri_variable_manager->addEntry(idx);

  if(flag)
  {
    buildTriBaseAndAdjLocal(idx);
    buildOppFaces(idx);
  }

  return idx;
}


/*---------------------------------------------------------------------------*/
TSimplexID SimplexMesh::addTriangle(const TInt AIndexPoint_0,
                       const TInt AIndexPoint_1,
                       const TInt AIndexPoint_2,
                        bool flag)
{
  return addTriangle(simplicesNode::SimplicesNode(this, AIndexPoint_0),
                     simplicesNode::SimplicesNode(this, AIndexPoint_1),
                     simplicesNode::SimplicesNode(this, AIndexPoint_2),
                     flag);
}
/*****************************************************************************/
std::vector<TInt> SimplexMesh::deleteTriangle(const TInt ATriIndex)
{
  TSimplexID ATriangleIndex = std::abs(ATriIndex);
  std::vector<TInt> nodesTri;
  TInt border = std::numeric_limits<int>::min();
  if(ATriangleIndex == 0 || ATriangleIndex > m_tri_ids.capacity())
  {
    //TODO
    //Exception
  }
  else
  {
    if(m_tri_ids[ATriangleIndex] != 0)
    {
      /*reconstruction of the m_tri_adj */
      std::vector<std::vector<TSimplexID>> directNeighborSimplex{};
      for(unsigned int localNode = 0 ; localNode < 3 ; localNode++)
      {
        TInt nodeA = m_tri_nodes[ATriangleIndex][localNode];
        TInt nodeB = m_tri_nodes[ATriangleIndex][(localNode + 1) % 3];
        TInt nodeC = m_tri_nodes[ATriangleIndex][(localNode + 2) % 3];
        std::vector<TSimplexID> clockWiseTri = SimplicesTriangle(this, ATriangleIndex).findclockWiseTrianglesbyShell(nodeA, nodeB, nodeC);
        directNeighborSimplex.push_back(clockWiseTri);
      }
      SimplicesTriangle currentTriangle = SimplicesTriangle(this, ATriangleIndex);

      for(auto const clockTri : directNeighborSimplex)
      {
        if(clockTri.size() > 0)
        {
          if(clockTri.size() == 1)
          {
            TSimplexID triangleA = clockTri.back();
            std::vector<TInt> otherNode = SimplicesTriangle(this, triangleA).otherNodesInTriangle(currentTriangle);
            if(otherNode.size() == 1)
            {
              TInt node = otherNode.front();
              TInt localNode = SimplicesTriangle(this, triangleA).getLocalNode(node);
              m_tri_adj[triangleA][localNode] = border;
            }
            else
            {
              //TODO exception
              std::cout << "otherNode.size() != 1" << std::endl;
            }
          }
          else
          {
            //////////////////////
            TSimplexID triangleA = clockTri.front();
            TSimplexID triangleB = clockTri.back();

            std::vector<TInt> otherNode = SimplicesTriangle(this, triangleB).otherNodesInTriangle(currentTriangle);
            if(otherNode.size() == 1)
            {
              TInt node = otherNode.front();
              TInt localNode = SimplicesTriangle(this, triangleB).getLocalNode(node);
              m_tri_adj[triangleB][localNode] = triangleA;
            }
            else
            {
              //TODO exception
              std::cout << "otherNode.size() != 1" << std::endl;
            }
            //////////////////////
          }
        }
      }
      //not useful
      m_tri_adj[ATriangleIndex][0] = border;
      m_tri_adj[ATriangleIndex][1] = border;
      m_tri_adj[ATriangleIndex][2] = border;
      //
      int border = std::numeric_limits<int>::min();
      nodesTri = std::vector<TSimplexID>{m_tri_nodes[ATriangleIndex][0], m_tri_nodes[ATriangleIndex][1], m_tri_nodes[ATriangleIndex][2]};
      //On recontruit les vector adj des tetra adj a ce triangle.
      TSimplexID tetra0 = m_tri_nodes[ATriangleIndex][3];
      TSimplexID tetra1 = m_tri_adj[ATriangleIndex][3];

      std::vector<TInt> nodesTet0NotInNodesTri;
      std::vector<TInt> nodesTet1NotInNodesTri;

      if(tetra0 != border)
      {
        nodesTet0NotInNodesTri = SimplicesCell(this, tetra0).getOtherNodeInSimplex(nodesTri);
      }
      if(tetra1 != border)
      {
        nodesTet1NotInNodesTri = SimplicesCell(this, tetra1).getOtherNodeInSimplex(nodesTri);
      }


      if(nodesTet0NotInNodesTri.size() == 1)
      {
        if(m_tet_ids[tetra0] != 0)
        {
          TInt nodeAdjTet0 = nodesTet0NotInNodesTri[0];
          TInt nodeAdjTet0Local = SimplicesCell(this, tetra0).getLocalNode(nodeAdjTet0);
          if(nodeAdjTet0Local < 4)
          {
              m_tet_adj[tetra0][nodeAdjTet0Local] = tetra1;
          }
          else
          {
            //TODO exception
            std::cout << "nodeAdjTet0Local > 4" << std::endl;
          }
        }
        else
        {
          //TODO exeption
          std::cout << "nodesTet0NotInNodesTri.size() == 1" << std::endl;
        }
      }
      if(nodesTet1NotInNodesTri.size() == 1 )
      {
        if(m_tet_ids[tetra1] != 0)
        {
          TInt nodeAdjTet1 = nodesTet1NotInNodesTri[0];
          TInt nodeAdjTet1Local = SimplicesCell(this, tetra1).getLocalNode(nodeAdjTet1);
          if(nodeAdjTet1Local < 4)
          {
              m_tet_adj[tetra1][nodeAdjTet1Local] = tetra0;
          }
          else
          {
            //TODO exception
          }
        }
        else
        {
          //TODO exeption
        }
      }
    }
  }
  m_tri_ids.unselect(ATriangleIndex);
  return std::move(nodesTri);
}
/*---------------------------------------------------------------------------*/
TSimplexID SimplexMesh::addTetraedre(const simplicesNode::SimplicesNode&& ANode0,
                                    const simplicesNode::SimplicesNode&& ANode1,
                                    const simplicesNode::SimplicesNode&& ANode2,
                                    const simplicesNode::SimplicesNode&& ANode3,
                                    const bool buildAdjInfo)
{
  TInt border = std::numeric_limits<int>::min();
  std::vector<TInt> borders{border, border, border, border};

  if(m_tet_ids.top() + 1 >=  m_tet_ids.capacity())
  {
    //resize of the m_tri_node vector..
    const TInt previousCapacityBitVector = m_tet_ids.capacity() + 1;
    m_tet_nodes.resize(previousCapacityBitVector*2, borders);
    m_tet_adj.resize(previousCapacityBitVector*2, borders);
  }

  const TInt idx = m_tet_ids.getFreeIndex();
  m_tet_ids.assign(idx);

  /*Building the array of the tri_node before push it back to the m_tri_nodes*/
  std::vector<TSimplexID> tetra{ANode0.getGlobalNode(), ANode1.getGlobalNode(), ANode2.getGlobalNode(), ANode3.getGlobalNode()};

  //TODO//
  //m_tet_variable_manager->addEntry(idx);

  /*rebuild the base vector after ading aTetra*/
  m_tet_nodes[idx] = tetra;

  std::vector<TSimplexID> ptrAdj{border, border, border, border};
  m_tet_adj[idx] = ptrAdj;

  reorientTetra(idx); /*voir si on peut reorienté que le tetra crée*/

  buildBaseLocal(idx);

  if(buildAdjInfo)
  {
    buildTetBaseAndAdjLocal(idx);
  }
  m_tet_variable_manager->addEntry(idx);
  return idx;
}
/*---------------------------------------------------------------------------*/
void SimplexMesh::reorientTetra(const TSimplexID & tetIndx) /*voir si on peut reorienté que le tetra crée*/
{
  if(m_tet_ids[tetIndx] != 0)
  {
    simplicesCell::SimplicesCell(this, tetIndx).reorientTet();
  }
}
/*---------------------------------------------------------------------------*/
void SimplexMesh::buildBaseLocal(const TSimplexID& tetIndx)
{
  int errorId = std::numeric_limits<int>::min();
  unsigned int sizeNodeInTet = 4;

  for(unsigned int nodeLocal = 0; nodeLocal < sizeNodeInTet; nodeLocal++)
  {
    TInt nodeGlobal = SimplicesCell(this, tetIndx).getNode(nodeLocal).getGlobalNode();
    if(m_base[nodeGlobal] == errorId)
    {
      m_base[nodeGlobal] = tetIndx;
    }
  }
}
/*---------------------------------------------------------------------------*/
void SimplexMesh::buildTriBaseAndAdjLocal(const TSimplexID & triIndx)
{

  TInt border = std::numeric_limits<int>::min();
  const std::vector<TInt> nodes = m_tri_nodes[triIndx];

  const SimplicesNode node0 = SimplicesNode(this, nodes[0]);
  const SimplicesNode node1 = SimplicesNode(this, nodes[1]);
  const SimplicesNode node2 = SimplicesNode(this, nodes[2]);

  std::vector<TSimplexID>&& ballNode0 = node0.ballOf();
  std::vector<TSimplexID>&& ballNode1 = node1.ballOf();
  std::vector<TSimplexID>&& ballNode2 = node2.ballOf();

  std::vector<std::vector<TSimplexID>> ball012{ballNode0, ballNode1, ballNode2};
  std::vector<TSimplexID> adjSimplex = intersectionSimplex(ball012);

  if(adjSimplex.size() > 2)
  {
    //TODO assertion
    std::cout << "adjSimplex.size() > 2 ---> " << adjSimplex.size() << std::endl;
  }

  for(auto const simplex : adjSimplex)
  {
    SimplicesCell currentSimplex = SimplicesCell(this, simplex);
    std::vector<TInt> otherNode = currentSimplex.getOtherNodeInSimplex(nodes);
    if(otherNode.size() == 1)
    {
      TInt localNode =  currentSimplex.getLocalNode(otherNode.front());
      std::vector<TInt> orderedNodes = currentSimplex.getOrderedFace(localNode);
      bool flag = false;

      for(unsigned int nodeIdx = 0 ; nodeIdx < 3 ; nodeIdx++)
      {
        if(nodes[0] == orderedNodes[nodeIdx] && nodes[1] == orderedNodes[(nodeIdx + 1) % 3])
        {
          flag = true;
          m_tri_nodes[triIndx][3] = simplex;
          m_tet_adj[simplex][localNode] = -triIndx;
          break;
        }
      }
      if(flag == false)
      {
        m_tri_adj[triIndx][3] = simplex;
        m_tet_adj[simplex][localNode] = -triIndx;
      }
    }
    else
    {
      //TODO exception
      std::cout << "otherNode.size() != 1 -> " << otherNode.size() << std::endl;
    }
  }
}
/*---------------------------------------------------------------------------*/
void SimplexMesh::buildTetBaseAndAdjLocal(const TSimplexID & tetIndx)
{
  SimplicesCell currentCell = SimplicesCell(this, tetIndx);
  TInt border = std::numeric_limits<int>::min();

  TInt nodeA = m_tet_nodes[tetIndx][0];
  TInt nodeB = m_tet_nodes[tetIndx][1];
  TInt nodeC = m_tet_nodes[tetIndx][2];
  TInt nodeD = m_tet_nodes[tetIndx][3];

  std::vector<TInt> nodes{nodeA, nodeB, nodeC, nodeD};

  std::vector<TSimplexID> ballA = SimplicesNode(this, nodeA).ballOf();
  std::vector<TSimplexID> ballB = SimplicesNode(this, nodeB).ballOf();
  std::vector<TSimplexID> ballC = SimplicesNode(this, nodeC).ballOf();
  std::vector<TSimplexID> ballD = SimplicesNode(this, nodeD).ballOf();

  std::vector<std::vector<TSimplexID>> balls{ballA, ballB, ballC, ballD};
  for(unsigned int ballNode = 0 ; ballNode < balls.size() ; ballNode++)
  {
    std::vector<TSimplexID> ball0 = balls[ballNode];
    std::vector<TSimplexID> ball1 = balls[(ballNode + 1) % balls.size()];
    std::vector<TSimplexID> ball2 = balls[(ballNode + 2) % balls.size()];

    std::vector<std::vector<TSimplexID>> ball012{ball0, ball1, ball2};
    std::vector<TSimplexID>&& intersectionBall012 = intersectionSimplex(ball012);
    std::sort (intersectionBall012.begin(), intersectionBall012.end());

    if(intersectionBall012.size() != 0)
    {
      TSimplexID simplex = intersectionBall012.front();
      if(simplex != tetIndx)
      {
        TInt nodeA = nodes[ballNode];
        TInt nodeB = nodes[(ballNode + 1) % balls.size()];
        TInt nodeC = nodes[(ballNode + 2) % balls.size()];
        const std::vector<TInt>  vecNode{nodeA, nodeB, nodeC};

        if(simplex >= 0)
        {
          SimplicesCell cell = SimplicesCell(this, simplex);
          std::vector<TInt> nodesCurrentSimplex = currentCell.getOtherNodeInSimplex(vecNode);
          std::vector<TInt> nodesSimplex = cell.getOtherNodeInSimplex(vecNode);
          if(nodesCurrentSimplex.size() == 1 && nodesSimplex.size() == 1)
          {
            TInt localNodeCurrentSimplex = currentCell.getLocalNode(nodesCurrentSimplex.front());
            TInt localNodeSimplex = cell.getLocalNode(nodesSimplex.front());

            m_tet_adj[currentCell.simplexId()][localNodeCurrentSimplex] = cell.simplexId();
            m_tet_adj[cell.simplexId()][localNodeSimplex] = currentCell.simplexId();
          }
          else
          {
            std::cout << "!(nodesCurrentSimplex.size() == 1 && nodesSimplex.size() == 1)" << std::endl;
          }
        }
        else
        {
          SimplicesTriangle cell = SimplicesTriangle(this, simplex);
          std::vector<TInt> nodesCurrentSimplex = currentCell.getOtherNodeInSimplex(vecNode);
          if(nodesCurrentSimplex.size() == 1)
          {
            TInt globalNode = nodesCurrentSimplex.front();
            TInt localNodeCurrentSimplex = currentCell.getLocalNode(globalNode);
            m_tet_adj[currentCell.simplexId()][localNodeCurrentSimplex] = - cell.simplexId();
          }

          if(m_tri_nodes[-simplex][3] == border && m_tri_adj[-simplex][3] != border)
          {
            m_tri_nodes[-simplex][3] = currentCell.simplexId();
          }
          else if(m_tri_nodes[-simplex][3] != border && m_tri_adj[-simplex][3] == border)
          {
            m_tri_adj[-simplex][3] = currentCell.simplexId();
          }
          else if(m_tri_nodes[-simplex][3] == border && m_tri_adj[-simplex][3] == border)
          {
            std::cout << "m_tri_nodes[-simplex][3] == border && m_tri_adj[-simplex][3] == border" << std::endl;
          }
          else
          {
            std::cout << "m_tri_nodes[-simplex][3] != border && m_tri_adj[-simplex][3] != border" << std::endl;
          }
        }
      }
    }
  }
}
/******************************************************************************/
void SimplexMesh::buildAdjInfoGlobal()
{
    std::vector<std::vector<TSimplexID>> nodesToSimplex(m_node_ids.size());
    for(unsigned int tetIter = 0; tetIter < m_tet_ids.capacity() ; tetIter++)
    {
      if(m_tet_ids[tetIter])
      {
        SimplicesCell cell = SimplicesCell(this, tetIter);
        std::vector<TInt>&& nodes = cell.getNodes();
        nodesToSimplex[nodes[0]].push_back(tetIter);
        nodesToSimplex[nodes[1]].push_back(tetIter);
        nodesToSimplex[nodes[2]].push_back(tetIter);
        nodesToSimplex[nodes[3]].push_back(tetIter);
      }
    }

    for(unsigned int tetIter = 0 ; tetIter < m_tet_ids.capacity() ; tetIter++)
    {
      if(m_tet_ids[tetIter])
      {
        TInt border = std::numeric_limits<int>::min();
        std::vector<TSimplexID> ptrAdj{border, border, border, border};

        SimplicesCell cell = SimplicesCell(this, tetIter);
        std::vector<TInt>&& nodes = cell.getNodes();
        std::vector<TSimplexID> ballNode0 = nodesToSimplex[nodes[0]];
        std::vector<TSimplexID> ballNode1 = nodesToSimplex[nodes[1]];
        std::vector<TSimplexID> ballNode2 = nodesToSimplex[nodes[2]];
        std::vector<TSimplexID> ballNode3 = nodesToSimplex[nodes[3]];

        SimplicesNode ANode0 = SimplicesNode(this, nodes[0]);
        SimplicesNode ANode1 = SimplicesNode(this, nodes[1]);
        SimplicesNode ANode2 = SimplicesNode(this, nodes[2]);
        SimplicesNode ANode3 = SimplicesNode(this, nodes[3]);

        std::vector<std::vector<TSimplexID>> ballsNode{ballNode0, ballNode1, ballNode2, ballNode3};
        std::vector<SimplicesNode>           Anodes{ANode0, ANode1, ANode2, ANode3};
        TInt sizeNodeInTet = 4;


        for(unsigned int ballIter = 0; ballIter < sizeNodeInTet; ballIter++)
        {
          std::vector<TSimplexID> ball = ballsNode[ballIter];
          SimplicesNode          Anode = Anodes[ballIter];

          if(ball.size() != 0)
          {
            for(auto const & simplexInBall : ball)
            {
              std::vector<TInt>&&  vecNode = SimplicesCell(this, tetIter).intersectionNodes(SimplicesCell(this, simplexInBall));
              if(vecNode.size() == 3)
              {
                std::vector<TInt> nodes              = SimplicesCell(this, tetIter).getOtherNodeInSimplex(vecNode);
                std::vector<TInt> nodesSimplexInBall = SimplicesCell(this, simplexInBall).getOtherNodeInSimplex(vecNode);

                if(nodes.size() == 1 && nodesSimplexInBall.size() == 1)
                {
                  TInt localNode              = SimplicesCell(this, tetIter).getLocalNode(nodes[0]);
                  TInt localNodeSimplexInBall = SimplicesCell(this, simplexInBall).getLocalNode(nodesSimplexInBall[0]);
                  if(m_tet_adj[simplexInBall][localNodeSimplexInBall] < 0)
                  {
                    if(m_tet_adj[simplexInBall][localNodeSimplexInBall] == border)
                    {
                      m_tet_adj[tetIter][localNode] = simplexInBall;
                      m_tet_adj[simplexInBall][localNodeSimplexInBall] = tetIter;
                    }
                    else
                    {
                      TSimplexID triangleID =  - m_tet_adj[simplexInBall][localNodeSimplexInBall];
                      ptrAdj[localNode] = m_tet_adj[simplexInBall][localNodeSimplexInBall];
                      bool flagAdj = (m_tri_adj[triangleID][3] == border) ? true : false;
                      (flagAdj) ? m_tri_adj[triangleID][3] = tetIter : m_tri_nodes[triangleID][3] = tetIter;
                    }
                  }
                  else
                  {
                      //ptrAdj[localNode] = simplexInBall;
                      //m_tet_adj[simplexInBall][localNodeSimplexInBall] = tetIter;
                  }
                }
              }
            }
            //m_tet_adj[tetIter] = ptrAdj;
          }
        }
      }
    }
}
/******************************************************************************/
void SimplexMesh::buildOppFaces(const TSimplexID idxTri)
{
  TInt triIdx = -idxTri;
  TInt border = std::numeric_limits<int>::min();
  TInt node0 = m_tri_nodes[-triIdx][0];
  TInt node1 = m_tri_nodes[-triIdx][1];
  TInt node2 = m_tri_nodes[-triIdx][2];

  std::vector<TInt> nodes{node0, node1, node2};
  for(TInt edge = 0 ; edge < nodes.size() ; edge++)
  {
    TInt nodeA = nodes[edge];
    TInt nodeB = nodes[(edge + 1) % nodes.size()];
    TInt nodeC = nodes[(edge + 2) % nodes.size()];
    std::vector<TSimplexID> clockWiseOrdererdTriangles = SimplicesTriangle(this, -idxTri).buildclockWiseTrianglesbyShell(nodeA, nodeB, nodeC);

    /*reconstruction of the adj component*/
    if(clockWiseOrdererdTriangles.size()  > 0)
    {
      TSimplexID simplexA = clockWiseOrdererdTriangles.front();
      TSimplexID simplexB = clockWiseOrdererdTriangles.back();

      SimplicesTriangle triangle  = SimplicesTriangle(this, -triIdx);
      SimplicesTriangle triangleA = SimplicesTriangle(this, -simplexA);
      SimplicesTriangle triangleB = SimplicesTriangle(this, -simplexB);

      //build m_tet_adj for triangle to triangle A
      std::vector<TInt> gobalNodeVec = triangle.otherNodesInTriangle(triangleA);
      if(gobalNodeVec.size() == 1)
      {
        TInt localNode = triangle.getLocalNode(gobalNodeVec.front());
        m_tri_adj[-triIdx][localNode] = triangleA.simplexId();
      }
      else
      {
        std::cout << "globalNodeVec.size() != 1" << std::endl;
      }


      //build m_tet_adj for triangleB to triangle
      gobalNodeVec = triangleB.otherNodesInTriangle(triangle);
      if(gobalNodeVec.size() == 1)
      {
        TInt localNode = triangleB.getLocalNode(gobalNodeVec.front());
        m_tri_adj[-simplexB][localNode] = triangle.simplexId();
      }
      else
      {
        std::cout << "globalNodeVec.size() != 1" << std::endl;
      }

    }
  }
}
/*---------------------------------------------------------------------------*/
TSimplexID SimplexMesh::addTetraedre(const TInt AIndexPoint0,
                                    const TInt AIndexPoint1,
                                    const TInt AIndexPoint2,
                                    const TInt AIndexPoint3,
                                    const bool buildAdjInfo)
{
  return addTetraedre(simplicesNode::SimplicesNode(this, AIndexPoint0),
                      simplicesNode::SimplicesNode(this, AIndexPoint1),
                      simplicesNode::SimplicesNode(this, AIndexPoint2),
                      simplicesNode::SimplicesNode(this, AIndexPoint3),
                      buildAdjInfo);
}
/******************************************************************************/
bool SimplexMesh::checkMeshLocal(const SimplicesNode node)
{
  TInt border = std::numeric_limits<int>::min();
  const std::vector<TSimplexID>&& ball = node.ballOf();
  for(auto const tet : ball)
  {
    SimplicesCell currentSimplex = SimplicesCell(this, tet);
    TInt currentLocalNode = currentSimplex.getLocalNode(node.getGlobalNode());
    for(unsigned int localNode = 0 ; localNode < 4 ; localNode++)
    {
      if(currentLocalNode != localNode)
      {
        TSimplexID adjSimplexId = currentSimplex.oppositeTetraIdx(localNode);
        if(adjSimplexId != border)
        {
          SimplicesCell adjSimplex = SimplicesCell(this, adjSimplexId);
          const std::vector<TInt>&& intersectionNode = currentSimplex.intersectionNodes(adjSimplex);
          if(intersectionNode.size() == 3 /*&& (intersectionNode)*/)
          {
            const std::vector<TInt>&& nodeAdjCell= adjSimplex.getOtherNodeInSimplex(intersectionNode);
            if(nodeAdjCell.size() == 1)
            {
              TInt localNodeAdjSimplex = adjSimplex.getLocalNode(nodeAdjCell.front());
              if(adjSimplex.oppositeTetraIdx(localNodeAdjSimplex) != tet)
              {
                return false;
              }
            }
            else
            {
              return false;
            }
          }
          else
          {
            return false;
          }
        }
      }
    }
  }

  return true;
}
/******************************************************************************/

bool SimplexMesh::checkMesh()
{
  TInt border = std::numeric_limits<int>::min();
  bool flagRes = true;;
/*  for(unsigned int iter = 0; iter < m_tet_ids.capacity(); iter++)
  {
    if(m_tet_ids[iter] != 0)
    {
      std::vector<TSimplexID>&& neighborTet = SimplicesCell(this, iter).adjacentTetra();
      unsigned int nodeLocal = 0;
      for(auto const & tet : neighborTet)
      {
        const unsigned int nodeGlobal = SimplicesCell(this, iter).getNode(nodeLocal).getGlobalNode();
        if(tet != border)
        {
          flagRes = false;
          std::vector<TInt>&& intersectedNodes = SimplicesCell(this, iter).intersectionNodes(SimplicesCell(this, tet));
          if(intersectedNodes.size() == 3)
          {
            std::vector<TInt>&& adjNodes = SimplicesCell(this, tet).getNodes();
            adjNodes.erase(std::remove_if(adjNodes.begin(), adjNodes.end(), [&](TInt node)
            {
            bool flag = false;
            if(std::find(intersectedNodes.begin(), intersectedNodes.end(), node) != intersectedNodes.end())
            {
              flag = true;
            }
            return flag;
            }), adjNodes.end());
            if(adjNodes.size() == 1)
            {
              TSimplexID checkSimplex = SimplicesCell(this, tet).oppositeTetraIdx(SimplicesNode(this, adjNodes.front()));
              if(checkSimplex != iter)
              {
                std::cout << "checkSimplex != iter" << std::endl;
                std::cout << "checkSimplex : " << checkSimplex << "  tet : " << iter  << std::endl;
              }
              else
              {
                flagRes = true;
              }
            }
            else
            {
              std::cout << "adjNodes.size() != 1" << std::endl;
            }
          }
          else
          {
            std::cout << "Simplices [" << iter << ";" << tet << "] are adjacent but don't have 3 intersected node" << std::endl;
          }


          if(flagRes == false)
          {
            break;
          }
        }
        else
        {
          for(unsigned int tet = 0; tet < m_tet_ids.capacity(); tet++)
          {
            if(tet != iter)
            {
              if(m_tet_ids[tet] == 1)
              {
                SimplicesCell cell(this, iter);
                SimplicesCell cellToCompare(this, tet);
                const std::vector<TInt>&& intersectionNode = cell.intersectionNodes(cellToCompare);
                if(intersectionNode.size() == 3)
                {
                  if(std::find(intersectionNode.begin(), intersectionNode.end(), nodeGlobal) == intersectionNode.end())
                  {
                    std::cout << "Simplices [" << tet << ";" << iter << "] are adjacent but adjcent vector is wrong" << std::endl;
                  }
                }
              }
            }
          }
        }
        nodeLocal++;
      }
    }
    if(flagRes == false)
    {
      break;
    }
  }
*/
  unsigned int tetraNbr = getNbTetra();
  unsigned int cpt = 0;
  for(unsigned int cell0 = 0; cell0 < m_tet_ids.capacity(); cell0++)
  {
    if(m_tet_ids[cell0] == 1)
    {
      cpt++;
      std::cout << "avancement : " << (float)cpt / (float)tetraNbr * 100.0f<< " %" << std::endl;
      for(unsigned int cell1 = 0; cell1 < m_tet_ids.capacity(); cell1++)
      {
        if(m_tet_ids[cell1] == 1)
        {
          if(cell0 != cell1)
          {
              SimplicesCell cell(this, cell0);
              SimplicesCell cellToCompare(this, cell1);
              const std::vector<TInt>&& intersectionNode = cell.intersectionNodes(cellToCompare);
              if(intersectionNode.size() == 3)
              {
                const std::vector<TInt>&& compNodeCell0 = cell.getOtherNodeInSimplex(intersectionNode);
                const std::vector<TInt>&& compNodeCell1 = cellToCompare.getOtherNodeInSimplex(intersectionNode);
                if(compNodeCell0.size() == 1 && compNodeCell1.size() == 1)
                {
                  SimplicesNode node0 = SimplicesNode(this, compNodeCell0.back());
                  SimplicesNode node1 = SimplicesNode(this, compNodeCell1.back());

                  TInt localNode0 = cell.getLocalNode(node0.getGlobalNode());
                  TInt localNode1 = cellToCompare.getLocalNode(node1.getGlobalNode());
                  if(m_tet_adj[cell0][localNode0] != cell1 || m_tet_adj[cell1][localNode1] != cell0)
                  {
                      flagRes = false;
                      std::cout << "Simplices [" << cell0 << ";" << cell1 << "] are adjacent but adjcent vector is wrong" << std::endl;
                  }
                }
              }
            }
          }
        }
        if(!flagRes)
        {
          break;
        }
      }
    }

  return flagRes;
}
/******************************************************************************/
bool SimplexMesh::doCellExist(const TSimplexID simplex) const
{
  if(simplex >= 0 || simplex < m_tet_ids.capacity())
  {
    return (m_tet_ids[simplex] != 0);
  }
  return false;
}
/******************************************************************************/
bool SimplexMesh::doNodeExist(const TInt node) const
{
  if(node >= 0 || node < m_node_ids.capacity())
  {
    return (m_node_ids[node] == 1);
  }
  return false;
}
/******************************************************************************/
void SimplexMesh::buildSimplexHull()
{
  TSimplexID border = std::numeric_limits<int>::min();
  struct borderSimplex{
    TSimplexID simplexId;
    TInt localNode;
  };
  std::vector<borderSimplex> borderSimplices;

  /*this loop find all the border tetraedron with the local node pointing the border of the mesh*/
  for(unsigned int tet = 0 ; tet < m_tet_ids.capacity() ; tet++)
  {
    if(m_tet_ids[tet] != 0)
    {
      for(unsigned int localNode = 0 ; localNode < 4 ; localNode++)
      {
        TSimplexID oppositeSimplex = m_tet_adj[tet][localNode];//SimplicesCell(this, tet).oppositeTetraIdx(localNode);
        if(oppositeSimplex == border)
        {
          borderSimplex bS;
          bS.simplexId = tet;
          bS.localNode = localNode;
          borderSimplices.push_back(bS);
        }
      }
    }
  }

  /*Building the hull with the borderSimplices*/
  for(auto const borderSimplex : borderSimplices)
  {
    const std::vector<TInt>&& faceNodes = SimplicesCell(this, borderSimplex.simplexId).getOrderedFace(borderSimplex.localNode);
    TSimplexID triIdx = addTriangle(faceNodes[0], faceNodes[1], faceNodes[2]);
  }

  //for(auto const tet : m_tet_ids){deleteTetra(tet);}
  Variable<int>* BND_VERTEX_COLOR = getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  Variable<int>* BND_CURVE_COLOR = getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR = getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  if(!(BND_VERTEX_COLOR == nullptr || BND_CURVE_COLOR == nullptr || BND_SURFACE_COLOR == nullptr))
  {
    gmds::Variable<int>* BND_TRIANGLES = newVariable<int, SimplicesTriangle>("BND_TRIANGLES");
    for(unsigned int tri = 1 ; tri < m_tri_ids.capacity() ; tri++)
    {
      if(m_tri_ids[tri] == 1)
      {
        std::vector<TInt> nodes = SimplicesTriangle(this, tri).getNodes();
        //0 --> corner | 1 --> ridge | 2--> surface | 3 --> volume node
        unsigned int dim0 = ((*BND_VERTEX_COLOR)[nodes[0]] != 0)? CORNER : ((*BND_CURVE_COLOR)[nodes[0]] != 0)?  RIDGE : ((*BND_SURFACE_COLOR)[nodes[0]] != 0)? SURFACE : VOLUME;
        unsigned int dim1 = ((*BND_VERTEX_COLOR)[nodes[1]] != 0)? CORNER : ((*BND_CURVE_COLOR)[nodes[1]] != 0)?  RIDGE : ((*BND_SURFACE_COLOR)[nodes[1]] != 0)? SURFACE : VOLUME;
        unsigned int dim2 = ((*BND_VERTEX_COLOR)[nodes[2]] != 0)? CORNER : ((*BND_CURVE_COLOR)[nodes[2]] != 0)?  RIDGE : ((*BND_SURFACE_COLOR)[nodes[2]] != 0)? SURFACE : VOLUME;
        int index;

        unsigned int dimMax = std::max(dim0, std::max(dim1, dim2));
        unsigned int dimMin = std::min(dim0, std::min(dim1, dim2));

        if(dimMax == SURFACE)
        {
          TInt surfaceNode = (dimMax == dim0)? nodes[0] : (dimMax == dim1)? nodes[1] : nodes[2];
          index = (*BND_SURFACE_COLOR)[surfaceNode];
          (*BND_TRIANGLES)[tri] = index;
        }
        else
        {
          if(dimMin == CORNER)
          {
            SimplicesTriangle triangle = SimplicesTriangle(this, tri);
            TInt cornerNode = (dimMin == dim0)? nodes[0] : (dimMin == dim1)? nodes[1] : nodes[2];
            TInt localNode   = triangle.getLocalNode(cornerNode);
            /*look for the adjacent surface triangle ---> to change if some triangle will be added on the volume mesh*/
            TSimplexID adjTriangle = triangle.neighborTriangle(localNode);
            std::vector<TInt> nodesAdjTriangle = SimplicesTriangle(this, adjTriangle).otherNodesInTriangle(triangle);
            if(nodesAdjTriangle.size() == 1)
            {
              index = (*BND_SURFACE_COLOR)[nodesAdjTriangle.front()];

              if(index != 0)
              {
                (*BND_TRIANGLES)[tri] = index;
              }
            }
          }
        }
      }
    }

    for(unsigned int tri = 1 ; tri < m_tri_ids.capacity() ; tri++)
    {
      if(m_tri_ids[tri] == 1)
      {
        if((*BND_TRIANGLES)[tri] == 0)
        {
          gmds::BitVector cyclingCheck(m_tri_ids.capacity());
          cyclingCheck.assign(tri);
          TInt idx = 0;
          do {
            idx = findRemainTriangleIdx(tri, cyclingCheck);
          } while(idx == 0 || idx > 150);


          (*BND_TRIANGLES)[tri] = idx;
        }
      }
    }
    //filling edgeTianglesIndices
    Variable<int>* BND_CURVE_COLOR = getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
    for(unsigned int nodeIdx = 0 ; nodeIdx < m_node_ids.capacity() ; nodeIdx++)
    {
      if(m_node_ids[nodeIdx] != 0)
      {
        if((*BND_CURVE_COLOR)[nodeIdx] != 0)
        {
          unsigned int indiceNode      = (*BND_CURVE_COLOR)[nodeIdx];
          std::vector<TSimplexID> ball = SimplicesNode(this, nodeIdx).ballOf();
          std::cout << "ball -> " << ball.size() << std::endl;
          std::set<unsigned int> trianglesIndices{};
          if(edgeTianglesIndices.find(indiceNode) == edgeTianglesIndices.end())
          {
            for(auto const simplex : ball)
            {
              if(simplex < 0 && simplex != border)
              {
                unsigned int triangleIndice = (*BND_TRIANGLES)[-simplex];
                trianglesIndices.insert(triangleIndice);
              }
            }

            if(trianglesIndices.size() == 2)
            {
              std::set<unsigned int>::iterator it = trianglesIndices.begin();
              unsigned int triangleId0 = *it;
              std::advance(it, 1);
              unsigned int triangleId1 = *it;
              std::pair<unsigned int, unsigned int> pairTrianglesIndedices = std::make_pair(triangleId0, triangleId1);
              edgeTianglesIndices[indiceNode] = pairTrianglesIndedices;
            }
            else
            {
              std::cout << "trianglesIndices.size() != 2 --> " << trianglesIndices.size() << " for node : "<< nodeIdx << " of index : " << indiceNode <<  std::endl;
            }
          }
        }
      }
    }
  }
}
/******************************************************************************/
TInt SimplexMesh::findRemainTriangleIdx(const TInt tri, gmds::BitVector& cyclingCheck)
{
  TInt max = std::numeric_limits<TInt>::max();
  unsigned int sizeFace = 3;
  std::vector<TInt> nodes = SimplicesTriangle(this, tri).getNodes();
  Variable<int>* BND_CURVE_COLOR = getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_TRIANGLES   = getVariable<int,SimplicesTriangle>("BND_TRIANGLES");


  for(unsigned int local = 0 ; local < sizeFace ; local++)
  {
    TInt nodeA = nodes[local];
    TInt nodeB = nodes[(local + 1) % sizeFace];

    if((*BND_CURVE_COLOR)[nodeA] != 0 && (*BND_CURVE_COLOR)[nodeB] != 0)
    {
      if((*BND_CURVE_COLOR)[nodeA] != (*BND_CURVE_COLOR)[nodeB])
      {
        TInt nodeC_Local = (local + 2) % sizeFace;
        TInt adjTriangle = m_tri_adj[tri][nodeC_Local];
        if(cyclingCheck[adjTriangle] == 0)
        {
          cyclingCheck.assign(adjTriangle);
          if((*BND_TRIANGLES)[adjTriangle] != 0)
          {
            return (*BND_TRIANGLES)[adjTriangle];
          }
          else
          {
            TInt ii = findRemainTriangleIdx(adjTriangle, cyclingCheck);
            return ii;
          }
        }
      }
    }
  }
}
/******************************************************************************/
std::vector<TInt> SimplexMesh::deleteTetra(const TInt ATetraIndex)
{
  int errorId = std::numeric_limits<int>::min();
  std::vector<TInt> nodes;
  if(ATetraIndex < 0 || ATetraIndex > m_tet_ids.capacity())
  {
    /*TODO execption*/
    std::cout << "ATetraIndex < 0 || ATetraIndex > m_tet_ids.capacity()" << std::endl;
  }
  else
  {
      if(m_tet_ids[ATetraIndex] != 0)
      {
        nodes.resize(4);
        nodes[0] = m_tet_nodes[ATetraIndex][0];
        nodes[1] = m_tet_nodes[ATetraIndex][1];
        nodes[2] = m_tet_nodes[ATetraIndex][2];
        nodes[3] = m_tet_nodes[ATetraIndex][3];

        //Réatribution de m_base[node] pour tout node du tétra ATetraIndex
        for(auto const & node : nodes)
        {
          if(m_base[node] == ATetraIndex)
          {
            //CHange the direction of the base
            std::vector<TSimplexID> ball{};
            std::vector<TSimplexID> && ballTMP = SimplicesNode(this, node).ballOf();
            for(auto const simplexInBall : ballTMP)
            {
              if(simplexInBall >= 0)
              {
                ball.push_back(simplexInBall);
              }
            }
            if(ball.size() >= 2)
            {
              m_base[node] = (ball[0] == ATetraIndex)? ball[1]: ball[0];
            }
            else
            {
              m_base[node] = errorId;
            }
          }
        }

        std::vector<TSimplexID> adjSimplex;
        adjSimplex.resize(4);
        adjSimplex[0] = m_tet_adj[ATetraIndex][0];
        adjSimplex[1] = m_tet_adj[ATetraIndex][1];
        adjSimplex[2] = m_tet_adj[ATetraIndex][2];
        adjSimplex[3] = m_tet_adj[ATetraIndex][3];

        for(auto const & simplexAdj : adjSimplex)
        {
          if(simplexAdj != errorId )
          {
            if(simplexAdj >= 0)
            {
              if(m_tet_ids[simplexAdj] != 0)
              {
                std::vector<TInt>&&  vecNode = SimplicesCell(this, ATetraIndex).intersectionNodes(SimplicesCell(this, simplexAdj));

                if(vecNode.size() == 3)
                {
                  std::vector<TInt> nodes = SimplicesCell(this, simplexAdj).getOtherNodeInSimplex(vecNode);
                  if(nodes.size() == 1)
                  {
                    TInt localNode = SimplicesCell(this, simplexAdj).getLocalNode(SimplicesNode(this, nodes[0]).getGlobalNode());
                    m_tet_adj[simplexAdj][localNode] = errorId;
                  }
                }
                else
                {
                  /*todo Exception*/
                  std::cout << "vecNode.size() != 3 " << vecNode.size()  << std::endl;
                }
              }
            }
            else
            {
              bool deleteFromTriNodesFlag = (m_tri_nodes[-simplexAdj][3] == ATetraIndex)? true : false;
              bool deleteFromTriAdjFlag   = (m_tri_adj[-simplexAdj][3] == ATetraIndex)? true : false;

              if(deleteFromTriNodesFlag && !deleteFromTriAdjFlag)
              {
                m_tri_nodes[-simplexAdj][3] = errorId;
              }
              else if(deleteFromTriAdjFlag && !deleteFromTriNodesFlag)
              {
                m_tri_adj[-simplexAdj][3] = errorId;
              }
            }
          }
        }
        m_tet_ids.unselect(ATetraIndex);
      }
  }
  return std::move(nodes);
}
/*****************************************************************************/
SimplexMesh& SimplexMesh::buildRobustLayerMesh(const unsigned int nbrLayer)
{
  deleteNode(632);
  return *this;
  //INITIALISATION
  std::vector<TSimplexID>    tets;
  std::vector<TInt>          nodeslayer;
  unsigned int               sizeNodeInTet = 4;
  int errorId = std::numeric_limits<int>::min();
  for(TInt tetIter = 0; tetIter < m_tet_ids.capacity(); tetIter++)
  {
    if(m_tet_ids[tetIter] == 1)
    {
      if(SimplicesCell(this, tetIter).inBorder())
      {
        float lengthFromTheCenter = 1.5;
        unsigned int flagCpt = 0;
        std::vector<TInt> nodes;
        for(unsigned int localNode = 0; localNode < sizeNodeInTet; localNode++)
        {
          nodes.push_back(SimplicesCell(this, tetIter).getNode(localNode).getGlobalNode());
        }
        for(auto const & node : nodes)
        {
          gmds::math::Point point = m_coords[node];
          float nodeFromCenter = sqrt(point.X()*point.X()* + point.Y()*point.Y() + point.Z()*point.Z());
          if(nodeFromCenter <= lengthFromTheCenter)
          {
            flagCpt = flagCpt + 1;
          }
        }
        if(flagCpt == 4)
        {
          tets.push_back(tetIter);
        }
      }
    }
  }
  std::vector<std::vector<TInt>> nodesInLayer;
  for(auto const & tet : tets)
  {
    std::vector<TInt>&& nodes = SimplicesCell(this, tet).getNodes();
    for(auto node : nodes)
    {
      std::vector<TSimplexID>&& ball = SimplicesNode(this, node).ballOf(true);
      if(std::find(ball.begin(), ball.end(), errorId) == ball.end())
      {
        nodes.erase(std::find(nodes.begin(), nodes.end(), node), nodes.end());
      }
    }
    nodesInLayer.push_back(nodes);
  }

  std::set<TInt> setNodesLayer;
  for(auto const & nodes : nodesInLayer)
  {
    setNodesLayer.insert(nodes.begin(), nodes.end());
  }
  nodeslayer.resize(setNodesLayer.size());
  std::move(setNodesLayer.begin(), setNodesLayer.end(), nodeslayer.begin());
  //FIN INITIALISATION


  std::vector<TSimplexID> markedSimplex{};
  std::vector<TSimplexID> markedSimplexInPreviousLayer{};
  CriterionRAIS criterionRAIS(new VolumeCriterion());
  float sizeLayer = 1.0f;
  std::set<TSimplexID> testNodes;
  //gmds::Variable<math::Vector3d>* var = newVariable<math::Vector3d, SimplicesNode>("normal");
  for(unsigned int layer = 0; layer < nbrLayer; layer++)
  {
    std::vector<TInt> futurNodes{};
    for(auto const node : nodeslayer)
    {
      std::vector<std::vector<TInt>> orderedFaces;
      std::vector<TSimplexID>&& ball = SimplicesNode(this, node).ballOf();
      ball.erase(std::remove_if(ball.begin(), ball.end(),
      [&](TSimplexID simplex)
      {
          bool flag = false;
          if(std::find(markedSimplexInPreviousLayer.begin(), markedSimplexInPreviousLayer.end(), simplex) != markedSimplexInPreviousLayer.end())
          {
            flag = true;
          }
          return flag;
      }), ball.end());

      if(ball.size() != 0)
      {
        math::Vector3d normal(0.0f, 0.0f, 0.0f);
        for(auto const & tet : ball)
        {
          std::vector<TInt>&& nodes = SimplicesCell(this, tet).getNodes();
          intersection(nodes, nodeslayer);
          if(nodes.size() == 3)
          {
            math::Vector3d normalTmp = SimplicesCell(this, tet).normalOfFace(nodes); // return the interior normal of the cell
            normal += normalTmp;
          }
        }
        normal.normalize();
        //TODO check the tetra marked and fill the marked vector
        std::vector<TSimplexID>&& initialCavity = SimplicesNode(this, node).ballOf();
        gmds::math::Point pointToAdd = m_coords[node] + sizeLayer * gmds::math::Point(normal[0], normal[1], normal[2]);
        float lengthFromCenter = std::sqrt(pointToAdd[0] * pointToAdd[0] + pointToAdd[1] * pointToAdd[1] + pointToAdd[2] * pointToAdd[2]);
        const TInt nodeToInsert = addNode(pointToAdd);
        futurNodes.push_back(nodeToInsert);
        //PointInsertion pi(this, SimplicesNode(this, nodeToInsert), criterionRAIS, initialCavity, markedSimplex);
        std::vector<TSimplexID>&& ballCurrentNode = SimplicesNode(this, node).ballOf();
        std::move(ballCurrentNode.begin(), ballCurrentNode.end(), std::back_inserter(markedSimplex));

      }
    }
    nodeslayer.clear();
    std::copy(markedSimplex.begin(), markedSimplex.end(), std::back_inserter(markedSimplexInPreviousLayer));
    std::move(futurNodes.begin(), futurNodes.end(), std::back_inserter(nodeslayer));
  }

  bool flag = false;
  if(flag)
  {
    for(unsigned int tetIter = 0; tetIter < m_tet_ids.capacity(); tetIter++)
    {
      if(m_tet_ids[tetIter] != 0)
      {
        if(std::find(testNodes.begin(), testNodes.end(), tetIter) == testNodes.end())
        {
          deleteTetra(tetIter);
        }
      }
    }
  }
  else
  {
    for(unsigned int tetIter = 0; tetIter < m_tet_ids.capacity(); tetIter++)
    {
      if(m_tet_ids[tetIter] != 0)
      {
        if(std::find(markedSimplex.begin(), markedSimplex.end(), tetIter) == markedSimplex.end())
        {
          deleteTetra(tetIter);
        }
      }
    }
  }

}
/*****************************************************************************/
SimplexMesh& SimplexMesh::buildRobustLayerMeshOrderedNode(const unsigned int nbrLayer)
{
  //INITIALISATION
  std::vector<TSimplexID>    tets;
  std::vector<TInt>          nodeslayer;
  unsigned int               sizeNodeInTet = 4;
  int                        border        = std::numeric_limits<int>::min();
  for(TInt node = 0; node < m_node_ids.capacity(); node++)
  {
    if(m_node_ids[node] == 1)
    {
      std::vector<TSimplexID>&& ball = SimplicesNode(this, node).ballOf(true);
      if(std::find(ball.begin(), ball.end(), border) != ball.end())
      {
        nodeslayer.push_back(node);
      }
    }
  }

  std::vector<math::Vector3d> orderedNormal{};
  std::vector<math::Vector3d> UnorderedNormal{};
  std::vector<unsigned int>   UnordererdIndx{};
  std::vector<TSimplexID> markedSimplex{};
  std::vector<TSimplexID> markedSimplexInPreviousLayer{};
  CriterionRAIS criterionRAIS(new VolumeCriterion());
  float sizeLayer = 0.5f;
  std::set<TSimplexID> testNodes;
  struct DataNode
  {
    double val;
    int    node;
    math::Vector3d normal;
  };

  //gmds::Variable<math::Vector3d>* var = newVariable<math::Vector3d, SimplicesNode>("normal");
  for(unsigned int layer = 0; layer < nbrLayer; layer++)
  {
    std::vector<TInt> futurNodes{};
    std::vector<DataNode> dataNodes{};

    for(auto const node : nodeslayer)
    {
      std::vector<std::vector<TInt>> orderedFaces;
      std::vector<TSimplexID>&& ball = SimplicesNode(this, node).ballOf();
      ball.erase(std::remove_if(ball.begin(), ball.end(),
      [&](TSimplexID simplex)
      {
          bool flag = false;
          if(std::find(markedSimplexInPreviousLayer.begin(), markedSimplexInPreviousLayer.end(), simplex) != markedSimplexInPreviousLayer.end())
          {
            flag = true;
          }
          return flag;
      }), ball.end());

      if(ball.size() != 0)
      {
        std::vector<Vector3d> normalFaces{};
        math::Vector3d normal(0.0, 0.0, 0.0);
        std::vector<math::Vector3d> normals;
        for(auto const & tet : ball)
        {
          std::vector<TInt>&& nodes = SimplicesCell(this, tet).getNodes();
          intersection(nodes, nodeslayer);
          std::vector<std::vector<TInt>> FacesNodes{};
          if(std::find(nodes.begin(), nodes.end(), node) != nodes.end() )
          {
            if(nodes.size() == 4 )
            {
              nodes.erase(std::remove(nodes.begin(), nodes.end(), node), nodes.end());
              std::vector<TInt> faceA{node, nodes[0], nodes[1]};
              std::vector<TInt> faceB{node, nodes[1], nodes[2]};
              std::vector<TInt> faceC{node, nodes[0], nodes[2]};
              FacesNodes.push_back(std::move(faceA));
              FacesNodes.push_back(std::move(faceB));
              FacesNodes.push_back(std::move(faceC));

            }
            else if(nodes.size() == 3)
            {
              FacesNodes.push_back(std::move(nodes));
            }

            if(FacesNodes.size() != 0)
            {
              std::vector<math::Vector3d>&& normalFaces = SimplicesCell(this, tet).normalOfFaces(FacesNodes);
              for(auto & normalFace : normalFaces)
              {
                  normalFace.normalize();
                  normal += normalFace;
              }
              normal.normalize();
              std::move(normalFaces.begin(), normalFaces.end(), std::back_inserter(normals));
            }
          }
        }

        std::vector<math::Vector3d> vectorToDelete{};
        double epsilonNeg = -0.998;
        double epsilonPos =  0.998;
        normals.erase(std::remove_if(normals.begin(), normals.end(),
        [&](math::Vector3d& vec)
        {
          bool flag = false;
          for(auto const & vecComp : normals)
          {
            double dotComp = vecComp.dot(vec);
            if(dotComp <= epsilonNeg)
            {
                vectorToDelete.push_back(vecComp);
                flag = true;
            }
          }
          return flag;
        }), normals.end());


        normals.erase(std::remove_if(normals.begin(), normals.end(),
        [&](math::Vector3d& vec)
        {
          bool flag = false;
          for(auto const & vecComp : vectorToDelete)
          {
            double dotComp = vecComp.dot(vec);
            if(dotComp >= epsilonPos)
            {
                vectorToDelete.push_back(vecComp);
                flag = true;
            }
          }
          return flag;
        }) , normals.end());

        Eigen::MatrixXd matNormalByFace(3, normals.size());
        Eigen::MatrixXd oneMatrix       = Eigen::MatrixXd::Ones(normals.size(), normals.size());
        for(unsigned int column = 0; column < normals.size(); column++)
        {
          matNormalByFace(0,column) = normals[column].X();
          matNormalByFace(1,column) = normals[column].Y();
          matNormalByFace(2,column) = normals[column].Z();
        }

        Eigen::MatrixXd covMat = (matNormalByFace.transpose() * matNormalByFace + oneMatrix)* 0.5; //coincidence si les coeff sont > 0 normalement il faut normaliser 0.5*(covMat + 1)
        /*std::cout << "covMat" << std::endl;
        std::cout << covMat << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;*/


        double val       = 1.0;
        double exponent  = 2.0;
        for(unsigned int i = 0; i < normals.size() - 1; i++)
        {
          for(unsigned int j = 1; j < normals.size(); j++)
          {
            val *= covMat(i,j);//std::pow(covMat(i,j), exponent);
          }
        }

        DataNode dataNode;
        dataNode.val    = val;
        dataNode.node   = node;
        dataNode.normal = normal;

        dataNodes.push_back(dataNode);

      }
    }

    std::sort(dataNodes.begin(), dataNodes.end(),
    [&](const DataNode& dataNode0, const DataNode& dataNode1)
    {
      return dataNode0.val < dataNode1.val;
    });

    for(auto const & dataNode : dataNodes)
    {
      std::vector<TSimplexID>&& initialCavity = SimplicesNode(this, dataNode.node).ballOf();
      initialCavity.erase(std::remove_if(initialCavity.begin(), initialCavity.end(),
      [&](TSimplexID simplex)
      {
          bool flag = false;
          if(std::find(markedSimplexInPreviousLayer.begin(), markedSimplexInPreviousLayer.end(), simplex) != markedSimplexInPreviousLayer.end())
          {
            flag = true;
          }
          return flag;
      }), initialCavity.end());
      gmds::math::Point pointToAdd = m_coords[dataNode.node] + sizeLayer * gmds::math::Point(dataNode.normal[0], dataNode.normal[1], dataNode.normal[2]);
      float lengthFromCenter = std::sqrt(pointToAdd[0] * pointToAdd[0] + pointToAdd[1] * pointToAdd[1] + pointToAdd[2] * pointToAdd[2]);
      const TInt nodeToInsert = addNode(pointToAdd);
      futurNodes.push_back(nodeToInsert);
      //PointInsertion pi(this, SimplicesNode(this, nodeToInsert), criterionRAIS, initialCavity, markedSimplex);
      std::vector<TSimplexID>&& ballCurrentNode = SimplicesNode(this, dataNode.node).ballOf();
      std::move(ballCurrentNode.begin(), ballCurrentNode.end(), std::back_inserter(markedSimplex));
    }


    nodeslayer.clear();
    std::copy(markedSimplex.begin(), markedSimplex.end(), std::back_inserter(markedSimplexInPreviousLayer));
    std::move(futurNodes.begin(), futurNodes.end(), std::back_inserter(nodeslayer));


    int i = 0;
    for(auto const dataNode : dataNodes)
    {
      if(dataNode.val == 1)
      {
        std::cout << "dataNode.node : " << dataNode.node << std::endl;
        deleteNode(dataNode.node);
        if(i == 200)
        {
          break;
        }
      }
      i++;
    }
  }

  return *this;
}
/*****************************************************************************/
SimplexMesh& SimplexMesh::buildRobustLayerMeshOrderedNode01(const unsigned int nbrLayer)
{
  //INITIALISATION
  std::vector<TSimplexID>    tets;
  std::vector<TInt>          nodeslayer;
  unsigned int               sizeNodeInTet = 4;
  int                        border        = std::numeric_limits<int>::min();
  for(TInt node = 0; node < m_node_ids.capacity(); node++)
  {
    if(m_node_ids[node] == 1)
    {
      std::vector<TSimplexID>&& ball = SimplicesNode(this, node).ballOf(true);
      if(std::find(ball.begin(), ball.end(), border) != ball.end())
      {
        nodeslayer.push_back(node);
      }
    }
  }

  std::vector<math::Vector3d> orderedNormal{};
  std::vector<math::Vector3d> UnorderedNormal{};
  std::vector<unsigned int>   UnordererdIndx{};
  std::vector<TSimplexID> markedSimplex{};
  std::vector<TSimplexID> markedSimplexInPreviousLayer{};
  CriterionRAIS criterionRAIS(new VolumeCriterion());
  float sizeLayer = 0.5f;
  std::set<TSimplexID> testNodes;
  struct DataNode
  {
    double meanCourvature;
    double gaussianCourvature;
    int    node;
    math::Vector3d normal;
  };

  gmds::Variable<double>* varMeanCurvature     = newVariable<double, SimplicesNode>("meanCurvature");
  gmds::Variable<double>* varGaussianCurvature = newVariable<double, SimplicesNode>("gaussianCurvature");
  for(unsigned int layer = 0; layer < nbrLayer; layer++)
  {
    std::vector<TInt> futurNodes{};
    std::vector<DataNode> dataNodes{};

    for(auto const node : nodeslayer)
    {
      std::vector<TSimplexID>&& ball = SimplicesNode(this, node).ballOf();
      ball.erase(std::remove_if(ball.begin(), ball.end(),
      [&](TSimplexID simplex)
      {
          bool flag = false;
          if(std::find(markedSimplexInPreviousLayer.begin(), markedSimplexInPreviousLayer.end(), simplex) != markedSimplexInPreviousLayer.end())
          {
            flag = true;
          }
          return flag;
      }), ball.end());

      if(ball.size() != 0)
      {
        std::vector<Vector3d> normalFaces{};
        math::Vector3d normal(0.0, 0.0, 0.0);
        std::vector<math::Vector3d> normals;
        for(auto const & tet : ball)
        {
          std::vector<std::vector<TInt>> FacesNodes{};
          std::vector<TInt>&& nodes = SimplicesCell(this, tet).getNodes();
          intersection(nodes, nodeslayer);
          if(std::find(nodes.begin(), nodes.end(), node) != nodes.end() )
          {
            if(nodes.size() == 4 )
            {
              nodes.erase(std::remove(nodes.begin(), nodes.end(), node), nodes.end());
              std::vector<TInt> faceA{node, nodes[0], nodes[1]};
              std::vector<TInt> faceB{node, nodes[1], nodes[2]};
              std::vector<TInt> faceC{node, nodes[0], nodes[2]};
              FacesNodes.push_back(std::move(faceA));
              FacesNodes.push_back(std::move(faceB));
              FacesNodes.push_back(std::move(faceC));
            }
            else if(nodes.size() == 3)
            {
              FacesNodes.push_back(std::move(nodes));
            }

            if(FacesNodes.size() != 0)
            {
              std::vector<math::Vector3d>&& normalFaces = SimplicesCell(this, tet).normalOfFaces(FacesNodes);
              for(auto & normalFace : normalFaces)
              {
                  normalFace.normalize();
                  normal += normalFace;
              }
            }
          }
        }

        //Construction de la base local au node courant pour construire le path de Monge..
        Eigen::Vector3d zLocal = Eigen::Vector3d(normal.X(), normal.Y(), normal.Z());
        Eigen::Vector3d uLocal = zLocal.unitOrthogonal();
        Eigen::Vector3d vLocal = zLocal.cross(uLocal);
        zLocal.normalize();
        uLocal.normalize();
        vLocal.normalize();

        //Matrice de la base local;
        Eigen::Matrix3d P;
        P.col(2) = zLocal;
        P.col(0) = vLocal;
        P.col(1) = uLocal;

        std::vector<TInt>&& neihborsOfNode = SimplicesNode(this, node).commonNodeInBall(nodeslayer);
        std::vector<Eigen::Vector3d> nodesInLocalBase{};
        //calculer les coordonnées des vosins dans la base local
        for(auto const neighbor : neihborsOfNode)
        {
            if(neighbor != node)
            {
              math::Point recenterNeighborPoint = m_coords[neighbor] - m_coords[node];
              Eigen::Vector3d recenterVec       = Eigen::Vector3d(recenterNeighborPoint.X(), recenterNeighborPoint.Y(), recenterNeighborPoint.Z());
              Eigen::Vector3d localNode         = P.inverse() * recenterVec;
              nodesInLocalBase.push_back(localNode);
            }
        }

        //Moindre carré pour déterminer les coefficient (a, b, c) du patch de monge
        Eigen::MatrixXf A = Eigen::MatrixXf(nodesInLocalBase.size(),3);
        Eigen::VectorXf B = Eigen::VectorXf(nodesInLocalBase.size());
        unsigned int rowCpt = 0;
        for(auto const & nodeInBaseLocal : nodesInLocalBase)
        {
          B.row(rowCpt) << nodeInBaseLocal[2];
          A.row(rowCpt) << nodeInBaseLocal[0]*nodeInBaseLocal[0]/2,nodeInBaseLocal[1]*nodeInBaseLocal[1]/2,nodeInBaseLocal[0]*nodeInBaseLocal[1];
          rowCpt++;
        }

        //Shape operator
        Eigen::Vector3f  a_tmp = (A.transpose() * A).ldlt().solve(A.transpose() * B);
        float k1 = a_tmp[0];
        float k2 = a_tmp[1];
        double gaussianCurvature = (double)(k1*k2);
        double meanCurvature = (double)((k1+k2) /2.0);

        varMeanCurvature->set(node, meanCurvature);
        varGaussianCurvature->set(node, gaussianCurvature);

      }
    }
/*
    std::sort(dataNodes.begin(), dataNodes.end(),
    [&](const DataNode& dataNode0, const DataNode& dataNode1)
    {
      return dataNode0.val > dataNode1.val;
    });

    for(auto const & dataNode : dataNodes)
    {
      std::vector<TSimplexID>&& initialCavity = SimplicesNode(this, dataNode.node).ballOf();
      initialCavity.erase(std::remove_if(initialCavity.begin(), initialCavity.end(),
      [&](TSimplexID simplex)
      {
          bool flag = false;
          if(std::find(markedSimplexInPreviousLayer.begin(), markedSimplexInPreviousLayer.end(), simplex) != markedSimplexInPreviousLayer.end())
          {
            flag = true;
          }
          return flag;
      }), initialCavity.end());
      gmds::math::Point pointToAdd = m_coords[dataNode.node] + sizeLayer * gmds::math::Point(dataNode.normal[0], dataNode.normal[1], dataNode.normal[2]);
      float lengthFromCenter = std::sqrt(pointToAdd[0] * pointToAdd[0] + pointToAdd[1] * pointToAdd[1] + pointToAdd[2] * pointToAdd[2]);
      const TInt nodeToInsert = addNode(pointToAdd);
      futurNodes.push_back(nodeToInsert);
      PointInsertion pi(this, SimplicesNode(this, nodeToInsert), criterionRAIS, initialCavity, markedSimplex);
      std::vector<TSimplexID>&& ballCurrentNode = SimplicesNode(this, dataNode.node).ballOf();
      std::move(ballCurrentNode.begin(), ballCurrentNode.end(), std::back_inserter(markedSimplex));
    }


    nodeslayer.clear();
    std::copy(markedSimplex.begin(), markedSimplex.end(), std::back_inserter(markedSimplexInPreviousLayer));
    std::move(futurNodes.begin(), futurNodes.end(), std::back_inserter(nodeslayer));


    int i = 0;
    for(auto const dataNode : dataNodes)
    {
      std::cout << std::endl;
      std::cout << std::endl;
      if(dataNode.val == 1)
      {
        std::cout << "dataNode.node : " << dataNode.node << std::endl;
        deleteNode(dataNode.node);
        if(i == 200)
        {
          break;
        }
      }
      i++;
    }*/
  }

  return *this;
}
/*****************************************************************************/
std::vector<TInt> SimplexMesh::buildQuadFaceFromCoordPt0(const std::vector<math::Point>& ptCoords, std::vector<TSimplexID>& markedSimplex)
{
  std::vector<TInt> nodes{};
  if(ptCoords.size() == 4)
  {
    //BUILDING THE TWO BASE fACE OF THE FUTUR QUAND COMPOSED BY TET
    //the ptCoord vector is ordered..
    int errorId = std::numeric_limits<int>::min();
    CriterionRAIS criterionRAIS(new VolumeCriterion());

    for(auto const & ptCoord : ptCoords)
    {
      TInt nodeId = addNode(ptCoord);
      nodes.push_back(nodeId);
      std::vector<TSimplexID>&& ball = SimplicesNode(this, nodeId).ballOf();
      if(ball.size() == 0) // check si le node existait deja et était déjà inséré
      {
        //PointInsertion(this, SimplicesNode(this, nodeId), criterionRAIS);
      }
    }



    unsigned int nodesSizeBase = 4;
    float epsilon = 10E-5f;
    float weight = 1.0f/8.0f;
    math::Point centerQuad = weight * (ptCoords[0] + ptCoords[1] + ptCoords[2] + ptCoords[3] + ptCoords[4] + ptCoords[5] + ptCoords[6] + ptCoords[7]);
    weight = 0.25f;
    math::Point centerFace = weight *  (ptCoords[0] + ptCoords[1] + ptCoords[2] + ptCoords[3]);
    weight = 0.5f;
    std::vector<TSimplexID> FaceCavity{};

    bool flag = false;
    for(unsigned int node = 0; node < nodesSizeBase; node++)
    {
      //initialCavity is the only tet in the direction of the extremity of the current edge
      std::vector<TSimplexID>&& tets = SimplicesNode(this, nodes[node]).ballOf();
      if(tets.size() != 0)
      {
        const math::Point& nextNode  = SimplicesNode(this, nodes[(node + 1) % nodesSizeBase]).getCoords();
        math::Vector3d normal =  (SimplicesNode(this, nodes[(node + 1) % nodesSizeBase]).getCoords() - SimplicesNode(this, nodes[node]).getCoords());
        normal.normalize();
        flag = false;
        while(flag != true)
        {
          epsilon = epsilon * 0.1f;
          const math::Point directNeighbor = SimplicesNode(this, nodes[node]).getCoords() + epsilon * normal;
          TSimplexID tetraContainingDirectNeighbor;
          for(auto & tet : tets)
          {
            if(SimplicesCell(this,tet).isPointInSimplex(directNeighbor))
            {
              flag = true;
              //PointInsertion
              std::vector<TSimplexID> initialCavity{tet};
              //PointInsertion pi(this, SimplicesNode(this, nodes[(node + 1) % nodesSizeBase]), criterionRAIS, initialCavity, markedSimplex);
            }
          }
        }
      }
      else
      {
        //TODO
        //exception
        std::cout << SimplicesNode(this, nodes[node]) << "has no tet in his ball" << std::endl;
        return std::move(nodes);
      }
    }



    //initialisation de la cavité pour construire une face "parfaite"...
    for(unsigned int node = 0; node < nodesSizeBase; node++)
    {
      math::Point centerEdge = weight * (ptCoords[(node + 1) % nodesSizeBase] + ptCoords[node]);
      math::Vector3d normalEdge =  centerFace - centerEdge;
      normalEdge.normalize();
      flag = false;
      while(flag != true)
      {
        epsilon = epsilon * 0.1f;
        math::Point NodeNeighborEdge = centerEdge + epsilon * normalEdge;
        std::vector<TSimplexID> && shell = SimplicesNode(this, nodes[(node + 1) % nodesSizeBase]).shell(SimplicesNode(this, nodes[node]));
        for(auto & tet : shell)
        {
          if(SimplicesCell(this, tet).isPointInSimplex(NodeNeighborEdge))
          {
            flag = true;
            FaceCavity.push_back(tet);
            break;
          }
        }
      }
    }

    //on retire les tetra marqué dans les quad générés précedemnent
    FaceCavity.erase(std::remove_if(FaceCavity.begin(), FaceCavity.end(),[=](TSimplexID tetInInitialCavity)
    {
      bool flag = false;
      if(std::find(markedSimplex.begin(), markedSimplex.end(), tetInInitialCavity) != markedSimplex.end())
      {
        flag = true;
      }
      return flag;
    }),FaceCavity.end());

    //PointInsertion pi(this, SimplicesNode(this, nodes[0]), criterionRAIS, FaceCavity);
    //The base face (0, 1, 2, 3) of the quad is now build

  }
  else
  {
    //TODO exception
    std::cout << "ptCoords size != 4" << std::endl;
    return std::move(nodes);
  }


  return std::move(nodes);
}
/*****************************************************************************/
void SimplexMesh::buildQuadFaceFromNode(const std::vector<TInt>& nodesFace, std::vector<TSimplexID>& markedSimplex)
{
  if(nodesFace.size() == 4)
  {
    if(SimplicesNode(this, nodesFace[0]).shell(SimplicesNode(this, nodesFace[1])).size() > 0 &&
        SimplicesNode(this, nodesFace[1]).shell(SimplicesNode(this, nodesFace[2])).size() > 0 &&
          SimplicesNode(this, nodesFace[2]).shell(SimplicesNode(this, nodesFace[3])).size() > 0 &&
            SimplicesNode(this, nodesFace[3]).shell(SimplicesNode(this, nodesFace[0])).size() > 0 )

            {
              std::vector<TSimplexID>&& ballNode0 = SimplicesNode(this, nodesFace[0]).ballOf();
              std::vector<TSimplexID>&& ballNode1 = SimplicesNode(this, nodesFace[1]).ballOf();
              std::vector<TSimplexID>&& ballNode2 = SimplicesNode(this, nodesFace[2]).ballOf();
              std::vector<TSimplexID>&& ballNode3 = SimplicesNode(this, nodesFace[3]).ballOf();

              std::vector<std::vector<TSimplexID>> ball012{ballNode0, ballNode1, ballNode2};
              std::vector<std::vector<TSimplexID>> ball032{ballNode0, ballNode3, ballNode2};

              std::vector<std::vector<TSimplexID>> ball013{ballNode0, ballNode1, ballNode3};
              std::vector<std::vector<TSimplexID>> ball123{ballNode1, ballNode2, ballNode3};


              std::vector<TSimplexID>&& intersectionBall012 = intersectionSimplex(ball012);
              std::vector<TSimplexID>&& intersectionBall032 = intersectionSimplex(ball032);

              std::vector<TSimplexID>&& intersectionBall013 = intersectionSimplex(ball013);
              std::vector<TSimplexID>&& intersectionBall123 = intersectionSimplex(ball123);

              if((intersectionBall012.size() == 0 && intersectionBall032.size() == 0)
                  || (intersectionBall013.size() == 0 && intersectionBall123.size() == 0))

                  {

                    //BUILDING THE TWO BASE fACE OF THE FUTUR QUAND COMPOSED BY TET
                    //the ptCoord vector is ordered..
                    int errorId = std::numeric_limits<int>::min();
                    CriterionRAIS criterionRAIS(new VolumeCriterion());

                    unsigned int nodesSizeBase = 4;
                    float epsilon = 10E-5f;
                    std::vector<math::Point> ptCoords{
                      SimplicesNode(this, nodesFace[0]).getCoords(),
                      SimplicesNode(this, nodesFace[1]).getCoords(),
                      SimplicesNode(this, nodesFace[2]).getCoords(),
                      SimplicesNode(this, nodesFace[3]).getCoords()
                    };

                    float weight  = 0.25f;
                    math::Point centerFace = weight *  (ptCoords[0] + ptCoords[1] + ptCoords[2] + ptCoords[3]);
                    std::vector<TSimplexID> FaceCavity{};

                    bool flag = false;
                    weight = 0.5f;
                    //initialisation de la cavité pour construire une face "parfaite"...
                    for(unsigned int node = 0; node < nodesSizeBase; node++)
                    {
                      math::Point centerEdge = weight * (ptCoords[(node + 1) % nodesSizeBase] + ptCoords[node]);
                      math::Vector3d normalEdge =  centerFace - centerEdge;
                      normalEdge.normalize();
                      flag = false;
                      epsilon = 10E-5f;
                      std::vector<TSimplexID> && shell = SimplicesNode(this, nodesFace[(node + 1) % nodesSizeBase]).shell(SimplicesNode(this, nodesFace[node]));

                      while(flag != true)
                      {
                        math::Point NodeNeighborEdge = centerEdge + epsilon * normalEdge;
                        for(auto & tet : shell)
                        {
                          if(SimplicesCell(this,tet).isPointInSimplex(NodeNeighborEdge))
                          {
                            flag = true;
                            FaceCavity.push_back(tet);
                            break;
                          }
                        }
                        epsilon += 10E-5;
                      }
                    }

                    //on retire les tetra marqué dans les quad générés précedemnent
                    FaceCavity.erase(std::remove_if(FaceCavity.begin(), FaceCavity.end(),[=](TSimplexID tetInInitialCavity)
                    {
                      bool flag = false;
                      if(std::find(markedSimplex.begin(), markedSimplex.end(), tetInInitialCavity) != markedSimplex.end())
                      {
                        flag = true;
                      }
                      return flag;
                    }),FaceCavity.end());

                    //PointInsertion pi(this, SimplicesNode(this, nodesFace[0]), criterionRAIS, FaceCavity);
                    //The base face (0, 1, 2, 3) of the quad is now build
                    }
                  }
    }
    else
    {
      //TODO exception
      std::cout << "ptCoords size != 4" << std::endl;
    }
}
/*****************************************************************************/
std::vector<TInt> SimplexMesh::buildQuadFaceFromCoordPt(const std::vector<math::Point>& ptCoords, std::vector<TSimplexID>& markedSimplex)
{
  std::vector<TInt> nodes{};
  if(ptCoords.size() == 4)
  {
    //BUILDING THE TWO BASE fACE OF THE FUTUR QUAND COMPOSED BY TET
    //the ptCoord vector is ordered..
    int errorId = std::numeric_limits<int>::min();
    CriterionRAIS criterionRAIS(new VolumeCriterion());

    for(auto const & ptCoord : ptCoords)
    {
      TInt nodeId = addNode(ptCoord);
      nodes.push_back(nodeId);
      std::vector<TSimplexID>&& ball = SimplicesNode(this, nodeId).ballOf();
      if(ball.size() == 0) // check si le node existait deja et était déjà inséré
      {
        //PointInsertion(this, SimplicesNode(this, nodeId), criterionRAIS);
      }
    }

    unsigned int nodesSizeBase = 4;

    bool flag = false;
    for(unsigned int node = 0; node < nodesSizeBase; node++)
    {
      EdgeInsertion(this, SimplicesNode(this, nodes[node]), SimplicesNode(this, nodes[(node + 1) % nodesSizeBase]), criterionRAIS, markedSimplex);
    }
    buildQuadFaceFromNode(nodes, markedSimplex);
  }
  else
  {
    //TODO exception
    std::cout << "ptCoords size != 4" << std::endl;
    return std::move(nodes);
  }


  return std::move(nodes);
}
/*****************************************************************************/
SimplexMesh& SimplexMesh::buildPrismFromFace0(const std::vector<TInt> & nodesDown, const std::vector<math::Point> & nodesUpCoord, std::vector<TSimplexID>& markedSimplex, const int debug)
{

  if(nodesDown.size() == 3 && nodesUpCoord.size() == 3 )
  {
    const std::vector<TInt>& orientedNodeDown = nodesDown;
    const std::vector<math::Point>& orientedNodeUp   = nodesUpCoord;

    std::vector<TInt> nodesUp{};
    //now : build the prims of the faces build before
    //build the first tet of the prism {0, 1, 2, 4}
    std::vector<TSimplexID>&& ballNode0 = SimplicesNode(this, orientedNodeDown[0]).ballOf();
    std::vector<TSimplexID>&& ballNode1 = SimplicesNode(this, orientedNodeDown[1]).ballOf();
    std::vector<TSimplexID>&& ballNode2 = SimplicesNode(this, orientedNodeDown[2]).ballOf();

    std::vector<std::vector<TSimplexID>> ball012{};

    ball012.push_back(ballNode0);
    ball012.push_back(ballNode1);
    ball012.push_back(ballNode2);


    std::vector<TSimplexID>&& intersectionBall012 = intersectionSimplex(ball012);
    std::vector<TInt> node012{orientedNodeDown[0], orientedNodeDown[1], orientedNodeDown[2]};
    TInt errorId = std::numeric_limits<int>::min();
    TSimplexID startingSimplex = errorId;
    math::Vector3d vec = orientedNodeUp[0] - SimplicesNode(this, orientedNodeDown[0]).getCoords();
    vec.normalize();


    for(auto const & tet : intersectionBall012)
    {
      std::vector<TInt> nodesInSimplex_NOT012 = SimplicesCell(this, tet).getOtherNodeInSimplex(node012);
      if(nodesInSimplex_NOT012.size() == 1)
      {
        math::Vector3d vecComp   = (SimplicesNode(this, nodesInSimplex_NOT012[0]).getCoords() - SimplicesNode(this, orientedNodeDown[0]).getCoords());
        float dotComp = vec.dot(vecComp);
        if(dotComp > 0.0f)
        {
          startingSimplex = tet;
          break;
        }
      }
      else
      {
        //TODO
        std::cout << "nodesInSimplex.size() != 1" << std::endl;
        return *this;
      }
    }


    if(startingSimplex != errorId)
    {

        //create the first tetra {0, 1, 2 , 4}
        std::vector<TSimplexID> firstInitCavity{startingSimplex};
        CriterionRAIS criterionRAIS(new VolumeCriterion());

        TInt nodeId = addNode(orientedNodeUp[0]);
        nodesUp.push_back(nodeId);

        if(debug == 3 && startingSimplex != errorId)
        {
          //std::cout << SimplicesNode(this, nodeId) << std::endl;
          return *this;
        }
        //PointInsertion(this, SimplicesNode(this, nodesUp.back()), criterionRAIS, firstInitCavity, markedSimplex);

        //find the first marked tetra
        TInt firstNodeDownLocal = 0;
        TInt node0             = orientedNodeDown[0];
        TInt node1             = orientedNodeDown[1];
        TInt node2             = orientedNodeDown[2];
        TInt node4             = nodesUp.back();

        std::vector<TSimplexID>&& ballNode0       = SimplicesNode(this, node0).ballOf();
        std::vector<TSimplexID>&& ballNode1       = SimplicesNode(this, node1).ballOf();
        std::vector<TSimplexID>&& ballNode2       = SimplicesNode(this, node2).ballOf();
        std::vector<TSimplexID>&& ballNode4       = SimplicesNode(this, node4).ballOf();

        std::vector<std::vector<TSimplexID>> ballNodes;

        ballNodes.push_back(ballNode0);
        ballNodes.push_back(ballNode1);
        ballNodes.push_back(ballNode2);
        ballNodes.push_back(ballNode4);

        std::vector<TSimplexID> inter = intersectionSimplex(ballNodes);
        if(inter.size() == 1)
        {
          markedSimplex.push_back(inter.front());
          //first tet marked is inter[0]
          //now iterate over the first tetra build by chosing the tetra oposite to the exNode from where it was build
          std::vector<TInt> nodesDownLocalToUseForBuildThePrismLocal{1, 2}; // {1, 2, 2, 3, 0}
          for(const auto & nodeDownLocal : nodesDownLocalToUseForBuildThePrismLocal)
          {
            firstInitCavity.clear();
            TInt nodeDownPreviousLocal = nodeDownLocal - 1;
            TInt nodeDownPrevious      = orientedNodeDown[nodeDownPreviousLocal];
            TInt nodeUp                = addNode(orientedNodeUp[nodeDownLocal]);

            nodesUp.push_back(nodeUp);
            TSimplexID tetraOpp   = getOppositeCell(nodeDownPrevious, markedSimplex.back());
            firstInitCavity.push_back(tetraOpp);

            //PointInsertion(this, SimplicesNode(this, nodesUp.back()), criterionRAIS, firstInitCavity, markedSimplex);

            tetraOpp = getOppositeCell(nodeDownPrevious, markedSimplex.back());
            markedSimplex.push_back(tetraOpp);
          }
        }
        else
        {
          std::cout << "firstMarkedSimplex.size() != 1" << std::endl;
          return *this;
        }
    }
    else
    {
      //TODO exception
      std::cout << "starting Simplex == errorId" << std::endl;
    }
  }
  else
  {
    //TODO
    std::cout << "nodesDownAndNodeUp[0].size() == 3 && nodesDownAndNodeUp[1].size() == 3 " << std::endl;
    return *this;
  }

  return *this;
}
/****************************************************************************/
SimplexMesh& SimplexMesh::buildPrismFromFace(const std::vector<TInt> & nodesDown, const std::vector<math::Point> & nodesUpCoord, std::vector<TSimplexID>& markedSimplex, const int debug)
{

  if(nodesDown.size() == 3 && nodesUpCoord.size() == 3 )
  {
    const std::vector<TInt>& orientedNodeDown = nodesDown;
    const std::vector<math::Point>& orientedNodeUp   = nodesUpCoord;

    std::vector<TInt> nodesUp{};
    //now : build the prims of the faces build before
    //build the first tet of the prism {0, 1, 2, 4}
    std::vector<TSimplexID>&& ballNode0 = SimplicesNode(this, orientedNodeDown[0]).ballOf();
    std::vector<TSimplexID>&& ballNode1 = SimplicesNode(this, orientedNodeDown[1]).ballOf();
    std::vector<TSimplexID>&& ballNode2 = SimplicesNode(this, orientedNodeDown[2]).ballOf();

    std::vector<std::vector<TSimplexID>> ball012{};

    ball012.push_back(ballNode0);
    ball012.push_back(ballNode1);
    ball012.push_back(ballNode2);


    std::vector<TSimplexID>&& intersectionBall012 = intersectionSimplex(ball012);
    std::vector<TInt> node012{orientedNodeDown[0], orientedNodeDown[1], orientedNodeDown[2]};
    TInt errorId = std::numeric_limits<int>::min();
    TSimplexID startingSimplex = errorId;
    math::Vector3d vec = orientedNodeUp[0] - SimplicesNode(this, orientedNodeDown[0]).getCoords();
    vec.normalize();


    for(auto const & tet : intersectionBall012)
    {
      std::vector<TInt> nodesInSimplex_NOT012 = SimplicesCell(this, tet).getOtherNodeInSimplex(node012);
      if(nodesInSimplex_NOT012.size() == 1)
      {
        math::Vector3d vecComp   = (SimplicesNode(this, nodesInSimplex_NOT012[0]).getCoords() - SimplicesNode(this, orientedNodeDown[0]).getCoords());
        float dotComp = vec.dot(vecComp);
        if(dotComp > 0.0f)
        {
          startingSimplex = tet;
          break;
        }
      }
      else
      {
        //TODO
        std::cout << "nodesInSimplex.size() != 1" << std::endl;
        return *this;
      }
    }


    if(startingSimplex != errorId)
    {

        //create the first tetra {0, 1, 2 , 4}
        std::vector<TSimplexID> firstInitCavity{startingSimplex};
        CriterionRAIS criterionRAIS(new VolumeCriterion());

        TInt nodeId = addNode(orientedNodeUp[0]);
        nodesUp.push_back(nodeId);

        if(debug == 3 && startingSimplex != errorId)
        {
          //std::cout << SimplicesNode(this, nodeId) << std::endl;
          return *this;
        }
        //PointInsertion(this, SimplicesNode(this, nodesUp.back()), criterionRAIS, firstInitCavity, markedSimplex);

        //find the first marked tetra
        TInt firstNodeDownLocal = 0;
        TInt node0             = orientedNodeDown[0];
        TInt node1             = orientedNodeDown[1];
        TInt node2             = orientedNodeDown[2];
        TInt node4             = nodesUp.back();

        std::vector<TSimplexID>&& ballNode0       = SimplicesNode(this, node0).ballOf();
        std::vector<TSimplexID>&& ballNode1       = SimplicesNode(this, node1).ballOf();
        std::vector<TSimplexID>&& ballNode2       = SimplicesNode(this, node2).ballOf();
        std::vector<TSimplexID>&& ballNode4       = SimplicesNode(this, node4).ballOf();

        std::vector<std::vector<TSimplexID>> ballNodes;

        ballNodes.push_back(ballNode0);
        ballNodes.push_back(ballNode1);
        ballNodes.push_back(ballNode2);
        ballNodes.push_back(ballNode4);

        std::vector<TSimplexID> inter = intersectionSimplex(ballNodes);
        if(inter.size() == 1)
        {
          markedSimplex.push_back(inter.front());
          //first tet marked is inter[0]
          //now iterate over the first tetra build by chosing the tetra oposite to the exNode from where it was build
          std::vector<TInt> nodesDownLocalToUseForBuildThePrismLocal{1, 2}; // {1, 2, 2, 3, 0}
          for(const auto & nodeDownLocal : nodesDownLocalToUseForBuildThePrismLocal)
          {
            firstInitCavity.clear();
            TInt nodeDownPreviousLocal = nodeDownLocal - 1;
            TInt nodeDownPrevious      = orientedNodeDown[nodeDownPreviousLocal];
            TInt nodeUp                = addNode(orientedNodeUp[nodeDownLocal]);

            nodesUp.push_back(nodeUp);
            TSimplexID tetraOpp   = getOppositeCell(nodeDownPrevious, markedSimplex.back());
            firstInitCavity.push_back(tetraOpp);

            //PointInsertion(this, SimplicesNode(this, nodesUp.back()), criterionRAIS, firstInitCavity, markedSimplex);

            tetraOpp = getOppositeCell(nodeDownPrevious, markedSimplex.back());
            markedSimplex.push_back(tetraOpp);
          }
        }
        else
        {
          std::cout << "firstMarkedSimplex.size() != 1" << std::endl;
          return *this;
        }
    }
    else
    {
      //TODO exception
      std::cout << "starting Simplex == errorId" << std::endl;
    }
  }
  else
  {
    //TODO
    std::cout << "nodesDownAndNodeUp[0].size() == 3 && nodesDownAndNodeUp[1].size() == 3 " << std::endl;
    return *this;
  }

  return *this;
}
/*****************************************************************************/
SimplexMesh& SimplexMesh::buildQuadFromNodes0(const std::vector<math::Point>& ptCoords, std::vector<TSimplexID>& marquedSimplex)
{
  if(ptCoords.size() == 8)
  {
    std::vector<math::Point> nodeDownCoord{ptCoords[0], ptCoords[1], ptCoords[2], ptCoords[3]};
    CriterionRAIS criterionRAIS(new VolumeCriterion());
    static int j = 0;

    std::cout << "SimplexMesh::buildQuadFromNodes " << j <<std::endl;
    std::vector<TInt>&& nodesDown = buildQuadFaceFromCoordPt(nodeDownCoord, marquedSimplex);

    std::vector<std::vector<TInt>> nodeDownAndUp0{};
    std::vector<TInt> nodesDownFace0           {nodesDown[0], nodesDown[1], nodesDown[2]};
    std::vector<math::Point> nodeUpCordFace0   {ptCoords[4], ptCoords[5], ptCoords[6]};
    std::cout << "buildPrismFromFace 0" <<std::endl;
    buildPrismFromFace(nodesDownFace0, nodeUpCordFace0, marquedSimplex, j);

    if(j == 3)
    {
      return *this;
    }
    j++;

    std::vector<TInt> nodesDownFace1{nodesDown[0], nodesDown[3], nodesDown[2]};
    std::vector<math::Point> nodeUpCordFace1{ptCoords[4], ptCoords[7], ptCoords[6]};
    std::vector<std::vector<TInt>> nodeDownAndUp1{};
    std::cout << "buildPrismFromFace 1" <<std::endl;
    buildPrismFromFace(nodesDownFace1, nodeUpCordFace1, marquedSimplex);
  }
  else
  {
    std::cout << "ptCoords.size() == 8" << std::endl;
  }
  return *this;
}
/*****************************************************************************/
SimplexMesh& SimplexMesh::buildQuadFromNodes(const std::vector<math::Point>& ptCoords, std::vector<TSimplexID>& markedSimplex)
{
  if(ptCoords.size() == 8)
  {
    std::cout << "SimplexMesh::buildQuadFromNodes" << std::endl;
    std::vector<math::Point> nodeDownCoord{ptCoords[0], ptCoords[1], ptCoords[2], ptCoords[3]};
    std::vector<math::Point> nodeUpCoord{ptCoords[4], ptCoords[5], ptCoords[6], ptCoords[7]};

    CriterionRAIS criterionRAIS(new VolumeCriterion());

    std::cout << "buildQuadFaceFromCoordPt0" << std::endl;
    std::vector<TInt>&& nodesDown    = buildQuadFaceFromCoordPt(nodeDownCoord, markedSimplex);
    std::cout << "buildQuadFaceFromCoordPt1" << std::endl;
    std::vector<TInt>&& nodesUp      = buildQuadFaceFromCoordPt(nodeUpCoord, markedSimplex);


    //return *this;

    EdgeInsertion(this, SimplicesNode(this, nodesDown[0]), SimplicesNode(this, nodesUp[0]), criterionRAIS, markedSimplex);
    EdgeInsertion(this, SimplicesNode(this, nodesDown[1]), SimplicesNode(this, nodesUp[1]), criterionRAIS, markedSimplex);
    EdgeInsertion(this, SimplicesNode(this, nodesDown[2]), SimplicesNode(this, nodesUp[2]), criterionRAIS, markedSimplex);
    EdgeInsertion(this, SimplicesNode(this, nodesDown[3]), SimplicesNode(this, nodesUp[3]), criterionRAIS, markedSimplex);

    return *this;

    std::vector<TInt> nodeSide0{nodesDown[0], nodesDown[1], nodesUp[1], nodesUp[0]};
    std::vector<TInt> nodeSide1{nodesDown[1], nodesUp[1], nodesUp[2], nodesDown[2]};
    std::vector<TInt> nodeSide2{nodesDown[2], nodesDown[3], nodesUp[3], nodesUp[2]};
    std::vector<TInt> nodeSide3{nodesDown[3], nodesUp[3], nodesUp[0], nodesDown[0]};


    buildQuadFaceFromNode(nodeSide0, markedSimplex);
    buildQuadFaceFromNode(nodeSide1, markedSimplex);
    buildQuadFaceFromNode(nodeSide2, markedSimplex);
    buildQuadFaceFromNode(nodeSide3, markedSimplex);

  }
  else
  {
    std::cout << "ptCoords.size() == 8" << std::endl;
  }
  return *this;
}

/*****************************************************************************/
TSimplexID SimplexMesh::getOppositeFace(const TInt ANodeID, const TSimplexID ATriID)
{
  size_t localIndex             = -1;
  const size_t indexBufferSize  =  3;
  const std::vector<TSimplexID> buffer = m_tri_nodes[ATriID];
  if(ATriID >= 0)
  {
    for(size_t i = 0; i< indexBufferSize ; i++)
    {
      TSimplexID node = buffer[i];
      if(node == ANodeID)
      {
        localIndex = i;
        break;
      }
    }

    if(localIndex == -1)
    {
      //si il n'existe ni triangle ni tetra opposé cherché
      return std::numeric_limits<int>::min();
    }
    else
    {
      return m_tri_adj[ATriID][localIndex];
    }
  }
}
/*---------------------------------------------------------------------------*/

void SimplexMesh::reorientAllTetra()
{
  for(TSimplexID Tet = 0; Tet < m_tet_nodes.size(); Tet++)
  {
    if(m_tet_ids[Tet] != 0)
    {
      simplicesCell::SimplicesCell(this, Tet).reorientTet();
    }
  }
}
/*---------------------------------------------------------------------------*/
TSimplexID SimplexMesh::buildOppositeFacesVector(const TSimplexID currentTriangle, const TInt currentNode)
{
  const size_t TriangleNodeSize = 3;
  const int errorId             = std::numeric_limits<int>::min();
  TInt       nodesTri[2]        = {errorId, errorId};
  bool flag                     = false;
  TSimplexID adjTriangle        = errorId;
  TInt iter                     = 0;

  for(TInt indexTri = 0; indexTri < TriangleNodeSize; indexTri++)
  {
    if(m_tri_nodes[currentTriangle][indexTri] != currentNode)
    {
      nodesTri[iter] = m_tri_nodes[currentTriangle][indexTri];
      iter++;
    }
  }

  for(TSimplexID compTri = 1; compTri < m_tri_nodes.size(); compTri++)
  {
    if(currentTriangle != compTri && m_tri_ids[compTri] != 0)
    {
      const std::vector<TSimplexID> vec = m_tri_nodes[compTri];
      const TInt arrayComp[3] = {vec[0], vec[1], vec[2]};
      flag = containNodes<2,3>(nodesTri, &arrayComp[0]);
      if(flag)
      {
        adjTriangle = compTri;
        break;
      }
    }
  }
  return adjTriangle;
}
/*---------------------------------------------------------------------------*/
bool SimplexMesh::checkSimplexContenaing(const gmds::math::Point& pt, /*std::vector<*/TSimplexID/*>*/& tetraContainingPt)
{
  bool flag = false;
  for(TSimplexID tet = 0; tet < m_tet_ids.capacity(); tet++)
  {
    if(m_tet_ids[tet] != 0)
    {
        if(SimplicesCell(this, tet).isPointInSimplex(pt))
        {
          tetraContainingPt = tet;
          flag = true;
        }
    }
  }
  return flag;
}
/*---------------------------------------------------------------------------*/
bool SimplexMesh::simplexContaining(const SimplicesNode& node, TSimplexID& tetraContainingNode)
{
  return simplexContaining(node.getCoords(), tetraContainingNode);
}
/*---------------------------------------------------------------------------*/
bool SimplexMesh::simplexContaining(const Point& pt, TSimplexID& tetraContainingPt)
{

  std::vector<bool>  markedSimplex(m_tet_ids.capacity(), false);
  for(TInt idx = 0; idx <  markedSimplex.size(); idx++)
  {
    if(m_tet_ids[idx] == 1)
    {
      markedSimplex[idx] = true;
    }
  }


  TSimplexID Tet          =  firstBitTo1(markedSimplex);
  TSimplexID previousTet  =  Tet;
  TSimplexID noTet        =  std::numeric_limits<int>::min();
  bool       flag         =  false;


  while(previousTet != noTet)
  {
    if(m_tet_ids[previousTet] != 0)
    {
        markedSimplex[previousTet] = false;
        flag               = isInSimplex(Tet, pt);
        if(flag)
        {
            flag       = true;
            break;
        }
        else if(Tet  == noTet || markedSimplex[Tet] == false)
        {
          Tet          =  firstBitTo1(markedSimplex);
        }
        previousTet = Tet;
    }
  }

  if(flag)
  {
      tetraContainingPt = Tet;
  }
  else
  {
    tetraContainingPt = noTet;
  }

  return flag;

}
/*---------------------------------------------------------------------------*/
bool SimplexMesh::checkSimplicesContenaing(const gmds::math::Point& pt, std::vector<TSimplexID>& tetraContainingPt, TSimplexID simplexToCheckFirst)
{
  TSimplexID border = std::numeric_limits<int>::min();
  TSimplexID nextTet = border;
  bool flag = false;
  double u = 0.0 ; double v = 0.0 ; double w = 0.0 ; double t = 0.0 ;
  closestSimplex closestSimplex;
  std::vector<math::Orientation::Sign> uvwt;
  uvwt.reserve(4);



  TSimplexID simplexNextToPt = (simplexToCheckFirst == border)?m_octree->findSimplexNextTo(pt):simplexToCheckFirst;
  BitVector cyclingCheck(getBitVectorTet().capacity());
  if(simplexNextToPt != border)
  {
    nextTet = nextSimplexToCheckOrientation(simplexNextToPt, pt, uvwt, cyclingCheck);
    if(nextTet >= 0)
    {
      tetraContainingPt.push_back(nextTet);
    }
  }
  else
  {
    gmds::BitVector m_tet_ids = getBitVectorTet();
    TSimplexID firstValidTet = 0;
    for(TSimplexID tet = 0; tet < m_tet_ids.capacity(); tet++)
    {
      if(m_tet_ids[tet] != 0)
      {
          firstValidTet = tet;
          break;
      }
    }

    for(TSimplexID tet = firstValidTet; tet < m_tet_ids.capacity(); tet++)
    {
      if(m_tet_ids[tet] != 0)
      {
        if(cyclingCheck[tet] == 0)
        {
          nextTet = nextSimplexToCheckOrientation(tet, pt, uvwt, cyclingCheck);
          if(nextTet < 0)
          {
            continue;
          }
          else
          {
            tetraContainingPt.push_back(nextTet);
            break;
          }
        }
      }
    }
  }




  /*if the node did not find a simplex contenaing it, we find the minimal projection
    of the node on the surface in order to find the initial cavity*/
  if(tetraContainingPt.size() == 0)
  {
    double distance = std::numeric_limits<double>::max();
    tetraContainingPt.resize(1);

    for(TSimplexID tri = 1; tri < m_tri_ids.capacity(); tri++)
    {
      if(m_tri_ids[tri] != 0)
      {
        double minProj = 0.0;
        math::Point projectedPoint;

        SimplicesTriangle triangle = SimplicesTriangle(this, tri);
        std::vector<TInt> nodes = triangle.getNodes();
        math::Point pt0 = SimplicesNode(this, nodes[0]).getCoords();
        math::Point pt1 = SimplicesNode(this, nodes[1]).getCoords();
        math::Point pt2 = SimplicesNode(this, nodes[2]).getCoords();

        if(pointInTriangle(pt, pt0, pt1, pt2, minProj, projectedPoint))
        {
          minProj = std::fabs(minProj);
          if(minProj < distance)
          {
            tetraContainingPt[0] = m_tri_nodes[tri][3];
            distance = minProj;
          }
        }
      }
    }
    /*Now the simplex have been found, we will calculate the Orientation::Sign in order to find the right cavity*/
    SimplicesCell currentSimplex(this,tetraContainingPt.front());
    for(unsigned int faceLocal = 0 ; faceLocal < 4 ; faceLocal++)
    {
      const std::vector<TInt>&& faceNodes = currentSimplex.getOrderedFace(faceLocal);
      uvwt[faceLocal] = currentSimplex.orientation(faceLocal, pt, true);
    }
  }

  /*Now uvwt is found for every node (surface & volume) we can correct the initial cavity with the Orientation::Sign*/
  if(tetraContainingPt.size() != 0)
  {
    flag = true;
    std::vector<int> neighborFace{};
    SimplicesCell currentSimplex(this,tetraContainingPt.front());
    for(unsigned int faceLocal = 0 ; faceLocal < 4 ; faceLocal++)
    {
      math::Orientation::Sign sign = uvwt[faceLocal];
      if(sign == math::Orientation::Sign::ZERO)
      {
        neighborFace.push_back(faceLocal);
      }
    }

    if(neighborFace.size() != 0)
    {
      if(neighborFace.size() == 1)
      {
        TSimplexID oppositeSimplex = currentSimplex.oppositeTetraIdx(neighborFace.front());
        if(oppositeSimplex != border & oppositeSimplex >= 0)
        {
          tetraContainingPt.push_back(oppositeSimplex);
        }
      }
      else if(neighborFace.size() == 2)
      {
        std::vector<TInt> nodesEdge{};
        gmds::BitVector bitNeighborFace(4);
        for(auto const & face : neighborFace){bitNeighborFace.assign(face);}
        for(unsigned int faceComp = 0 ; faceComp < 4 ; faceComp++)
        {
          if(bitNeighborFace[faceComp] == 0)
          {
            nodesEdge.push_back(faceComp);
          }
        }
        if(nodesEdge.size() == 2)
        {
          const SimplicesNode nodeA = currentSimplex.getNode(nodesEdge.front());
          const SimplicesNode nodeB = currentSimplex.getNode(nodesEdge.back());
          const std::vector<TSimplexID>&& shell = nodeA.shell(nodeB);
          tetraContainingPt.clear();
          std::copy_if(shell.begin(), shell.end(), std::back_inserter(tetraContainingPt), [&](const TSimplexID simplex)
          {
            return (simplex >= 0);
          });
          //tetraContainingPt = nodeA.shell(nodeB);
        }
        else
        {
          //TODO assertion
          std::cout << "nodesEdge.size() == 2 not valid" << std::endl;
        }
      }
      else if(neighborFace.size() == 3)
      {
        TInt nodeBallLocal;
        gmds::BitVector bitNeighborFace(4);
        for(auto const & face : neighborFace){bitNeighborFace.assign(face);}
        for(unsigned int faceComp = 0 ; faceComp < 4 ; faceComp++)
        {
          if(bitNeighborFace[faceComp] == 0)
          {
            nodeBallLocal = faceComp;
            break;
          }
        }
        const SimplicesNode nodeA = currentSimplex.getNode(nodeBallLocal);
        const std::vector<TSimplexID>&& ball = nodeA.ballOf();
        tetraContainingPt.clear();
        std::copy_if(ball.begin(), ball.end(), std::back_inserter(tetraContainingPt), [&](const TSimplexID simplex)
        {
          return (simplex >= 0);
        });
      }
    }
  }

  return flag;
}
/*---------------------------------------------------------------------------*/
TSimplexID SimplexMesh::nextSimplexToCheck(const TSimplexID currentSimplex, const math::Point& pt, double& u, double& v, double& w, double& t, BitVector& cyclingCheck, closestSimplex& closestSimplexInfo)
{
  TSimplexID cycling = std::numeric_limits<int>::min();
  TSimplexID border  = cycling;
  TSimplexID nextTet = currentSimplex;

  if(cyclingCheck[currentSimplex] == 0)
  {
    const double epsilon = 0.0;//-10E-4;
    if(u >= epsilon && v >= epsilon && w >= epsilon && t >= epsilon)
    {
      return currentSimplex;
    }
    else if(u < epsilon)
    {
      nextTet = SimplicesCell(this, currentSimplex).oppositeTetraIdx(0);
    }
    else if(v < epsilon)
    {
      nextTet = SimplicesCell(this, currentSimplex).oppositeTetraIdx(1);
    }
    else if(w < epsilon)
    {
      nextTet = SimplicesCell(this, currentSimplex).oppositeTetraIdx(2);
    }
    else
    {
      nextTet = SimplicesCell(this, currentSimplex).oppositeTetraIdx(3);
    }

    if(nextTet == border)
    {
      return border;
    }
    cyclingCheck.assign(currentSimplex);
  }
  else
  {
    return cycling;
  }

  u = SimplicesCell(this, nextTet).signedBarycentric(0, pt);
  v = SimplicesCell(this, nextTet).signedBarycentric(1, pt);
  w = SimplicesCell(this, nextTet).signedBarycentric(2, pt);
  t = SimplicesCell(this, nextTet).signedBarycentric(3, pt);

  std::vector<TInt> nodes = SimplicesCell(this, nextTet).getNodes();
  unsigned int proximFace = 0;
  double minLenght = std::numeric_limits<double>::max();
  double distance = 0.0;
  for(unsigned int face = 0 ; face < 4 ; face++)
  {
    const math::Point p0 = SimplicesNode(this,nodes[face]).getCoords();
    const math::Point p1 = SimplicesNode(this,nodes[(face+1) % 4]).getCoords();
    const math::Point p2 = SimplicesNode(this,nodes[(face+2) % 4]).getCoords();
    math::Point projectedPoint;
    if(pointInTriangle(pt, p0, p1, p2, distance, projectedPoint))
    {
      math::Vector3d u = p1 - p0;
      math::Vector3d v = p2 - p0;
      math::Vector3d n = u.cross(v);
      n.normalize();

      math::Vector3d w = pt - p0;
      double distancePtTriangle = std::fabs(w.dot(n));
      if(distancePtTriangle < minLenght)
      {
        proximFace = face;
        minLenght = distancePtTriangle;
      }
    }

  }

  /*
  const std::vector<TInt>&& nodes = SimplicesCell(this, nextTet).getNodes();
  const math::Point && nodeACoord = SimplicesNode(this, nodes[0]).getCoords();
  const math::Point && nodeBCoord = SimplicesNode(this, nodes[1]).getCoords();
  const math::Point && nodeCCoord = SimplicesNode(this, nodes[2]).getCoords();
  const math::Point && nodeDCoord = SimplicesNode(this, nodes[3]).getCoords();

  const math::Point faceBarA = 1.0 / 3.0 * (nodeACoord + nodeBCoord + nodeCCoord);
  const math::Point faceBarB = 1.0 / 3.0 * (nodeACoord + nodeBCoord + nodeDCoord);
  const math::Point faceBarC = 1.0 / 3.0 * (nodeACoord + nodeDCoord + nodeCCoord);
  const math::Point faceBarD = 1.0 / 3.0 * (nodeCCoord + nodeBCoord + nodeDCoord);

  const math::Vector3d A_pt = faceBarA - pt;
  const math::Vector3d B_pt = faceBarB - pt;
  const math::Vector3d C_pt = faceBarC - pt;
  const math::Vector3d D_pt = faceBarD - pt;

  const double lenghtA = A_pt.norm();
  const double lenghtB = B_pt.norm();
  const double lenghtC = C_pt.norm();
  const double lenghtD = D_pt.norm();

  const double minLenght = std::min(lenghtA, std::min(lenghtB,std::min(lenghtC, lenghtD)));
  */
  /*
  const math::Vector3d A_pt = nodeACoord - pt;
  const math::Vector3d B_pt = nodeBCoord - pt;
  const math::Vector3d C_pt = nodeCCoord - pt;
  const math::Vector3d D_pt = nodeDCoord - pt;

  const double lenghtA = A_pt.norm();
  const double lenghtB = B_pt.norm();
  const double lenghtC = C_pt.norm();
  const double lenghtD = D_pt.norm();

  const double minLenght = std::min(lenghtA, std::min(lenghtB,std::min(lenghtC, lenghtD)));
  */

  if(minLenght < closestSimplexInfo.minLenghtPtToNodeOfSimplex)
  {
    closestSimplexInfo.closeSimplex = nextTet;
    closestSimplexInfo.minLenghtPtToNodeOfSimplex = minLenght;
    closestSimplexInfo.u = u ; closestSimplexInfo.v = v ; closestSimplexInfo.w = w ; closestSimplexInfo.t = t ;
  }

  return nextSimplexToCheck(nextTet, pt, u, v, w, t, cyclingCheck, closestSimplexInfo);
}
/*---------------------------------------------------------------------------*/
TSimplexID SimplexMesh::nextSimplexToCheckOrientation(const TSimplexID currentSimplex, const math::Point& pt, std::vector<math::Orientation::Sign>& uvwt, BitVector& cyclingCheck)
{
  TSimplexID cycling = std::numeric_limits<int>::min();
  TSimplexID border  = cycling;
  TSimplexID nextTet = currentSimplex;
  uvwt[0] = SimplicesCell(this, currentSimplex).orientation(0, pt);
  uvwt[1] = SimplicesCell(this, currentSimplex).orientation(1, pt);
  uvwt[2] = SimplicesCell(this, currentSimplex).orientation(2, pt);
  uvwt[3] = SimplicesCell(this, currentSimplex).orientation(3, pt);

  if(cyclingCheck[currentSimplex] == 0)
  {
    if(uvwt[0] > 0  && uvwt[1]  > 0 && uvwt[2]  > 0 && uvwt[3]  > 0)
    {
      return currentSimplex;
    }
    else if(uvwt[0] < 0)
    {
      nextTet = SimplicesCell(this, currentSimplex).oppositeTetraIdx(0);
    }
    else if(uvwt[1] < 0)
    {
      nextTet = SimplicesCell(this, currentSimplex).oppositeTetraIdx(1);
    }
    else if(uvwt[2] < 0)
    {
      nextTet = SimplicesCell(this, currentSimplex).oppositeTetraIdx(2);
    }
    else if(uvwt[3] < 0)
    {
      nextTet = SimplicesCell(this, currentSimplex).oppositeTetraIdx(3);
    }

    cyclingCheck.assign(currentSimplex);
    if(nextTet < 0)
    {
      return nextTet;
    }
  }
  else
  {
    return cycling;
  }



  return nextSimplexToCheckOrientation(nextTet, pt, uvwt, cyclingCheck);
}
/*---------------------------------------------------------------------------*/
TSimplexID SimplexMesh::firstBitTo1(const std::vector<bool>& vec)
{
  int errorId     = std::numeric_limits<int>::min();
  bool flag       = false;
  TSimplexID res  = 0;

  for(TSimplexID idx  = 0; idx < vec.size(); idx++)
  {
    if(vec[idx]  == true)
    {
      flag = true;
      res = idx;
      break;
    }
  }

  return (flag == true )? res : errorId;
}
/*---------------------------------------------------------------------------*/
bool SimplexMesh::isInSimplex(TSimplexID& tetra, const gmds::math::Point& pt)
{
  bool flag = false;
  TSimplexID previousTet = tetra;
  int errorId = std::numeric_limits<int>::min();
  double epsilon = 10E-35;
  if(tetra <= m_tet_ids.size())
  {
    if(m_tet_ids[tetra] != 0)
    {

        double u = SimplicesCell(this, tetra).signedBarycentricNormalized(0,pt);
        double v = SimplicesCell(this, tetra).signedBarycentricNormalized(1,pt);
        double w = SimplicesCell(this, tetra).signedBarycentricNormalized(2,pt);
        double t = SimplicesCell(this, tetra).signedBarycentricNormalized(3,pt);


        /*TODO à optimiser...*/
        if(u >= 0.0)
        {
          if(v >= 0.0)
          {
            if(w >= 0.0)
            {
              if(t >= 0.0)
              {
                flag = true;
              }
              else if(t < 0.0)
              {
                tetra = m_tet_adj[tetra][3];
                flag  = false;
              }
              /*else if(t == 0.0)
              {
                //tetra = m_tet_adj[tetra][3];
                flag  = true;
              }*/
            }
            /*else if(w == 0.0)
            {
              //tetra = m_tet_adj[tetra][2];
              flag = true;
            }*/
            else if(w < 0.0)
            {
              tetra = m_tet_adj[tetra][2];
              flag = false;
            }
          }
          /*else if(v == 0.0)
          {
            //tetra = m_tet_adj[tetra][1];
            flag = true;
          }*/
          else if(v < 0.0)
          {
            tetra = m_tet_adj[tetra][1];
            flag = false;
          }
        }
        /*else if(u == 0.0)
        {
          //tetra = m_tet_adj[tetra][0];
          flag = true;
        }*/
        else if(u < 0.0)
        {
          tetra = m_tet_adj[tetra][0];
          flag = false;
        }
    }

    tetra = (tetra >= 0 || tetra == errorId)? tetra :
            (m_tri_nodes[-tetra][3] != previousTet || m_tri_nodes[-tetra][3] == errorId)? m_tri_nodes[-tetra][3] :
             m_tri_adj[-tetra][3];
  }

  return flag;
}
/*---------------------------------------------------------------------------*/
TInt SimplexMesh::nextNode()
{
  TInt res = 0;
  while(m_node_ids[nodeIndx] != 1 )
  {
    nodeIndx++;
    if(nodeIndx > m_node_ids.capacity())
    {
      break;
    }
  }

  res = nodeIndx;
  nodeIndx++;
  return res;
}
/*---------------------------------------------------------------------------*/
TInt SimplexMesh::getFirstNode()
{
  nodeIndx = 0;
  return nextNode();
}
/*---------------------------------------------------------------------------*/
TInt SimplexMesh::nextTet()
{
  TInt res = 0;
  while(m_tet_ids[tetIndx] != 1 )
  {
    tetIndx++;
    if(tetIndx > m_tet_ids.capacity())
    {
      break;
    }
  }

  res = tetIndx;
  tetIndx++;
  return res;
}
/*---------------------------------------------------------------------------*/
TInt SimplexMesh::getFirstTet()
{
  tetIndx = 0;
  return nextTet();
}
/*---------------------------------------------------------------------------*/

TInt SimplexMesh::nextTri()
{
  TInt res = 0;
  while(m_tri_ids[triIndx] != 1 )
  {
    triIndx++;
    if(triIndx > m_tri_ids.capacity())
    {
      break;
    }
  }

  res = triIndx;
  triIndx++;
  return res;
}
/*---------------------------------------------------------------------------*/
TInt SimplexMesh::getFirstTri()
{
  triIndx = 0;
  return nextTri();
}
/*****************************************************************************/
std::set<TSimplexID> SimplexMesh::simplexInBorder()
{
  std::set<TSimplexID> s;
  unsigned int simplexNodeSize = 4;
  for (TSimplexID tetIndx = 0 ; tetIndx < m_tet_nodes.size() ; tetIndx++)
  {
    if(m_tet_ids[tetIndx] != 0)
    {
      for(TInt node = 0; node < simplexNodeSize; node++)
      {
          if(SimplicesCell(this, tetIndx).inBorder())
          {
            s.insert(tetIndx);
            break;
          }
      }
    }
  }


  return std::move(s);
}
/******************************************************************************/
std::vector<TSimplexID> SimplexMesh::intersectionSimplex(std::vector<std::vector<TSimplexID>> & vectorOfBalls)
{

  for(auto & vector : vectorOfBalls)
  {
    std::sort(vector.begin(), vector.end());
  }

  std::vector<TSimplexID> intersectionVec = vectorOfBalls[0];
  std::vector<TSimplexID> tempBufferSimplex;
  for(int idxVector = 0 ; idxVector < vectorOfBalls.size() - 1 ; idxVector++)
  {
    std::set_intersection(intersectionVec.begin(), intersectionVec.end(),
                          vectorOfBalls[idxVector + 1].begin(), vectorOfBalls[idxVector + 1].end(),
                          std::back_inserter(tempBufferSimplex));
    intersectionVec.clear();
    intersectionVec = tempBufferSimplex;
    tempBufferSimplex.clear();
  }

  return std::move(intersectionVec);
}
/******************************************************************************/
gmds::BitVector SimplexMesh::FindTetInHexa(const std::vector<std::vector<TInt>>& nodesQuads) const
{
  gmds::BitVector tetInBitVector;

}
/******************************************************************************/
void SimplexMesh::buildQuadEdgeFromNodes(const std::vector<TInt>& nodes, const CriterionRAIS& criterion)
{
  unsigned int quadNodeSize = 4;
  if(nodes.size() == quadNodeSize)
  {
    for(unsigned int nodeIter = 0; nodeIter < nodes.size() ; nodeIter++)
    {
      std::vector<TSimplexID>mS{};
      EdgeInsertion(this, SimplicesNode(this, nodes[nodeIter]), SimplicesNode(this, nodes[(nodeIter + 1) % quadNodeSize]), criterion, mS);
    }
  }
  else
  {
    //TODO exception
    std::cout << "nodes.size() != quadNodeSize" << std::endl;
  }
}
/******************************************************************************/
void SimplexMesh::buildHexaEdgeFromNodes(const std::vector<TInt>& nodes, const operators::CriterionRAIS& criterion)
{
  unsigned int hexaNodeSize = 8;
  std::vector<std::vector<TInt>> nodesQuad{
    {nodes[0], nodes[1], nodes[2], nodes[3]},
    {nodes[4], nodes[5], nodes[6], nodes[7]},
    {nodes[0], nodes[1], nodes[5], nodes[4]},
    {nodes[1], nodes[5], nodes[6], nodes[2]},
    {nodes[2], nodes[3], nodes[7], nodes[6]},
    {nodes[3], nodes[7], nodes[4], nodes[0]}
  };

  if(nodes.size() == hexaNodeSize)
  {
    for(auto const & nodeQuad : nodesQuad)
    {
      buildQuadEdgeFromNodes(nodeQuad, criterion);
      //return;
    }
  }
  else
  {
    //TODO exception
    std::cout << "nodes.size() != hexaNodeSize" << std::endl;
  }

  return;
}
/******************************************************************************/
void SimplexMesh::buildFacesFromEdges(const std::vector<TInt>& nodes, const operators::CriterionRAIS& criterion)
{
  unsigned int sizeFaceEdge = 3;
  if(nodes.size() == sizeFaceEdge)
  {
      std::vector<TSimplexID>&& ballA = SimplicesNode(this, nodes[0]).ballOf();
      std::vector<TSimplexID>&& ballB = SimplicesNode(this, nodes[1]).ballOf();
      std::vector<TSimplexID>&& ballC = SimplicesNode(this, nodes[2]).ballOf();

      std::vector<std::vector<TSimplexID>>&& ballsABC{
        ballA, ballB, ballC
      };

      std::vector<TSimplexID>&& intersectionBallABC = intersectionSimplex(ballsABC);
      if(intersectionBallABC.size() == 0) // no face existance
      {
        std::vector<TSimplexID>&& shellAB = SimplicesNode(this, nodes[0]).shell(SimplicesNode(this, nodes[1]));
        std::vector<TSimplexID>&& shellBC = SimplicesNode(this, nodes[1]).shell(SimplicesNode(this, nodes[2]));

        if(shellAB.size() != 0 && shellBC.size() != 0)
        {
          math::Point centerEdge = 0.5 * (SimplicesNode(this, nodes[0]).getCoords() + SimplicesNode(this, nodes[1]).getCoords());
          math::Point centerFace = 1.0/3.0 * (SimplicesNode(this, nodes[0]).getCoords() + SimplicesNode(this, nodes[1]).getCoords() + SimplicesNode(this, nodes[2]).getCoords());
          math::Vector3d normalEdge =  centerFace - centerEdge;
          normalEdge.normalize();
          bool flag = false;
          double epsilon = 10E-5f;
          std::vector<TSimplexID> && shell = SimplicesNode(this, nodes[0]).shell(SimplicesNode(this, nodes[1]));
          std::vector<TSimplexID> cavity{};

          while(flag != true)
          {
            math::Point NodeNeighborEdge = centerEdge + epsilon * normalEdge;
            for(auto & tet : shell)
            {
              if(SimplicesCell(this,tet).isPointInSimplex(NodeNeighborEdge))
              {
                flag = true;
                cavity.push_back(tet);
                break;
              }
            }
            epsilon *= 0.1;
          }

          if(cavity.size() == 1)
          {
            //PointInsertion pi(this, SimplicesNode(this, nodes[2]), criterion, cavity);
          }
          else
          {
            //TODO
            std::cout << "cavity.size() == 1" << std::endl;
          }
        }
      }
  }
  else
  {
    //TODO exception
    std::cout << "nodes.size() != sizeFaceEdge" << std::endl;
  }
  return;
}
/******************************************************************************/
void SimplexMesh::buildHexaFaceFromNodes(const std::vector<TInt>& nodes, const operators::CriterionRAIS& criterion)
{
  unsigned int quadNodeSize = 4;
  if(nodes.size() == quadNodeSize)
  {
    std::vector<TSimplexID>&& ball0 = SimplicesNode(this, nodes[0]).ballOf();
    std::vector<TSimplexID>&& ball1 = SimplicesNode(this, nodes[1]).ballOf();
    std::vector<TSimplexID>&& ball2 = SimplicesNode(this, nodes[2]).ballOf();
    std::vector<TSimplexID>&& ball3 = SimplicesNode(this, nodes[3]).ballOf();

    std::vector<std::vector<TSimplexID>> ball012{
      ball0, ball1, ball2
    };

    std::vector<std::vector<TSimplexID>> ball032{
      ball0, ball3, ball2
    };

    std::vector<std::vector<TSimplexID>> ball013{
      ball0, ball1, ball3
    };

    std::vector<std::vector<TSimplexID>> ball123{
      ball1, ball2, ball3
    };


    std::vector<TSimplexID>&& intersection012 = intersectionSimplex(ball012);
    std::vector<TSimplexID>&& intersection032 = intersectionSimplex(ball032);
    std::vector<TSimplexID>&& intersection013 = intersectionSimplex(ball013);
    std::vector<TSimplexID>&& intersection123 = intersectionSimplex(ball123);

    if((intersection012.size() > 0 && intersection032.size() > 0)||
       (intersection013.size() > 0 && intersection123.size() > 0)) //the faces are already build
    {
      //the faces are already build
    }
    else
    {
      if(intersection012.size() == 0 && intersection032.size() == 0)
      {
        if(intersection013.size() == 0 && intersection123.size() == 0)
        {
          std::vector<TSimplexID> nodeFace0{nodes[0], nodes[1], nodes[2]};
          std::vector<TSimplexID> nodeFace1{nodes[0], nodes[3], nodes[2]};
          buildFacesFromEdges(nodeFace0, criterion);
          buildFacesFromEdges(nodeFace1, criterion);
        }
      }
      else if(intersection012.size() == 0 && intersection032.size() > 0)
      {
        std::vector<TSimplexID> nodeFace{nodes[0], nodes[1], nodes[2]};
        buildFacesFromEdges(nodeFace, criterion);
      }
      else if(intersection012.size() > 0 && intersection032.size() == 0)
      {
        std::vector<TSimplexID> nodeFace{nodes[0], nodes[3], nodes[2]};
        buildFacesFromEdges(nodeFace, criterion);
      }

      else if(intersection013.size() == 0 && intersection123.size() > 0)
      {
        std::vector<TSimplexID> nodeFace{nodes[0], nodes[1], nodes[3]};
        buildFacesFromEdges(nodeFace, criterion);
      }
      else if(intersection013.size() > 0 && intersection123.size() == 0)
      {
        std::vector<TSimplexID> nodeFace{nodes[1], nodes[2], nodes[3]};
        buildFacesFromEdges(nodeFace, criterion);
      }
      else
      {

      }
    }
  }
  else
  {
    //TODO Exeption
    std::cout << "nodes.size() == quadNodeSize" << std::endl;
  }
}
/******************************************************************************/
void SimplexMesh::buildHexaedre(const std::vector<TInt>& nodes, const operators::CriterionRAIS& criterion)
{
  if(nodes.size() == 8)
  {
    buildHexaEdgeFromNodes(nodes, criterion);
    // TODO check if all the edge exist !!
    std::vector<std::vector<TInt>> nodesQuads{
      {nodes[0], nodes[1], nodes[2], nodes[3]},
      {nodes[4], nodes[5], nodes[6], nodes[7]},
      {nodes[0], nodes[1], nodes[5], nodes[4]}, //bas
      {nodes[1], nodes[5], nodes[6], nodes[2]}, //droite
      {nodes[2], nodes[3], nodes[7], nodes[6]}, //haut
      {nodes[3], nodes[7], nodes[4], nodes[0]} // gauche
    };

    for(auto const & nodesQuad : nodesQuads)
    {
      buildHexaFaceFromNodes(nodesQuad, criterion);
    }
  }
  else
  {
    std::cout << "nodes.size() != 8" << std::endl;
  }
}
/******************************************************************************/
bool SimplexMesh::MollerTriangleIntersection(const simplicesNode::SimplicesNode& nodeA, const math::Vector3d& dir, const std::vector<TInt>& triangleNodeId, math::Vector3d& tuv)
{
  if(triangleNodeId.size() == 3)
  {
    //Real moller intersection
    math::Vector3d ray_direction = dir;
    math::Vector3d ray_origin    = nodeA.getCoords();

    //on regarde maintenant si il y'a intersection... Möller..
    math::Vector3d E1 = SimplicesNode(this, triangleNodeId[1]).getCoords() - SimplicesNode(this, triangleNodeId[0]).getCoords();
    math::Vector3d E2 = SimplicesNode(this, triangleNodeId[2]).getCoords() - SimplicesNode(this, triangleNodeId[0]).getCoords();
    math::Vector3d T = ray_origin - SimplicesNode(this, triangleNodeId[0]).getCoords();


    math::Vector3d DxE2 = ray_direction.cross(E2);
    double den = DxE2.dot(E1);
    if (den > -10E-20 && den < 10E-20) // ray and triangle are parallele
    {
      return false;
    }

    math::Vector3d TxE1 = T.cross(E1);
    double compo1 = TxE1.dot(E2);
    double compo2 = DxE2.dot(T);
    double compo3 = TxE1.dot(ray_direction);
    tuv = (1.0/den) * math::Vector3d(compo1,compo2,compo3);


    //condition d'intersection..
    if(tuv[0] < 0.0)
    {
      std::cout << "tuv[0] < 0" << std::endl;
      return false;
    }

    double epsilon = 10E-14;

    if(tuv[1] < -epsilon  || tuv[1] > 1.0 + epsilon )
    {
      std::cout << "tuv[1] : " << tuv[1] << std::endl;
      std::cout << std::endl;
      return false;
    }
    else if(tuv[2] < -epsilon || tuv[2] > 1.0 + epsilon )
    {
      std::cout << "tuv[2] : " << tuv[2] << std::endl;
      std::cout << std::endl;
      return false;
    }
    else
    {
      return true;
    }
  }
  else
  {
    //todo exception
    std::cout << "triangleNodeId.size() != 3" << std::endl;
  }
  return false;
}
/******************************************************************************/
void SimplexMesh::hexaBuildPerformance(const std::vector<std::vector<TInt>>& orderedNodesQuads, double& edgePerformance, double& hexaPerformance)
{
  unsigned int hexaWishNbr = orderedNodesQuads.size();
  unsigned int edgeWishNbr = hexaWishNbr * 12;
  std::vector<double> res{};
  unsigned int hexaBuild = 0;
  unsigned int edgeBuild = 0;

  if(orderedNodesQuads.size())
  {
    for(auto const nodes : orderedNodesQuads)
    {
      unsigned int facenbr = 0;
      if(nodes.size() != 8)
      {
        //TODO Exception
        std::cout << "nodes.size() != 8 in hexaBuildPerformance" << std::endl;
        return ;
      }
      if(m_node_ids[nodes[0]] != 0 && m_node_ids[nodes[1]] != 0 && m_node_ids[nodes[2]] != 0 && m_node_ids[nodes[3]] != 0 &&
         m_node_ids[nodes[4]] != 0 && m_node_ids[nodes[5]] != 0 && m_node_ids[nodes[6]] != 0 && m_node_ids[nodes[7]] != 0 )
         {
           std::vector<TInt> nodesFace0 {nodes[0], nodes[1], nodes[2], nodes[3]}; // front face
           std::vector<TInt> nodesFace1 {nodes[4], nodes[5], nodes[6], nodes[7]}; // back face
           std::vector<TInt> nodesFace2 {nodes[0], nodes[1], nodes[5], nodes[4]};
           std::vector<TInt> nodesFace3 {nodes[1], nodes[5], nodes[6], nodes[2]};
           std::vector<TInt> nodesFace4 {nodes[2], nodes[6], nodes[7], nodes[3]};
           std::vector<TInt> nodesFace5 {nodes[3], nodes[7], nodes[4], nodes[0]};

           for(unsigned int nodeIter = 0; nodeIter < nodesFace0.size() ; nodeIter++) //premiere face
           {
             std::vector<TSimplexID>&& shell = SimplicesNode(this, nodesFace0[nodeIter]).shell(SimplicesNode(this, nodesFace0[(nodeIter + 1) % (nodesFace0.size())]));
             if(shell.size())
             {
               edgeBuild++;
             }
           }

           for(unsigned int nodeIter = 0; nodeIter < nodesFace1.size() ; nodeIter++) //seconde face
           {
             std::vector<TSimplexID>&& shell = SimplicesNode(this, nodesFace1[nodeIter]).shell(SimplicesNode(this, nodesFace1[(nodeIter + 1) % nodesFace1.size() ]));
             if(shell.size())
             {
               edgeBuild++;
             }
           }


           std::vector<TSimplexID>&& shell04 = SimplicesNode(this, nodes[0]).shell(SimplicesNode(this, nodes[4]));
           std::vector<TSimplexID>&& shell15 = SimplicesNode(this, nodes[1]).shell(SimplicesNode(this, nodes[5]));
           std::vector<TSimplexID>&& shell26 = SimplicesNode(this, nodes[2]).shell(SimplicesNode(this, nodes[6]));
           std::vector<TSimplexID>&& shell37 = SimplicesNode(this, nodes[3]).shell(SimplicesNode(this, nodes[7]));

           (shell04.size())?edgeBuild++ : 1;
           (shell15.size())?edgeBuild++ : 1;
           (shell26.size())?edgeBuild++ : 1;
           (shell37.size())?edgeBuild++ : 1;

           if(isFaceBuild(nodesFace0) && isFaceBuild(nodesFace1) && isFaceBuild(nodesFace2) && isFaceBuild(nodesFace3) && isFaceBuild(nodesFace4) && isFaceBuild(nodesFace5))
           {
             hexaBuild++;
           }
         }
    }
    edgePerformance = static_cast<double>(edgeBuild) / static_cast<double>(edgeWishNbr) * 100.0;
    hexaPerformance = static_cast<double>(hexaBuild) / static_cast<double>(hexaWishNbr) * 100.0;

  }
}
/******************************************************************************/
bool SimplexMesh::isFaceBuild(const std::vector<TInt>& nodes)
{
  bool flag = false;
  if(nodes.size() == 4)
  {
    std::vector<TSimplexID> && ball0 = SimplicesNode(this, nodes[0]).ballOf();
    std::vector<TSimplexID> && ball1 = SimplicesNode(this, nodes[1]).ballOf();
    std::vector<TSimplexID> && ball2 = SimplicesNode(this, nodes[2]).ballOf();
    std::vector<TSimplexID> && ball3 = SimplicesNode(this, nodes[3]).ballOf();

    std::vector<std::vector<TSimplexID>> bal012{{ball0}, {ball1}, {ball2}};
    std::vector<std::vector<TSimplexID>> bal230{{ball2}, {ball3}, {ball0}};
    std::vector<TSimplexID>&& inter012  = intersectionSimplex(bal012);
    std::vector<TSimplexID>&& inter230  = intersectionSimplex(bal230);

    std::vector<std::vector<TSimplexID>> bal013{{ball0}, {ball1}, {ball3}};
    std::vector<std::vector<TSimplexID>> bal321{{ball3}, {ball2}, {ball1}};
    std::vector<TSimplexID>&& inter013  = intersectionSimplex(bal013);
    std::vector<TSimplexID>&& inter321  = intersectionSimplex(bal321);

    if((inter012.size() && inter230.size()) || (inter013.size() && inter321.size()))
    {
      flag = true;
    }
  }
  else
  {
      std::cout << "nodes.size() != 4" << std::endl;
  }
  return flag;
}
/******************************************************************************/
void SimplexMesh::rebuildCavity(CavityOperator::CavityIO& cavityIO, const TInt nodeToConnect)
{
  const std::vector<std::vector<TInt>>&         pointsToConnect          = cavityIO.getNodesToReconnect();
  const std::vector<TSimplexID>&                extSimplexBorder         = cavityIO.getOppositesSimplex();
  const std::vector<TSimplexID>&                trianglesConnectedToP    = cavityIO.getTrianglesConnectedToPInCavity();
  const std::vector<std::vector<int>>&          localsInfoForTriangles   = cavityIO.getTriangleReconstructionInfo();
  const std::vector<std::vector<TSimplexID>>&   oppositeTriangles        = cavityIO.getoppositeTriangle();
  const std::vector<std::vector<TInt>>&         trianglesIndices         = cavityIO.getTriangleIndices();

  unsigned int sizeFace = 3;
  std::vector<TSimplexID> tetIdToBuildInfo{};
  std::vector<TSimplexID> triangles{};

  std::vector<TInt> nodesBall{};
  std::vector<TSimplexID> ballNodes{};
  TInt border = std::numeric_limits<int>::min();
  unsigned int nbrTetToAdd = extSimplexBorder.size();
  gmds::BitVector triConnectedToP(triCapacity());

  gmds::Variable<int>* BND_TRIANGLES     = getVariable<int, SimplicesTriangle>("BND_TRIANGLES");

  for(auto const simplex : trianglesConnectedToP)
  {
    if(triConnectedToP[-simplex] == 0)
    {
      triConnectedToP.assign(-simplex);
    }
  }
  auto my_make = [=](const TInt a, const TInt b){
    return (a < b)? std::make_pair(a, b) : std::make_pair(b, a);
  };

  std::multimap<std::pair<TInt, TInt>, TSimplexID> hash_face; // multimap  for the reconnection of the tetraedron
  std::map<TInt, std::vector<TSimplexID>> hash_edge; // map  for the reconnection of the surface triangle

  unsigned int faceIter = 0;

  for(auto const& face : pointsToConnect)
  {

    TSimplexID neighborSimplex         = extSimplexBorder[faceIter];

    if(neighborSimplex != border)
    {
      if(neighborSimplex < 0)
      {
        if(triConnectedToP[-neighborSimplex] == 1)
        {
          faceIter++;
          tetIdToBuildInfo.push_back(border);
          continue;
        }
      }
    }

    TSimplexID simplex = addTetraedre(face[0], face[1], face[2], nodeToConnect, false);
    SimplicesCell cell0(this, simplex);
    std::vector<TSimplexID> vecAdj{border, border, border, border};
    m_tet_adj[simplex] = vecAdj;
    m_tet_adj[simplex][3] = neighborSimplex;
    tetIdToBuildInfo.push_back(simplex);

    std::pair<TInt, TInt> key_1 = my_make(face[0], face[1]);
    std::pair<TInt, TInt> key_2 = my_make(face[1], face[2]);
    std::pair<TInt, TInt> key_3 = my_make(face[2], face[0]);

    hash_face.insert(std::pair<std::pair<TInt, TInt>, TSimplexID>(key_1, simplex));
    hash_face.insert(std::pair<std::pair<TInt, TInt>, TSimplexID>(key_2, simplex));
    hash_face.insert(std::pair<std::pair<TInt, TInt>, TSimplexID>(key_3, simplex));
    if(neighborSimplex != border)
    {
      if(neighborSimplex >= 0)
      {
        std::vector<TInt>&& nodesNeighborTet = SimplicesCell(this, neighborSimplex).getOtherNodeInSimplex(face);
        if(nodesNeighborTet.size() == 1)
        {
          TInt nodeLocalNeighborSimplex = SimplicesCell(this, neighborSimplex).getLocalNode(nodesNeighborTet.front());
          m_tet_adj[neighborSimplex][nodeLocalNeighborSimplex] = simplex;
        }
      }
      else
      {
        m_tri_nodes[-neighborSimplex][3] = simplex;
      }
    }

    for(unsigned int lN = 0 ; lN < sizeFace ; lN++)
    {
      unsigned int localNodeforTriangle  = localsInfoForTriangles[faceIter][lN];
      TSimplexID oppositeTriangle        = -oppositeTriangles[faceIter][lN];

      if(localNodeforTriangle != border)
      {
        TInt nodeB = face[(localNodeforTriangle + 1) % sizeFace];
        TInt nodeA = face[(localNodeforTriangle + 2) % sizeFace];

        TInt triangle = -addTriangle(nodeA, nodeB, nodeToConnect, false);
        triangles.push_back(triangle);
        (*BND_TRIANGLES)[-triangle] = trianglesIndices[faceIter][lN];

        m_tri_nodes[-triangle][3] = simplex;
        m_tet_adj[simplex][localNodeforTriangle] = triangle;
        std::vector<TInt> oppoNode = SimplicesTriangle(this, oppositeTriangle).getOtherNodeInSimplex(std::vector<TInt>{nodeA, nodeB});

        if(oppoNode.size() == 1)
        {
          unsigned int localnode = SimplicesTriangle(this, oppositeTriangle).getLocalNode(oppoNode.front());
          if(localnode < 3)
          {
              m_tri_adj[oppositeTriangle][localnode] = -triangle;
              m_tri_adj[-triangle][2] = oppositeTriangle;
          }
          else
          {
              std::cout << "localnode >= 3" << std::endl;
          }
        }
        else
        {
          std::cout << "oppoNode.size() != 1 --> "<<  oppoNode.size() << std::endl;
        }

        hash_edge[nodeA].push_back(-triangle);
        hash_edge[nodeB].push_back(-triangle);
      }
      else
      {
        triangles.push_back(border);
      }
    }
    faceIter++;
  }

  /////Triangle reconnection//////
  ////////////////////////////////
  for(auto const & pair : hash_edge)
  {
    TInt node = pair.first;
    std::vector<TInt> commonEdge{node, nodeToConnect};
    if(pair.second.size() == 2)
    {
      TSimplexID idxTriangleA = pair.second.front();
      TSimplexID idxTriangleB = pair.second.back();

      SimplicesTriangle triangleA = SimplicesTriangle(this, -idxTriangleA);
      SimplicesTriangle triangleB = SimplicesTriangle(this, -idxTriangleB);

      std::vector<TInt> nodesTriangleA = triangleA.getNodes();
      TInt nodeA = nodesTriangleA[0];
      TInt nodeB = nodesTriangleA[1];
      TInt localNodeA = (nodeA != node)? 0 : 1;


      std::vector<TInt> nodesTriangleB = triangleB.getNodes();
      nodeA = nodesTriangleB[0];
      nodeB = nodesTriangleB[1];
      TInt localNodeB = (nodeA != node)? 0 : 1;

      m_tri_adj[idxTriangleA][localNodeA] = idxTriangleB;
      m_tri_adj[idxTriangleB][localNodeB] = idxTriangleA;
    }
    else
    {
      std::cout << "pair.second.size() != 2" << std::endl;
    }
  }

  ////////////////////////////////
  std::vector<std::pair<TInt, TInt>> pairToDelete{};
  std::pair <std::multimap<std::pair<TInt, TInt>, TSimplexID>::iterator, std::multimap<std::pair<TInt, TInt>, TSimplexID>::iterator> pair_it;
  for(auto const & data : hash_face)
  {
    pair_it = hash_face.equal_range(data.first);
    int count = std::distance(pair_it.first, pair_it.second);
    if(count == 1)
    {
      pairToDelete.push_back(data.first);
    }
  }

  for(auto const & pairToDel : pairToDelete)
  {
    hash_face.erase(pairToDel);
  }

  std::vector<std::vector<TInt>> pointsToConnectCopy = pointsToConnect;
  std::vector<TSimplexID> tetIdToBuildInfoCopy = tetIdToBuildInfo;
  while(!hash_face.empty())
  {
    if(tetIdToBuildInfoCopy.size() == 0 && pointsToConnectCopy.size() == 0)
    {
      return;
    }

    TSimplexID simplexId   = tetIdToBuildInfoCopy.back();
    std::vector<TInt> face = pointsToConnectCopy.back();

    tetIdToBuildInfoCopy.pop_back();
    pointsToConnectCopy.pop_back();

    if(simplexId == border){continue;}
    if(face.size() == 0)
    {
      return;
    }

    for(unsigned int node = 0 ; node < face.size(); node++)
    {
      TInt node_0 = face[node];
      TInt node_1 = face[(node + 1) % face.size()];
      TInt node_2 = face[(node + 2) % face.size()];

      std::pair<TInt, TInt> pairNode = my_make(node_1, node_2);
      std::pair <std::multimap<std::pair<TInt, TInt>, TSimplexID>::iterator, std::multimap<std::pair<TInt, TInt>, TSimplexID>::iterator> pair_it;
      pair_it = hash_face.equal_range(pairNode);
      int count = std::distance(pair_it.first, pair_it.second);
      if(count == 2)
      {
        std::multimap<std::pair<TInt, TInt>, TSimplexID>::iterator it = pair_it.first;
        TSimplexID neighborSimplex = (it->second != simplexId) ? it->second : (++it)->second;
        TInt localNode0  = SimplicesCell(this, simplexId).getLocalNode(node_0);
        std::vector<TInt> intersectionFace{node_1, node_2, nodeToConnect};
        std::vector<TInt>&& nodesNeighborTet      = SimplicesCell(this, neighborSimplex).getOtherNodeInSimplex(intersectionFace);
        if(nodesNeighborTet.size() == 1)
        {
          TInt nodeLocalNeighborSimplex = SimplicesCell(this, neighborSimplex).getLocalNode(nodesNeighborTet.front());
          m_tet_adj[simplexId][localNode0] = neighborSimplex;
          m_tet_adj[neighborSimplex][nodeLocalNeighborSimplex] = simplexId;
          hash_face.erase(pair_it.first, pair_it.second);
        }
      }
    }
  }

  //correction of the cavity if plane tetraedron is was build by reinsertion of node for example
  for(unsigned int cellIdx = 0 ; cellIdx < tetIdToBuildInfo.size() ; cellIdx++)
  {
    const TSimplexID simplex      =  tetIdToBuildInfo[cellIdx];
    if(simplex != border)
    {
      const std::vector<TInt>nodes  =  std::vector<TInt>{m_tet_nodes[simplex][0], m_tet_nodes[simplex][1], m_tet_nodes[simplex][2]};
      for(unsigned int localNode = 0 ; localNode < nodes.size() ; localNode++)
      {
        TInt node = nodes[localNode];
        if(node == nodeToConnect)
        {
          TSimplexID oppositeSimplex  = m_tet_adj[simplex][3];
          TSimplexID neighborSimplex = m_tet_adj[simplex][localNode];

          if(neighborSimplex < 0){continue;}
          if(oppositeSimplex >= 0)
          {
            std::vector<TInt> othersNode0 = SimplicesCell(this, oppositeSimplex).getOtherNodeInSimplex(nodes);
            std::vector<TInt> othersNode1 = SimplicesCell(this, neighborSimplex).getOtherNodeInSimplex(nodes);

            if(othersNode0.size() == 1 && othersNode1.size() == 1)
            {

              TInt nodeLocalOppositeSimplex = SimplicesCell(this, oppositeSimplex).getLocalNode(othersNode0.front());
              m_tet_adj[oppositeSimplex][nodeLocalOppositeSimplex]   = neighborSimplex;

              TInt nodeLocalNeighborSimplex = SimplicesCell(this, neighborSimplex).getLocalNode(othersNode1.front());
              m_tet_adj[neighborSimplex][nodeLocalNeighborSimplex]   = oppositeSimplex;
            }
            else
            {
              std::cout << "othersNode0.size() != 1 || othersNode1.size() == 1" << std::endl;
            }
          }
          else
          {
            if(oppositeSimplex != border)
            {
              m_tri_nodes[-oppositeSimplex][3] = neighborSimplex;
              std::vector<TInt> othersNode = SimplicesCell(this, neighborSimplex).getOtherNodeInSimplex(nodes);
              if(othersNode.size() == 1)
              {
                TInt nodeLocalNeighborSimplex = SimplicesCell(this, neighborSimplex).getLocalNode(othersNode.front());
                m_tet_adj[neighborSimplex][nodeLocalNeighborSimplex]   = oppositeSimplex;
              }
              else
              {
                std::cout << "othersNode.size() != 1" << std::endl;
              }
            }
          }
          m_tet_ids.unselect(simplex);
          break;
        }
      }
    }
  }

  for(unsigned int triangleIdx = 0 ; triangleIdx < triangles.size() ; triangleIdx++)
  {
    const TSimplexID simplex       =  triangles[triangleIdx];
    if(simplex != border)
    {
      const std::vector<TInt> nodes{ m_tri_nodes[-simplex][0],  m_tri_nodes[-simplex][1]};

      for(unsigned int localNode = 0 ; localNode < 2 ; localNode++)
      {
        TInt node = nodes[localNode];
        if(node == nodeToConnect)
        {
          TSimplexID oppositeTriangle = m_tri_adj[-simplex][2];
          TSimplexID neighborTriangle = m_tri_adj[-simplex][localNode];

          std::vector<TInt> othersNode0 = SimplicesTriangle(this, oppositeTriangle).getOtherNodeInSimplex(nodes);
          std::vector<TInt> othersNode1 = SimplicesTriangle(this, neighborTriangle).getOtherNodeInSimplex(nodes);

          if(othersNode0.size() == 1 && othersNode1.size() == 1)
          {
            TInt nodeLocalOppositeTriangle = SimplicesTriangle(this, oppositeTriangle).getLocalNode(othersNode0.front());
            m_tri_adj[oppositeTriangle][nodeLocalOppositeTriangle] = neighborTriangle;

            TInt nodeLocalNeighborTriangle = SimplicesTriangle(this, neighborTriangle).getLocalNode(othersNode1.front());
            m_tri_adj[neighborTriangle][nodeLocalNeighborTriangle] = oppositeTriangle;
          }
          else
          {
            std::cout << "othersNode0.size() != 1 || othersNode1.size() == 1" << std::endl;
          }
          m_tri_ids.unselect(-simplex);
        }
      }
    }
  }
}
/******************************************************************************/
void SimplexMesh::fillBNDVariable()
{
  Variable<int>* BND_NODES_TOPO     = getVariable<int,SimplicesNode>("BND_NODES_TOPO");
  Variable<int>* BND_NODES_INDICES  = getVariable<int,SimplicesNode>("BND_NODES_INDICES");

  Variable<int>* BND_VERTEX_COLOR   = getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  Variable<int>* BND_CURVE_COLOR    = getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR  = getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  for(unsigned int node = 0 ; node < m_node_ids.capacity() ; node++)
  {
    if(m_node_ids[node] == 1)
    {
      if      ((*BND_VERTEX_COLOR)[node]  != 0){BND_NODES_TOPO->set(node, topo::CORNER);  BND_NODES_INDICES->set(node, (*BND_VERTEX_COLOR)[node]);}
      else if ((*BND_CURVE_COLOR) [node]  != 0){BND_NODES_TOPO->set(node, topo::RIDGE);   BND_NODES_INDICES->set(node, (*BND_CURVE_COLOR)[node]);}
      else if ((*BND_SURFACE_COLOR)[node] != 0){BND_NODES_TOPO->set(node, topo::SURFACE); BND_NODES_INDICES->set(node, (*BND_SURFACE_COLOR)[node]);}
    }
  }
}
/******************************************************************************/
void SimplexMesh::pointsLabeling(const std::vector<math::Point> &points, std::vector<int>& pointsLabeling, std::vector<int>& topoIndex, std::vector<TSimplexID>& nearNodes)
{
  pointsLabeling.reserve(points.size());
  topoIndex.reserve(points.size());
  nearNodes.reserve(points.size());

  int border = std::numeric_limits<int>::min();
  Variable<int>* BND_VERTEX_COLOR = getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  Variable<int>* BND_CURVE_COLOR = getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR = getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  int cpt = 0;
  for(auto const & point : points)
  {
    std::vector<TSimplexID> tetraContainingPt{};
    SimplicesNode::nodeNeighborInfo AnodeInfo;
    bool flag = checkSimplicesContenaing(point, tetraContainingPt);

    if(flag)
    {
      TSimplexID simplex = tetraContainingPt.front();
      bool flag = false;
      for(auto const currentSimplex : tetraContainingPt)
      {
        SimplicesCell cell = SimplicesCell(this, currentSimplex);
        for(unsigned int face = 0 ; face < 4 ; face++)
        {
          const TSimplexID oppositeSimplex = cell.oppositeTetraIdx(face);
          if(oppositeSimplex < 0)
          {
            simplex = currentSimplex;
            flag = true;
            break;
          }
        }
        if(flag)
        {
          break;
        }
      }

      //on trie les tétraèdre dans tetraContenaingPt de sorte a voir ceux qui sont aux bord du maillage en premier si ils existent
      /*std::sort(tetraContainingPt.begin(), tetraContainingPt.end(), [=](const TSimplexID simplexA, const TSimplexID simplexB)
      {
        std::cout << "simplexA ---> " << simplexA << std::endl;
        const SimplicesCell cellA = SimplicesCell(this, simplexA);
        bool hasBorderA = false;
        std::cout << cellA << std::endl;
        const SimplicesCell cellB = SimplicesCell(this, simplexA);
        std::cout << cellB << std::endl;

        for(unsigned int face = 0 ; face < 4 ; face++)
        {
          std::cout << "face -> " << face << std::endl;
          const TSimplexID oppositeSimplexA = cellA.oppositeTetraIdx(face);
          std::cout << "oppositeSimplexA -> " << oppositeSimplexA << std::endl;
          if(oppositeSimplexA < 0){hasBorderA = true; break;}
        }
        return hasBorderA;
      });*/

      //TSimplexID simplex = tetraContainingPt.front();
      nearNodes.push_back(simplex);
      SimplicesCell cell = SimplicesCell(this, simplex);
      double minLenght = std::numeric_limits<double>::max();
      unsigned int face = 0;

      std::vector<TInt> facesId{};
      double volTmp = std::numeric_limits<double>::max();
      double volTot = 0.0;
      for(unsigned int faceId = 0 ; faceId < 4 ; faceId++)
      {
          double bar = cell.signedBarycentricNormalized(faceId, point) ;
          volTot += bar;
          if(bar < 0.04)
          {
            facesId.push_back(faceId);
          }

          if(bar < volTmp)
          {
            face = faceId;
            volTmp = bar;
          }
      }


      unsigned int dimMax;
      unsigned int dimMin;
      TInt node0 ;
      TInt node1 ;
      TInt node2 ;
      unsigned int dim0;
      unsigned int dim1;
      unsigned int dim2;

      if(facesId.size() > 1)
      {
        std::vector<TInt> intersectionVec = cell.intersectionFaces(facesId);
        if(intersectionVec.size() == 1)
        {
          TInt node = intersectionVec.front();
          unsigned int dim = ((*BND_VERTEX_COLOR)[node] != 0)? CORNER : ((*BND_CURVE_COLOR)[node] != 0)?  RIDGE : ((*BND_SURFACE_COLOR)[node] != 0)? SURFACE : VOLUME;
          pointsLabeling.push_back(dim);
          if(dim == CORNER){topoIndex.push_back((*BND_VERTEX_COLOR)[node]);}
          else if(dim == RIDGE){topoIndex.push_back((*BND_CURVE_COLOR)[node]);}
          else if(dim == SURFACE){topoIndex.push_back((*BND_SURFACE_COLOR)[node]);}
          else {topoIndex.push_back(0);}
          continue;
        }
        else //intersectionVec.size() == 2
        {
          TInt node0 = intersectionVec.front();
          TInt node1 = intersectionVec.back();

          dim0 = ((*BND_VERTEX_COLOR)[node0] != 0)? CORNER : ((*BND_CURVE_COLOR)[node0] != 0)?  RIDGE : ((*BND_SURFACE_COLOR)[node0] != 0)? SURFACE : VOLUME;
          dim1 = ((*BND_VERTEX_COLOR)[node1] != 0)? CORNER : ((*BND_CURVE_COLOR)[node1] != 0)?  RIDGE : ((*BND_SURFACE_COLOR)[node1] != 0)? SURFACE : VOLUME;

          dimMax = std::max(dim0, dim1);
          dimMin = std::min(dim0, dim1);
        }
      }
      else
      {
        std::vector<TInt> orderedNodesFace = cell.getOrderedFace(face);

        node0 = orderedNodesFace[0];
        node1 = orderedNodesFace[1];
        node2 = orderedNodesFace[2];

        dim0 = ((*BND_VERTEX_COLOR)[node0] != 0)? CORNER : ((*BND_CURVE_COLOR)[node0] != 0)?  RIDGE : ((*BND_SURFACE_COLOR)[node0] != 0)? SURFACE : VOLUME;
        dim1 = ((*BND_VERTEX_COLOR)[node1] != 0)? CORNER : ((*BND_CURVE_COLOR)[node1] != 0)?  RIDGE : ((*BND_SURFACE_COLOR)[node1] != 0)? SURFACE : VOLUME;
        dim2 = ((*BND_VERTEX_COLOR)[node2] != 0)? CORNER : ((*BND_CURVE_COLOR)[node2] != 0)?  RIDGE : ((*BND_SURFACE_COLOR)[node2] != 0)? SURFACE : VOLUME;

        dimMax = std::max(dim0, std::max(dim1, dim2));
        dimMin = std::min(dim0, std::min(dim1, dim2));
      }


      if(dimMax == VOLUME)
      {
        pointsLabeling.push_back(VOLUME);
      }
      else
      {
        if(dimMax == SURFACE)
        {
          if(facesId.size() == 1)
          {
            if((dim0 == SURFACE && (dim0 == dim1)))
            {
              if((*BND_SURFACE_COLOR)[node0] != (*BND_SURFACE_COLOR)[node1] )
              {
                pointsLabeling.push_back(VOLUME);
              }
              else
              {
                pointsLabeling.push_back(SURFACE);
              }
            }
            else if((dim1 == SURFACE && (dim1 == dim2)))
            {
              if((*BND_SURFACE_COLOR)[node1] != (*BND_SURFACE_COLOR)[node2] )
              {
                pointsLabeling.push_back(VOLUME);
              }
              else
              {
                pointsLabeling.push_back(SURFACE);
              }
            }
            else if((dim2 == SURFACE && (dim2 == dim0)))
            {
              if((*BND_SURFACE_COLOR)[node2] != (*BND_SURFACE_COLOR)[node0] )
              {
                pointsLabeling.push_back(VOLUME);
              }
              else
              {
                pointsLabeling.push_back(SURFACE);
              }
            }
            else
            {
              pointsLabeling.push_back(SURFACE);
            }
          }
          else
          {
            pointsLabeling.push_back(SURFACE);
          }
        }
        else if(dimMax == RIDGE)
        {
          pointsLabeling.push_back(RIDGE);
        }
      }


      // build topo
      int dim = pointsLabeling.back();
      if(dim == VOLUME)
      {
        topoIndex.push_back(0);
      }
      else
      {
        if(dim == SURFACE)
        {
          TSimplexID oppositeSimplex = cell.oppositeTetraIdx(face);
          if(oppositeSimplex < 0 && oppositeSimplex != border)
          {
            Variable<int>* BND_TRIANGLES = getVariable<int,SimplicesTriangle>("BND_TRIANGLES");
            topoIndex.push_back((*BND_TRIANGLES)[-oppositeSimplex]);
          }
          else
          {
            struct faceNodeInfo
            {
              double lenght;
              int face;
            };

            std::vector<faceNodeInfo> facesInfo{};
            facesInfo.reserve(4);
            for(unsigned int faceId = 0 ; faceId < 4 ; faceId++)
            {
              faceNodeInfo info;
              double lenght;
              std::vector<TInt> orderedNodesFace = cell.getOrderedFace(faceId);
              const math::Point p0 = SimplicesNode(this, orderedNodesFace[0]).getCoords();
              const math::Point p1 = SimplicesNode(this, orderedNodesFace[1]).getCoords();
              const math::Point p2 = SimplicesNode(this, orderedNodesFace[2]).getCoords();
              math::Point projectedPoint;
              pointInTriangle(point, p0, p1, p2, lenght, projectedPoint);
              info.lenght = lenght;
              info.face   = faceId;
              facesInfo.push_back(info);
            }
            std::sort(facesInfo.begin(), facesInfo.end(), [&](const faceNodeInfo faceInfoA, const faceNodeInfo faceInfoB){
              return (faceInfoA.lenght < faceInfoB.lenght);
            });

            //face 0 && face1
            std::vector<TInt> orderedFaces = cell.getOrderedFace(facesInfo.front().face);
            TInt node0 = orderedFaces[0];
            TInt node1 = orderedFaces[1];
            TInt node2 = orderedFaces[2];
            TInt node = ((*BND_SURFACE_COLOR)[node0] != 0)? node0 : ((*BND_SURFACE_COLOR)[node1] != 0)?  node1 : ((*BND_SURFACE_COLOR)[node2] != 0)? SURFACE : node2;
            if(node != border)
            {
              topoIndex.push_back((*BND_SURFACE_COLOR)[node]);
            }
            else
            {
              //TODO exception
              std::cout << "node == border" << std::endl;
            }
          }
        }
        else
        {
          struct faceNodeInfo
          {
            double lenght;
            int face;
          };

          std::vector<faceNodeInfo> facesInfo{};
          facesInfo.reserve(4);
          for(unsigned int faceId = 0 ; faceId < 4 ; faceId++)
          {
            faceNodeInfo info;
            double lenght;
            std::vector<TInt> orderedNodesFace = cell.getOrderedFace(faceId);
            const math::Point p0 = SimplicesNode(this, orderedNodesFace[0]).getCoords();
            const math::Point p1 = SimplicesNode(this, orderedNodesFace[1]).getCoords();
            const math::Point p2 = SimplicesNode(this, orderedNodesFace[2]).getCoords();
            math::Point projectedPoint;
            pointInTriangle(point, p0, p1, p2, lenght, projectedPoint);
            info.lenght = lenght;
            info.face   = faceId;
            facesInfo.push_back(info);
          }
          std::sort(facesInfo.begin(), facesInfo.end(), [&](const faceNodeInfo faceInfoA, const faceNodeInfo faceInfoB){
            return (faceInfoA.lenght < faceInfoB.lenght);
          });

          if(dim == RIDGE)
          {
            //face0 & face1
            std::vector<TInt> face0 = cell.getOrderedFace(facesInfo[0].face);
            std::vector<TInt> face1 = cell.getOrderedFace(facesInfo[1].face);

            std::vector<TInt> intersectionNodes{};
            std::sort(face0.begin(), face0.end());
            std::sort(face1.begin(), face1.end());
            std::set_intersection(face0.begin(), face0.end(), face1.begin(), face1.end(), std::back_inserter(intersectionNodes));
            if(intersectionNodes.size() == 2)
            {
              TInt node0 = intersectionNodes[0];
              TInt node1 = intersectionNodes[1];
              unsigned int dim0 = ((*BND_VERTEX_COLOR)[node0] != 0)? CORNER : ((*BND_CURVE_COLOR)[node0] != 0)?  RIDGE : ((*BND_SURFACE_COLOR)[node0] != 0)? SURFACE : VOLUME;
              unsigned int dim1 = ((*BND_VERTEX_COLOR)[node1] != 0)? CORNER : ((*BND_CURVE_COLOR)[node1] != 0)?  RIDGE : ((*BND_SURFACE_COLOR)[node1] != 0)? SURFACE : VOLUME;

              unsigned int dimMax = std::max(dim0, dim1);
              if(dimMax == RIDGE)
              {
                TInt node = (dim0 >= dim1) ? node0 : node1;
                topoIndex.push_back((*BND_CURVE_COLOR)[node]);
              }
              else
              {
                std::cout << "dimMax != RIDGE" << std::endl;
              }
            }
            else
            {
              std::cout << "intersectionNodes.size() != 2" << std::endl;
            }
          }
          else // dim == CORNER
          {
              //face0 & face1 & face2
              std::vector<TInt> face0 = cell.getOrderedFace(facesInfo[0].face);
              std::vector<TInt> face1 = cell.getOrderedFace(facesInfo[1].face);
              std::vector<TInt> face2 = cell.getOrderedFace(facesInfo[2].face);
              std::vector<TInt> intersectionNodesFace01{};
              std::vector<TInt> intersectionNodes{};
              std::sort(face0.begin(), face0.end());
              std::sort(face1.begin(), face1.end());
              std::sort(face2.begin(), face2.end());
              std::set_intersection(face0.begin(), face0.end(), face1.begin(), face1.end(), std::back_inserter(intersectionNodesFace01));
              std::set_intersection(intersectionNodesFace01.begin(), intersectionNodesFace01.end(), face2.begin(), face2.end(), std::back_inserter(intersectionNodes));
              if(intersectionNodes.size() == 1)
              {
                TInt node = intersectionNodes.front();
                unsigned int dim = ((*BND_VERTEX_COLOR)[node] != 0)? CORNER : ((*BND_CURVE_COLOR)[node] != 0)?  RIDGE : ((*BND_SURFACE_COLOR)[node] != 0)? SURFACE : VOLUME;
                if(dim == CORNER)
                {
                  topoIndex.push_back((*BND_VERTEX_COLOR)[node]);
                }
                else
                {
                  //TODO exception
                  std::cout << "dim != CORNER" << std::endl;
                }
              }
              else
              {
                std::cout << "intersectionNodes.size() != 1" << std::endl;
              }
          }
        }
      }
    }
    else
    {
      //TODO exception
      std::cout << "tetraContenaingPt is empty " << std::endl;
      break;
    }
  }
}
/******************************************************************************/
void SimplexMesh::edgesRemove(const gmds::BitVector& nodeBitVector)
{
  gmds::ISimplexMeshIOService ioService(this);
  unsigned int cpt = 0;
  TInt border = std::numeric_limits<TInt>::min();
  Variable<int>* BND_VERTEX_COLOR = getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  Variable<int>* BND_CURVE_COLOR = getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR = getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
  std::vector<TInt> nodesNotConnected{};
  CriterionRAIS criterionRAIS(new VolumeCriterion());

  //Dimension node boundary vertex --> 0 / node on curve --> 1 / node on surface --> 2 / node on the mesh --> 3
  for(TInt node = 0; node < m_node_ids.capacity(); node++)
  {
    //std::cout << "node -> " << node << std::endl;
    //Collapsing Ni --> Nj
    int dim_Ni = 0;
    int index_Ni = 0;
    if(m_node_ids[node] != 0)
    {
      //check if the node is one of the boundaries of the mesh
      if(nodeBitVector[node] == 0)
      {
        if((*BND_VERTEX_COLOR)[node] != 0){dim_Ni = 0; index_Ni = (*BND_VERTEX_COLOR)[node];}
        else if((*BND_CURVE_COLOR)[node] != 0){dim_Ni = 1; index_Ni = (*BND_CURVE_COLOR)[node];}
        else if((*BND_SURFACE_COLOR)[node] != 0){dim_Ni = 2; index_Ni = (*BND_SURFACE_COLOR)[node];}
        else{dim_Ni = 4;}

        //if(dim_Ni == 4){continue;}

        if(dim_Ni != 0)
        {
          SimplicesNode simpliceNode(this, node);
          std::vector<TSimplexID>&& ball = simpliceNode.ballOf();

          if(ball.size() == 0)
          {
            continue;
          }

          std::set<TInt> nodes;
          for(auto const & tet : ball)
          {
            if(tet >= 0)
            {
              SimplicesCell simpliceCell(this, tet);
              std::vector<TInt>&& cellNodes = simpliceCell.getNodes();
              for(auto const cellNode : cellNodes)
              {
                if(m_node_ids[cellNode] != 0 && cellNode != node)
                {
                    nodes.insert(cellNode);
                }
              }
            }
          }

          CriterionRAIS criterionRAIS(new VolumeCriterion());
          //2 choix --> on trie cette fois si les longueur e ne gardant que les node donnant nodeToInsert.checkStar(ball, criterionRAIS) vrai
          auto cmp = [](const double& lenghtNodeA, const double& lenghtNodeB){return lenghtNodeA < lenghtNodeB;};
          std::map<double, dataNode, decltype(cmp)> adjacentNodesInfo{cmp};
          std::map<double, dataNode, decltype(cmp)> adjacentNodesInfoTmp{cmp};
          std::vector< std::pair<double, dataNode>> concattedAdjNodes;
          if(nodes.size() != 0)
          {
            const math::Point&& pt = simpliceNode.getCoords();
            for(auto const cellNode : nodes)
            {
              SimplicesNode nodeTet(this, cellNode);
              math::Vector3d vec = nodeTet.getCoords() - pt;
              double lenght = vec.norm();
              dataNode nodeInfo;
              nodeInfo.node = cellNode;
              nodeInfo.lenght = lenght;
              if      ((*BND_VERTEX_COLOR)[cellNode] != 0){continue;/*nodeInfo.dim_Nj = 0; nodeInfo.index_Nj = (*BND_VERTEX_COLOR)[cellNode];*/}
              else if ((*BND_CURVE_COLOR)[cellNode] != 0){nodeInfo.dim_Nj = 1; nodeInfo.index_Nj = (*BND_CURVE_COLOR)[cellNode];}
              else if ((*BND_SURFACE_COLOR)[cellNode] != 0){nodeInfo.dim_Nj = 2; nodeInfo.index_Nj = (*BND_SURFACE_COLOR)[cellNode];}
              else    {nodeInfo.dim_Nj = 4;}
              //if(nodeBitVector[cellNode] == 1)
              {
                std::vector<TSimplexID> ballOnlyTet{};
                std::copy_if(ball.begin(), ball.end(),std::back_inserter(ballOnlyTet),[&](TInt simplex){
                  return (simplex >= 0);
                });

                //if(nodeTet.checkStar(ballOnlyTet, criterionRAIS))
                {
                  adjacentNodesInfo.insert( std::pair<double, dataNode>(lenght, nodeInfo));
                }
              }
              /*else
              {
                adjacentNodesInfoTmp.insert( std::pair<double, dataNode>(lenght, nodeInfo));
              }*/
            }
            //concatenation of adjacentNodesInfo & adjacentNodesInfoTMP
            std::copy(adjacentNodesInfo.begin(), adjacentNodesInfo.end(), std::back_inserter(concattedAdjNodes));
            //std::copy(adjacentNodesInfoTmp.begin(), adjacentNodesInfoTmp.end(), std::back_inserter(concattedAdjNodes));
          }
          bool cavityFound  = false;
          for(auto const & dataInfo : concattedAdjNodes)
          {
            dataNode data = dataInfo.second;
            if(data.dim_Nj <= dim_Ni)
            {
              if(data.dim_Nj == dim_Ni)
              {
                if(data.index_Nj != index_Ni)
                {
                  continue;
                }
              }
              SimplicesNode nodeToInsert(this, data.node);
              bool status = false;
              //std::cout << "NODE TO INSERT : " << data.node << " | dimension : " << data.dim_Nj << " | index : " << data.index_Nj<<  std::endl;
              //std::cout << "FROM NODE      : " << node << " | dimension : " << dim_Ni << " | index : " << index_Ni <<  std::endl;
              if(dim_Ni == 4 && (data.dim_Nj == 1 || data.dim_Nj == 2)){break;}
              if(dim_Ni == 2 && (data.dim_Nj == 1 || data.dim_Nj == 0)){break;}
              std::vector<TSimplexID> deletedSimplex{};
              PointInsertion(this, nodeToInsert, criterionRAIS, status, ball, nodeBitVector, deletedSimplex);
              if(status)
              {

                //std::cout << "iteration --->   " << cpt << std::endl;
                //std::vector<TSimplexID> v =  SimplicesNode(this, 584).ballOf();
                //std::cout << "SIZE BALL OF 584---> " << v.size() << std::endl;

                //if(cpt > 540)
                {
                  //gmds::VTKWriter vtkWriter(&ioService);
                  //vtkWriter.setCellOptions(gmds::N|gmds::F|gmds::R);
                  //vtkWriter.setDataOptions(gmds::N|gmds::F|gmds::R);
                  //vtkWriter.write("ModeleCAD5_iteration_" + std::to_string(cpt) + ".vtk");
                  //vtkWriter.write("ModeleCAD5_iteration_BIS_" + std::to_string(cpt) + ".vtk");
                  //vtkWriter.write("ModeleCAD5_iteration_" + std::to_string(cpt) + ".vtk");
                  //vtkWriter.write("ModeleCAD5_CAV_" + std::to_string(cpt) + ".vtk");
                  /*SimplicesTriangle triangleTest = SimplicesTriangle(this,6563);
                  std::cout << "triangleTest ---> " << triangleTest << std::endl;
                  if(m_tri_adj[6563][2] == border)
                  {
                    std::cout << "triangleTest -> " << triangleTest << std::endl;
                    std::cout << "ITERATION -> " << cpt << std::endl;
                    return;
                  }*/
                  //if(cpt == 549){return;}
                }
                //std::cout << std::endl;
                //std::cout << std::endl;
                cpt++;

                break;
              }
            }
          }
        }
      }
    }
  }
}
/******************************************************************************/
void SimplexMesh::buildEdges(const std::vector<std::vector<TInt>>& AEdges, const gmds::BitVector& nodeBitVector)
{
  const TInt border = std::numeric_limits<TInt>::min();
  const gmds::BitVector& bitNodes = getBitVectorNodes();
  unsigned int cpt_EdgeAlreadyBuild = 0;
  double cpt_EdgeBuilt = 0.0;
  double sizeEdge = AEdges.size() / 2.0;
  for(auto const & edge : AEdges)
  {
    if(edge.size() == 2)
    {
      TInt nodeAidx = edge.front();
      TInt nodeBidx = edge.back();
      std::cout << "bitNodes.capacity() --> " << bitNodes.capacity() << std::endl;
      std::cout << "nodeAidx --> " << nodeAidx << std::endl;
      std::cout << "nodeBidx --> " << nodeBidx << std::endl;

      if(nodeAidx != border && nodeBidx != border)
      {
        if(bitNodes[nodeAidx] != 0 && bitNodes[nodeBidx] != 0)
        {
          SimplicesNode nodeA = SimplicesNode(this, nodeAidx);
          SimplicesNode nodeB = SimplicesNode(this, nodeBidx);
          if(nodeA.shell(nodeB).size() == 0)
          {

            std::cout << nodeA << std::endl;
            std::cout << nodeB << std::endl;
            std::vector<TSimplexID> cavityTmp{};
            std::vector<TSimplexID> ballA = nodeA.ballOf();
            math::Vector3d vec = nodeB.getCoords() - nodeA.getCoords();
            math::Vector3d u = math::Vector3d(vec.Y(), - vec.X(), 0.0);
            math::Vector3d v = vec.cross(u);
            u.normalize(); v.normalize();

            math::Point pt1 = nodeA.getCoords();
            math::Point pt2 = pt1 + u;
            math::Point pt3 = pt1 + v;
            std::vector<TInt> nodes{};

            std::cout << "ballA.size() -> " << ballA.size() << std::endl;
            if(ballA.size() == 0){continue;}
            for(auto const & simplex : ballA)
            {
              std::cout << "simplex -> " << simplex << std::endl;
              if(simplex >= 0){nodes=  SimplicesCell(this, simplex).getNodes();}
              else if(simplex < 0 && simplex != border){ nodes = SimplicesTriangle(this, -simplex).getNodes();}

              for(auto const node : nodes)
              {
                const math::Point pt0 = SimplicesNode(this,node).getCoords();
                math::Orientation::Sign sign = math::Orientation::orient3d(pt0, pt1, pt2, pt3);
                if(sign == 1)
                {
                  cavityTmp.push_back(simplex);
                  break;
                }
              }
            }
            CriterionRAIS criterionRAIS(new VolumeCriterion());
            bool status = false;
            std::vector<TSimplexID> deletedSimplex{};
            std::vector<TSimplexID> cavity{};
            if(cavityTmp.size() != 0)
            {
                for(auto const & simplex : cavityTmp)
                {
                  if(simplex >= 0)
                  {
                    cavity.push_back(cavityTmp.front());
                    break;
                  }
                }
                std::cout << "PointInsertion" << std::endl;
                PointInsertion(this, nodeB, criterionRAIS, status, cavity, nodeBitVector, deletedSimplex);
            }
            if(!status)
            {
              std::cout << "edge [" << nodeAidx << " ; " << nodeBidx << "] was not built " << std::endl;
            }
            else
            {
              cpt_EdgeBuilt = cpt_EdgeBuilt + 1.0;;
              std::cout << "edge [" << nodeAidx << " ; " << nodeBidx << "] was built " << std::endl;
            }
            std::cout << std::endl;
            std::cout << std::endl;
          }
          else
          {
            cpt_EdgeAlreadyBuild = cpt_EdgeAlreadyBuild + 1.0;
            std::cout << "edge [" << nodeAidx << " ; " << nodeBidx << "] is already built ! " << std::endl;
          }
        }
      }
    }
  }

  std::cout << "EDGE BUILT = " << cpt_EdgeBuilt / sizeEdge << std::endl;;
  std::cout << "EDGE ALREADY BUILT = " << cpt_EdgeAlreadyBuild / sizeEdge << std::endl;
}
/******************************************************************************/
bool SimplexMesh::isHexaEdgeBuild(const std::vector<std::vector<TInt>>& ANodesFaces)
{
  for(auto const & face : ANodesFaces)
  {
    if(!isFaceBuild(face))
    {
      return false;
    }
  }
  return true;
}
/******************************************************************************/
std::vector<TSimplexID> SimplexMesh::hex2tet(const std::vector<TInt>& ANodesHex)
{
  //use is hexa build before..
  //intialisation of the tet..
  std::vector<TSimplexID> res{};
  if(ANodesHex.size() == 8)
  {
    TInt border = std::numeric_limits<TInt>::min();
    TSimplexID initTet = border;
    std::vector<TSimplexID> && ball0Tmp = SimplicesNode(this, ANodesHex[0]).ballOf();
    std::vector<TSimplexID> && ball1Tmp = SimplicesNode(this, ANodesHex[1]).ballOf();
    std::vector<TSimplexID> && ball2Tmp = SimplicesNode(this, ANodesHex[2]).ballOf();
    std::vector<TSimplexID> && ball3Tmp = SimplicesNode(this, ANodesHex[3]).ballOf();

    std::vector<TSimplexID> ball0; std::vector<TSimplexID> ball2;
    std::vector<TSimplexID> ball1; std::vector<TSimplexID> ball3;

    std::copy_if(ball0Tmp.begin(), ball0Tmp.end(), std::back_inserter(ball0), [&](TSimplexID simplex){return (simplex >= 0);});
    std::copy_if(ball1Tmp.begin(), ball1Tmp.end(), std::back_inserter(ball1), [&](TSimplexID simplex){return (simplex >= 0);});
    std::copy_if(ball2Tmp.begin(), ball2Tmp.end(), std::back_inserter(ball2), [&](TSimplexID simplex){return (simplex >= 0);});
    std::copy_if(ball3Tmp.begin(), ball3Tmp.end(), std::back_inserter(ball3), [&](TSimplexID simplex){return (simplex >= 0);});

    std::vector<std::vector<TSimplexID>> balls{{ball0}, {ball1}, {ball2}};
    std::vector<TSimplexID>&& inter  = intersectionSimplex(balls);
    if(!(inter.size()))
    {
      balls.pop_back();
      balls.push_back(ball3);
      inter  = intersectionSimplex(balls);
      if(!(inter.size()))
      {
          std::cout << "!(inter.size()" << std::endl;
          return res;
      }
    }
    //test with hexa's local node 4
    std::vector<TSimplexID> ball4 = SimplicesNode(this, ANodesHex[4]).ballOf();
    balls.push_back(ball4) ;
    inter = intersectionSimplex(balls);
    if(!inter.size())
    {
      balls.pop_back();
      std::vector<TSimplexID> ball5 = SimplicesNode(this, ANodesHex[5]).ballOf();
      balls.push_back(ball5) ;
      inter = intersectionSimplex(balls);
      if(!inter.size())
      {
        balls.pop_back();
        std::vector<TSimplexID> ball6 = SimplicesNode(this, ANodesHex[6]).ballOf();
        balls.push_back(ball6) ;
        inter = intersectionSimplex(balls);
        if(!inter.size())
        {
          balls.pop_back();
          std::vector<TSimplexID> ball7 = SimplicesNode(this, ANodesHex[7]).ballOf();
          balls.push_back(ball7) ;
          inter = intersectionSimplex(balls);
          if(!inter.size())
          {
            std::cout << "inter.size() for second face is null" << std::endl;
            return res;
          }
        }
      }
    }

    if(inter.size() != 1 )
    {
      std::cout << "inter.size() != 1 " << std::endl;
      return res;
    }

    initTet = inter.front();
    std::vector<TSimplexID> to_do{initTet};
    gmds::BitVector markedTet(getBitVectorTet().capacity());

    //create a cube boxe to look if a node is inside this box
    const BitVector& nodesVector = getBitVectorNodes();
    double epsilon = 10E-5;
    double min = std::numeric_limits<double>::min();
    double max = std::numeric_limits<double>::max();

    std::vector<double> nodes{max, min, max, min, max, min}; // --> xmin, xmax, ymin, ymax, zmin, zmax

    for(auto const nodeIds : ANodesHex)
    {
      if(nodesVector[nodeIds] != 0)
      {
        math::Point pt = SimplicesNode(this, nodeIds).getCoords();
        double x = pt.X();
        double y = pt.Y();
        double z = pt.Z();

        nodes[0] = (x < nodes[0])? x : nodes[0];
        nodes[1] = (x > nodes[1])? x : nodes[1];

        nodes[2] = (y < nodes[2])? y : nodes[2];
        nodes[3] = (y > nodes[3])? y : nodes[3];

        nodes[4] = (z < nodes[4])? z : nodes[4];
        nodes[5] = (z > nodes[5])? z : nodes[5];
      }
    }

    double xmin = nodes[0] - epsilon; double xmax = nodes[1] + epsilon;
    double ymin = nodes[2] - epsilon; double ymax = nodes[3] + epsilon;
    double zmin = nodes[4] - epsilon; double zmax = nodes[5] + epsilon;



    while(!to_do.empty())
    {
      TSimplexID t = to_do.back();
      res.push_back(t);
      markedTet.assign(t);
      to_do.pop_back();
      std::vector<TSimplexID> adjSimplex = SimplicesCell(this, t).adjacentTetra();
      for(auto const simplex : adjSimplex)
      {
        if(simplex >= 0)
        {
          if(markedTet[simplex] == 0)
          {
            markedTet.assign(simplex);
            std::vector<TInt> nodes = SimplicesCell(this, simplex).getNodes();
            unsigned int cpt = 0;
            for(unsigned int idx = 0 ; idx < nodes.size() ; idx++)
            {
              TInt node = nodes[idx];
              if(!(node == ANodesHex[0] || node == ANodesHex[1] || node == ANodesHex[2] || node == ANodesHex[3] ||
                node == ANodesHex[4] || node == ANodesHex[5] || node == ANodesHex[6] || node == ANodesHex[7]))
                {
                  math::Point pt = SimplicesNode(this, node).getCoords();
                  double x = pt.X() ; double y = pt.Y() ; double z = pt.Z() ;
                  if(x >= xmin && x <= xmax)
                  {
                    if(y >= ymin && y <= ymax)
                    {
                      if(z >= zmin && z <= zmax)
                      {
                        to_do.push_back(simplex);
                        break;
                      }
                    }
                  }
                }
                else
                {
                  cpt++;
                }
            }

            if(cpt == 4)
            {
              to_do.push_back(simplex);
            }
          }
        }
      }
    }
  }
  else
  {
    std::cout << "ANodesHex.size() != 4" << std::endl;
  }

  return res;
}
/******************************************************************************/
bool SimplexMesh::pointInTriangle(const math::Point& query_point,
                     const math::Point& triangle_vertex_0,
                     const math::Point& triangle_vertex_1,
                     const math::Point& triangle_vertex_2,
                     double& distance,
                     math::Point& projectedPoint)
{
  double epsilon = 10E-4;
    // u=P2−P1
    math::Vector3d u = triangle_vertex_1 - triangle_vertex_0;
    // v=P3−P1
    math::Vector3d v = triangle_vertex_2 - triangle_vertex_0;
    // n=u×v
    math::Vector3d n = u.cross(v);
    // w=P−P1
    math::Vector3d w = query_point - triangle_vertex_0;
    // Barycentric coordinates of the projection P′of P onto T:
    // γ=[(u×w)⋅n]/n²
    double gamma = u.cross(w).dot(n) / n.dot(n);
    // β=[(w×v)⋅n]/n²
    double beta = w.cross(v).dot(n) / n.dot(n);
    double alpha = 1 - gamma - beta;

    projectedPoint = alpha*triangle_vertex_0 + beta*triangle_vertex_1 + gamma*triangle_vertex_2;
    // distance from the face
    math::Vector3d normal = n;
    normal.normalize();
    distance = w.dot(normal);
    // The point P′ lies inside T if:
    return ((-epsilon <= alpha) && (alpha <= 1.0 + epsilon) &&
            (-epsilon <= beta)  && (beta  <= 1.0 + epsilon) &&
            (-epsilon <= gamma) && (gamma <= 1.0 + epsilon));
}
/******************************************************************************/
std::vector<TSimplexID> SimplexMesh::cavityIntersectedByEdge(SimplicesNode& simpliceNodeA, SimplicesNode& simpliceNodeB)
{
  std::vector<TSimplexID> v{};
  return v;
}
/****************************************************************************/
std::vector<VariableItf*> SimplexMesh::getAllVariables(ECellType AType) const
{
	switch (AType){
	case GMDS_NODE :
		return m_node_variable_manager->getAllVariables();
		break;

  case GMDS_TETRA :
    return m_tet_variable_manager->getAllVariables();
    break;

  case GMDS_TRIANGLE :
    return m_tri_variable_manager->getAllVariables();
    break;

	default:
		throw GMDSException("Unmanaged type of value -> impossible to access to a variable");
	}
}
