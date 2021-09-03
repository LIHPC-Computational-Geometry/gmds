/*---------------------------------------------------------------------------*/
#include <gmds/sofiane/SimplexMesh.h>
#include <gmds/sofiane/ICriterion.h>
#include <gmds/sofiane/PointInsertion.h>
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
  m_node_ids = gmds::BitVector(initSize);
  m_tet_ids = gmds::BitVector(initSize);
  m_tri_ids = gmds::BitVector(initSize);

  //init the tet tri and node size with the same size of ids ...
  m_coords.resize(initSize);
  m_tet_nodes.resize(initSize);
  m_tri_nodes.resize(initSize);
  TInt idx0 = m_tri_ids.getFreeIndex();
  m_tri_ids.assign(idx0);

  m_node_variable_manager = new VariableManager();
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
  m_triangle_variable_manager = simplexMesh.m_triangle_variable_manager;
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
    m_triangle_variable_manager = simplexMesh.m_triangle_variable_manager;
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
  m_triangle_variable_manager = simplexMesh.m_triangle_variable_manager;

  /*Clear all the simplexMesh instance data*/
  simplexMesh.clear();


  /*creer une fonction qui clear le vector dans variable manager ?*/
  /*simplexMesh.m_node_variable_manager;
  simplexMesh.m_tet_variable_manager;
  simplexMesh.m_triangle_variable_manager;
  */
}
/*****************************************************************************/
SimplexMesh::SimplexMesh(const StructuredGrid& structuredGrid)
{
  clear();
  /*je ne veux pas mettre de auto car on ne sait plus trop ce qu'on manipule ... a voir..*/
  const std::unordered_map<math::Vector3i, math::Point, StructuredGrid::hash_function> & grid = structuredGrid.getUnorderedMap();
  /*Clear all the simplexMesh instance data*/
  for(auto const & pair : grid)
  {
      /*TODO faire un check pour voir si il ya pas de doublon (peut etre mettre ce check directement dans structuredgrid)*/
      addNode(pair.second);
  }
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
  //a voir si c'est cette objet qui s'occupe de detruire ces intance ou non*/
  //ou rajouter un if == nullptr ...
  for(int i = 0; i< m_tet_nodes.size(); i++)
  {
    if(m_tet_nodes[i] != nullptr && m_tet_ids[i] != 0)
    {
        delete[] m_tet_nodes[i];
    }

  }

  for(int i = 0; i< m_tet_adj.size(); i++)
  {
    if(m_tet_nodes[i] != nullptr && m_tet_ids[i] != 0)
    {
        delete[] m_tet_adj[i];
    }
  }

  for(int i = 0; i< m_tri_nodes.size(); i++)
  {
    if(m_tri_nodes[i] != nullptr && m_tri_ids[i] != 0)
    {
        delete[] m_tri_nodes[i];
    }
  }

  for(int i = 0; i< m_tri_adj.size(); i++)
  {
    if(m_tri_adj[i] != nullptr && m_tri_ids[i] != 0)
    {
        delete[] m_tri_adj[i];
    }
  }



  if(m_node_variable_manager != nullptr)
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
  }

}
/*---------------------------------------------------------------------------*/
SimplexMesh::SimplexMesh(gmds::BitVector node_ids,
            std::vector<gmds::math::Point> coords,
            std::vector<TSimplexID> base,
            gmds::BitVector tet_ids,
            std::vector<TSimplexID* > tet_nodes,
            std::vector<TSimplexID* > tet_adj,
            gmds::BitVector tri_ids,
            std::vector<TSimplexID* > tri_nodes,
            std::vector<TSimplexID* > tri_adj)
{
  m_node_ids    = node_ids;
  m_coords      = coords;
  m_base        = base;

  m_tet_ids     = tet_ids;
  m_tet_nodes   = tet_nodes;
  m_tet_adj     = tet_adj;

  m_tri_ids     = tri_ids;
  m_tri_nodes   = tri_nodes;
  m_tri_adj     = tri_adj;
}

/*---------------------------------------------------------------------------*/
SimplexMesh::SimplexMesh(
            std::vector<gmds::math::Point> coords,
            std::vector<TSimplexID* > tet_nodes,
            std::vector<TSimplexID* > tet_adj,
            std::vector<TSimplexID* > tri_nodes,
            std::vector<TSimplexID* > tri_adj)
{
  m_coords      = coords;

  m_tet_nodes   = tet_nodes;
  m_tet_adj     = tet_adj;

  m_node_ids = gmds::BitVector(m_coords.size());
  m_node_ids.fillAll();
  m_tet_ids = gmds::BitVector(m_tet_nodes.size());
  m_tet_ids.fillAll();
  buildBaseVector();

  m_tri_ids = gmds::BitVector(m_tri_nodes.size());
  m_tri_ids.fillAll();
  m_tri_nodes   = tri_nodes;
  m_tri_adj     = tri_adj;

}

/*---------------------------------------------------------------------------*/
SimplexMesh::SimplexMesh(
            std::vector<gmds::math::Point> coords,
            std::vector<TSimplexID* > tet_nodes,
            std::vector<TSimplexID* > tri_nodes,
            const bool flag)
{
  m_coords      = coords;

  m_tet_nodes   = tet_nodes;

  m_node_ids = gmds::BitVector(m_coords.size());
  m_node_ids.fillAll();
  m_tet_ids = gmds::BitVector(m_tet_nodes.size());
  m_tet_ids.fillAll();
  buildBaseVector();

  m_tri_nodes.resize(tri_nodes.size() + 1);
  m_tri_nodes[0] = 0;
  m_tri_ids = gmds::BitVector(m_tri_nodes.size());
  m_tri_ids.fillAll();

  for(TInt i = 0; i < tri_nodes.size(); i++)
  {
      m_tri_nodes[i + 1]   = tri_nodes[i];
  }

  buildAdjTetVector();
  buildOppFacesVector();

  if(flag)
  {
    reorientAllTetra();
  }

}


/*---------------------------------------------------------------------------*/
/*use this function evry time update is done for the tet_node or the node...*/
void SimplexMesh::buildBaseVector()
{
  int errorId = std::numeric_limits<int>::min();
  //m_base.clear();
  bool noTet = true;
  for(TInt iterTet = 0 ; iterTet < m_tet_ids.capacity() ; iterTet++)
  {
    if(m_tet_ids[iterTet] == 1)
    {
      noTet = false;
      break;
    }
  }

  if(!noTet)
  {

    m_base.resize(m_coords.size());
    for(auto & value : m_base)
    {
      value = errorId;
    }
    std::vector<std::set<TSimplexID>> IndexNodes2Tet;
    std::map<TInt, std::vector<TSimplexID>> IndexNodes2TetMap;
    IndexNodes2Tet.resize(m_base.size());
    std::map<TSimplexID, std::vector<int>> Tet2IndexNodes;
    size_t iter = 0;


    /*fill the current mapping between the tet & index node*/
    for(auto const & indexTet : m_tet_nodes)
    {
        if(m_tet_ids[iter] != 0) /*continuer a utiliser cette condition ou le bitvector de tet_node ...*/
        {
            std::vector<int> tet{*indexTet, *(indexTet + 1), *(indexTet + 2), *(indexTet + 3)};
            Tet2IndexNodes.insert(std::pair<TSimplexID, std::vector<int>>(iter, tet));
        }
        iter++;
    }
    iter = 0;

    /*build the base vector with the previous mapping*/
    for(auto const & tet2node : Tet2IndexNodes)
    {
      for(auto const & indexNode : tet2node.second)
      {
          if(indexNode != errorId)
          {
              IndexNodes2Tet[indexNode].insert(tet2node.first);
          }
      }
    }



    /*pushing the IndexNodes2Tet set to the IndexNodes2TetMap (vector have a better performance..)*/
    for(int node = 0; node <IndexNodes2Tet.size(); node++)
    {
      if(m_node_ids[node] != 0)
      {
          std::vector<TSimplexID> vTet;
          vTet.insert(vTet.begin(), IndexNodes2Tet[node].begin(), IndexNodes2Tet[node].end());
          IndexNodes2TetMap[node] = vTet;
      }
    }

    for(auto const & val : IndexNodes2TetMap)
    {
      if(m_node_ids[val.first] != 0)
      {
          std::vector<TSimplexID> Tets = val.second;
          if(Tets.size() != 0)
          {
            for(auto const & valComp : IndexNodes2TetMap)
            {
              if(val.first != valComp.first)
              {
                const std::vector<TSimplexID> & Tets2Comp = valComp.second;
                setComparator(Tets, Tets2Comp);
                //if(Tets.size() == 1)
                {
                  m_base[val.first] = *(Tets.begin());
                  break;
                }
                /*else
                {
                  m_base[val.first] = *(Tets.begin());
                }*/
              }
            }
          }
      }
    }
  }
}
/*****************************************************************************/
void SimplexMesh::buildAdjTetVector()
{
  m_tet_adj.clear();
  m_tet_adj.resize(m_tet_nodes.size());
  TSimplexID AdjTet_or_tri  = std::numeric_limits<int>::min();


  for(size_t TetIndex = 0; TetIndex < m_tet_nodes.size(); TetIndex++)
  {
    if(m_tet_ids[TetIndex] != 0)
    {
      TSimplexID* adjTmp = new TSimplexID[4]{AdjTet_or_tri, AdjTet_or_tri, AdjTet_or_tri, AdjTet_or_tri};

      for(size_t nodeIndex = 0; nodeIndex < 4; nodeIndex++)
      {

        AdjTet_or_tri           = adjacentTet_or_tri(TetIndex, m_tet_nodes[TetIndex][nodeIndex]);
        adjTmp[nodeIndex]       = AdjTet_or_tri;
      }
      m_tet_adj[TetIndex]       = adjTmp;
    }
  }
}

/*****************************************************************************/
TSimplexID SimplexMesh::adjacentTet_or_tri(const TInt currentTet, const TInt ANode)
{
  const size_t TetNodeIndexSize = 4;
  TSimplexID* currentTetNodes = new TSimplexID[3]{-1, -1, -1};
  size_t iter = 0;
  for(size_t nodeLocal = 0; nodeLocal < TetNodeIndexSize; nodeLocal++)
  {
    if(ANode != m_tet_nodes[currentTet][nodeLocal])
    {
      currentTetNodes[iter] = m_tet_nodes[currentTet][nodeLocal];
      iter++;
    }
  }


  bool flag          = false;
  TSimplexID AdjTet_or_tri  = std::numeric_limits<int>::min();

  /*Check if a triangle is adjacent to the curerent tet.*/
  for(size_t triComp = 1 /*start at 1*/; triComp < m_tri_nodes.size(); triComp++)
  {
    if(m_tri_ids[triComp] != 0)
    {
      const TSimplexID* arrayComp    = m_tri_nodes[triComp];

      flag = containNodes<3,3>(currentTetNodes, arrayComp);
      if(flag)
      {
        /*the adjacent triangle of a tet must hyave index < 0 in the adjTet vector...to recognize if it's a triangle or not*/
        AdjTet_or_tri = -triComp;
        return AdjTet_or_tri;
      }
    }
  }


  for(size_t TetComp = 0; TetComp < m_tet_nodes.size(); TetComp++)
  {
    if(m_tet_ids[TetComp] != 0)
    {
      if(TetComp != currentTet)
      {
        const TSimplexID* arrayComp    = m_tet_nodes[TetComp];

        flag = containNodes<3,4>(currentTetNodes, arrayComp);
        if(flag)
        {
          AdjTet_or_tri = TetComp;
          return AdjTet_or_tri;
        }

      }
    }
  }

  delete [] currentTetNodes;
  currentTetNodes = nullptr;

  return AdjTet_or_tri;
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
TInt SimplexMesh::getOppositeCell(const TSimplexID ANodeID, const TSimplexID ATetID)
{
  const int errorId               = std::numeric_limits<int>::min();
  size_t localIndex               = errorId;
  const size_t indexBufferSize    =  4;
  const TSimplexID* indexBuffer   = m_tet_nodes[ATetID];


  localIndex = SimplicesCell(this, ATetID).getLocalNode(ANodeID);
  //si la valeur retournée est < 0 --> c'est la valeur negative de l'index d'un triangle dans m_tri_nodes..
  if(localIndex == errorId)
  {
    //si il n'existe ni triangle ni tetra opposé cherché
    return errorId;
  }
  else
  {
    if(m_tet_adj[ATetID] != nullptr)
    {
        return m_tet_adj[ATetID][localIndex];
    }
    else
    {
      return errorId;
    }

  }

}

/*---------------------------------------------------------------------------*/
TInt SimplexMesh::addNode(const math::Point& pt)
{
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
/*****************************************************************************/
TInt SimplexMesh::addNode(const Point&& pt)
{
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
bool SimplexMesh::deleteNode(const SimplicesNode& simpliceNode)
{
  return deleteNode(simpliceNode.getGlobalNode());
}
/*****************************************************************************/
bool SimplexMesh::deleteNode(const TInt indexNode)
{
  bool flag = false;
  if(indexNode > m_node_ids.capacity() || indexNode < 0)
  {
    /*TODO mettre une exeption a la place d'un bool*/
    /*Node index  <0 || > m_node_ids.capacity()*/
    flag = false;
  }
  else
  {
    const std::vector<TSimplexID>&& ballOfNode = SimplicesNode(this, indexNode).ballOf();
    const std::vector<TSimplexID>&& triangleAdjToNode = SimplicesNode(this, indexNode).linksTri();

    /*Dans un premier temps on delete les tetra liée au node indexNode*/
    for(const auto& simplex : ballOfNode)
    {
      if(m_tet_ids[simplex] != 0)
      {
          deleteTetra(simplex);
      }
    }

    /*Dans un second temps on supprime les triangles liée au node indexNode*/
    for(const auto & triangle : triangleAdjToNode)
    {
      if(m_tri_ids[triangle] != 0)
      {
        deleteTriangle(-triangle);
      }
    }

    m_node_ids.unselect(indexNode);
    m_node_variable_manager->removeEntry(indexNode);
    flag = true;
  }
  return flag;
}
/*---------------------------------------------------------------------------*/
TSimplexID SimplexMesh::addTriangle(const simplicesNode::SimplicesNode& ANode0,
                       const simplicesNode::SimplicesNode& ANode1,
                       const simplicesNode::SimplicesNode& ANode2)
{
  int errorId = std::numeric_limits<int>::min();
  TSimplexID idx = errorId;
  if(m_tri_ids.top() + 1 >= m_tri_ids.capacity())
  {
    //resize of the m_tri_node vector..
    const TInt previousCapacityBitVector = m_tri_ids.capacity() + 1;
    m_tri_nodes.resize(previousCapacityBitVector*2);
    m_tri_adj.resize(previousCapacityBitVector*2);
  }

  idx = m_tri_ids.getFreeIndex();
  m_tri_ids.assign(idx);


  /*Building the array of the tri_node before push it back to the m_tri_nodes*/
  TSimplexID* triangle = new TSimplexID[4]{ANode0.getGlobalNode(), ANode1.getGlobalNode(), ANode2.getGlobalNode(), errorId};
  TSimplexID* triangleAdj = new TSimplexID[4]{errorId, errorId, errorId, errorId};

  m_tri_nodes[idx] = triangle;
  m_tri_adj[idx]   = triangleAdj;


  //TODO//
  //m_triangle_variable_manager->addEntry(idx);
  std::vector<TSimplexID>&& ballNode0 = ANode0.ballOf();
  std::vector<TSimplexID>&& ballNode1 = ANode1.ballOf();
  std::vector<TSimplexID>&& ballNode2 = ANode2.ballOf();
  std::vector<TSimplexID> adjSimplex;
  for(auto const & simplexBall0 : ballNode0)
  {
    for(auto const & simplexBall1 : ballNode1)
    {
      for(auto const & simplexBall2 : ballNode2)
      {
        if(simplexBall0 == simplexBall1 && simplexBall0 == simplexBall2)
        {
          adjSimplex.push_back(simplexBall0);
        }
      }
    }
  }

  std::vector<TInt> nodes{ANode0.getGlobalNode(), ANode1.getGlobalNode(), ANode2.getGlobalNode()};
  for(auto const & simplex : adjSimplex)
  {
    std::vector<TInt> otherNode = SimplicesCell(this, simplex).getOtherNodeInSimplex(nodes);
    if(otherNode.size() == 1)
    {
      TInt localNode =  SimplicesCell(this, simplex).getLocalNode(otherNode[0]);

      if(localNode < /*4*/5)
      {
          m_tet_adj[simplex][localNode] = - idx;
      }
      else
      {
        //TODO exception
      }
    }
    else
    {
      //TODO exception
    }
  }

  TSimplexID otherSimplex;
  if(adjSimplex.size() == 0)
  {
    triangle[3]  = errorId;
    triangleAdj[3] = errorId;
  }
  else if(adjSimplex.size() == 1)
  {
    triangle[3] = adjSimplex[0];
    triangleAdj[3] = errorId;
  }
  else if(adjSimplex.size() > 1)
  {
    triangle[3] = adjSimplex[0];
    triangleAdj[3] = adjSimplex[1];
  }

  buildOppFaces( ANode0, ANode1, ANode2);

  return idx;
}


/*---------------------------------------------------------------------------*/
TSimplexID SimplexMesh::addTriangle(const TInt AIndexPoint_0,
                       const TInt AIndexPoint_1,
                       const TInt AIndexPoint_2)
{
  return addTriangle(simplicesNode::SimplicesNode(this, AIndexPoint_0),
                     simplicesNode::SimplicesNode(this, AIndexPoint_1),
                     simplicesNode::SimplicesNode(this, AIndexPoint_2));
}
/*****************************************************************************/
std::vector<TInt> SimplexMesh::deleteTriangle(const TInt ATriIndex)
{
  TSimplexID ATriangleIndex = std::abs(ATriIndex);
  std::vector<TInt> nodesTri;

  if(ATriangleIndex == 0 || ATriangleIndex > m_tri_ids.capacity())
  {
    //TODO
    //Exception
  }
  else
  {
    if(m_tri_ids[ATriangleIndex] != 0)
    {
      int errorId = std::numeric_limits<int>::min();
      nodesTri = std::vector<TSimplexID>{m_tri_nodes[ATriangleIndex][0], m_tri_nodes[ATriangleIndex][1], m_tri_nodes[ATriangleIndex][2]};
      //On recontruit les vector adj des tetra adj a ce triangle.
      TSimplexID tetra0 = m_tri_nodes[ATriangleIndex][3];
      TSimplexID tetra1 = m_tri_adj[ATriangleIndex][3];

      std::vector<TInt> nodesTet0NotInNodesTri;
      std::vector<TInt> nodesTet1NotInNodesTri;

      if(tetra0 != errorId)
      {
        nodesTet0NotInNodesTri = SimplicesCell(this, tetra0).getOtherNodeInSimplex(nodesTri);
      }
      if(tetra1 != errorId)
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
          }
        }
        else
        {
          //TODO exeption
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

      if(m_node_ids[nodesTri[0]] != 0 || m_node_ids[nodesTri[1]] != 0 || m_node_ids[nodesTri[2]] != 0 )
      {
        SimplicesNode Anode0 = SimplicesNode(this, nodesTri[0]);
        SimplicesNode Anode1 = SimplicesNode(this, nodesTri[1]);
        SimplicesNode Anode2 = SimplicesNode(this, nodesTri[2]);
        buildOppFaces(Anode0, Anode1, Anode2);
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
                                    const simplicesNode::SimplicesNode&& ANode3)
{
  if(m_tet_ids.top() + 1 >=  m_tet_ids.capacity())
  {
    //resize of the m_tri_node vector..
    const TInt previousCapacityBitVector = m_tet_ids.capacity() + 1;
    m_tet_nodes.resize(previousCapacityBitVector*2, nullptr);
    m_tet_adj.resize(previousCapacityBitVector*2, nullptr);
  }
  const TInt idx = m_tet_ids.getFreeIndex();

  m_tet_ids.assign(idx);
  /*Building the array of the tri_node before push it back to the m_tri_nodes*/
  TSimplexID* tetra = new TSimplexID[4]{ANode0.getGlobalNode(), ANode1.getGlobalNode(), ANode2.getGlobalNode(), ANode3.getGlobalNode()};
  m_tet_nodes[idx] = tetra;
  //TODO//
  //m_tet_variable_manager->addEntry(idx);

  /*rebuild the base vector after ading aTetra*/
  reorientTetra(idx); /*voir si on peut reorienté que le tetra crée*/
  buildBaseLocal(idx);
  buildBaseAndAdjLocal(idx,
                      ANode0,
                      ANode1,
                      ANode2,
                      ANode3);
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

void SimplexMesh::buildBaseAndAdjLocal(const TSimplexID & tetIndx,
                                      const simplicesNode::SimplicesNode& ANode0,
                                      const simplicesNode::SimplicesNode& ANode1,
                                      const simplicesNode::SimplicesNode& ANode2,
                                      const simplicesNode::SimplicesNode& ANode3)
{
  //std::cout << "adding : " << tetIndx << std::endl;
  /*if(ANode0.getGlobalNode() == 1831 || ANode0.getGlobalNode() == 1826 || ANode0.getGlobalNode() == 1827)
  {
    std::cout << "ANode0 : " << ANode0.getGlobalNode() << std::endl;
    std::cout << "ANode1 : " << ANode1.getGlobalNode() << std::endl;
    std::cout << "ANode2 : " << ANode2.getGlobalNode() << std::endl;
    std::cout << "ANode3 : " << ANode3.getGlobalNode() << std::endl;
    std::cout << "tet index : " << tetIndx << std::endl;
    std::cout << std::endl;

  }*/
  std::vector<TSimplexID> ballNode0;//  = ANode0.ballOf();
  std::vector<TSimplexID> ballNode1;//  = ANode1.ballOf();
  std::vector<TSimplexID> ballNode2;//  = ANode2.ballOf();
  std::vector<TSimplexID> ballNode3;//  = ANode3.ballOf();
  unsigned int tetIter = 0;
  for(tetIter = 0; tetIter < m_tet_ids.capacity(); tetIter++)
  {
    if(m_tet_ids[tetIter] != 0)
    {
      if(tetIter != tetIndx)
      {
        if(SimplicesCell(this, tetIter).containNode(ANode0) )
        {
            ballNode0.push_back(tetIter);
        }
        if(SimplicesCell(this, tetIter).containNode(ANode1))
        {
            ballNode1.push_back(tetIter);
        }
        if(SimplicesCell(this, tetIter).containNode(ANode2))
        {
            ballNode2.push_back(tetIter);
        }
        if(SimplicesCell(this, tetIter).containNode(ANode3))
        {
            ballNode3.push_back(tetIter);
        }
      }
    }
  }

  ballNode0.erase(std::find(ballNode0.begin(), ballNode0.end(), tetIndx), ballNode0.end());
  ballNode1.erase(std::find(ballNode1.begin(), ballNode1.end(), tetIndx), ballNode1.end());
  ballNode2.erase(std::find(ballNode2.begin(), ballNode2.end(), tetIndx), ballNode2.end());
  ballNode3.erase(std::find(ballNode3.begin(), ballNode3.end(), tetIndx), ballNode3.end());



  TInt errorId       = std::numeric_limits<int>::min();
  std::vector<std::vector<TSimplexID>> ballsNode{ballNode0, ballNode1, ballNode2, ballNode3};
  std::vector<SimplicesNode>           Anodes{ANode0, ANode1, ANode2, ANode3};
  TSimplexID* ptrAdj = new TSimplexID[4]{errorId, errorId, errorId, errorId};
  TInt sizeNodeInTet = 4;

  for(unsigned int ballIter = 0; ballIter < sizeNodeInTet; ballIter++)
  {
    std::vector<TSimplexID> ball = ballsNode[ballIter];
    SimplicesNode          Anode = Anodes[ballIter];
    /*if(tetIndx == 485)
    {
      std::cout << "Anode : " << Anode.getGlobalNode() << std::endl;
      std::cout << "ball.size() : " << ball.size() << std::endl;
      for(auto const & tet : ball)
      {
        std::cout << "tet in ball : " << tet << std::endl;
      }
      std::cout << std::endl;
    }*/
    if(ball.size() != 0)
    {
      for(auto const & simplexInBall : ball)
      {
        std::vector<TInt>&&  vecNode = SimplicesCell(this, tetIndx).intersectionNodes(SimplicesCell(this, simplexInBall));
        if(vecNode.size() == 3)
        {
          std::vector<TInt> nodes              = SimplicesCell(this, tetIndx).getOtherNodeInSimplex(vecNode);
          std::vector<TInt> nodesSimplexInBall = SimplicesCell(this, simplexInBall).getOtherNodeInSimplex(vecNode);

          if(nodes.size() == 1 && nodesSimplexInBall.size() == 1)
          {
            TInt localNode              = SimplicesCell(this, tetIndx).getLocalNode(nodes[0]);
            TInt localNodeSimplexInBall = SimplicesCell(this, simplexInBall).getLocalNode(nodesSimplexInBall[0]);

            /*if(tetIndx == 485)
            {
              std::cout << "ball.size() : " << ball.size() << std::endl;
              std::cout << "simplexInBall : " << simplexInBall << std::endl;
              std::cout << "localNodeSimplexInBall : " << localNodeSimplexInBall << std::endl;
            }*/
            /*if(simplexInBall == 35)
            {
              std::cout << "ballIter : " << ballIter <<std::endl;
              for(auto const & node : vecNode)
              {
                std::cout << "node : "<< node << std::endl;
              }
              std::cout << std::endl;
              std::cout << "ball Size " << ball.size() << std::endl;
              std::cout << "simplex 35 adj to : " << tetIndx << " by node : " << localNodeSimplexInBall<< std::endl;
              std::cout << std::endl;
            }*/



            if(m_tet_adj[simplexInBall][localNodeSimplexInBall] < 0)
            {
              if(m_tet_adj[simplexInBall][localNodeSimplexInBall] == errorId)
              {
                ptrAdj[localNode] = simplexInBall;
                m_tet_adj[simplexInBall][localNodeSimplexInBall] = tetIndx;
              }
              else
              {
                TSimplexID triangleID =  - m_tet_adj[simplexInBall][localNodeSimplexInBall];
                ptrAdj[localNode] = m_tet_adj[simplexInBall][localNodeSimplexInBall];
                bool flagAdj = (m_tri_adj[triangleID][3] == errorId) ? true : false;
                (flagAdj) ? m_tri_adj[triangleID][3] = tetIndx : m_tri_nodes[triangleID][3] = tetIndx;
              }
            }
          }
        }
      }
    }
  }
  m_tet_adj[tetIndx] = ptrAdj;


  /*
  TInt errorId       = std::numeric_limits<int>::min();
  std::vector<SimplicesNode>           Anodes{ANode0, ANode1, ANode2, ANode3};
  TSimplexID* ptrAdj = new TSimplexID[4]{errorId, errorId, errorId, errorId};
  TInt sizeNodeInTet = 4;

  for(unsigned int nodeIter = 0; nodeIter < sizeNodeInTet; nodeIter++)
  {
    SimplicesNode          Anode = Anodes[nodeIter];
    if(Anode.ballOf().size() == 0)
    {
      m_base[Anode.getGlobalNode()] = tetIndx;
    }

    std::vector<SimplicesNode> nodesComp {Anodes[(nodeIter + 1) % 4], Anodes[(nodeIter + 2) % 4], Anodes[(nodeIter + 3) % 4]};

    const std::vector<TSimplexID>&& ballNode0  = nodesComp[0].ballOf();
    const std::vector<TSimplexID>&& ballNode1  = nodesComp[1].ballOf();
    const std::vector<TSimplexID>&& ballNode2  = nodesComp[2].ballOf();
    std::cout << "Anode" << std::endl;
    std::cout << Anode << std::endl;
    std::cout << nodesComp[0] << std::endl;
    std::cout << "ballNode0" << std::endl;
    for(auto const & tet : ballNode0)
    {
      std::cout << tet << std::endl;
    }
    std::cout << nodesComp[1] << std::endl;
    std::cout << "ballNode1" << std::endl;
    for(auto const & tet : ballNode1)
    {
      std::cout << tet << std::endl;
    }
    std::cout << nodesComp[2] << std::endl;
    std::cout << "ballNode2" << std::endl;
    for(auto const & tet : ballNode2)
    {
      std::cout << tet << std::endl;
    }
    std::vector<TSimplexID> ballConcat;
    for(auto const & ball0 : ballNode0)
    {
      if(std::find(ballNode1.begin(), ballNode1.end(), ball0) != ballNode1.end())
      {
        if(std::find(ballNode2.begin(), ballNode2.end(), ball0) != ballNode2.end())
        {
            ballConcat.push_back(ball0);
        }
      }
    }
    std::cout << "ballConcat" << std::endl;
    for(auto const & tet : ballConcat)
    {
      std::cout << tet << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    ballConcat.erase(std::remove(ballConcat.begin(), ballConcat.end(), tetIndx), ballConcat.end());


    if(Anode.ballOf().size() == 0)
    {
      TInt localNode = SimplicesCell(this, tetIndx).getLocalNode(Anode.getGlobalNode());
      ptrAdj[localNode] = errorId;
    }
    if(ballConcat.size() == 1)
    {
      if(ballConcat[0] >= 0)
      {
        TSimplexID simplex = ballConcat[0];
        if(simplex >= 0)
        {
          std::vector<TInt> vecNode{nodesComp[0].getGlobalNode(), nodesComp[1].getGlobalNode(), nodesComp[2].getGlobalNode()};
          std::vector<TInt> nodesSimplex = SimplicesCell(this, simplex).getOtherNodeInSimplex(vecNode);

          TInt localNode              = SimplicesCell(this, tetIndx).getLocalNode(Anode.getGlobalNode());
          TInt localNodeSimplex       = SimplicesCell(this, simplex).getLocalNode(nodesSimplex[0]);

          ptrAdj[localNode] = simplex;
          m_tet_adj[simplex][localNodeSimplex] = tetIndx;
        }
        else
        {*/
          // triangle
          /*std::cout << "cogno " << std::endl;
          TInt localNode              = SimplicesCell(this, tetIndx).getLocalNode(Anode.getGlobalNode());
          TSimplexID triangleID =  ballConcat[0];
          ptrAdj[localNode] = triangleID;
          bool flagAdj = (m_tri_adj[-triangleID][3] == errorId) ? true : false;
          (flagAdj == true) ? m_tri_adj[-triangleID][3] = tetIndx : m_tri_nodes[-triangleID][3] = tetIndx;*/
        /*}
      }
    }
    else //size == 2
    {*/
      /*std::cout << "tetIdx : "<< tetIndx << std::endl;
      std::cout << "ballConcat.size() : "<< ballConcat.size() << std::endl;
      std::cout << "triangleID : " <<  ballConcat[1] << " " <<  ballConcat[0] << std::endl;

      TSimplexID triangleID = (ballConcat[0] < 0)? ballConcat[0] : ballConcat[1];
      TInt localNode        = SimplicesCell(this, tetIndx).getLocalNode(Anode.getGlobalNode());
      ptrAdj[localNode]     = triangleID;
      bool flagAdj          = (m_tri_adj[-triangleID][3] == errorId) ? true : false;
      std::cout << "OKOK" << std::endl;
      (flagAdj ) ? m_tri_adj[-triangleID][3] = tetIndx : m_tri_nodes[-triangleID][3] = tetIndx;*/
    /*}
  }*/
  //m_tet_adj[tetIndx] = ptrAdj;
}
/******************************************************************************/
void SimplexMesh::buildAdjTet(const TSimplexID & currentTet , const simplicesNode::SimplicesNode& ANode0,
                                                        const simplicesNode::SimplicesNode& ANode1,
                                                        const simplicesNode::SimplicesNode& ANode2,
                                                        const simplicesNode::SimplicesNode& ANode3)
{
  int errorId = std::numeric_limits<int>::min();

  const std::vector<TSimplexID>&& ballNode0  = ANode0.ballOf();

  const std::vector<TSimplexID>&& ballNode1  = ANode1.ballOf();

  const std::vector<TSimplexID>&& ballNode2  = ANode2.ballOf();

  const std::vector<TSimplexID>&& ballNode3  = ANode3.ballOf();

  std::vector<TSimplexID> concatenationBall;


  concatenationBall.insert(concatenationBall.end(), ballNode0.begin(), ballNode0.end());
  concatenationBall.insert(concatenationBall.end(), ballNode1.begin(), ballNode1.end());
  concatenationBall.insert(concatenationBall.end(), ballNode2.begin(), ballNode2.end());
  concatenationBall.insert(concatenationBall.end(), ballNode3.begin(), ballNode3.end());

  TSimplexID* nodePtr    = m_tet_nodes[currentTet];
  TSimplexID* faceTet0   =  nodePtr;
  TSimplexID* faceTet1   = (nodePtr + 1);
  TSimplexID* faceTet2   = (nodePtr + 2);
  TSimplexID* faceTet3   = (nodePtr + 3);
  const std::vector<TSimplexID*> faces{faceTet0, faceTet1, faceTet2, faceTet3};
  TSimplexID* adjData = (TSimplexID*)malloc(4*sizeof(TSimplexID));
  if(ballNode0.size() != 0)
  {
    for(int faceIter = 0; faceIter < faces.size(); faceIter++)
    {
      int nodelocalAdjToface = (faceIter + 3 ) % 4;

      TSimplexID* face = faces[faceIter];

      for(auto const & simplex : concatenationBall)
      {
        if(containNodes<3,4>(face, m_tet_nodes[simplex]))
        {
          adjData[nodelocalAdjToface] = simplex;
        }
        else
        {
          adjData[nodelocalAdjToface] = errorId;
        }
      }
    }
  }
  m_tet_adj[currentTet] = adjData;
}
/******************************************************************************/
void SimplexMesh::buildOppFaces(const simplicesNode::SimplicesNode& ANode0,
                                const simplicesNode::SimplicesNode& ANode1,
                                const simplicesNode::SimplicesNode& ANode2)
{
  std::vector<TSimplexID> shellTriangleOrdererdWithHole01 = ANode0.shellTriangleOrdererdWithHole(ANode1);
  std::vector<TSimplexID> shellTriangleOrdererdWithHole12 = ANode1.shellTriangleOrdererdWithHole(ANode2);
  std::vector<TSimplexID> shellTriangleOrdererdWithHole20 = ANode2.shellTriangleOrdererdWithHole(ANode0);

  std::vector<TInt> nodes01{ANode0.getGlobalNode(), ANode1.getGlobalNode()};
  std::vector<TInt> nodes12{ANode1.getGlobalNode(), ANode2.getGlobalNode()};
  std::vector<TInt> nodes20{ANode2.getGlobalNode(), ANode0.getGlobalNode()};
  for(unsigned int triangleIter = 0; triangleIter < shellTriangleOrdererdWithHole01.size(); triangleIter++)
  {
    unsigned int triangle = - shellTriangleOrdererdWithHole01[triangleIter];
    if(m_tri_ids[triangle] != 0)
    {
      std::vector<TInt> nodeExt = SimplicesTriangle(this, triangle).getOtherNodeInSimplex(nodes01);
      if(nodeExt.size() == 1)
      {
        TInt localNode = SimplicesTriangle(this,triangle).getLocalNode(nodeExt[0]);
        unsigned int triangleNext = - shellTriangleOrdererdWithHole01[(triangleIter + 1) % shellTriangleOrdererdWithHole01.size()];
        if(triangle != triangleNext)
        {
            m_tri_adj[triangle][localNode] = triangleNext;
        }
      }
      else
      {
        //TODO exception
      }
    }
    else
    {
      //TODO exception
    }

  }
  for(unsigned int triangleIter = 0; triangleIter < shellTriangleOrdererdWithHole12.size(); triangleIter++)
  {
    unsigned int triangle = - shellTriangleOrdererdWithHole12[triangleIter];
    if(m_tri_ids[triangle] != 0)
    {
      std::vector<TInt> nodeExt = SimplicesTriangle(this, triangle).getOtherNodeInSimplex(nodes12);
      if(nodeExt.size() == 1)
      {
        TInt localNode = SimplicesTriangle(this,triangle).getLocalNode(nodeExt[0]);
        unsigned int triangleNext = - shellTriangleOrdererdWithHole12[(triangleIter + 1) % shellTriangleOrdererdWithHole12.size()];
        if(triangle != triangleNext)
        {
            m_tri_adj[triangle][localNode] = triangleNext;
        }
      }
      else
      {
        //TODO exception
      }
    }
    else
    {
      //TODO exception
    }
  }

  for(unsigned int triangleIter = 0; triangleIter < shellTriangleOrdererdWithHole20.size(); triangleIter++)
  {
    unsigned int triangle = - shellTriangleOrdererdWithHole20[triangleIter];
    if(m_tri_ids[triangle] != 0)
    {
      std::vector<TInt> nodeExt = SimplicesTriangle(this, triangle).getOtherNodeInSimplex(nodes20);
      if(nodeExt.size() == 1)
      {
        TInt localNode = SimplicesTriangle(this,triangle).getLocalNode(nodeExt[0]);
        unsigned int triangleNext = - shellTriangleOrdererdWithHole20[(triangleIter + 1) % shellTriangleOrdererdWithHole20.size()];
        if(triangle != triangleNext)
        {
            m_tri_adj[triangle][localNode] = triangleNext;
        }
      }
      else
      {
        //TODO exception
      }
    }
    else
    {
      //TODO exception
    }
  }
}
/*---------------------------------------------------------------------------*/
TSimplexID SimplexMesh::addTetraedre(const TInt AIndexPoint0,
                                    const TInt AIndexPoint1,
                                    const TInt AIndexPoint2,
                                    const TInt AIndexPoint3)
{
  return addTetraedre(simplicesNode::SimplicesNode(this, AIndexPoint0),
                      simplicesNode::SimplicesNode(this, AIndexPoint1),
                      simplicesNode::SimplicesNode(this, AIndexPoint2),
                      simplicesNode::SimplicesNode(this, AIndexPoint3));
}
/******************************************************************************/
std::vector<TInt> SimplexMesh::deleteTetra(const TInt ATetraIndex)
{
  int errorId = std::numeric_limits<int>::min();
  std::vector<TInt> nodes;

  if(ATetraIndex < 0 || ATetraIndex > m_tet_ids.capacity())
  {
    /*TODO execption*/
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
            std::vector<TSimplexID> && ball = SimplicesNode(this, node).ballOf();
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
          if(simplexAdj != errorId)
          {
            if(simplexAdj >= 0)
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

              //on supprime le triangle si il est adjacent a rien...à discuter avec Franck
              if(m_tri_nodes[-simplexAdj][3] == errorId && m_tri_adj[-simplexAdj][3] == errorId)
              {
                deleteTriangle(simplexAdj);
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
  float sizeLayer = 0.2f;
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
        PointInsertion pi(this, SimplicesNode(this, nodeToInsert), criterionRAIS, initialCavity, markedSimplex);
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
TSimplexID SimplexMesh::getOppositeFace(const TInt ANodeID, const TSimplexID ATriID)
{
  size_t localIndex             = -1;
  const size_t indexBufferSize  =  3;
  const TSimplexID* indexBuffer = m_tri_nodes[ATriID];
  if(ATriID >= 0)
  {
    for(size_t i = 0; i< indexBufferSize ; i++)
    {
      if(*(indexBuffer + i) == ANodeID)
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
void SimplexMesh::buildOppFacesVector()
{
  const size_t TriangleNodeSize = 3;
  const int errorId             = std::numeric_limits<int>::min();

  m_tri_adj.clear();
  m_tri_adj.resize(m_tri_nodes.size());


  for(TSimplexID currentTri = 1; currentTri < m_tri_nodes.size(); currentTri++)
  {
    if(m_tri_ids[currentTri] != 0)
    {
      TSimplexID* trianglesOpp      = new TSimplexID[3]{errorId, errorId, errorId};;
      for(TInt nodeTri = 0; nodeTri < TriangleNodeSize; nodeTri++)
      {
          trianglesOpp[nodeTri] = buildOppositeFacesVector(currentTri, m_tri_nodes[currentTri][nodeTri]);
      }
      m_tri_adj[currentTri]    = trianglesOpp;
    }
  }

  /*Building of the adjacente simplex of tirangles (m_tri_nodes[3] & m_tri_adj[3]...)*/
  std::vector<TSimplexID> v;
  for(TSimplexID currentTri = 1 /*start at 1*/; currentTri < m_tri_nodes.size(); currentTri++)
  {
    if(m_tri_ids[currentTri] != 0)
    {
      v = buildAdjTriVector(currentTri);
      /*TODO rajouter une relation d'ordre grâce aux normales..*/
      m_tri_nodes[currentTri][3] = v[0];
      m_tri_adj  [currentTri][3] = v[1];
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

      const TSimplexID* arrayComp = m_tri_nodes[compTri];
      flag = containNodes<2,3>(nodesTri, arrayComp);
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
std::vector<TSimplexID> SimplexMesh::buildAdjTriVector(const TSimplexID currentTriangle)
{
  const int errorId = std::numeric_limits<int>::min();
  std::vector<TSimplexID> v{errorId, errorId};
  int iter = 0;

  TSimplexID* nodeTri = m_tri_nodes[currentTriangle];

  for(TSimplexID Tet = 0; Tet < m_tet_nodes.size(); Tet++)
  {
    if(m_tet_ids[Tet] != 0)
    {
      TSimplexID* arrayComp = m_tet_nodes[Tet];
      bool flag                  = containNodes<3,4>(nodeTri, arrayComp);
      if(flag && iter < 2)
      {
        v[iter]             = Tet;
        iter++;
      }
    }
  }

  return v;
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

      /*std::cout << "u : " << u << std::endl;
      std::cout << "v : " << v << std::endl;
      std::cout << "w : " << w << std::endl;
      std::cout << "t : " << t << std::endl;
      std::cout << "tetra : " << tetra << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;*/

  }


  tetra = (tetra >= 0 || tetra == errorId)? tetra :
          (m_tri_nodes[-tetra][3] != previousTet || m_tri_nodes[-tetra][3] == errorId)? m_tri_nodes[-tetra][3] :
           m_tri_adj[-tetra][3];

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
template<typename T>
void SimplexMesh::intersection(std::vector<std::vector<T>>& buffer0, const std::vector<T>& buffer1)
{
  for(auto & simplices : buffer0)
  {
      simplices.erase(std::remove_if(simplices.begin(), simplices.end(), [&](T node)
      {
        bool flag = false;
        if(std::find(buffer1.begin(), buffer1.end(), node) == buffer1.end())
        {
          flag = true;
        };
        return flag;
      }), simplices.end());
  }
}
/******************************************************************************/
template<typename T>
void SimplexMesh::intersection(std::vector<T>& buffer0, const std::vector<T>& buffer1)
{
  buffer0.erase(std::remove_if(buffer0.begin(), buffer0.end(), [&](T node)
                {
                  bool flag = false;
                  if(std::find(buffer1.begin(), buffer1.end(), node) == buffer1.end())
                  {
                    flag = true;
                  }
                  return flag;
                }
              ), buffer0.end());
}
/******************************************************************************/
std::vector<TSimplexID> SimplexMesh::intersectionSimplex(std::vector<std::vector<TSimplexID>> & vectorOfBalls)
{
  std::vector<TSimplexID> intersectionVec;
/*
  for(auto & vector : vectorOfBalls)
  {
    std::sort(vector.begin(), vector.end());
  }

  for(int idxVector = 0 ; idxVector < vectorOfBalls.size() -1 ; idxVector++)
  {
    std::set_intersection(vectorOfBalls[idxVector].begin(), vectorOfBalls[idxVector].end(),
                          vectorOfBalls[idxVector + 1].begin(), vectorOfBalls[idxVector + 1].end(),
                          std::back_inserter(intersectionSec));
  }

  for(const auto & simplex : intersectionSec)
  {
    intersectionVec.push_back(simplex);
  }
*/
  return std::move(intersectionVec);
}
/******************************************************************************/
void SimplexMesh::saveAtSimplexMeshFormat(const std::string && destination, const std::string&& name)
{
  int errorId = std::numeric_limits<int>::min();
  std::ofstream myfile;
  const std::string SimplexMeshFormat = destination + name + ".simplexMesh";
  myfile.open (SimplexMeshFormat);
  if(!myfile.is_open())
  {
    /*TODO exception*/
    std::cout << "impossible d'ouvrir le fichier : " << SimplexMeshFormat << std::endl;
  }

  //m_node_ids
  myfile << "m_node_ids" << "\n";
  myfile << m_node_ids.size() << "\n";
  for(int iter = 0; iter < m_node_ids.size() ; iter ++)
  {
    myfile << m_node_ids[iter] << "\n";
  }

  //m_coords
  myfile << "m_coords" << "\n";
  myfile << m_coords.size() << "\n";
  for(auto const & coord : m_coords)
  {
    myfile << coord[0] << " " << coord[1] << " " << coord[2] << "\n";
  }

  //m_base
  myfile << "m_base" << "\n";
  myfile << m_base.size() << "\n";
  for(auto const & base : m_base)
  {
    myfile << base << "\n";
  }
  /****************************************************************************/
  //m_node_ids
  myfile << "m_tet_ids" << "\n";
  myfile << m_tet_ids.size() << "\n";
  for(int iter = 0; iter < m_tet_ids.size() ; iter ++)
  {
    myfile << m_tet_ids[iter] << "\n";
  }

  //m_tet_nodes
  myfile << "m_tet_nodes" << "\n";
  myfile << m_tet_nodes.size() << "\n";
  for(int iter = 0; iter < m_tet_nodes.size() ; iter ++)
  {
    if(m_tet_nodes[iter] != nullptr)
    {
      myfile << *m_tet_nodes[iter] << " " << *(m_tet_nodes[iter] + 1) << " " << *(m_tet_nodes[iter] + 2) << " "<< *(m_tet_nodes[iter] + 3) << "\n";
    }
    else
    {
      myfile << errorId << " " << errorId << " " << errorId << " "<< errorId << "\n";
    }
  }

  //m_tet_adj
  myfile << "m_tet_adj" << "\n";
  myfile << m_tet_adj.size() << "\n";
  for(int iter = 0; iter < m_tet_adj.size() ; iter ++)
  {
    if(m_tet_adj[iter] != nullptr)
    {
        myfile << *m_tet_adj[iter] << " " << *(m_tet_adj[iter] + 1) << " " << *(m_tet_adj[iter] + 2) << " "<< *(m_tet_adj[iter] + 3) << "\n";
    }
    else
    {
      myfile << errorId << " " << errorId << " " << errorId << " "<< errorId << "\n";
    }
  }
  /****************************************************************************/

  //m_node_ids
  myfile << "m_tri_ids" << "\n";
  myfile << m_tri_ids.size() << "\n";
  for(int iter = 0; iter < m_tri_ids.size() ; iter ++)
  {
    myfile << m_tri_ids[iter] << "\n";
  }

  //m_tri_nodes
  myfile << "m_tri_nodes" << "\n";
  myfile << m_tri_nodes.size() << "\n";
  for(int iter = 0; iter < m_tri_nodes.size() ; iter ++)
  {
    if(m_tri_nodes[iter] != nullptr)
    {
      myfile << *m_tri_nodes[iter] << " " << *(m_tri_nodes[iter] + 1) << " " << *(m_tri_nodes[iter] + 2) << *(m_tri_nodes[iter] + 3) << "\n";
    }
    else
    {
      myfile << errorId << " " << errorId << " " << errorId << " "<< errorId << "\n";
    }
  }

  //m_tri_adj
  myfile << "m_tri_adj" << "\n";
  myfile << m_tri_adj.size() << "\n";
  for(int iter = 0; iter < m_tri_adj.size() ; iter ++)
  {
    if(m_tri_adj[iter] != nullptr)
    {
      myfile << *m_tri_adj[iter] << " " << *(m_tri_adj[iter] + 1) << " " << *(m_tri_adj[iter] + 2) << *(m_tri_adj[iter] + 3) << "\n";
    }
    else
    {
      myfile << errorId << " " << errorId << " " << errorId << " "<< errorId << "\n";
    }
  }
  /****************************************************************************/
  myfile.close();
}
/****************************************************************************/
std::vector<VariableItf*> SimplexMesh::getAllVariables(ECellType AType) const
{
	switch (AType){
	case GMDS_NODE:{
					   return m_node_variable_manager->getAllVariables();
	}
		break;
	default:
		throw GMDSException("Unmanaged type of value -> impossible to access to a variable");
	}
}
/******************************************************************************/
void SimplexMesh::loadSimplexMeshFormat(const std::string && file)
{
  //TODO
  //mettre une condition sur l'extension verifier que c'est bien .simplexMesh
  if(1)
  {
    std::ifstream myfile;
    myfile.open(file);
    if(!myfile.is_open())
    {
      /*TODO exceptio*/
      std::cout << "impossible d'ouvrir le fichier : " << file << std::endl;
    }
    else
    {
      std::string line;
      getline (myfile,line);
      {
          if(line == "m_node_ids")
          {
            getline (myfile,line);
            int nbr_line = std::stoi(line);
            m_node_ids.clear();
            m_node_ids.resize(nbr_line);
            for(int iter = 0; iter < nbr_line ; iter++)
            {
              getline (myfile,line);
              int flag = std::stoi(line);
              if(flag == 1)
              {
                m_node_ids.assign(iter);
              }
            }
          }
          std::cout << "M_NODE_IDS PASSED" << std::endl;
          getline (myfile,line);
          if(line == "m_coords")
          {
            getline (myfile,line);
            int nbr_line = std::stoi(line);
            m_coords.clear();
            m_coords.resize(nbr_line);
            for(int iter = 0; iter < nbr_line ; iter++)
            {
              gmds::math::Point pt;
              getline (myfile,line,' ');
              pt[0] = std::stoi(line);
              getline (myfile,line,' ');
              pt[1] = std::stoi(line);
              getline (myfile,line);
              pt[2] = std::stoi(line);
              m_coords[iter] = pt;
            }
          }
          std::cout << "M_COORDS PASSED" << std::endl;

          getline (myfile,line);
          if(line == "m_base")
          {
            getline (myfile,line);
            int nbr_line = std::stoi(line);
            m_base.clear();
            m_base.resize(nbr_line);
            for(int iter = 0; iter < nbr_line ; iter++)
            {
              getline (myfile,line);
              int base = std::stoi(line);
              m_base[iter] = base;
            }
          }
          std::cout << "M_BASE PASSED" << std::endl;
          getline (myfile,line);

          if(line == "m_tet_ids")
          {
            getline (myfile,line);
            int nbr_line = std::stoi(line);
            m_tet_ids.clear();
            m_tet_ids.resize(nbr_line);
            for(int iter = 0; iter < nbr_line ; iter++)
            {
              getline (myfile,line);
              int flag = std::stoi(line);
              if(flag == 1)
              {
                m_tet_ids.assign(iter);
              }
            }
          }
          std::cout << "M_TET_IDS PASSED" << std::endl;
          getline (myfile,line);

          if(line == "m_tet_nodes")
          {
            getline (myfile,line);
            int nbr_line = std::stoi(line);

            m_tet_nodes.clear();
            m_tet_nodes.resize(nbr_line);
            for(int iter = 0; iter < nbr_line ; iter++)
            {
              getline (myfile,line,' ');
              TSimplexID nodeA = std::stoi(line);
              getline (myfile,line,' ');
              TSimplexID nodeB = std::stoi(line);
              getline (myfile,line,' ');
              TSimplexID nodeC = std::stoi(line);
              getline (myfile,line);
              TSimplexID nodeD = std::stoi(line);

              TSimplexID* bufferNode = new TSimplexID[4]{nodeA, nodeB, nodeC, nodeD};
              m_tet_nodes[iter] = bufferNode;
            }

          }
          std::cout << "M_TET_NODE PASSED" << std::endl;
          getline (myfile,line);
          if(line == "m_tet_adj")
          {
            getline (myfile,line);
            int nbr_line = std::stoi(line);

            m_tet_adj.clear();
            m_tet_adj.resize(nbr_line);
            for(int iter = 0; iter < nbr_line ; iter++)
            {
              getline (myfile,line,' ');
              TSimplexID nodeA = std::stoi(line);
              getline (myfile,line,' ');
              TSimplexID nodeB = std::stoi(line);
              getline (myfile,line,' ');
              TSimplexID nodeC = std::stoi(line);
              getline (myfile,line);
              TSimplexID nodeD = std::stoi(line);

              TSimplexID* bufferNode = new TSimplexID[4]{nodeA, nodeB, nodeC, nodeD};
              m_tet_adj[iter] = bufferNode;
            }
          }
          std::cout << "M_TET_ADJ PASSED" << std::endl;
          /*getline (myfile,line);
          if(line == "m_tri_ids")
          {
            getline (myfile,line);
            int nbr_line = std::stoi(line);
            m_tri_ids.clear();
            m_tri_ids.resize(nbr_line);
            for(int iter = 0; iter < nbr_line ; iter++)
            {
              getline (myfile,line);
              int flag = std::stoi(line);
              if(flag == 1)
              {
                m_tri_ids.assign(iter);
              }
            }
          }
          std::cout << "M_TRi_IDS PASSED" << std::endl;
          getline (myfile,line);

          if(line == "m_tri_nodes")
          {
            getline (myfile,line);
            int nbr_line = std::stoi(line);

            m_tri_nodes.clear();
            m_tri_nodes.resize(nbr_line);
            for(int iter = 0; iter < nbr_line ; iter++)
            {
              getline (myfile,line,' ');
              TSimplexID nodeA = std::stoi(line);
              getline (myfile,line,' ');
              TSimplexID nodeB = std::stoi(line);
              getline (myfile,line);
              TSimplexID nodeC = std::stoi(line);
              std::cout << nodeA << " " << nodeB << " " << nodeC <<std::endl;
              TSimplexID* bufferNode = new TSimplexID[4]{nodeA, nodeB, nodeC};
              m_tri_nodes[iter] = bufferNode;
            }
          }
          else if(line == "m_tri_adj")
          {
            //std::cout << line << std::endl;
          }*/
          std::cout << "ALL PASSED" << std::endl;
      }
    }
  }
}
