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

  myfile.open("data_AdaptationMesh.txt");
  if(myfile.is_open())
  {
    myfile << "DATA ADAPTAION FOR VISUALIZATION" << std::endl;
  }
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
  if(simplices.size() != 0)
  {
    for(unsigned int tet = 0 ; tet < m_tet_ids.capacity(); tet++)
    {
      if(std::find(simplices.begin(), simplices.end(), tet) == simplices.end())
      {
        m_tet_ids.unselect(tet);
      }
    }

    for(unsigned int tri = 1 ; tri < m_tri_ids.capacity(); tri++)
    {
      if(std::find(simplices.begin(), simplices.end(), -tri) == simplices.end())
      {
        m_tri_ids.unselect(tri);
      }
    }
  }
  else
  {
    for(unsigned int tet = 0 ; tet < m_tet_ids.capacity(); tet++)
    {
      m_tet_ids.unselect(tet);
    }
  }
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
  double epsilon = 10E-4;

  //check if the point pt is on any other node already in the mesh;
  if(tetraContainingPt.size() == 0)
  {
    flag = checkSimplicesContenaing(pt, tetraContainingPt, simplexToCheckFirst);
  }


  if(flag)
  {
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

  if(m_tet_nodes.size() == 0)
  {
    for(unsigned int nodeId = 0 ; nodeId < m_node_ids.capacity() ; nodeId++)
    {
      if(m_node_ids[nodeId] != 0)
      {
        math::Vector3d VecBetweenPt = (SimplicesNode(this, nodeId).getCoords() - pt);
        double lenght = VecBetweenPt.norm();
        if(lenght < epsilon)//pt is on nodeOfTet
        {
          alreadyAdd = true;
          return nodeId;
        }
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
  m_node_variable_manager->addEntry(idx);
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
void SimplexMesh::moveNodeCoord(const TInt node, const math::Point& newCoord)
{
  if(m_node_ids[node] != 0)
  {
    m_coords[node] = newCoord;
  }
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
  unsigned int sizeEdge = 3;
  unsigned int sizeFace = 4;

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
      SimplicesTriangle currentTriangle = SimplicesTriangle(this, ATriangleIndex);
      std::vector<TSimplexID> neighborTriangles = currentTriangle.adjacentTriangle();
      for(auto const neighborTriangle : neighborTriangles)
      {
        //if neighborTriangle == border, neighborTriangle has already been destroy it's not a bug
        if(neighborTriangle != border)
        {
          bool flag = false;
          for(unsigned int edge = 0 ; edge < sizeEdge ; edge++)
          {
            if(m_tri_adj[neighborTriangle][edge] == ATriangleIndex)
            {
              flag = true;
              m_tri_adj[neighborTriangle][edge] = border;
              break;
            }
          }
          if(!flag)
          {
            /*todo Exception*/
            /*gmds::ISimplexMeshIOService ioServiceMesh(this);
            std::vector<TSimplexID> v{ATetraIndex, simplexAdj};
            deleteAllSimplicesBut(v);
            gmds::VTKWriter vtkWriterMeshER(&ioServiceMesh);
            vtkWriterMeshER.setCellOptions(gmds::N|gmds::R|gmds::F);
            vtkWriterMeshER.setDataOptions(gmds::N|gmds::R|gmds::F);
            vtkWriterMeshER.write("ADJACENT_PROBLEM.vtk");*/
            std::cout << "Triangle "<< ATriangleIndex << " is not adjacent to triangle " << neighborTriangle << std::endl;
            throw gmds::GMDSException("Triangle adjacent relation problem !");
          }
        }
      }
      //On recontruit les relation du tetra voisin du triangle supprimé
      TSimplexID neigborTetra = m_tri_nodes[ATriangleIndex][3];
      if(neigborTetra != border)
      {
        bool flag = false;
        for(unsigned int face = 0 ; face < sizeFace ; face++)
        {
          if(m_tet_adj[neigborTetra][face] == -ATriangleIndex)
          {
            flag = true;
            m_tet_adj[neigborTetra][face] = border;
            break;
          }
        }
        if(!flag)
        {
          /*todo Exception*/
          /*gmds::ISimplexMeshIOService ioServiceMesh(this);
          std::vector<TSimplexID> v{ATetraIndex, simplexAdj};
          deleteAllSimplicesBut(v);
          gmds::VTKWriter vtkWriterMeshER(&ioServiceMesh);
          vtkWriterMeshER.setCellOptions(gmds::N|gmds::R|gmds::F);
          vtkWriterMeshER.setDataOptions(gmds::N|gmds::R|gmds::F);
          vtkWriterMeshER.write("ADJACENT_PROBLEM.vtk");*/
          std::cout << "Triangle "<< ATriangleIndex << " is not adjacent to tetra " << neigborTetra << std::endl;
          throw gmds::GMDSException("Triangle adjacent relation problem !");
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
      if(m_tet_ids[tetIter] != 0)
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
void SimplexMesh::initializeEdgeStructure()
{
  Variable<TInt>* BND_CURVE_COLOR = nullptr;
  Variable<TInt>* BND_VERTEX_COLOR = nullptr;
  try {
    BND_CURVE_COLOR = getVariable<TInt,SimplicesNode>("BND_CURVE_COLOR");
    BND_VERTEX_COLOR = getVariable<TInt,SimplicesNode>("BND_VERTEX_COLOR");
  }catch (GMDSException e) {
    throw gmds::GMDSException(e);
  }

  std::multimap<TInt, std::vector<TInt>> labels2Nodes{};
  std::multimap<TInt, std::vector<TInt>> corner2label{};

  for(TInt nodeIdx = 0 ; nodeIdx < m_node_ids.capacity() ; nodeIdx++)
  {
    if(m_node_ids[nodeIdx] != 0)
    {
      TInt labelCurve = (*BND_CURVE_COLOR)[nodeIdx];
      TInt labelCorner = (*BND_VERTEX_COLOR)[nodeIdx];
      if(labelCurve != 0)
      {
        auto iter = labels2Nodes.find(labelCurve);
        if(iter != labels2Nodes.end())
        {
          iter->second.push_back(nodeIdx);
        }
        else
        {
          std::vector<TInt> v{nodeIdx};
          labels2Nodes.insert(std::make_pair(labelCurve,v));
        }
      }
      else if(labelCorner != 0)
      {
        const std::vector<TInt> nodes = SimplicesNode(this, nodeIdx).neighborNodes();
        auto iter = corner2label.find(nodeIdx);
        std::vector<TInt> v{};
        if(iter == corner2label.end())
        {
          for(auto const node : nodes)
          {
            if(m_node_ids[node] != 0)
            {
              if((*BND_CURVE_COLOR)[node] != 0)
              {
                v.push_back(node);
              }
            }
            else
            {
              throw gmds::GMDSException("m_node_ids[node] == 0");
            }
          }
          corner2label.insert(std::make_pair(nodeIdx, v));
        }
      }
    }
  }

  gmds::BitVector nodeCycling(m_node_ids.capacity());
  for(auto const data : labels2Nodes)
  {
    for(auto const node0 : data.second)
    {
        const SimplicesNode sNode0 = SimplicesNode(this, node0);
        for(auto const node1 : data.second)
        {
          if(nodeCycling[node1] == 0)
          {
            if(node0 != node1)
            {
              std::pair<TInt, TInt> p{std::min(node0, node1), std::max(node0, node1)};
              const SimplicesNode sNode1 = SimplicesNode(this, node1);
              if(sNode0.shell(sNode1).size() != 0)
              {
                m_edgesStructure.insert(std::make_pair(data.first, p));
              }
            }
          }
        }
        nodeCycling.assign(node0);
    }
  }

  //Curve nodes have been connected together
  //Connection to corner node with curve node
  for(auto const c2n : corner2label)
  {
    TInt corner = c2n.first;
    for(auto const curveNode : c2n.second)
    {
      TInt labelCurve = (*BND_CURVE_COLOR)[curveNode];
      std::pair<TInt, TInt> p{std::min(corner, curveNode), std::max(corner, curveNode)};
      m_edgesStructure.insert(std::make_pair(labelCurve, p));
    }
  }


  /*for(auto const data : m_edgesStructure)
  {
    std::cout << "data --> " << data.first << " | [" << data.second.first << " : " << data.second.second << "]" << std::endl;
  }
  std::cout << "INITIALISATION EDGE STRUCTURE END" << std::endl;*/

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
void SimplexMesh::fillHexahedron(const TInt ANode0, const TInt ANode1, const TInt ANode2, const TInt ANode3,
                    const TInt ANode4, const TInt ANode5, const TInt ANode6, const TInt ANode7)
{
  std::vector<TInt> nodes{ANode0, ANode1, ANode2, ANode3, ANode4, ANode5, ANode6, ANode7};
  m_hexahedronData.push_back(nodes);
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
void SimplexMesh::checkMesh()
{
  TInt border = std::numeric_limits<int>::min();
  unsigned int tetraNbr = getNbTetra();
  unsigned int cpt = 0;
  unsigned int sizeFace = 3;

  for(unsigned int cell0 = 0; cell0 < m_tet_ids.capacity(); cell0++)
  {
    if(m_tet_ids[cell0] == 1)
    {
      const SimplicesCell c0 = SimplicesCell(this, cell0);
      for(unsigned int face = 0; face < sizeFace; face++)
      {
        TSimplexID oppositeCell = m_tet_adj[cell0][face];
        if(oppositeCell == border)
        {
          std::cout << "c0 -> " << c0 << std::endl;
          throw gmds::GMDSException("oppositeCell == border");
        }

        if(oppositeCell >= 0)
        {
          const SimplicesCell c1 = SimplicesCell(this, oppositeCell);
          std::vector<TInt> nodesFace =  c0.getOrderedFace(face);
          if(nodesFace.size() != 3)
          {
            throw gmds::GMDSException("nodeFace.size() != 3");
          }

          std::vector<TInt> otherNode = c1.getOtherNodeInSimplex(nodesFace);
          if(otherNode.size() != 1)
          {
            throw gmds::GMDSException("otherNode.size() != 1");
          }

          if(m_tet_adj[oppositeCell][c1.getLocalNode(otherNode.front())] != cell0)
          {
            std::cout << "cell0 -> " << c0 << std::endl;
            std::cout << "opposite Cell -> " << c1 << std::endl;
            throw gmds::GMDSException("m_tet_adj[oppositeCell][otherNode.size()] != cell0");
          }
        }
        else
        {

        }
      }
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
                    TSimplexID adjSimplex = m_tet_adj[cell0][localNode0];
                    if(adjSimplex >= 0)
                    {
                      const SimplicesCell adjCell = SimplicesCell(this, adjSimplex);
                      //2 plane place tetrahedron can exist so we have to ensure that one of then is realy adjacent to cell
                      if(adjCell.intersectionNodes(cell).size() != 3)
                      {
                        std::cout << "Simplices [" << cell0 << ";" << cell1 << "] are adjacent but adjcent vector is wrong" << std::endl;
                        std::cout << "cell0 -> " << cell << std::endl;
                        std::cout << "cell1 -> " << cellToCompare << std::endl;
                        std::vector<TSimplexID> v{static_cast<int>(cell0), static_cast<int>(cell1)};
                        std::vector<TSimplexID> adjSimplices = cell.adjacentTetra();
                        for(auto const simplex : adjSimplices)
                        {
                          if(simplex >= 0)
                          {
                            v.push_back(simplex);
                          }
                        }
                        deleteAllSimplicesBut(v);
                        gmds::ISimplexMeshIOService ioServiceMesh(this);
                        gmds::VTKWriter vtkWriterMeshER(&ioServiceMesh);
                        vtkWriterMeshER.setCellOptions(gmds::N|gmds::R|gmds::F);
                        vtkWriterMeshER.setDataOptions(gmds::N|gmds::R|gmds::F);
                        vtkWriterMeshER.write("ADJACENT_PROBLEM.vtk");
                        throw gmds::GMDSException("WRONG MESH");
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
}
/******************************************************************************/
void SimplexMesh::checkMeshCavity(const std::vector<TSimplexID>& cavity)
{
  TInt border = std::numeric_limits<int>::min();
  unsigned int sizeFACE = 4;
  unsigned int sizeEDGE = 3;

  for(auto const s : cavity)
  {
    if(s >= 0)
    {
      const SimplicesCell c0 = SimplicesCell(this, s);
      for(unsigned int face = 0; face < sizeFACE; face++)
      {
        TSimplexID oppositeCell = m_tet_adj[s][face];
        if(oppositeCell == border)
        {
          std::cout << "c0 -> " << c0 << std::endl;
          throw gmds::GMDSException("oppositeCell == border");
        }

        if(oppositeCell >= 0)
        {
          const SimplicesCell c1 = SimplicesCell(this, oppositeCell);
          std::vector<TInt> nodesFace =  c0.getOrderedFace(face);
          if(nodesFace.size() != 3)
          {
            throw gmds::GMDSException("nodeFace.size() != 3");
          }

          std::vector<TInt> otherNode = c1.getOtherNodeInSimplex(nodesFace);
          if(otherNode.size() != 1)
          {
            throw gmds::GMDSException("otherNode.size() != 1");
          }

          if(m_tet_adj[oppositeCell][c1.getLocalNode(otherNode.front())] != s)
          {
            std::cout << "TET WARNING" << std::endl;
            std::cout << "simplex -> " << c0 << std::endl;
            std::cout << "opposite Cell -> " << c1 << std::endl;
            throw gmds::GMDSException("m_tet_adj[oppositeCell][otherNode.size()] != simplex");
          }
        }
      }
    }
    else
    {
      //tetra neigbor of the current triangle s
      TSimplexID oppositeCell = m_tri_nodes[-s][3];
      if(oppositeCell < 0)
      {
        std::cout << "TRIANGLE WARNING" << std::endl;
        throw gmds::GMDSException("oppositeCell < 0");
      }
      const SimplicesCell c1 = SimplicesCell(this, oppositeCell);
      std::vector<TInt> nodesFace = SimplicesTriangle(this, s).getNodes();
      std::vector<TInt> otherNode = c1.getOtherNodeInSimplex(nodesFace);
      if(m_tet_adj[oppositeCell][c1.getLocalNode(otherNode.front())] != s)
      {
        std::cout << "TRIANGLE WARNING" << std::endl;
        std::cout << "triangle -> " << s << std::endl;
        std::cout << "opposite Cell -> " << c1 << std::endl;
        throw gmds::GMDSException("m_tet_adj[oppositeCell][otherNode.size()] != simplex");
      }

      //triagle neigbor
      for(unsigned int edge = 0; edge < sizeEDGE; edge++)
      {
        TSimplexID oppositeCell = m_tri_adj[-s][edge];
        if(oppositeCell == border)
        {
          std::cout << "TRIANGLE WARNING" << std::endl;
          std::cout << "s -> " << s << std::endl;
          throw gmds::GMDSException("oppositeCell == border");
        }

        std::vector<TInt> nodesEdge =  SimplicesTriangle(this, s).getOppositeEdge(edge);
        const SimplicesTriangle c1 = SimplicesTriangle(this, oppositeCell);

        std::vector<TInt> otherNode = c1.getOtherNodeInSimplex(nodesEdge);
        if(otherNode.size() != 1)
        {
          std::cout << "TRIANGLE WARNING" << std::endl;
          throw gmds::GMDSException("otherNode.size() != 1");
        }

        if(m_tri_adj[oppositeCell][c1.getLocalNode(otherNode.front())] != -s)
        {
          std::cout << "TRIANGLE WARNING" << std::endl;
          std::cout << "triangle -> " << s << std::endl;
          std::cout << "oppositeCell -> " << oppositeCell << std::endl;
          std::cout << "m_tri_adj[oppositeCell][c1.getLocalNode(otherNode.front())] -> " << m_tri_adj[oppositeCell][c1.getLocalNode(otherNode.front())] << std::endl;
          throw gmds::GMDSException("m_tet_adj[oppositeCell][otherNode.size()] != simplex");
        }
      }
    }
  }
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

  Variable<int>* BND_VERTEX_COLOR  = nullptr;
  Variable<int>* BND_CURVE_COLOR   = nullptr;
  Variable<int>* BND_SURFACE_COLOR = nullptr;

  try{
    BND_VERTEX_COLOR  = getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
    BND_CURVE_COLOR   = getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
    BND_SURFACE_COLOR = getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
  } catch(GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  /*this loop find all the border tetraedron with the local node pointing the border of the mesh*/
  unsigned int sizeFace = 4;
  for(unsigned int tet = 0 ; tet < m_tet_ids.capacity() ; tet++)
  {
    if(m_tet_ids[tet] != 0)
    {
      for(unsigned int localNode = 0 ; localNode < sizeFace ; localNode++)
      {
        TSimplexID oppositeSimplex = m_tet_adj[tet][localNode];
        if(oppositeSimplex == border)
        {
          TInt node0 = m_tet_nodes[tet][(localNode + 1) % sizeFace];
          TInt node1 = m_tet_nodes[tet][(localNode + 2) % sizeFace];
          TInt node2 = m_tet_nodes[tet][(localNode + 3) % sizeFace];
          //This condition won't build a null triangle on a curve (if a plate tetra exist in the mesh 's border)
          if((*BND_CURVE_COLOR)[node0] != 0 && (*BND_CURVE_COLOR)[node0] == (*BND_CURVE_COLOR)[node1]  && (*BND_CURVE_COLOR)[node0] == (*BND_CURVE_COLOR)[node2])
          {
            //continue;
          }
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


  gmds::Variable<int>* BND_TRIANGLES = newVariable<int, SimplicesTriangle>("BND_TRIANGLES");
  std::vector<TInt> notIndexedCornerTriangle{};
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

      if(dimMax != VOLUME)
      {
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
            //look for the adjacent surface triangle ---> to change if some triangle will be added on the volume mesh
            TSimplexID adjTriangle = triangle.neighborTriangle(localNode);
            std::vector<TInt> nodesAdjTriangle = SimplicesTriangle(this, adjTriangle).otherNodesInTriangle(triangle);
            if(nodesAdjTriangle.size() == 1)
            {
              index = (*BND_SURFACE_COLOR)[nodesAdjTriangle.front()];

              if(index != 0)
              {
                (*BND_TRIANGLES)[tri] = index;
              }
              else
              {
                notIndexedCornerTriangle.push_back(tri);
              }
            }
          }
        }
      }
      else
      {
        gmds::ISimplexMeshIOService ioServiceMesh(this);
        std::vector<TSimplexID> v{static_cast<int>(tri)};
        deleteAllTrianglesBut(v);
        gmds::VTKWriter vtkWriterMeshER(&ioServiceMesh);
        vtkWriterMeshER.setCellOptions(gmds::N|gmds::F);
        vtkWriterMeshER.setDataOptions(gmds::N|gmds::F);
        vtkWriterMeshER.write("TRIANGLE_INSIDE_VOLUME.vtk");
        std::cout << SimplicesTriangle(this, tri) << std::endl;
        throw gmds::GMDSException("dimMax == VOLUME, triangle is inside the volume");
      }
    }
  }

  for(auto const triangle : notIndexedCornerTriangle)
  {
    int index;
    TSimplexID tri = triangle;
    unsigned int sizeFace = 3;
    gmds::BitVector cyclingCheck(m_tri_ids.capacity());
    cyclingCheck.assign(tri);
    for(;;)
    {
      if((*BND_TRIANGLES)[tri] != 0)
      {
        index = (*BND_TRIANGLES)[tri];
        break;
      }
      std::vector<TInt> nodes{m_tri_nodes[tri][0], m_tri_nodes[tri][1], m_tri_nodes[tri][2]};
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
              tri = adjTriangle;
              break;
            }
          }
        }
      }
    }
    (*BND_TRIANGLES)[triangle] = index;
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
        } while(idx == 0 || idx > 200);

        (*BND_TRIANGLES)[tri] = idx;
      }
    }
  }

  //filling edgeTianglesIndices & cornerSurfaceConnexion & cornerEdgeConnexion map
  auto my_make = [=](const TInt a, const TInt b){
    return (a < b)? std::make_pair(a, b) : std::make_pair(b, a);
  };
  for(unsigned int nodeIdx = 0 ; nodeIdx < m_node_ids.capacity() ; nodeIdx++)
  {
    if(m_node_ids[nodeIdx] != 0)
    {
      if((*BND_CURVE_COLOR)[nodeIdx] != 0)
      {
        unsigned int indiceNode      = (*BND_CURVE_COLOR)[nodeIdx];
        std::vector<TSimplexID> ball = SimplicesNode(this, nodeIdx).ballOf();
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
            std::pair<unsigned int, unsigned int> pairTrianglesIndedices = my_make(triangleId0, triangleId1);
            edgeTianglesIndices[indiceNode] = pairTrianglesIndedices;
          }
          else
          {
            std::cout << "trianglesIndices.size() != 2 --> " << trianglesIndices.size() << " for node : "<< nodeIdx << " of index : " << indiceNode <<  std::endl;
            //throw GMDSException("problem with the mesh file (the edge is maybe not realy one)");
          }
        }
      }
      else if((*BND_VERTEX_COLOR)[nodeIdx] != 0)
      {
        unsigned int indiceNode      = (*BND_VERTEX_COLOR)[nodeIdx];
        std::vector<TSimplexID> ball = SimplicesNode(this, nodeIdx).ballOf();
        std::set<unsigned int> trianglesIndices{};
        if(cornerSurfaceConnexion.find(indiceNode) == cornerSurfaceConnexion.end())
        {
          for(auto const simplex : ball)
          {
            if(simplex < 0 && simplex != border)
            {
              unsigned int triangleIndice = (*BND_TRIANGLES)[-simplex];
              trianglesIndices.insert(triangleIndice);
            }
          }

          std::vector<unsigned int> trianglesIndicesVec{};
          std::copy(trianglesIndices.begin(), trianglesIndices.end(), std::back_inserter(trianglesIndicesVec));
          cornerSurfaceConnexion[indiceNode] = trianglesIndicesVec;
        }
        else
        {
          std::cout << "It exists many corner with the same index | corner node -> "<<  nodeIdx << " of index -> " <<(*BND_VERTEX_COLOR)[nodeIdx]  << std::endl;
          throw gmds::GMDSException("");
        }
      }
    }
  }

  for(unsigned int nodeIdx = 0 ; nodeIdx < m_node_ids.capacity() ; nodeIdx++)
  {
    if(m_node_ids[nodeIdx] != 0)
    {
      if((*BND_VERTEX_COLOR)[nodeIdx] != 0)
      {
        std::vector<unsigned int> indicesEdges{};
        unsigned int indiceNode      = (*BND_VERTEX_COLOR)[nodeIdx];
        std::vector<unsigned int> indicesSurfaces = cornerSurfaceConnexion.at(indiceNode);

        for(auto const & edgeTriangle : edgeTianglesIndices)
        {
          unsigned int edgeIdx = edgeTriangle.first;
          unsigned int surfaceIdx0 = edgeTriangle.second.first;
          unsigned int surfaceIdx1 = edgeTriangle.second.second;

          unsigned int SubSurfaceConnexionCpt = 0;
          for(auto const indiceSurface : indicesSurfaces)
          {
            if(indiceSurface == surfaceIdx0 || indiceSurface == surfaceIdx1)
            {
              SubSurfaceConnexionCpt++;
            }
          }

          if(SubSurfaceConnexionCpt == 2)
          {
            indicesEdges.push_back(edgeIdx);
          }
        }
        cornerEdgeConnexion[indiceNode] = indicesEdges;
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
  unsigned int sizeFace = 4;
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
                bool flag = false;
                for(unsigned int face = 0 ; face < sizeFace ; face++)
                {
                  if(m_tet_adj[simplexAdj][face] == ATetraIndex)
                  {
                    flag = true;
                    m_tet_adj[simplexAdj][face] = errorId;
                    break;
                  }
                }
                if(!flag)
                {
                  gmds::ISimplexMeshIOService ioServiceMesh(this);
                  std::vector<TSimplexID> v{ATetraIndex, simplexAdj};
                  deleteAllSimplicesBut(v);
                  gmds::VTKWriter vtkWriterMeshER(&ioServiceMesh);
                  vtkWriterMeshER.setCellOptions(gmds::N|gmds::R|gmds::F);
                  vtkWriterMeshER.setDataOptions(gmds::N|gmds::R|gmds::F);
                  vtkWriterMeshER.write("ADJACENT_PROBLEM.vtk");
                  std::cout << ATetraIndex << " is not adjacent to " << simplexAdj << std::endl;
                  throw gmds::GMDSException("adjacent relation problem !");
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
        math::Vector3d normal({0.0f, 0.0f, 0.0f});
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
        math::Vector3d normal({0.0, 0.0, 0.0});
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
        math::Vector3d normal({0.0, 0.0, 0.0});
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
					math::Vector3d vtemp=m_coords[neighbor] - m_coords[node];
              math::Point recenterNeighborPoint(vtemp.X(),vtemp.Y(),vtemp.Z());
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
  //uvwt.reserve(4);


  if(m_tet_ids.size() == 0){return false;}
  if(m_tri_ids.size() == 0){return false;}
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
  //std::cout << "  currentSimplex -> " << SimplicesCell(this,currentSimplex) << std::endl;
 TSimplexID cycling = std::numeric_limits<int>::min();
 TSimplexID border  = cycling;
 TSimplexID nextTet = currentSimplex;

 uvwt.clear();
 uvwt.push_back(SimplicesCell(this, currentSimplex).orientation(0, pt));
 uvwt.push_back(SimplicesCell(this, currentSimplex).orientation(1, pt));
 uvwt.push_back(SimplicesCell(this, currentSimplex).orientation(2, pt));
 uvwt.push_back(SimplicesCell(this, currentSimplex).orientation(3, pt));
 //std::cout << "  uvwt -> " << uvwt[0] << " | " << uvwt[1] << " | " << uvwt[2] << " | " << uvwt[3] << " | " << std::endl;

 std::vector<TInt> nodes = SimplicesCell(this, currentSimplex).getNodes();

 if(cyclingCheck[currentSimplex] == 0)
 {
   if(uvwt[0] > 0  && uvwt[1]  > 0 && uvwt[2]  > 0 && uvwt[3]  > 0)
   {
     return currentSimplex;
   }

   unsigned int cpt = 0;
   std::vector<unsigned int> adjTetra{};
   std::vector<unsigned int> adjTri{};
   SimplicesCell cell = SimplicesCell(this, currentSimplex);
   for(unsigned int idx = 0; idx < uvwt.size() ; idx++)
   {
     math::Orientation::Sign bar = uvwt[idx];
     if(bar < 0)
     {
       TSimplexID adjCell = cell.oppositeTetraIdx(idx);
       if(adjCell >= 0)
       {
         adjTetra.push_back(adjCell);
       }
       else
       {
         adjTri.push_back(adjCell);
       }
     }
   }

   if(adjTetra.size() != 0)
   {
     nextTet = adjTetra.front();
   }
   else
   {
     if(adjTri.size() != 0)
     {
       nextTet = adjTri.front();
     }
     else
     {
       //std::cout << "adjTri.size() == 0" << std::endl;
     }
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
    math::Vector3d ray_origin    = math::vec(nodeA.getCoords());

    //on regarde maintenant si il y'a intersection... Möller..
    math::Vector3d E1 = SimplicesNode(this, triangleNodeId[1]).getCoords() - SimplicesNode(this, triangleNodeId[0]).getCoords();
    math::Vector3d E2 = SimplicesNode(this, triangleNodeId[2]).getCoords() - SimplicesNode(this, triangleNodeId[0]).getCoords();
    math::Vector3d T = ray_origin - math::vec(SimplicesNode(this, triangleNodeId[0]).getCoords());


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
    tuv = (1.0/den) * math::Vector3d({compo1, compo2, compo3});


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
void SimplexMesh::setBase(const TInt node, const TSimplexID simplex)
{
  if(this != nullptr)
  {
    if(m_base.size() >= node -1)
    {
      m_base[node] = simplex;
    }
    else
    {
      std::cout << "m_base.size() >= node -1 -> " << (m_base.size() >= node -1) << std::endl;
      std::cout << "node -> " << node - 1<< std::endl;
      std::cout << "m_base.size() -> " << m_base.size() << std::endl;
      throw gmds::GMDSException("m_base.size() < node -1");
    }
  }
}
/******************************************************************************/
void SimplexMesh::getEdgeSizeInfo(double& meanEdges, double& maxEdge, double& minEdge)
{
  auto my_make = [=](const TInt a, const TInt b){
    return (a < b)? std::make_pair(a, b) : std::make_pair(b, a);
  };

  meanEdges = 0.0;
  std::map<TInt, TInt> edges{};
  std::set<double> sizeEdges{};
  for(unsigned int node = 0 ; node < m_node_ids.capacity() ; node++)
  {
    if(m_node_ids[node] != 0)
    {
      const SimplicesNode N(this, node);
      std::vector<TInt> nodes = N.neighborNodes();
      for(auto const n : nodes)
      {
        std::pair<TInt, TInt> edge = my_make(node, n);
        edges.insert(edge);
      }
    }
  }

  for(auto const edge : edges)
  {
    const SimplicesNode N0(this, edge.first);
    const SimplicesNode N1(this, edge.second);

    const math::Point coord0 = N0.getCoords();
    const math::Point coord1 = N1.getCoords();

    const math::Vector3d vec = coord0 - coord1;
    const double norm = vec.norm();
    sizeEdges.insert(norm);
  }

  for(auto const sizeEdge : sizeEdges)
  {
    meanEdges += sizeEdge;
  }

  meanEdges /= sizeEdges.size();
  maxEdge = *(sizeEdges.begin());
  minEdge = *(--sizeEdges.end());

  myfile.open("data_AdaptationMesh.txt", std::ios_base::app);
  if (myfile.is_open())
  {
    myfile << "meanValue -> " << meanEdges << std::endl;
    myfile << "maxEdge -> " << maxEdge << std::endl;
    myfile << "minEdge -> " << minEdge << std::endl;
    std::cout << std::endl;
    myfile.close();
  }
  else
  {
    throw gmds::GMDSException("Unable to open file");
  }
}
/******************************************************************************/
void SimplexMesh::getEdgeSizeInfowithMetric(double& meanEdges, double& minEdge, double& maxEdge)
{
  Variable<Eigen::Matrix3d>* metric = nullptr;
  try{
    metric = getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

    auto my_make = [=](const TInt a, const TInt b){
    return (a < b)? std::make_pair(a, b) : std::make_pair(b, a);
  };

  meanEdges = 0.0;
  std::map<TInt, TInt> edges{};
  std::set<double> sizeEdges{};
  for(unsigned int node = 0 ; node < m_node_ids.capacity() ; node++)
  {
    if(m_node_ids[node] != 0)
    {
      const SimplicesNode N(this, node);
      std::vector<TInt> nodes = N.neighborNodes();
      for(auto const n : nodes)
      {
        std::pair<TInt, TInt> edge = my_make(node, n);
        edges.insert(edge);
      }
    }
  }

  unsigned int edgesAboveSqrt2CPT = 0;
  unsigned int edgesUnderSqrt_2CPT = 0;
  unsigned int goodSizeEdge_CPT = 0;

  for(auto const edge : edges)
  {
    const SimplicesNode N0(this, edge.first);
    const SimplicesNode N1(this, edge.second);

    const math::Point coord0 = N0.getCoords();
    const math::Point coord1 = N1.getCoords();

    Metric<Eigen::Matrix3d> M0 = Metric<Eigen::Matrix3d>((*metric)[edge.first]);
    Metric<Eigen::Matrix3d> M1 = Metric<Eigen::Matrix3d>((*metric)[edge.second]);

    double metricLenght        = M0.metricDist(vec(coord0), vec(coord1), M1);
    if(metricLenght > sqrt(2))
    {
      edgesAboveSqrt2CPT++;
    }
    else if(metricLenght < sqrt(2) * 0.5)
    {
      edgesUnderSqrt_2CPT++;
    }
    else
    {
      goodSizeEdge_CPT++;
    }
    sizeEdges.insert(metricLenght);
  }

  for(auto const sizeEdge : sizeEdges)
  {
    meanEdges += sizeEdge;
  }

  meanEdges /= sizeEdges.size();
  minEdge = *(sizeEdges.begin());
  maxEdge = *(--sizeEdges.end());

  myfile.open("data_AdaptationMesh.txt", std::ios_base::app);
  if (myfile.is_open())
  {
    myfile << "meanValue -> " << meanEdges << std::endl;
    myfile << "maxEdge -> " << maxEdge << std::endl;
    myfile << "minEdge -> " << minEdge << std::endl;
    myfile << "edgesUnderSqrt_2CPT -> " << edgesUnderSqrt_2CPT / static_cast<double>(edges.size())<< std::endl;
    myfile << "edgesAboveSqrt2CPT -> " << edgesAboveSqrt2CPT / static_cast<double>(edges.size())<< std::endl;
    myfile << "edgeWithWrongSize -> " << (edgesUnderSqrt_2CPT + edgesAboveSqrt2CPT) / static_cast<double>(edges.size()) << std::endl;
    myfile << "goodSizeEdge_CPT -> " << goodSizeEdge_CPT / static_cast<double>(edges.size()) << std::endl;
    myfile.close();
  }
  else
  {
    throw gmds::GMDSException("Unable to open file");
  }
}
/******************************************************************************/
Eigen::Matrix3d SimplexMesh::getAnalyticMetric(const Point& pt) const
{
  double epsilon = 0.01;
  Eigen::Matrix3d m = Eigen::MatrixXd::Identity(3, 3);

  //CONSTANT ISOTROPE METRIC
  double metricX = 0.05*(1.0 - pt.X()) + 0.1*pt.X();
  double metricY = 0.05*(1.0 - pt.X()) + 0.1*pt.X();
  double metricZ = 0.05*(1.0 - pt.X()) + 0.1*pt.X();

  /*if(pt.Y() <= 0.5)
  {
    metricX = 0.1;
    metricY = 0.1;
    metricZ = 0.1;
  }
  else
  {
    metricX = 0.2;
    metricY = 0.2;
    metricZ = 0.2;
  }*/


  /*metricX = 0.01*(1.0 - pt.Y()) + 0.1*pt.Y();
  metricY = 0.01*(1.0 - pt.Y()) + 0.1*pt.Y();
  metricZ = 0.01*(1.0 - pt.Y()) + 0.1*pt.Y();*/

  //double metricX = std::atan(80.0* (std::pow(pt.X(), 4) /** std::pow(pt.Y(), 4)*/));
  //double metricY = std::atan(80.0* (std::pow(pt.X(), 4) /** std::pow(pt.Y(), 4)*/));
  //double metricZ = std::atan(80.0* (std::pow(pt.X(), 4) /** std::pow(pt.Y(), 4)*/));

  /*double metricX = (1.0 - std::exp(-1.0*(std::pow(pt.X() - 0.5, 2) + std::pow(pt.Y() - 0.5, 2)))) + epsilon;
  double metricY = (1.0 - std::exp(-1.0*(std::pow(pt.X() - 0.5, 2) + std::pow(pt.Y() - 0.5, 2)))) + epsilon;
  double metricZ = (1.0 - std::exp(-1.0*(std::pow(pt.X() - 0.5, 2) + std::pow(pt.Y() - 0.5, 2)))) + epsilon;*/


  m(0,0) = 1.0 / (metricX*metricX);
  m(1,1) = 1.0 / (metricY*metricY);
  m(2,2) = 1.0 / (metricZ*metricZ);
  return m;

}
/******************************************************************************/
void SimplexMesh::setAnalyticMetric(const TInt node)
{
  //std::cout << "NODE  -> " << node << std::endl;
  Variable<Eigen::Matrix3d>* metric = nullptr;
  try{
    metric = getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  //analytic isotrope metric here !
  gmds::math::Point pt = m_coords[node];
  double epsilon = 0.01;

  //CONSTANT ISOTROPE METRIC
  (*metric)[node] =  Eigen::MatrixXd::Identity(3, 3);
  double metricX = 0.05*(1.0 - pt.X()) + 0.1*pt.X();
  double metricY = 0.05*(1.0 - pt.X()) + 0.1*pt.X();
  double metricZ = 0.05*(1.0 - pt.X()) + 0.1*pt.X();
  /*std::cout << "pt -> " << pt << std::endl;
  std::cout << "metricX -> " << metricX << std::endl;
  std::cout << "metricY -> " << metricY << std::endl;
  std::cout << "metricZ -> " << metricZ << std::endl;*/



  /*if(pt.Y() <= 0.5)
  {
    metricX = 0.1;
    metricY = 0.1;
    metricZ = 0.1;
  }
  else
  {
    metricX = 0.2;
    metricY = 0.2;
    metricZ = 0.2;
  }*/

  //double metricX = std::atan(80.0* (std::pow(pt.X(), 4) /** std::pow(pt.Y(), 4)*/));
  //double metricY = std::atan(80.0* (std::pow(pt.X(), 4) /** std::pow(pt.Y(), 4)*/));
  //double metricZ = std::atan(80.0* (std::pow(pt.X(), 4) /** std::pow(pt.Y(), 4)*/));

  /*double metricX = (1.0 - std::exp(-1.0*(std::pow(pt.X() - 0.5, 2) + std::pow(pt.Y() - 0.5, 2)))) + epsilon;
  double metricY = (1.0 - std::exp(-1.0*(std::pow(pt.X() - 0.5, 2) + std::pow(pt.Y() - 0.5, 2)))) + epsilon;
  double metricZ = (1.0 - std::exp(-1.0*(std::pow(pt.X() - 0.5, 2) + std::pow(pt.Y() - 0.5, 2)))) + epsilon;*/

  (*metric)[node](0,0) = 1.0 / (metricX*metricX);
  (*metric)[node](1,1) = 1.0 / (metricY*metricY);
  (*metric)[node](2,2) = 1.0 / (metricZ*metricZ);
  //std::cout << "*metric -> " << (*metric)[node] << std::endl;
}
/******************************************************************************/
void SimplexMesh::rebuildCav(CavityOperator::CavityIO& cavityIO, const std::vector<std::vector<TInt>>& deleted_Tet, const std::vector<std::vector<TInt>>& deleted_Tri,  const TInt nodeToConnect, std::vector<TSimplexID>& createdCells)
{
  gmds::Variable<int>* BND_SURFACE_COLOR = nullptr;
  gmds::Variable<int>* BND_CURVE_COLOR = nullptr;
  gmds::Variable<int>* BND_VERTEX_COLOR = nullptr;
  gmds::Variable<int>* BND_TRIANGLES   = nullptr;

  try{
    BND_SURFACE_COLOR    = getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
    BND_CURVE_COLOR      = getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
    BND_VERTEX_COLOR     = getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
    BND_TRIANGLES        = getVariable<int, SimplicesTriangle>("BND_TRIANGLES");
  } catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }
  const std::vector<std::vector<TInt>>&               pointsToConnect              = cavityIO.getNodesToReconnect();
  const std::vector<TSimplexID>&                      extSimplexBorder             = cavityIO.getOppositesSimplex();
  const std::vector<std::vector<TInt>>&               pointsToConnect_tri          = cavityIO.getNodesToReconnect_Tri();
  const std::vector<TSimplexID>&                      extSimplexBorder_tri         = cavityIO.getOppositesSimplex_Tri();
  const std::vector<unsigned int>&                    trianglesIndices             = cavityIO.getIndices_Tri();
  const std::map<std::pair<TInt, TInt>, TSimplexID> & mapForReinsertion            = cavityIO.getReinsertionData();
  const std::map<TInt, TSimplexID> &            mapForReinsertionSurface           = cavityIO.getReinsertionSurfaceData();

  std::map<std::pair<TInt, TInt>, std::pair<unsigned int, TSimplexID>>  mapForTrianglesColor         = cavityIO.getTrianglesColor();

  TInt border = std::numeric_limits<int>::min();
  unsigned int sizeFace = 3;
  //Rebuild volume
  std::multimap<std::pair<TInt, TInt>, TSimplexID> hash_face; // multimap  for the reconnection of the tetraedron
  std::multimap<std::pair<TInt, TInt>, TSimplexID> hash_face_Copy; // multimap  for the reconnection of the tetraedron

  auto my_make = [=](const TInt a, const TInt b){
    return (a < b)? std::make_pair(a, b) : std::make_pair(b, a);
  };

  std::vector<TSimplexID> tets{};
  std::vector<TSimplexID> trianglesCreated;
  std::vector<TSimplexID> tetsCreated;
  for(unsigned int idx = 0 ; idx < pointsToConnect.size() ; idx++)
  {
    //std::cout << "idx -> " << idx << std::endl;
    //extracting cell data for reconstruction
    std::vector<TInt> nodes       = pointsToConnect[idx];
    TSimplexID oppositeSimplex    = extSimplexBorder[idx];

    //sort the edge's face of the simpelx created for the adjacent reconstruction
    TInt nodeA = nodes[0]; TInt nodeB = nodes[1]; TInt nodeC = nodes[2];

    std::pair<TInt, TInt> p0 = my_make(nodeA, nodeB);
    std::pair<TInt, TInt> p1 = my_make(nodeB, nodeC);
    std::pair<TInt, TInt> p2 = my_make(nodeC, nodeA);

    TSimplexID cell = addTetraedre(nodeA, nodeB, nodeC, nodeToConnect, false);
    createdCells.push_back(cell);
    tetsCreated.push_back(cell);

    tets.push_back(cell);
    hash_face.insert(std::pair<std::pair<TInt, TInt>, TSimplexID>(p0, cell));
    hash_face.insert(std::pair<std::pair<TInt, TInt>, TSimplexID>(p1, cell));
    hash_face.insert(std::pair<std::pair<TInt, TInt>, TSimplexID>(p2, cell));

    std::vector<TSimplexID> vecAdj{border, border, border, border};
    m_tet_adj[cell] = vecAdj;
    m_tet_adj[cell][3] = oppositeSimplex;

    if(oppositeSimplex != border)
    {
      if(oppositeSimplex >= 0)
      {
        std::vector<TInt>&& nodesNeighborTet = SimplicesCell(this, oppositeSimplex).getOtherNodeInSimplex(nodes);
        if(nodesNeighborTet.size() == 1)
        {
          TInt nodeLocalNeighborSimplex = SimplicesCell(this, oppositeSimplex).getLocalNode(nodesNeighborTet.front());
          m_tet_adj[oppositeSimplex][nodeLocalNeighborSimplex] = cell;
        }
        else
        {
          throw gmds::GMDSException("nodesNeighborTet.size() == 1");
        }
      }
      else
      {
        m_tri_nodes[-oppositeSimplex][3] = cell;
      }
    }
  }

  hash_face_Copy = hash_face;
  while(!hash_face.empty())
  {
    std::multimap<std::pair<TInt, TInt>, TSimplexID>::iterator it0 = hash_face.begin();
    std::pair<std::pair<TInt, TInt>, TSimplexID> p = *it0;
    std::pair<TInt, TInt> pairNode = my_make(p.first.first, p.first.second); //std::pair(nodeA, nodeB)
    std::pair <std::multimap<std::pair<TInt, TInt>, TSimplexID>::iterator, std::multimap<std::pair<TInt, TInt>, TSimplexID>::iterator> pair_it;
    pair_it = hash_face.equal_range(pairNode);
    int count = std::distance(pair_it.first, pair_it.second);

    if(count > 2)
    {
      std::cout << "hash_face filling algo" << std::endl;
      std::cout << "p.first.first -> " << p.first.first << " | " << "p.first.second -> " << p.first.second << std::endl;

      deleteAllSimplicesBut(tetsCreated);
      gmds::ISimplexMeshIOService ioServiceMesh(this);
      gmds::VTKWriter vtkWriterMesh(&ioServiceMesh);
      vtkWriterMesh.setCellOptions(gmds::N|gmds::R|gmds::F);
      vtkWriterMesh.setDataOptions(gmds::N|gmds::R|gmds::F);
      vtkWriterMesh.write("DEBUG_MESH.vtk");
      std::cout << "count -> " << count << std::endl;
      throw gmds::GMDSException("count > 2");
    }
    else if(count == 2)
    {
      std::multimap<std::pair<TInt, TInt>, TSimplexID>::iterator it1 = hash_face.begin();
      ++it1;
      TSimplexID simplexA = it0->second;
      TSimplexID simplexB = it1->second;

      const SimplicesCell cA = SimplicesCell(this, simplexA);
      const SimplicesCell cB = SimplicesCell(this, simplexB);

      std::vector<TInt> intersectionFaces = cA.intersectionNodes(cB);
      if(intersectionFaces.size() != 3)
      {
        std::cout << "cellA -> " << cA << std::endl;
        std::cout << "cellB -> " << cB << std::endl;
        gmds::ISimplexMeshIOService ioServiceMesh(this);
        std::vector<TSimplexID> v{simplexA, simplexB};
        deleteAllSimplicesBut(v);
        gmds::VTKWriter vtkWriterMeshER(&ioServiceMesh);
        vtkWriterMeshER.setCellOptions(gmds::N|gmds::R|gmds::F);
        vtkWriterMeshER.setDataOptions(gmds::N|gmds::R|gmds::F);
        vtkWriterMeshER.write("INTERSECTION_FACE_WRONG.vtk");
        std::cout << "intersectionFaces.size() = "  << intersectionFaces.size() << std::endl;
        throw gmds::GMDSException("intersectionFaces.size() != 3");
      }

      std::vector<TInt> otherNode_A = cA.getOtherNodeInSimplex(intersectionFaces);
      std::vector<TInt> otherNode_B = cB.getOtherNodeInSimplex(intersectionFaces);

      if(otherNode_A.size() != 1 || otherNode_B.size() != 1)
      {
        throw gmds::GMDSException("otherNode_A != 1 || otherNode_B != 1");
      }

      TInt localNode_A = cA.getLocalNode(otherNode_A.front());
      TInt localNode_B = cB.getLocalNode(otherNode_B.front());

      m_tet_adj[simplexA][localNode_A] = simplexB;
      m_tet_adj[simplexB][localNode_B] = simplexA;

      hash_face.erase(it0);
      hash_face.erase(it1);
    }
    else if(count == 1)
    {
      //reinsertion node
      const TSimplexID simplexA = p.second;
      TInt nodeA = p.first.first ; TInt nodeB = p.first.second ; TInt nodeC = nodeToConnect ;
      std::vector<TInt> nodesFace{nodeA, nodeB, nodeC};

      std::pair<TInt, TInt> p_data = my_make(nodeA, nodeB);

      TSimplexID simplexAdj = border;
      try{
        simplexAdj = mapForReinsertion.at(p_data);
      } catch(std::out_of_range& oor)
      {
        //if the surface is being reconstruct it's normal look in cavityOperator enlargement function
        //mapForReinsertion is not fill if oppositeCell is a triangle
        //if not throw exception
        if((*BND_SURFACE_COLOR)[nodeToConnect] == 0 && (*BND_CURVE_COLOR)[nodeToConnect] == 0 && (*BND_VERTEX_COLOR)[nodeToConnect] == 0)
        {
          std::cout << "nodeToConnect -> " << nodeToConnect << std::endl;
          std::cout << "p_data -> " << p_data.first << " | " << p_data.second <<std::endl;
          std::cout << "simplexA ->" << simplexA << std::endl;

          std::vector<TSimplexID> simplices{simplexA};
          deleteAllSimplicesBut(simplices);

          gmds::ISimplexMeshIOService ioServiceMesh(this);
          gmds::VTKWriter vtkWriterMesh(&ioServiceMesh);
          vtkWriterMesh.setCellOptions(gmds::N|gmds::R|gmds::F);
          vtkWriterMesh.setDataOptions(gmds::N|gmds::R|gmds::F);
          vtkWriterMesh.write("DEBUG_MESH.vtk");

          std::cerr << "Out of Range error: " << oor.what() << '\n';
        }
        else
        {
          hash_face.erase(it0);
          continue;
        }
      }

      const SimplicesCell cellAdj = SimplicesCell(this, simplexAdj);
      const SimplicesCell cell    = SimplicesCell(this, simplexA);

      std::vector<TInt> otherNode_Adj = cellAdj.getOtherNodeInSimplex(nodesFace);
      std::vector<TInt> otherNode_A   = cell.getOtherNodeInSimplex(nodesFace);
      if(otherNode_Adj.size() != 1 || otherNode_A.size() != 1)
      {
        throw gmds::GMDSException(" ! (otherNode_Adj.size() != 1 || otherNode_A.size() != 1)");
      }

      TInt localNode_Adj    = cellAdj.getLocalNode(otherNode_Adj.front());
      TInt localNode_cell   = cell.getLocalNode(otherNode_A.front());

      m_tet_adj[simplexA][localNode_cell] = simplexAdj;
      m_tet_adj[simplexAdj][localNode_Adj] = simplexA;

      hash_face.erase(it0);
    }
  }

  //REBULDING THE SURFACE MESH
  std::multimap<std::pair<TInt, TInt>, TSimplexID> hash_edge{};
  while(!mapForTrianglesColor.empty())
  {
    std::multimap<std::pair<TInt, TInt>, std::pair<unsigned int, TSimplexID>>::iterator it0 = mapForTrianglesColor.begin();
    std::pair<std::pair<TInt, TInt>, std::pair<unsigned int, TSimplexID>> p = *it0;
    TInt nodeA  = p.first.first;
    TInt nodeB  = p.first.second;

    unsigned int indice         = p.second.first;
    TSimplexID oppositeTriangle = p.second.second;
    std::pair<TInt, TInt> pairNode = my_make(nodeA, nodeB);

    TSimplexID T = border;
    std::multimap<std::pair<TInt, TInt>, TSimplexID>::iterator it = hash_face_Copy.find(pairNode);
    if(it == hash_face_Copy.end())
    {
      std::cout << "node to connect -> " << nodeToConnect << std::endl;
      std::cout << "pairNode.first -> " << pairNode.first << std::endl;
      std::cout << "pairNode.second -> " << pairNode.second << std::endl;

      deleteAllSimplicesBut(tetsCreated);
      gmds::ISimplexMeshIOService ioServiceMesh(this);
      gmds::VTKWriter vtkWriterMesh(&ioServiceMesh);
      vtkWriterMesh.setCellOptions(gmds::N|gmds::R|gmds::F);
      vtkWriterMesh.setDataOptions(gmds::N|gmds::R|gmds::F);
      vtkWriterMesh.write("DEBUG_MESH.vtk");

      throw gmds::GMDSException("hash_face_Copy.find(pairNode) == hash_face_Copy.end()");
    }
    else
    {
      T = it->second;
    }

    //sort the node in order to build a correct (not invert) triangle for the normal
    std::vector<TInt> nodes = SimplicesTriangle(this, oppositeTriangle).getNodes();
    std::vector<TInt> orderedNodes{};
    for(unsigned int ln = 0 ; ln < sizeFace ; ++ln)
    {
      if(nodes[ln] != nodeA && nodes[ln] != nodeB)
      {
        orderedNodes.push_back(nodes[(ln + 1) % sizeFace]);
        orderedNodes.push_back(nodes[(ln + 2) % sizeFace]);
        break;
      }
    }
    if(orderedNodes.size() != 2)
    {
      throw gmds::GMDSException("orderedNodes.size() != 2");
    }

    TSimplexID t = addTriangle(orderedNodes.back(), orderedNodes.front(), nodeToConnect, false);
    createdCells.push_back(-t);
    std::vector<TInt> v{nodeA, nodeB};
    std::vector<TInt> otherNodeInOppTriangle = SimplicesTriangle(this, oppositeTriangle).getOtherNodeInSimplex(v);
    if(otherNodeInOppTriangle.size() != 1)
    {
      throw gmds::GMDSException("otherNodeInOppTriangle.size() != 1");
    }
    TInt lNode_oppositeTriangle = SimplicesTriangle(this, oppositeTriangle).getLocalNode(otherNodeInOppTriangle.front());
    m_tri_adj[oppositeTriangle][lNode_oppositeTriangle] = t;
    m_tri_adj[t][2] = oppositeTriangle;
    (*BND_TRIANGLES)[t] = indice;

    v = std::vector<TInt>{nodeA, nodeB, nodeToConnect};
    std::vector<TInt> otherNodesIn_T = SimplicesCell(this, T).getOtherNodeInSimplex(v);
    if(otherNodesIn_T.size() != 1)
    {
      throw gmds::GMDSException("otherNodesIn_T.size() != 1");
    }
    TInt lNode_T = SimplicesCell(this, T).getLocalNode(otherNodesIn_T.front());
    m_tri_nodes[t][3] = T;
    m_tet_adj[T][lNode_T] = -t;


    std::pair<TInt, TInt> p0 = my_make(nodeA, nodeToConnect);
    std::pair<TInt, TInt> p1 = my_make(nodeB, nodeToConnect);
    hash_edge.insert(std::pair<std::pair<TInt, TInt>, TSimplexID>(p0, t));
    hash_edge.insert(std::pair<std::pair<TInt, TInt>, TSimplexID>(p1, t));

    mapForTrianglesColor.erase(it0);
  }

  //TRIANGLE RECOMBINAISONBETWEEN THE PREVIOUS CREATED TRIANGLE
  while(!hash_edge.empty())
  {
    std::multimap<std::pair<TInt, TInt>, TSimplexID>::iterator it0 = hash_edge.begin();
    std::pair<std::pair<TInt, TInt>, TSimplexID> p = *it0;
    std::pair<TInt, TInt> pairNode = my_make(p.first.first, p.first.second); //std::pair(nodeA, nodeB)
    std::pair <std::multimap<std::pair<TInt, TInt>, TSimplexID>::iterator, std::multimap<std::pair<TInt, TInt>, TSimplexID>::iterator> pair_it;
    pair_it = hash_edge.equal_range(pairNode);
    int count = std::distance(pair_it.first, pair_it.second);

    if(count > 2)
    {
      std::cout << "hash_edge filling algo" << std::endl;
      std::cout << "p.first.first -> " << p.first.first << " | " << "p.first.second -> " << p.first.second << std::endl;
      gmds::ISimplexMeshIOService ioServiceMesh(this);
      gmds::VTKWriter vtkWriterMesh(&ioServiceMesh);
      vtkWriterMesh.setCellOptions(gmds::N|gmds::R|gmds::F);
      vtkWriterMesh.setDataOptions(gmds::N|gmds::R|gmds::F);
      vtkWriterMesh.write("DEBUG_MESH.vtk");
      std::cout << "count -> " << count << std::endl;
      throw gmds::GMDSException("count > 2");
    }
    else if(count == 2)
    {
      std::multimap<std::pair<TInt, TInt>, TSimplexID>::iterator it1 = hash_edge.begin();
      ++it1;
      TSimplexID simplexA = it0->second;
      TSimplexID simplexB = it1->second;

      const SimplicesTriangle cA = SimplicesTriangle(this, simplexA);
      const SimplicesTriangle cB = SimplicesTriangle(this, simplexB);

      std::vector<TInt> intersectionEdge = cA.intersectionNodes(cB);
      if(intersectionEdge.size() != 2)
      {
        throw gmds::GMDSException("intersectionEdge.size() != 2");
      }

      std::vector<TInt> otherNode_A = cA.getOtherNodeInSimplex(intersectionEdge);
      std::vector<TInt> otherNode_B = cB.getOtherNodeInSimplex(intersectionEdge);

      if(otherNode_A.size() != 1 || otherNode_B.size() != 1)
      {
        throw gmds::GMDSException("otherNode_A != 1 || otherNode_B != 1");
      }

      TInt localNode_A = cA.getLocalNode(otherNode_A.front());
      TInt localNode_B = cB.getLocalNode(otherNode_B.front());

      m_tri_adj[simplexA][localNode_A] = simplexB;
      m_tri_adj[simplexB][localNode_B] = simplexA;

      hash_edge.erase(it0);
      hash_edge.erase(it1);
    }
    else if(count ==1)
    {
      //reinsertion node
      const TSimplexID simplexA = p.second;
      TInt nodeA = (p.first.first == nodeToConnect)? p.first.second: p.first.first ;
      std::vector<TInt> nodesEdge{p.first.first, p.first.second};

      TSimplexID simplexAdj = border;
      try{
        simplexAdj = mapForReinsertionSurface.at(nodeA);
      } catch(std::out_of_range& oor)
      {
        std::cout << "nodeToConnect -> " << nodeToConnect << std::endl;
        std::cout << "nodeA -> " << nodeA <<std::endl;
        std::cout << "simplexA ->" << simplexA << std::endl;

        std::vector<TSimplexID> simplices{simplexA};
        deleteAllSimplicesBut(simplices);

        gmds::ISimplexMeshIOService ioServiceMesh(this);
        gmds::VTKWriter vtkWriterMesh(&ioServiceMesh);
        vtkWriterMesh.setCellOptions(gmds::N|gmds::R|gmds::F);
        vtkWriterMesh.setDataOptions(gmds::N|gmds::R|gmds::F);
        vtkWriterMesh.write("DEBUG_MESH.vtk");

        std::cerr << "Out of Range error: " << oor.what() << '\n';
        throw gmds::GMDSException("mapForReinsertionSurface.at(nodeA) does not exist");
      }

      const SimplicesTriangle triangleAdj = SimplicesTriangle(this, simplexAdj);
      const SimplicesTriangle triangle    = SimplicesTriangle(this, simplexA);

      std::vector<TInt> otherNode_Adj = triangleAdj.getOtherNodeInSimplex(nodesEdge);
      std::vector<TInt> otherNode_A   = triangle.getOtherNodeInSimplex(nodesEdge);
      if(otherNode_Adj.size() != 1 || otherNode_A.size() != 1)
      {
        throw gmds::GMDSException(" ! (otherNode_Adj.size() != 1 || otherNode_A.size() != 1)");
      }

      TInt localNode_Adj    = triangleAdj.getLocalNode(otherNode_Adj.front());
      TInt localNode_cell   = triangle.getLocalNode(otherNode_A.front());

      m_tri_adj[simplexA][localNode_cell] = simplexAdj;
      m_tri_adj[simplexAdj][localNode_Adj] = simplexA;

      hash_edge.erase(it0);
    }
  }
}
/******************************************************************************/
void SimplexMesh::rebuildCavity(CavityOperator::CavityIO& cavityIO, const std::vector<std::vector<TInt>>& deleted_Tet, const std::vector<std::vector<TInt>>& deleted_Tri,  const TInt nodeToConnect, std::vector<TSimplexID>& createdCells)
{
  Variable<int>* BND_CURVE_COLOR    = getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
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
    TSimplexID neighborSimplex = extSimplexBorder[faceIter];

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
    createdCells.push_back(simplex);

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

        if(nodeA != nodeToConnect && nodeB != nodeToConnect)
        {
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
              throw GMDSException("localnode >= 3");
            }
          }
          else
          {
            std::cout << "nodeA -> " << nodeA << std::endl;
            std::cout << "nodeB -> " << nodeB << std::endl;
            std::cout << std::endl;
            std::cout << "triangle : " << triangle << std::endl;
            std::cout << SimplicesTriangle(this, triangle) << std::endl;
            std::cout << "oppositeTriangle : " << oppositeTriangle << std::endl;
            std::cout << SimplicesTriangle(this, oppositeTriangle) << std::endl;
            std::cout << std::endl;
            std::vector<TSimplexID> cellInfo{};
            std::vector<TSimplexID> trianglesInfo{std::abs(triangle), std::abs(oppositeTriangle)};
            deleteAllTrianglesBut(trianglesInfo);
            deleteAllSimplicesBut(cellInfo);

            gmds::ISimplexMeshIOService ioServiceMesh(this);
            gmds::VTKWriter vtkWriterMesh(&ioServiceMesh);
            vtkWriterMesh.setCellOptions(gmds::N|gmds::R|gmds::F);
            vtkWriterMesh.setDataOptions(gmds::N|gmds::R|gmds::F);
            vtkWriterMesh.write("DEBUG_MESH.vtk");
            std::cout << "oppoNode.size() -> " << oppoNode.size() << std::endl;
            for(auto const node : oppoNode)
            {
              std::cout << "  node -> " << node << std::endl;
            }
            throw GMDSException("oppoNode.size() != 1");
          }

          hash_edge[nodeA].push_back(-triangle);
          hash_edge[nodeB].push_back(-triangle);
        }
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
    if(pair.second.size() == 0)
    {
      throw GMDSException("pair.second.size() == 0");
    }
    else if(pair.second.size() <= 2)
    {
      if(pair.second.size() == 1)
      {
        TSimplexID idxTriangleA = pair.second.front();
        SimplicesTriangle triangleA = SimplicesTriangle(this, idxTriangleA);
        std::vector<TInt> nodesTriangleA = triangleA.getNodes();
        TInt localNode = (nodesTriangleA[0] != node)? 0 : 1;

        /*const std::vector<TSimplexID> shell = SimplicesNode(this, node).shell(SimplicesNode(this, nodeToConnect));
        for(auto const triangle : shell)
        {
          if(triangle < 0 && triangle != idxTriangleA)
          {
            m_tri_adj[idxTriangleA][localNode] = -triangle;
            std::vector<TInt> nodeShell{node, nodeToConnect};
            std::vector<TInt> otherNodes = SimplicesTriangle(this, -triangle).getOtherNodeInSimplex(nodeShell);
            if(otherNodes.size() == 1)
            {
              TInt lNode = SimplicesTriangle(this, -triangle).getLocalNode(otherNodes.front());
              m_tri_adj[-triangle][lNode] = idxTriangleA;
              break;
            }
            else
            {
              throw gmds::GMDSException("otherNodes.size() == 1 in hash EDGE");
            }
          }
        }*/

        for(auto const triangles : oppositeTriangles)
        {
          for(auto const triangle : triangles)
          {
            if(triangle != border)
            {
              std::vector<TInt> interNodes = triangleA.intersectionNodes(SimplicesTriangle(this,triangle));
              if(interNodes.size() == 2)
              {
                std::vector<TInt> otherNodesA = triangleA.otherNodesInTriangle(SimplicesTriangle(this,triangle));
                std::vector<TInt> otherNodes = SimplicesTriangle(this,triangle).otherNodesInTriangle(triangleA);

                if(otherNodesA.size() == 1 && otherNodes.size() == 1)
                {
                  const TInt lNodeA = triangleA.getLocalNode(otherNodesA.front());
                  const TInt lNode = SimplicesTriangle(this,triangle).getLocalNode(otherNodes.front());
                  m_tri_adj[idxTriangleA][lNodeA] = -triangle;
                  m_tri_adj[-triangle][lNode] = idxTriangleA;
                }
                else
                {
                  throw gmds::GMDSException("!(otherNodesA.size() == 1 && otherNodes.size() == 1)");
                }
              }
            }
          }
        }
        if(m_tri_adj[idxTriangleA][localNode] == border)
        {
          std::cout << "idxTriangleA -> " << SimplicesTriangle(this, idxTriangleA) << std::endl;
          std::cout << "localNode    -> " << localNode << std::endl;
          gmds::ISimplexMeshIOService ioServiceMesh(this);
          gmds::VTKWriter vtkWriterMesh(&ioServiceMesh);
          vtkWriterMesh.setCellOptions(gmds::N|gmds::R|gmds::F);
          vtkWriterMesh.setDataOptions(gmds::N|gmds::R|gmds::F);
          vtkWriterMesh.write("DEBUG_MESH.vtk");


          SimplexMesh simplexMeshDebug = SimplexMesh();
          std::vector<TInt> new_Nodes_Already_Added{};
          std::map<TInt, TInt> m{};
          for(auto const nodes : deleted_Tet)
          {
            std::vector<TInt> new_Nodes{};
            for(auto const node : nodes)
            {
              if(std::find(new_Nodes_Already_Added.begin(), new_Nodes_Already_Added.end(), node) == new_Nodes_Already_Added.end())
              {
                new_Nodes_Already_Added.push_back(node);
                const math::Point& pt = SimplicesNode(this, node).getCoords();
                new_Nodes.push_back(simplexMeshDebug.addNode(pt));
                std::pair<TInt, TInt> p{node, new_Nodes.back()};
                m.insert(p);
                if(node == nodeToConnect)
                {
                  std::cout << "nex node to connect in debug mesh is -> " << new_Nodes.back() << std::endl;
                }
              }
              else
              {
                new_Nodes.push_back(m[node]);
              }
            }
            simplexMeshDebug.addTetraedre(new_Nodes[0], new_Nodes[1], new_Nodes[2], new_Nodes[3]);
          }

          for(auto const nodes : deleted_Tri)
          {
            std::vector<TInt> new_Nodes{};
            for(auto const node : nodes)
            {
              const math::Point& pt = SimplicesNode(this, node).getCoords();
              new_Nodes.push_back(simplexMeshDebug.addNode(pt));
            }
            simplexMeshDebug.addTriangle(new_Nodes[0], new_Nodes[1], new_Nodes[2]);
          }
          gmds::ISimplexMeshIOService ioServiceMeshCavity(&simplexMeshDebug);
          gmds::VTKWriter vtkWriterMeshCavity(&ioServiceMeshCavity);
          vtkWriterMeshCavity.setCellOptions(gmds::N|gmds::R|gmds::F);
          vtkWriterMeshCavity.setDataOptions(gmds::N|gmds::R|gmds::F);
          vtkWriterMeshCavity.write("DEBUG_MESH_CAVITY.vtk");

          throw GMDSException("m_tri_adj[idxTriangleA][localNodeA] == border");
        }
      }
      else if(pair.second.size() == 2)
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
    }
    else
    {
      SimplexMesh simplexMeshDebug = SimplexMesh();
      std::vector<TInt> new_Nodes_Already_Added{};
      std::map<TInt, TInt> m{};
      for(auto const nodes : deleted_Tet)
      {
        std::vector<TInt> new_Nodes{};
        for(auto const node : nodes)
        {
          if(std::find(new_Nodes_Already_Added.begin(), new_Nodes_Already_Added.end(), node) == new_Nodes_Already_Added.end())
          {
            new_Nodes_Already_Added.push_back(node);
            const math::Point& pt = SimplicesNode(this, node).getCoords();
            new_Nodes.push_back(simplexMeshDebug.addNode(pt));
            std::pair<TInt, TInt> p{node, new_Nodes.back()};
            m.insert(p);
            if(node == nodeToConnect)
            {
              std::cout << "nex node to connect in debug mesh is -> " << new_Nodes.back() << std::endl;
            }
          }
          else
          {
            new_Nodes.push_back(m[node]);
          }
        }
        simplexMeshDebug.addTetraedre(new_Nodes[0], new_Nodes[1], new_Nodes[2], new_Nodes[3]);
      }

      for(auto const nodes : deleted_Tri)
      {
        std::vector<TInt> new_Nodes{};
        for(auto const node : nodes)
        {
          const math::Point& pt = SimplicesNode(this, node).getCoords();
          new_Nodes.push_back(simplexMeshDebug.addNode(pt));
        }
        simplexMeshDebug.addTriangle(new_Nodes[0], new_Nodes[1], new_Nodes[2]);
      }
      gmds::ISimplexMeshIOService ioServiceMesh(&simplexMeshDebug);
      gmds::VTKWriter vtkWriterMesh(&ioServiceMesh);
      vtkWriterMesh.setCellOptions(gmds::N|gmds::R|gmds::F);
      vtkWriterMesh.setDataOptions(gmds::N|gmds::R|gmds::F);
      vtkWriterMesh.write("DEBUG_MESH_pair.vtk");

      throw GMDSException("pair.second.size() > 2");
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
      break;
      //return;
    }

    TSimplexID simplexId   = tetIdToBuildInfoCopy.back();
    std::vector<TInt> face = pointsToConnectCopy.back();

    tetIdToBuildInfoCopy.pop_back();
    pointsToConnectCopy.pop_back();

    if(simplexId == border){continue;}
    if(face.size() == 0)
    {
      throw gmds::GMDSException("face.size() == 0");
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

  //correction of the cavity if plane tetraedron  was build by reinsertion of node for example
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
              throw GMDSException("othersNode0.size() != 1 || othersNode1.size() == 1");
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
                throw GMDSException("othersNode.size() != 1");
              }
            }
          }
          m_tet_ids.unselect(simplex);
          createdCells.erase(std::find(createdCells.begin(), createdCells.end(), simplex));
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
            throw GMDSException("othersNode0.size() != 1 || othersNode1.size() == 1");
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
bool SimplexMesh::edgeRemove(const TInt nodeA, const TInt nodeB)
{
  Variable<int>* BND_VERTEX_COLOR  = nullptr;
  Variable<int>* BND_CURVE_COLOR   = nullptr;
  Variable<int>* BND_SURFACE_COLOR = nullptr;

  try{
    BND_VERTEX_COLOR  = getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
    BND_CURVE_COLOR   = getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
    BND_SURFACE_COLOR = getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
  } catch(gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  const SimplicesNode sNodeA(this, nodeA);
  const SimplicesNode sNodeB(this, nodeB);

  unsigned int nodeADim = ((*BND_VERTEX_COLOR)[nodeA] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[nodeA] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[nodeA] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
  unsigned int nodeALabel = ((*BND_VERTEX_COLOR)[nodeA] != 0)?(*BND_VERTEX_COLOR)[nodeA]:((*BND_CURVE_COLOR)[nodeA] != 0)?(*BND_CURVE_COLOR)[nodeA]:((*BND_SURFACE_COLOR)[nodeA] != 0)?(*BND_SURFACE_COLOR)[nodeA]:0;
  unsigned int nodeBDim = ((*BND_VERTEX_COLOR)[nodeB] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[nodeB] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[nodeB] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
  unsigned int nodeBLabel = ((*BND_VERTEX_COLOR)[nodeB] != 0)?(*BND_VERTEX_COLOR)[nodeB]:((*BND_CURVE_COLOR)[nodeB] != 0)?(*BND_CURVE_COLOR)[nodeB]:((*BND_SURFACE_COLOR)[nodeB] != 0)?(*BND_SURFACE_COLOR)[nodeB]:0;

  bool status = false;
  CriterionRAIS criterionRAIS(new VolumeCriterion());
  const gmds::BitVector markedNodes{};
  std::vector<TSimplexID> deletedSimplex{};
  std::vector<TInt> deletedNodes{};
  const std::multimap<TInt, TInt>& facesAlreadyBuilt{};
  std::vector<TSimplexID> cellsCreated{};
  std::vector<TSimplexID> ball{};

  TInt nodeToReinsert = 0;

  if(nodeADim == topo::CORNER && nodeBDim == topo::CORNER)
  {
    return status;
  }
  if((nodeADim == SimplexMesh::topo::CORNER))
  {
    if((nodeBDim == SimplexMesh::topo::RIDGE))
    {
      ball = sNodeB.ballOf();
      nodeToReinsert = nodeA;
    }
  }
  if((nodeADim == SimplexMesh::topo::RIDGE))
  {
    if((nodeBDim == SimplexMesh::topo::CORNER))
    {
      ball = sNodeA.ballOf();
      nodeToReinsert = nodeB;
    }
  }
  if((nodeADim == SimplexMesh::topo::RIDGE))
  {
    if((nodeBDim == SimplexMesh::topo::VOLUME))
    {
      ball = sNodeB.ballOf();
      nodeToReinsert = nodeA;
    }
  }
  if((nodeADim == SimplexMesh::topo::VOLUME))
  {
    if((nodeBDim == SimplexMesh::topo::RIDGE))
    {
      ball = sNodeA.ballOf();
      nodeToReinsert = nodeB;
    }
  }
  if((nodeADim == SimplexMesh::topo::SURFACE))
  {
    if((nodeBDim == SimplexMesh::topo::VOLUME))
    {
      ball = sNodeB.ballOf();
      nodeToReinsert = nodeA;
    }
  }
  if((nodeADim == SimplexMesh::topo::VOLUME))
  {
    if((nodeBDim == SimplexMesh::topo::SURFACE))
    {
      ball = sNodeA.ballOf();
      nodeToReinsert = nodeB;
    }
  }//////
  if((nodeADim == SimplexMesh::topo::RIDGE))
  {
    if((nodeBDim == SimplexMesh::topo::SURFACE))
    {
      ball = sNodeB.ballOf();
      nodeToReinsert = nodeA;
    }
  }
  if((nodeADim == SimplexMesh::topo::SURFACE))
  {
    if((nodeBDim == SimplexMesh::topo::RIDGE))
    {
      ball = sNodeA.ballOf();
      nodeToReinsert = nodeB;
    }
  }
  //
  if(ball.size() == 0)
  {
    //the ball can be whatever
    ball = sNodeA.ballOf();
    nodeToReinsert = nodeB;
  }
  //

  if(ball.size() == 0){return false;}
  PointInsertion pi(this, SimplicesNode(this, nodeToReinsert), criterionRAIS, status, ball, markedNodes, deletedNodes, facesAlreadyBuilt, cellsCreated);
  return status;
}
/******************************************************************************/
unsigned int SimplexMesh::edgesRemove(const gmds::BitVector& nodeBitVector, std::vector<TInt>& deletedNodes)
{
  gmds::ISimplexMeshIOService ioService(this);
  static int cpt = 0;
  unsigned int edgesRemovedNbr = 0;
  TInt border = std::numeric_limits<TInt>::min();
  Variable<int>* BND_VERTEX_COLOR = getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  Variable<int>* BND_CURVE_COLOR = getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR = getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
  std::vector<TInt> nodesNotConnected{};
  CriterionRAIS criterionRAIS(new VolumeCriterion());
  unsigned int nodeNotInserted = 0;
  gmds::BitVector surfaceNodesAdded(m_node_ids.capacity());
  for(TInt node = 0; node < m_node_ids.capacity(); node++)
  {
    if(nodeBitVector[node] == 1)
    {
      if((*BND_VERTEX_COLOR)[node] != 0 || (*BND_CURVE_COLOR)[node] != 0 || (*BND_SURFACE_COLOR)[node] != 0)
      {
        surfaceNodesAdded.assign(node);
      }
    }
  }

  //Dimension node boundary vertex --> 0 / node on curve --> 1 / node on surface --> 2 / node on the mesh --> 3
  for(TInt node = 0; node < m_node_ids.capacity(); node++)
  {
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

        //if(dim_Ni != 0)
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
                adjacentNodesInfo.insert( std::pair<double, dataNode>(lenght, nodeInfo));
              }
            }
            //concatenation of adjacentNodesInfo & adjacentNodesInfoTMP
            std::copy(adjacentNodesInfo.begin(), adjacentNodesInfo.end(), std::back_inserter(concattedAdjNodes));
          }

          bool status = false;
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
              //std::cout << "node being inserted --> " << data.node << " Of dimension -> " <<  data.dim_Nj << " and label -> " << data.index_Nj << std::endl;
              //std::cout << "from node --> " << node << " Of dimension -> " <<  dim_Ni << " and label -> " << index_Ni << std::endl;

              //if(dim_Ni == 4 && (data.dim_Nj == 0 || data.dim_Nj == 1 || data.dim_Nj == 2)){continue;}
              const std::multimap<TInt, TInt> facesAlreadyBuilt{};
              std::vector<TSimplexID> cellsCreated{};
              PointInsertion(this, nodeToInsert, criterionRAIS, status, ball, surfaceNodesAdded, deletedNodes, facesAlreadyBuilt, cellsCreated);
              if(status)
              {
                edgesRemovedNbr++;
                break;
              }
            }
          }
          if(!status)
          {
            nodeNotInserted++;
          }
        }
      }
    }
  }

  return edgesRemovedNbr;
}
/******************************************************************************/
std::vector<TSimplexID> SimplexMesh::initializeCavityWith(const TInt nodeAidx, const TInt nodeBidx)
{
  std::vector<TSimplexID> cavity{};
  const TInt border = std::numeric_limits<TInt>::min();
  SimplicesNode nodeA = SimplicesNode(this, nodeAidx);
  SimplicesNode nodeB = SimplicesNode(this, nodeBidx);
  std::vector<TSimplexID> ballA = nodeA.ballOf();
  std::vector<TSimplexID> ballB = nodeB.ballOf();

  if(ballA.size() == 0 || ballB.size() == 0)
  {
    return cavity;
  }

  std::vector<double> distance{};
  for(auto const & simplex : ballA)
  {
    std::vector<TInt> nodes{};
    std::vector<TInt> currentNode{nodeA.getGlobalNode()};
    if(simplex >= 0){nodes = SimplicesCell(this, simplex).getOtherNodeInSimplex(currentNode);}
    else{ distance.push_back(border) ; continue;}
    math::Point p0 = SimplicesNode(this, nodes[0]).getCoords();
    math::Point p1 = SimplicesNode(this, nodes[1]).getCoords();
    math::Point p2 = SimplicesNode(this, nodes[2]).getCoords();
    math::Vector3d vec = nodeB.getCoords() - nodeA.getCoords();
    math::Plane P(p0, p1, p2);
    math::Ray r(nodeA.getCoords(), vec);
    math::Point p;

    //math::Point barFace = 1.0 / 3.0 * (p0 + p1 + p2);
    //math::Vector3d vecTest = barFace - nodeA.getCoords();
    if(r.intersect3D(P,p) /*&& vecTest.dot(vec) >= 0.0*/)
    {
      math::Vector3d v = p - nodeA.getCoords();
      distance.push_back(v.norm());
    }
    else
    {
      distance.push_back(border);
    }
  }

  TSimplexID idSimplex = 0;
  double distanceMax = std::numeric_limits<double>::max();
  for(unsigned int id = 0 ; id < distance.size() ; id++)
  {
    double d = distance[id];
    if(d >= 0.0 && d < distanceMax)
    {
      idSimplex = ballA[id];
      distanceMax = d;
    }
  }
  cavity.push_back(idSimplex);
  return cavity;

}
/******************************************************************************/
unsigned int SimplexMesh::buildEdges(const std::multimap<TInt, TInt>& AEdges, const gmds::BitVector& nodeBitVector)
{
  const TInt border = std::numeric_limits<TInt>::min();
  const gmds::BitVector& bitNodes = getBitVectorNodes();
  unsigned int cpt_EdgeAlreadyBuild = 0;
  double cpt_EdgeBuilt = 0.0;
  double sizeEdge = AEdges.size();
  std::vector<TInt> allDeletedNodes{};
  gmds::ISimplexMeshIOService ioService(this);

  unsigned int cpt = 0;
  for(auto const & edge : AEdges)
  {
    TInt nodeAidx = edge.first;
    TInt nodeBidx = edge.second;

    /*std::cout << "nodeAidx -> " << nodeAidx << std::endl;
    std::cout << "nodeBidx -> " << nodeBidx << std::endl;
    if(nodeAidx == 23517 && nodeBidx == 72468)
    {
      gmds::VTKWriter vtkWriter(&ioService);
      vtkWriter.setCellOptions(gmds::N|gmds::R|gmds::F);
      vtkWriter.setDataOptions(gmds::N|gmds::R|gmds::F);
      vtkWriter.write("TEST.vtk");
    }*/
    if(nodeAidx != border && nodeBidx != border)
    {
      if(bitNodes[nodeAidx] != 0 && bitNodes[nodeBidx] != 0)
      {
        SimplicesNode nodeA = SimplicesNode(this, nodeAidx);
        SimplicesNode nodeB = SimplicesNode(this, nodeBidx);

        if(nodeA.shell(nodeB).size() == 0)
        {
          std::vector<TInt> cavity = initializeCavityWith(nodeA.getGlobalNode(), nodeB.getGlobalNode());
          //for(auto const & tet : cavity){std::cout << "tet in cav -> " << tet << std::endl;}
          CriterionRAIS criterionRAIS(new VolumeCriterion());
          bool status = false;
          std::vector<TInt> deletedNodes{};
          const std::multimap<TInt, TInt> facesAlreadyBuilt{};
          std::vector<TSimplexID> createdCells{};
          //std::cout << SimplicesNode(this, edge.first) << std::endl;
          //std::cout << SimplicesNode(this, edge.second) << std::endl;
          //if(cpt > 500)
          {

            //gmds::VTKWriter vtkWriterTEST(&ioService);
            //vtkWriterTEST.setCellOptions(gmds::N|gmds::R|gmds::F);
            //vtkWriterTEST.setDataOptions(gmds::N|gmds::R|gmds::F);
            //vtkWriterTEST.write("testEdge4_CPT_" + std::to_string(cpt) +  ".vtk");
            //std::cout << "CPT -> " << cpt << std::endl;
          }
          PointInsertion(this, nodeB, criterionRAIS, status, cavity, nodeBitVector, deletedNodes, facesAlreadyBuilt, createdCells);
          //cpt++;
          if(!status)
          {
            //std::cout << "edge [" << nodeAidx << " ; " << nodeBidx << "] was not built " << std::endl;
          }
          else
          {
            cpt_EdgeBuilt = cpt_EdgeBuilt + 1.0;
          }
        }
        else
        {
          cpt_EdgeAlreadyBuild = cpt_EdgeAlreadyBuild + 1.0;
        }
      }
    }
  }

  return cpt_EdgeBuilt;
}
/******************************************************************************/
bool SimplexMesh::isHexaEdgeBuild(const std::vector<std::vector<TInt>>& ANodesFaces)
{
  for(auto const & face : ANodesFaces)
  {
    //std::cout <<  "   face -> " << face[0] << " | " << face[1] << " | " << face[2] << " | " << face[3] << std::endl;
    if(!isFaceBuild(face))
    {
      return false;
    }
  }
  return true;
}
/******************************************************************************/
void SimplexMesh::whatFaceIsBuilt(const std::vector<std::vector<TInt>>& ANodesFaces, std::multimap<TInt, TInt>& facesAlreadyBuilt)
{
  std::vector<bool> res{};
  if(ANodesFaces.size() == 6)
  {
    for(unsigned int idx = 0 ; idx < ANodesFaces.size() ; idx++)
    {
      std::vector<TInt> face = ANodesFaces[idx];
      if(isFaceBuild(face))
      {
        TInt node0 = face[0] ; TInt node1 = face[1] ; TInt node2 = face[2] ; TInt node3 = face[3] ;
        std::vector<TInt> subFace0{node0, node1, node2};
        std::vector<TInt> subFace1{node2, node3, node0};
        std::sort(subFace0.begin(), subFace0.end());
        std::sort(subFace1.begin(), subFace1.end());

        std::pair<TInt, TInt> p0{subFace0[1], subFace0[2]};
        std::pair<TInt, TInt> p1{subFace1[1], subFace1[2]};
        //facesAlreadyBuilt.insert(std::make_pair(subFace0[0], p0));
        //facesAlreadyBuilt.insert(std::make_pair(subFace1[0], p1));
      }
    }
  }
  else
  {
    std::cout << "ANodesFaces.size() != 6" << std::endl;
    return;
  }
}
/******************************************************************************/
static int faceAdded = 0;
bool SimplexMesh::buildFace(const std::vector<TInt>& nodes, const gmds::BitVector& nodeAdded, const std::multimap<TInt, TInt>& facesAlreadyBuilt)
{
  TInt n0 = nodes[0] ; TInt n1 = nodes[1] ; TInt n2 = nodes[2] ; TInt n3 = nodes[3] ;
  TInt n4 = nodes[4] ; TInt n5 = nodes[5] ; TInt n6 = nodes[6] ; TInt n7 = nodes[7] ;

  std::vector<TInt> face0{n0, n1, n2, n3}; std::vector<TInt> face3{n2, n6, n7, n3};
  std::vector<TInt> face1{n4, n5, n6, n7}; std::vector<TInt> face4{n1, n5, n6, n2};
  std::vector<TInt> face2{n0, n1, n5, n4}; std::vector<TInt> face5{n0, n4, n7, n3};
  std::vector<std::vector<TInt>> faces{face0, face1, face2, face3, face4, face5};

  Variable<int>* BND_VERTEX_COLOR = getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  Variable<int>* BND_CURVE_COLOR = getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR = getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  for(unsigned int idx = 0 ; idx < faces.size() ; idx++)
  {
    std::vector<TInt> face = faces[idx];
    const SimplicesNode node0 = SimplicesNode(this, face[0]);
    const SimplicesNode node1 = SimplicesNode(this, face[1]);
    const SimplicesNode node2 = SimplicesNode(this, face[2]);
    const SimplicesNode node3 = SimplicesNode(this, face[3]);

    std::vector<TSimplexID> e01 = node0.shell(node1);
    std::vector<TSimplexID> e12 = node1.shell(node2);
    std::vector<TSimplexID> e23 = node2.shell(node3);
    std::vector<TSimplexID> e30 = node3.shell(node0);

    if(e01.size() != 0 && e12.size() != 0 && e23.size() != 0 && e30.size() != 0)
    {
      //check for the diagonale
      std::vector<TSimplexID> e02 = node0.shell(node2);
      std::vector<TSimplexID> e13 = node1.shell(node3);


      if(e02.size() == 0 && e13.size() == 0)
      {
        for(unsigned int idx = 0 ; idx < face.size() - 2 ; idx++)
        {
          TInt nodeA = face[idx];
          TInt nodeB = face[idx + 2];
          std::vector<TSimplexID> cavity = initializeCavityWith(nodeA, nodeB);
          CriterionRAIS criterionRAIS(new VolumeCriterion());
          std::vector<TSimplexID> deletedNodes{};
          std::vector<TSimplexID> createdCells{};
          bool status = false;
          PointInsertion(this, SimplicesNode(this, nodeB), criterionRAIS, status, cavity, nodeAdded, deletedNodes, facesAlreadyBuilt, createdCells);
          if(!status)
          {
            std::vector<TSimplexID> cavity = initializeCavityWith(nodeB, nodeA);
            PointInsertion(this, SimplicesNode(this, nodeA), criterionRAIS, status, cavity, nodeAdded, deletedNodes, facesAlreadyBuilt, createdCells);
          }
        }
      }

      e02 = node0.shell(node2);
      e13 = node1.shell(node3);

      if(e02.size() != 0 && e13.size() != 0)
      {
        //flat tetrahedron exist
        continue;
      }
      else
      {
        for(unsigned int node = 0 ; node < 2 ; node++)
        {
          TInt nodeAidx = face[node];
          TInt nodeBidx = face[(node + 2) % 4];
          const SimplicesNode nodeA = SimplicesNode(this, nodeAidx);
          const SimplicesNode nodeB = SimplicesNode(this, nodeBidx);
          const std::vector<TSimplexID> shellAB = nodeA.shell(nodeB);

          if(nodeA.shell(nodeB).size())
          {
            TInt nodeCidx = face[(node + 1) % 4];
            TInt nodeDidx = face[(node + 3) % 4];

            std::vector<TInt> otherNodes{nodeCidx, nodeDidx};
            std::vector<TInt> localNodeConnected{0,1};
            for(unsigned int idx = 0 ; idx < otherNodes.size() ; idx++)
            {
              const SimplicesNode currentNode = SimplicesNode(this, otherNodes[idx]);
              for(auto const simplex : shellAB)
              {
                if(simplex >= 0)
                {
                  if(SimplicesCell(this, simplex).containNode(currentNode))
                  {
                    localNodeConnected.erase(std::remove(localNodeConnected.begin(), localNodeConnected.end(), idx), localNodeConnected.end());
                  }
                }
              }
            }


            if(localNodeConnected.size() == 0)
            {
              //face already built
              break;
            }
            else if(localNodeConnected.size() == 1)
            {
              TInt localNodeFace = localNodeConnected.front();
              const SimplicesNode nodeToInsert = SimplicesNode(this, otherNodes[localNodeFace]);
              std::vector<TSimplexID> cavity = shellAB;
              CriterionRAIS criterionRAIS(new VolumeCriterion());
              std::vector<TSimplexID> deletedNodes{};
              std::vector<TSimplexID> createdCells{};
              bool status = false;
              PointInsertion(this, nodeToInsert, criterionRAIS, status, cavity, nodeAdded, deletedNodes, facesAlreadyBuilt, createdCells);
              if(status == true)
              {
                break;
              }
            }
            else
            {
              //none of the simplex in shellAB is connected to their opposite node
              //in the current face being built
              std::vector<TSimplexID> cavity = initializeCavityWith(nodeAidx, nodeBidx);
              CriterionRAIS criterionRAIS(new VolumeCriterion());
              std::vector<TSimplexID> deletedNodes{};
              std::vector<TSimplexID> createdCells{};
              bool status = false;
              PointInsertion(this, SimplicesNode(this, nodeBidx), criterionRAIS, status, cavity, nodeAdded, deletedNodes, facesAlreadyBuilt, createdCells);
              if(!status)
              {
                std::vector<TSimplexID> cavity = initializeCavityWith(nodeBidx, nodeAidx);
                PointInsertion(this, SimplicesNode(this, nodeAidx), criterionRAIS, status, cavity, nodeAdded, deletedNodes, facesAlreadyBuilt, createdCells);
                if(status == true)
                {
                  break;
                }
              }
            }
          }
        }
      }
    }
    else
    {
      //std::cout << "One of the edge was not built : e01 | e12 | e23 | e30 => " << node0.getGlobalNode() << "-" << node1.getGlobalNode() << " | "<< node1.getGlobalNode() << "-" << node2.getGlobalNode() << " | "<< node2.getGlobalNode() << "-" << node3.getGlobalNode() << " | "<< node3.getGlobalNode() << "-" << node0.getGlobalNode() << std::endl;
      return false;
    }
  }
  return isHexaEdgeBuild(faces);
}
/******************************************************************************/
/*static int faceAdded = 0;
bool SimplexMesh::buildFace(const std::vector<TInt>& nodes, const gmds::BitVector& nodeAdded)
{
  bool status = false;
  TInt n0 = nodes[0] ; TInt n1 = nodes[1] ; TInt n2 = nodes[2] ; TInt n3 = nodes[3] ;
  TInt n4 = nodes[4] ; TInt n5 = nodes[5] ; TInt n6 = nodes[6] ; TInt n7 = nodes[7] ;

  std::vector<TInt> face0{n0, n1, n2, n3}; std::vector<TInt> face3{n2, n6, n7, n3};
  std::vector<TInt> face1{n4, n5, n6, n7}; std::vector<TInt> face4{n1, n5, n6, n2};
  std::vector<TInt> face2{n0, n1, n5, n4}; std::vector<TInt> face5{n0, n4, n7, n3};
  std::vector<std::vector<TInt>> faces{face0, face1, face2, face3, face4, face5};

  Variable<int>* BND_VERTEX_COLOR = getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  Variable<int>* BND_CURVE_COLOR = getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR = getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  for(unsigned int idx = 0 ; idx < faces.size() ; idx++)
  {
    std::cout << "face id local -> " << idx << std::endl;
    std::vector<TInt> face = faces[idx];
    const SimplicesNode node0 = SimplicesNode(this, face[0]);
    const SimplicesNode node1 = SimplicesNode(this, face[1]);
    const SimplicesNode node2 = SimplicesNode(this, face[2]);
    const SimplicesNode node3 = SimplicesNode(this, face[3]);

    std::cout << node0.getGlobalNode() << std::endl;
    std::cout << node1.getGlobalNode() << std::endl;
    std::cout << node2.getGlobalNode() << std::endl;
    std::cout << node3.getGlobalNode() << std::endl;
    std::cout << face[0] << " ----> " << (*BND_SURFACE_COLOR)[face[0]] << " || " << (*BND_CURVE_COLOR)[face[0]] << std::endl;
    std::cout << face[1] << " ----> " << (*BND_SURFACE_COLOR)[face[1]] << " || " << (*BND_CURVE_COLOR)[face[1]] << std::endl;
    std::cout << face[2] << " ----> " << (*BND_SURFACE_COLOR)[face[2]] << " || " << (*BND_CURVE_COLOR)[face[2]] << std::endl;
    std::cout << face[3] << " ----> " << (*BND_SURFACE_COLOR)[face[3]] << " || " << (*BND_CURVE_COLOR)[face[3]] << std::endl;

    std::vector<TSimplexID> e01 = node0.shell(node1);
    std::vector<TSimplexID> e12 = node1.shell(node2);
    std::vector<TSimplexID> e23 = node2.shell(node3);
    std::vector<TSimplexID> e30 = node3.shell(node0);

    std::cout << "size Edge -> " << e01.size() << " : " << e12.size() << " : " <<  e23.size() << " : " <<  e30.size() << std::endl;
    if(e01.size() != 0 && e12.size() != 0 && e23.size() != 0 && e30.size() != 0)
    {
      //check for the diagonale
      std::vector<TSimplexID> e02 = node0.shell(node2);
      std::vector<TSimplexID> e13 = node1.shell(node3);

      std::cout << " nodes --> " << node0.getGlobalNode()  << " : " << node2.getGlobalNode() << " || e02.size() --> " << e02.size() << std::endl;
      std::cout << " nodes --> " << node1.getGlobalNode()  << " : " << node3.getGlobalNode() << " || e13.size() --> " << e13.size() << std::endl;

      if(e02.size() == 0 && e13.size() == 0)
      {
        std::cout << "diag not built" << std::endl;
        std::vector<SimplicesNode> nodes{node0, node2, node1, node3};
        for(unsigned int node = 0 ; node < 3 ; node = node + 2)
        {
          TInt nodeA = nodes[node].getGlobalNode();
          TInt nodeB = nodes[node + 1].getGlobalNode();
          std::vector<TSimplexID> cavity = initializeCavityWith(nodeA, nodeB);
          std::cout << "cavity --> " << cavity.front() << std::endl;
          CriterionRAIS criterionRAIS(new VolumeCriterion());
          std::vector<TSimplexID> deletedSimplex{};
          std::cout << "PointInsertion" << std::endl;
          status = false;
          PointInsertion(this, SimplicesNode(this, nodeB), criterionRAIS, status, cavity, nodeAdded, deletedSimplex);
          if(status == true)
          {
            std::cout << "faceAdded --> "<< faceAdded++ << std::endl;
            break;
          }
          else
          {
            std::vector<TSimplexID> cavity = initializeCavityWith(nodeB, nodeA);
            PointInsertion(this, SimplicesNode(this, nodeA), criterionRAIS, status, cavity, nodeAdded, deletedSimplex);
            if(status == true)
            {
              std::cout << "faceAdded --> "<< faceAdded++ << std::endl;
              break;
            }
          }
        }
        if(status == false)
        {
          std::cout << "diag can not be built" << std::endl;
          return false;
        }
      }
    }
    else
    {
      std::cout << "One of the edge is was not built" << std::endl;
      return false;
    }
  }
  std::cout << "faceAdded nbr --> " << faceAdded << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  return true;
}*/
/******************************************************************************/
std::vector<TSimplexID> SimplexMesh::hex2tet(const std::vector<TInt>& ANodesHex)
{
  //use is hexa build before..
  //intialisation of the tet..
  std::vector<TSimplexID> res{};
  if(ANodesHex.size() == 8)
  {
    TInt n0 = ANodesHex[0] ; TInt n1 = ANodesHex[1] ; TInt n2 = ANodesHex[2] ; TInt n3 = ANodesHex[3] ;
    TInt n4 = ANodesHex[4] ; TInt n5 = ANodesHex[5] ; TInt n6 = ANodesHex[6] ; TInt n7 = ANodesHex[7] ;

    std::vector<TInt> face0{n0, n3, n2, n1}; std::vector<TInt> face3{n2, n3, n7, n6};
    std::vector<TInt> face1{n4, n5, n6, n7}; std::vector<TInt> face4{n1, n2, n6, n5};
    std::vector<TInt> face2{n0, n1, n5, n4}; std::vector<TInt> face5{n0, n4, n7, n3};

    /*std::vector<TInt> face0{n0, n1, n2, n3}; std::vector<TInt> face3{n2, n6, n7, n3};
    std::vector<TInt> face1{n7, n6, n5, n4}; std::vector<TInt> face4{n1, n5, n6, n2};
    std::vector<TInt> face2{n4, n5, n1, n0}; std::vector<TInt> face5{n0, n3, n7, n4};*/

    std::vector<std::vector<TInt>> faces{face0, face1, face2, face3, face4, face5};
    if(!isHexaEdgeBuild(faces)) //check if hex faces are built
    {
      return res;
    }

    struct faceDataInfo{
      std::vector<math::Vector3d> normals;
      math::Point barycentreFace;
    };
    std::vector<faceDataInfo> facesDataInfo;
    for(auto const & face : faces)
    {
      faceDataInfo faceInfo;
      math::Point p0 = SimplicesNode(this, face[0]).getCoords();
      math::Point p1 = SimplicesNode(this, face[1]).getCoords();
      math::Point p2 = SimplicesNode(this, face[2]).getCoords();
      math::Point p3 = SimplicesNode(this, face[3]).getCoords();

      if(SimplicesNode(this, face[0]).shell(SimplicesNode(this, face[2])).size() > 0)
      {
        faceInfo.barycentreFace = 0.5 * (p0 + p2);
      }
      else if(SimplicesNode(this, face[1]).shell(SimplicesNode(this, face[3])).size() > 0)
      {
        faceInfo.barycentreFace = 0.5 * (p1 + p3);
      }
      else
      {
        std::cout << "faceInfo.barycentreFace shell size null" << std::endl;
        return res;
      }


      math::Vector3d vec0B = p0 - faceInfo.barycentreFace;
      math::Vector3d vec10 = p1 - p0;
      math::Vector3d vec0  = vec0B.cross(vec10); vec0.normalize();

      math::Vector3d vec1B = p1 - faceInfo.barycentreFace;
      math::Vector3d vec21 = p2 - p1;
      math::Vector3d vec1  = vec1B.cross(vec21); vec1.normalize();

      math::Vector3d vec2B = p2 - faceInfo.barycentreFace;
      math::Vector3d vec32 = p3 - p2;
      math::Vector3d vec2  = vec2B.cross(vec32); vec2.normalize();

      math::Vector3d vec3B = p3 - faceInfo.barycentreFace;
      math::Vector3d vec03 = p0 - p3;
      math::Vector3d vec3  = vec3B.cross(vec03); vec3.normalize();

      faceInfo.normals.push_back(vec0); faceInfo.normals.push_back(vec1);
      faceInfo.normals.push_back(vec2); faceInfo.normals.push_back(vec3);
      facesDataInfo.push_back(faceInfo);
    }

    TInt border = std::numeric_limits<TInt>::min();
    TSimplexID initTet = border;
    std::vector<TSimplexID> && ball0Tmp = SimplicesNode(this, n0).ballOf();
    std::vector<TSimplexID> && ball1Tmp = SimplicesNode(this, n1).ballOf();
    std::vector<TSimplexID> && ball2Tmp = SimplicesNode(this, n2).ballOf();
    std::vector<TSimplexID> && ball3Tmp = SimplicesNode(this, n3).ballOf();

    std::vector<TSimplexID> ball0; std::vector<TSimplexID> ball2;
    std::vector<TSimplexID> ball1; std::vector<TSimplexID> ball3;

    std::copy_if(ball0Tmp.begin(), ball0Tmp.end(), std::back_inserter(ball0), [&](TSimplexID simplex){return (simplex >= 0);});
    std::copy_if(ball1Tmp.begin(), ball1Tmp.end(), std::back_inserter(ball1), [&](TSimplexID simplex){return (simplex >= 0);});
    std::copy_if(ball2Tmp.begin(), ball2Tmp.end(), std::back_inserter(ball2), [&](TSimplexID simplex){return (simplex >= 0);});
    std::copy_if(ball3Tmp.begin(), ball3Tmp.end(), std::back_inserter(ball3), [&](TSimplexID simplex){return (simplex >= 0);});

    std::vector<std::vector<TSimplexID>> balls{{ball0}, {ball1}, {ball2}, {ball3}};
    std::vector<TSimplexID>&& inter  = intersectionSimplex(balls);
    if(inter.size()){initTet = inter.front();}

    if(initTet == border)
    {
      std::vector<std::vector<TSimplexID>> balls012{{ball0}, {ball1}, {ball2}};
      std::vector<std::vector<TSimplexID>> balls013{{ball0}, {ball1}, {ball3}};
      std::vector<TSimplexID>&& inter012  = intersectionSimplex(balls012);
      std::vector<TSimplexID>&& inter013  = intersectionSimplex(balls013);
      std::vector<TInt> nodes{};
      if(inter012.size()){inter = inter012; nodes.push_back(n0) ; nodes.push_back(n1) ; nodes.push_back(n2);}
      else if(inter013.size()){inter = inter013; nodes.push_back(n0) ; nodes.push_back(n1) ; nodes.push_back(n3);}
      else{return res;}

      for(auto const simplex : inter)
      {
        if(simplex >= 0)
        {
          std::vector<TInt> otherNodes = SimplicesCell(this, simplex).getOtherNodeInSimplex(nodes);
          const TInt currentNode       = otherNodes.front();
          if(currentNode == n0 || currentNode == n1 || currentNode == n2 || currentNode == n3 ||
            currentNode == n4 || currentNode == n5 || currentNode == n6 || currentNode == n7)
          {
            initTet = simplex;
            break;
          }
          else
          {
            const math::Point pt = SimplicesNode(this,currentNode).getCoords();
            bool insideCell = true;
            for(auto const & faceDataInfo : facesDataInfo)
            {
              math::Vector3d vec    = pt - faceDataInfo.barycentreFace;
              for(auto const & normalSubFace : faceDataInfo.normals)
              {
                math::Vector3d normal = normalSubFace;
                if(normal.dot(vec) > 0.0)
                {
                  insideCell = false;
                  break;
                }
              }
            }
            if(insideCell)
            {
              initTet = simplex;
              break;
            }
          }
        }
      }
    }


    if(initTet == border)
    {
      std::cout << "initTet == border " << std::endl;
      return res;
    }

    std::vector<TSimplexID> to_do{initTet};
    gmds::BitVector markedTet(getBitVectorTet().capacity());


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
              if(!(node == n0 || node == n1 || node == n2 || node == n3 ||
                node == n4 || node == n5 || node == n6 || node == n7))
                {
                  math::Point pt = SimplicesNode(this, node).getCoords();
                  bool insideCell = true;
                  unsigned int faceNbr = 0;
                  for(auto const & faceDataInfo : facesDataInfo)
                  {
                    math::Vector3d vec = pt - faceDataInfo.barycentreFace;
                    for(auto const & normalSubFace : faceDataInfo.normals)
                    {
                      math::Vector3d normal = normalSubFace;

                      if(normal.dot(vec) > 0.0)
                      {
                        insideCell = false;
                        break;
                      }
                    }
                    if(insideCell == false)
                    {
                      break;
                    }
                    faceNbr++;
                  }
                  if(insideCell)
                  {
                    cpt++;
                  }
                  else
                  {
                    break;
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
double SimplexMesh::subSurfaceFactor(const std::vector<std::vector<TInt>>& faces)
{
  double res = std::numeric_limits<double>::max();
  for(auto const & face : faces)
  {
    if(face.size() == 3)
    {
      //compute face's normal
      math::Vector3d n0 = SimplicesNode(this, face[0]).getNormal();
      math::Vector3d n1 = SimplicesNode(this, face[1]).getNormal();
      math::Vector3d n2 = SimplicesNode(this, face[2]).getNormal();
      n0.normalize() ; n1.normalize() ; n2.normalize() ;

      math::Vector3d normalFace = (n0 + n1).normalize() + n2 / 3.0;
      normalFace.normalize();
      double epsilon = std::min(n0.dot(normalFace), std::min(n1.dot(normalFace), n2.dot(normalFace)));
      if(epsilon < res)
      {
        res = epsilon;
      }
    }
    else
    {
      std::cout << "face.size() != 3" << std::endl;
      return res;
    }
  }

  return res;
}
/******************************************************************************/
double SimplexMesh::computeQualityEdge(const gmds::math::Point& pA, const gmds::math::Point& pB,
                                       const Eigen::Matrix3d& MA, const Eigen::Matrix3d& MB) const

{
  //https://pages.saclay.inria.fr/frederic.alauzet/download/George_Mesh%20generation%20and%20mesh%20adaptivity%20theory%20and%20techniques.pdf
  //compute the qualité of simplex with each simplex's node metric and the quality of the simplex will be the worst quality previously computed
  math::Vector3d vec = pB - pA;
  Eigen::Vector3d edge = Eigen::Vector3d(vec.X(), vec.Y(), vec.Z());
  double edgeInMetric = 0.5 * sqrt(edge.dot(MA * edge)) + 0.5 * sqrt(edge.dot(MB * edge));

  return edgeInMetric;
}
/******************************************************************************/
double SimplexMesh::computeQualityEdge(const TInt nodeA, const TInt nodeB)

{
  Variable<Eigen::Matrix3d>* metric = nullptr;
  try{
    metric = getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }
  //compute the mean of the metric for the current simplex in order to approximate a quality factor using mean metric
  const gmds::math::Point pA = SimplicesNode(this, nodeA).getCoords();
  const gmds::math::Point pB = SimplicesNode(this, nodeB).getCoords();
  Eigen::Matrix3d MA = (*metric)[nodeA];
  Eigen::Matrix3d MB = (*metric)[nodeB];

  return computeQualityEdge(pA, pB, MA, MB);
}
/******************************************************************************/
double SimplexMesh::computeQualityElement(const gmds::math::Point& pA, const gmds::math::Point& pB, const gmds::math::Point& pC, const gmds::math::Point& pD,
                                          const Eigen::Matrix3d& MA, const Eigen::Matrix3d& MB, const Eigen::Matrix3d& MC, const Eigen::Matrix3d& MD) const

{
  double alpha = 36.0 * sqrt(2.0);
  gmds::math::Point points[4] = {pA, pB, pC, pD};
  gmds::math::Tetrahedron tet(points);
  const double volumeK  = tet.getVolume();

  //https://pages.saclay.inria.fr/frederic.alauzet/download/George_Mesh%20generation%20and%20mesh%20adaptivity%20theory%20and%20techniques.pdf
  //compute the qualité of simplex with each simplex's node metric and the quality of the simplex will be the worst quality previously computed
  std::vector<double> qualities{};
  std::vector<Eigen::Matrix3d> metrics{MA, MB, MC, MD};
  math::Vector3d vec01 = pB - pA;
  math::Vector3d vec02 = pC - pA;
  math::Vector3d vec12 = pC - pB;
  math::Vector3d vec03 = pD - pA;
  math::Vector3d vec13 = pD - pB;
  math::Vector3d vec23 = pD - pC;

  Eigen::Vector3d e01 = Eigen::Vector3d(vec01.X(), vec01.Y(), vec01.Z());
  Eigen::Vector3d e02 = Eigen::Vector3d(vec02.X(), vec02.Y(), vec02.Z());
  Eigen::Vector3d e12 = Eigen::Vector3d(vec12.X(), vec12.Y(), vec12.Z());
  Eigen::Vector3d e03 = Eigen::Vector3d(vec03.X(), vec03.Y(), vec03.Z());
  Eigen::Vector3d e13 = Eigen::Vector3d(vec13.X(), vec13.Y(), vec13.Z());
  Eigen::Vector3d e23 = Eigen::Vector3d(vec23.X(), vec23.Y(), vec23.Z());

  std::vector<Eigen::Vector3d> edges{e01,e02,e12,e03,e13,e23};
  for(auto const metric : metrics)
  {
    double den = 0.0;
    for(auto const edge : edges)
    {
      //den += (edge.dot(metric * edge) * edge.dot(metric * edge));
      den += edge.dot(metric * edge);
    }
    if(den == 0.0)
    {
      throw gmds::GMDSException("den == 0.0");
    }
    else
    {
      qualities.push_back((alpha * sqrt(metric.determinant()) * volumeK) / den);
    }
  }
  return std::min(std::min(qualities[0], qualities[1]), std::min(qualities[2], qualities[3]));
}
/******************************************************************************/
double SimplexMesh::computeQualityElement(const TInt nodeA, const TInt nodeB, const TInt nodeC, const TInt nodeD)
{
  Variable<Eigen::Matrix3d>* metric = nullptr;
  try{
    metric = getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }
  //compute the mean of the metric for the current simplex in order to approximate a quality factor using mean metric
  const gmds::math::Point pA = SimplicesNode(this, nodeA).getCoords();
  const gmds::math::Point pB = SimplicesNode(this, nodeB).getCoords();
  const gmds::math::Point pC = SimplicesNode(this, nodeC).getCoords();
  const gmds::math::Point pD = SimplicesNode(this, nodeD).getCoords();


  Eigen::Matrix3d MA = (*metric)[nodeA];
  Eigen::Matrix3d MB = (*metric)[nodeB];
  Eigen::Matrix3d MC = (*metric)[nodeC];
  Eigen::Matrix3d MD = (*metric)[nodeD];

  return computeQualityElement(pA, pB, pC, pD, MA, MB, MC, MD);
}
/******************************************************************************/
double SimplexMesh::computeQualityElement(const TSimplexID simplex)
{
  if(simplex >= 0)
  {
    if(m_tet_ids[simplex] != 0)
    {
      std::vector<TInt> nodes = SimplicesCell(this, simplex).getNodes();
      return computeQualityElement(nodes[0], nodes[1], nodes[2], nodes[3]);
    }
    else
    {
      std::cout << "simplex -> " << simplex << std::endl;
      std::cout << SimplicesCell(this, simplex) << std::endl;
      throw gmds::GMDSException("m_tet_ids[simplex] != 0");
    }
  }
  else
  {
    throw gmds::GMDSException("simplex < 0");
  }
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
    //std::cout << "ABG -> " << alpha << " | " << beta << " | " << gamma << std::endl;
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
