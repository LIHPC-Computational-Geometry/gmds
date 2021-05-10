/*----------------------------------------------------------------------------*/
#include <gmds/hybridMeshAdapt/SimplicesCell.h>
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace hybrid;
using namespace simplicesNode;
using namespace simplicesCell;
/*----------------------------------------------------------------------------*/
SimplicesCell::SimplicesCell(SimplexMesh* simplexMesh, const TSimplexID simplexId)
{
    m_simplex_mesh = simplexMesh;
    m_simplexId    = simplexId;

    if(m_simplex_mesh->m_tet_ids[m_simplexId] == 0)
    {
      /*TODO exeption le node d'existe pas ... le creer avant de poouvoir l'utiliser*/
      std::cout << "Creer la cellule " << m_simplexId <<  " avant de l'utiliser !!" << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
SimplicesCell::SimplicesCell(const SimplicesCell& simplicesCell)
{
  //TODO
}
/*----------------------------------------------------------------------------*/
SimplicesCell::SimplicesCell(SimplicesCell&& simplicesCell)
{
  //TODO
}
/*----------------------------------------------------------------------------*/
SimplicesCell::~SimplicesCell()
{
}
/*----------------------------------------------------------------------------*/
std::vector<TSimplexID> SimplicesCell::getNodes() const
{
  std::vector<TSimplexID> v{};
  if(m_simplex_mesh != nullptr)
  {
    v.resize(4);
    v[0] = m_simplex_mesh->m_tet_nodes[m_simplexId][0];
    v[1] = m_simplex_mesh->m_tet_nodes[m_simplexId][1];
    v[2] = m_simplex_mesh->m_tet_nodes[m_simplexId][2];
    v[3] = m_simplex_mesh->m_tet_nodes[m_simplexId][3];
  }

  return std::move(v);
}
/*----------------------------------------------------------------------------*/
double SimplicesCell::getVolumeOfCell()
{
  if(m_simplexId > m_simplex_mesh->m_tet_nodes.size() || m_simplexId < 0)
  {
    return 0.0;
  }
  if(m_simplex_mesh != nullptr)
  {
    const TInt Id1 =  m_simplex_mesh->m_tet_nodes[m_simplexId][0];
    const TInt Id2 =  m_simplex_mesh->m_tet_nodes[m_simplexId][1];
    const TInt Id3 =  m_simplex_mesh->m_tet_nodes[m_simplexId][2];
    const TInt Id4 =  m_simplex_mesh->m_tet_nodes[m_simplexId][3];

    const math::Point p1 = m_simplex_mesh->m_coords[Id1];
    const math::Point p2 = m_simplex_mesh->m_coords[Id2];
    const math::Point p3 = m_simplex_mesh->m_coords[Id3];
    const math::Point p4 = m_simplex_mesh->m_coords[Id4];

    const math::VectorND<3, double> e1(p2 - p1);
    const math::VectorND<3, double> e2(p3 - p1);
    const math::VectorND<3, double> e3(p4 - p1);

    const Eigen::Vector3d v1(e1[0], e1[1], e1[2]);
    const Eigen::Vector3d v2(e2[0], e2[1], e2[2]);
    const Eigen::Vector3d v3(e3[0], e3[1], e3[2]);

    Eigen::Matrix3d m;
    m.col(0) = v1; m.col(1) = v2; m.col(2) = v3;

    return 1.0 / 6.0 *(v3.dot(v1.cross(v2)));
  }
  else
  {
      return 0.0;
  }
}
/*----------------------------------------------------------------------------*/
std::vector<TSimplexID> SimplicesCell::neighborTetra(const TInt indexNodeGlobal, bool boundariesAccepted) const
{
  std::vector<TSimplexID> v{};
  TInt sizeTabIndex = 4;
  TSimplexID errorId = std::numeric_limits<int>::min();
  unsigned int cmptBoundaries = 0;
  if(m_simplex_mesh != nullptr)
  {
    for(int i = 0; i < sizeTabIndex; i++)
    {
      TSimplexID indexGlobalInTetra = m_simplex_mesh->m_tet_nodes[m_simplexId][i];
      if(indexGlobalInTetra != indexNodeGlobal)
      {
        TInt tetraAdj =  m_simplex_mesh->getOppositeCell(indexGlobalInTetra, m_simplexId);
        if(tetraAdj == errorId)
        {
          cmptBoundaries++;
        }
        if(tetraAdj  != errorId)
        {
          //cette boucle if permet de gerer en plus les triangle adjacent aux faces...
          if(tetraAdj < 0)
          {
            TSimplexID triAdj = tetraAdj;
            bool flag0 = (std::find(v.begin(), v.end(), m_simplex_mesh->m_tri_nodes[-triAdj][3]) == v.end()
                          && m_simplex_mesh->m_tri_nodes[-triAdj][3] != m_simplexId);

            TSimplexID tetraAdj0  = m_simplex_mesh->m_tri_nodes[-triAdj][3];

            bool flag1 = (std::find(v.begin(), v.end(), m_simplex_mesh->m_tri_adj[-triAdj][3]) == v.end()
                          && m_simplexId != m_simplex_mesh->m_tri_adj[-triAdj][3]);

            TSimplexID tetraAdj1  = m_simplex_mesh->m_tri_adj[-triAdj][3];

            if(flag0 ==true && tetraAdj0 != errorId){v.push_back(tetraAdj0);}
            if(flag1 ==true && tetraAdj1 != errorId){v.push_back(tetraAdj1);}
          }
          else //tetraAdj >= 0
          {
              v.push_back(tetraAdj);
          }
        }
      }
    }
  }
  if(boundariesAccepted && cmptBoundaries > 0)
  {
    v.push_back(errorId);
  }
  return std::move(v);
}
/******************************************************************************/
std::vector<TSimplexID> SimplicesCell::neighborTri(const TInt indexNodeGlobal) const
{
  std::vector<TSimplexID> v{};
  TInt sizeTabIndex = 4;
  TSimplexID errorId = std::numeric_limits<int>::min();

  if(m_simplex_mesh != nullptr)
  {
    for(int i = 0; i < sizeTabIndex; i++)
    {
      TSimplexID indexGlobalInTetra = m_simplex_mesh->m_tet_nodes[m_simplexId][i];

      if(indexGlobalInTetra != indexNodeGlobal)
      {
        TInt triAdj =  m_simplex_mesh->getOppositeCell(indexGlobalInTetra, m_simplexId);
        if(triAdj  != errorId)
        {
          //cette boucle if permet de gerer en plus les triangle adjacent aux faces...
          if(triAdj < 0 && triAdj != 0)
          {
            v.push_back(triAdj);
          }
        }
      }
    }
  }

  return std::move(v);
}
/******************************************************************************/
std::vector<TSimplexID> SimplicesCell::neighborTri(const SimplicesNode& simpliceNode) const
{
  neighborTri(simpliceNode.getGlobalNode());
}
/*----------------------------------------------------------------------------*/
std::vector<TSimplexID> SimplicesCell::neighborTri() const
{

}
/*----------------------------------------------------------------------------*/
bool SimplicesCell::containNode(const simplicesNode::SimplicesNode& simplicesNode) const
{
  bool flag = false;
  TInt sizeTabIndexCell = 4;


  if(m_simplex_mesh != nullptr)
  {
    //changer le 4 pour quand il y'aura autre chose que simplement des tetra....

    for(int i = 0; i < sizeTabIndexCell ; i++)
    {
      flag = (simplicesNode.getGlobalNode() == m_simplex_mesh->m_tet_nodes[m_simplexId][i])?true:false;
      if(flag)
      {
        return flag;
      }
    }
  }

  return flag;
}
/*----------------------------------------------------------------------------*/
bool SimplicesCell::containNode(const simplicesNode::SimplicesNode& simplicesNode, TInt& indexLocal) const
{
  bool flag = false;
  TInt sizeTabIndexCell = 4;


  if(m_simplex_mesh != nullptr)
  {
    //changer le 4 pour quand il y'aura autre chose que simplement des tetra....

    for(int i = 0; i < sizeTabIndexCell ; i++)
    {
      flag = (simplicesNode.getGlobalNode() == m_simplex_mesh->m_tet_nodes[m_simplexId][i])?true:false;
      if(flag)
      {
        indexLocal = i;
        return flag;
      }
    }
  }

  return flag;
}
/*----------------------------------------------------------------------------*/
void SimplicesCell::reorientTet()
{
  if(m_simplex_mesh != nullptr)
  {
    if(m_simplex_mesh->m_tet_ids[m_simplexId] != 0)
    {
      double signedVolumeOfTet = getVolumeOfCell();
      if(signedVolumeOfTet < 0.0)
      {
        TInt indexTmp = m_simplex_mesh->m_tet_nodes[m_simplexId][2];
        m_simplex_mesh->m_tet_nodes[m_simplexId][2] = m_simplex_mesh->m_tet_nodes[m_simplexId][1];
        m_simplex_mesh->m_tet_nodes[m_simplexId][1] = indexTmp;
      }
    }
  }
}
/*----------------------------------------------------------------------------*/
double SimplicesCell::signedBarycentricNormalized(const TInt index, const gmds::math::Point& pt)
{
  double signedBarNormalized = 0.0;
  double volumeCell          = std::fabs(getVolumeOfCell());
  if(m_simplex_mesh->m_tet_ids[m_simplexId] != 0)
  {
    if(!(volumeCell <= 0))
    {
        signedBarNormalized  = signedBarycentric(index, pt) / volumeCell;
    }
    else
    {
      /*THROW*/
    }
  }

  return signedBarNormalized;
}
/*----------------------------------------------------------------------------*/
double SimplicesCell::signedBarycentric(const TInt index, const gmds::math::Point& pt)
{
  /*extract the orientation of the face seen by index*/
  int FacesOrientation[4][3] = {
    {1, 3, 2},
    {0, 2, 3},
    {1, 0, 3},
    {0, 1, 2},
  };

  double signedBar = 0.0;
  if(!(index < 0 || index > 3))
  {
    const TInt Node0 = m_simplex_mesh->m_tet_nodes[m_simplexId][FacesOrientation[index][0]];
    const TInt Node1 = m_simplex_mesh->m_tet_nodes[m_simplexId][FacesOrientation[index][1]];
    const TInt Node2 = m_simplex_mesh->m_tet_nodes[m_simplexId][FacesOrientation[index][2]];

    const math::Point pt0 = m_simplex_mesh->m_coords[Node0];
    const math::Point pt1 = m_simplex_mesh->m_coords[Node1];
    const math::Point pt2 = m_simplex_mesh->m_coords[Node2];

    const math::VectorND<3, double> e0(pt1 - pt0);
    const math::VectorND<3, double> e1(pt2 - pt0);

    const Eigen::Vector3d v1(e0[0], e0[1], e0[2]);
    const Eigen::Vector3d v2(e1[0], e1[1], e1[2]);

    const Eigen::Vector3d normalTet = v1.cross(v2);

    const Eigen::Vector3d v3(pt[0] - pt0[0], pt[1] - pt0[1], pt[2] - pt0[2]);
    signedBar =  normalTet.dot(v3) * 1.0 / 6.0;

    /*std::cout << "v3 --> " << v3[0] <<" " << v3[1] <<" " << v3[2] <<" " << std::endl;
    std::cout << "normalTet --> " << normalTet[0] <<" " << normalTet[1] <<" " << normalTet[2] <<" " << std::endl;
    std::cout << std::endl;*/
  }
  else
  {
    /*TODO THROW*/
  }
  return signedBar;
}
/*----------------------------------------------------------------------------*/
math::Vector3d SimplicesCell::normalOfFace(const std::vector<TInt>& nodes)
{
  math::Vector3d normal(0.0, 0.0, 0.0);
  if(nodes.size() != 3)
  {
    //TODO exception
    std::cout << "cell " << m_simplexId << " doesnt have a face of "<< nodes.size() << std::endl;
    return normal;
  }
  else
  {
    if(containNodes(nodes))
    {
      TInt localNode0 = getLocalNode(nodes[0]);
      TInt localNode1 = getLocalNode(nodes[1]);
      TInt localNode2 = getLocalNode(nodes[2]);

      math::Point pt0 = m_simplex_mesh->m_coords[m_simplex_mesh->m_tet_nodes[m_simplexId][localNode0]];
      math::Point pt1 = m_simplex_mesh->m_coords[m_simplex_mesh->m_tet_nodes[m_simplexId][localNode1]];
      math::Point pt2 = m_simplex_mesh->m_coords[m_simplex_mesh->m_tet_nodes[m_simplexId][localNode2]];

      math::Vector3d vec0 = math::Vector3d(pt2.X() - pt0.X(), pt2.Y() - pt0.Y(), pt2.Z() - pt0.Z());
      math::Vector3d vec1 = math::Vector3d(pt1.X() - pt0.X(), pt1.Y() - pt0.Y(), pt1.Z() - pt0.Z());

      normal = vec0.cross(vec1);
      std::vector<TInt>&& otherNode = getOtherNodeInSimplex(nodes);
      if(otherNode.size() != 1)
      {
        //TODO ECXEPTION
        std::cout << "cell " << m_simplexId << " doesnt have more or less than 4 nodes" << std::endl;
        return normal;
      }
      else
      {
        TInt localNode4 = getLocalNode(otherNode[0]);
        if(localNode4 == localNode0 || localNode4 == localNode1 || localNode4 == localNode2)
        {
          //TODO ECXEPTION
          std::cout << "cell " << m_simplexId << " doesnt have the same local Node" << std::endl;
          return normal;
        }
        else
        {
          math::Point ptTest = m_simplex_mesh->m_coords[m_simplex_mesh->m_tet_nodes[m_simplexId][localNode4]];
          math::Vector3d vecTest = math::Vector3d(ptTest.X() - pt0.X(), ptTest.Y() - pt0.Y(), ptTest.Z() - pt0.Z());
          if(normal.dot(vecTest) < 0.0)
          {
            normal = -normal;
          }
        }
      }
    }
    else
    {
      std::cout << "cell " << m_simplexId << " don't contain all of nodes" << std::endl;
      return normal;
    }
  }
}
/*----------------------------------------------------------------------------*/
bool SimplicesCell::correspondance(const TInt localIndex, const TInt generalIndex) const
{
  bool flag = false;
  if(m_simplex_mesh->m_node_ids[generalIndex] == true)
  {
    flag = (m_simplex_mesh->m_tet_nodes[m_simplexId][localIndex] == generalIndex);
  }
  else
  {
    flag = false;
  }
  return flag;
}
/*----------------------------------------------------------------------------*/
TInt SimplicesCell::getLocalNode(const TInt generalIndex) const
{
  int errorId = std::numeric_limits<int>::min();
  TInt localNode = errorId;
  size_t sizeSimplexNodes = 4;
  for(TInt iter = 0; iter < sizeSimplexNodes ; iter++)
  {
    if(m_simplex_mesh->m_tet_nodes[m_simplexId][iter] == generalIndex)
    {
      localNode = iter;
      break;
    }
  }

  return localNode;
}
/*----------------------------------------------------------------------------*/
TSimplexID SimplicesCell::oppositeTetraIdx(const simplicesNode::SimplicesNode& simplicesNode) const
{
  size_t sizeLocalIndex = 4;
  int errorId = std::numeric_limits<int>::min();
  TSimplexID oppositeSimplex;

  if(containNode(simplicesNode) && m_simplex_mesh->m_node_ids[simplicesNode.getGlobalNode()])
  {
    for(TInt localIndex = 0; localIndex < sizeLocalIndex; localIndex++)
    {
      if(correspondance(localIndex, simplicesNode.getGlobalNode()))
      {
          oppositeSimplex = m_simplex_mesh->m_tet_adj[m_simplexId][localIndex];
      }
    }
  }

  if(oppositeSimplex < 0 && oppositeSimplex != errorId)
  {
    oppositeSimplex = (m_simplexId == m_simplex_mesh->m_tri_nodes[-oppositeSimplex][3])?m_simplex_mesh-> m_tri_adj[-oppositeSimplex][3]: m_simplex_mesh->m_tri_nodes[-oppositeSimplex][3];
  }
  return oppositeSimplex;
}
/******************************************************************************/
TSimplexID SimplicesCell::oppositeTetraIdx(const TInt indexLocal) const
{
  TSimplexID oppositeSimplex;
  if(!(indexLocal < 0 || indexLocal > 3))
  {
      oppositeSimplex = m_simplex_mesh->m_tet_adj[m_simplexId][indexLocal];
  }
  else
  {
    /*TODO exeption*/
  }

  return oppositeSimplex;

}
/******************************************************************************/
std::vector<TSimplexID> SimplicesCell::oppositeTetraVectorPrivated(const SimplicesNode& simplicesNode) const
{
  std::vector<TSimplexID> v{};
  v.resize(3);
  if(m_simplex_mesh != nullptr)
  {
    unsigned int nodeLocalSize = 4;
    unsigned int cpt = 0;
    for(unsigned int nodeLocal = 0; nodeLocal < nodeLocalSize ; nodeLocal++)
    {
      if(!correspondance(nodeLocal,simplicesNode.getGlobalNode()))
      {
        v[cpt] = m_simplex_mesh->m_tet_adj[m_simplexId][nodeLocal];
        cpt++;
      }
    }
  }

  return std::move(v);

}
/******************************************************************************/
std::vector<TSimplexID> SimplicesCell::adjacentTetra()
{
  std::vector<TSimplexID> v{};
  v.resize(4);
  unsigned int nodeLocalSize = 4;

  for(unsigned int nodeLocal = 0; nodeLocal < nodeLocalSize; nodeLocal++)
  {
    v[nodeLocal] = m_simplex_mesh->m_tet_adj[m_simplexId][nodeLocal];
  }

  return std::move(v);
}
/******************************************************************************/
std::vector<TInt> SimplicesCell::intersectionNodes(const SimplicesCell& simplicesCell)
{
  std::vector<TInt> v{};
  if(m_simplex_mesh != nullptr)
  {
    std::vector<TInt> nodesT0{
      getNode(0).getGlobalNode(),
      getNode(1).getGlobalNode(),
      getNode(2).getGlobalNode(),
      getNode(3).getGlobalNode()};

    for(unsigned int nodeLocalT0 = 0; nodeLocalT0 < nodesT0.size(); nodeLocalT0++)
    {
      for(unsigned int nodeLocalT1 = 0; nodeLocalT1 < nodesT0.size(); nodeLocalT1++)
      {
        if(nodesT0[nodeLocalT0] == simplicesCell.getNode(nodeLocalT1).getGlobalNode())
        {
          v.push_back(nodesT0[nodeLocalT0]);
        }
      }
    }
  }
  else
  {
      /*TODO exeption*/
  }

  return std::move(v);
}
/*----------------------------------------------------------------------------*/
SimplicesNode  SimplicesCell::getNode(const TInt indexLocal) const
{
  return SimplicesNode(m_simplex_mesh, m_simplex_mesh->m_tet_nodes[m_simplexId][indexLocal]);
}
/*----------------------------------------------------------------------------*/
std::vector<unsigned int> SimplicesCell::nodes() const
{
    std::vector<unsigned int> v{};
    v.resize(4);

    v[0] = (unsigned int)m_simplex_mesh->m_tet_nodes[m_simplexId][0];
    v[1] = (unsigned int)m_simplex_mesh->m_tet_nodes[m_simplexId][1];
    v[2] = (unsigned int)m_simplex_mesh->m_tet_nodes[m_simplexId][2];
    v[3] = (unsigned int)m_simplex_mesh->m_tet_nodes[m_simplexId][3];

    return std::move(v);
}
/*----------------------------------------------------------------------------*/
bool SimplicesCell::inBorder() const
{
  bool flag = false;
  unsigned int simplexNodeSize = 4;
  int border = std::numeric_limits<int>::min();

  for(TInt node = 0; node < simplexNodeSize; node++)
  {
    if(m_simplex_mesh->m_tet_adj[m_simplexId][node] == border)
    {
      flag = true;
      break;
    }
  }
  return flag;
}
/*----------------------------------------------------------------------------*/
bool SimplicesCell::containNodes(const std::vector<TInt>& simplicesNode) const
{
  std::vector<simplicesNode::SimplicesNode> nodes{};
  for(const auto & node : simplicesNode)
  {
    nodes.push_back(SimplicesNode(m_simplex_mesh, node));
  }

  return containNodes(nodes);
}
/*----------------------------------------------------------------------------*/
bool SimplicesCell::containNodes(const std::vector<simplicesNode::SimplicesNode>& simplicesNode) const
{
  bool flag = true;
  for(auto const & node : simplicesNode)
  {
    if(containNode(node) == false)
    {
      flag = false;
      break;
    }
  }
  return flag;
}
/*----------------------------------------------------------------------------*/
bool SimplicesCell::containAtLeast(TInt N, const std::vector<simplicesNode::SimplicesNode>& simplicesNode) const
{
  bool flag = false;
  size_t cmpt = 0;

  for(auto const & node : simplicesNode)
  {
    if(containNode(node) == true)
    {
      cmpt++;
    }
  }

  return (cmpt >= N)? true : false;
}
/*----------------------------------------------------------------------------*/
double SimplicesCell::signedBarycentric(const gmds::math::Point& pt, const std::vector<SimplicesNode>& nodes)
{
  if(nodes.size() > 3)
  {
    /*TODO*/
    //exception
    std::cout << "nodes size is higher than 3 can not identifiate the face " << std::endl;
  }
  else if(!containAtLeast(nodes.size(), nodes))
  {
    /*TODO*/
    //exception
    std::cout << "one of the node is not on part of the tetra : " << m_simplexId << std::endl;
  }

  TInt   indexLocal = 0;
  size_t sizeNodeInSimplex = 4;
  std::vector<TInt> nodeGlobal{nodes[0].getGlobalNode(), nodes[1].getGlobalNode(), nodes[2].getGlobalNode()};

  for(TInt node = 0  ; node < sizeNodeInSimplex ; node++)
  {
    TInt globalNode = m_simplex_mesh->m_tet_nodes[m_simplexId][node];
    if(std::find(nodeGlobal.begin(), nodeGlobal.end(), globalNode) == nodeGlobal.end())
    {
      indexLocal = node;
      break;
    }
  }

  return signedBarycentric(indexLocal, pt);
}
/*----------------------------------------------------------------------------*/
std::vector<SimplicesNode> SimplicesCell::removeNodeNotInSimplex(const std::vector<SimplicesNode>& nodes) const
{
  std::vector<TInt> cptIndx;
  std::vector<SimplicesNode> resNodes;

  for(TInt idxNode = 0; idxNode < nodes.size() ; idxNode++)
  {
    const SimplicesNode & node = nodes[idxNode];
    if(containNode(node) == true)
    {
      cptIndx.push_back(idxNode);
    }
  }

  for(auto const & idx : cptIndx)
  {
    resNodes.push_back(nodes[idx]);
  }

  return std::move(resNodes);
}
/*----------------------------------------------------------------------------*/
SimplicesNode SimplicesCell::removeNodesFromVec(std::vector<SimplicesNode>& nodes) const
{
  SimplicesNode res;
  if(m_simplex_mesh != nullptr)
  {
    if(m_simplex_mesh->m_tet_ids[m_simplexId] != 0)
    {
      if(containNodes(nodes) == false)
      {
        /*TODO ... exception*/
        std::cout << "One Or more node is not on the tetra : " << *this << std::endl;
      }
      SimplicesNode nodeA(m_simplex_mesh, m_simplex_mesh->m_tet_nodes[m_simplexId][0]);
      SimplicesNode nodeB(m_simplex_mesh, m_simplex_mesh->m_tet_nodes[m_simplexId][1]);
      SimplicesNode nodeC(m_simplex_mesh, m_simplex_mesh->m_tet_nodes[m_simplexId][2]);
      SimplicesNode nodeD(m_simplex_mesh, m_simplex_mesh->m_tet_nodes[m_simplexId][3]);

      std::vector<SimplicesNode> vecNodes{nodeA, nodeB, nodeC, nodeD};

      for(auto const nodeInVec : vecNodes)
      {
        if(std::find(nodes.begin(), nodes.end(), nodeInVec) == nodes.end())
        {
          res = nodeInVec;
          break;
        }
      }
    }
    else
    {
      /*TODO exception*/
    }
  }

  return std::move(res);
}
/*----------------------------------------------------------------------------*/
std::vector<TInt> SimplicesCell::getOtherNodeInSimplex(const std::vector<TSimplexID>& v) const
{
  std::vector<TInt> res{};
  if(m_simplex_mesh != nullptr)
  {
    std::vector<TInt> nodeInTet{m_simplex_mesh->m_tet_nodes[m_simplexId][0],
                                m_simplex_mesh->m_tet_nodes[m_simplexId][1],
                                m_simplex_mesh->m_tet_nodes[m_simplexId][2],
                                m_simplex_mesh->m_tet_nodes[m_simplexId][3]};

    for(auto const & nodeIndx : nodeInTet)
    {
      if(std::find(v.begin(), v.end(), nodeIndx) == v.end())
      {
        res.push_back(nodeIndx);
      }
    }

    return std::move(res);
  }


}
