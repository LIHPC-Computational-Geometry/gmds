/*----------------------------------------------------------------------------*/
#include <gmds/hybridMeshAdapt/SimplicesCell.h>
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
/*----------------------------------------------------------------------------*/
#include <chrono>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace hybrid;
using namespace simplicesNode;
using namespace simplicesCell;
using namespace operators;
/*----------------------------------------------------------------------------*/
SimplicesCell::SimplicesCell(SimplexMesh* simplexMesh, const TSimplexID simplexId)
{
    m_simplex_mesh = simplexMesh;
    m_simplexId    = simplexId;
    math::Orientation::initialize();
    if(m_simplex_mesh->m_tet_ids[m_simplexId] == 0)
    {
      /*TODO exeption le node d'existe pas ... le creer avant de poouvoir l'utiliser*/
      std::cout << "Creer la cellule " << m_simplexId <<  " avant de l'utiliser !!" << std::endl;
      throw gmds::GMDSException("cell doesnt exist") ;
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
  math::Orientation::finalize();
}
/*----------------------------------------------------------------------------*/
std::vector<TInt> SimplicesCell::getNodes() const
{
  std::vector<TInt> v{};
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
double SimplicesCell::getVolumeOfCell() const
{
  if(m_simplexId > m_simplex_mesh->m_tet_nodes.size() || m_simplexId < 0)
  {
    //todo exception
    return 0.0;
  }
  if(m_simplex_mesh != nullptr)
  {
    /*const TInt Id1 =  m_simplex_mesh->m_tet_nodes[m_simplexId][0];
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

    //Eigen::Matrix3d m;
    //m.col(0) = v1; m.col(1) = v2; m.col(2) = v3;

    return 1.0 / 6.0 *(v3.dot(v1.cross(v2)));*/


    const TInt Node0 = m_simplex_mesh->m_tet_nodes[m_simplexId][0]; //1
    const TInt Node1 = m_simplex_mesh->m_tet_nodes[m_simplexId][1]; //3
    const TInt Node2 = m_simplex_mesh->m_tet_nodes[m_simplexId][2]; //2
    const TInt Node3 = m_simplex_mesh->m_tet_nodes[m_simplexId][3];

    const math::Point pt0 = m_simplex_mesh->m_coords[Node0];
    const math::Point pt1 = m_simplex_mesh->m_coords[Node1];
    const math::Point pt2 = m_simplex_mesh->m_coords[Node2];
    const math::Point pt3 = m_simplex_mesh->m_coords[Node3];

    const math::VectorND<3, double> e0(pt1 - pt0);
    const math::VectorND<3, double> e1(pt2 - pt0);
    const math::VectorND<3, double> e2(pt3 - pt0);

    const Eigen::Vector3d v1(e0[0], e0[1], e0[2]);
    const Eigen::Vector3d v2(e1[0], e1[1], e1[2]);
    const Eigen::Vector3d v3(e2[0], e2[1], e2[2]);

    Eigen::Matrix3d mat;
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;

    return 1.0/ 6.0 * mat.determinant();
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
        //TSimplexID tetraAdj =  m_simplex_mesh->getOppositeCell(indexGlobalInTetra, m_simplexId);
        TSimplexID tetraAdj =  oppositeTetraIdx(i);

        if(tetraAdj == errorId)
        {
          cmptBoundaries++;
        }
        if(tetraAdj  != errorId)
        {
          //cette boucle if permet de gerer en plus les triangle adjacent aux faces...
          if(tetraAdj < 0)
          {
            if(m_simplex_mesh->m_tri_ids[-tetraAdj] != 0)
            {
                v.push_back(tetraAdj);
            }
          }
          else //tetraAdj >= 0
          {
            if(m_simplex_mesh->m_tet_ids[tetraAdj] != 0)
            {
                v.push_back(tetraAdj);
            }
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
  std::vector<TSimplexID> res{};
  neighborTri(simpliceNode.getGlobalNode());

  return res;
}
/*----------------------------------------------------------------------------*/
std::vector<TSimplexID> SimplicesCell::neighborTri() const
{
  std::vector<TSimplexID> res{};
  return res;
}
/*----------------------------------------------------------------------------*/
bool SimplicesCell::containNode(const simplicesNode::SimplicesNode& simplicesNode) const
{
  bool flag = false;
  TInt sizeTabIndexCell = 4;


  if(m_simplex_mesh != nullptr)
  {
    for(int i = 0; i < sizeTabIndexCell ; i++)
    {
      flag = (simplicesNode.getGlobalNode() == m_simplex_mesh->m_tet_nodes[m_simplexId][i])?true:false;
      if(flag)
      {
        break;;
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
std::vector<double> SimplicesCell::uvwt(const gmds::math::Point& pt) const
{
  std::vector<double> uvwt;
  unsigned int sizeFACE = 4;
  for(unsigned int face = 0 ; face < sizeFACE ; face++)
    uvwt.push_back(signedBarycentricNormalized(face, pt));

  return uvwt;
}
/*----------------------------------------------------------------------------*/
double SimplicesCell::signedBarycentricNormalized(const TInt index, const gmds::math::Point& pt) const
{
  double signedBarNormalized = 0.0;
  double volumeCell          = std::fabs(getVolumeOfCell());
  if(m_simplex_mesh->m_tet_ids[m_simplexId] != 0)
  {
    if((volumeCell != 0.0))
    {
        signedBarNormalized  = signedBarycentric(index, pt) / volumeCell;
    }
  }

  return signedBarNormalized;
}
/*----------------------------------------------------------------------------*/
double SimplicesCell::signedBarycentric(const TInt index, const gmds::math::Point& pt) const
{
  /*extract the orientation of the face seen by index*/
  double signedBar = 0.0;
  if(!(index < 0 || index > 3))
  {
    const TInt Node0 = m_simplex_mesh->m_tet_nodes[m_simplexId][FacesOrientation[index][0]]; //1
    const TInt Node1 = m_simplex_mesh->m_tet_nodes[m_simplexId][FacesOrientation[index][1]]; //3
    const TInt Node2 = m_simplex_mesh->m_tet_nodes[m_simplexId][FacesOrientation[index][2]]; //2

    const math::Point pt0 = m_simplex_mesh->m_coords[Node0];
    const math::Point pt1 = m_simplex_mesh->m_coords[Node1];
    const math::Point pt2 = m_simplex_mesh->m_coords[Node2];

    const math::VectorND<3, double> e0(pt1 - pt0);
    const math::VectorND<3, double> e1(pt2 - pt0);

    const Eigen::Vector3d v1(e0[0], e0[1], e0[2]);
    const Eigen::Vector3d v2(e1[0], e1[1], e1[2]);

    Eigen::Vector3d normalTet = v1.cross(v2);
    //normalTet.normalize();

    Eigen::Vector3d v3(pt[0] - pt0[0], pt[1] - pt0[1], pt[2] - pt0[2]);
    signedBar =  normalTet.dot(v3) * 1.0 / 6.0;
  }
  else
  {
    /*TODO THROW*/
  }
  return signedBar;
}

/*----------------------------------------------------------------------------*/
bool SimplicesCell::isInCell(const gmds::math::Point& pt, bool inverseOrientation) const
{
  math::Orientation::Sign sign;
  unsigned int sizeCell = 4;
  for(unsigned int face = 0 ; face < sizeCell ; face++)
  {
    const TInt Node0 = m_simplex_mesh->m_tet_nodes[m_simplexId][FacesOrientation[face][0]];
    const TInt Node1 = m_simplex_mesh->m_tet_nodes[m_simplexId][FacesOrientation[face][1]];
    const TInt Node2 = m_simplex_mesh->m_tet_nodes[m_simplexId][FacesOrientation[face][2]];


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

    if(sign == math::Orientation::Sign::NEGATIVE)
    {
      return false;
    }
  }

  return true;
}
/*----------------------------------------------------------------------------*/
bool SimplicesCell::isCellClose(const gmds::math::Point& pt, std::vector<double>& UVWT, double epsilon) const
{
  UVWT =  uvwt(pt);

  for(double c : UVWT)
    if(c < - epsilon)
      return false;

  return true;
}
/*----------------------------------------------------------------------------*/
math::Orientation::Sign SimplicesCell::orientation(const TInt faceIdx, const gmds::math::Point& pt, bool inverseOrientation) const
{
  math::Orientation::Sign sign;
  if(!(faceIdx < 0 || faceIdx > 3))
  {
    const TInt Node0 = m_simplex_mesh->m_tet_nodes[m_simplexId][FacesOrientation[faceIdx][0]];
    const TInt Node1 = m_simplex_mesh->m_tet_nodes[m_simplexId][FacesOrientation[faceIdx][1]];
    const TInt Node2 = m_simplex_mesh->m_tet_nodes[m_simplexId][FacesOrientation[faceIdx][2]];


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
  }
  /*TODO assertion if faceIdx <0 | faceIdx >3*/
  return sign;
}
/*----------------------------------------------------------------------------*/
double SimplicesCell::dihedralAngle(const unsigned int localNode0, const unsigned int localNode1) const
{
    if((localNode0 < 0 || localNode0 > 3) || (localNode1 < 0 || localNode1 > 3))
    {
      std::cout << "localNode0 -> " << localNode0 << std::endl;
      std::cout << "localNode1 -> " << localNode1 << std::endl;
      throw gmds::GMDSException("(localNode0 < 0 || localNode0 > 3) || (localNode1 < 0 || localNode1 > 3)");
    }
    if(localNode0 == localNode1)
    {
      throw gmds::GMDSException("localNode0 == localNode1");
    }

    std::vector<unsigned int> localsNode{0,1,2,3};
    localsNode.erase(std::remove_if(localsNode.begin(), localsNode.end(), [=](unsigned int lN){
      return (lN == localNode0 || lN == localNode1);
    }), localsNode.end());

    const math::Point pt0 = getNode(localNode0).getCoords();
    const math::Point pt1 = getNode(localNode1).getCoords();
    const math::Point pt2 = getNode(localsNode.front()).getCoords();
    const math::Point pt3 = getNode(localsNode.back()).getCoords();

    return gmds::math::dihedralAngle(pt0, pt1, pt2, pt3);
}
/*----------------------------------------------------------------------------*/
std::set<double> SimplicesCell::minAndmaxDihedralAngle() const
{
  std::set<double> s{};
  s.insert(dihedralAngle(0, 1));
  s.insert(dihedralAngle(0, 2));
  s.insert(dihedralAngle(0, 3));
  s.insert(dihedralAngle(1, 2));
  s.insert(dihedralAngle(1, 3));
  s.insert(dihedralAngle(2, 3));

  return s;
}
/*----------------------------------------------------------------------------*/
bool SimplicesCell::isSliver() const
{
  //a checker si cela convient
  double epsilon = 0.1;
  double pi      = M_PI;
  const std::set<double> s = minAndmaxDihedralAngle();

  return  (*(s.begin()) < epsilon) || (*(--s.end()) > (M_PI - epsilon));
}
/*----------------------------------------------------------------------------*/
math::Vector3d SimplicesCell::normalOfFace(const std::vector<TInt>& nodes) const
{
  math::Vector3d normal({0.0, 0.0, 0.0});
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

      math::Vector3d vec0 = math::Vector3d({pt2.X() - pt0.X(), pt2.Y() - pt0.Y(), pt2.Z() - pt0.Z()});
      math::Vector3d vec1 = math::Vector3d({pt1.X() - pt0.X(), pt1.Y() - pt0.Y(), pt1.Z() - pt0.Z()});

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
          math::Vector3d vecTest = math::Vector3d({ptTest.X() - pt0.X(), ptTest.Y() - pt0.Y(), ptTest.Z() - pt0.Z()});
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

  return normal;
}
/*----------------------------------------------------------------------------*/
std::vector<math::Vector3d> SimplicesCell::normalOfFaces(const std::vector<std::vector<TInt>>& FacesNodes)
{
  std::vector<math::Vector3d> normalFaces;
  if(m_simplex_mesh != nullptr)
  {
    for(auto const & FaceNodes : FacesNodes)
    {
      math::Vector3d normalFace = normalOfFace(FaceNodes);
      normalFaces.push_back(normalFace);
    }
  }

  return std::move(normalFaces);
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
  int border = std::numeric_limits<int>::min();
  TSimplexID oppositeSimplex = border;

  if(containNode(simplicesNode) && m_simplex_mesh->m_node_ids[simplicesNode.getGlobalNode()] != 0)
  {
    for(TInt localIndex = 0; localIndex < sizeLocalIndex; localIndex++)
    {
      if(correspondance(localIndex, simplicesNode.getGlobalNode()))
      {
        oppositeSimplex = m_simplex_mesh->m_tet_adj[m_simplexId][localIndex];
        if(oppositeSimplex != border)
        {
          oppositeSimplex = (m_simplex_mesh->m_tet_ids[oppositeSimplex] != 0 )? oppositeSimplex : border;
        }
      }
    }
  }

  if(oppositeSimplex < 0 && oppositeSimplex != border)
  {
    oppositeSimplex = (m_simplexId == m_simplex_mesh->m_tri_nodes[-oppositeSimplex][3])?m_simplex_mesh-> m_tri_adj[-oppositeSimplex][3]: m_simplex_mesh->m_tri_nodes[-oppositeSimplex][3];
  }
  return oppositeSimplex;
}
/******************************************************************************/
TSimplexID SimplicesCell::oppositeTetraIdx(const TInt indexLocal) const
{
  TSimplexID oppositeSimplex;
  TSimplexID border = std::numeric_limits<TSimplexID>::min();
  oppositeSimplex = border;

  if(!(indexLocal < 0 || indexLocal > 3))
  {
    oppositeSimplex = m_simplex_mesh->m_tet_adj[m_simplexId][indexLocal];
    if(oppositeSimplex != border)
    {
      if(oppositeSimplex >= 0)
      {
        if(m_simplex_mesh->m_tet_ids[oppositeSimplex] != 1)
        {
          oppositeSimplex = border;
        }
      }
      else
      {
        if(m_simplex_mesh->m_tri_ids[-oppositeSimplex] != 1)
        {
          oppositeSimplex = border;
        }
      }
    }
  }
  else
  {
    /*TODO exeption*/
  }

  return oppositeSimplex;

}
/******************************************************************************/
std::vector<TInt> SimplicesCell::getOrderedFace(const TInt indexFace) const
{
  std::vector<TInt> v{};
  if(!(indexFace < 0 || indexFace > 3))
  {
    v.reserve(3);
    v.push_back(m_simplex_mesh->m_tet_nodes[m_simplexId][FacesOrientation[indexFace][0]]);
    v.push_back(m_simplex_mesh->m_tet_nodes[m_simplexId][FacesOrientation[indexFace][1]]);
    v.push_back(m_simplex_mesh->m_tet_nodes[m_simplexId][FacesOrientation[indexFace][2]]);
  }
  else
  {
    /*TODO assertion...*/
  }

  return std::move(v);
}
/******************************************************************************/
std::vector<TInt> SimplicesCell::visibleFaces(const math::Point& coordNode) const
{
  CriterionRAIS criterionRAIS(new VolumeCriterion());
  const unsigned int sizeFace = 4;
  std::vector<TInt> res{};
  for(unsigned int face = 0 ; face < sizeFace ; face++)
  {
    bool isfaceVisible = !criterionRAIS.execute(m_simplex_mesh, m_simplexId, face, coordNode);
    if(isfaceVisible)
    {
      res.push_back(face);
    }
  }
  return res;
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
std::vector<TSimplexID> SimplicesCell::adjacentTetra() const
{
  TInt border = std::numeric_limits<int>::min();
  std::vector<TSimplexID> v{border, border, border, border};
  unsigned int nodeLocalSize = 4;

  v[0] = m_simplex_mesh->m_tet_adj[m_simplexId][0];
  v[1] = m_simplex_mesh->m_tet_adj[m_simplexId][1];
  v[2] = m_simplex_mesh->m_tet_adj[m_simplexId][2];
  v[3] = m_simplex_mesh->m_tet_adj[m_simplexId][3];

  return std::move(v);
}
/******************************************************************************/
std::vector<TSimplexID> SimplicesCell::directConnectedSimplex () const
{
  gmds::BitVector tetAlreadyIn(m_simplex_mesh->getBitVectorTet().capacity());
  gmds::BitVector triAlreadyIn(m_simplex_mesh->getBitVectorTri().capacity());
  std::vector<TSimplexID> res{};

  std::vector<TInt> nodes = getNodes();
  for(auto const n : nodes)
  {
    std::vector<TSimplexID> ball = SimplicesNode(m_simplex_mesh, n).ballOf();
    for(auto const simplex : ball)
    {
      if(simplex >= 0)
      {
        if(tetAlreadyIn[simplex] == 0)
        {
          res.push_back(simplex);
          tetAlreadyIn.assign(simplex);
        }
      }
      else
      {
        if(triAlreadyIn[-simplex] == 0)
        {
          res.push_back(simplex);
          triAlreadyIn.assign(-simplex);
        }
      }
    }
  }
  return res;
}
/******************************************************************************/
std::vector<TInt> SimplicesCell::intersectionFaces(const std::vector<TInt>& localFaces) const
{
  std::vector<TInt> res{};
  if(localFaces.size() > 0)
  {
    if(localFaces.size() == 1)
    {
      res = getOrderedFace(localFaces.front());
    }
    else
    {
      std::vector<std::vector<TInt>> facesNodes{};
      for(auto const face : localFaces)
      {
        std::vector<TInt> orderedFaces = getOrderedFace(face);
        std::sort(orderedFaces.begin(), orderedFaces.end());
        facesNodes.push_back(orderedFaces);
      }

      if(localFaces.size() == 2)
      {
        std::set_intersection(facesNodes.front().begin(), facesNodes.front().end(), facesNodes.back().begin(), facesNodes.back().end(), std::back_inserter(res));
      }
      else if(localFaces.size() == 3)
      {
        std::vector<TInt> tmp;
        std::set_intersection(facesNodes[0].begin(), facesNodes[0].end(), facesNodes[1].begin(), facesNodes[1].end(), std::back_inserter(tmp));
        std::set_intersection(tmp.begin(), tmp.end(), facesNodes[2].begin(), facesNodes[2].end(), std::back_inserter(res));
      }
    }
  }
  return res;
}
/******************************************************************************/
std::vector<TInt> SimplicesCell::intersectionNodes(const SimplicesCell& simplicesCell) const
{
  std::vector<TInt> v{};
  if(m_simplex_mesh != nullptr)
  {
    //std::cout << "std::vector<TInt> nodesT0" << std::endl;
    std::vector<TInt> nodesT0{
      getNode(0).getGlobalNode(),
      getNode(1).getGlobalNode(),
      getNode(2).getGlobalNode(),
      getNode(3).getGlobalNode()
    };

    std::vector<TInt> nodesT1{
      simplicesCell.getNode(0).getGlobalNode(),
      simplicesCell.getNode(1).getGlobalNode(),
      simplicesCell.getNode(2).getGlobalNode(),
      simplicesCell.getNode(3).getGlobalNode()
    };



    for(unsigned int nodeLocalT0 = 0; nodeLocalT0 < nodesT0.size(); nodeLocalT0++)
    {
      for(unsigned int nodeLocalT1 = 0; nodeLocalT1 < nodesT1.size(); nodeLocalT1++)
      {
        if(nodesT0[nodeLocalT0] == nodesT1[nodeLocalT1])
        {
          v.push_back(nodesT0[nodeLocalT0]);
          break;
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
/******************************************************************************/
void SimplicesCell::intersectionSimplexFacesForUnbuildStruct(const SimplicesCell& simplicesCell, std::vector<TSimplexID>& intersectionNodes)
{
  if(m_simplex_mesh != nullptr)
  {
    std::vector<TInt> nodesT0{
      getNode(0).getGlobalNode(),
      getNode(1).getGlobalNode(),
      getNode(2).getGlobalNode(),
      getNode(3).getGlobalNode()
    };

    std::vector<TInt> nodesT1{
      simplicesCell.getNode(0).getGlobalNode(),
      simplicesCell.getNode(1).getGlobalNode(),
      simplicesCell.getNode(2).getGlobalNode(),
      simplicesCell.getNode(3).getGlobalNode()
    };

    unsigned int cpt = 0;
    for(unsigned int nodeLocalT0 = 0; nodeLocalT0 < nodesT0.size(); nodeLocalT0++)
    {
      unsigned int sizeV = intersectionNodes.size();
      for(unsigned int nodeLocalT1 = 0; nodeLocalT1 < nodesT1.size(); nodeLocalT1++)
      {
        if(nodesT0[nodeLocalT0] == nodesT1[nodeLocalT1])
        {
          intersectionNodes.push_back(nodesT0[nodeLocalT0]);
          break;
        }
      }

      if(intersectionNodes.size() == sizeV)
      {
        cpt++;
        if(cpt > 1)
        {
          break;
        }
      }
    }
  }
  else
  {
      /*TODO exeption*/
  }

}
/*----------------------------------------------------------------------------*/
bool SimplicesCell::intersectionSimplexFaces(const SimplicesCell& simplicesCell, std::vector<TInt>& simplicesNodeLocal)
{
  if(m_simplex_mesh != nullptr)
  {
    unsigned int nodeLocalNbr = 4;
    std::vector<TInt> nodesT0{
      getNode(0).getGlobalNode(),
      getNode(1).getGlobalNode(),
      getNode(2).getGlobalNode(),
      getNode(3).getGlobalNode()
    };

    for(unsigned int nodeLocal = 0; nodeLocal < nodeLocalNbr ; nodeLocal++)
    {
      TSimplexID oppositeCell  = oppositeTetraIdx(nodeLocal);
      if(oppositeCell == simplicesCell.simplexId())
      {
        simplicesNodeLocal.reserve(2);
        simplicesNodeLocal.push_back(nodeLocal);
        for(unsigned int nodeLocalBis = 0 ; nodeLocalBis < nodeLocalNbr ; nodeLocalBis++)
        {
          if(simplicesCell.oppositeTetraIdx(nodeLocalBis) == m_simplexId)
          {
            simplicesNodeLocal.push_back(nodeLocalBis);
            return true;
          }
        }
      }
    }
  }
  else
  {
      /*TODO exeption*/
  }
  return false;
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
double SimplicesCell::signedBarycentric(const gmds::math::Point& pt, const std::vector<SimplicesNode>& nodes) const
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

  }
  return std::move(res);

}
/*----------------------------------------------------------------------------*/
bool SimplicesCell::isPointInSimplex(const math::Point& pt) const
{
  bool flag = false;
  double u = signedBarycentric(0,pt);
  double v = signedBarycentric(1,pt);
  double w = signedBarycentric(2,pt);
  double t = signedBarycentric(3,pt);


  double epsilon = 0.0;
  if(u >= epsilon && v >= epsilon && w >= epsilon && t >= epsilon)
  {
    flag = true;
  }
  return flag;
}
/*----------------------------------------------------------------------------*/
bool SimplicesCell::isPointInSimplices(const math::Point& pt, double& u, double& v, double& w, double& t) const
{
  bool flag = false;
  u = signedBarycentric(0,pt);
  v = signedBarycentric(1,pt);
  w = signedBarycentric(2,pt);
  t = signedBarycentric(3,pt);

  double epsilon = 0.0;
  if(u >= epsilon && v >= epsilon && w >= epsilon && t >= epsilon)
  {
    flag = true;
  }
  return flag;
}
/*----------------------------------------------------------------------------*/
unsigned int SimplicesCell::checkFaceNbrVisibility(std::vector<std::vector<TInt>>& facesId, const SimplicesNode & simpliceNode)
{
  facesId.clear();
  unsigned int nbrFaceNotVisible = 0;
  std::vector<TInt> nodes = getNodes();
  if(nodes.size() == 4)
  {
    std::vector<std::vector<TInt>> faces{{nodes[0], nodes[1], nodes[2]}, {nodes[0], nodes[1], nodes[3]},
                                         {nodes[0], nodes[2], nodes[3]}, {nodes[1], nodes[2], nodes[3]}};
    std::vector<math::Vector3d> normals = normalOfFaces(faces);

    for(unsigned int iter = 0; iter < faces.size(); iter ++)
    {
      TInt nodeIdFace        = faces[iter][0];
      SimplicesNode nodeFace = SimplicesNode(m_simplex_mesh, nodeIdFace);
      math::Vector3d vec     = simpliceNode.getCoords() - nodeFace.getCoords();
      if(vec.dot(normals[iter]) < 0.0)
      {
        facesId.push_back(faces[iter]);
        nbrFaceNotVisible++;
      }
    }
  }
  else
  {
    //TODO exception
    std::cout << "nodes.size() != 4" << std::endl;
  }

  return nbrFaceNotVisible;
}
/*----------------------------------------------------------------------------*/
std::vector<double> SimplicesCell::getBestNodeJacobianNormalise()
{
  std::vector<double> ans {};
  std::vector<TInt> nodes = getNodes();
  for(unsigned int lN = 0 ; lN < 4 ; lN++)
  {
    unsigned int lN1 = (lN+1) % 4;
    unsigned int lN2 = (lN+2) % 4;
    unsigned int lN3 = (lN+3) % 4;

    //compute the jacobian for each node of the tetra
    TInt A = nodes[lN1] ; TInt B = nodes[lN2] ; TInt C = nodes[lN3];
    TInt P = nodes[lN];
    const math::Point pt = SimplicesNode(m_simplex_mesh, P).getCoords();
    const math::Point ptA = SimplicesNode(m_simplex_mesh, A).getCoords();
    const math::Point ptB = SimplicesNode(m_simplex_mesh, B).getCoords();
    const math::Point ptC = SimplicesNode(m_simplex_mesh, C).getCoords();

    const math::Vector3d v0 = ptA - pt;
    const math::Vector3d v1 = ptB - pt;
    const math::Vector3d v2 = ptC - pt;
    double num = std::fabs((v0.dot(v1.cross(v2))));
    double den = v0.norm()*v1.norm()*v2.norm();
    ans.push_back(num / den);
  }
  return ans;
}
/*----------------------------------------------------------------------------*/
std::vector<TInt> SimplicesCell::commonNode(SimplicesCell cell1)
{
  std::vector<TInt> ans{};
  std::vector<TInt> nodes0 = getNodes();
  std::vector<TInt> nodes1 = cell1.getNodes();
  std::vector<TInt>::iterator it;
  std::vector<TInt>::iterator end;
  std::vector<TInt> inter(nodes0.size() + nodes1.size());

  end = std::set_intersection(
        nodes0.begin(), nodes0.end(),
        nodes1.begin(), nodes1.end(),
        inter.begin());

  std::copy(inter.begin(), end, std::back_inserter(ans));
  return ans;

}
/*----------------------------------------------------------------------------*/
