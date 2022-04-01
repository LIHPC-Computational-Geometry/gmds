/******************************************************************************/
#include<chrono>
/******************************************************************************/
#include<gmds/hybridMeshAdapt/SimplicesNode.h>
/******************************************************************************/
#include<gmds/hybridMeshAdapt/SimplexMesh.h> //pour la forward declaration
/******************************************************************************/
using namespace gmds;
using namespace hybrid;
using namespace simplicesNode;
using namespace simplicesCell;
using namespace simplicesTriangle;
using namespace operators;
using namespace math;
/******************************************************************************/
SimplicesNode::SimplicesNode(SimplexMesh* simplexMesh, const TInt indexPoint)
{
  m_simplex_mesh = simplexMesh;
  m_indexPoint   = indexPoint;
  if(m_simplex_mesh->m_node_ids[m_indexPoint] != 1)
  {
    /*TODO exeption le node d'existe pas ... le creer avant de poouvoir l'utiliser*/
    std::cout << "Creer le node " << m_indexPoint <<  " avant de l'utiliser !!" << std::endl;
    GMDSException("ERREUR NODE INEXISTANTE");
  }
}
/******************************************************************************/
SimplicesNode::~SimplicesNode()
{

}
/******************************************************************************/
bool SimplicesNode::checkExistance()
{
  bool flag = false;
  if(m_simplex_mesh->m_node_ids[m_indexPoint] == 1)
  {
    flag = true;
  }
  return flag;
}
/******************************************************************************/
std::vector<TSimplexID> SimplicesNode::ballOf(bool boundariesAccepted) const
{
  //std::cout << "ballOf(bool boundariesAccepted)" << std::endl;
  //std::cout << "m_indexPoint -> " << m_indexPoint << std::endl;
  std::vector<TSimplexID> v{};
  std::vector<TSimplexID> to_do{};
  int border = std::numeric_limits<int>::min();
  bool addBorder = false;

  if(m_simplex_mesh->m_base.size() != 0)
  {
    //permet de voir si base ne contient pas un tetra deja supprimé (dans le futu penser a ameliorer build base vector pour evité de faire ce genre de chose....)
    if(m_simplex_mesh->m_base[m_indexPoint] != border)
    {
      if(m_simplex_mesh->m_tet_ids[m_simplex_mesh->m_base[m_indexPoint]] != 0)
      {
        BitVector flagTet = m_simplex_mesh->m_tet_ids;
        BitVector flagTri = m_simplex_mesh->m_tri_ids;

        to_do.push_back(m_simplex_mesh->m_base[m_indexPoint]);
        //std::cout << "m_simplex_mesh->m_base[m_indexPoint] -> " << m_simplex_mesh->m_base[m_indexPoint] << std::endl;
        //std::cout << "flagTet.capacity() -> " << flagTet.capacity() << std::endl;
        //flagTet.unselect(m_simplex_mesh->m_base[m_indexPoint]);
        while(!to_do.empty())
        {
          TSimplexID tetraId =  to_do.back();
          //std::cout << "tetraId --> " << tetraId << std::endl;
          if(tetraId != border)
          {
            /*contain tet(>= 0) and triangle( < 0)*/
            std::vector<TSimplexID> neighborTetras;
            if(tetraId >=0)
            {
              if(flagTet[tetraId] == 0){to_do.pop_back(); continue;}
              neighborTetras = SimplicesCell(m_simplex_mesh, tetraId).neighborTetra(m_indexPoint, boundariesAccepted);
            }
            else
            {
              if(flagTri[-tetraId] == 0){to_do.pop_back(); continue;}
              neighborTetras = SimplicesTriangle(m_simplex_mesh, -tetraId).neighborTetra();
            }
            //std::cout << "tetraId -> "<< tetraId << std::endl;
            to_do.pop_back();

            if(tetraId >= 0) /*tetraId is a Tet*/
            {
              if(m_simplex_mesh->m_tet_ids[tetraId] != 0)
              {
                if(flagTet[tetraId] != 0)
                {
                  v.push_back(tetraId);
                  flagTet.unselect(tetraId);
                }
              }
            }
            else /*tetraId is a tri*/
            {
              if(m_simplex_mesh->m_tri_ids[-tetraId] != 0)
              {
                if(flagTri[-tetraId] != 0)
                {
                  v.push_back(tetraId);
                  flagTri.unselect(-tetraId);
                }
              }
            }


            for(auto const & val : neighborTetras)
            {
              if(val != border)
              {
                if(val >= 0)
                {
                  if(flagTet[val] != 0)
                  {
                    to_do.push_back(val);
                  }
                }
                else
                {
                  if(flagTri[-val] != 0)
                  {
                    to_do.push_back(val);
                  }
                }
              }
              else
              {
                to_do.push_back(val);
              }
            }
          }
          else if(tetraId == border && boundariesAccepted == true)
          {
            addBorder = true;
            to_do.pop_back();
          }
          else
          {
            to_do.pop_back();
          }
        }

        if(addBorder == true)
        {
            v.push_back(border); //in order to not use find when looking for border just use the last item ...
        }

      }
    }
  }

  return std::move(v);

}
/********************************************************************************/
std::vector<TSimplexID> SimplicesNode::linksTri() const
{
  std::vector<TSimplexID> v{};
  std::vector<TSimplexID> to_do{};
  if(m_simplex_mesh->m_base.size() != 0)
  {
    std::vector<TSimplexID>&& neighborTetras = ballOf();
    for(auto const& tetIdx : neighborTetras)
    {
      std::vector<TSimplexID>&& v_tmp = SimplicesCell(m_simplex_mesh, tetIdx).neighborTri(this->getGlobalNode());
      v.resize(v_tmp.size() + v.size());
      std::move(std::begin(v_tmp), std::end(v_tmp), std::back_inserter(v));
    }

  }

    // A OPTIMISER!!!!!! (on a des duplication de simplex dans le tableau v ... donc on retire ceux qui existe deja dans ce tableau..)
    std::vector<TSimplexID> v_tmp = v;

    for(auto& currentVal : v)
    {
      int iter = 0;
      for(auto& val : v_tmp)
      {
        if(val == currentVal)
        {
          if(iter != 0)
          {
            val = std::numeric_limits<int>::min();
          }
          iter++;
        }
      }
    }

    std::vector<TSimplexID> res;
    for(auto const& val_tmp : v_tmp)
    {
      if(val_tmp != std::numeric_limits<int>::min() && val_tmp)
      {
        res.push_back(val_tmp);
      }
    }
    return std::move(res);

}
/******************************************************************************/
std::vector<TSimplexID> SimplicesNode::shell(const SimplicesNode& simplicesNode) const
{
  std::vector<TSimplexID> v{};
  std::vector<TSimplexID>&& ballOf_vector0 = ballOf();
  std::vector<TSimplexID>&& ballOf_vector1 = simplicesNode.ballOf();
  std::sort(ballOf_vector0.begin(), ballOf_vector0.end());
  std::sort(ballOf_vector1.begin(), ballOf_vector1.end());
  std::set_intersection(ballOf_vector0.begin(), ballOf_vector0.end(), ballOf_vector1.begin(), ballOf_vector1.end(), std::back_inserter(v));

  return std::move(v);
}
/******************************************************************************/
Point SimplicesNode::getCoords() const
{
  Point p;
  if(m_simplex_mesh->m_node_ids[m_indexPoint] == true)
  {
    p = m_simplex_mesh->m_coords[m_indexPoint];
  }
  return std::move(p);
}
/******************************************************************************/
math::Vector3d SimplicesNode::getNormal() const
{
  unsigned int border = std::numeric_limits<TInt>::min();
  std::vector<TSimplexID> ball = ballOf();
  math::Vector3d interpolateNormal{};
  unsigned int triangleNbr = 0;
  for(auto const simplex : ball)
  {
    if(simplex < 0 && simplex != border)
    {
      triangleNbr++;
      math::Vector3d n = SimplicesTriangle(m_simplex_mesh, -simplex).getNormal();
      n.normalize();
      interpolateNormal = interpolateNormal + n;
    }
  }

  if(triangleNbr)
  {
      interpolateNormal = interpolateNormal / triangleNbr;
  }
  else
  {
    std::cout << "triangleNbr == 0" << std::endl;
  }

  return interpolateNormal;
}
/******************************************************************************/
bool SimplicesNode::isInBorder() const
{
  bool boundariesAccepted = true;
  bool flag               = false;
  bool boundary           = std::numeric_limits<int>::min();
  std::vector<TSimplexID>&& ball = ballOf(boundariesAccepted);

  if(std::find(ball.begin(), ball.end(), boundary) != ball.end())
  {
    flag = true;
  }

  return flag;
}
/******************************************************************************/
std::vector<TSimplexID>  SimplicesNode::shellTriangleOrdererdWithHole(const SimplicesNode& simplicesNode) const
{
  std::vector<TSimplexID> ordererdTriangles{};
  std::vector<TSimplexID>&& shellTet = this->shell(simplicesNode);

  if(shellTet.size() > 0)
  {
    unsigned int     nbrTetraShell = shellTet.size();
    std::vector<TSimplexID>&& triangleShell = this->triangleUnordererdShell(simplicesNode);
    unsigned int nbrTriangleShell = triangleShell.size();
    unsigned int nbrSimplex = nbrTetraShell + nbrTriangleShell;
    int errorId = std::numeric_limits<int>::min();
    bool hole = false;
    TSimplexID firstTet = shellTet[0];
    TSimplexID currentSimplex = firstTet; // can be positive or neg
    std::vector<TSimplexID> markedSimplex{firstTet};

    std::vector<TInt> nodes{this->getGlobalNode(), simplicesNode.getGlobalNode()};
    for(;;)
    {
      std::vector<TSimplexID> nodesExt;
      if(currentSimplex >= 0)
      {
        if(m_simplex_mesh->m_tet_ids[currentSimplex] != 0)
        {
          nodesExt = SimplicesCell(m_simplex_mesh, currentSimplex).getOtherNodeInSimplex(nodes);
          if(nodesExt.size() != 2)
          {
            //TODO ecxeption
          }
          else
          {
            TSimplexID simplex = m_simplex_mesh->m_tet_adj[currentSimplex][SimplicesCell(m_simplex_mesh, currentSimplex).getLocalNode(nodesExt[0])];
            if(std::find(markedSimplex.begin(), markedSimplex.end(), simplex) != markedSimplex.end())
            {
              simplex = m_simplex_mesh->m_tet_adj[currentSimplex][SimplicesCell(m_simplex_mesh, currentSimplex).getLocalNode(nodesExt[1])];
            }
            currentSimplex = simplex;
            if(currentSimplex != errorId && std::find(markedSimplex.begin(), markedSimplex.end(), currentSimplex) == markedSimplex.end())
            {
              markedSimplex.push_back(currentSimplex);
            }
          }
        }
        else
        {
          //TODO exception
        }
      }
      else if(currentSimplex < 0)
      {
        if(currentSimplex == errorId)
        {
          hole = true;
        }
        else
        {
          TSimplexID simplex = m_simplex_mesh->m_tri_nodes[-currentSimplex][3];
          if(std::find(markedSimplex.begin(), markedSimplex.end(), simplex) != markedSimplex.end())
          {
            simplex = m_simplex_mesh->m_tri_adj[-currentSimplex][3];
          }
          currentSimplex = simplex;
          if(currentSimplex != errorId && std::find(markedSimplex.begin(), markedSimplex.end(), currentSimplex) == markedSimplex.end())
          {
            markedSimplex.push_back(currentSimplex);
          }
        }
      }
      if((markedSimplex.size() >= nbrSimplex) || hole == true)
      {
        break;
      }
    }

    if(hole)
    {
      TSimplexID firstSimplex = markedSimplex.back();
      markedSimplex.clear();
      currentSimplex = firstSimplex;

      for(unsigned int simplexIter = 0; simplexIter < nbrSimplex; simplexIter++)
      {
        if(currentSimplex >=0)
        {
          std::vector<TSimplexID>&& nodesExt = SimplicesCell(m_simplex_mesh, currentSimplex).getOtherNodeInSimplex(nodes);
          if(std::find(markedSimplex.begin(), markedSimplex.end(), currentSimplex) == markedSimplex.end())
          {
            markedSimplex.push_back(currentSimplex);
          }
          TSimplexID simplex = m_simplex_mesh->m_tet_adj[currentSimplex][SimplicesCell(m_simplex_mesh, currentSimplex).getLocalNode(nodesExt[0])];
          if(simplex == errorId || std::find(markedSimplex.begin(), markedSimplex.end(), simplex) != markedSimplex.end())
          {
            simplex = m_simplex_mesh->m_tet_adj[currentSimplex][SimplicesCell(m_simplex_mesh, currentSimplex).getLocalNode(nodesExt[1])];
          }
          currentSimplex = simplex;
        }
        else
        {
          if(std::find(markedSimplex.begin(), markedSimplex.end(), currentSimplex) == markedSimplex.end())
          {
            markedSimplex.push_back(currentSimplex);
          }
          TSimplexID simplex = m_simplex_mesh->m_tri_nodes[-currentSimplex][3];
          if(simplex == errorId || std::find(markedSimplex.begin(), markedSimplex.end(), simplex) != markedSimplex.end())
          {
            simplex = m_simplex_mesh->m_tri_adj[-currentSimplex][3];
          }
          currentSimplex = simplex;
        }
      }
      //le dernier simplex marker dans cette situation est le simplex qui voit le bord;
    }

    for(auto const & simplexMarked : markedSimplex)
    {
      if(simplexMarked < 0  && simplexMarked != errorId && std::find(ordererdTriangles.begin(), ordererdTriangles.end(), simplexMarked) == ordererdTriangles.end())
      {
        ordererdTriangles.push_back(simplexMarked);
      }
    }
  }


  return std::move(ordererdTriangles);
}
/******************************************************************************/
std::vector<hybrid::TSimplexID> SimplicesNode::triangleUnordererdShell(const SimplicesNode& simplicesNode) const
{
  std::vector<TSimplexID>&& linkTriA = this->linksTri();
  std::vector<TSimplexID>&& linkTriB = simplicesNode.linksTri();
  std::vector<TSimplexID> v{};

  for(auto const & triangleA : linkTriA)
  {
    if(std::find(linkTriB.begin(), linkTriB.end() ,triangleA) != linkTriB.end())
    {
      v.push_back(triangleA);
    }
  }

  return std::move(v);
}
/******************************************************************************/
std::vector<TInt>  SimplicesNode::commonNodeInBall(const std::vector<TInt>& nodeVec) const
{
  std::vector<TInt> v{};
  std::set<TInt> setV{};
  std::vector<TSimplexID>&& ball = ballOf();
  std::vector<std::vector<TInt>> nodeOfTet{};
  for(auto const & tet : ball)
  {
    std::vector<TInt> && nodes = SimplicesCell(m_simplex_mesh, tet).getNodes();
    nodeOfTet.push_back(nodes);
  }
  m_simplex_mesh->intersection(nodeOfTet, nodeVec);
  for(auto const & nodes : nodeOfTet)
  {
    for(auto const & node : nodes)
    {
      setV.insert(node);
    }
  }

  std::move(setV.begin(), setV.end(), std::back_inserter(v));

  return v;
}
/******************************************************************************/
std::vector<TSimplexID> SimplicesNode::directSimplices(const math::Vector3d& vector) const
{
  double epsilon                    = 0.1;
  bool simplexFound                 = false;
  int errorId                       = std::numeric_limits<int>::min();
  std::vector<TSimplexID> && ball   = ballOf();
  std::vector<TSimplexID> dirSimplex{};
  math::Vector3d dir = vector;
  dir.normalize();

  while(simplexFound != true)
  {
    math::Point pt = getCoords() + epsilon * dir;
    std::vector<TSimplexID> tetraContainingPt{};
    nodeNeighborInfo AnodeInfo;
    if(m_simplex_mesh->checkSimplicesContenaing(pt, tetraContainingPt))
    {
      dirSimplex = tetraContainingPt;
      simplexFound = true;
    }
    epsilon *= 0.1;
    if(epsilon < 10E-20)
    {
      break;
    }
  }

  return std::move(dirSimplex);

}
/******************************************************************************/
TSimplexID SimplicesNode::directSimplex(const math::Vector3d& vector) const
{
  double epsilon                    = 1.0;
  bool simplexFound                 = false;
  int errorId                       = std::numeric_limits<int>::min();
  std::vector<TSimplexID> && ball   = ballOf();
  TSimplexID directSimplexCrossingVector = errorId;
  math::Vector3d dir = vector;
  dir.normalize();

  while(simplexFound != true)
  {
    math::Point pt = getCoords() + epsilon * dir;
    for(auto const & simplex : ball)
    {
      if(SimplicesCell(m_simplex_mesh, simplex).isPointInSimplex(pt))
      {
        directSimplexCrossingVector = simplex;
        simplexFound = true;
        break;
      }
    }
    epsilon *= 0.1;
    if(epsilon < 10E-20)
    {
      break;
    }
  }

  return directSimplexCrossingVector;

}
/******************************************************************************/
std::vector<TSimplexID> SimplicesNode::simplicesIntersectedbyRay(const SimplicesNode& simpliceNodeB) const
{
  std::vector<TSimplexID> v;
  math::Vector3d normalAB = simpliceNodeB.getCoords() - getCoords();
  normalAB.normalize();
  TSimplexID simplexDirect = directSimplex(normalAB);
  int errorId = std::numeric_limits<int>::min();

  if(simplexDirect != errorId)
  {
    v.push_back(simplexDirect);
    TSimplexID simplexAdj;// = SimplicesCell(m_simplex_mesh, simplexDirect).oppositeTetraIdx(*this);
    //todo
  }

  return std::move(v);
}
/******************************************************************************/
std::vector<TSimplexID> SimplicesNode::simplicesIntersectedbyRayAndCheckVisibility(std::vector<std::vector<TInt>>& facesId, const SimplicesNode& simpliceNodeB) const
{
  std::vector<TSimplexID> v;
  math::Vector3d normalAB = simpliceNodeB.getCoords() - getCoords();
  normalAB.normalize();
  TSimplexID simplexDirect = directSimplex(normalAB);
  int errorId = std::numeric_limits<int>::min();

  if(simplexDirect != errorId)
  {
    v.push_back(simplexDirect);
    TSimplexID simplexAdj;// = SimplicesCell(m_simplex_mesh, simplexDirect).oppositeTetraIdx(*this);
    if(simplexAdj != errorId)
    {
      v.push_back(simplexAdj); // par copie


      for(;;)
      {
        unsigned int nbrFaceNotVisible = SimplicesCell(m_simplex_mesh, simplexAdj).checkFaceNbrVisibility(facesId, simpliceNodeB);

        if(facesId.size() == nbrFaceNotVisible)
        {
          if(nbrFaceNotVisible > 1)
          {
            break;
          }
          else if(nbrFaceNotVisible == 1)
          {
            std::vector<TInt> face          = facesId[0];
            TSimplexID simplexTMP           = simplexAdj;
            std::vector<TSimplexID> && ball = SimplicesNode(m_simplex_mesh, face[0]).ballOf();
            for(auto const & tet : ball)
            {
              if(tet != simplexAdj && SimplicesCell(m_simplex_mesh, tet).containNodes(face))
              {
                v.push_back(tet);
                simplexAdj = tet;
                break;
              }
            }

          }
          else
          {
            //std::cout << "nbrFaceNotVisible < 1 || nbrFaceNotVisible != 1" << std::endl;
            break;
          }

          if(SimplicesCell(m_simplex_mesh, simplexAdj).containNode(simpliceNodeB))
          {
            break;
          }

        }
      }
    }
    else
    {
      std::cout << "simplexAdj != errorId" << std::endl;
    }
  }

  return std::move(v);
}
/******************************************************************************/
bool SimplicesNode::checkStar(const std::vector<TSimplexID>& cavity, const CriterionRAIS& criterion) const
{
  int border = std::numeric_limits<int>::min();
  struct tetAndNode{
    TSimplexID tet;
    std::vector<TInt> nodesLocal;
  };
  std::vector<tetAndNode> borderSimplex{};
  std::vector<TSimplexID> trianglesConnectedToP{};
  std::vector<TSimplexID> trianglesNotConnectedToP{};
  CavityOperator::CavityIO cavIO = CavityOperator::CavityIO(m_simplex_mesh, cavity, trianglesConnectedToP, trianglesNotConnectedToP);

  for(auto const & tet : cavity)
  {
    std::vector<TInt> nodesLocal{};
    if(cavIO.simplexInborderOfCavity(tet, nodesLocal))
    {
      tetAndNode tmpData;
      tmpData.tet = tet;
      tmpData.nodesLocal = nodesLocal;
      borderSimplex.push_back(tmpData);
    }
  }

  for(auto const& dataBorderTet : borderSimplex)
  {
    if(dataBorderTet.tet != border)
    {
      for(auto const& nodeLocal : dataBorderTet.nodesLocal)
      {
        const SimplicesCell cell = SimplicesCell(m_simplex_mesh, dataBorderTet.tet);
        const TInt nodeGlobal = cell.getNode(nodeLocal).getGlobalNode();
        const std::vector<TInt> nodeVec{nodeGlobal};
        const std::vector<TInt>&& faceNodes = cell.getOtherNodeInSimplex(nodeVec);
        if(std::find(faceNodes.begin(), faceNodes.end(), m_indexPoint) != faceNodes.end())
        {
            //break;
            continue;
        }

        //bool flag = criterion.execute(m_simplex_mesh, dataBorderTet.tet, nodeLocal, getCoords());
        bool flag = false;
        math::Orientation::Sign orientation = SimplicesCell(m_simplex_mesh, dataBorderTet.tet).orientation(nodeLocal, getCoords());
        flag = (orientation < 0)? true:false;
        if(flag == true)
        {
          return false;
        }
      }
    }
  }
  return true;
}
/******************************************************************************/
bool SimplicesNode::isFaceVisible(math::Vector3d& normalOfFace, math::Vector3d& vecFacePt) const
{
  bool flag = true;
  vecFacePt.normalize();
  normalOfFace.normalize();
  double epsilonTheta = M_PI / 2.0  - M_PI* 3.0 / 180.0;
  double epsilonLenght = cos(epsilonTheta);
  if(vecFacePt.dot(normalOfFace) >= epsilonLenght)
  {

  }
  else
  {
    flag = false;
  }

  return flag;
}
/******************************************************************************/
bool SimplicesNode::isAttachToSimplex() const
{
  TSimplexID border  = std::numeric_limits<int>::min();
  TSimplexID simplex = m_simplex_mesh->m_base[m_indexPoint];
  return (simplex != border);
}
/******************************************************************************/
std::vector<TInt> SimplicesNode::complentaryNodeShell(const SimplicesNode& simpliceNode) const
{
  std::set<TInt> s{};
  std::vector<TInt> v{};
  std::vector<TInt>nodeShell{this->getGlobalNode(), simpliceNode.getGlobalNode()};
  std::vector<TSimplexID> shell = this->shell(simpliceNode);
  for(auto const simplice : shell)
  {
    if(simplice >= 0)
    {
      std::vector<TInt> v = SimplicesCell(m_simplex_mesh, simplice).getOtherNodeInSimplex(nodeShell);
      copy(v.begin(), v.end(),inserter(s, s.begin()));
    }
  }
  std::copy(s.begin(), s.end(), std::back_inserter(v));
  return v;
}
/******************************************************************************/
std::vector<TInt> SimplicesNode::directNeighboorNodeId() const
{
  std::vector<TInt> res{};
  std::set<TInt> s{};
  const std::vector<TInt>&& ball = ballOf();

  for(auto const simplex : ball)
  {
    if(simplex > 0)
    {
      SimplicesCell cell(m_simplex_mesh, simplex);
      const std::vector<TInt>&& directNodes = cell.getNodes();
      std::copy(directNodes.begin(), directNodes.end(), std::inserter(s, s.begin()));
    }
  }

  std::copy_if(s.begin(), s.end(), std::back_inserter(res), [&](const TInt currentNode){
    return (currentNode != m_indexPoint);
  });

  return res;
}
/******************************************************************************/
void SimplicesNode::detectType(const nodeNeighborInfo& nodeInfo) const
{
  TInt border = std::numeric_limits<int>::min();
  if(ballOf(true).back() == border)
  {
    Variable<int>* BND_VERTEX_COLOR = m_simplex_mesh->getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
    Variable<int>* BND_CURVE_COLOR = m_simplex_mesh->getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
    Variable<int>* BND_SURFACE_COLOR = m_simplex_mesh->getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
    if(BND_VERTEX_COLOR != nullptr && BND_CURVE_COLOR != nullptr && BND_SURFACE_COLOR != nullptr)
    {
      if(nodeInfo.nodes.size() == 1)
      {
        TInt node = nodeInfo.nodes.front();
        if((*BND_VERTEX_COLOR)[node] != 0)
        {
          (* BND_VERTEX_COLOR)[m_indexPoint] = (*BND_VERTEX_COLOR)[node];
        }
        else if((*BND_CURVE_COLOR)[node] != 0)
        {
          (*BND_CURVE_COLOR)[m_indexPoint] = (*BND_CURVE_COLOR)[node];
        }
        else if((*BND_SURFACE_COLOR)[node] != 0)
        {
          (*BND_SURFACE_COLOR)[m_indexPoint] = (*BND_SURFACE_COLOR)[node];
        }
      }
      else if(nodeInfo.nodes.size() == 2)
      {
        TInt nodeA = nodeInfo.nodes.front();
        TInt nodeB = nodeInfo.nodes.back();

        struct proprNode{
            int isOn;
            int number;
        };
        //node on vertexBoundary -> 0 / node on curveBoundary -> 1 / node on surfaceBoundary -> 2 / node on cell ->3
        proprNode proprnodeA;
        proprNode proprnodeB;

        proprnodeA.isOn = ((* BND_VERTEX_COLOR)[nodeA] != 0)? 0 : ((* BND_CURVE_COLOR)[nodeA] != 0 )? 1 : ((* BND_SURFACE_COLOR)[nodeA] != 0)? 2 : 3;
        proprnodeB.isOn = ((* BND_VERTEX_COLOR)[nodeB] != 0)? 0 : ((* BND_CURVE_COLOR)[nodeB] != 0 )? 1 : ((* BND_SURFACE_COLOR)[nodeB] != 0)? 2 : 3;

        if(proprnodeA.isOn == 0){proprnodeA.number = (* BND_VERTEX_COLOR)[nodeA];}
        if(proprnodeA.isOn == 1){proprnodeA.number = (* BND_CURVE_COLOR)[nodeA];}
        if(proprnodeA.isOn == 2){proprnodeA.number = (* BND_SURFACE_COLOR)[nodeA];}

        if(proprnodeB.isOn == 0){proprnodeB.number = (* BND_VERTEX_COLOR)[nodeB];}
        if(proprnodeB.isOn == 1){proprnodeB.number = (* BND_CURVE_COLOR)[nodeB];}
        if(proprnodeB.isOn == 2){proprnodeB.number = (* BND_SURFACE_COLOR)[nodeB];}

        if(proprnodeA.isOn == proprnodeB.isOn)
        {
          if(proprnodeA.isOn == 0)      {(* BND_VERTEX_COLOR)[m_indexPoint] = proprnodeA.number;}
          else if(proprnodeA.isOn == 1) {(* BND_CURVE_COLOR)[m_indexPoint] = proprnodeA.number;}
          else if(proprnodeA.isOn == 2) {(* BND_SURFACE_COLOR)[m_indexPoint] = proprnodeA.number;}
        }
        else if(proprnodeA.isOn < proprnodeB.isOn)
        {
          if(proprnodeB.isOn == 0)      {(* BND_VERTEX_COLOR)[m_indexPoint] = proprnodeB.number;}
          else if(proprnodeB.isOn == 1) {(* BND_CURVE_COLOR)[m_indexPoint] = proprnodeB.number;}
          else if(proprnodeB.isOn == 2) {(* BND_SURFACE_COLOR)[m_indexPoint] = proprnodeB.number;}
        }
        else
        {
          if(proprnodeA.isOn == 0)      {(* BND_VERTEX_COLOR)[m_indexPoint] = proprnodeA.number;}
          else if(proprnodeA.isOn == 1) {(* BND_CURVE_COLOR)[m_indexPoint] = proprnodeA.number;}
          else if(proprnodeA.isOn == 2) {(* BND_SURFACE_COLOR)[m_indexPoint] = proprnodeA.number;}
        }
      }
      else if(nodeInfo.nodes.size() == 3)
      {
        TInt nodeA = nodeInfo.nodes.front();
        TInt nodeB = nodeInfo.nodes[1];
        TInt nodeC = nodeInfo.nodes.back();



        struct proprNode{
            int isOn;
            int number;
        };
        //node on vertexBoundary -> 0 / node on curveBoundary -> 1 / node on surfaceBoundary -> 2 / node on cell ->3
        proprNode proprnodeA;
        proprNode proprnodeB;
        proprNode proprnodeC;

        proprnodeA.isOn = ((* BND_VERTEX_COLOR)[nodeA] != 0)? 0 : ((* BND_CURVE_COLOR)[nodeA] != 0 )? 1 : ((* BND_SURFACE_COLOR)[nodeA] != 0)? 2 : 3;
        proprnodeB.isOn = ((* BND_VERTEX_COLOR)[nodeB] != 0)? 0 : ((* BND_CURVE_COLOR)[nodeB] != 0 )? 1 : ((* BND_SURFACE_COLOR)[nodeB] != 0)? 2 : 3;
        proprnodeC.isOn = ((* BND_VERTEX_COLOR)[nodeC] != 0)? 0 : ((* BND_CURVE_COLOR)[nodeC] != 0 )? 1 : ((* BND_SURFACE_COLOR)[nodeC] != 0)? 2 : 3;

        if(proprnodeA.isOn == 0){proprnodeA.number = (* BND_VERTEX_COLOR)[nodeA];}
        if(proprnodeA.isOn == 1){proprnodeA.number = (* BND_CURVE_COLOR)[nodeA];}
        if(proprnodeA.isOn == 2){proprnodeA.number = (* BND_SURFACE_COLOR)[nodeA];}

        if(proprnodeB.isOn == 0){proprnodeB.number = (* BND_VERTEX_COLOR)[nodeB];}
        if(proprnodeB.isOn == 1){proprnodeB.number = (* BND_CURVE_COLOR)[nodeB];}
        if(proprnodeB.isOn == 2){proprnodeB.number = (* BND_SURFACE_COLOR)[nodeB];}

        if(proprnodeC.isOn == 0){proprnodeC.number = (* BND_VERTEX_COLOR)[nodeC];}
        if(proprnodeC.isOn == 1){proprnodeC.number = (* BND_CURVE_COLOR)[nodeC];}
        if(proprnodeC.isOn == 2){proprnodeC.number = (* BND_SURFACE_COLOR)[nodeC];}

        if(proprnodeA.isOn == proprnodeB.isOn && proprnodeA.isOn == proprnodeC.isOn)
        {
          if(proprnodeA.isOn == 0)      {(* BND_VERTEX_COLOR)[m_indexPoint] = proprnodeA.number;}
          else if(proprnodeA.isOn == 1) {(* BND_CURVE_COLOR)[m_indexPoint] = proprnodeA.number;}
          else if(proprnodeA.isOn == 2) {(* BND_SURFACE_COLOR)[m_indexPoint] = proprnodeA.number;}
        }
        else
        {

          if(proprnodeB.isOn < proprnodeC.isOn)
          {
            proprnodeB.isOn = proprnodeC.isOn;
            proprnodeB.number = proprnodeC.number;
          }

          if(proprnodeA.isOn < proprnodeB.isOn)
          {
            if(proprnodeB.isOn == 0)      {(* BND_VERTEX_COLOR)[m_indexPoint] = proprnodeB.number;}
            else if(proprnodeB.isOn == 1) {(* BND_CURVE_COLOR)[m_indexPoint] = proprnodeB.number;}
            else if(proprnodeB.isOn == 2) {(* BND_SURFACE_COLOR)[m_indexPoint] = proprnodeB.number;}
          }
          else
          {
            if(proprnodeA.isOn == 0)      {(* BND_VERTEX_COLOR)[m_indexPoint] = proprnodeA.number;}
            else if(proprnodeA.isOn == 1) {(* BND_CURVE_COLOR)[m_indexPoint] = proprnodeA.number;}
            else if(proprnodeA.isOn == 2) {(* BND_SURFACE_COLOR)[m_indexPoint] = proprnodeA.number;}
          }
        }
      }
    }
  }
}
/******************************************************************************/
