#include<gmds/sofiane/SimplicesNode.h>
/******************************************************************************/
#include<gmds/sofiane/SimplexMesh.h> //pour la forward declaration
/******************************************************************************/
using namespace gmds;
using namespace hybrid;
using namespace simplicesNode;
using namespace simplicesCell;
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
  }
}
/******************************************************************************/
SimplicesNode::~SimplicesNode()
{

}
/******************************************************************************/
std::vector<TSimplexID> SimplicesNode::ballOf(bool boundariesAccepted) const
{
  std::vector<TSimplexID> v{};
  std::vector<TSimplexID> to_do{};
  int errorId = std::numeric_limits<int>::min();
  bool addBorder = false;

  if(m_simplex_mesh->m_base.size() != 0)
  {
    //permet de voir si base ne contient pas un tetra deja supprimé (dans le futu penser a ameliorer build base vector pour evité de faire ce genre de chose....)
    if(m_simplex_mesh->m_base[m_indexPoint] != errorId)
    {
      if(m_simplex_mesh->m_tet_ids[m_simplex_mesh->m_base[m_indexPoint]] != 0)
      {
        to_do.push_back(m_simplex_mesh->m_base[m_indexPoint]);
        while(!to_do.empty())
        {
          TSimplexID tetraId =  to_do.back();
          /*affectation semantic (&&) et operateur de copie semantique (dans neighborTetra il y'a std::move(v)) */
          if(tetraId != errorId )
          {
            SimplicesCell(m_simplex_mesh, tetraId);
            std::vector<TSimplexID>&& neighborTetras = SimplicesCell(m_simplex_mesh, tetraId).neighborTetra(m_indexPoint, boundariesAccepted);

            to_do.pop_back();
            if(std::find(v.begin(), v.end(), tetraId) == v.end())
            {
              v.push_back(tetraId);
            }

            for(auto const & val : neighborTetras)
            {
              if(find(v.begin(), v.end(), val) == v.end())
              {
                to_do.push_back(val);
              }
            }
          }
          else if(tetraId == errorId && boundariesAccepted == true)
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
            v.push_back(errorId); //in order to not use find when looking for border just use the last item ...
        }

      }
    }

  }
  return std::move(v);
}
/******************************************************************************/
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
  const std::vector<TSimplexID>&& ballOf_vector = ballOf();
  std::vector<TSimplexID> v;

  for(const auto & val : ballOf_vector)
  {
    if(SimplicesCell(m_simplex_mesh, val).containNode(simplicesNode))
    {
      v.push_back(val);
    }
  }

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
std::vector<hybrid::TSimplexID> SimplicesNode::triangleUnordererdShell (const SimplicesNode& simplicesNode) const
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
