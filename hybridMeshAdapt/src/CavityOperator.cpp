/******************************************************************/
#include <gmds/hybridMeshAdapt/CavityOperator.h>
#include <gmds/hybridMeshAdapt/ICriterion.h>
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
}
/******************************************************************/
CavityOperator::CavityIO::CavityIO(SimplexMesh* simplexMesh, std::vector<TSimplexID>&& cavityIn):
m_simplex_mesh(simplexMesh)
{
  m_cavityIn = cavityIn;
  cavityIn.clear();
}
/******************************************************************/
CavityOperator::CavityIO::CavityIO(CavityIO&& cavityIO)
{
  m_simplex_mesh = cavityIO.m_simplex_mesh;
  m_cavityIn  = cavityIO.m_cavityIn;

  cavityIO.m_cavityIn.clear();
  cavityIO.m_simplex_mesh = nullptr;
}
/******************************************************************/
CavityOperator::CavityIO& CavityOperator::CavityIO::CavityIO::operator=(CavityIO&& cavityIO)
{
  if(this != &cavityIO)
  {
    m_cavityIn.clear();
    /*A VOIR SI JAMAIS IL YA UN PROBLEME LORS D'UNE COPIE*/
    m_simplex_mesh = nullptr;

    m_cavityIn  = cavityIO.m_cavityIn;
    m_simplex_mesh = cavityIO.m_simplex_mesh;

    cavityIO.m_cavityIn.clear();
    cavityIO.m_simplex_mesh = nullptr;
  }
  return *this;
}
/******************************************************************/
std::vector<TInt>& CavityOperator::CavityIO::CavityIO::nodeInCavity()
{
  m_nodeInCavity.clear();
  unsigned int sizeindexLocalTet = 4;

  for(const auto & simplexIn : simplexInCavityUnmodifiable())
  {
    for(unsigned int nodeLocal = 0; nodeLocal < sizeindexLocalTet; nodeLocal++)
    {
      bool flag = true;
      SimplicesNode simpliceNode  = SimplicesCell(m_simplex_mesh, simplexIn).getNode(nodeLocal);
      if(std::find(m_nodeInCavity.begin(), m_nodeInCavity.end(), simpliceNode.getGlobalNode()) == m_nodeInCavity.end())
      {
        bool boundariesAccepted = true;
        std::vector<TSimplexID>&& ballOfNode = simpliceNode.ballOf(boundariesAccepted);

        for(auto const & simplexInBall : ballOfNode)
        {
            if(std::find(simplexInCavityUnmodifiable().begin(), simplexInCavityUnmodifiable().end(), simplexInBall) == simplexInCavityUnmodifiable().end())
            {
              flag = false;
              break;
            }
        }

        if(flag)
        {
          m_nodeInCavity.push_back(simpliceNode.getGlobalNode());
        }
      }
    }
  }

  return m_nodeInCavity;
}
/******************************************************************//*
std::vector<std::vector<TInt>> CavityOperator::CavityIO::CavityIO::pointToConnect()
{
  std::vector<std::vector<TInt>> v{};
  const std::vector<TInt>& nodesInCavity = getNodeInCavity();

  for(const auto & simplexInCavity : simplexInCavityUnmodifiable())
  {
    std::cout << SimplicesCell(m_simplex_mesh, simplexInCavity) << std::endl;
    std::vector<TInt> pointToConnect =
    {
      SimplicesCell(m_simplex_mesh,simplexInCavity).getNode(0).getGlobalNode(),
      SimplicesCell(m_simplex_mesh,simplexInCavity).getNode(1).getGlobalNode(),
      SimplicesCell(m_simplex_mesh,simplexInCavity).getNode(2).getGlobalNode(),
      SimplicesCell(m_simplex_mesh,simplexInCavity).getNode(3).getGlobalNode()
    };

    if(nodesInCavity.size() != 0 )
    {
      for(const auto & nodeInCavity : nodesInCavity)
      {
        pointToConnect.erase(std::remove(pointToConnect.begin(), pointToConnect.end(), nodeInCavity), pointToConnect.end());
      }
      if(pointToConnect.size() == 3)
      {
          v.push_back(pointToConnect);
      }
    }
    else // ca pue ca a revoir
    {
      v.push_back(std::vector<TInt>{pointToConnect[0], pointToConnect[1], pointToConnect[2]});
      v.push_back(std::vector<TInt>{pointToConnect[1], pointToConnect[2], pointToConnect[3]});
      v.push_back(std::vector<TInt>{pointToConnect[2], pointToConnect[3], pointToConnect[0]});
      v.push_back(std::vector<TInt>{pointToConnect[3], pointToConnect[0], pointToConnect[1]});
    }
  }


  return std::move(v);
}*/
/******************************************************************************/
std::vector<std::vector<TInt>> CavityOperator::CavityIO::CavityIO::pointToConnect()
{
  unsigned int sizeNodeInSimplex = 4;
  std::vector<std::vector<TInt>> v{};
  const std::vector<TInt>& nodesInCavity = getNodeInCavity();
  TInt errorId = std::numeric_limits<int>::min();

  for(const auto & simplexInCavity : simplexInCavityUnmodifiable())
  {
    std::vector<TInt> pointToConnect =
    {
      SimplicesCell(m_simplex_mesh,simplexInCavity).getNode(0).getGlobalNode(),
      SimplicesCell(m_simplex_mesh,simplexInCavity).getNode(1).getGlobalNode(),
      SimplicesCell(m_simplex_mesh,simplexInCavity).getNode(2).getGlobalNode(),
      SimplicesCell(m_simplex_mesh,simplexInCavity).getNode(3).getGlobalNode()
    };

    for(TInt node = 0; node < sizeNodeInSimplex; node++)
    {
        TSimplexID simplexOpp = m_simplex_mesh->getOppositeCell(pointToConnect[node], simplexInCavity);
        if(simplexOpp < 0 && simplexOpp != errorId) /* triangle was founded*/
        {
          std::vector<TSimplexID>&& simplexAdjtoTri = SimplicesTriangle(m_simplex_mesh, -simplexOpp).neighborSimplex();
          simplexAdjtoTri.erase(std::find(simplexAdjtoTri.begin(), simplexAdjtoTri.end(), simplexInCavity),simplexAdjtoTri.end());
          simplexOpp = simplexAdjtoTri[0];
        }

        if(std::find(simplexInCavityUnmodifiable().cbegin(), simplexInCavityUnmodifiable().cend(), simplexOpp) == simplexInCavityUnmodifiable().cend())
        {
          v.push_back(std::vector<TInt>{pointToConnect[(node +1) % 4], pointToConnect[(node + 2)% 4], pointToConnect[(node + 3) % 4]});
        }
    }
  }
  return std::move(v);
}
/******************************************************************************/
CavityOperator::CavityIO CavityOperator::cavity(std::vector<TSimplexID>& initCavity, const simplicesNode::SimplicesNode& nodeToInsert, const CriterionRAIS& criterion, const std::vector<TSimplexID> markedSimplex)
{
  return std::move(cavity(initCavity, nodeToInsert, nodeToInsert, criterion, markedSimplex));
}
/**********************************************************************/
CavityOperator::CavityIO CavityOperator::cavity(std::vector<TSimplexID>& initCavity, const simplicesNode::SimplicesNode& previousNode, const simplicesNode::SimplicesNode& nextNode, const CriterionRAIS& criterion, const std::vector<TSimplexID> markedSimplex)
{
  int errorId =  std::numeric_limits<int>::min();
  gmds::math::Point nextPt = nextNode.getCoords();
  std::vector<TSimplexID> bords;
  size_t sizeLocalNode = 4;
  for(TInt idx = 0 ; idx < initCavity.size() ; idx++)
  {
    TSimplexID simplexId= initCavity[idx];
    TInt nodeIndexLocal = -1;

    if(simplicesCell::SimplicesCell(m_simplex_mesh, simplexId).containNode(previousNode, nodeIndexLocal))
    {
      if(criterion.execute(m_simplex_mesh, simplexId, nodeIndexLocal, nextPt))
      {
        TSimplexID nextSimplexToAdd = simplicesCell::SimplicesCell(m_simplex_mesh, simplexId).oppositeTetraIdx(simplicesCell::SimplicesCell(m_simplex_mesh, simplexId).getNode(nodeIndexLocal));
        if(std::find(initCavity.begin(), initCavity.end(), nextSimplexToAdd) == initCavity.end() && nextSimplexToAdd != errorId &&
           std::find(markedSimplex.begin(), markedSimplex.end(), nextSimplexToAdd) == markedSimplex.end())
        {
          initCavity.resize(initCavity.size() + 1);
          initCavity[initCavity.size() - 1] = nextSimplexToAdd;
        }
      }
    }
    else
    {
      for(nodeIndexLocal = 0; nodeIndexLocal < sizeLocalNode; nodeIndexLocal++)
      {
        if(criterion.execute(m_simplex_mesh, simplexId, nodeIndexLocal, nextPt))
        {
            TSimplexID nextSimplexToAdd = simplicesCell::SimplicesCell(m_simplex_mesh, simplexId).oppositeTetraIdx(simplicesCell::SimplicesCell(m_simplex_mesh, simplexId).getNode(nodeIndexLocal));

            if(std::find(initCavity.begin(), initCavity.end(), nextSimplexToAdd) == initCavity.end() && nextSimplexToAdd != errorId &&
               std::find(markedSimplex.begin(), markedSimplex.end(), nextSimplexToAdd) == markedSimplex.end())
            {
              initCavity.resize(initCavity.size() + 1);
              initCavity[initCavity.size() - 1] = nextSimplexToAdd;
            }
        }
      }
    }
  }

  return std::move(CavityIO(m_simplex_mesh, std::move(initCavity)));
}
