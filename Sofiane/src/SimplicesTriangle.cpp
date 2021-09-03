#include "gmds/sofiane/SimplicesTriangle.h"
#include "gmds/sofiane/SimplexMesh.h"
/*---------------------------------------------------------------------------*/
using namespace gmds;
using namespace math;
/*---------------------------------------------------------------------------*/
namespace gmds
{
  namespace hybrid
  {
    namespace simplicesTriangle
    {
      SimplicesTriangle::SimplicesTriangle()
      {

      }

      SimplicesTriangle::~SimplicesTriangle()
      {

      }

      SimplicesTriangle::SimplicesTriangle(SimplexMesh* simplexMesh, const TSimplexID indexTriangle)
      {
        m_simplex_mesh = simplexMesh;
        m_simplexId = std::abs(indexTriangle);
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

      std::vector<TSimplexID> SimplicesTriangle::neighborTriangle (const TInt indexNodeGlobal)
      {
        if(m_simplex_mesh != nullptr)
        {

        }
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
      std::vector<unsigned int> SimplicesTriangle::nodes()
      {
        std::vector<unsigned int> v{};
        v.resize(3);
        v[0] = (unsigned int)m_simplex_mesh->m_tri_nodes[m_simplexId][0];
        v[1] = (unsigned int)m_simplex_mesh->m_tri_nodes[m_simplexId][1];
        v[2] = (unsigned int)m_simplex_mesh->m_tri_nodes[m_simplexId][2];

        return std::move(v);
      }
      /*---------------------------------------------------------------------------*/
      std::vector<TSimplexID> SimplicesTriangle::neighborSimplex()
      {
        std::vector<TSimplexID> v{};
        v.resize(2);
        v[0] = m_simplex_mesh->m_tri_nodes[m_simplexId][3];
        v[1] = m_simplex_mesh->m_tri_adj[m_simplexId][3];
        return std::move(v);
      }

      TInt SimplicesTriangle::getLocalNode(const TInt globalNode)
      {
        int errorId    = std::numeric_limits<int>::min();
        TInt localNode = errorId;
        size_t sizeTriangleBuffer = 3;
        std::vector<TSimplexID> nodesGlobal{m_simplex_mesh->m_tri_nodes[m_simplexId][0], m_simplex_mesh->m_tri_nodes[m_simplexId][1], m_simplex_mesh->m_tri_nodes[m_simplexId][2]};

        for(TInt nodeLocal = 0; nodeLocal < sizeTriangleBuffer; nodeLocal++)
        {
          TInt nodeGlobal = nodesGlobal[nodeLocal];
          if(nodeGlobal == globalNode)
          {
            localNode = nodeLocal;
          }
        }

        return localNode;
      }


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
    }
  }
}
