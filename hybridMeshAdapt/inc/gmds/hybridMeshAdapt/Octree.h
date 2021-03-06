#ifndef OCTREE_H
#define OCTREE_H

// gmds Headers
#include <gmds/math/Point.h>
#include <gmds/hybridMeshAdapt/CommonInfo.h>
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
#include <gmds/hybridMeshAdapt/SimplicesNode.h>
////////////////////////////////////////////////////////////////////////////////
// STL Headers
#include <vector>

namespace gmds
{
  namespace hybrid
  {
    class Octree
    {
    public:
      Octree(double xmin, double xmax,
             double ymin, double ymax,
             double zmin, double zmax,
             SimplexMesh* simplexMesh,
             const unsigned int numbersMaxSimplices);

      Octree(SimplexMesh* simplexMesh,
             const unsigned int numbersMaxSimplices);

      ~Octree();

      void initialize();

      void preprocess();

      bool belongToOc(const math::Point& pt);

      TSimplexID findSimplexNextTo(const math::Point& pt);

      std::vector<TInt> findNodesNextTo(const math::Point& pt);

    private:
      double m_xmin, m_xmax, m_ymin, m_ymax, m_zmin, m_zmax;
      unsigned int m_numbersMaxSimplices;
      SimplexMesh* m_simplexMesh;// = nullptr;

      std::vector<TInt> m_nodes;
      std::vector<Octree*> m_ocs{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    };
  }
}

#endif //OCTREE_H
