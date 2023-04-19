#ifndef OCTREE_H
#define OCTREE_H
////////////////////////////////////////////////////////////////////////////////
#include <gmds/math/Point.h>
#include <gmds/math/Hexahedron.h>
////////////////////////////////////////////////////////////////////////////////
#include <gmds/hybridMeshAdapt/CommonInfo.h>
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
#include <gmds/hybridMeshAdapt/SimplicesNode.h>
////////////////////////////////////////////////////////////////////////////////
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
////////////////////////////////////////////////////////////////////////////////
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
////////////////////////////////////////////////////////////////////////////////
// STL Headers
#include <vector>
#include <unordered_set>

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
             const unsigned int numbersMaxSimplices,
             std::string strName = "Octree");

      Octree(SimplexMesh* simplexMesh,
              unsigned int subdivision,
              double xmin, double xmax,
              double ymin, double ymax,
              double zmin, double zmax);

      void addNode(TInt node, std::unordered_set<TInt>& seen);

      ~Octree();

      void initialize();

      void preprocess();

      bool belongToOc(const math::Point& pt) const ;


      TSimplexID findSimplexNextTo(const math::Point& pt);

      std::vector<TInt> findNodesNextTo(const math::Point& pt);

      std::vector<TSimplexID> findSimplicesInOc(const math::Point& pt);

      std::vector<TSimplexID> findTriangleInOc(const math::Point& pt);

      void writeOctree(Mesh& mesh, std::vector<std::vector<Node>>& nodes) const ;

      void setRootOctree(Octree* rootOc){m_rootOc = rootOc;}

      void setParentOctree(Octree* parentOc){m_parentOc = parentOc;}

      std::vector<double> getBorderOctree() const;

      std::vector<int> getSimplices() const {return m_simplices;}

      std::vector<int> getTriangles() const {return m_triangles;}

    private:
      double m_xmin, m_xmax, m_ymin, m_ymax, m_zmin, m_zmax;
      unsigned int m_numbersMaxSimplices;
      SimplexMesh* m_simplexMesh;// = nullptr;

      std::vector<TInt> m_nodes;
      std::vector<TInt> m_simplices;
      std::vector<TInt> m_triangles;
      std::vector<Octree*> m_ocs{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
      Octree* m_rootOc = nullptr;
      Octree* m_parentOc = nullptr;
    };
  }
}

#endif //OCTREE_H
