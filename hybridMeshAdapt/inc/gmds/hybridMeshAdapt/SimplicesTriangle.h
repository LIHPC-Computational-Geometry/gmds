#ifndef SIMPLICESTRIANGLE_H_
#define SIMPLICESTRIANGLE_H_
/*----------------------------------------------------------------------------*/
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
#include <gmds/hybridMeshAdapt/CommonInfo.h>
#include <gmds/hybridMeshAdapt/SimplicesNode.h>
#include "gmds/hybridMeshAdapt/SimplicesCell.h"
/*----------------------------------------------------------------------------*/

namespace gmds
{
  namespace hybrid
  {
    class SimplexMesh;

    namespace simplicesTriangle
    {

      class SimplicesTriangle
      {
      public:
        /*------------------------------------------------------------------------*/
        /** \brief Constructor.
         */
        SimplicesTriangle();

        ~SimplicesTriangle();

        SimplicesTriangle(SimplexMesh* simplexMesh, const TSimplexID indexTrianngle);

        double getAreaOfCell                      ();

        TSimplexID neighborTriangle               (const TInt indexNodeLocal) const ;

        std::vector<TSimplexID> adjacentTriangle  () const ;

        bool containNode                          (const gmds::hybrid::simplicesNode::SimplicesNode& simplicesNode);

        bool containNode                          (const std::vector<TInt>& nodes);

        std::vector<TInt> intersectionNodes             (const SimplicesTriangle& triangles);

        /*Reorient the Tetra (the normal is out of the Tetra ) if the normal is inside the Tetra*/
        void reorientTriangle                     ();

        std::vector<unsigned int> nodes           () const;

        std::vector<TInt> getNodes                () const;

        /*triangles's normal are oriented outside the mesh*/
        math::Vector3d    getNormal               ()  const;

        std::vector<TSimplexID> neighborSimplex   () const ;

        std::vector<TInt> otherNodesInTriangle    (const SimplicesTriangle& simpliceTriangle) const ;

        std::vector<TInt> getOtherNodeInSimplex   (const std::vector<TInt>& v) const;

        TInt getLocalNode                         (const TInt generalIndex) const;

        std::vector<TSimplexID> neighborTetra     () const ;

        std::vector<TInt> getOppositeEdge(const unsigned int localNode) const ;

        /*nodeA & node B form the edge we are revolve around*/
        std::vector<TSimplexID> buildclockWiseTrianglesbyShell(const TInt nodeA, const TInt nodeB, const TInt nodeC) const ;

        std::vector<TSimplexID> findclockWiseTrianglesbyShell(const TInt nodeA, const TInt nodeB, const TInt nodeC) const ;

        TSimplexID simplexId() const {return m_simplexId;}

        math::Orientation::Sign orientation(const gmds::math::Point& pt, bool inverseOrientation = false) const;

        friend std::ostream&  operator<<(std::ostream& os, const SimplicesTriangle& simplicesTriangle)
        {
          std::vector<TInt > nodes = simplicesTriangle.getNodes();
          std::vector<TSimplexID> adjSimplices = simplicesTriangle.adjacentTriangle();
          const TInt Node0 = nodes[0];
          const TInt Node1 = nodes[1];
          const TInt Node2 = nodes[2];

          std::vector<TSimplexID> adjSimplex = simplicesTriangle.neighborSimplex();

          os << "simplex ID  = " << - simplicesTriangle.simplexId() << "| Node : " << Node0 << " " << Node1 << " " << Node2  << " direct normal Simplex ->" << adjSimplex.front() << std::endl;
          os << "simplex Adj = " <<  " : " << adjSimplices[0] << " " << adjSimplices[1] << " " << adjSimplices[2] << " indirect normal Simplex -> " << adjSimplex.back() << std::endl;

          return os;
        }

      private:
        SimplexMesh* m_simplex_mesh          = nullptr;

        TSimplexID  m_simplexId             = -1;
      };
    }
  }
}


#endif //SIMPLICESTRIANGLE_H_
