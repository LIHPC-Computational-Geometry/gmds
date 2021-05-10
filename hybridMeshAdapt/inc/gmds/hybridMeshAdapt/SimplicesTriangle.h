#ifndef SIMPLICESTRIANGLE_H_
#define SIMPLICESTRIANGLE_H_
/*----------------------------------------------------------------------------*/
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
#include <gmds/hybridMeshAdapt/CommonInfo.h>
#include <gmds/hybridMeshAdapt/SimplicesNode.h>
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

        std::vector<TSimplexID> neighborTriangle  (const TInt indexNodeGlobal);

        bool containNode                          (const gmds::hybrid::simplicesNode::SimplicesNode& simplicesNode);

        /*Reorient the Tetra (the normal is out of the Tetra ) if the normal is inside the Tetra*/
        void reorientTriangle                     ();

        std::vector<unsigned int> nodes           ();

        std::vector<TSimplexID> neighborSimplex   ();

        TInt getLocalNode                         (const TInt globalNode);

        std::vector<TInt> getOtherNodeInSimplex   (const std::vector<TInt>& v) const;

        TInt getLocalNode                         (const TInt generalIndex) const;
        
      private:
        SimplexMesh* m_simplex_mesh          = nullptr;

        TSimplexID  m_simplexId             = -1;

      };
    }
  }
}


#endif //SIMPLICESTRIANGLE_H_
