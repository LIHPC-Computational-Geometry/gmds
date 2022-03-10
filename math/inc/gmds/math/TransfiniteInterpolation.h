/*----------------------------------------------------------------------------*/
#ifndef GMDS_TRANSFINITEINTERPOLATION_H
#define GMDS_TRANSFINITEINTERPOLATION_H
/*----------------------------------------------------------------------------*/
// gmds file headers
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Point.h>
#include "GMDSMath_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
    namespace math{
/*----------------------------------------------------------------------------*/
/** \class TransfiniteInterpolation
 *  \brief This class gathers algorithms to compute transfinite interpolation
 *         for surfacic and  blocks
 */
/*----------------------------------------------------------------------------*/
        class GMDSMath_API TransfiniteInterpolation
        {

        public:

            /*------------------------------------------------------------------------*/
            /** \brief  Compute a transfinite interpolation from a grid of points.
             * \param[in/out] AGrid the points we compute It must be initialized with values
             *                      on the grid boundary
             */
            static bool compute(std::vector<std::vector< Point > > & AGrid);


/**
  Transfinite volume meshes

                     a0   s0 s1  f0  s0 s1 s5 s4
   s7        s6      a1   s1 s2  f1  s1 s2 s6 s5
    *-------*       a2   s3 s2   f2  s3 s2 s6 s7
     |\s4    |\      a3   s0 s3  f3  s0 s3 s7 s4
     | *-------* s5  a4   s4 s5  f4  s0 s1 s2 s3
     | |   s2| |     a5   s5 s6  f5  s4 s5 s6 s7
  s3 *-|-----* |     a6   s7 s6
      \|      \|     a7   s4 s7
       *-------*     a8   s0 s4
       s0       s1    a9   s1 s5
               a10  s2 s6
                      a11  s3 s7

  - onlvolume right now
  - the faces have to be meshed previously by a grid
*/

            static bool compute(std::vector<std::vector<std::vector<Point>>>& AGrid);

        };
/*----------------------------------------------------------------------------*/
    } // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif //GMDS_TRANSFINITEINTERPOLATION_H
/*----------------------------------------------------------------------------*/


