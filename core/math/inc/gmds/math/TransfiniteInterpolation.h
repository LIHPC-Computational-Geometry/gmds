/*----------------------------------------------------------------------------*/
#ifndef GMDS_TRANSFINITEINTERPOLATION_H
#define GMDS_TRANSFINITEINTERPOLATION_H
/*----------------------------------------------------------------------------*/
// gmds file headers
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Point.h>
#include <gmds/utils/Array.h>
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
        class GMDSMath_API TransfiniteInterpolation {
        public:
            /** \brief  Compute a transfinite interpolation from a grid of points.
             * \param[in/out] AGrid the points we compute It must be initialized with values
             *                      on the grid boundary */
            static bool compute(std::vector<std::vector< Point > > & AGrid);

            /** \brief  Compute a transfinite interpolation from a grid of points.
             * \param[in/out] AGrid the points we compute It must be initialized with values
             *                      on the grid boundary */
            static bool computeQuad(Array2D<Point>& AGrid);

            /** \brief  Compute a transfinite interpolation from a triangular grid of
             *          points.
             * \param[in/out] AGrid the points we compute It must be initialized with values
             *                      on the grid boundary */
            static bool computeTri(TriArray<Point>& AGrid);

            /** Transfinite volume meshes
             *   - the faces have to be meshed previously by a grid
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


