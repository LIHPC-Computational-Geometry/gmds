//
// Created by rochec on 01/12/23.
//

#ifndef GMDS_TRANSFINITEINTERPOLATION_3D_H
#define GMDS_TRANSFINITEINTERPOLATION_3D_H

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
/** \class TransfiniteInterpolation
 *  \brief This class gathers algorithms to compute transfinite interpolation
 *         for surfacic and  blocks
 */
/*----------------------------------------------------------------------------*/
class GMDSMath_API TransfiniteInterpolation_3D {
 public:
	/** \brief Compute a transfinite interpolation from a volumic grid of points. The faces have
	 		* to be meshed previously by a grid
	 		* \param[in/out] AGrid the points we compute It must be initialized with values
         *                      on the faces of the grid boundary
	 */

	static bool computeHex(Array3D<math::Point>& AGrid);

};
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/

#endif     // GMDS_TRANSFINITEINTERPOLATION_3D_H
