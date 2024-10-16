//
// Created by rochec on 21/07/2022.
//

#ifndef GMDS_SMOOTHLINESWEEPINGYAO_H
#define GMDS_SMOOTHLINESWEEPINGYAO_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_AERO_export.h"
#include <gmds/aero/AbstractSmoothLineSweeping_2D.h>
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  Yao's Line Sweeping algorithm
 *  \brief
 */
class LIB_GMDS_AERO_API SmoothLineSweepingYao: public AbstractSmoothLineSweeping_2D {
 public:

	/*------------------------------------------------------------------------*/
	/** \brief Default constructor
	*	 @param AB a block
	*	 @param Anb_max_it max iteration number
	*	 @param Atheta damping parameter, theta is in [0,1]
	*/
	SmoothLineSweepingYao(Blocking2D::Block* AB, int Anb_max_it, double Atheta);

	/*------------------------------------------------------------------------*/

 protected:
	/*-------------------------------------------------------------------*/
	/** @brief Compute the new position of the node
	 * @param i index of the node in the block
	 * @param j index of the node in the block
	 */
	math::Point ComputeNewPosition(int i, int j) override;
	/*-------------------------------------------------------------------*/

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_SMOOTHLINESWEEPINGYAO_H
