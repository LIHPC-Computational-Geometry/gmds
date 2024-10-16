//
// Created by rochec on 22/07/2022.
//

#ifndef GMDS_SMOOTHLINESWEEPINGORTHO_H
#define GMDS_SMOOTHLINESWEEPINGORTHO_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_AERO_export.h"
#include <gmds/aero/AbstractSmoothLineSweeping_2D.h>
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  Yao's Line Sweeping algorithm + orthogonalization
 *  \brief
 */
class LIB_GMDS_AERO_API SmoothLineSweepingOrtho: public AbstractSmoothLineSweeping_2D {
 public:

	/*------------------------------------------------------------------------*/
	/** \brief Default constructor
	*	 @param AB a block
	*	 @param Anb_max_it max iteration number
	*	 @param Atheta damping parameter, theta is in [0,1]
	 */
	SmoothLineSweepingOrtho(Blocking2D::Block* AB, int Anb_max_it, double Atheta);

	/*------------------------------------------------------------------------*/

 protected:
	/*-------------------------------------------------------------------*/
	/** @brief Compute the new position of the node
	 * @param i index of the node in the block
	 * @param j index of the node in the block
	 */
	math::Point ComputeNewPosition(int i, int j) override;
	/*-------------------------------------------------------------------*/

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Compute the orthogonal node
	 * @param i index of the node in the block
	 * @param j index of the node in the block
	 */
	math::Point ComputeOrtho(int i, int j);
	/*-------------------------------------------------------------------*/

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_SMOOTHLINESWEEPINGORTHO_H
