//
// Created by rochec on 28/11/23.
//

#ifndef GMDS_CONTROLPOINTSSMOOTHING_2D_H
#define GMDS_CONTROLPOINTSSMOOTHING_2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_AERO_export.h"
#include <gmds/ig/Blocking2D.h>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_AERO_API ControlPointsSmoothing_2D {
 public:
	/*--------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS
	} STATUS;

	/*-------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param ACtrlPoints Control Points of a 2D Blocking
	 */
	ControlPointsSmoothing_2D(Blocking2D* ACtrlPoints);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/
 private:

 private:
	/** Control Points we work on */
	Blocking2D *m_CtrlPts;
};
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/

#endif     // GMDS_CONTROLPOINTSSMOOTHING_2D_H
