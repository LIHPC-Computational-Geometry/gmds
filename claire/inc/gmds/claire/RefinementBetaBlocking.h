//
// Created by rochec on 08/08/2022.
//

#ifndef GMDS_REFINEMENTBETABLOCKING_H
#define GMDS_REFINEMENTBETABLOCKING_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include "gmds/ig/Mesh.h"
#include <gmds/ig/Blocking2D.h>
#include <gmds/claire/Params.h>
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API RefinementBetaBlocking
{
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
         *  @param[in] ABlocking2D the blocking 2D
         *  @param[in] Aparams_aero parameters for aero algorithm
         *
	 */
	RefinementBetaBlocking(Blocking2D* ABlocking2D, ParamsAero Aparams_aero);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/

 private:

 private:
	/** Blocking */
	Blocking2D *m_blocking;
	/** Params pour l'aéro */
	ParamsAero m_params_aero;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_REFINEMENTBETABLOCKING_H
