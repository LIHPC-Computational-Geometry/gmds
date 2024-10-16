//
// Created by rochec on 29/05/24.
//

#ifndef GMDS_REFINEMENTBETABLOCK3D_H
#define GMDS_REFINEMENTBETABLOCK3D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_AERO_export.h"
#include "gmds/ig/Mesh.h"
#include <gmds/aero/Blocking3D.h>
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_AERO_API RefinementBetaBlock3D
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
         *  @param[in] ABlock block to refine
         *  @param[in] Abf block face to give the refinement direction
         *  @param[in] Asize_first_edge target size of the first edge
         *
	 */
	RefinementBetaBlock3D(Blocking3D::Block* ABlock,
	                      Blocking3D::BlockFace* Abf,
	                      double Asize_first_edge);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/

 private:

 private:
	/** Block to refine */
	Blocking3D::Block* m_Block;
	/** Block face for the refinement direction */
	Blocking3D::BlockFace* m_bf;
	/** Target size first edge */
	double m_size_first_edge;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_REFINEMENTBETABLOCK3D_H
