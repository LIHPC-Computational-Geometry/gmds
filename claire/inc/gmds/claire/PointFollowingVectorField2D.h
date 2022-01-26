//
// Created by rochec on 26/01/2022.
//

#ifndef GMDS_POINTFOLLOWINGVECTORFIELD2D_H
#define GMDS_POINTFOLLOWINGVECTORFIELD2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/claire/DistanceMap.h>
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API PointFollowingVectorField2D
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
         *  @param AMesh the mesh where we work on
	 */
	PointFollowingVectorField2D(Mesh *AMesh);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();

 private:


 private:
	/** mesh we work on */
	Mesh *m_mesh;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_POINTFOLLOWINGVECTORFIELD2D_H
