//
// Created by rochec on 21/01/2022.
//

#ifndef GMDS_GRADIENTCOMPUTATION_3D_H
#define GMDS_GRADIENTCOMPUTATION_3D_H

/*----------------------------------------------------------------------------*/
#include "GMDSAero_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/aero/DistanceMap.h>
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class GMDSAero_API GradientComputation_3D
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
	GradientComputation_3D(Mesh *AMesh, Variable<double>* Adistance);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Donne le vecteur gradient sur une région
	 */
	math::Vector3d computeGradientOnSimpleRegion(TCellID region_id);
	/*-------------------------------------------------------------------*/

 private:
	/** mesh we work on */
	Mesh *m_mesh;
	/** champ scalaire dont on souhaite obtenir le gradient */
	Variable<double>* m_distance;
	/** champ de gradient défini sur chaque face */
	Variable<math::Vector3d>* m_gradient3D;
};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_GRADIENTCOMPUTATION_3D_H
