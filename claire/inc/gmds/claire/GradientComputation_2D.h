//
// Created by rochec on 21/01/2022.
//

#ifndef GMDS_GRADIENTCOMPUTATION_2D_H
#define GMDS_GRADIENTCOMPUTATION_2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/claire/DistanceMap.h>
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API GradientComputation_2D
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
	GradientComputation_2D(Mesh *AMesh, Variable<double>* Adistance);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Retourne le vecteur orthogonal (rotation de pi/2 dans le sens trigo)
	 * au vecteur v0v1
	 */
	math::Vector3d getNormalVector(TCellID n0_id, TCellID n1_id);
	/*-------------------------------------------------------------------*/
	/** @brief Donne le vecteur gradient sur une face
	 */
	math::Vector3d computeGradientOnSimpleFace(TCellID face_id);
	/*-------------------------------------------------------------------*/

 private:
	/** mesh we work on */
	Mesh *m_mesh;
	/** champ scalaire dont on souhaite obtenir le gradient */
	Variable<double>* m_distance;
	/** champ de gradient défini sur chaque face */
	Variable<math::Vector3d>* m_gradient2D;
};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_GRADIENTCOMPUTATION_2D_H
