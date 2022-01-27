//
// Created by rochec on 27/01/2022.
//

#ifndef GMDS_LEASTSQUARESGRADIENTCOMPUTATION2D_H
#define GMDS_LEASTSQUARESGRADIENTCOMPUTATION2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/claire/DistanceMap.h>
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API LeastSquaresGradientComputation2D {
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
	LeastSquaresGradientComputation2D(Mesh *AMesh, Variable<double>* Adistance);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();

 private:

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

#endif     // GMDS_LEASTSQUARESGRADIENTCOMPUTATION2D_H
