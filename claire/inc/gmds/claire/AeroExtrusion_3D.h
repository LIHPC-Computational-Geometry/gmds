//
// Created by rochec on 18/11/2022.
//

#ifndef GMDS_AEROEXTRUSION_3D_H
#define GMDS_AEROEXTRUSION_3D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include "gmds/ig/Mesh.h"
#include <gmds/claire/AeroException.h>
#include <gmds/claire/Front_3D.h>
#include <gmds/claire/Params.h>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API AeroExtrusion_3D
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
         *  @param[in] AMeshT the triangular mesh where we work on
         *  @param[in] AMeshH the quad mesh to generate
         *  @param[in] Aparams_aero parameters for aero algorithm
         *  @param[in] A_DistanceField distance field for extrusion
         *  @param[in] A_VectorField vector field for extrusion
         *
	 */
	AeroExtrusion_3D(Mesh *AMeshT, Mesh *AMeshH, ParamsAero Aparams_aero, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/

 private:

 private:
	/** triangular mesh we work on */
	Mesh *m_meshT;
	/** quad mesh to generate */
	Mesh *m_meshH;
	/** Params pour l'aéro */
	ParamsAero m_params_aero;
	/** Distance Field for extrusion */
	Variable<double>* m_DistanceField;
	/** Vector Field for extrusion */
	Variable<math::Vector3d>* m_VectorField;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_AEROEXTRUSION_3D_H
