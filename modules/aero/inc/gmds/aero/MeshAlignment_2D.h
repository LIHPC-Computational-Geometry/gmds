//
// Created by rochec on 13/03/2023.
//

#ifndef GMDS_MESHALIGNMENT_2D_H
#define GMDS_MESHALIGNMENT_2D_H

/*----------------------------------------------------------------------------*/
#include "GMDSAero_export.h"
#include "gmds/ig/Mesh.h"
#include <gmds/aero/FastLocalize.h>
namespace gmds {
/*----------------------------------------------------------------------------*/
class GMDSAero_API MeshAlignment_2D
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
         *  @param[in] AMeshTri the mesh on which we know the vector field
         *  @param[in] A_VectorField vector field for extrusion
         *  @param[in] AMeshQuad the mesh on which we want to compute the deviation
         *
	 */
	MeshAlignment_2D(Mesh *AMeshTri, Variable<math::Vector3d>* A_VectorField, Mesh *AMeshQuad);
	/*-------------------------------------------------------------------*/
	/** @brief Default destructor.
	 */
	virtual ~MeshAlignment_2D();
	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/
	/** @brief Return the variable which contains the angle deviation
	 */
	Variable<double>* getVariableDeviation();
	/*-------------------------------------------------------------------*/
 private:
	/*-------------------------------------------------------------------*/
	/** @brief Compute the interpolated value of the vector field at point p
	 * in the triangle tri_id.
	 	* \param[in] tri_id id of the triangle where point p is in m_meshTri
   	* \param[in] p point p
		*
		* \return Interpolated vector
	 */
	math::Vector3d ComputeInterpolatedVector(TCellID tri_id, const math::Point& p);
	/*-------------------------------------------------------------------*/
	/** @brief Compute the min angle between a vector and a cross.
	 	* \param[in] v_ref vector
   	* \param[in] v_cross vector of the cross
		*
		* \return Min angle
	 */
	double minAngleBtwVectorandCross(const math::Vector3d& v, const math::Vector3d& v_cross);
	/*-------------------------------------------------------------------*/
	/** @brief Compute the max angle deviation at a single node n of m_meshQuad.
	 	* \param[in] n_id node id
		*
		* \return Max angle deviation
	 */
	double maxAngleDeviationAtNode(TCellID n_id, const math::Vector3d& v_ref);
	/*-------------------------------------------------------------------*/
 private:
	/** Mesh on which we know the vector field */
	Mesh *m_meshTri;
	/** Mesh on which we want to compute the deviation */
	Mesh *m_meshQuad;
	/** The vector field we want check if we are aligned */
	Variable<math::Vector3d>* m_VectorField;
	/** k-d tree */
	FastLocalize m_fl;
	/** The deviation computed from the vector field */
	Variable<double>* m_var_deviation;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_MESHALIGNMENT_2D_H
