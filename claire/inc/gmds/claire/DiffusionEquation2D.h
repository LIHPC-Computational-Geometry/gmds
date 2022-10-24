//
// Created by rochec on 24/10/2022.
//

#ifndef GMDS_DIFFUSIONEQUATION2D_H
#define GMDS_DIFFUSIONEQUATION2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <Eigen/Sparse>
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API DiffusionEquation2D
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
	DiffusionEquation2D(Mesh *AMesh, int AmarkFrontNodes_int, int AmarkFrontNodes_out);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();

 private:

	/*-------------------------------------------------------------------*/
	/** @brief Basis function of the reference triangle.
    	*  @param[in] hat_i summit of the reference triangle (1, 2 or 3)
	 */
	double phi_hat(int hat_i, math::Point X_hat);
	/*-------------------------------------------------------------------*/
	/** @brief Grad of basis function of the reference triangle.
    	*  @param[in] hat_i summit of the reference triangle (1, 2 or 3)
	 */
	Eigen::Vector2d grad_phi_hat(int hat_i);
	/*-------------------------------------------------------------------*/
	/** @brief Transformation of the point X_hat in the triangle K.
	 	*  @param[in] triK_id id of the triangle K
    	*  @param[in] X_hat point in the ref triangle K_hat
	 */
	Eigen::Vector2d FK(TCellID triK_id, math::Point X_hat);
	/*-------------------------------------------------------------------*/
	/** @brief Jacobian of the function FK.
	 	*  @param[in] triK_id id of the triangle K
	 */
	Eigen::Matrix2d JFK(TCellID triK_id);
	/*-------------------------------------------------------------------*/
	/** @brief Absolute value of the determinant of the Jacobian of the function FK.
	 	*  @param[in] triK_id id of the triangle K
	 */
	double absdetJFK(TCellID triK_id);
	/*-------------------------------------------------------------------*/
	/** @brief Gradient of the function phi_i (gradphi_i) where i is the hati-ème
	 * summit of K.
	 	*  @param[in] hat_i summit of the reference triangle (1, 2 or 3)
	 	*  @param[in] triK_id id of the triangle K
	 */
	Eigen::Vector2d grad_phi(int i_hat, TCellID triK_id);
	/*-------------------------------------------------------------------*/
 private:
	/** mesh we work on */
	Mesh *m_mesh;
	/** mark on the interior nodes, first boundary */
	int m_markNodes_int;
	/** mark on the exterior nodes, second boundary */
	int m_markNodes_out;
};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_DIFFUSIONEQUATION2D_H
