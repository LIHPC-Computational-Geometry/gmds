//
// Created by rochec on 24/10/2022.
//

#ifndef GMDS_DIFFUSIONEQUATION2D_H
#define GMDS_DIFFUSIONEQUATION2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_AERO_export.h"
#include <gmds/ig/Mesh.h>
#include <Eigen/Sparse>
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_AERO_API DiffusionEquation2D
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
	DiffusionEquation2D(Mesh *AMesh, TInt AmarkFrontNodes_int, TInt AmarkFrontNodes_out, Variable<double>* Adistance);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();

 private:

	/*-------------------------------------------------------------------*/
	/** @brief Basis function of the reference triangle.
    	*  @param[in] hat_i summit of the reference triangle (1, 2 or 3)
	 */
	static double phi_hat(int hat_i, math::Point X_hat);
	/*-------------------------------------------------------------------*/
	/** @brief Grad of basis function of the reference triangle.
    	*  @param[in] hat_i summit of the reference triangle (1, 2 or 3)
	 */
	static Eigen::Vector2d grad_phi_hat(int hat_i);
	/*-------------------------------------------------------------------*/
	/** @brief Transformation of the point X_hat in the triangle K.
	 	*  @param[in] triK_id id of the triangle K
    	*  @param[in] X_hat point in the ref triangle K_hat
	 */
	Eigen::Vector2d FK(TCellID triK_id, const math::Point& X_hat);
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
	/** @brief Quadrature points and weights for the Summit Formula.
	 	*  @param[in]
	 	*  @param[in]
	 */
	static void quadraturePointsAndWeightsSummitpointFormula(std::vector<double> &weights, std::vector<math::Point> &points);
	/*-------------------------------------------------------------------*/

	/*-------------------------------------------------------------------*/
	/** @brief Assemble the stiffness and mass matrix for the resolution.
	 */
	void assembleMassAndStiffnessMatrices();
	/*-------------------------------------------------------------------*/
	/** @brief Apply the Boundary Conditions on the system matrix.
	 */
	void applyBCToSystemMatrix(Eigen::SparseMatrix<double,Eigen::RowMajor> & systemMatrix);
	/*-------------------------------------------------------------------*/
	/** @brief Apply the Dirichlet Boundary Conditions on the second member.
	 */
	void applyBCToSecondMember(Eigen::SparseVector<double> & secondMember);
	/*-------------------------------------------------------------------*/

	/*-------------------------------------------------------------------*/
	/** @brief Initialisation for the time scheme.
	 */
	void initialisation();
	/*-------------------------------------------------------------------*/
	/** @brief One iteration of the time scheme.
	 */
	void oneTimeStep();
	/*-------------------------------------------------------------------*/
 private:
	/** mesh we work on */
	Mesh *m_mesh;
	/** mark on the interior nodes, first boundary */
	TInt m_markNodes_int;
	/** mark on the exterior nodes, second boundary */
	TInt m_markNodes_out;
	/** carte des distances par rapport au front concerné */
	Variable<double>* m_distance;

	/** We store the correspondance from node ids to local index in m_nodes */
	std::map<gmds::TCellID, int> m_id_local_index;

	/** Stiffness Matrix */
	Eigen::SparseMatrix<double,Eigen::RowMajor> m_stiffness;
	/** Mass Matrix */
	Eigen::SparseMatrix<double,Eigen::RowMajor> m_mass;

	/** Time */
	double m_t;
	/** Time step */
	double m_dt;
	/** Iteration */
	int m_it;
	/** Iteration max */
	int m_it_max;
	/** Diffusion coefficient */
	double m_sigma;

	/** Solveur */
	Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> > m_solver;

	/** Initial conditions */
	Eigen::SparseVector<double> m_sol_0;
	/** Solution at time n */
	Eigen::SparseVector<double> m_sol_n;
};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_DIFFUSIONEQUATION2D_H
