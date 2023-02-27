//
// Created by rochec on 04/02/2022.
//

#ifndef GMDS_ADVECTEDPOINTRK4_2D_H
#define GMDS_ADVECTEDPOINTRK4_2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/claire/FastLocalize.h>
#include <string>
#include <map>
#include <fstream>
#include <Eigen/Sparse>
#include <Eigen/Eigen>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API AdvectedPointRK4_2D
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
	AdvectedPointRK4_2D(Mesh *AMesh, FastLocalize *A_fl, math::Point A_Pstart, double A_d0, Variable<double>* A_distance, Variable<math::Vector3d>* A_gradient2D);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/
	/** @brief Donne le point final
	 */
	math::Point getPend();

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Return true if the point M is in the triangle
	 */
	bool isInTriangle(TCellID face_id, math::Point &M);
	/*-------------------------------------------------------------------*/
	/** @brief Return in which triangle M is
	 */
	TCellID inWhichTriangle(math::Point &M, TCellID f0_id);
	/*-------------------------------------------------------------------*/
	/** @brief Compute the minimal edge's lenght
	 */
	double minEdgeLenght();
	/*-------------------------------------------------------------------*/
	/** @brief Retourne la matrice inverse de A sur la face face_id
	 */
	Eigen::Matrix3d getInvMatrixA(TCellID face_id);
	/*-------------------------------------------------------------------*/
	/** @brief Interpolation de la distance au point M
	 */
	double interpolationDistance(TCellID face_id, Eigen::Matrix3d &Mat_A_Inv, math::Point M);
	/*-------------------------------------------------------------------*/
	/** @brief Interpolation du gradient au point M
	 */
	math::Vector3d interpolationGradient(TCellID face_id, Eigen::Matrix3d &Mat_A_Inv, math::Point M);
	/*-------------------------------------------------------------------*/
	/** @brief Applique le schéma Runge Kutta d'ordre 4 pour résoudre dx/dt = grad
	 */
	math::Point RungeKutta4(math::Point &yn, math::Vector3d grad_yn, double dt);
	/*-------------------------------------------------------------------*/
	/** @brief Write the discrete path in a vtk field
	 */
	void writeDiscretePathInVTK();
	/*-------------------------------------------------------------------*/

 private:
	/** mesh we work on */
	Mesh *m_mesh;
	/** k-d tree */
	FastLocalize* m_fl;
	/** starting point */
	math::Point m_Pstart ;
	/** ending point */
	math::Point m_Pend ;
	/** the distance */
	double m_d0 ;
	/** carte des distances */
	Variable<double> *m_distance;
	/** the vector field to follow */
	Variable<math::Vector3d>* m_gradient2D ;
	/** liste des points intermédiaires */
	std::vector<math::Vector3d> m_discrete_path;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_ADVECTEDPOINTRK4_2D_H
