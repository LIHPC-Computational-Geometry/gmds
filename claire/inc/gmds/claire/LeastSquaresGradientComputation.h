//
// Created by rochec on 27/01/2022.
//

#ifndef GMDS_LEASTSQUARESGRADIENTCOMPUTATION_H
#define GMDS_LEASTSQUARESGRADIENTCOMPUTATION_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/claire/DistanceMap.h>
#include <string>
#include <map>
#include <fstream>
#include <Eigen/Sparse>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API LeastSquaresGradientComputation
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
	LeastSquaresGradientComputation(Mesh *AMesh, Variable<double>* Adistance);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Construit les matrices pour la résolution du système
	 */
	void buildMatrix(TCellID n_id, Eigen::SparseMatrix<double> &M, Eigen::VectorXd &b);
	/*-------------------------------------------------------------------*/
	/** @brief Donne le vecteur gradient sur un sommet
	 */
	math::Vector3d computeGradientOnSimpleVertex(TCellID node_id);
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

#endif     // GMDS_LEASTSQUARESGRADIENTCOMPUTATION_H
