//
// Created by rochec on 04/02/2022.
//

#ifndef GMDS_RESLU_H
#define GMDS_RESLU_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <string>
#include <map>
#include <fstream>
#include <Eigen/Sparse>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API ResLU
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
         *  @param
	 */
	ResLU(Eigen::SparseMatrix<double> A_A, Eigen::VectorXd A_b);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();

 public:
	/*-------------------------------------------------------------------*/
	/** @brief Donne le vecteur solution de Ax=b
	 */
	Eigen::VectorXd getSolution();
	/*-------------------------------------------------------------------*/
	/** @brief Donne la matrice LU
	 */
	Eigen::SparseMatrix<double> getLU();
	/*-------------------------------------------------------------------*/

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Décompose la matrice A en L et U
	 */
	void computeLU();
	/*-------------------------------------------------------------------*/
	/** @brief Résout Ax=b
	 */
	void solveLU();
	/*-------------------------------------------------------------------*/
 private:
	/** Matrice du système à résoudre */
	Eigen::SparseMatrix<double> m_A;
	/** Matrice de la décomposition LU */
	Eigen::SparseMatrix<double> m_LU;
	/** Vecteur Ax=b */
	Eigen::VectorXd m_b;
	/** Vecteur solution du système Ax=b */
	Eigen::VectorXd m_x;
	/** Taille de la matrice A */
	int m_n;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_RESLU_H
