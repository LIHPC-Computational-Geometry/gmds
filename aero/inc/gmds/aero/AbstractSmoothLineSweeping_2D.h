//
// Created by rochec on 20/07/2022.
//

#ifndef GMDS_ABSTRACTSMOOTHLINESWEEPING_2D_H
#define GMDS_ABSTRACTSMOOTHLINESWEEPING_2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_AERO_export.h"
#include <gmds/ig/Blocking2D.h>
#include <gmds/math/Point.h>
#include <gmds/utils/Array.h>
#include <string>
#include <map>
#include <fstream>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** \class  AbstractSmoothLineSweeping_2D
 *  \brief  Classe abstraite pour l'algorithme de lissage.
 */
class LIB_GMDS_AERO_API AbstractSmoothLineSweeping_2D{

 public:
	/*--------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS
	} STATUS;
	/*--------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param
	 */
	AbstractSmoothLineSweeping_2D(Blocking2D::Block* AB, int Anb_max_it, double Atheta);
	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/

 protected:
	/*-------------------------------------------------------------------*/
	/** @brief Compute the new position of the node
	 */
	virtual math::Point ComputeNewPosition(int i, int j)=0;
	/*-------------------------------------------------------------------*/

	/*-------------------------------------------------------------------*/
	/** @brief One step of the smoothing algorithm
	 */
	void One_Step_Smoothing();
	/*-------------------------------------------------------------------*/
	/** @brief Boundary smoothing algorithm
	 */
	void BoundarySlipping();
	/*-------------------------------------------------------------------*/
	/** @brief Compute the L2 norm of the relative error between the pos
	 * of the nodes on the bloc and the old coords
	 */
	double L2_norm_relative_error();
	/*-------------------------------------------------------------------*/
	/** @brief Update new coords
	 */
	void Update_new_coords();
	/*-------------------------------------------------------------------*/

 protected:
	/** block we work on */
	Blocking2D::Block* m_B;
	/** discretization I */
	int m_Nx;
	/** discretization J */
	int m_Ny;
	/** old points */
	Array2D<math::Point> m_P_new;
	/** nb max iterations */
	int m_nb_max_iterations;
	/** tolerance */
	double m_tol;
	/** damping parameter */
	double m_theta;

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_ABSTRACTSMOOTHLINESWEEPING_2D_H
