//
// Created by rochec on 20/07/2022.
//

#ifndef GMDS_ABSTRACTSMOOTHLINESWEEPING_2D_H
#define GMDS_ABSTRACTSMOOTHLINESWEEPING_2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Blocking2D.h>
#include <string>
#include <map>
#include <fstream>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** \class  AbstractSmoothLineSweeping_2D
 *  \brief  Classe abstraite pour l'algorithme de lissage.
 */
class LIB_GMDS_CLAIRE_API AbstractSmoothLineSweeping_2D{

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

	/*-------------------------------------------------------------------*/
	/** @brief Find the point at alpha of a branch between 3 points
	 */
	math::Point WeightedPointOnBranch(const math::Point A, const math::Point B, const math::Point C, double alpha);
	/*-------------------------------------------------------------------*/
	/** @brief Compute the L2 norm of the relative error between the pos
	 * of the nodes on the bloc and the old coords
	 */
	double L2_norm_relative_error(Array2D<TCoord>* old_coord_x, Array2D<TCoord>* old_coord_y);
	/*-------------------------------------------------------------------*/
	/** @brief Update old coords
	 */
	void Update_old_coords(Array2D<TCoord>* old_coord_x, Array2D<TCoord>* old_coord_y);
	/*-------------------------------------------------------------------*/

 protected:
	/** block we work on */
	Blocking2D::Block* m_B;
	/** discretization I */
	int m_Nx;
	/** discretization J */
	int m_Ny;
	/** old coord x */
	Array2D<TCoord> m_old_coord_x;
	/** old coord y */
	Array2D<TCoord> m_old_coord_y;
	/** nb max iterations */
	int m_nb_max_iterations;
	/** tolerance */
	int m_tol;
	/** damping parameter */
	double m_theta;

	typedef struct {
		unsigned int val[3][3];
	} stencil;
	/** stencil for each node id */
	std::map<TCellID, stencil> m_stencil;


};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_ABSTRACTSMOOTHLINESWEEPING_2D_H
