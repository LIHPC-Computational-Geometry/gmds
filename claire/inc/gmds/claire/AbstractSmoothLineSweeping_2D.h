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
	AbstractSmoothLineSweeping_2D(Blocking2D::Block* B, int nb_max_it);
	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/

 protected:
	/*-------------------------------------------------------------------*/
	/** @brief Find the middle of a branch between 3 points
	 */
	math::Point FindMidBranche(const math::Point A, const math::Point B, const math::Point C);
	/*-------------------------------------------------------------------*/

 protected:
	/** block we work on */
	Blocking2D::Block* m_B;
	/** nb max iterations */
	int m_nb_max_iterations;

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
