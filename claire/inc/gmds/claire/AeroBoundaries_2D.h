//
// Created by rochec on 23/03/2022.
//

#ifndef GMDS_AEROBOUNDARIES_2D_H
#define GMDS_AEROBOUNDARIES_2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/claire/AbstractAeroBoundaries.h>
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  AeroBoundaries_2D
 *  \brief  Caractéristiques des frontières du maillage 2D.
 */
class LIB_GMDS_CLAIRE_API AeroBoundaries_2D: public AbstractAeroBoundaries {
 public:

	/*------------------------------------------------------------------------*/
	/** \brief Default constructor
	 */
	AeroBoundaries_2D(Mesh *AMesh);

	/*------------------------------------------------------------------------*/
	/** \brief Function to be called for mesh generation
	 */
	AbstractAeroBoundaries::STATUS execute();
	/*------------------------------------------------------------------------*/

 protected:
	/*-------------------------------------------------------------------*/
	/** \brief Marque les noeuds sur les bords dans la marque
	 * m_markBoundaryNodes
	 */
	void MarkBoundariesNodes();
	/*-------------------------------------------------------------------*/
	/** \brief Colorie les bords. Une couleur par bord connexe, 0 pour
	 * l'intérieur du domaine. Couleurs stockées dans m_var_color_bords.
	 */
	void ColoriageBordsConnexes();
	/*-------------------------------------------------------------------*/

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_AEROBOUNDARIES_2D_H
