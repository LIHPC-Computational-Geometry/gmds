//
// Created by rochec on 24/03/2022.
//

#ifndef GMDS_AEROBOUNDARIES_3D_H
#define GMDS_AEROBOUNDARIES_3D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/claire/AbstractAeroBoundaries.h>
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  AeroBoundaries_3D
 *  \brief  Caractéristiques des frontières du maillage 3D.
 */
class LIB_GMDS_CLAIRE_API AeroBoundaries_3D: public AbstractAeroBoundaries {
 public:

	/*------------------------------------------------------------------------*/
	/** \brief Default constructor
	 */
	explicit AeroBoundaries_3D(Mesh *AMesh);

	/*------------------------------------------------------------------------*/

 protected:
	/*-------------------------------------------------------------------*/
	/** \brief Marque les noeuds sur les bords dans la marque
	 * m_markBoundaryNodes
	 */
	void MarkBoundariesNodes() override;
	/*-------------------------------------------------------------------*/
	/** \brief Identifie de quelle couleur est le front Amont dans la
	 * variable m_var_color_bords.
	 */
	void WhichColorIsAmont() override;
	/*-------------------------------------------------------------------*/

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_AEROBOUNDARIES_3D_H
