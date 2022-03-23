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
	virtual AbstractAeroBoundaries::STATUS execute();
	/*------------------------------------------------------------------------*/

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_AEROBOUNDARIES_2D_H
