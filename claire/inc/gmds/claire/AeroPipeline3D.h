//
// Created by rochec on 09/02/2022.
//

#ifndef GMDS_AEROPIPELINE3D_H
#define GMDS_AEROPIPELINE3D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/claire/AbstractAeroPipeline.h>
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  AeroPipeline2D
 *  \brief  Pipeline de génération de maillages 3D pour l'aéro.
 */
class LIB_GMDS_CLAIRE_API AeroPipeline3D: public AbstractAeroPipeline {
 public:

	/*------------------------------------------------------------------------*/
	/** \brief Default constructor
	 */
	AeroPipeline3D();
	/*------------------------------------------------------------------------*/
	/** \brief Function to be called for mesh generation
	 */
	virtual void execute();
	/*-------------------------------------------------------------------*/

 private:


};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_AEROPIPELINE3D_H
