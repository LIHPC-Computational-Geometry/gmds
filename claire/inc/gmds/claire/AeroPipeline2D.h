//
// Created by rochec on 09/02/2022.
//

#ifndef GMDS_AEROPIPELINE2D_H
#define GMDS_AEROPIPELINE2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/claire/AbstractAeroPipeline.h>
#include <gmds/claire/Params.h>
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  AeroPipeline2D
 *  \brief  Pipeline de génération de maillages 2D pour l'aéro.
 */
class LIB_GMDS_CLAIRE_API AeroPipeline2D: public AbstractAeroPipeline {
 public:

	/*------------------------------------------------------------------------*/
	/** \brief Default constructor
	 */
	AeroPipeline2D(ParamsAero Aparams);

	/*------------------------------------------------------------------------*/
	/** \brief Function to be called for mesh generation
	 */
	virtual void execute();
	/*------------------------------------------------------------------------*/

 private:
	/*----------------------------------------------------------------------------*/
	/** @brief Lecture du fichier de maillage au format .vtk
	 */
	void LectureMaillage();
	/*----------------------------------------------------------------------------*/
	/** @brief Initialisation des marques sur les fronts
	 */
	void InitialisationFronts();
	/*----------------------------------------------------------------------------*/

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_AEROPIPELINE2D_H
