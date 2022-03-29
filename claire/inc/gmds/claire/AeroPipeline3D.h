//
// Created by rochec on 09/02/2022.
//

#ifndef GMDS_AEROPIPELINE3D_H
#define GMDS_AEROPIPELINE3D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/claire/AbstractAeroPipeline.h>
#include <gmds/claire/AeroBoundaries_3D.h>
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
	AeroPipeline3D(ParamsAero Aparams);
	/*------------------------------------------------------------------------*/
	/** \brief Function to be called for mesh generation
	 */
	virtual AbstractAeroPipeline::STATUS execute();
	/*------------------------------------------------------------------------*/

 private:
	/*------------------------------------------------------------------------*/
	/** \brief Function to read the initial mesh
	 */
	void LectureMaillage();
	/*------------------------------------------------------------------------*/
	/** \brief Function to write the hex mesh and the tetra mesh
	 */
	void EcritureMaillage();
	/*------------------------------------------------------------------------*/

 protected:
	/** Tetra mesh we work on */
	//Mesh m_mTetra;
	/** Hexa mesh generated */
	//Mesh m_mHexa;
	/** Données des bords */
	AeroBoundaries_3D* m_Bnd ;


};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_AEROPIPELINE3D_H
