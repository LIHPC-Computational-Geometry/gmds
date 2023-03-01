//
// Created by rochec on 09/02/2022.
//

#ifndef GMDS_AEROPIPELINE_3D_H
#define GMDS_AEROPIPELINE_3D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/claire/AbstractAeroPipeline.h>
#include <gmds/claire/AeroBoundaries_3D.h>
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  AeroPipeline_2D
 *  \brief  Pipeline de génération de maillages 3D pour l'aéro.
 */
class LIB_GMDS_CLAIRE_API AeroPipeline_3D : public AbstractAeroPipeline {
 public:

	/*------------------------------------------------------------------------*/
	/** \brief Default constructor
	 */
	explicit AeroPipeline_3D(std::string Aparams);
	/*------------------------------------------------------------------------*/
	/** \brief Function to be called for mesh generation
	 */
	AbstractAeroPipeline::STATUS execute() override;
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
	/** \brief Generate the blocking of the geometry surface. It can be read
	 * 	in a .vtk file for instance.
	 */
	void GeometrySurfaceBlockingGeneration();
	/*------------------------------------------------------------------------*/
	/** \brief
	 */
	void PreTraitementMeshTet();
	/*------------------------------------------------------------------------*/

 protected:
	/** Données des bords */
	AeroBoundaries_3D* m_Bnd ;


};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_AEROPIPELINE_3D_H
