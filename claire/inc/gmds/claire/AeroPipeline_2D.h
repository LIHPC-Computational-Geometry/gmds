//
// Created by rochec on 09/02/2022.
//

#ifndef GMDS_AEROPIPELINE_2D_H
#define GMDS_AEROPIPELINE_2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/claire/AbstractAeroPipeline.h>
#include <gmds/claire/AeroBoundaries_2D.h>
#include <gmds/ig/Blocking2D.h>
#include <gmds/claire/Params.h>
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  AeroPipeline_2D
 *  \brief  Pipeline de génération de maillages 2D pour l'aéro.
 */
class LIB_GMDS_CLAIRE_API AeroPipeline_2D : public AbstractAeroPipeline {
 public:

	/*------------------------------------------------------------------------*/
	/** \brief Default constructor
	 */
	AeroPipeline_2D(std::string Aparams);

	/*------------------------------------------------------------------------*/
	/** \brief Function to be called for mesh generation
	 */
	virtual AbstractAeroPipeline::STATUS execute();
	/*------------------------------------------------------------------------*/

 private:
	/*----------------------------------------------------------------------------*/
	/** @brief Lecture du fichier de maillage au format .vtk
	 */
	void LectureMaillage();
	/*----------------------------------------------------------------------------*/
	/** @brief Ecritures des fichiers de maillages triangulaires et quad au
	 * format .vtk
	 */
	void EcritureMaillage();
	/*----------------------------------------------------------------------------*/
	/** @brief Créé les sommets des blocs sur le bord de couleur color pour le
	 * maillage quad généré.
	 */
	void DiscretisationParoi(int color);
	/*----------------------------------------------------------------------------*/
	/** @brief Convertisseur m_mQuad rempli des blocs à la structure Blocking2D.
	 */
	void ConvertisseurMeshToBlocking();
	/*----------------------------------------------------------------------------*/
	/** @brief Update the second linker at the node n2 with de datas of the first
	 * linker.
	 * @param linker_1 the reference linker
	 * @param n_1 the node in the mesh on the linker_1
	 * @param linker_2 the linker to update
	 * @param n_2 the node in the second mesh to initialize the classification
	 * @return void
	 */
	void UpdateLinker(cad::GeomMeshLinker* linker_1, Node n_1, cad::GeomMeshLinker* linker_2, Node n_2);
	/*----------------------------------------------------------------------------*/
	/** @brief Update the linker for the last layer (on amont boundary)
	 * @return void
	 */
	void UpdateLinkerLastLayer();
	/*----------------------------------------------------------------------------*/
	/** @brief Transmet la classification des noeuds de blocs aux noeuds intérieurs
	 * aux blocs
	 * @param
	 * @return void
	 */
	void BlockingClassification();
	/*----------------------------------------------------------------------------*/
	/** @brief Compute the vector field used for the extrusion.
	 * @param
	 * @return void
	 */
	void ComputeVectorFieldForExtrusion();
	/*----------------------------------------------------------------------------*/
	/** @brief Mesh refinement.
	 * @param
	 * @return void
	 */
	void MeshRefinement();
	/*----------------------------------------------------------------------------*/
	/** @brief Check the mesh alignement with a vector field.
	 * @param
	 * @return void
	 */
	void MeshAlignement();
	/*----------------------------------------------------------------------------*/
 protected:
	/** blocking 2D */
	Blocking2D m_Blocking2D;
	/** Données des bords */
	AeroBoundaries_2D* m_Bnd;
	/** Linker Blocking à la géométrie */
	cad::GeomMeshLinker* m_linker_BG;

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_AEROPIPELINE_2D_H
