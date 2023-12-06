//
// Created by rochec on 09/02/2022.
//

#ifndef GMDS_AEROPIPELINE_3D_H
#define GMDS_AEROPIPELINE_3D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/claire/AbstractAeroPipeline.h>
#include <gmds/claire/AeroBoundaries_3D.h>
#include <gmds/claire/Blocking3D.h>
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
	explicit AeroPipeline_3D(std::string Aparams, std::string &Aworking_dir);
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
	/** \brief Split the tets with the 4 nodes on the wall into 3 different tets.
	 */
	void PreTraitementMeshTet();
	/*------------------------------------------------------------------------*/
	/** \brief  Classify the surface blocking on the geometry.
	 */
	void SurfaceBlockingClassification();
	/*------------------------------------------------------------------------*/
	/** @brief Compute the vector field used for the extrusion.
	 */
	void ComputeVectorFieldForExtrusion();
	/*------------------------------------------------------------------------*/
	/** \brief  Classify all the entities (nodes, edges, faces) of the blocking
	 * (stored in a mesh) on the geometry.
	 */
	void BlockingGeometricClassification();
	/*------------------------------------------------------------------------*/
	/** @brief Init the Blocking3D and the Control Points structure from the
	 * hex mesh built using the AeroPipeline3D algorithm. Initialize the
	 * geometric classification of the block corners too in m_linker_BG.
	 */
	void initBlocking3DfromMesh();
	/*----------------------------------------------------------------------------*/
	/** @brief Update the Variable "Couche" for each inner nodes of the blocking
	 * m_Blocking3D.
	 */
	void updateLayerValues();
	/*----------------------------------------------------------------------------*/
	/** @brief Compute the position of the inner nodes of each block of the
	 * m_Blocking3D, using the parametric space [0,1]x[0,1]x[0,1] for each block
	 * considered as Bezier Hex defined by the control points stored in
	 * m_CtrlPts
	 */
	void computeBlockNodesPositionsFromCtrlPoints();
	/*----------------------------------------------------------------------------*/
	/** @brief Compute the position of the boundary control points stored in
	 * m_CtrlPts, in order to interpolate the physical boundaries of the geometry.
	 */
	void computeControlPointstoInterpolateBoundaries();
	/*----------------------------------------------------------------------------*/

 protected:
	/** Données des bords */
	AeroBoundaries_3D* m_Bnd ;
	/** blocking 3D */
	Blocking3D m_Blocking3D;
	/** control points of the 3D blocking */
	Blocking3D m_CtrlPts;
	/** Linker Blocking à la géométrie */
	cad::GeomMeshLinker* m_linker_BG;

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_AEROPIPELINE_3D_H
