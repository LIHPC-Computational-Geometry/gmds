//
// Created by rochec on 18/11/2022.
//

#ifndef GMDS_AEROEXTRUSION_3D_H
#define GMDS_AEROEXTRUSION_3D_H

/*----------------------------------------------------------------------------*/
#include "GMDSAero_export.h"
#include "gmds/ig/Mesh.h"
#include <gmds/aero/AeroException.h>
#include <gmds/aero/Front_3D.h>
#include <gmds/aero/Params.h>
#include <gmds/aero/FastLocalize.h>
namespace gmds {
/*----------------------------------------------------------------------------*/
class GMDSAero_API AeroExtrusion_3D
{

 public:
	/*--------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS
	} STATUS;

	/*-------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param[in] AMeshT the triangular mesh where we work on
         *  @param[in] AMeshH the quad mesh to generate
         *  @param[in] Aparams_aero parameters for aero algorithm
         *  @param[in] A_DistanceField distance field for extrusion
         *  @param[in] A_VectorField vector field for extrusion
         *
	 */
	AeroExtrusion_3D(Mesh *AMeshT, Mesh *AMeshH, ParamsAero& Aparams_aero, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Calcule la position idéale du prochain noeud pour chaque
	 	* noeud du front.
	 	* \param[in] AFront the front
   	* \param[in] A_distance the distance variable on the first mesh
		* \param[in] dist_cible the distance wanted for the first layer
		*
		* \return  a map with (TCellID, TCellID) for ideal positions
	 */
	std::map<TCellID, TCellID> ComputeIdealPositions(Front_3D AFront, double dist_cible, Variable<double>* A_distance, Variable<math::Vector3d>* A_vectors);
	/*-------------------------------------------------------------------*/
	/** @brief Construit une couche de mailles à partir d'un front. Ici,
	 	* des mailles peuvent être fusionnées ou insérées.
	 	* \param[in] Front_IN the front before
   	* \param[in] A_distance the distance variable on the first mesh
		* \param[in] dist_cible the distance wanted for the first layer
	 	* \param[in] A_vectors le champ de vecteurs à utiliser
		*
		* \return  the front computed
	 */
	Front_3D ComputeLayer(Front_3D Front_IN, Variable<double>* A_distance, double dist_cible, Variable<math::Vector3d>* A_vectors);
	/*-------------------------------------------------------------------*/
	/** @brief Construit la première couche de blocs. Pour cette couche,
	 	* les conditions sont particulières.
	 	* \param[in] Front_IN front faces and nodes
   	* \param[in] A_distance the distance variable on the first mesh
	 	* \param[in] A_vectors le champ de vecteurs à utiliser
		*
		* \return  the first front computed
	 */
	Front_3D Compute1stLayer(Front_3D Front_IN, Variable<double>* A_distance, Variable<math::Vector3d>* A_vectors);
	/*-------------------------------------------------------------------*/
	/** @brief
	 	* \param[in] Front_IN the front before
   	* \param[in] A_distance the distance variable on the first mesh
		*
		* \return
	 */
	double ComputeMaxDistOnFront(Front_3D Front_IN, Variable<double>* A_distance);
	/*-------------------------------------------------------------------*/
	/** @brief Init the next Front from a Front_IN.
	 	* \param[in] Front_IN front en entrée
		*
		* \return
	 */
	Front_3D InitFrontOUT(Front_3D &Front_IN);
	/*-------------------------------------------------------------------*/
	/** @brief Return, in a map, the singular nodes of the front, and the
	 	* template to apply.
	 	* \param[in] AFront the front
   	* \param[in] front_edges_classification the front edges classification
		*
		* \return  a map with (TCellID, int) for the template to apply on each singular node
	 */
	std::map<TCellID, int> getSingularNodes(Front_3D &AFront, Variable<int>* front_edges_classification);
	/*-------------------------------------------------------------------*/
	/** @brief Return, in a map, the singular edges of the front, and the
	 	* template to apply.
	 	* \param[in] AFront the front
   	* \param[in] front_edges_classification the front edges classification
   	* \param[in] mark_singEdgesTreated mark the already treated edges
		*
		* \return  a map with (TCellID, int) for the template to apply on each singular edge
	 */
	std::map<TCellID, int> getSingularEdges(Front_3D &AFront, Variable<int>* front_edges_classification, TInt mark_singEdgesTreated);
	/*-------------------------------------------------------------------*/
	/** @brief Write a file for debug
		*
		* \return
	 */
	void WriteVTKforDebug();
	/*-------------------------------------------------------------------*/
 private:
	/** triangular mesh we work on */
	Mesh *m_meshT;
	/** k-d tree */
	FastLocalize m_fl;
	/** quad mesh to generate */
	Mesh *m_meshH;
	/** Params pour l'aéro */
	ParamsAero m_params_aero;
	/** Distance Field for extrusion */
	Variable<double>* m_DistanceField;
	/** Vector Field for extrusion */
	Variable<math::Vector3d>* m_VectorField;
	/** Compteur d'hexa */
	int m_iteration;
	/** Variable on patterns */
	Variable<int>* m_Patterns;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_AEROEXTRUSION_3D_H
