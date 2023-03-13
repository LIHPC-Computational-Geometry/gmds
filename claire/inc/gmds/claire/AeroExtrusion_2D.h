//
// Created by rochec on 14/04/2022.
//

#ifndef GMDS_AEROEXTRUSION_2D_H
#define GMDS_AEROEXTRUSION_2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include "gmds/ig/Mesh.h"
#include <gmds/claire/AeroException.h>
#include <gmds/claire/Front.h>
#include <gmds/claire/Params.h>
#include <gmds/claire/FastLocalize.h>
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API AeroExtrusion_2D
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
         *  @param[in] AMeshQ the quad mesh to generate
         *  @param[in] Aparams_aero parameters for aero algorithm
         *  @param[in] A_VectorField vector field for extrusion
         *
	 */
	AeroExtrusion_2D(Mesh *AMeshT, Mesh *AMeshQ, ParamsAero& Aparams_aero, Variable<math::Vector3d>* A_VectorField);

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
	std::map<TCellID, TCellID> ComputeIdealPositions(Front AFront, double dist_cible, Variable<double>* A_distance, Variable<math::Vector3d>* A_vectors);
	/*-------------------------------------------------------------------*/
	/** @brief Construit la première couche de blocs. Pour cette couche,
	 	* les conditions sont particulières.
   	* \param[in] A_distance the distance variable on the first mesh
		* \param[in] dist_cible the distance wanted for the first layer
	 	* \param[in] A_vectors le champ de vecteurs à utiliser
		*
		* \return  the first front computed
	 */
	Front Compute1stLayer(Variable<double>* A_distance, double dist_cible, Variable<math::Vector3d>* A_vectors);
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
	Front ComputeLayer(Front Front_IN, Variable<double>* A_distance, double dist_cible, Variable<math::Vector3d>* A_vectors);
	/*-------------------------------------------------------------------*/
	/** @brief Retourne un noeud singulier du front.
	 	* \param[in] Front_IN the front
   	* \param[out] n_id l'id du noeud singulier
		* \param[out] type type de la singularité
		*
		* \return
	 */
	void getSingularNode(Front Front_IN, TCellID &node_id, int &type);
	/*-------------------------------------------------------------------*/
	/** @brief Insertion d'un bloc.
	 	* \param[in] Front_IN the front
   	* \param[in] n_id l'id du noeud concerné
		*
		* \return
	 */
	void Insertion(Front &Front_IN, TCellID n_id, Variable<double>* A_distance, double dist_cible, Variable<math::Vector3d>* A_vectors);
	/*-------------------------------------------------------------------*/
	/** @brief Insertion de deux bloc.
	 	* \param[in] Front_IN the front
   	* \param[in] n_id l'id du noeud auquel 2 blocs sont insérés
		*
		* \return
	 */
	void Insertion_Double(Front &Front_IN, TCellID n_id, Variable<double>* A_distance, double dist_cible, Variable<math::Vector3d>* A_vectors);
	/*-------------------------------------------------------------------*/
	/** @brief Fusion de deux blocs.
	 	* \param[in] Front_IN the front
   	* \param[in] n_id l'id du noeud concerné
		*
		* \return
	 */
	void Fusion(Front &Front_IN, TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief Créé un quad normal sur la couche à partir d'une arête
	 	* \param[in] e_id l'arête concernée
	 	* \param[in] Front_IN front en entrée
		*
		* \return
	 */
	void CreateNormalQuad(TCellID e_id, Front &Front_IN);
	/*-------------------------------------------------------------------*/
 private:
	/** triangular mesh we work on */
	Mesh *m_meshT;
	/** k-d tree */
	FastLocalize m_fl;
	/** quad mesh to generate */
	Mesh *m_meshQ;
	/** Params pour l'aéro */
	ParamsAero m_params_aero;
	/** Vector Field for extrusion */
	Variable<math::Vector3d>* m_VectorField;
	/** Compteur d'hexa */
	int m_iteration;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_AEROEXTRUSION_2D_H
