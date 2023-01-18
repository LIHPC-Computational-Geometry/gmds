//
// Created by rochec on 18/11/2022.
//

#ifndef GMDS_AEROEXTRUSION_3D_H
#define GMDS_AEROEXTRUSION_3D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include "gmds/ig/Mesh.h"
#include <gmds/claire/AeroException.h>
#include <gmds/claire/Front_3D.h>
#include <gmds/claire/Params.h>
#include <gmds/claire/FastLocalize.h>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API AeroExtrusion_3D
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
	AeroExtrusion_3D(Mesh *AMeshT, Mesh *AMeshH, ParamsAero Aparams_aero, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/

	struct face_info
	{
		TCellID     					f_id;
		std::map<TCellID, TCellID> next_ideal_nodes;	// (f_node_id, next_ideal_n_id)
		std::map<TCellID, TCellID> next_nodes;			// (f_node_id, next_n_id)
	};
	struct edge_info
	{
		TCellID     					e_id;
		int 								singularity_type;
		std::map<TCellID, bool>		CORNER_n_face_created;
		std::map<TCellID, bool>		END_n_face_created;
		std::map<std::pair<TCellID, TCellID>, TCellID> CORNER_next_nodes;			// ((f_adj_id, e_node_id), next_n_id)
		std::map<TCellID, TCellID> diag_next_node;
	};

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
	 	* \param[in] A_Front_IN front faces and nodes
   	* \param[in] A_distance the distance variable on the first mesh
	 	* \param[in] A_vectors le champ de vecteurs à utiliser
		*
		* \return  the first front computed
	 */
	Front_3D Compute1stLayer(Front_3D A_Front_IN, Variable<double>* A_distance, Variable<math::Vector3d>* A_vectors);
	/*-------------------------------------------------------------------*/
	/** @brief Init the faces struct info for the Front, according to the
	 	* ideal position of each node computed.
	 	*
	 	* \param[in] Front front en entrée
	 	* \param[in] map_new_nodes the map of new nodes
		*
		* \return
	 */
	void InitFaceStructInfo(Front_3D &Front, std::map<TCellID, TCellID> map_new_nodes);
	/*-------------------------------------------------------------------*/
	/** @brief Init the edges struct info for the Front, according to the edges
	 	* classification.
	 	* \param[in] Front front en entrée
		*
		* \return
	 */
	void InitEdgeStructInfo(Front_3D &Front);
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
	std::map<TCellID, int> getSingularEdges(Front_3D &AFront, Variable<int>* front_edges_classification, int mark_singEdgesTreated);
	/*-------------------------------------------------------------------*/
	/** @brief Hexa insertion at the node n_id.
	 	* \param[in] AFront the front
   	* \param[in] n_id the node
   	* \param[in] map_new_nodes map with the ideal next nodes of the front
		*
		* \return  the id of the hexa
	 */
	TCellID TemplateNode3Corner(Front_3D &AFront, TCellID n_id, std::map<TCellID, TCellID> map_new_nodes, double dc);
	/*-------------------------------------------------------------------*/
	/** @brief Hexa insertion at the node n_id.
	 	* \param[in] AFront the front
   	* \param[in] n_id the node
   	* \param[in] mark_edgesTreated mark the edges that can't insert or collapse an hex now
   	* \param[in] mark_facesTreated mark the faces that can't insert an hex now
		*
		* \return  the id of the new hexa
	 */
	TCellID TemplateNode2Corner1End(Front_3D &AFront, TCellID n_id, double dc, int mark_edgesTreated, int mark_facesTreated);
	/*-------------------------------------------------------------------*/
	/** @brief Hexa insertion at the node n_id.
	 	* \param[in] AFront the front
   	* \param[in] n_id the node
   	* \param[in] mark_edgesTreated mark the edges that can't insert or collapse an hex now
   	* \param[in] mark_facesTreated mark the faces that can't insert an hex now
		*
		* \return  the id of the new hexa
	 */
	TCellID TemplateNode1Corner2End(Front_3D &AFront, TCellID n_id, double dc, int mark_edgesTreated, int mark_facesTreated);
	/*-------------------------------------------------------------------*/
	/** @brief Hexa insertion at the node n_id.
	 	* \param[in] AFront the front
   	* \param[in] n_id the node
   	* \param[in] mark_edgesTreated mark the edges that can't insert or collapse an hex now
   	* \param[in] mark_facesTreated mark the faces that can't insert an hex now
		*
		* \return  the id of the new hexa
	 */
	TCellID TemplateNode3End(Front_3D &AFront, TCellID n_id, int mark_edgesTreated, int mark_facesTreated);
	/*-------------------------------------------------------------------*/
	/** @brief Advancing front template on edge classified as corner.
	 	* \param[in] AFront the front
   	* \param[in] e_id the node
		*
		* \return  the id of the hexa
	 */
	TCellID TemplateEdgeCorner(Front_3D &AFront, TCellID e_id, double dc);
	/*-------------------------------------------------------------------*/
	/** @brief Advancing front template on edge classified as end.
	 	* \param[in] AFront the front
   	* \param[in] e_id the node
   	* \param[in] mark_edgesTreated mark the edges that can't insert or collapse an hex now
		* \param[in] mark_facesTreated mark the faces that can't insert an hex now
		*
		* \return  the id of the hexa
	 */
	TCellID TemplateEdgeEnd(Front_3D &AFront, TCellID e_id, double dc, int mark_edgesTreated, int mark_facesTreated);
	/*-------------------------------------------------------------------*/
	/** @brief Créé un hex normal sur la couche à partir d'une face
	 	* \param[in] f_id la face concernée
	 	* \param[in] Front_IN front en entrée
	 	* \param[in] map_new_nodes map with the ideal next nodes of the front
		*
		* \return
	 */
	void TemplateFace(TCellID f_id, Front_3D &Front_IN, std::map<TCellID, TCellID> map_new_nodes);
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
	/** Infos sur les noeuds auxquels se connecter pour chaque face du front */
	std::map<TCellID, face_info> m_FaceInfo;
	/** Infos sur les noeuds auxquels se connecter pour chaque edge du front */
	std::map<TCellID, edge_info> m_EdgeInfo;
	/** Compteur d'hexa */
	int m_iteration;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_AEROEXTRUSION_3D_H
