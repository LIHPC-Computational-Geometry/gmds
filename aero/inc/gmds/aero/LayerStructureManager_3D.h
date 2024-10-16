//
// Created by rochec on 27/04/23.
//

#ifndef GMDS_LAYERSTRUCTUREMANAGER_3D_H
#define GMDS_LAYERSTRUCTUREMANAGER_3D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_AERO_export.h"
#include "gmds/ig/Mesh.h"
#include <gmds/aero/Front_3D.h>
#include <gmds/aero/Params.h>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_AERO_API LayerStructureManager_3D
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
    *  @param[in] AMeshH the quad mesh to generate
    *  @param[in] AFront the front
    *  @param[in] A_map_new_nodes the map of the TCellID of the ideal next node
    *  				for each node of the front AFront
    *
	 */
	LayerStructureManager_3D(Mesh *AMeshH, Front_3D *AFront, std::map<TCellID, TCellID> &A_map_new_nodes);
	/*-------------------------------------------------------------------*/
	/** @brief Destructor.
    *  @param
	 */
	~LayerStructureManager_3D();
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

		std::map<TCellID, bool>		REVERSAL_n_faces_created;
		std::map<std::pair<TCellID, TCellID>, TCellID> REVERSAL_diag_nodes;
		std::map<std::pair<TCellID, TCellID>, TCellID> REVERSAL_adj_nodes;
		std::map<TCellID, TCellID> REVERSAL_medium_node;
	};

 public:
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] f_id face id
	 *
	 * \return
	 */
	void setFaceTreated(TCellID f_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] e_id edge id
	 *
	 * \return
	 */
	void setEdgeTreated(TCellID e_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] f_id face id
	 *
	 * \return
	 */
	bool isFaceTreated(TCellID f_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] e_id edge id
	 *
	 * \return
	 */
	bool isEdgeTreated(TCellID e_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] f_id face id
	 * \param[in] n_id face node id
	 * \param[in] ne_id next node id
	 *
	 * \return
	 */
	void setFaceNextNode(TCellID f_id, TCellID n_id, TCellID ne_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] f_id face id
	 * \param[in] n_id next node id
	 *
	 * \return
	 */
	TCellID getFaceNextNode(TCellID f_id, TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] f_id face id
	 * \param[in] n_id next node id
	 *
	 * \return
	 */
	TCellID getFaceIdealNextNode(TCellID f_id, TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] e_id corner edge id
	 * \param[in] n_id node id
	 *
	 * \return
	 */
	void setCornerFaceCreated(TCellID e_id, TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] e_id end edge id
	 * \param[in] n_id node id
	 *
	 * \return
	 */
	void setEndFaceCreated(TCellID e_id, TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] n_id node id
	 *
	 * \return
	 */
	void setReversalFacesCreated(TCellID e_id, TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] n_id node id
	 *
	 * \return
	 */
	bool isCornerFaceCreated(TCellID e_id, TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] n_id node id
	 *
	 * \return
	 */
	bool isEndFaceCreated(TCellID e_id, TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] n_id node id
	 *
	 * \return
	 */
	bool areReversalFacesCreated(TCellID e_id, TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] e_id corner edge id
	 * \param[in] f_id adjacent node on f side
	 * \param[in] n_id edge node id
	 * \param[in] ne_id next node id
	 *
	 * \return
	 */
	void setCornerNextAdjNode(TCellID e_id, TCellID f_id, TCellID n_id, TCellID ne_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] e_id corner edge id
	 * \param[in] f_id adjacent node on f side
	 * \param[in] n_id edge node id
	 *
	 * \return
	 */
	TCellID getCornerNextAdjNode(TCellID e_id, TCellID f_id, TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] e_id corner edge id
	 * \param[in] n_id edge node id
	 * \param[in] ne_id next node id
	 *
	 * \return
	 */
	void setNextDiagNode(TCellID e_id, TCellID n_id, TCellID ne_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] e_id corner edge id
	 * \param[in] n_id edge node id
	 *
	 * \return
	 */
	TCellID getNextDiagNode(TCellID e_id, TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] e_id corner edge id
	 * \param[in] f_id adjacent node on f side
	 * \param[in] n_id edge node id
	 * \param[in] ne_id next node id
	 *
	 * \return
	 */
	void setReversalNextAdjNode(TCellID e_id, TCellID f_id, TCellID n_id, TCellID ne_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] e_id corner edge id
	 * \param[in] f_id adjacent node on f side
	 * \param[in] n_id edge node id
	 *
	 * \return
	 */
	TCellID getReversalNextAdjNode(TCellID e_id, TCellID f_id, TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] e_id corner edge id
	 * \param[in] f_id adjacent node on f side
	 * \param[in] n_id edge node id
	 * \param[in] ne_id next node id
	 *
	 * \return
	 */
	void setReversalNextDiagNode(TCellID e_id, TCellID f_id, TCellID n_id, TCellID ne_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] e_id corner edge id
	 * \param[in] f_id adjacent node on f side
	 * \param[in] n_id edge node id
	 *
	 * \return
	 */
	TCellID getReversalNextDiagNode(TCellID e_id, TCellID f_id, TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] e_id corner edge id
	 * \param[in] n_id edge node id
	 * \param[in] ne_id next node id
	 *
	 * \return
	 */
	void setReversalNextMediumNode(TCellID e_id, TCellID n_id, TCellID ne_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 * \param[in] e_id corner edge id
	 * \param[in] n_id edge node id
	 *
	 * \return
	 */
	TCellID getReversalNextMediumNode(TCellID e_id, TCellID n_id);
	/*-------------------------------------------------------------------*/
 private:
	/*-------------------------------------------------------------------*/
	/** @brief Init the faces struct info for the Front, according to the
	 	* ideal position of each node computed.
	 	*
	 	* \param[in] Front front en entrée
	 	* \param[in] map_new_nodes the map of new nodes
		*
		* \return
	 */
	void InitFaceStructInfo(std::map<TCellID, TCellID> map_new_nodes);
	/*-------------------------------------------------------------------*/
	/** @brief Init the edges struct info for the Front, according to the edges
	 	* classification.
	 	* \param[in] Front front en entrée
		*
		* \return
	 */
	void InitEdgeStructInfo();
	/*-------------------------------------------------------------------*/
 private:
	/** quad mesh to generate */
	Mesh *m_meshH;
	/** */
	Front_3D* m_Front;
	/** Infos sur les noeuds auxquels se connecter pour chaque face du front */
	std::map<TCellID, face_info> m_FaceInfo;
	/** Infos sur les noeuds auxquels se connecter pour chaque edge du front */
	std::map<TCellID, edge_info> m_EdgeInfo;
	/** Mark the treated edges */
	TInt m_mark_edgesTreated;
	/** Mark the treated faces */
	TInt m_mark_facesTreated;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_LAYERSTRUCTUREMANAGER_3D_H
