//
// Created by rochec on 18/11/2022.
//

#ifndef GMDS_FRONT_3D_H
#define GMDS_FRONT_3D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_AERO_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/aero/FastLocalize.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
/** \class Front_3D
     *  \brief A Front 3D is a set of adjacent faces and nodes. This front is
     *  closed.
 */
class LIB_GMDS_AERO_API Front_3D {
 public:
	/*--------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS
	} STATUS;

	/*-------------------------------------------------------------------*/
	/** \brief Default constructor
	 */
	Front_3D(int front_id, std::vector<TCellID>& nodes_Id, std::vector<TCellID>& faces_Id);
	/*-------------------------------------------------------------------*/

 public:
	/*-------------------------------------------------------------------*/
	/** @brief Set the front id
	 * 	@param[in] layer_id nouvel id du front
	 */
	void setFrontID(int layer_id);
	/*-------------------------------------------------------------------*/
	/** @brief Give the front id
	 */
	int getFrontID();
	/*-------------------------------------------------------------------*/
	/** @brief Returns the front nodes id
	 */
	std::vector<TCellID> getNodes();
	/*-------------------------------------------------------------------*/
	/** @brief Returns the front faces id
	 */
	std::vector<TCellID> getFaces();
	/*-------------------------------------------------------------------*/
	/** @brief Returns the front edges id
	 */
	std::vector<TCellID> getEdges();
	/*-------------------------------------------------------------------*/
	/** @brief Add the node id to the front
	 * 	@param n_id id of the node to add
	 */
	void addNodeId(TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief Add the face id to the front
	 * 	@param f_id id of the face to add
	 */
	void addFaceId(TCellID f_id);
	/*-------------------------------------------------------------------*/
	/** @brief Returns the ordered front edges around a node of the front.
	 * 	@param m	the mesh
	 *		@param n_id id of the node
	 *
	 *		\return  a vector of the edges id ordered
	 */
	std::vector<TCellID> orderedFrontEdgesAroundNode(Mesh *m, TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief Returns the two front faces around an edge.
	 * 	@param m	the mesh
	 *		@param e_id id of the edge
	 *
	 *		\return  a vector of the two faces id
	 */
	std::vector<TCellID> edgeFacesOnFront(Mesh *m, TCellID e_id);
	/*-------------------------------------------------------------------*/
	/** @brief Returns the outgoing normal to a face on the front.
	 * 	@param mesh_H	the hex mesh that knows the front
	 * 	@param mesh_T	the tet mesh used as a support in our algorithm
	 * 	@param fl the kd tree of the mesh_T
	 *		@param f_id id of the face
	 *
	 *		\return  a 3d vector
	 */
	math::Vector3d outgoingNormal(Mesh *mesh_H, Mesh *mesh_T, FastLocalize* fl, Variable<math::Vector3d>* A_VectorField, TCellID f_id);
	/*-------------------------------------------------------------------*/
	/** @brief Returns the face of the front adjacent to the face f_id
	 * 	and the edge e_id.
	 * 	@param m	the mesh
	 *		@param f_id id of the face
	 *		@param e_id id of the edge
	 *
	 *		\return  a 3d vector
	 */
	TCellID adjacentFaceOnFront(Mesh *m, TCellID f_id, TCellID e_id);
	/*-------------------------------------------------------------------*/
	/** @brief Compute the edges of the front and store them in m_edgesId.
	 * 	@param m	the mesh
	 *
	 *		\return
	 */
	void ComputeEdgesOnFront(Mesh *m);
	/*-------------------------------------------------------------------*/

 private:
	/** Id du front, de la couche */
	int m_FrontID;
	/** Liste d'id des noeuds du front */
	std::vector<TCellID> m_nodesId;
	/** Liste d'id des faces du front */
	std::vector<TCellID> m_facesId;
	/** Liste d'id des arêtes du front */
	std::vector<TCellID> m_edgesId;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_FRONT_3D_H
