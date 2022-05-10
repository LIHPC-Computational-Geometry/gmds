//
// Created by Paul Bourmaud on 12/04/2022.
//

#ifndef GMDS_Tools_H
#define GMDS_Tools_H
/*----------------------------------------------------------------------------*/
#include "gmds/ig/Mesh.h"
#include "gmds/ig/MeshDoctor.h"
#include "GMDSIgAlgo_export.h"
#include "gmds/paul/Grid.h"

/*----------------------------------------------------------------------------*/


namespace gmds{
//Tools
class Tools{
 public:
	Tools(GridBuilderAround *AGrid);

	/** \brief  function for get the value of activate, if the face was delete or note
  *
  *  \return an int for the value of activate (false = 0 and true=1)
	 */
	int getValueActivateFace(const int faceIDChecked);

	/** \brief  Information about a edge "existing" between 2 nodes
	 *
	 * \param i1 id first node
	 * \param i2 id second node
	 * \param faceID id face Common
  *  \return an bool
  */
	bool checkExistEdge(const int i1, const int i2,const int faceID);

	/** \brief  Check if a face existing between 2 nodes
	 *
	 * \param i1 id first node
	 * \param i2 id second node
  *
  *  \return an bool
	 */
	bool checkCommonFace(const int i1, const int i2);


	/** \brief  Verifie that 2 nodes follow each other
	 *
	 * \param i1 id first node
	 * \param i2 id second node
	 * \param faceID id face Common
  *  \return an bool
	 */

	bool checkFollowIdNode(const int i1, const int i2,const int faceID);


	/** \brief  get list of nodes for a specify face
	 *
	 * \param faceID id face
	 *
		*\return an vector of nodes
	 */

	std::vector<Node> getListNodesOfFace(const int faceID);

	/** \brief  get list of faces for a specify node
	 *
	 * \param nodeID id node
	 *
		*\return an vector of faces
	 */

	std::vector<Face> getListFacesOfNode(const int nodeID);

	/** \brief  get face(s) common between 2 nodes
	 *
	 * \param i1 id first node
	 * \param i2 id second node
	 *
		*\return an vector of faces
	 */

	std::vector<Face> getFacesCommon(const int i1, const int i2);

	/** \brief  get only one face common between 2 nodes
	 *
	 * \param i1 id first node
	 * \param i2 id second node
	 *
		*\return an int (the id of face)
	 */

	int getIdOneCommonFace(const int i1, const int i2);

	/** \brief  get previous node id in the list of nodes for a face
	 *
	 * \param idNode id node
	 * \param idFaceNode id face
	 *
		*\return an int (the id of the node)
	 */

	int getIdPreviousNode(const int idNode, const int idFaceNode);

	/** \brief  get next node id in the list of nodes for a face
	 *
	 * \param idNode id node
	 * \param idFaceNode id face
	 *
		*\return an int (the id of the node)
	 */

	int getIdNextNode(const int idNode, const int idFaceNode);

	/** \brief  get the opposite edge(s) of a edge
	 *
	 * \param i1 id first node
	 * \param i2 id second node
	 *
		*\return an vector of vector of int (vec<vec<int>>)
	 */

	std::vector<std::vector<int>> getOtherNodes(const int i1, const int i2);

	/** \brief  get all the faces from a edge (simulated by the id of 2 nodes)
	 *
	 * \param i1 id first node
	 * \param i2 id second node
	 *
		*\return an vector of face (vec<face>)
	 */

	std::vector<Face> getAllFacesChain(const int i1, const int i2);

	/** \brief  get all the edge (pair of id of node) from a edge (simulated by the id of 2 nodes)
	 *
	 * \param i1 id first node
	 * \param i2 id second node
	 *
		*\return an vector of vector of int (vec<vec<int>>)
	 */

	std::vector<std::vector<int>> getAllNodesChain(const int i1, const int i2);

	//en cours d'implé, par sur  necessaire pour le moment
	void suppElementVector(const std::vector<int> vectorInt);

	/** \brief get all nodes with a common edge for a specifiy node
	 *
	 * @param i1 id node
	 * @return an vector of nodes
	 */
	std::vector<Node> nodesAroundANode(const int i1);

	/** \brief get opposite nodes
	 *
	 * @param i1 id first node (node target)
	 * @param i2 id second node (node to know the edge)
	 * @return an vector of nodes
	 */
	std::vector<Node> getOppositeNodes(const int i1,const int i2);

	/** \brief get list nodes on the same line of the first element
	 *
	 * @param i1 id first node (node target)
	 * @param i2 id second node (other node to know the edge)
	 * @return vector int
	 */
	std::vector<int> getListFirstNodesChain(const int i1, const int i2);

	/** \brief get list nodes on the same line of the second element
	 *
	 * @param i1 id first node (other node to know the edge)
	 * @param i2 id second node (node target)
	 * @return vector int
	 */

	std::vector<int> getListSecondNodesChain(const int i1, const int i2);

	/** \brief create middle point between 2 points
	 *
	 * @param i1 first point id
	 * @param i2  second point id
	 * @return
	 */
	Node createMiddleNode(Node i1,Node i2);

	/** @brief function to create au new middle ligne of nodes
	 *
	 * @param i1 first point id
	 * @param i2 second point id
	 */
	void createAllMiddlePoint(Node i1, Node i2);

	/** @brief get a pair of edges for a specify face and a edge reference (for know the direction of the edges)
	 *
	 * @param i1 first node of the edge reference
	 * @param i2 second node of  the edge reference
	 * @param faceID face target id
	 * @return two edges (vector of vector of int, with only 2 int by vector for the pair)
	 */
	std::vector<std::vector<int>> getPairNodesFace(const int i1, const int i2, const int faceID);

	/** @brief return one boundary edge of a chain faces
	 *
	 * @param firstNodeId
	 * @param secondNodeId
	 * @return vector of nodes (the edge)
	 */
	std::vector<Node> getBoundaryEdge(Node firstNodeId, Node secondNodeId);

	//premier test de la fonction pour couper les aretes
	void joinFaceToNodes(Node i1,Node i2);

	/** @brief create new edges and new faces during the action executeCutEdge
	 *
	 * @param i1 first node of the edge to cut
	 * @param i2 second node of the edge to cut
	 */
	void createEdge(Node i1, Node i2);

	/** @brief get the face common of 2 edges
	 *
	 * @param i1 first node first edge
	 * @param i2 second node first edge
	 * @param i3 first node second edge
	 * @param i4 second node second edge
	 * @return face
	 */
	Face getFaceEdges(const int i1,const int i2,const int i3,const int i4);

	/** @brief get all nodes in the same way, use for checkVertical function
	 *
	 * @param i1 first node
	 * @param i2 second node
	 * @return vector of int (id nodes)
	 */
	std::vector<int> getAllNodesSameWay(Node i1, Node i2);

	/** @brief check if a edge is vertical
	 *
	 * @param i1 first node
	 * @param i2 second node
	 * @return true if vertical, false if horizontal
	 */
	bool checkVertical(Node i1, Node i2 );

	/** get all boundary nodes for a mesh
	 *
	 * @param AMesh the mesh
	 * @return map nodes
	 */
	std::map<TCellID ,TCellID> getBoundaryNodes (Mesh *AMesh);

	double calcRangePoints(Node node1, Node node2);

	gmds::GridBuilderAround g_grid;

};
}

#endif     // GMDS_Tools_H
