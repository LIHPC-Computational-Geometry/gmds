//
// Created by Paul Bourmaud on 12/04/2022.
//

#ifndef GMDS_Tools_H
#define GMDS_Tools_H
/*----------------------------------------------------------------------------*/
#include "gmds/ig/Mesh.h"
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

	void suppElementVector(const std::vector<int> vectorInt);

	gmds::GridBuilderAround g_grid;

};
}

#endif     // GMDS_Tools_H
