//
// Created by rochec on 21/12/2022.
//

#ifndef GMDS_NODENEIGHBOURHOODONFRONT_3D_H
#define GMDS_NODENEIGHBOURHOODONFRONT_3D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/claire/Front_3D.h>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API NodeNeighbourhoodOnFront_3D
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
         *  @param[in] AMesh the triangular mesh where we work on
         *  @param[in] AFront front 3D
         *  @param[in] An_id node id
         *
	 */
	NodeNeighbourhoodOnFront_3D(Mesh *AMesh, Front_3D *AFront, TCellID An_id);

	/*-------------------------------------------------------------------*/
	/** @brief Destructor.
         *  @param
	 */
	~NodeNeighbourhoodOnFront_3D();
	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/
	/** @brief Returns the front edges ordered
	 */
	std::vector<TCellID> getOrderedEdges();
	/*-------------------------------------------------------------------*/
	/** @brief Returns the front faces ordered
	 */
	std::vector<TCellID> getOrderedFaces();
	/*-------------------------------------------------------------------*/
	/** @brief Two adjacent faces to e_id on the front
	 *
	 * 	\return a vector sized two, with the 2 ids of the adjacent faces
	 * 	on the front
	 */
	std::vector<TCellID> adjFacesToEdge(TCellID e_id);
	/*-------------------------------------------------------------------*/
	/** @brief The next edge of the face f around the node m_n
	 *
	 * 	\return
	 */
	TCellID nextEdgeOfFace(TCellID f_id, TCellID e_id);
	/*-------------------------------------------------------------------*/
	/** @brief Adjacent face to Edge_1 on the front, in the side of Edge_2
	 *
	 * 	\return a face id
	 */
	TCellID adjFaceToEdge1InEdge2SideAvoidingEdge3(TCellID e1_id, TCellID e2_id, TCellID e3_id);
	/*-------------------------------------------------------------------*/
 private:
	/** @brief Returns the ordered front edges around a node of the front.
	 *
	 *		\return  fill the vector of the edges id ordered m_orderedEdges
	 */
	void orderedFrontEdgesAroundNode();
	/*-------------------------------------------------------------------*/
	/** @brief Returns the ordered front edges around a node of the front.
	 *
	 *		\return  fill the vector of the faces id ordered m_orderedFaces
	 */
	void orderedFrontFacesAroundNode();
	/*-------------------------------------------------------------------*/
 private:
	/** the quad mesh we work on */
	Mesh *m_mesh;
	/** the front */
	Front_3D *m_Front;
	/** Node */
	TCellID m_n_id;
	/** Ordered edges around the node */
	std::vector<TCellID> m_orderedEdges;
	/** Ordered faces around the node */
	std::vector<TCellID> m_orderedFaces;
};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_NODENEIGHBOURHOODONFRONT_3D_H
