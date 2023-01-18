//
// Created by rochec on 21/03/2022.
//

#ifndef GMDS_UTILS_H
#define GMDS_UTILS_H
/*----------------------------------------------------------------------------*/
// gmds file headers
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include "gmds/ig/Mesh.h"
#include "gmds/ig/Blocking2D.h"
#include <gmds/claire/FastLocalize.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace math {
/*----------------------------------------------------------------------------*/
/** \class Utils
 *  \brief
 **/
class LIB_GMDS_CLAIRE_API Utils {
 public:
	/*------------------------------------------------------------------------*/
	/** \brief  Compute the distance between two nodes given the ids
         *
         * \param[in] AMesh the mesh
         * \param[in] n0_id first node id
         * \param[in] n1_id second node id
         *
         * \return  the distance between the nodes n0 and n1
	 */
	static double distFromNodeIds(Mesh *AMesh, TCellID n0_id, TCellID n1_id);

	/*------------------------------------------------------------------------*/
	/** \brief  Return the common edge between 2 points if it exists, NullID
	 		* otherwise (ou si la connectivité n'est pas renseignée)
         *
         * \param[in] n0_id first node id
         * \param[in] n1_id second node id
         *
         * \return  the edge between the nodes n0 and n1
	 */
	static TCellID CommonEdge(Mesh *AMesh, TCellID n0_id, TCellID n1_id);

	/*------------------------------------------------------------------------*/
	/** \brief  Return the common face between 4 points if it exists, NullID
	 		* otherwise (ou si la connectivité n'est pas renseignée)
         *
         * \param[in] n0_id first node id
         * \param[in] n1_id second node id
         * \param[in] n2_id thirst node id
			* \param[in] n3_id fourth node id
         *
         * \return  the face between the nodes n0, n1, n2 and n3
	 */
	static TCellID CommonFace(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id, TCellID n3_id);

	/*------------------------------------------------------------------------*/
	/** \brief  Return the common face between 3 points if it exists, NullID
	 		* otherwise (ou si la connectivité n'est pas renseignée)
         *
         * \param[in] n0_id first node id
         * \param[in] n1_id second node id
         * \param[in] n2_id thirst node id
         *
         * \return  the face between the nodes n0, n1 and n2
	 */
	static TCellID CommonFace3Nodes(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id);

	/*------------------------------------------------------------------------*/
	/** \brief  Retire les noeuds qui ne sont connectés à rien dans le maillage
         *
         * \param[in] AMesh the mesh
         *
         * \return  Nothing
	 */
	static void MeshCleaner(Mesh *AMesh);
	/*------------------------------------------------------------------------*/
	/** @brief Donne le vecteur des noeuds adjacents à n dans le maillage m.
	 	*
		* \param[in] m the mesh
	 	* \param[in] n the node
		*
		* \return  a vector of Nodes
	 */
	static std::vector<Node> AdjacentNodes(Mesh* m, Node n);
	/*----------------------------------------------------------------------------*/
	/** @brief Analyse la qualité d'un maillage composé de quad.
	 	*
		* \param[in] m the mesh
		*
		* \return  nothing
	 */
	static void AnalyseQuadMeshQuality(Mesh* m);
	/*----------------------------------------------------------------------------*/
	/** @brief Build a mesh 2D from a Blocking2D.
	 	*
		* \param[in] blocking2D the blocking
		*
		* \return  the mesh
	 */
	static void BuildMesh2DFromBlocking2D(Blocking2D* blocking2D, Mesh* m, const int mark_block_nodes, const int mark_first_layer_nodes, const int mark_farfield_nodes);
	/*----------------------------------------------------------------------------*/
	/** @brief Return the point at position alpha of the branch. alpha = 0.5 returns
	 * the mid point on the branch.
	 	*
		* \param[in] A, B, C the three ordered points of the branch
		*
		* \return  the point at position alpha from the point A
	 */
	static math::Point WeightedPointOnBranch(const math::Point A, const math::Point B, const math::Point C, double alpha);
	/*----------------------------------------------------------------------------*/
	/** @brief Return true if the point M is in the triangle T
	 	*
		* \param[in] T1, T2, T3 the point of the triangle
	 	* \param[in] M, the point to check if it's in the triangle
		*
		* \return  the point at position alpha from the point A
	 */
	static bool isInTriangle(const math::Point T1, const math::Point T2, const math::Point T3, const math::Point M);
	/*----------------------------------------------------------------------------*/
	/** @brief Return true if the points M1 and M2 are on the same side of the plan
	 	* defined by the points T1/T2/T3
	 	*
		* \param[in] T1, T2, T3 the point of the plan
	 	* \param[in] M1, M2 the points to check the side of
		*
		* \return  true if M1 and M2 are on the same side, false otherwise
	 */
	static bool sameSide(const math::Point T1, const math::Point T2, const math::Point T3, const math::Point M1, const math::Point M2);
	/*----------------------------------------------------------------------------*/
	/** @brief Return true if the point M is in the tetra T
	 	*
		* \param[in] T1, T2, T3, T4 the point of the tetra
	 	* \param[in] M the point to check if it's in the tetra
		*
		* \return  true if M is in the tetra, false otherwise
	 */
	static bool isInTetra(const math::Point T1, const math::Point T2, const math::Point T3, const math::Point T4, const math::Point M);
	/*----------------------------------------------------------------------------*/
	/** @brief Linear 2D interpolation with 3 points
	 	*
		* \param[in] P1, P2, P3 the 3 points were the value is known
	 	* \param[in] c1, c2, c3 the 3 values à the points P1, P2 and P3
	 	* \param[in] M the position we want to interpolate the value
		*
		* \return  the interpolated value
	 */
	static double linearInterpolation2D3Pt(const math::Point P1, const math::Point P2, const math::Point P3, const math::Point M, const double c1, const double c2, const double c3);
	/*----------------------------------------------------------------------------*/
	/** @brief Reavel the curved block edges
	 	*
		* \param[in] blocking2D the blocking
		*
		* \return  the mesh
	 */
	static void CurveBlockEdgesReavel(Blocking2D* blocking2D, Mesh* m);
	/*----------------------------------------------------------------------------*/
	/** \brief  Create the edge defined by the two nodes to the mesh, and create the
	 		* connectivities N<->E.
         *
         * \param[in] AMesh the mesh
         * \param[in] n0_id first node id
         * \param[in] n1_id second node id
         *
         * \return  the id of the new face
	 */
	static TCellID GetOrCreateEdgeAndConnectivitiesN2E(Mesh *AMesh, TCellID n0_id, TCellID n1_id);
	/*----------------------------------------------------------------------------*/
	/** \brief  Create the quad defined by the four nodes to the mesh and create the
	 		* four edges. Add all the connectivities N<->F, N<->E and E<->F.
         *
         * \param[in] AMesh the mesh
         * \param[in] n0_id first node id
         * \param[in] n1_id second node id
         * \param[in] n2_id thirst node id
			* \param[in] n3_id fourth node id
         *
         * \return  the id of the new face
	 */
	static TCellID GetOrCreateQuadAndConnectivities(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id, TCellID n3_id);
	/*----------------------------------------------------------------------------*/
	/** @brief Create the hexa from the 8 nodes, and all the connectivities
	 	*
	 	* \param[in] AMesh the mesh
  		* \param[in] n0 first node
  		* \param[in] n1 second node
  		* \param[in] n2 thirst node
  		* \param[in] n3 fourth node
  		* \param[in] n4
  		* \param[in] n5
  		* \param[in] n6
  		* \param[in] n7 last node
		*
		* \return the id of the region (hexa) created
	 */
	static TCellID CreateHexaNConnectivities(Mesh *Amesh, Node n0, Node n1, Node n2, Node n3, Node n4, Node n5, Node n6, Node n7);
	/*-------------------------------------------------------------------*/
	/** @brief Check if the hex mesh is valid
	 	*
	 	* \param[in] AMesh the mesh
		*
		* \return true if the mesh is valid, based on the values checked
	 */
	static bool isThisHexMeshValid(Mesh *Amesh);
	/*-------------------------------------------------------------------*/
	/** @brief Compute the position of a point advected along a constant
	 	* vector field, using the AdvectedPointRK4 algorithm.
	 	*
	 	* \param[in] AMesh the mesh
	 	* \param[in] fl tree of the mesh AMesh
	 	* \param[in] M strating point
	 	* \param[in] dist_cible target distance in the distance field
	 	* \param[in] A_distance the distance field
	 	* \param[in] A_vector the vector
		*
		* \return the point of the final position of the point M
	 */
	static math::Point AdvectedPointRK4_UniqVector_3D(Mesh *Amesh, FastLocalize *fl, math::Point M, double dist_cible, Variable<double>* A_distance, math::Vector3d A_vector);
	/*-------------------------------------------------------------------*/
	/** @brief Return the 2 faces adjacents to the edge e_id in the hexa
	 	*
	 	* \param[in] AMesh the mesh
	 	* \param[in] e_id the edge id
	 	* \param[in] r_id the region id
		*
		* \return the two faces adjacents to the edge e_id in the hexa
	 */
	static std::vector<Face> getFacesAdjToEdgeInHexa(Mesh *Amesh, TCellID e_id, TCellID r_id);
	/*-------------------------------------------------------------------*/
	/** @brief Return the opposite edge of e in a quad face
	 	*
	 	* \param[in] AMesh the mesh
	 	* \param[in] e_id the edge id
	 	* \param[in] f_id the face id
		*
		* \return the opposite edge of e in a quad face
	 */
	static Edge oppositeEdgeInFace(Mesh *Amesh, TCellID e_id, TCellID f_id);
	/*-------------------------------------------------------------------*/
	/** @brief Compute the minimal edge's lenght
	 	*
	 	* \param[in] AMesh the mesh
		*
		* \return
	 */
	static double minEdgeLenght(Mesh *Amesh);
	/*-------------------------------------------------------------------*/

};
/*----------------------------------------------------------------------------*/
}  // namespace math
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/

#endif     // GMDS_UTILS_H
