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
	/** \brief  Retire les noeuds qui ne sont connectés à rien dans le maillage
         *
         * \param[in] AMesh the mesh
         *
         * \return  Nothing
	 */
	static void MeshCleaner(Mesh *AMesh);
	/*------------------------------------------------------------------------*/
	/** @brief Donne le vecteur des noeuds adjacents à n dans le maillage m.
	 */
	static std::vector<Node> AdjacentNodes(Mesh* m, Node n);
	/*----------------------------------------------------------------------------*/

};
/*----------------------------------------------------------------------------*/
}  // namespace math
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/

#endif     // GMDS_UTILS_H
