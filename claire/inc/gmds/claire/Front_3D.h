//
// Created by rochec on 18/11/2022.
//

#ifndef GMDS_FRONT_3D_H
#define GMDS_FRONT_3D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
/** \class Front_3D
     *  \brief A Front 3D is a set of adjacent faces and nodes. This front is
     *  closed.
 */
class LIB_GMDS_CLAIRE_API Front_3D {
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
	Front_3D(int front_id, std::vector<TCellID> nodes_Id, std::vector<TCellID> faces_Id);
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

 private:
	/** Id du front, de la couche */
	int m_FrontID;
	/** Liste d'id des noeuds du front */
	std::vector<TCellID> m_nodesId;
	/** Liste d'id des arêtes du front */
	std::vector<TCellID> m_facesId;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_FRONT_3D_H
