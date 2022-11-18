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
	Front_3D();
	/*-------------------------------------------------------------------*/

 public:
	/*-------------------------------------------------------------------*/
	/** @brief Change l'id du front
	 * 	@param[in] layer_id nouvel id du front
	 */
	void setFrontID(int layer_id);
	/*-------------------------------------------------------------------*/
	/** @brief Retourne l'id du front
	 */
	int getFrontID();
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
