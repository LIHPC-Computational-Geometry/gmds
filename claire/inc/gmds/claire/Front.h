//
// Created by rochec on 20/04/2022.
//

#ifndef GMDS_FRONT_H
#define GMDS_FRONT_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <iostream>
#include <string>
#include <map>
#include <fstream>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API Front {
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
	Front();
	/*-------------------------------------------------------------------*/

 public:
	/*-------------------------------------------------------------------*/
	/** @brief Retourne les id des noeuds du front dans un vecteur
	 */
	std::vector<TCellID> getNodes();
	/*-------------------------------------------------------------------*/
	/** @brief Ajoute l'id du noeud au front
	 * 	@param n_id id du noeud à ajouter
	 */
	void addNodeId(TCellID n_id);
	/*-------------------------------------------------------------------*/
 private:
	/** Liste d'id des noeuds du front */
	std::vector<TCellID> m_nodesId;
	/** Tas des couples (node id, type du noeud) */
	std::map<std::vector<TCellID>, int> m_NodeType;


};
/*----------------------------------------------------------------------------*/
}     // namespace gmds


#endif     // GMDS_FRONT_H
