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
	/** @brief Retourne les id des arêtes du front dans un vecteur
	 */
	std::vector<TCellID> getEdges();
	/*-------------------------------------------------------------------*/
	/** @brief Ajoute l'id du noeud au front
	 * 	@param n_id id du noeud à ajouter
	 */
	void addNodeId(TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief Ajoute l'id de l'arête au front
	 * 	@param e_id id de l'arête à ajouter
	 */
	void addEdgeId(TCellID e_id);
	/*-------------------------------------------------------------------*/
	/** @brief Initialise le Front à partir de l'ID de la couche souhaitée.
	 * 	@param[in] layer_id id de la couche des noeuds
	 */
	void initializeFromLayerId(Mesh *m, int layer_id);
	/*-------------------------------------------------------------------*/
 private:
	/** Liste d'id des noeuds du front */
	std::vector<TCellID> m_nodesId;
	/** Liste d'id des arêtes du front */
	std::vector<TCellID> m_edgesId;
	/** Tas des couples (node id, type du noeud) */
	std::map<std::vector<TCellID>, int> m_NodeType;


};
/*----------------------------------------------------------------------------*/
}     // namespace gmds


#endif     // GMDS_FRONT_H
