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

	struct multiple_info{
		TCellID     					n_id;
		TCellID 							ideal_node_id;
		int 								nodeType;
		std::map<TCellID, TCellID> next_nodes;	// (e_id, next_n_id)
	};

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
	/** @brief Retourne les id des noeuds voisins sur le front
	 * 	@param[in] n_id id du noeud concerné
	 */
	std::vector<TCellID> getNeighbors(TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief Retourne le type du noeud d'id n_id
	 * 	@param[in] n_id id du noeud concerné
	 */
	int getNodeType(TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief Retourne le noeud dans la couche suivante pour la
	 * 	construction du quad.
	 * 	@param[in] n_id id du noeud concerné
	 * 	@param[in] n_neighbors_id id du noeud voisin qui fait
	 * 	partie du quad
	 */
	TCellID getNextNode(TCellID n_id, TCellID n_neighbors_id);
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
	 * 	@param[in] m the mesh we work on
	 * 	@param[in] layer_id id de la couche des noeuds
	 */
	void initializeFromLayerId(Mesh *m, int layer_id);
	/*-------------------------------------------------------------------*/
	/** @brief Initialise la map m_NodeInfo avec une map déjà existante.
	 * 	Cette map fait correspondre à chaque noeud du front un noeud.
	 * 	@param[in] m the mesh we work on
	 * 	@param[in] map_idealpositions map contenant les ID des noeuds
	 * 	créés à la position idéale.
	 */
	void initializeNodeType(Mesh *m, std::map<TCellID, TCellID> map_idealpositions);
	/*-------------------------------------------------------------------*/
	/** @brief Initialise la map m_NodeNeighbors des noeuds voisins à un noeud
	 * 	sur un front.
	 * 	@param[in] m the mesh we work on
	 */
	void initializeNodeNeighbors(Mesh *m);
	/*-------------------------------------------------------------------*/
 private:
	/** Liste d'id des noeuds du front */
	std::vector<TCellID> m_nodesId;
	/** Liste d'id des arêtes du front */
	std::vector<TCellID> m_edgesId;
	/** Infos sur les noeuds multiples et contractés */
	std::map<TCellID, multiple_info> m_NodeInfo;
	/** Noeuds voisins sur le front */
	std::map<TCellID, std::vector<TCellID>> m_NodeNeighbors;


};
/*----------------------------------------------------------------------------*/
}     // namespace gmds


#endif     // GMDS_FRONT_H
