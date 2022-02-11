//
// Created by rochec on 10/02/2022.
//

#ifndef GMDS_ABSTRACTLEVELSET_H
#define GMDS_ABSTRACTLEVELSET_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/claire/DistanceMap.h>
#include <string>
#include <map>
#include <fstream>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** \class  AbstractLevelSet
 *  \brief  Classe abstraite pour l'algorithme de level set.
 */
class LIB_GMDS_CLAIRE_API AbstractLevelSet{

 public:
	/*--------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS
	} STATUS;
	/*--------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param
	 */
	AbstractLevelSet(Mesh *AMesh, int AmarkFrontNodes, Variable<double>* Adistance);
	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/
	/** @brief Obtient la distance pour un noeud
	 */
	void getValue(TCellID n_id, double &v0);
	/*-------------------------------------------------------------------*/

 protected:
	/*-------------------------------------------------------------------*/
	/** @brief Initialize the DistanceMap
	 */
	void initialisationDistances();
	/*-------------------------------------------------------------------*/
	/** @brief Défini la distance pour un noeud
	 */
	void setValue(TCellID n_id, double v0);
	/*-------------------------------------------------------------------*/
	/** @brief Récupère le voisinage du noeud
	 */
	virtual std::vector<Node> getNeighbors(Node n)=0;
	/*-------------------------------------------------------------------*/

 protected:
	/** mesh we work on */
	Mesh* m_mesh;
	/** ids of the nodes on the front to advance */
	int m_markFrontNodes;
	/** Tas des couple (distance provisoire, liste d'id) */
	DistanceMap m_DistanceMap;
	/** carte des distances par rapport au front concerné */
	Variable<double>* m_distance;


};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_ABSTRACTLEVELSET_H
