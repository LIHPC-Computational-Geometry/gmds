//
// Created by rochec on 13/01/2022.
//

#ifndef GMDS_LEVELSET2D_H
#define GMDS_LEVELSET2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/claire/DistanceMap.h>
/*----------------------------------------------------------------------------*/
#include <string>
#include <map>
/*----------------------------------------------------------------------------*/
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API LevelSet2D {
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
         *  @param AMesh the mesh where we work on
	 */
	LevelSet2D(Mesh *AMesh, int AmarkFrontNodes);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Initialize the DistanceMap
	 */
	void initialisationDistances();
	/*-------------------------------------------------------------------*/
	/** @brief Défini la distance pour un noeud
	 */
	void setValue(TCellID n_id, double v0);
	/*-------------------------------------------------------------------*/

 public:
	/*-------------------------------------------------------------------*/
	/** @brief Obtient la distance pour un noeud
	 */
	void getValue(TCellID n_id, double &v0);
	/*-------------------------------------------------------------------*/

 private:
	/** mesh we work on */
	Mesh *m_mesh;
	/** ids of the nodes on the front to advance */
	int m_markFrontNodes;
	/** Tas des couple (distance provisoire, liste d'id) */
	DistanceMap m_DistanceMap;
	/** carte des distances par rapport au front concerné */
	Variable<double>* m_distance;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds


#endif     // GMDS_LEVELSET2D_H
