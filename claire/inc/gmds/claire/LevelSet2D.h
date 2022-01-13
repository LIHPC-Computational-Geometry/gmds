//
// Created by rochec on 13/01/2022.
//

#ifndef GMDS_LEVELSET2D_H
#define GMDS_LEVELSET2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
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

	// OPTION 1
	LevelSet2D(Mesh *AMesh, std::vector<TCellID> Afront_nodes_Ids);
	// OPTION 2
	//LevelSet2D(Mesh *AMesh, Variable<int> *Afront_nodes_Ids);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();

 private:
	/** mesh we work on */
	Mesh *m_mesh;
	/** ids of the nodes on the front to advance */
	// OPTION 1
	std::vector<TCellID> m_front_nodes_Ids;
	// OPTION 2
	//Variable<int>* m_front_nodes_Ids;
	/** carte des distances par rapport au front concerné */
	Variable<double>* m_distance;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds


#endif     // GMDS_LEVELSET2D_H
