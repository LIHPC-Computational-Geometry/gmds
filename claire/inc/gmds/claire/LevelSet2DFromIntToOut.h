//
// Created by rochec on 19/01/2022.
//

#ifndef GMDS_LEVELSET2DFROMINTTOOUT_H
#define GMDS_LEVELSET2DFROMINTTOOUT_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include "gmds/ig/Mesh.h"
#include "DistanceMap.h"
/*----------------------------------------------------------------------------*/
#include <string>
#include <map>
/*----------------------------------------------------------------------------*/
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API LevelSet2DFromIntToOut {
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
	LevelSet2DFromIntToOut(Mesh *AMesh, int AmarkFrontNodesInt, int AmarkFrontNodesOut);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();

 private:
	/** mesh we work on */
	Mesh *m_mesh;
	/** ids of the nodes on the interior front to advance */
	int m_markFrontNodesInt;
	/** ids of the nodes on the exterior front to advance */
	int m_markFrontNodesOut;
	/** carte des distances par rapport au front concerné */
	Variable<double>* m_distance;


};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_LEVELSET2DFROMINTTOOUT_H
