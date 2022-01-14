//
// Created by rochec on 14/01/2022.
//

#ifndef GMDS_DISTANCEMAP_H
#define GMDS_DISTANCEMAP_H

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
class LIB_GMDS_CLAIRE_API DistanceMap {
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

	DistanceMap();

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Initialize the distances field
	 */

	/*-------------------------------------------------------------------*/

 private:
	/** Tas des couples (distance provisoire, liste d'ids) */
	std::map<double, std::vector<TCellID>> m_map;


};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_DISTANCEMAP_H
