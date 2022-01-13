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
	LevelSet2D(Mesh *AMesh);
	//LevelSet2D(Mesh *AMesh,const Variable<int> *AFrontIds);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();

 private:
	/** mesh we work on */
	Mesh *m_mesh;
	/** carte des distances par rapport au front concerné */
	Variable<double>* m_distance;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds


#endif     // GMDS_LEVELSET2D_H
