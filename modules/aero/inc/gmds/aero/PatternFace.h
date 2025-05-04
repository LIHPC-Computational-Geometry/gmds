//
// Created by rochec on 28/04/23.
//

#ifndef GMDS_PATTERNFACE_H
#define GMDS_PATTERNFACE_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_AERO_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/aero/Front_3D.h>
#include <gmds/aero/LayerStructureManager_3D.h>
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_AERO_API PatternFace
{
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
    *  @param AFront the front we build a layer on
    *  @param Af_id face id of the front on which the pattern is applied
    *  @param AStructManager manager of the structure of the layer
	 */
	PatternFace(Mesh *AMesh, Front_3D *AFront, TCellID Af_id, LayerStructureManager_3D *AStructManager);
	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/
 private:
	/*-------------------------------------------------------------------*/
	/** @brief
	 */
	void computeNewHex();
	/*-------------------------------------------------------------------*/
 private:
	/** mesh we work on */
	Mesh* m_mesh;
	/** front we want to build a layer on */
	Front_3D* m_Front;
	/** face id of the front where the pattern is applied */
	TCellID m_f_id;
	/** structure manager for the patterns */
	LayerStructureManager_3D* m_StructManager;
	/** the new hex created */
	std::vector<Region> m_hex;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_PATTERNFACE_H
