//
// Created by rochec on 28/04/23.
//

#ifndef GMDS_ABSTRACTPATTERNEDGE_H
#define GMDS_ABSTRACTPATTERNEDGE_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_AERO_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/aero/Front_3D.h>
#include <gmds/aero/LayerStructureManager_3D.h>
#include <gmds/aero/FastLocalize.h>
#include <string>
#include <map>
#include <fstream>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** \class  AbstractPatternEdge
 *  \brief
 */
class LIB_GMDS_AERO_API AbstractPatternEdge{

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
	AbstractPatternEdge(Mesh *AMesh, Front_3D *AFront, TCellID Ae_id, LayerStructureManager_3D *AStructManager,
	                    Mesh *AMeshT, FastLocalize *Afl,
	                    double dc, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField);
	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/
 public:
	/*-------------------------------------------------------------------*/
	/** @brief
	 *
	 * \return
	 */
	std::vector<Region> getNewHex();
	/*-------------------------------------------------------------------*/
 protected:
	/*-------------------------------------------------------------------*/
	/** @brief
	 */
	virtual void computeNewHex()=0;
	/*-------------------------------------------------------------------*/
 protected:
	/** mesh we work on */
	Mesh* m_mesh;
	/** front we want to build a layer on */
	Front_3D* m_Front;
	/** */
	TCellID m_e_id;
	/** structure manager for the patterns */
	LayerStructureManager_3D* m_StructManager;
	/** */
	Mesh* m_meshT;
	/** */
	FastLocalize* m_fl;
	/** Target distance */
	double m_dc;
	/** Distance Field */
	Variable<double>* m_DistanceField;
	/** Vector Field */
	Variable<math::Vector3d>* m_VectorField;
	/** */
	std::vector<Region> m_hex;

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_ABSTRACTPATTERNEDGE_H
