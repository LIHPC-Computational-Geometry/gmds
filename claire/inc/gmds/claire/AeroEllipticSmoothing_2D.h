//
// Created by rochec on 23/08/2022.
//

#ifndef GMDS_AEROELLIPTICSMOOTHING_2D_H
#define GMDS_AEROELLIPTICSMOOTHING_2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/cadfac/FACManager.h>
#include <string>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API AeroEllipticSmoothing_2D
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
	 */
	AeroEllipticSmoothing_2D(Mesh *AMesh, Variable<int>* Alayer, cad::FACManager* Amanager, cad::GeomMeshLinker* Alinker);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Boundary smoothing algorithm
	 */
	void BoundarySlipping();
	/*-------------------------------------------------------------------*/

 private:
	/** mesh we work on */
	Mesh *m_mesh;
	/** layer ids */
	Variable<int>* m_layer;
	/** Manager */
	cad::FACManager* m_manager;
	/** linker to the geometry */
	cad::GeomMeshLinker* m_linker;
	/** slipping nodes */
	std::vector<TCellID> m_slippingNodesId;
};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_AEROELLIPTICSMOOTHING_2D_H
