//
// Created by rochec on 26/09/23.
//

#ifndef GMDS_MFEMMESHWRITER_H
#define GMDS_MFEMMESHWRITER_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_AERO_export.h"
#include "gmds/ig/Mesh.h"
#include <gmds/cadfac/FACManager.h>
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_AERO_API MFEMMeshWriter
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
         *  @param[in] AMesh the mesh to write
         *  @param[in] AFileName the name of the written file
         *
	 */
	MFEMMeshWriter(Mesh *AMeshT, std::string AFileName);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/
 private:
	/*-------------------------------------------------------------------*/
	/** @brief Ordering nodes in m_ordering_nodes tables
	 */
	void orderingNodes();
	/*-------------------------------------------------------------------*/
	/** @brief Write the header
	 */
	void header();
	/*-------------------------------------------------------------------*/
	/** @brief Write the elements (faces in 2D)
	 */
	void writeElements();
	/*-------------------------------------------------------------------*/
	/** @brief Write the elements (regions in 3D)
	 */
	void writeElements3D();
	/*-------------------------------------------------------------------*/
	/** @brief Write boundary elements (segments in 2D)
	 */
	void writeBnd();
	/*-------------------------------------------------------------------*/
	/** @brief Write boundary elements (faces in 3D)
	 */
	void writeBnd3D();
	/*-------------------------------------------------------------------*/
	/** @brief Write the nodes
	 */
	void writeNodes();
	/*-------------------------------------------------------------------*/
	/** @brief Write the nodes
	 */
	void writeNodes3D();
	/*-------------------------------------------------------------------*/

 private:
	/** the mesh to write */
	Mesh *m_mesh;
	/** the name of the written file */
	std::string m_filename;
	/** */
	std::ofstream m_stream;
	/** map to order nodes from 0 to nbr_nodes-1 */
	std::map<TCellID,int> m_ordering_nodes;
	/** Manager */
	cad::FACManager* m_manager;
	/** Linker input mesh to geometry */
	cad::GeomMeshLinker* m_linker;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_MFEMMESHWRITER_H
