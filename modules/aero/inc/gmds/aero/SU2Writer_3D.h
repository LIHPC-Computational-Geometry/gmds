//
// Created by rochec on 11/01/24.
//

#ifndef GMDS_SU2WRITER_3D_H
#define GMDS_SU2WRITER_3D_H

/*----------------------------------------------------------------------------*/
#include "GMDSAero_export.h"
#include "gmds/ig/Mesh.h"
#include <gmds/aero/Front.h>
#include <gmds/aero/Params.h>
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class GMDSAero_API SU2Writer_3D
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
	SU2Writer_3D(Mesh *AMeshT, std::string AFileName);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/

 private:
	/** the mesh to write */
	Mesh *m_mesh;
	/** the name of the written file */
	std::string m_filename;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_SU2WRITER_3D_H
