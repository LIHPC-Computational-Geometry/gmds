//
// Created by rochec on 19/05/2022.
//

#ifndef GMDS_SU2WRITER_H
#define GMDS_SU2WRITER_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include "gmds/ig/Mesh.h"
#include <gmds/claire/AeroException.h>
#include <gmds/claire/Front.h>
#include <gmds/claire/Params.h>
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API SU2Writer
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
	SU2Writer(Mesh *AMeshT, std::string AFileName);

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

#endif     // GMDS_SU2WRITER_H
