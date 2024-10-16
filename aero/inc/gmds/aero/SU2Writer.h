//
// Created by rochec on 19/05/2022.
//

#ifndef GMDS_SU2WRITER_H
#define GMDS_SU2WRITER_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_AERO_export.h"
#include "gmds/ig/Mesh.h"
#include <gmds/aero/AeroException.h>
#include <gmds/aero/Front.h>
#include <gmds/aero/Params.h>
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_AERO_API SU2Writer
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
         *  @param[in] Ax_lim_inout limit between inlet and outlet
         *
	 */
	SU2Writer(Mesh *AMeshT, std::string AFileName, double Ax_lim_inout);

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
	/** limit between inlet and outlet */
	double m_x_lim_inout;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_SU2WRITER_H
