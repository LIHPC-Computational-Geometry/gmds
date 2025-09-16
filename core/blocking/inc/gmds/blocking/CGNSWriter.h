//
// Created by calderans on 22/03/2022.
//
/*----------------------------------------------------------------------------*/
#ifndef GMDS_CGNSWRITER_H
#define GMDS_CGNSWRITER_H
/*----------------------------------------------------------------------------*/
#include <sstream>
/*----------------------------------------------------------------------------*/
#include <cgnslib.h>
/*----------------------------------------------------------------------------*/
#include "GMDSBlocking_export.h"
#include "gmds/ig/Blocking2D.h"
/*----------------------------------------------------------------------------*/
namespace gmds {

namespace blocking {

class GMDSBlocking_API CGNSWriter
{
 public:
	/** @brief Constructor
		 *
		 * @param AMeshService an implementation of an io service to write data
		 * 						  into a mesh
	 */
	CGNSWriter(Blocking2D *ABlocking);

	CGNSWriter(Mesh *AMesh);

	/*------------------------------------------------------------------------*/
	/** \brief Destructor. */
	virtual ~CGNSWriter();

	void write(const std::string &AFileName, const std::string &AWorkingDir);

	void writeBoundaryCondition(int &id_bc, cgsize_t *pts, int id_zone, char ABCtype[32], int AEdgeID) const;

 protected:
	void initialize(const std::string &AFileName);

	void writeZones();

	void writeTri();

	void finalize(const std::string &AWorkingDir) const;

	Blocking2D *m_blocks;
	Mesh *m_mesh;

	int m_cellDim;
	int m_physdim;

	int m_indexFile;
	int m_indexBase;
	int m_indexZone;

	std::ofstream *m_stream;
};
/*----------------------------------------------------------------------------*/
} // namespace blocking
/*----------------------------------------------------------------------------*/
}// namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_CGNSWRITER_H
