//
// Created by calderans on 29/06/23.
//

#ifndef GMDS_CGNSWRITER3D_H
#define GMDS_CGNSWRITER3D_H
/*----------------------------------------------------------------------------*/
#include <sstream>
/*----------------------------------------------------------------------------*/
#include <cgnslib.h>
/*----------------------------------------------------------------------------*/
#include "Blocking3D.h"
#include "GMDSAero_export.h"
// #include "gmds/ig/Blocking2D.h"
/*----------------------------------------------------------------------------*/
namespace gmds {

namespace aero {

class GMDSAero_API CGNSWriter3D
{
 public:
	/** @brief Constructor
		 *
		 * @param AMeshService an implementation of an io service to write data
		 * 						  into a mesh
	 */
	CGNSWriter3D(Blocking3D *ABlocking);

	CGNSWriter3D(Mesh *AMesh);

	CGNSWriter3D();

	/*------------------------------------------------------------------------*/
	/** \brief Destructor. */
	virtual ~CGNSWriter3D();

	void write(const std::string &AInFileName, const std::string &AOutFileName, const std::string &AWorkingDir);

	void writeBoundaryCondition(int &id_bc, cgsize_t *pts, int id_zone, char ABCtype[32], int AEdgeID) const;


 protected:
	void initialize(const std::string &AOutFileName, const std::string &dir);

	void writeZones();

	void writeTri();

	void finalize(const std::string &AWorkingDir) const;

	void _getIndicesIdAndVal(const int *ipnts1, const int *ipnts2, bool *filtre, int &ind, int &val);



 protected:
	Mesh *m_blocks;
	Mesh *m_mesh;

	Variable<std::vector<TCellID>>* m_block_grid;
	Variable<int>* VarDiscrI;
	Variable<int>* VarDiscrJ;
	Variable<int>* VarDiscrK;

	Variable<int>* axis;

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
#endif     // GMDS_CGNSWRITER3D_H
