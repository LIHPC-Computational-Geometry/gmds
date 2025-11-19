//
// Created by calderans on 29/06/23.
//

#ifndef GMDS_CGNSWRITERND_H
#define GMDS_CGNSWRITERND_H
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

class GMDSAero_API CGNSWriterND
{
 public:
	/** @brief Constructor
		 *
		 * @param AMeshService an implementation of an io service to write data
		 * 						  into a mesh
	 */
	CGNSWriterND(Blocking3D *ABlocking, int ADim);

	CGNSWriterND(Mesh *AMesh, int ADim);

	CGNSWriterND();

	/*------------------------------------------------------------------------*/
	/** \brief Destructor. */
	virtual ~CGNSWriterND();

	void write(const std::string &AInFileName, const std::string &AOutFileName, const std::string &AWorkingDir);

 protected:
	void initialize(const std::string &AOutFileName, const std::string &dir);

	void writeZones();

	void writeTri();

	void finalize(const std::string &AWorkingDir) const;

	void _getIndicesIdAndVal(const int *ipnts1, const int *ipnts2, bool *filtre, int &ind, int &val);

	void writeConnections3D(const Region& Ablock, int iFace, int& index_tf, const std::vector<Variable<int>*>& zone_vars);
	void writeConnections2D(const Face& Ablock, int iEdge, int& index_tf, const std::vector<Variable<int>*>& zone_vars);

	void writeBoundaryCondition3D(int &num_bc, const Region& Ablock, int iFace, const std::vector<Variable<int>*>& bc_vars) const;
	void writeBoundaryCondition2D(int &num_bc, const Face& Ablock, int iEdge, const std::vector<Variable<int>*>& bc_vars) const;

 protected:
	Mesh *m_blocks;
	Mesh *m_mesh;

	Variable<std::vector<TCellID>>* m_block_grid;
	Variable<int>* VarDiscrI;
	Variable<int>* VarDiscrJ;
	Variable<int>* VarDiscrK;

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
#endif     // GMDS_CGNSWRITERND_H
