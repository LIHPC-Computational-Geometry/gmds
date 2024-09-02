/*----------------------------------------------------------------------------*/
#ifndef GMDS_LIMA_WRITER_H
#define GMDS_LIMA_WRITER_H
/*----------------------------------------------------------------------------*/
#include <map>
#include <sstream>
/*----------------------------------------------------------------------------*/
// headers of GMDS files
#include <gmds/ig/Mesh.h>
#include "GMDSIo_export.h"
/*----------------------------------------------------------------------------*/
#include <Lima/lima++.h>
/*----------------------------------------------------------------------------*/
namespace gmds {

class GMDSIo_API LimaWriter
{
 public:
	/*------------------------------------------------------------------------*/
	/** \brief  Constructor.
	 *
	 *  \param AMesh the mesh we want to write into a file.
	 */
	LimaWriter(Mesh &AMesh);

	/*------------------------------------------------------------------------*/
	/** \brief  Destructor.	*/
	virtual ~LimaWriter();

	/*------------------------------------------------------------------------*/
	/** \brief  Set the mesh length unit. It is the conversion factor from meters
	 */
	void setLengthUnit(double AUnit);

	/*------------------------------------------------------------------------*/
	/** \brief  Write the content of mesh_ into the file named AFileName.
	 */
	void write(const std::string &AFileName, MeshModel AModel, int ACompact = false);

 private:

	void writeClassic(const std::string &AFileName, MeshModel AModel, int ACompact = false);

	void writeNodes(Lima::Maillage& ALimaMesh);
	void writeEdges(Lima::Maillage& ALimaMesh);
	void writeFaces(Lima::Maillage& ALimaMesh);
	void writeRegions(Lima::Maillage& ALimaMesh);

//	void writeClouds();
//	void writeLines();
//	void writeSurfaces();
//	void writeVolumes();

 private:
	/* a mesh */
	gmds::Mesh &m_mesh;

	/* length unit */
	double m_length_unit;

	/* connection between original nodes ID and lima nodes */
	std::vector<Lima::Noeud> m_nodes_connection;
};
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_LIMA_WRITER_H