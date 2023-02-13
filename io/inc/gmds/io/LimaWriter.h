/*----------------------------------------------------------------------------*/
#ifndef GMDS_LIMA_WRITER_API_H
#define GMDS_LIMA_WRITER_API_H
/*----------------------------------------------------------------------------*/
#include <map>
#include <sstream>
/*----------------------------------------------------------------------------*/
// headers of GMDS files
#include <gmds/ig/Mesh.h>
#include "GMDSIo_export.h"
/*----------------------------------------------------------------------------*/
#include <Lima/malipp2.h>
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

	/*------------------------------------------------------------------------*/
	/** \brief  Activate the zlib compression.
	 */
	void activateZlibCompression();

 private:
	void writeNodes();
	void writeEdges();
	void writeFaces();
	void writeRegions();

	void writeClouds();
	void writeLines();
	void writeSurfaces();
	void writeVolumes();

	void writeNodesAttributes();
	void writeEdgesAttributes();
	void writeFacesAttributes();
	void writeRegionsAttributes();

	void writeCloudsAttributes();
	void writeLinesAttributes();
	void writeSurfacesAttributes();
	void writeVolumesAttributes();

	/**@brief Check if the node ids are continuously ordered and
	 * initialized the MaliPPWriter2 accordingly.
	 */
	void checkContinuousNodes();
	/**@brief Check if the edge ids are continuously ordered and
* initialized the MaliPPWriter2 accordingly.
	 */
	void checkContinuousEdges();
	/**@brief Check if the face ids are continuously ordered and
* initialized the MaliPPWriter2 accordingly.
	 */
	void checkContinuousFaces();
	/**@brief Check if the region ids are continuously ordered and
	 * initialized the MaliPPWriter2 accordingly.
	 */
	void checkContinuousRegions();
 private:
	/* a mesh */
	gmds::Mesh &m_mesh;

	/* length unit */
	double m_length_unit;

	Lima::MaliPPWriter2 *m_writer;
};
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_LIMAWRITER_API_H