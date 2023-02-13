/*----------------------------------------------------------------------------*/
#ifndef GMDS_IO_LIMA_READER_H_
#define GMDS_IO_LIMA_READER_H_
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include "GMDSIo_export.h"
/*----------------------------------------------------------------------------*/
#include <Lima/lima++.h>
/*----------------------------------------------------------------------------*/
namespace gmds{

class LimaReader{

 public:
	/*------------------------------------------------------------------------*/
	/** @brief  Constructor.
    *  @param AMesh the mesh in which we want to copy the content of a Lima
     *  	   file.
	 */
	LimaReader(Mesh&AMesh);

	/*------------------------------------------------------------------------*/
	/** \brief  Destructor.	*/
	virtual ~LimaReader();

	/*------------------------------------------------------------------------*/
	/** \brief  Set the mesh length unit. It is the conversion factor from meters
	 */
	double getLengthUnit();

	/*------------------------------------------------------------------------*/
	/** \brief  Read the content of the file named outputName_ and write it in
     *   		mesh_.
	 */
	void read(const std::string &AFileName, MeshModel AModel);

 private:
	void readNodes(Lima::Maillage &ALimaMesh);
	void readEdges(Lima::Maillage &ALimaMesh);
	void readFaces(Lima::Maillage &ALimaMesh);
	void readRegions(Lima::Maillage &ALimaMesh);

	//the input mesh we work with
	Mesh& m_mesh;
	/* connection between original nodes ID and gmds nodes */
	std::vector<TCellID> m_lima2gmds_node_ids;

	/* length unit */
	double m_lengthUnit;
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_IO_LIMA_READER_H
/*----------------------------------------------------------------------------*/