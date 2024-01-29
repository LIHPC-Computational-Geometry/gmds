
#ifndef GMDS_ELLIPTICMORPH_H
#define GMDS_ELLIPTICMORPH_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_MORPHMESH_export.h"
#include "gmds/ig/Mesh.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace morphmesh {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_MORPHMESH_API EllipticMorph
{

 public:
	/*------------------------------------------------------------------------*/
	/** @brief Constructor.
	 */
	EllipticMorph(const std::string& AFilename, Mesh* AMesh);

	/*------------------------------------------------------------------------*/
	/** @brief Destructor.
	 */
	~EllipticMorph();

	/*------------------------------------------------------------------------*/
	/** @brief Constructor.
	 */
	void execute();

 private:

	/*Surement a retirer plus tard*/
	void markLockedCells();

	/** @brief End the elliptic morph algorithm by writing the mesh result in a .vtk or .mli2 file
	 */
	void finalize();

 private:

	/* Locked nodes */
	int m_locked;
	/* Nodes on the surface of the mesh */
	std::vector<TCellID> m_surfNodes;
	std::vector<TCellID> m_locked_faces;
	std::vector<Node> m_morphed;
	std::vector<Node> m_lockedIn;
	std::vector<Node> m_lockedOut;

	std::string m_inputName;
	std::string m_outputName;

	std::vector<std::string> m_to_morph;
	std::vector<std::string> m_ext_lock;
	std::vector<std::string> m_int_lock;

	std::vector<std::vector<double>> m_ellipses;


	/* The mesh to deform */
	Mesh* m_mesh;
};
/*----------------------------------------------------------------------------*/
}// namespace morphmesh
/*----------------------------------------------------------------------------*/
}// namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_ELLIPTICMORPH_H
