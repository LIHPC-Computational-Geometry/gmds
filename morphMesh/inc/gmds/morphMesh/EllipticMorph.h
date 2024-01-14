
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
	EllipticMorph(Mesh* AMesh);

	/*------------------------------------------------------------------------*/
	/** @brief Destructor.
	 */
	~EllipticMorph();

	/*------------------------------------------------------------------------*/
	/** @brief Constructor.
	 */
	void execute(const std::vector<std::vector<double>>& AEllipses);

 private:

	/*Surement a retirer plus tard*/
	void markLockedCells();


 private:

	/* Locked nodes */
	int m_locked;
	/* Nodes on the surface of the mesh */
	std::vector<TCellID> m_surfNodes;
	std::vector<TCellID> m_locked_faces;


	/* The mesh to deform */
	Mesh* m_mesh;
};
/*----------------------------------------------------------------------------*/
}// namespace morphmesh
/*----------------------------------------------------------------------------*/
}// namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_ELLIPTICMORPH_H
