/*----------------------------------------------------------------------------*/
#ifndef GMDS_MORPHMESH_H
#define GMDS_MORPHMESH_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_MORPHMESH_export.h"
#include "gmds/ig/Mesh.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace morphmesh {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_MORPHMESH_API MorphMesh
{

 public:
	/*------------------------------------------------------------------------*/
	/** @brief Constructor.
	 */
	MorphMesh(Mesh* AMesh, const std::vector<math::Point> &APoints, double ARadius);

	/*------------------------------------------------------------------------*/
	/** @brief Destructor.
	 */
	~MorphMesh();

	/*------------------------------------------------------------------------*/
	/** @brief Constructor.
	 */
	 void execute();

  private:

	 /*Surement a retirer plus tard*/
	 void markLockedCells();

	 math::Point computeLocalOrigin(const math::Point& AP);

	 double findHomothetyRatio(const math::Point& AP, const Node& ANode);

  private:

	 /* Locked nodes */
	 int m_locked;
	 /* List of "homothetic points" that will define the mesh deformation */
	 std::vector<math::Point> m_targets;
	 /* The radius of an area that will be modified by a homothetic point*/
	 double m_radius;
	 /* Nodes on the surface of the mesh */
	 std::vector<TCellID> m_surfNodes;

	 /* The mesh to deform */
	 Mesh* m_mesh;
};
/*----------------------------------------------------------------------------*/
}// namespace morphmesh
/*----------------------------------------------------------------------------*/
}// namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_MORPHMESH_H
