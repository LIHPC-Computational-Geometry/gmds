//
// Created by chenyt on 10/07/24.
//

#ifndef GMDS_CROSSFIELD_H
#define GMDS_CROSSFIELD_H

#include "LIB_GMDS_MEDIALAXIS_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/math/Matrix.h>
#include <gmds/math/Vector.h>
#include <gmds/math/Point.h>
#include <gmds/math/Cross2D.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace medialaxis{
/*----------------------------------------------------------------------------*/
/** \class  dummy
 *  \brief  dummy class.
 */
class CrossField
{
 public:
	/*-------------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param
	 */
	explicit CrossField(Mesh &AMesh);

	/*-------------------------------------------------------------------------*/
	/** @brief Default destructor.
         *  @param
	 */
	virtual ~CrossField();

	/*----------------------------------------------------------------------------*/
	/** @brief Find the connected components of the boundary and attach an ID to each of them.
         *  @param AN a reference, updated to the number of connected components
	 */
	void setBoundaryConnectedComponents(int &AN);

	/*----------------------------------------------------------------------------*/
	/** @brief Set the boundaries connexions.
         *  @param AV a vector of couples of connecting nodes
	 */
	void setBoundariesConnexions(std::vector<std::vector<Node>> AV);

	/*----------------------------------------------------------------------------*/
	/** @brief Return a set of couples of nodes to link to connect boundaries.
         *  @param AV a vector of couples of connected boundary points given by the medial axis
	 */
	std::vector<std::vector<Node>> connectedBoundaryNodes(std::vector<std::vector<math::Point>> &AV);

	/*----------------------------------------------------------------------------*/
	/** @brief Find an edges path connecting two nodes in the mesh, using Dijkstra algorithm.
         *  @param AN1, AN2 two nodes
	 */
	std::vector<Edge> findPath(Node &AN1, Node &AN2);

	/*-------------------------------------------------------------------------*/
	/** @brief Signed area of a triangle.
         *  @param AP1, AP2, AP3 three points forming a triangle
	 */
	double signedArea(const math::Point &AP1, const math::Point &AP2, const math::Point &AP3);

	/*-------------------------------------------------------------------------*/
	/** @brief Positively orientate the two nodes on an edge with respect to the orientation given by a triangle.
         *  @param AE, AF an edge and a triangle face
	 */
	std::vector<Node> orientateEdge(const Edge AE, const Face AF);

	/*-------------------------------------------------------------------------*/
	/** @brief Return the sign of the orientation of the the given edge in the given triangle.
         *  @param AE, AF an edge and a triangle face
	 */
	int orientation(const Edge AE, const Face AF);

	/*-------------------------------------------------------------------------*/
	/** @brief Compute the three barycentric coordinates of a point in a triangle.
         *  @param AP, AF a point and a triangle face
	 */
	std::vector<double> barycentricCoordinates(const math::Point AP, const Face AF);

	/*----------------------------------------------------------------------------*/
	/** @brief Return the Id of the face of the mesh containing the given point.
         *  @param AP a point
	 */
	TCellID locateInMesh(const math::Point AP);

	/*----------------------------------------------------------------------------*/
	/** @brief Set the singularity indexes.
         *  @param ANodes, AIndexes a vector of singular nodes and the vector of the corresponding indexes
	 */
	void setSingularitiesIndexes(const std::vector<Node> ANodes, const std::vector<double> AIndexes);

	/*----------------------------------------------------------------------------*/
	/** @brief Check if an edge is an interior edge.
         *  @param AE an edge
	 */
	bool isInterior(const Edge &AE);

	/*----------------------------------------------------------------------------*/
	/** @brief Check if a node is an interior node.
         *  @param AN a node
	 */
	bool isInterior(const Node &AN);

	/*----------------------------------------------------------------------------*/
	/** @brief Number of interior nodes.
         *  @param
	 */
	int NbInteriorNodes();

	/*----------------------------------------------------------------------------*/
	/** @brief Number of interior edges.
         *  @param
	 */
	int NbInteriorEdges();

	/*----------------------------------------------------------------------------*/
	/** @brief Return the complex number associated to a point.
         *  @param AP a point
	 */
	std::complex<double> point2complex(math::Point &AP);

	/*----------------------------------------------------------------------------*/
	/** @brief Initialize crosses reference angle.
         *  @param
	 */
	void initializeCrossReferenceAngle();

	/*----------------------------------------------------------------------------*/
	/** @brief Initialize adjustment angles on boundary edges.
         *  @param
	 */
	void initializeAdjustmentAngles();

	/*----------------------------------------------------------------------------*/
	/** @brief Compute the adjustment angles for each interior edge.
         *  @param
	 */
	void computeAdjustmentAngles();

	/*----------------------------------------------------------------------------*/
	/** @brief Propagate crosses from boundary.
         *  @param
	 */
	void propagateRefAngleFromBoundary();

	/*----------------------------------------------------------------------------*/
	/** @brief Build crosses from reference angles.
         *  @param
	 */
	void buildCrossesFromRefAngles();

	/*----------------------------------------------------------------------------*/
	/** @brief Build the cross field.
         *  @param
	 */
	void buildCrossField();

 private:
	Mesh* m_mesh;
	gmds::Variable<math::Cross2D>* m_cross_field_2D;
	std::vector<std::vector<Node>> m_boundaries_connexions;
};
/*----------------------------------------------------------------------------*/
}  // end namespace medialaxis
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/

#endif     // GMDS_CROSSFIELD_H
