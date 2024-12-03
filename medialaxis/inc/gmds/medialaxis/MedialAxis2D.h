#ifndef GMDS_MEDIALAXIS_MEDIALAXIS_H
#define GMDS_MEDIALAXIS_MEDIALAXIS_H
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_MEDIALAXIS_export.h"
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace medialaxis{
/*----------------------------------------------------------------------------*/
/** \class  dummy
 *  \brief  dummy class.
 */
class LIB_GMDS_MEDIALAXIS_API MedialAxis2D
{

 public:
	/*-------------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param
	 */
	explicit MedialAxis2D();

	/*-------------------------------------------------------------------------*/
	/** @brief Default destructor.
         *  @param
	 */
	virtual ~MedialAxis2D();

	/*-------------------------------------------------------------------------*/
	/** \brief Add a medial point.
	      *  @param a point.
	 */
	Mesh getMeshRepresentation();

	/*-------------------------------------------------------------------------*/
	/** \brief Add a medial point.
	      *  @param a point.
	 */
	Node newMedPoint(const math::Point &APnt);

	/*-------------------------------------------------------------------------*/
	/** \brief Get a medial point.
	      *  @param a medial point Id.
	 */
	Node getMedPoint(const TCellID APointID);

	/*-------------------------------------------------------------------------*/
	/** \brief Add a medial edge.
	      *  @param two medial points id.
	 */
	Edge newMedEdge(const TCellID &AN1, const TCellID &AN2);

	/*-------------------------------------------------------------------------*/
	/** \brief Update the singularities attributes
	      *  @param
	 */
	void setSingularities();

	/*-------------------------------------------------------------------------*/
	/** \brief Get m_singular_nodes
	      *  @param
	 */
	std::vector<Node> getSingularNodes();

	/*-------------------------------------------------------------------------*/
	/** \brief Get m_singularity_indexes
	      *  @param
	 */
	std::vector<double> getSingularityIndexes();

	/*-------------------------------------------------------------------------*/
	/** \brief Return the minimum medial edge length
	      *  @param
	 */
	double minMedEdgeLength();

	/*-------------------------------------------------------------------------*/
	/** \brief Return the maximal medial edge length
	      *  @param
	 */
	double maxMedEdgeLength();

	/*-------------------------------------------------------------------------*/
	/** \brief Return the mean medial edge length
	      *  @param
	 */
	double meanMedEdgeLength();

	/*-------------------------------------------------------------------------*/
	/** \brief Return the number of medial points
	      *  @param
	 */
	int getNbMedPoints();

	/*-------------------------------------------------------------------------*/
	/** \brief Update the medial axis connectivity
	      *  @param A medial point id
	 */
	void updateConnectivity();

	/*-------------------------------------------------------------------------*/
	/** \brief Return groups of medial points satisfying the following criterion: two medial points connected by an edge of length
	 * inferior to the given tolerance are in the same group
	      *  @param ATol a tolerance
	 */
	std::vector<std::vector<Node>> medialPointsGroups(const double &ATol);

	/*-------------------------------------------------------------------------*/
	/** \brief Return the group ID of the given medial point
	      *  @param ANode a node
	 */
	int medialPoint2Group(const Node &ANode);

	/*-------------------------------------------------------------------------*/
	/** \brief Refine the medial axis by adding medial points on too long edges (edges of length > tolerance)
	      *  @param ATol a tolerance
	 */
	void refine(const double& ATol);

	/*-------------------------------------------------------------------------*/
	/** \brief Set the ID of the primal triangle
	      *  @param AMedPointID, ATriID, a medial point ID and its primal triangle ID
	 */
	void setPrimalTriangleID(TCellID AMedPointID, TCellID ATriID);

	/*-------------------------------------------------------------------------*/
	/** \brief Returns the ID of the primal triangle
	      *  @param A medial point id
	 */
	TCellID primalTriangleID(TCellID AMedPointID);

	/*-------------------------------------------------------------------------*/
	/** \brief Set the medial point type
	      *  @param A medial point id
	 */
	void setMedialPointType();

	/*-------------------------------------------------------------------------*/
	/** \brief Get the medial point type
	      *  @param AId medial point id
	 */
	int getMedialPointType(TCellID AId);

	/*-------------------------------------------------------------------------*/
	/** \brief Set the medial edge type
	      *  @param
	 */
	void setMedialEdgeType();

	/*-------------------------------------------------------------------------*/
	/** \brief Attach to the medial point its corresponding boundary connected components
	      *  @param AId, AV a node ID and a vector of connected components IDs
	 */
	void setBoundaryConnectedComponentsIDs(TCellID AId, std::vector<int> &AV);

	/*-------------------------------------------------------------------------*/
	/** \brief Return the adjacency matrix M of the boundary connected components. If components
	 * I and J are connected by the medial axis, M(I,J) = ID of the medial point that connects them best
	      *  @param AN the number of boundary connected components
	 */
	Eigen::MatrixXi findOptimalBoundaryConnectedComponentsConnexion(int &AN);

	/*-------------------------------------------------------------------------*/
	/** \brief Set the cosine of the medial angle at medial points
	      *  @param
	 */
	void setCosMedialAngle();

	/*-------------------------------------------------------------------------*/
	/** \brief Smooth the values by placing each of them at the barycenter of its neighbours
	      *  @param
	 */
	void smoothCosMedialAngle();

	/*-------------------------------------------------------------------------*/
	/** \brief Attach to a medial point its corresponding touching points
	      *  @param A medial point ID, a points vector
	 */
	void setTouchingPoints(const TCellID &APointID, std::vector<math::Point> APoints);

	/*-------------------------------------------------------------------------*/
	/** \brief Attach to a medial point its corresponding tengancy nodes
	      *  @param A medial point ID, a vector of nodes
	 */
	void setTangentNodes(const TCellID &APointID, std::vector<Node> ANodes);

	/*-------------------------------------------------------------------------*/
	/** \brief Attach to a medial point its corresponding dual triangles
	      *  @param A medial point ID, a vector of faces
	 */
	void setDualTriangles(const TCellID &APointID, std::vector<Face> ATriangles);

	/*-------------------------------------------------------------------------*/
	/** \brief Get the touching points of a given medial point
	      *  @param A medial point ID
	 */
	std::vector<math::Point> getTouchingPoints(const TCellID &APointID);

	/*-------------------------------------------------------------------------*/
	/** \brief Get the tengancy nodes of a given medial point
	      *  @param A medial point ID
	 */
	std::vector<Node> getTangentNodes(const TCellID &APointID);

	/*-------------------------------------------------------------------------*/
	/** \brief Set the optimal change in angle of a cross following the curve formed by two adjacent medial radii
	      *  @param
	 */
	 void setFluxThroughMedialRadii();

	 /*-------------------------------------------------------------------------*/
	 /** \brief Return the orientation -1,1 or 0 of a curve formed by two vectors with respect to a third vector
	      *  @param Two vectors AU and VU forming a curve to orientate, a vector AE with respect to which to orientate the curve
	  */
	 double orientation(const math::Vector &AU, const math::Vector &AV, const math::Vector &AE);

	 /*-------------------------------------------------------------------------*/
	 /** \brief Set the flux residual across edges and intersection points
	      *  @param A medial edge ID, a value
	  */
	 void setFluxResidual();

	/*-------------------------------------------------------------------------*/
	/** \brief Set the medial radius at a given medial point
	      *  @param A medial point ID, a value
	 */
	void setMedialRadius(const TCellID &APointID, double AValue);

	/*-------------------------------------------------------------------------*/
	/** \brief Set the medial radius orthogonality default at a given medial point
	      *  @param A medial point ID, a value
	 */
	void setMedialRadiusOrthogonalityDefault(const TCellID &APointID, double AValue);

	/*-------------------------------------------------------------------------*/
	/** \brief Get the medial radius at a given medial point
	      *  @param A medial point ID
	 */
	double getMedialRadius(const TCellID &APointID);

	/*-------------------------------------------------------------------------*/
	/** \brief Get the medial radius orthogonality default at a given medial point
	      *  @param A medial point ID
	 */
	double getMedialRadiusOrthogonalityDefault(const TCellID &APointID);

	/*-------------------------------------------------------------------------*/
	/** \brief Place the singularities, avoiding placing them too close to the boundary, ie where the medial radius is too small,
	 * which is characterized by checking that the medial edge length is smaller than a given fraction of the medial radius
	      *  @param AMedRadiusFraction a medial radius fraction
	 */
	void placeSingularities(const double& AMedRadiusFraction);

	/*-------------------------------------------------------------------------*/
	/** \brief Remove singularity dipoles
	      *  @param
	 */
	void removeSingularityDipoles();

	/*-------------------------------------------------------------------------*/
	/** \brief Check if placed singularities are coherent, ie if the med radius orthogonality default is smaller than a tolerance
	      *  @param AOrthogonalityDefaultTol A tolerance
	 */
	void checkSingularities(const double& AOrthogonalityDefaultTol);

	/*-------------------------------------------------------------------------*/
	/** \brief  Check if singularities are close to an intersection point up to ANbPoints points. If it is the case,
	 * move the singularity to the IP.
	      *  @param ANbPoints a number of points
	 */
	void moveSingularitiesToIPs(const double& ANbPoints);

	/*-------------------------------------------------------------------------*/
	/** \brief Get a point's adjacent point with respect to one of its edges
	      *  @param A medial point id and a medial edge id.
	 */
	Node getNextPoint(const TCellID &APointID, const TCellID &AEdgeID);

	/*-------------------------------------------------------------------------*/
	/** \brief Get an edge's adjacent edge with respect to one of its points
	      *  @param A medial edge id and a medial point id.
	 */
	Edge getNextEdge(const TCellID &AEdgeID, const TCellID &APointID);

	/*-------------------------------------------------------------------------*/
	/** \brief Get the extremal point of a node's branch in the direction given by one of its adjacent edges
	      *  @param A medial point id and a medial edge id.
	 */
	Node getExtremPoint(const TCellID &APointID, const TCellID &AEdgeID);

	/*-------------------------------------------------------------------------*/
	/** \brief Set the type of each node's branch: 0 if the branch is interior, 1 if it is adjacent to the boundary, 2 if it is an intersection point
	      *  @param
	 */
	void setBranchTypeOnPoints();

	/*-------------------------------------------------------------------------*/
	/** \brief Set the type of each edge's branch: 0 if the branch is interior, 1 if it is adjacent to the boundary
	      *  @param
	 */
	void setBranchTypeOnEdges();

	/*-------------------------------------------------------------------------*/
	/** \brief Set the id of each branch, on their edges, and returns the number of branches
	      *  @param
	 */
	int setBranchIdOnEdges();

	/*-------------------------------------------------------------------------*/
	/** \brief Set the id of each branch, on their points
	      *  @param
	 */
	void setBranchIdOnPoints();

	/*-------------------------------------------------------------------------*/
	/** \brief Mark with 1 the edges whose branch corresponds to details, 0 the others
	      *  @param ATol a tolerance, that defines the typical size of a detail
	 */
	void identifyDetailsBranches(const double& ATol);

	/*-------------------------------------------------------------------------*/
	/** \brief Returns the value of corresponds_to_detail of a medial point
	      *  @param APointID a point ID
	 */
	int correspondsToDetail(const TCellID & APointID);

	/*-------------------------------------------------------------------------*/
	/** \brief Write the medial axis
	      *  @param
	 */
	void write(std::basic_string<char> AFileName);

 private:

	// Geometrical representation of the medial axis
	Mesh* m_mesh_representation;

	// Singularities identified on the medial axis
	std::vector<Node> m_singular_nodes;
	std::vector<double> m_singularity_indexes;
};
/*----------------------------------------------------------------------------*/
}  // end namespace medialaxis
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif  // GMDS_MEDIALAXIS_MEDIALAXIS_H