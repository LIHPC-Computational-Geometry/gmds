//
// Created by chenyt on 22/04/24.
//

#ifndef GMDS_MEDIALAXIS2DBUILDER_H
#define GMDS_MEDIALAXIS2DBUILDER_H
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_MEDIALAXIS_export.h"
#include "gmds/medialaxis/MedialAxis2D.h"
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace medialaxis{
/*----------------------------------------------------------------------------*/
/** \class  dummy
 *  \brief  dummy class.
 */
class LIB_GMDS_MEDIALAXIS_API MedialAxis2DBuilder{

 public:
	/*------------------------------------------------------------------------*/
	/** \brief friend class to access to protected/private data
	 */
	friend class Mesh;
	friend class MedialAxis2D;

	/*-------------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS
	} STATUS;


	/*-------------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param
	 */

	explicit MedialAxis2DBuilder(Mesh &AMesh);

 public:
	// Primal minimal Delaunay mesh
	Mesh* m_mesh;
	// Dual Voronoï medial axis
	MedialAxis2D *m_voronoi_medax;
	// Smoothed medial axis built from the Voronoï one
	MedialAxis2D * m_smoothed_medax;
	// Graph of the boundary connected components
	Mesh* m_boundary_connected_components_graph;
	// Set of couples of boundary points connecting boundary connected components
	std::vector<std::vector<math::Point>> m_connected_boundary_points;


	/*-------------------------------------------------------------------------*/
	/** @brief Default destructor.
         *  @param
	 */
	virtual ~MedialAxis2DBuilder();

	/*-------------------------------------------------------------------------*/
	/** \brief Get the medial axis
	 */
	MedialAxis2D * getMedialObject();

	/*-------------------------------------------------------------------------*/
	/** \brief Get the smoothed medial axis
	 */
	MedialAxis2D * getSmoothedMedialObject();

	/*-------------------------------------------------------------------------*/
	/** \brief Get the mesh
	 */
	Mesh* getMesh();

	/*-------------------------------------------------------------------------*/
	/** \brief Classify the edges
	 */
	void classEdges();

	/*-------------------------------------------------------------------------*/
	/** \brief Classify the faces
	 */
	void classFaces();

	/*-------------------------------------------------------------------------*/
	/** \brief Mark with 1 boundary edges belonging to a corner
	 */
	void markCorners();

	/*-------------------------------------------------------------------------*/
	/** \brief Attach an ID to each side of the geometry
	 */
	void setSideId();

	/*-------------------------------------------------------------------------*/
	/** \brief Mark with 1 Delaunay edges which correspond to internal constraints
	 */
	void markIntConstraints();

	/*-------------------------------------------------------------------------*/
	/** \brief Compute the inner product between medial radii and the corresponding boundary edges, which should be 0
	 * 	\param ABoundaryCurvatureTol a maximum distance away from 1 of the norm of the cosine of the boundary curvature angle
	 */
	void setMedialRadiiOrthogonalityDefault(const double &ABoundaryCurvatureTol);

	/*-------------------------------------------------------------------------*/
	/** \brief Attach to each medial point its corresponding touching points
	 */
	void setTouchingPoints();

	/*-------------------------------------------------------------------------*/
	/** \brief Set the medial radius on medial points
	 */
	void setMedialRadius();

	/*-------------------------------------------------------------------------*/
	/** \brief Mark the nodes of the input geometry considered as details
	 */
	void markGeometryDetailsOnPoints();

	/*-------------------------------------------------------------------------*/
	/** \brief Mark the edges of the input geometry considered as details
	 */
	void markGeometryDetailsOnEdges();

	/*-------------------------------------------------------------------------*/
	/** \brief Returns the maximal length of the Delaunay edges (gives a typical size for the domain)
	 */
	double maxDelEdgeLength();

	/*-------------------------------------------------------------------------*/
	/** \brief Build the Voronoï medial points and edges
	 */
	void buildVoronoiMedialPointsAndEdges();

	/*-------------------------------------------------------------------------*/
	/** \brief Build the Voronoï medial axis and set all its attributes
	 */
	void buildVoronoiMedialAxis();

	/*-------------------------------------------------------------------------*/
	/** \brief Build the smoothed medial axis from the groups of nodes of the Voronoï medial axis
	 * 	\param AGroups the set of medial points groups
	 */
	void buildSmoothedMedaxFromVoronoi(const std::vector<std::vector<Node>> AGroups);

	/*-------------------------------------------------------------------------*/
	/** \brief Identify and mark small details on the geometry using the medial axis
	 */
	void identifyDetails();

	/*-------------------------------------------------------------------------*/
	/** \brief Place singularities using the medial axis
	 */
	void placeSingularities();

	/*-------------------------------------------------------------------------*/
	/** \brief Clean the Delaunay mesh by deleting the isolated points
	 */
	void deleteIsolatedPoints();

	/*----------------------------------------------------------------------------*/
	/** @brief Find the connected components of the boundary and attach an ID to each of them.
         *  @param AN a reference, updated to the number of connected components
	 */
	void setBoundaryConnectedComponents(int &AN);

	/*----------------------------------------------------------------------------*/
	/** @brief Attach to each triangle its corresponding boundary connected components.
         *  @param
	 */
	void setBoundaryConnectedComponentsOnFaces();

	/*----------------------------------------------------------------------------*/
	/** @brief Attach to each medial point its corresponding boundary connected components.
         *  @param
	 */
	void setBoundaryConnectedComponentsOnMedPoints();

	/*----------------------------------------------------------------------------*/
	/** @brief Build the graph of the connected components from a components adjacency matrix.
         *  @param AM a components adjacency matrix
	 */
	void buildBoundaryComponentsGraph(Eigen::MatrixXi &AM);

	/*----------------------------------------------------------------------------*/
	/** @brief Return a list of edges forming a spanning tree of the connected components graph.
         *  @param
	 */
	std::vector<Edge> connectedComponentsSpanningTree();

	/*----------------------------------------------------------------------------*/
	/** @brief Return a list of couple of points, each couple being the optimal connexion between two boundary connected components.
         *  @param AST a spanning tree of the connected components graph
	 */
	std::vector<std::vector<math::Point>> connectedBoundaryPoints(std::vector<Edge> &AST);

	/*----------------------------------------------------------------------------*/
	/** @brief Connect the boundary connected components using the medial axis.
         *  @param
	 */
	void connectBoundaryConnectedComponents();

	/*-------------------------------------------------------------------------*/
	/** \brief
	 */
	MedialAxis2DBuilder::STATUS execute();

};
/*----------------------------------------------------------------------------*/
}  // end namespace medialaxis
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_MEDIALAXIS2DBUILDER_H
