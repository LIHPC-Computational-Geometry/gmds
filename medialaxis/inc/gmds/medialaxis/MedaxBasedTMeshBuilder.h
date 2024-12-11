#ifndef GMDS_MEDIALAXIS_MEDAXBASEDTMESHBUILDER_H
#define GMDS_MEDIALAXIS_MEDAXBASEDTMESHBUILDER_H
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
class LIB_GMDS_MEDIALAXIS_API MedaxBasedTMeshBuilder
{

 public:
	/*-------------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param
	 */
	explicit MedaxBasedTMeshBuilder(Mesh &AMedax, Mesh &AMinDel);

	/*-------------------------------------------------------------------------*/
	/** @brief Default destructor.
         *  @param
	 */
	virtual ~MedaxBasedTMeshBuilder()=default;

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
	/** \brief Set the sections IDs
	      *  @param
	 */
	void setSectionID();

	/*-------------------------------------------------------------------------*/
	/** \brief Compute the type of the medial sections
	      *  @param
	 */
	void computeSectionType();

	/*-------------------------------------------------------------------------*/
	/** \brief Returns the set of nodes forming the input section.
	      *  @param AID a section ID
	 */
	std::vector<Node> sectionNodes(int AID);

	/*-------------------------------------------------------------------------*/
	/** \brief Refine the topo rep by adding artificially singular points on the medial axis. YET TO FINALIZE
	      *  @param
	 */
	void refineByAddingSingularNodes();

	/*-------------------------------------------------------------------------*/
	/** \brief Build the nodes of the topological representation (IP, EP, singularities)
	      *  @param
	 */
	void buildTopoRepNodes();

	/*-------------------------------------------------------------------------*/
	/** \brief Build the edges of the topological representation (sections of the medial axis)
	      *  @param
	 */
	void buildTopoRepEdges();

	/*-------------------------------------------------------------------------*/
	/** \brief Write the topological representation
	      *  @param
	 */
	void writeTopoRep(std::basic_string<char> AFileName);

	/*-------------------------------------------------------------------------*/
	/** \brief Return 1 if the wings of the section wrap around the point, 0 if not
	      *  @param AN, AE a node and an edge of the topological representation
	 */
	 int wings(Node &AN, Edge &AE);

	 /*-------------------------------------------------------------------------*/
	 /** \brief Return 1 if the node is the basis of the edge, -1 if it is the end, 0 if non of them
	      *  @param AN, AE a node and an edge of the topological representation
	  */
	 int orientation(Node &AN, Edge &AE);

	 /*-------------------------------------------------------------------------*/
	 /** \brief Set the connectivity of the topological representation
	      *  @param
	  */
	 void setTopoRepConnectivity();

	  /*-------------------------------------------------------------------------*/
	  /** \brief  Return the given node's neighbouring node of the topological representation with respect to the given edge
	      *  @param AN, AE a node and an edge which contains the node
	   */
	  Node getNextSingularNode(Node &AN, Edge &AE);

	  /*-------------------------------------------------------------------------*/
	  /** \brief  Return the given edge's neighbouring edge of the topological representation with respect to the given node
	      *  @param AE, AN an edge and a node contained by th edge
	   */
	  Edge getNextMedialSection(Edge &AE, Node &AN);

	  /*-------------------------------------------------------------------------*/
	  /** \brief  Wander around the medial axis topological representation
	      *  @param
	   */
	  void browseTopoRep();

	  /*-------------------------------------------------------------------------*/
	  /** \brief  Takes as input two nodes of the minimal triangulation, and returns the nodes forming the shortest path to go from the first to the second, 
	   * following the boundary or the constraints.
	      *  @param
	   */
	  std::vector<Node> shortestPathAlongBoundaryOrConstraints(Node &AN1, Node &AN2);

	  /*-------------------------------------------------------------------------*/
	  /** \brief  Mark with 1 edges of the topological represetnation that generate triangles.
	      *  @param
	   */
	  void markEdgesGeneratingTriangles();

	  /*-------------------------------------------------------------------------*/
	  /** \brief  Build the nodes of the T-mesh from the nodes of the minimal triangulation.
	      *  @param
	   */
	  void buildTMeshNodesFromMinDelNodes();

	  /*-------------------------------------------------------------------------*/
	  /** \brief  Build the nodes of the medial axis based block decomposition
	      *  @param
	   */
	  void buildBlockDecompMedialAndBoundaryNodes();

	  /*-------------------------------------------------------------------------*/
	  /** \brief  Write the block decomposition
	      *  @param
	   */
	  void writeBlockDecomp(std::basic_string<char> AFileName);

	  /*-------------------------------------------------------------------------*/
	  /** \brief  
	      *  @param
	   */
	  void buildSection2MedialAndBoundaryNodesAdjacency();

	  /*-------------------------------------------------------------------------*/
	  /** \brief  Build middle nodes 
	      *  @param
	   */
	  void buildMiddleNodes();

	  /*-------------------------------------------------------------------------*/
	  /** \brief  Build blocks
	      *  @param
	   */
	  void buildBlocks();

	  /*-------------------------------------------------------------------------*/
	  /** \brief  Transforms the degenerate quads of the T-mesh into triangles.
	      *  @param
	   */
	  void transformDegenerateQuadsIntoTriangles();

	  /*-------------------------------------------------------------------------*/
	  /** \brief  Transforming triangles of the T-mesh into non-degenerate quads.
	      *  @param
	   */
	  void transformTrianglesIntoQuads();

	  /*-------------------------------------------------------------------------*/
	  /** \brief Sort adjacent edges with respect to their angle with the x-axis.
	      *  @param AN a node
	   */
	  std::vector<Edge> sortedAdjacentEdges(Node &AN);

	  /*-------------------------------------------------------------------------*/
	  /** \brief Take an edge and a node belonging to the edge. Return the edges closest to
	   * the given edge in the list of sorted edges of n.
	      *  @param AE, AN an edge and a node
	   */
	  std::vector<Edge> neighbouringEdges(Edge &AE, Node &AN);

	  /*-------------------------------------------------------------------------*/
	  /** \brief Take an section and a node belonging to the edge. Return the section closest to
	   * the given edge in the list of ordered sections of n.
	      *  @param ASection, AN a section and a node
	   */
	  std::vector<Edge> orderedNeigbourSections(Edge &ASection, Node &AN);

	  /*-------------------------------------------------------------------------*/
	  /** \brief  Set block decomposition connectivity
	      *  @param
	   */
	  void setBlockDecompConnectivity();

	  /*-------------------------------------------------------------------------*/
	  /** \brief  Mark with 1 edges separating different blocks
	      *  @param
	   */
	  void markBlocksSeparatingEdges();

	  /*-------------------------------------------------------------------------*/
	  /** \brief  Add big T-junctions (to have nice pictures)
	      *  @param
	   */
	  void addBigTJunctions();

	  /*-------------------------------------------------------------------------*/
	  /** \brief  Modify faces touching a constraint to ensure good connectivity  TO MODIFY
	      *  @param
	   */
	  void ensureConnectivityThroughConstraints();

	  /*----------------------------------------------------------------------------*/
	  /** @brief Find the connected components of the boundary and attach an ID to each of them.
	 	  *  @param
	  */
	  void setBoundaryConnectedComponents();

	  /*----------------------------------------------------------------------------*/
	  /** @brief Build the final T-mesh.
	 	  *  @param
	  */
	  void buildFinalTMesh();

	  /*----------------------------------------------------------------------------*/
	  /** @brief Write the final T-mesh.
	 	  *  @param
	  */
	  void writeFinalTMesh(std::basic_string<char> AFileName);

	  /*-------------------------------------------------------------------------*/
	  /** \brief  Set block connectivity of the T-mesh.
	      *  @param
	   */
	  void setFinalTMeshConnectivity();

	  /*-------------------------------------------------------------------------*/
	  /** \brief  Set block connectivity of the T-mesh.
	      *  @param
	   */
	  Mesh getFinalTMesh();

	  /*-------------------------------------------------------------------------*/
	  /** \brief  Mark the internal constraints on the final T-mesh.
	      *  @param
	   */
	  void markInternalConstraintsOnFinalTMesh();

 private:

	// A mesh representation of a medial axis, as built by the class MedialAxisBuilder
	Mesh* m_medax;

	// The minimal Delaunay triangulation from which the medial axis was created
	Mesh* m_min_delaunay;

	// Topological representation of the medial axis (used to build the quad)
	Mesh* m_topological_representation;
	// Nb of sections of the medial axis
	int m_nb_medial_sections;
	// medial axis based T-mesh
	Mesh* m_t_mesh;
	// Final T-mesh
	Mesh* m_final_t_mesh;

	// Connected components of the T-mesh
	std::vector<std::vector<Edge>> m_boundary_connected_components;
};
/*----------------------------------------------------------------------------*/
}  // end namespace medialaxis
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif  // GMDS_MEDIALAXIS_MEDAXBASEDTMESHBUILDER_H