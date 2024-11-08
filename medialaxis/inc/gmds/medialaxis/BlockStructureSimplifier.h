#ifndef GMDS_BLOCKSTRUCTURESIMPLIFIER_H
#define GMDS_BLOCKSTRUCTURESIMPLIFIER_H
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_MEDIALAXIS_export.h"
#include "gmds/medialaxis/MedialAxis2D.h"
#include "gmds/medialaxis/NonConformalHalfEdge.h"
#include "gmds/medialaxis/MedialAxisMath.h"
#include "gmds/medialaxis/QuantizationGraph.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKWriter.h"
#include "gmds/ig/MeshDoctor.h"
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/** \class  dummy
 *  \brief  dummy class.
 */
class LIB_GMDS_MEDIALAXIS_API BlockStructureSimplifier
{
 private:
	// Non conformal quad block decomposition
	Mesh* m_mesh;
	// Corresponding non-conformal half edges
	std::vector<NonConformalHalfEdge> m_half_edges;
	// Quantization graph
	QuantizationGraph* m_quantization_graph;
	// Half edges length
	std::vector<int> m_half_edges_lengths;
	// Mesh after removing half edges of length zero
	Mesh* m_simplified_mesh;

 public:

	/*-------------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param AMesh
	 */
	explicit BlockStructureSimplifier(Mesh &AMesh);

	/*-------------------------------------------------------------------------*/
	/** @brief Default destructor.
         *  @param
	 */
	virtual ~BlockStructureSimplifier()=default;

	/*-------------------------------------------------------------------------*/
	/** @brief Getters.
         *  @param
	 */
	std::vector<NonConformalHalfEdge> halfEdges(){return m_half_edges;}
	std::vector<int> halfEdgesLengths(){return m_half_edges_lengths;}

	/*-------------------------------------------------------------------------*/
	/** @brief Build the half edges.
         *  @param
	 */
	void buildHalfEdges();

	/*-------------------------------------------------------------------------*/
	/** @brief Build the quantization graph nodes.
         *  @param
	 */
	void buildQuantizationGraphNodes();

	/*-------------------------------------------------------------------------*/
	/** @brief Build the connected component of the given half edge.
         *  @param AHalfEdgeID an half edge ID
	 */
	void buildConnectedComponent(int AHalfEdgeID);

	/*-------------------------------------------------------------------------*/
	/** @brief Return the next of the next half edge of the given half edge.
         *  @param AHalfEdgeID an half edge ID
	 */
	int oppositeInQuad(int AHalfEdgeID);

	/*-------------------------------------------------------------------------*/
	/** @brief Build the connected component graph.
         *  @param
	 */
	void buildQuantizationGraph();

	/*-------------------------------------------------------------------------*/
	/** @brief Get the quantization graph.
         *  @param
	 */
	QuantizationGraph* getQuantizationGraph();

	/*-------------------------------------------------------------------------*/
	/** @brief Set half edges length.
         *  @param
	 */
	void setHalfEdgesLength();

	/*-------------------------------------------------------------------------*/
	/** @brief Return the set of sides of the input geometry, each side being a set of half-edge Ids.
         *  @param
	 */
	std::vector<std::vector<TCellID>> sides();

	/*-------------------------------------------------------------------------*/
	/** @brief Mark the T-junctions.
         *  @param
	 */
	void markTJunctions();

	/*-------------------------------------------------------------------------*/
	/** @brief Return true if the given node is on the boundary.
         *  @param
	 */
	bool isOnBoundary(Node AN);

	/*-------------------------------------------------------------------------*/
	/** @brief Return the set of half-edges that must be non zero, ie that are adjacent to a singularity.
         *  @param
	 */
	std::vector<TCellID> nonZeroHalfEdges();

	/*-------------------------------------------------------------------------*/
	/** @brief Return the set of groups of nodes, each group being a set of nodes connected by half-edges of length 0.
         *  @param
	 */
	std::vector<std::vector<Node>> groupsOfConfundedNodes();

	/*-------------------------------------------------------------------------*/
	/** @brief Return true if the given node is a corner of the boundary.
         *  @param
	 */
	bool isABoundaryCorner(Node AN);

	/*-------------------------------------------------------------------------*/
	/** @brief Compute the disappearance gain for each node.
         *  @param 
	 */
	void computeChosingGains();

	/*-------------------------------------------------------------------------*/
	/** @brief Return a representative for the given group of confunded nodes.
         *  @param AV
	 */
	Node representative(std::vector<Node> AV);

	/*-------------------------------------------------------------------------*/
	/** @brief Build the simplified mesh.
         *  @param 
	 */
	void buildSimplifiedMesh();

	/*-------------------------------------------------------------------------*/
	/** @brief Set the connectivity of the simplified mesh.
         *  @param 
	 */
	void setSimplifiedMeshConnectivity();

	/*-------------------------------------------------------------------------*/
	/** @brief Get the simplified mesh.
         *  @param 
	 */
	Mesh getSimplifiedMesh();

	/*-------------------------------------------------------------------------*/
	/** @brief Write the simplified mesh.
         *  @param 
	 */
	void writeSimplifiedMesh(std::basic_string<char> AFileName);

    /*-------------------------------------------------------------------------*/
	/** @brief Execute.
         *  @param 
	 */
	void execute();
};
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_BLOCKSTRUCTURESIMPLIFIER_H