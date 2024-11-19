//
// Created by chenyt on 26/09/24.
//

#ifndef GMDS_QUANTIZATIONSOLVER_H
#define GMDS_QUANTIZATIONSOLVER_H
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
class LIB_GMDS_MEDIALAXIS_API QuantizationSolver
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
	// Mesh after adding a singularity dipole
	Mesh* m_fixed_mesh;

 public:

	/*-------------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param AMesh
	 */
	explicit QuantizationSolver(Mesh &AMesh);

	/*-------------------------------------------------------------------------*/
	/** @brief Default destructor.
         *  @param
	 */
	virtual ~QuantizationSolver()=default;

	/*-------------------------------------------------------------------------*/
	/** @brief Getters.
         *  @param
	 */
	std::vector<NonConformalHalfEdge> halfEdges(){return m_half_edges;}
	std::vector<int> halfEdgesLengths(){return m_half_edges_lengths;}

	/*-------------------------------------------------------------------------*/
	/** @brief Find and mark the T-junctions of the input mesh. We associate to each t-junction the id of the face for which it is a non_conformal node. WARNING : GEOMETRY-DEPENDANT
         *  @param
	 */
	void findTJunctions();

	/*-------------------------------------------------------------------------*/
	/** @brief Regroups edges of the given face in groups of aligned edges.
         *  @param AF a face
	 */
	std::vector<std::vector<Edge>> alignedEdgesGroups(Face &AF);

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
	/** @brief Add on the fixed mesh the already existing nodes.
         *  @param 
	 */
	void addOldNodesOnFixedMesh();

	/*-------------------------------------------------------------------------*/
	/** @brief Write the fixed mesh.
         *  @param 
	 */
	void writeFixedMesh(std::basic_string<char> AFileName);

	/*-------------------------------------------------------------------------*/
	/** @brief Mark with 1 the faces affected when adding a singularity dipole in the face deliminted by the two given half edges, 
	 * ie the face itself and its neighbours not adjacent to AId1 or AId2.
         *  @param AId1, AId2 two half edges Ids delimiting a quad.
	 */
	void markAffectedFaces(int AId1, int AId2);

	/*-------------------------------------------------------------------------*/
	/** @brief Add on the fixed mesh the already existing faces not affected by the dipole.
         *  @param 
	 */
	void addOldFacesOnFixedMesh();

	/*-------------------------------------------------------------------------*/
	/** @brief Build faces of the fixed mesh oppposite to the given half-edge, taking into account eventual new nodes.
         *  @param AHalfEdgeId, AN1, AN2 a half-edge Id and two nodes of the fixed-mesh belonging to the half-edge.
	 */
	void buildOppFacesInFixedMesh(int AHalfEdgeId, Node AN1, Node AN2);

	/*-------------------------------------------------------------------------*/
	/** @brief Build the fixed mesh, where a singularity dipole has been added and oriented in the face delimited by two half edges whose ids are given.
         *  @param AId1, AId2 two half-edge ids. The dipole 3-5 is placed in the direction AId1->AId2.
	 */
	void buildFixedMesh(int AId1, int AId2);

	/*-------------------------------------------------------------------------*/
	/** @brief Get the fixed mesh.
         *  @param 
	 */
	Mesh getFixedMesh();

	/*-------------------------------------------------------------------------*/
	/** @brief Set the connectivity of the fixed mesh.
         *  @param 
	 */
	void setFixedMeshConnectivity();

	/*-------------------------------------------------------------------------*/
	/** @brief Build the complete quantization solution.
         *  @param 
	 */
	std::vector<std::vector<Node>> buildCompleteSolution();
};
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_QUANTIZATIONSOLVER_H
