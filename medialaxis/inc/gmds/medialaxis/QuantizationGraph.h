//
// Created by chenyt on 27/09/24.
//

#ifndef GMDS_QUANTIZATIONGRAPH_H
#define GMDS_QUANTIZATIONGRAPH_H
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_MEDIALAXIS_export.h"
#include "gmds/medialaxis/MedialAxisMath.h"
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/** \class  dummy
 *  \brief  dummy class.
 */
class LIB_GMDS_MEDIALAXIS_API QuantizationGraph
{
 private:
	// Mesh representation
	Mesh* m_mesh_representation;

	// Non zero conditions for the quantization solution :

	// Set of verticies with non zero condition
	std::vector<TCellID> m_non_zero_verticies;

	// Set of groups of verticies, with condition at least one non zero in each group
	std::vector<std::vector<TCellID>> m_non_zero_groups;

 public:

	/*-------------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param
	 */
	explicit QuantizationGraph();

	/*-------------------------------------------------------------------------*/
	/** @brief Default destructor.
         *  @param
	 */
	virtual ~QuantizationGraph()=default;

	/*-------------------------------------------------------------------------*/
	/** @brief Accessor to non zero verticies.
         *  @param
	 */
	void nonZeroVerticies(std::vector<TCellID> AV){m_non_zero_verticies = AV;}

	/*-------------------------------------------------------------------------*/
	/** @brief Accessor to non zero groups.
         *  @param
	 */
	void nonZeroGroups(std::vector<std::vector<TCellID>> AV){m_non_zero_groups = AV;}

	/*-------------------------------------------------------------------------*/
	/** @brief New node.
         *  @param
	 */
	Node newNode();

	/*-------------------------------------------------------------------------*/
	/** @brief New edge.
         *  @param AId1, AId2 two nodes Ids
	 */
	Edge newEdge(TCellID AId1, TCellID AId2);

	/*-------------------------------------------------------------------------*/
	/** @brief Display the edges of the graph.
         *  @param
	 */
	void display();

	/*-------------------------------------------------------------------------*/
	/** @brief Return 1 if the node has already been visited during the building of the edges, 0 if not.
         *  @param AID
	 */
	int alreadyVisited(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** @brief  Mark the given node as already visited.
         *  @param AID
	 */
	void markAsVisited(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** @brief Return 1 if the node has already been pushed during the building of the edges, 0 if not.
         *  @param AID
	 */
	int alreadyPushed(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** @brief  Mark the given node as already visited.
         *  @param AID
	 */
	void markAsPushed(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** @brief  Mark the given edge as being in a quad.
         *  @param AID an edge ID
	 */
	void markAsInQuad(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** \brief Update the graph connectivity.
	      *  @param
	 */
	void updateConnectivity();

	/*-------------------------------------------------------------------------*/
	/** \brief Get the leaving edges.
	      *  @param AN a node
	 */
	std::vector<Edge> getLeavingEdges(Node &AN);

	/*-------------------------------------------------------------------------*/
	/** \brief Get the entering edges.
	      *  @param AN a node
	 */
	std::vector<Edge> getEnteringEdges(Node &AN);

	/*-------------------------------------------------------------------------*/
	/** \brief Propagate length from a given root.
	      *  @param AID a node ID
	 */
	void propagateFromRoot(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** \brief Propagate length from roots.
	      *  @param
	 */
	void propagateFromRoots();

	/*-------------------------------------------------------------------------*/
	/** \brief Add 1 to the solstd::vector<TCellID> shortestPath;ution at every node of the cycle containing the given node.
	      *  @param AID a node ID
	 */
	void addOnCycle(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** \brief Propagate length from roots.
	      *  @param
	 */
	void addOnCycles();

	/*-------------------------------------------------------------------------*/
	/** \brief Build the complete quantization solution (minimizing the number of zeros).
	      *  @param
	 */
	void buildCompleteSolution();

	/*-------------------------------------------------------------------------*/
	/** @brief Display the quantization solution.
         *  @param
	 */
	void displaySolution();

	/*-------------------------------------------------------------------------*/
	/** @brief Get the value of the solution at the given node.
         *  @param AID a node ID
	 */
	int quantizationSolutionValue(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** @brief Return the shortest path joining the given node to a sink.
         *  @param AID a node ID
	 */
	std::vector<Node> shortestPathToASink(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** @brief Return the shortest path joining a source to the given node.
         *  @param AID a node ID
	 */
	std::vector<Node> shortestPathFromASource(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** @brief Return the shortest cycle containing the given node.
         *  @param AID a node ID
	 */
	std::vector<Node> shortestCycle(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** @brief Return the shortest elementary path containing the given node.
         *  @param AID a node ID
	 */
	std::vector<Node> shortestElementaryPath(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** @brief Return a cycle reachable from the given node.
         *  @param AID a node ID
	 */
	std::vector<Node> problematicCycle(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** @brief Get the edge corresponding to the two given nodes, in the right order.
         *  @param AN1, AN2 two nodes
	 */
	Edge getCorrespondingEdge(Node AN1, Node AN2);

	/*-------------------------------------------------------------------------*/
	/** @brief Return a quad edge approximatively in the middle of the given cycle. 
         *  @param AV a vector of nodes forming a cycle
	 */
	Edge middleQuadEdge(std::vector<Node> AV);

	/*-------------------------------------------------------------------------*/
	/** @brief Increase the quantization solution by 1 on the shortest elementary path containing the given node. The function return true if its succeeded
	 * increasing the solution, ie if the shortest elementary path containing the given node is not empty.
         *  @param AID a node ID
	 */
	bool increaseSolution(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** \brief Build minimal quantization solution (trying to minimize the number of blocks). The function returns couples of nodes 
	 * preventing from increasing the solution.
	      *  @param
	 */
	std::vector<std::vector<Node>> buildMinimalSolution();

	/*-------------------------------------------------------------------------*/
	/** \brief Check is for all half-edge e, l(e)=l(e.next.next). If not, artificially repair it. WARNING : this  solution doesn't garantee that the new solution is valid.
	      *  @param
	 */
	void roughlyRepairSolution();

	/*-------------------------------------------------------------------------*/
	/** \brief Check solution validity. The function returns couples of nodes 
	 * preventing from increasing the solution, in the cas where the solution is not valid.
	      *  @param
	 */
	std::vector<std::vector<Node>> checkSolutionValidity();



};
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/

#endif     // GMDS_QUANTIZATIONGRAPH_H
