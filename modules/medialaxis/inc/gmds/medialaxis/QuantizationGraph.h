//
// Created by chenyt on 27/09/24.
//

#ifndef GMDS_QUANTIZATIONGRAPH_H
#define GMDS_QUANTIZATIONGRAPH_H
/*----------------------------------------------------------------------------*/
#include "GMDSMedialaxis_export.h"
#include "gmds/medialaxis/MedialAxisMath.h"
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/** \class  dummy
 *  \brief  dummy class.
 */
class GMDSMedialaxis_API QuantizationGraph
{
 private:
	// Mesh representation
	Mesh* m_mesh_representation;

	// Non zero conditions for the quantization solution :

	// Set of verticies with non zero condition
	std::vector<TCellID> m_non_zero_verticies;

	// Set of groups of verticies, with condition at least one non zero in each group
	std::vector<std::vector<TCellID>> m_non_zero_groups;

	// Maximal number of iterations to find a path
	int m_max_it;

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
	/** @brief  Mark the given node as already pushed.
         *  @param AID
	 */
	void markAsPushed(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** @brief  Mark the given edge as being in a quad.
         *  @param AID an edge ID
	 */
	void markAsInQuad(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** @brief  Attach to the node the geometrical length of the half edge it represents.
         *  @param AID, ALength a node ID and a length
	 */
	void setGeometricalLength(TCellID AID, double ALength);

	/*-------------------------------------------------------------------------*/
	/** @brief  Mark the node as zero of the quantization solution.
         *  @param AID a node ID
	 */
	void markAsZero(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** \brief Update the graph connectivity.
	      *  @param
	 */
	void updateConnectivity();

	/*-------------------------------------------------------------------------*/
	/** \brief Get the leaving edges. ========= TO DELETE IN THE END =========
	      *  @param AN a node
	 */
	std::vector<Edge> getLeavingEdges(Node &AN);

	/*-------------------------------------------------------------------------*/
	/** \brief Get the entering edges. ========= TO DELETE IN THE END =========
	      *  @param AN a node
	 */
	std::vector<Edge> getEnteringEdges(Node &AN);

	/*-------------------------------------------------------------------------*/
	/** \brief Get the unique edge interior to a quad containing the given node.
	      *  @param AN a node
	 */
	Edge getIntEdge(Node &AN);

	/*-------------------------------------------------------------------------*/
	/** \brief Get set of edges between different quads containing the given node.
	      *  @param AN a node
	 */
	std::vector<Edge> getExtEdges(Node &AN);

	/*-------------------------------------------------------------------------*/
	/** \brief Given a node and an edge containing this node, return the other node of the edge.
	      *  @param AN, AE
	 */
	Node getOtherNode(Node AN, Edge AE);

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
	/** @brief Get the value of the flux through the given edge.
         *  @param AID a edge ID
	 */
	int fluxValue(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** @brief Return the shortest path joining the given node to a root, starting in the given direction (if 1 through the quad, if -1 else). 
         *  @param AID a node ID
	 */
	std::vector<Node> shortestHalfPath(TCellID AID, int ADirection);

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
	/** @brief Return a quad edge approximatively in the middle of the given cycle. The edge is given by a set of 2 nodes which orientates the edge.
         *  @param AV a vector of nodes forming a cycle
	 */
	std::vector<Node> middleQuadEdge(std::vector<Node> AV);

	/*-------------------------------------------------------------------------*/
	/** @brief Increase the quantization solution by 1 on the shortest elementary path containing the given node. The function return true if it succeeded
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

	/*-------------------------------------------------------------------------*/
	/** \brief Return the optimal (according to the input mesh size) path joining the given node to a root, starting in the given direction (if 1 through the quad, if -1 else). 
	      *  @param AID, ADirection a node ID and a direction.
	 */
	std::vector<Node> optimalHalfPath(TCellID AID, int ADirection);

	/*-------------------------------------------------------------------------*/
	/** \brief Return the optimal (according to the input mesh size) path joining two roots and containing the given node. 
	      *  @param AID a node ID.
	 */
	std::vector<Node> optimalPath(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** \brief Return the optimal (according to the input mesh size) cycle containing the given node. 
	      *  @param AID a node ID.
	 */
	std::vector<Node> optimalCycle(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** \brief Return the optimal (according to the input mesh size) elementary path containing the given node. NOT USED
	      *  @param AID a node ID.
	 */
	std::vector<Node> optimalElementaryPath(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** \brief Imroves the complete solution with respect to the target mesh size.
	      *  @param AMeshSize a target mesh size.
	 */
	void improveSolution(double AMeshSize);

};
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/

#endif     // GMDS_QUANTIZATIONGRAPH_H
