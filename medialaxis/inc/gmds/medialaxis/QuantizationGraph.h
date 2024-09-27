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
	/** \brief Add 1 to the solution at every node of the cycle containing the given node.
	      *  @param AID a node ID
	 */
	void addOnCycle(TCellID AID);

	/*-------------------------------------------------------------------------*/
	/** \brief Propagate length from roots.
	      *  @param
	 */
	void addOnCycles();

	/*-------------------------------------------------------------------------*/
	/** \brief Propagate length from roots.
	      *  @param
	 */
	void buildQuantizationSolution();

};
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/

#endif     // GMDS_QUANTIZATIONGRAPH_H
