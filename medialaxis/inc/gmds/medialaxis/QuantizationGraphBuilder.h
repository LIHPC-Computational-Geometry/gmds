//
// Created by chenyt on 26/09/24.
//

#ifndef GMDS_QUANTIZATIONGRAPHBUILDER_H
#define GMDS_QUANTIZATIONGRAPHBUILDER_H
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_MEDIALAXIS_export.h"
#include "gmds/medialaxis/MedialAxis2D.h"
#include "gmds/medialaxis/NonConformalHalfEdge.h"
#include "gmds/medialaxis/MedialAxisMath.h"
#include "gmds/medialaxis/QuantizationGraph.h"
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/** \class  dummy
 *  \brief  dummy class.
 */
class LIB_GMDS_MEDIALAXIS_API QuantizationGraphBuilder
{
 private:
	// Non conformal quad block decomposition
	Mesh* m_mesh;
	// Corresponding non-conformal half edges
	std::vector<NonConformalHalfEdge> m_half_edges;
	// Quantization graph
	QuantizationGraph* m_quantization_graph;

 public:

	/*-------------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param AMesh
	 */
	explicit QuantizationGraphBuilder(Mesh &AMesh);

	/*-------------------------------------------------------------------------*/
	/** @brief Default destructor.
         *  @param
	 */
	virtual ~QuantizationGraphBuilder()=default;

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
	void execute();

	/*-------------------------------------------------------------------------*/
	/** @brief Get the quantization graph.
         *  @param
	 */
	QuantizationGraph* getQuantizationGraph();
};
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_QUANTIZATIONGRAPHBUILDER_H
