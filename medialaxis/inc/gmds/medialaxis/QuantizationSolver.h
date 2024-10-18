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
	// Quantized mesh
	Mesh* m_quantized_mesh;
	// New nodes corresponding to half edges
	std::vector<std::vector<TCellID>> m_half_edges_to_new_nodes;

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
	/** @brief Build the group of half edges.
         *  @param AID, AV a half edge ID and a vector to mark already seen half edges
	 */
	std::vector<std::vector<int>> halfEdgesGroup(int AID, std::vector<int> &AV);

	/*-------------------------------------------------------------------------*/
	/** @brief Build the groups of half edges.
         *  @param
	 */
	std::vector<std::vector<std::vector<int>>> halfEdgesGroups();

	/*-------------------------------------------------------------------------*/
	/** @brief Set half edges length.
         *  @param
	 */
	void setHalfEdgesLength();

	/*-------------------------------------------------------------------------*/
	/** @brief Write the quantized mesh.
         *  @param
	 */
	void writeQuantizedMesh(std::basic_string<char> AFileName);

	/*-------------------------------------------------------------------------*/
	/** @brief Build the nodes of the quantized mesh that are on edges of the non-conformal mesh.
         *  @param
	 */
	void buildQuantizedMeshNodesOnEdges();

	/*-------------------------------------------------------------------------*/
	/** @brief Build the nodes of the quantized mesh that are internal to the faces of the non-conformal mesh.
         *  @param
	 */
	void buildQuantizedMeshInternalNodes();

	/*-------------------------------------------------------------------------*/
	/** @brief Build quantized mesh faces.
         *  @param
	 */
	void buildQuantizedMeshFaces();
};
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_QUANTIZATIONSOLVER_H
