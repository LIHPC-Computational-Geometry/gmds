#ifndef GMDS_CONFORMALMESHBUILDER_H
#define GMDS_CONFORMALMESHBUILDER_H
/*----------------------------------------------------------------------------*/
#include "GMDSMedialaxis_export.h"
#include "gmds/medialaxis/NonConformalHalfEdge.h"
#include "gmds/medialaxis/MedialAxisMath.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKWriter.h"
#include "gmds/ig/MeshDoctor.h"
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/** \class  dummy
 *  \brief  dummy class.
 */
class GMDSMedialaxis_API ConformalMeshBuilder
{
 private:
	// Non conformal quad block decomposition
	Mesh* m_mesh;
	// Corresponding non-conformal half edges
	std::vector<NonConformalHalfEdge> m_half_edges;
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
	explicit ConformalMeshBuilder(Mesh &AMesh, std::vector<NonConformalHalfEdge> AHalfEdges, std::vector<int> ALengths);

	/*-------------------------------------------------------------------------*/
	/** @brief Default destructor.
         *  @param
	 */
	virtual ~ConformalMeshBuilder()=default;

	/*-------------------------------------------------------------------------*/
	/** @brief Write the quantized mesh.
         *  @param
	 */
	void writeQuantizedMesh(std::basic_string<char> AFileName);

	/*-------------------------------------------------------------------------*/
	/** @brief Get the quantized mesh.
         *  @param 
	 */
	Mesh getQuantizedMesh();

	/*-------------------------------------------------------------------------*/
	/** @brief Set the connectivity of the quantized mesh.
         *  @param 
	 */
	void setQuantizedMeshConnectivity();

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
	/** @brief Execute.
         *  @param
	 */
	void execute();

};
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_CONFORMALMESHBUILDER_H