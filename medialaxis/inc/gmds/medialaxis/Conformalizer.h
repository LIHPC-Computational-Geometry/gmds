// An object that build a conformal mesh given a non conformal mesh, its list of non conformal half-edges, and the length of its edges.

#ifndef GMDS_CONFORMALIZER_H
#define GMDS_CONFORMALIZER_H
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_MEDIALAXIS_export.h"
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
class LIB_GMDS_MEDIALAXIS_API Conformalizer
{
 private:
	// Non conformal quad block decomposition
	Mesh* m_mesh;
	// Corresponding non-conformal half edges
	std::vector<NonConformalHalfEdge> m_half_edges;
	// Half edges length
	std::vector<int> m_half_edges_lengths;
	// Groups of nodes of the non-conformal mesh
	std::vector<std::vector<TCellID>> m_non_conformal_nodes_groups;
	// Conformal mesh
	Mesh* m_conformal_mesh;
    // New nodes corresponding to half edges
	std::vector<std::vector<TCellID>> m_half_edges_to_new_nodes;
    // Intermediate mesh
	Mesh* m_intermediate_mesh;
	// Groups of nodes of the intermediate mesh
	std::vector<std::vector<TCellID>> m_intermediate_nodes_groups;

 public:

	/*-------------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param AMesh
	 */
	explicit Conformalizer(Mesh &AMesh, std::vector<NonConformalHalfEdge> AHalfEdges, std::vector<int> ALengths);

	/*-------------------------------------------------------------------------*/
	/** @brief Default destructor.
         *  @param
	 */
	virtual ~Conformalizer()=default;

    /*-------------------------------------------------------------------------*/
	/** @brief Write the intermediate mesh.
         *  @param
	 */
	void writeConformalMesh(std::basic_string<char> AFileName);

	/*-------------------------------------------------------------------------*/
	/** @brief Write the intermediate mesh.
         *  @param
	 */
	void writeIntermediateMesh(std::basic_string<char> AFileName);

	/*-------------------------------------------------------------------------*/
	/** @brief Get the conformal mesh.
         *  @param 
	 */
	Mesh getConformalMesh();

	/*-------------------------------------------------------------------------*/
	/** @brief Set the connectivity of the conformal mesh.
         *  @param 
	 */
	void setConformalMeshConnectivity();

	/*-------------------------------------------------------------------------*/
	/** @brief Set the connectivity of the intermediate mesh.
         *  @param 
	 */
	void setIntermediateMeshConnectivity();

	/*-------------------------------------------------------------------------*/
	/** @brief Build the groups of nodes of the non-conformal mesh.
         *  @param 
	 */
	void buildNonConformalNodesGroups();

	/*-------------------------------------------------------------------------*/
	/** @brief Add old nodes on the conformal mesh, each node corresponding to a group of old nodes.
         *  @param 
	 */
	void addOldNodesOnConformalMesh();

    /*-------------------------------------------------------------------------*/
	/** @brief Adds on the conformal mesh nodes subdiving edges of the non conformal mesh, according to their length.
         *  @param 
	 */
	void addSubdividingNodes();

    /*-------------------------------------------------------------------------*/
	/** @brief Attach to each half-edge its corresponding new nodes.
         *  @param 
	 */
	void buildHalfEdges2newNodesConnectivity();

	/*-------------------------------------------------------------------------*/
	/** @brief Build intermediate mesh nodes.
         *  @param
	 */
	void buildIntermediateMeshNodes();

    /*-------------------------------------------------------------------------*/
	/** @brief Build intermediate mesh edges.
         *  @param
	 */
	void buildIntermediateMeshEdges();

	/*-------------------------------------------------------------------------*/
	/** @brief Build groups of intermediate nodes.
         *  @param
	 */
	void buildIntermediateNodesGroups();

	/*-------------------------------------------------------------------------*/
	/** @brief Builds one node per group of intermediate nodes, and ensure previous node/new node correspondance.
         *  @param
	 */
	void buildRepresentativeNodes();

    /*-------------------------------------------------------------------------*/
	/** @brief Build the internal nodes of the conformal mesh, and the face/new nodes connectivity.
         *  @param 
	 */
	void builIntNodesAndFaces2newNodesConnectivity();

    /*-------------------------------------------------------------------------*/
	/** @brief Build conformal mesh faces.
         *  @param
	 */
	void buildConformalMeshFaces();

	/*-------------------------------------------------------------------------*/
	/** @brief Deletes superfluous nodes.
         *  @param
	 */
	void deleteSuperfluousNodes();

	/*-------------------------------------------------------------------------*/
	/** @brief Execute the algorithm building the conformal mesh.
         *  @param
	 */
	void execute();

	/*-------------------------------------------------------------------------*/
	/** @brief Smooth the conformal mesh.
         *  @param
	 */
	void smooth();
};
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_CONFORMALIZER_H