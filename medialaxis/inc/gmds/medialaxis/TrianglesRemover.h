// An object that tries to remove the triangles of a given almost-quad mesh.

#ifndef GMDS_TRIANGLESREMOVER_H
#define GMDS_TRIANGLESREMOVER_H
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
class LIB_GMDS_MEDIALAXIS_API TrianglesRemover
{
 private:
	// Degenerated quad mesh to fix
	Mesh* m_degenerated_mesh;
    // T-mesh
	Mesh* m_t_mesh;
	// Final T-mesh
	Mesh* m_final_t_mesh;
    // Groups of nodes of the degenerate mesh : two nodes belong to the same group iff they are both triangle vertices of the same fan
    std::vector<std::vector<Node>> m_degenerated_nodes_groups;

 public:

	/*-------------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param AMesh
	 */
	explicit TrianglesRemover(Mesh &AMesh);

	/*-------------------------------------------------------------------------*/
	/** @brief Default destructor.
         *  @param
	 */
	virtual ~TrianglesRemover()=default;

    /*-------------------------------------------------------------------------*/
	/** @brief Write the degenerated mesh.
         *  @param
	 */
	void writeDegeneratedMesh(std::basic_string<char> AFileName);

	/*-------------------------------------------------------------------------*/
	/** @brief Write the T-mesh.
         *  @param
	 */
	void writeTMesh(std::basic_string<char> AFileName);

    /*-------------------------------------------------------------------------*/
	/** @brief Faces composing the fan.
         *  @param AFanId 
	 */
	std::vector<Face> fan(int AFanId);

	/*-------------------------------------------------------------------------*/
	/** @brief Return an ordered vector of nodes forming the outline of the shape formed by the given faces.
         *  @param AV a vector of faces forming the shape
	 */
	std::vector<Node> outline(std::vector<Face> &AV);

    /*-------------------------------------------------------------------------*/
	/** @brief Build the T-mesh.
         *  @param 
	 */
	void buildTMesh();

    /*-------------------------------------------------------------------------*/
	/** @brief Build the degenerated nodes groups.
         *  @param 
	 */
	void buildDegeneratedNodesGroups();

    /*-------------------------------------------------------------------------*/
	/** @brief Sub-outlines of the given fan outline of type 1.
         *  @param AOutline, AFan, AFanId a vector of nodes forming the outline of a fan, the fan, and the Id of the fan
	 */
	std::vector<std::vector<Node>> subOutlines1(std::vector<Node> AOutline, std::vector<Face> AFan, int AFanId);

    /*-------------------------------------------------------------------------*/
	/** @brief Sub-outlines of the given fan outline of type 2.
         *  @param AOutline, AFan, AFanId a vector of nodes forming the outline of a fan, the fan, and the Id of the fan
	 */
	std::vector<std::vector<Node>> subOutlines2(std::vector<Node> AOutline, std::vector<Face> AFan, int AFanId);

    /*-------------------------------------------------------------------------*/
	/** @brief Get the T-mesh.
         *  @param 
	 */
	Mesh getTMesh();

    /*-------------------------------------------------------------------------*/
	/** @brief Set the connectivity of the T-mesh.
         *  @param 
	 */
	void setTMeshConnectivity();

    /*-------------------------------------------------------------------------*/
	/** @brief Mark internal constraints on edges.
         *  @param 
	 */
	void markInternalConstraintsOnEdges();

	/*-------------------------------------------------------------------------*/
	/** @brief Transforms the unwanted triangles into quads.
         *  @param 
	 */
	void transformUnwantedTrianglesIntoQuads();

	/*-------------------------------------------------------------------------*/
	/** @brief Build the final T-mesh.
         *  @param 
	 */
	void buildFinalTMesh();

	/*-------------------------------------------------------------------------*/
	/** @brief Write the final T-mesh.
         *  @param
	 */
	void writeFinalTMesh(std::basic_string<char> AFileName);

	/*-------------------------------------------------------------------------*/
	/** @brief Get the final T-mesh.
         *  @param 
	 */
	Mesh getFinalTMesh();

    /*-------------------------------------------------------------------------*/
	/** @brief Set the connectivity of the final T-mesh.
         *  @param 
	 */
	void setFinalTMeshConnectivity();

	/*-------------------------------------------------------------------------*/
	/** @brief Mark internal constraints on edges of the final T-mesh.
         *  @param 
	 */
	void markInternalConstraintsOnFinalEdges();

};
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_TRIANGLESREMOVER_H