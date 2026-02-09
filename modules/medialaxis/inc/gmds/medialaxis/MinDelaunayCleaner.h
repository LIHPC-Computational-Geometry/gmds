#ifndef GMDS_MINDELAUNAYCLEANER_H
#define GMDS_MINDELAUNAYCLEANER_H
/*----------------------------------------------------------------------------*/
#include "GMDSMedialaxis_export.h"
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
class GMDSMedialaxis_API MinDelaunayCleaner
{
 private:
	// Minimal Delaunay triangulation to clean
	Mesh* m_mesh;
    // Cleaned mesh
	Mesh* m_cleaned_mesh;

 public:

	/*-------------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param AMesh
	 */
	explicit MinDelaunayCleaner(Mesh &AMesh);

	/*-------------------------------------------------------------------------*/
	/** @brief Default destructor.
         *  @param
	 */
	virtual ~MinDelaunayCleaner()=default;

    /*-------------------------------------------------------------------------*/
	/** @brief Compute the number of neighbours of each face.
         *  @param
	 */
	void setFacesTypes();

    /*-------------------------------------------------------------------------*/
	/** @brief Compute the medial points and radii.
         *  @param
	 */
	void computeMedialPointsAndRadii();

    /*-------------------------------------------------------------------------*/
	/** @brief Return a vector containing the faces of the medial branch of the input face.
         *  @param AF
	 */
	std::vector<Face> medialBranch(Face &AF);

    /*-------------------------------------------------------------------------*/
	/** @brief Set an ID to each medial branch.
         *  @param 
	 */
	void setBranchID();

    /*-------------------------------------------------------------------------*/
	/** @brief Mark small edges.
         *  @param 
	 */
	void markSmallEdges();

    /*-------------------------------------------------------------------------*/
	/** @brief Mark the faces to delete, those with too small area.
         *  @param 
	 */
	void markFacesToDelete();

    /*-------------------------------------------------------------------------*/
	/** @brief Build the cleaned mesh.
         *  @param 
	 */
	void buildCleanedMesh();

    /*-------------------------------------------------------------------------*/
	/** @brief Set the cleaned mesh connectivity.
         *  @param 
	 */
	void setCleanedMeshConnectivity();

    /*-------------------------------------------------------------------------*/
	/** @brief Get the cleaned mesh.
         *  @param 
	 */
	Mesh getCleanedMesh();

    /*-------------------------------------------------------------------------*/
	/** @brief Write the cleaned mesh.
         *  @param 
	 */
	void writeCleanedMesh(std::string AFileName);

};
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_MINDELAUNAYCLEANER_H