//
// Created by chenyt on 25/04/24.
//

#ifndef GMDS_MEDIALAXIS3DBUILDER_H
#define GMDS_MEDIALAXIS3DBUILDER_H
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_MEDIALAXIS_export.h"
#include "gmds/medialaxis/MedialAxis3D.h"
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace medialaxis{
/*----------------------------------------------------------------------------*/
/** \class  dummy
 *  \brief  dummy class.
 */
class LIB_GMDS_MEDIALAXIS_API MedialAxis3DBuilder{

 public:
	/*------------------------------------------------------------------------*/
	/** \brief friend class to access to protected/private data
	 */
	friend class Mesh;
	friend class MedialAxis3D;

	/*-------------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS
	} STATUS;


	/*-------------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param
	 */
	explicit MedialAxis3DBuilder(Mesh &AMesh);

 public:
	Mesh* m_mesh;
	MedialAxis3D* m_ma;


	/*-------------------------------------------------------------------------*/
	/** @brief Default destructor.
         *  @param
	 */
	virtual ~MedialAxis3DBuilder();

	/*-------------------------------------------------------------------------*/
	/** \brief Get the medial axis
	 */
	MedialAxis3D* getMedialObject();

	/*-------------------------------------------------------------------------*/
	/** \brief Get the mesh
	 */
	Mesh* getMesh();

	/*-------------------------------------------------------------------------*/
	/** \brief Check if an face of m_mesh is interior
	 */
	bool isAnInternalFace(const TCellID &AFaceID);

	/*-------------------------------------------------------------------------*/
	/** \brief Check if an edge of m_mesh is interior
	 */
	bool isAnInternalEdge(const TCellID &AEdgeID);

	/*-------------------------------------------------------------------------*/
	/** \brief Check if two tetras of m_mesh share a face
	 *     \param Two tetras' IDs
	 */
	bool shareFace(const TCellID &ATetraID1, const TCellID &ATetraID2);

	/*-------------------------------------------------------------------------*/
	/** \brief Check if the given tetrahedralization m_mesh is indeed of Delaunay.
	 * Warning, this function has quadratic complexity
	 *     \param
	 */
	bool isDelaunay();

	/*-------------------------------------------------------------------------*/
	/** \brief Returns the problematic nodes of the input minimal Delaunay mesh m_mesh, ie internal nodes (or Steiner points)
	 *     \param
	 */
	std::vector<Node> minDelProblematicPoints();

	/*-------------------------------------------------------------------------*/
	/** \brief Remove the problematic nodes returned by minDelProblematicPoints() (and all the cells that contain them)
	 *     \param
	 */
	void removeProblematicNodes();

	/*-------------------------------------------------------------------------*/
	/** \brief Remove the nearly flat tetras from m_mesh
	 *     \param
	 */
	void markNearlyFlatTetras();

	/*-------------------------------------------------------------------------*/
	/** \brief Deal with nearly flat tetras by moving their medial point to a reliable position
	 *     \param
	 */
	void dealWithNearlyFlatTetras();

	/*-------------------------------------------------------------------------*/
	/** \brief Computes the mesh maximal step
	 *     \param
	 */
	double meshStep();

	/*-------------------------------------------------------------------------*/
	/** \brief Regroups tetras whose circumcenters are close up to a tolerance
	 *     \param A tolerance ATol
	 */
	std::vector<std::vector<TCellID>> tetraFilter(const double &ATol);

	/*-------------------------------------------------------------------------*/
	/** \brief Builds the medial points
	 *     \param
	 */
	void buildMedialPoints();

	/*-------------------------------------------------------------------------*/
	/** \brief Builds the medial edges
	 *     \param
	 */
	void buildMedialEdges();

	/*-------------------------------------------------------------------------*/
	/** \brief Builds the medial faces
	 *     \param
	 */
	void buildMedialFaces();

	/*-------------------------------------------------------------------------*/
	/** \brief Places the problematic medial points (those coming from a Steiner tetra) at the
	 * barycenter of its non problematic neighbours
	 *     \param
	 */
	void repositionProblematicMedialPoints();

	/*-------------------------------------------------------------------------*/
	/** \brief Set the medial radius on medial points
	 */
	void setMedialRadius();

	/*-------------------------------------------------------------------------*/
	/** \brief Mark the nodes of the input geometry considered as details
	 */
	void markGeometryDetailsOnPoints();

	/*-------------------------------------------------------------------------*/
	/** \brief Mark the faces of the input geometry considered as details
	 */
	void markGeometryDetailsOnFaces();

	/*-------------------------------------------------------------------------*/
	/** \brief Returns the maximal length of the Delaunay edges (gives a typical size for the domain)
	 */
	double maxDelEdgeLength();

	/*-------------------------------------------------------------------------*/
	/** \brief Identify and mark small details on the geometry using the medial axis
	 */
	void identifyDetails();


	/*-------------------------------------------------------------------------*/
	/** \brief
	 */
	MedialAxis3DBuilder::STATUS execute();

	/*-------------------------------------------------------------------------*/
	/** \brief
	 */
	MedialAxis3DBuilder::STATUS executeFilter();

};
/*----------------------------------------------------------------------------*/
}  // end namespace medialaxis
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/

#endif     // GMDS_MEDIALAXIS3DBUILDER_H
