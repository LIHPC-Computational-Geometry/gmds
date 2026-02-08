//
// Created by chenyt on 25/04/24.
//

#ifndef GMDS_MEDIALAXIS3D_H
#define GMDS_MEDIALAXIS3D_H
/*----------------------------------------------------------------------------*/
#include "GMDSMedialaxis_export.h"
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace medialaxis{
/*----------------------------------------------------------------------------*/
/** \class  dummy
 *  \brief  dummy class.
 */
class GMDSMedialaxis_API MedialAxis3D{

 public:
	/*-------------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param
	 */
	explicit MedialAxis3D();

	/*-------------------------------------------------------------------------*/
	/** @brief Default destructor.
         *  @param
	 */
	virtual ~MedialAxis3D();

	/*-------------------------------------------------------------------------*/
	/** \brief Add a medial point.
	      *  @param a point.
	 */
	Node newMedPoint(const math::Point &APnt);

	/*-------------------------------------------------------------------------*/
	/** \brief Add a medial edge.
	      *  @param two medial points id.
	 */
	Edge newMedEdge(const TCellID &AN1, const TCellID &AN2);

	/*-------------------------------------------------------------------------*/
	/** \brief Add a medial face.
	      *  @param Anodes ordered vector of nodes IDs
	 */
	Face newMedFace(const std::vector<TCellID> &ANodes);

	/*-------------------------------------------------------------------------*/
	/** \brief Get a medial point.
	      *  @param a point ID.
	 */
	Node getMedPoint(const TCellID &AN);

	/*-------------------------------------------------------------------------*/
	/** \brief Delete medial point
	      *  @param a medial point ID
	 */
	void deleteMedPoint(const TCellID &AN);

	/*-------------------------------------------------------------------------*/
	/** \brief Delete medial edge
	      *  @param a medial edge ID
	 */
	void deleteMedEdge(const TCellID &AE);

	/*-------------------------------------------------------------------------*/
	/** \brief Delete medial face
	      *  @param a medial face ID
	 */
	void deleteMedFace(const TCellID &AF);

	/*-------------------------------------------------------------------------*/
	/** \brief Average the position of a point with respect to a vector of points
	      *  @param a medial point ID
	 */
	void averagePoint(const TCellID &AN, const std::vector<TCellID> ANodes);

	/*-------------------------------------------------------------------------*/
	/** \brief Update connectivity
	      *  @param
	 */
	void updateConnectivity();

	/*-------------------------------------------------------------------------*/
	/** \brief Set the medial point type
	      *  @param
	 */
	void setMedialPointType();

	/*-------------------------------------------------------------------------*/
	/** \brief Set the medial edge type
	      *  @param
	 */
	void setMedialEdgeType();

	/*-------------------------------------------------------------------------*/
	/** \brief Set the medial face type
	      *  @param
	 */
	void setMedialFaceType();

	/*-------------------------------------------------------------------------*/
	/** \brief Set the medial surfaces id on faces
	      *  @param
	 */
	int setMedialSurfaceIdOnFaces();

	/*-------------------------------------------------------------------------*/
	/** \brief Set the medial surfaces id on points
	      *  @param
	 */
	void setMedialSurfaceIdOnPoints();

	/*-------------------------------------------------------------------------*/
	/** \brief Mark the points and faces belonging to branches corresponding to details
	      *  @param
	 */
	void identifyDetailsSurfaces(const double &ATol);

	/*-------------------------------------------------------------------------*/
	/** \brief Returns the value of corresponds_to_detail of a medial point
	      *  @param APointID a point ID
	 */
	int correspondsToDetail(const TCellID & APointID);

	/*-------------------------------------------------------------------------*/
	/** \brief Set the medial surfaces type on faces
	      *  @param
	 */
	void setMedialSurfaceTypeOnFaces();

	/*-------------------------------------------------------------------------*/
	/** \brief Set the medial surfaces size on faces
	      *  @param
	 */
	void setMedialSurfaceSizeOnFaces();

	/*-------------------------------------------------------------------------*/
	/** \brief Set the medial radius at a given medial point
	      *  @param A medial point ID, a value
	 */
	void setMedialRadius(const TCellID &APointID, double AValue);

	/*-------------------------------------------------------------------------*/
	/** \brief Flatten the unwanted surfaces
	      *  @param
	 */
	void flattenUnwantedFaces();


	/*-------------------------------------------------------------------------*/
	/** \brief Write the medial axis
	      *  @param
	 */
	void write(std::string AFileName);

 private:

	Mesh* m_mesh_representation;
};
/*----------------------------------------------------------------------------*/
}  // end namespace medialaxis
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_MEDIALAXIS3D_H
