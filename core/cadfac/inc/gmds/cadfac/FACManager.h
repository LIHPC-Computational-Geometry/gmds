/*----------------------------------------------------------------------------*/
/*
 * FACManager.h
 *
 *  Created on: 1 juil. 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_FACETEDGEOMMANAGER_H_
#define GMDS_GEOM_FACETEDGEOMMANAGER_H_
/*----------------------------------------------------------------------------*/
#include <stack>
/*----------------------------------------------------------------------------*/
// gmds File headers
/*----------------------------------------------------------------------------*/
#include "gmds/utils/CommonTypes.h"
#include "gmds/utils/Exception.h"

#include "gmds/math/Vector.h"

#include "gmds/cad/GeomCurve.h"
#include "gmds/cad/GeomManager.h"
#include "gmds/cadfac/FACCurve.h"
#include "gmds/cadfac/FACPoint.h"
#include "gmds/cadfac/FACSurface.h"
#include "gmds/cadfac/FACVolume.h"

#include "GMDSCadFac_export.h"
#include "gmds/cad/GeomMeshLinker.h"
#include "gmds/ig/Edge.h"
#include "gmds/ig/Mesh.h"
#include "gmds/ig/MeshDoctor.h"
#include "gmds/ig/Node.h"

/*----------------------------------------------------------------------------*/
// avoid #include <gts.h> here by using forward declaration for
// the GTS data structure; also requires GNode of the glib
struct _GtsSurface;
typedef _GtsSurface GtsSurface;
struct _GNode;
typedef _GNode GNode;
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace cad {
/*----------------------------------------------------------------------------*/
/** \class  FACManager
 *  \brief  This class creates all the model entities and services relative to
 *  		the faceted representation of the geometry.
 *
 *
 */
/*----------------------------------------------------------------------------*/
class GMDSCadFac_API FACManager : public GeomManager
{

 public:
	/*------------------------------------------------------------------------*/
	/** \brief  Default constructor
	 */
	FACManager();

	/*------------------------------------------------------------------------*/
	/** \brief  Destructor
	 */
	virtual ~FACManager();

	/*------------------------------------------------------------------------*/
	/** @brief create a faceted model from a 3D tet mesh. There is no assumption
	 *         about AFromMesh, we get its boundary to build the model.
	 *
	 *         Warning, it means that we totally rebuild the model and we
	 *         consider that it is made of only one single volume!
	 *
	 *  @param AFromMesh 3D mesh used to build the model
	 */
	void initFrom3DMesh(Mesh *AFromMesh);

	/*------------------------------------------------------------------------*/
	/** @brief create a faceted model from a 3D tet mesh. There is no assumption
	 *         about AFromMesh, we get its boundary to build the model. A linker
	 *         given in arguments is filled up to keep the classification.
	 *
	 *         Warning, it means that we totally rebuild the model and we
	 *         consider that it is made of only one single volume!
	 *
	 *  @param AFromMesh 3D mesh used to build the model
	 */
	void initAndLinkFrom3DMesh(Mesh *AFromMesh, GeomMeshLinker *ALinker);

	/*------------------------------------------------------------------------*/
	/** @brief create a faceted model from a 2D mesh. There is no assumption
	 *         about AFromMesh, we get its boundary to build the model. A linker
	 *         given in arguments is filled up to keep the classification.
	 *
	 *         Warning, it means that we totally rebuild the model and we
	 *         consider that it is made of only one single volume!
	 *
	 *  @param AFromMesh 2D mesh used to build the model
	 */
	void initAndLinkFrom2DMesh(Mesh *AFromMesh, GeomMeshLinker *ALinker);

	/*------------------------------------------------------------------------*/
	/** \brief  reinitializes the complete model from data stored in mesh_.
	 */
	void updateFromMesh();

	/*------------------------------------------------------------------------*/
	/** @brief reinitialize the geom model
	 */
	void clear();
	/*------------------------------------------------------------------------*/
	/** \brief  creation of a geometric volume
	 */
	virtual GeomVolume *newVolume();

	/*------------------------------------------------------------------------*/
	/** \brief  creation of a geometric surface
	 */
	virtual GeomSurface *newSurface();

	/*------------------------------------------------------------------------*/
	/** \brief  creation of a geometric curve
	 */
	virtual GeomCurve *newCurve();

	/*------------------------------------------------------------------------*/
	/** \brief  creation of a geometric point
	 */
	virtual GeomPoint *newPoint();

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of points of the model.
	 *
	 *	\return the number of points.
	 */
	TInt getNbPoints() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of curves of the model.
	 *
	 *	\return the number of curves.
	 */
	TInt getNbCurves() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of surfaces of the model.
	 *
	 *	\return the number of surfaces.
	 */
	TInt getNbSurfaces() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of volumes of the model.
	 *
	 *	\return the number of volumes.
	 */
	TInt getNbVolumes() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the points of the model.
	 *
	 *  \param points the points of the model.
	 */
	void getPoints(std::vector<GeomPoint *> &points) const;
	std::vector<GeomPoint*> getPoints() const;
	/*------------------------------------------------------------------------*/
	/** \brief  Access to the curves of the model.
	 *
	 *  \param curves the curves of the model.
	 */
	void getCurves(std::vector<GeomCurve *> &curves) const;
	std::vector<GeomCurve *> getCurves() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the surface of the model.
	 *
	 *  \param surfaces the surfaces of the model.
	 */
	void getSurfaces(std::vector<GeomSurface *> &surfaces) const;
	std::vector<GeomSurface *> getSurfaces() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the volumes of the model.
	 *
	 *  \param volumes the volumes of the model.
	 */
	void getVolumes(std::vector<GeomVolume *> &volumes) const;
	std::vector<GeomVolume *> getVolumes() const;


    /*------------------------------------------------------------------------*/
    /** \brief  Gives access to the entity of id @AID and dimension @p ADim.
     *          Return NullPtr if it does not exist.
     *  \return A point
     */
    virtual GeomEntity* getEntity(TInt AID, TInt ADim);
	/*------------------------------------------------------------------------*/
	/** \brief  Gives access to the point of id @AID, Return NullPtr if it does
	 *          not exist.
	 *  \return A point
	 */
	virtual GeomPoint *getPoint(const TInt AID);
	/*------------------------------------------------------------------------*/
	/** \brief  Gives access to the curve of id @AID, Return NullPtr if it does
	 *          not exist.
	 *  \return A curve
	 */
	virtual GeomCurve *getCurve(const TInt AID);
	/*------------------------------------------------------------------------*/
	/** \brief  Gives access to the surface of id @AID, Return NullPtr if it
	 *          does not exist.
	 *  \return A surface
	 */
	virtual GeomSurface *getSurface(const TInt AID);
	/*------------------------------------------------------------------------*/
	/** \brief  Gives access to the volume of id @AID, Return NullPtr if it
	 *          does not exist.
	 *  \return A volume
	 */
	virtual GeomVolume *getVolume(const TInt AID);
	/*------------------------------------------------------------------------*/
	/** \brief  Gives access to the inner mesh view of the model.
	 *
	 *  Do not modify the referenced mesh content
	 */
	Mesh &getMeshView();
	const Mesh *getMeshView_ptr() const;
	/*------------------------------------------------------------------------*/
	/** \brief  Get the curve common to 2 points
	 *
	 *  \param return the id of the common curve, and -1 if it doesn't exist
	 */
	int getCommonCurve(GeomPoint *AP1, GeomPoint *AP2) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the surfaces common to 2 points
	 *
	 *  \param return the ids of the common surfaces (potentially empty)
	 */
	std::vector<int> getCommonSurfaces(GeomPoint *AP1, GeomPoint *AP2) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the surface common to 2 curves
	 *
	 *  \param return the id of the common surface, and -1 if it doesn't exist
	 */
	int getCommonSurface(GeomCurve *AC1, GeomCurve *AC2) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Write the mesh representation
	 */
	void write_surfaces(std::string AFilename) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Temporary method to get a FACSurface in order to have to geometry-to-tetmesh connectivity (Simon)
	 */
	FACSurface *getFACSurface(const TInt AID);

	/**@brief Build a gts tree data structure to accelerate information retrieval
    *        in the volume. It will consider the GTS surface formed of
    *        the triangular faces of the surfaces adjacent to the volume
    *
    *        Warning, we consider that it is made of only one single volume!
	 */
	void buildGTSTree(const Mesh *AFromMesh);

	/**@brief Check whether a point is inside the model
	 *        Warning, requires that the search tree was previously built
	 *
	 * @param APt a point
    * @return whether the point is inside the model
	 */
	bool is_in(gmds::math::Point APt) const;

 private:
	/*------------------------------------------------------------------------*/
	/** \brief  Build topological relations between geom entities. Must be done
	 *          in the init functions only!
	 */
	void buildTopologicalConnectivity(const int dim = 3);

	void OutwardNormal(Face *f);
	bool IsInsideFace(const int &, const math::Vector3d &, const math::Vector3d &, Face fi);

	math::Vector3d GetBarycentricCoefs(math::Vector3d &V1, math::Vector3d &V2, math::Vector3d &V3, math::Vector3d &I, math::Vector3d &normal);

	Mesh m_mesh;

	std::vector<FACVolume *> m_volumes;
	std::vector<FACSurface *> m_surfaces;
	std::vector<FACCurve *> m_curves;
	std::vector<FACPoint *> m_points;

	std::map<TInt, TInt> m_map_node_var_to_point;
	std::map<TInt, TInt> m_map_edge_var_to_curve;
	std::map<TInt, TInt> m_map_face_var_to_surf;

	// gmds and GTS data structure for fast retrieval
	Mesh m_mesh_is_inside;
	GNode* m_groot;
};
/*----------------------------------------------------------------------------*/
}     // namespace cad
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* gmds_GEOM_FACETEDGEOMMANAGER_H_ */
/*----------------------------------------------------------------------------*/
