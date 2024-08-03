/*----------------------------------------------------------------------------*/
/** \file    FacetedCurve.h
 *  \author  F. LEDOUX
 *  \date    30/05/2011
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_FACETEDCURVE_H_
#	define GMDS_GEOM_FACETEDCURVE_H_
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "GMDSCadFac_export.h"
#include "gmds/cad/GeomCurve.h"
#include "gmds/cadfac/FACPoint.h"
#include "gmds/ig/Edge.h"
#include "gmds/ig/Node.h"
#include "gmds/math/Point.h"
#include "gmds/math/Vector.h"
#include "gmds/utils/CommonTypes.h"
#include "gmds/utils/Exception.h"
/*----------------------------------------------------------------------------*/
#include <tuple>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace cad {
/*----------------------------------------------------------------------------*/
class FACSurface;
/*----------------------------------------------------------------------------*/
/** \class FACCurve
 *  \brief This class implements the curve services that are required by the
 *  	   mesh to the geometrical model.
 */
/*----------------------------------------------------------------------------*/
class GMDSCadFac_API FACCurve : public GeomCurve
{

 public:
	/*------------------------------------------------------------------------*/
	/** @brief  Default Constructor
	 *  @param AMeshSupport Mesh support for the faceted model
	 */
	FACCurve(Mesh *AMeshSupport);

	/*------------------------------------------------------------------------*/
	/** \brief  Constructor. A geometric curve is built as an ordered
	 * 			collection of points.
	 * 			The points vector must contain
	 * 		 	number of segments + 1 nodes : when the curve is a loop,
	 * 			the last point is also the first one.
	 *  \param AP1 		the first end point
	 *  \param AP2 		the second end point
	 *  \param ANodes 	the nodes used to discretize the curve
	 *  				(including end points).
	 *  \param AEdges 	the edges used to discretize the curve
	 *  				(including end points).
	 */
	FACCurve(Mesh *AMeshSupport, std::vector<TCellID> &ANodes, std::vector<TCellID> &AEdges, const std::string &AName = "Unknown curve");

	/*------------------------------------------------------------------------*/
	/** \brief  Constructor. A geometric curve is built as an ordered
	 * 			collection of points.
	 * 			The points vector must contain
	 * 		 	number of segments + 1 nodes : when the curve is a loop,
	 * 			the last point is also the first one.
	 *  \param AP1 		the first end point
	 *  \param AP2 		the second end point
	 *  \param APoints 	the node used to discretize the curve
	 *  				(including end points).
	 */
	FACCurve(Mesh *AMeshSupport,
	         FACPoint *AP1,
	         FACPoint *AP2,
	         std::vector<gmds::TCellID> &APoints,
	         std::vector<gmds::TCellID> &AEdges,
	         const std::string &AName = "Unknown curve");

	/*------------------------------------------------------------------------*/
	/** \brief  Destructor
	 */
	virtual ~FACCurve();

	/*------------------------------------------------------------------------*/
	/** \brief  Length of the curve
	 */
	double length() const override;

	/*------------------------------------------------------------------------*/
	/** \brief Move a point AP near the surface to the closest point on the
	 * 		   surface.
	 *  \param AP
	 */
	void project(math::Point &AP) const override;

	/*------------------------------------------------------------------------*/
	/** \brief Move a point AP near the surface to the closest point on the
	 * 		   surface. Also fills the tangent vector.
	 *  \param AP
	 *  \param AV the tangent vector
	 */
	void project(math::Point &AP, math::Vector3d &AV) const override;

	/*------------------------------------------------------------------------*/
	/** \brief Get the closest point from AP on the surface
	 *  \param AP a 3D point
	 *
	 *  \return the closest point of APoint on the surface
	 */
	math::Point closestPoint(const math::Point &AP) const override;

	/*------------------------------------------------------------------------*/
	/** \brief  computes the area of the entity.
	 */
	TCoord computeArea() const override;

	/*------------------------------------------------------------------------*/
	/** \brief  computes the bounding box
	 *
	 *	\param minXYZ The minimum coordinate of the bounding box.
	 *	\param maxXYZ The maximum coordinate of the bounding box.
	 */
	void computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const override;
	std::tuple<TCoord, TCoord, TCoord, TCoord, TCoord, TCoord> BBox() const override;

	/** \brief  Return whether the curve is a loop or not
	 *  \return a boolean
	 */
	bool isALoop() const override;

	/**@brief compute the curvature info of the curve
	 * @return the curvature type of this
	 */
	CurvatureInfo getCurvatureInfo() const override;

	/** @brief  Compute the tangent vector at the end point of a curve. The vector is oriented towards the curve
	 * 			If the curve is a loop, the starting point (param=0) is
	 * 			at param 1 too.
	 *  @param[in] AParam parameter equals to 0 or 1
	 *  @return a unit tangent vector at point of param @p AParam
	 */
	math::Vector3d computeTangent(const int AParam) const override;


	/** \brief  Compute the dihedral angle (max for each edge of the curve)
	 *
	 *  \return the dihedral angle
	 */
	TCoord computeDihedralAngle() const override;

	/*------------------------------------------------------------------------*/
	/** \brief  Returns a copy of the internal mesh representation nodes
	 *
	 *  \param ANodes a vector of mesh nodes
	 */
	void getMeshNodes(std::vector<Node> &ANodes) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Sets the internal mesh representation edges
	 *
	 *  \param AEdges a vector of mesh edges
	 */
	void setMeshEdges(const std::vector<Edge> &AEdges);

	/*------------------------------------------------------------------------*/
	/** \brief  Returns a copy of the internal mesh representation edges
	 *
	 *  \param AEdges a vector of mesh edges
	 */
	void getMeshEdges(std::vector<Edge> &AEdges) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Returns a triangulation of the curve
	 *
	 *  \param ASeg a triangulation
	 */
	virtual void getTriangulation(std::vector<gmds::math::Segment> &ASeg) const;

	virtual int id() const
	{
		return m_id;
	}

	/**@brief Accessor to the adjacent curves. Warning, there is no
	 *  assumption about the ordering
	 * @return curves that are adjacent to this point
	 */
	virtual std::vector<GeomPoint *> &points();
	/**@brief Accessor to the adjacent surfaces. Warning, there is no
	 *  assumption about the ordering
	 * @return surfaces that are adjacent to this point
	 */
	virtual std::vector<GeomSurface *> &surfaces();
	/**@brief Accessor to the adjacent volumes. Warning, there is no
	 *  assumption about the ordering
	 * @return volumes that are adjacent to this point
	 */
	virtual std::vector<GeomVolume *> &volumes();

	/**@brief Reset the global id counter to 1.
	 */
	static void resetIdCounter();

 private:
	Mesh *m_support;

	int m_id;
	static int m_next_id;

	std::vector<TCellID> m_mesh_nodes;
	std::vector<TCellID> m_mesh_edges;

	/** adjacent geometric points that bound this curves*/
	std::vector<GeomPoint *> m_adjacent_points;
	/** surfaces adjacent to this curve*/
	std::vector<GeomSurface *> m_adjacent_surfaces;
	/** volumes adjacent to this surface*/
	std::vector<GeomVolume *> m_adjacent_volumes;
};
/*----------------------------------------------------------------------------*/
}     // namespace cad
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GEOM_FACETEDCURVE_H_ */
/*----------------------------------------------------------------------------*/
