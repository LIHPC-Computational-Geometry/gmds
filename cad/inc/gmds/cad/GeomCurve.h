/*----------------------------------------------------------------------------*/
/** \file    GeomCurve.h
 *  \author  F. LEDOUX
 *  \date    09/21/2010
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOMCURVE_H
#define GMDS_GEOMCURVE_H
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "GMDSCad_export.h"
#include <gmds/cad/GeomEntity.h>
#include <gmds/math/Segment.h>
#include <gmds/utils/CommonTypes.h>
#include <gmds/utils/Exception.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace cad {
class GeomPoint;
class GeomSurface;
class GeomVolume;
/*----------------------------------------------------------------------------*/
/** \class GeomCurve
 *  \brief This class describe the services that are required by the
 *  	   mesh to the geometrical model. As a consequence, this interface only
 *  	   contains query methods.
 */
/*----------------------------------------------------------------------------*/
class GMDSCad_API GeomCurve : public GeomEntity
{
 public:
	/*---------------------------------------------------------------*/
	/** @brief  a curve can be convex, concave or smooth relatively to
	 *          a volume
	 */

	enum CurvatureInfo { Flat, Convex, Concave, SmoothConvex, SmoothConcave };
	/*---------------------------------------------------------------*/
	/** @brief  Constructor
	 */
	explicit GeomCurve(const std::string &AName = "Unknown curve") : GeomEntity(AName) {}

	/*------------------------------------------------------------------------*/
	/** @brief  provides the dimension of the geometrical entity.
	 */
	int dim() const override
	{
		return 1;
	}

	/*------------------------------------------------------------------------*/
	/** @brief  Length of the curve
	 */
	virtual double length() const = 0;

	/*------------------------------------------------------------------------*/
	/** \brief  Compute the dihedral angle (max of each edge of the curve)
	 *
	 *  \return the dihedral angle
	 */
	virtual TCoord computeDihedralAngle() const = 0;

	/** @brief  Compute the tangent vector at the end point of a curve. The vector is oriented towards the curve
	 * 			If the curve is a loop, the starting point (param=0) is
	 * 			at param 1 too.
	 *  @param[in] AParam parameter equals to 0 or 1
	 *  @return a unit tangent vector at point of param @p AParam
	 */
	virtual math::Vector3d computeTangent(const int AParam) const = 0;

	/*------------------------------------------------------------------------*/
	/** \brief Move a point AP near the surface to the closest point on the
	 * 		   surface.
	 *  \param AP
	 */
	void project(math::Point &AP) const override = 0;

	/** @brief  Return whether the curve is a loop or not
	 *  @return a boolean
	 */
	virtual bool isALoop() const = 0;

	/**@brief compute the curvature info of the curve
	 * @return the curvature type of this
	 */
	virtual CurvatureInfo getCurvatureInfo() const = 0;
	/**@brief Accessor to the adjacent curves. Warning, there is no
	 *  assumption about the ordering
	 * @return curves that are adjacent to this point
	 */
	virtual std::vector<GeomPoint *> &points() = 0;
	/**@brief Accessor to the adjacent surfaces. Warning, there is no
	 *  assumption about the ordering
	 * @return surfaces that are adjacent to this point
	 */
	virtual std::vector<GeomSurface *> &surfaces() = 0;
	/**@brief Accessor to the adjacent volumes. Warning, there is no
	 *  assumption about the ordering
	 * @return volumes that are adjacent to this point
	 */
	virtual std::vector<GeomVolume *> &volumes() = 0;
};
/*----------------------------------------------------------------------------*/
}     // namespace cad
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GEOMCURVE_H */
