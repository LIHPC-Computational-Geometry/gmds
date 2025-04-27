/*----------------------------------------------------------------------------*/
/** \file    GeomPoint.h
 *  \author  F. LEDOUX
 *  \date    02/08/2010
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_GEOMPOINT_H_
#define GMDS_GEOM_GEOMPOINT_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
/*----------------------------------------------------------------------------*/
#include "GMDSCad_export.h"
#include "gmds/cad/GeomEntity.h"
#include "gmds/utils/CommonTypes.h"
#include "gmds/utils/Exception.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace cad {
class GeomCurve;
class GeomSurface;
class GeomVolume;
/*----------------------------------------------------------------------------*/
/** \class GeomPoint
 *  \brief This class describe the services that are required by the
 *  	   mesh to the geometrical model. As a consequence, this interface only
 *  	   contains query methods.
 */
/*----------------------------------------------------------------------------*/
class GMDSCad_API GeomPoint : public GeomEntity
{
 public:
	/*------------------------------------------------------------------------*/
	/** \brief  Constructor
	 */
	explicit GeomPoint(const std::string &AName = "Unknown point") : GeomEntity(AName) {}

	/*------------------------------------------------------------------------*/
	/** \brief  provides the dimension of the geometrical entity.
	 */
	int dim() const override
	{
		return 0;
	}

	/*------------------------------------------------------------------------*/
	/** \brief  Access to X coordinate
	 *
	 *  \return value of the X coordinate
	 */
	virtual TCoord X() const = 0;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to Y coordinate
	 *
	 *  \return value of the Y coordinate
	 */
	virtual TCoord Y() const = 0;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to Z coordinate
	 *
	 *  \return value of the Z coordinate
	 */
	virtual TCoord Z() const = 0;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to X, Y and Z coordinates
	 *
	 *  \param  ACoordinates will receive the value of the X, Y and Z coordinates
	 */
	virtual void XYZ(TCoord ACoordinates[3]) const
	{
		ACoordinates[0] = X();
		ACoordinates[1] = Y();
		ACoordinates[2] = Z();
	};

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the point as a NumericPoint
	 *
	 *  \return a numeric point
	 */
	math::Point point() const
	{
		TCoord coordinates[3];
		XYZ(coordinates);
		return math::Point(coordinates[0], coordinates[1], coordinates[2]);
	};

	/*------------------------------------------------------------------------*/
	/** \brief Project the point AP unto the geometric entity.
	 *
	 *  \param AP the point to project
	 */
	void project(gmds::math::Point &AP) const override
	{
		AP.setXYZ(X(), Y(), Z());
	};

	/*------------------------------------------------------------------------*/
	/** \brief Get the closest point from AP on the surface
	 *  \param AP a 3D point
	 *
	 *  \return the closest point of APoint on the surface
	 */
	math::Point closestPoint(const math::Point &AP) const override{
		return math::Point(X(),Y(),Z());
	}

	int id() const override = 0;

	/**@brief Accessor to the adjacent curves. Warning, there is no
	 *  assumption about the ordering
	 * @return curves that are adjacent to this point
	 */
	virtual std::vector<GeomCurve *> &curves() = 0;
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
}     // end namespace cad
/*----------------------------------------------------------------------------*/
}     // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GEOM_GEOMPOINT_H_ */
