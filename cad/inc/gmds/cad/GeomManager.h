/*----------------------------------------------------------------------------*/
/** \file    GeomManager.h
 *  \author  F. LEDOUX
 *  \date    30/06/11
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_GEOMMANAGER_H_
#define GMDS_GEOM_GEOMMANAGER_H_
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "gmds/utils/Exception.h"
#include "gmds/utils/CommonTypes.h"
#include "gmds/cad/GeomPoint.h"
#include "gmds/cad/GeomCurve.h"
#include "gmds/cad/GeomSurface.h"
#include "gmds/cad/GeomVolume.h"
#include "GMDSCad_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace cad{
/*----------------------------------------------------------------------------*/
/** \class GeomManager
 *  \brief This interface gathers the factory methods required for every
 *  		geometric model, the access to all the geom entities stored in the
 *  		model as the responsability to delete geometric entities.
 *
 *  \param TBase the basic type used to store geometric coordinates.
 */
/*----------------------------------------------------------------------------*/
	class GMDSCad_API GeomManager {

public:

	/*------------------------------------------------------------------------*/
	/** \brief  creation of a geometric volume
	 */
	virtual GeomVolume* newVolume() =0;

	/*------------------------------------------------------------------------*/
	/** \brief  creation of a geometric surface
	 */
	virtual GeomSurface* newSurface() =0;

	/*------------------------------------------------------------------------*/
	/** \brief  creation of a geometric curve
	 */
	virtual GeomCurve* newCurve() =0;

	/*------------------------------------------------------------------------*/
	/** \brief  creation of a geometric point
	 */
	virtual GeomPoint* newPoint() =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of points of the model.
	 *
	 *	\return the number of points.
	 */
	virtual TInt getNbPoints() const =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of curves of the model.
	 *
	 *	\return the number of curves.
	 */
	virtual TInt getNbCurves() const =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of surfaces of the model.
	 *
	 *	\return the number of surfaces.
	 */
	virtual TInt getNbSurfaces() const =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of volumes of the model.
	 *
	 *	\return the number of volumes.
	 */
	virtual TInt getNbVolumes() const =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the points of the model.
	 *
	 *  \param points the points of the model.
	 */
	virtual void getPoints(std::vector<GeomPoint*>& points) const =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the curves of the model.
	 *
	 *  \param curves the curves of the model.
	 */
	virtual void getCurves(std::vector<GeomCurve*>& curves) const =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the surface of the model.
	 *
	 *  \param surfaces the surfaces of the model.
	 */
	virtual void getSurfaces(std::vector<GeomSurface*>& surfaces)const=0;

		/*------------------------------------------------------------------------*/
		/** \brief  Access to the volumes of the model.
         *  \param volumes the volumes of the model.
         */
		virtual void getVolumes(std::vector<GeomVolume*>& volumes) const =0;

        /*------------------------------------------------------------------------*/
        /** \brief  Gives access to the point of id @AID, Return NullPtr if it does
         *          not exist.
         *  \return A point
         */
        virtual GeomPoint* getPoint(const TInt AID)=0;
        /*------------------------------------------------------------------------*/
        /** \brief  Gives access to the curve of id @AID, Return NullPtr if it does
         *          not exist.
         *  \return A curve
         */
        virtual GeomCurve* getCurve(const TInt AID)=0;
        /*------------------------------------------------------------------------*/
        /** \brief  Gives access to the surface of id @AID, Return NullPtr if it
         *          does not exist.
         *  \return A surface
         */
        virtual GeomSurface* getSurface(const TInt AID)=0;
        /*------------------------------------------------------------------------*/
        /** \brief  Gives access to the volume of id @AID, Return NullPtr if it
         *          does not exist.
         *  \return A volume
         */
        virtual GeomVolume* getVolume(const TInt AID)=0;
		/*------------------------------------------------------------------------*/
		/** \brief  Get the curve common to 2 points
         *
         *  \param return the id of the common curve, and -1 if it doesn't exist
         */
		virtual int getCommonCurve(GeomPoint* AP1, GeomPoint* AP2) const =0;

		/*------------------------------------------------------------------------*/
		/** \brief  Get the surface common to 2 curves
         *
         *  \param return the id of the common surface, and -1 if it doesn't exist
         */
		virtual int getCommonSurface(GeomCurve* AC1, GeomCurve* AC2) const =0;
};
/*----------------------------------------------------------------------------*/
} // namespace cad
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GEOM_GEOMMANAGER_H_ */

