/*----------------------------------------------------------------------------*/
/** \file    GeomSurface.h
 *  \author  F. LEDOUX
 *  \date    09/21/2010
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_GEOMSURFACE_H_
#define GMDS_GEOM_GEOMSURFACE_H_
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "gmds/utils/Exception.h"
#include "gmds/utils/CommonTypes.h"
#include "gmds/cad/GeomEntity.h"
#include "GMDSCad_export.h"
#include "gmds/math/Point.h"
#include "gmds/math/Vector.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
    namespace cad{
        class GeomPoint;
        class GeomCurve;
        class GeomVolume;
/*----------------------------------------------------------------------------*/
/** \class Surface
 *  \brief This class describe the services that are required by the
 *  	   mesh to the geometrical model. As a consequence, this interface only
 *  	   contains query methods.
 */
/*----------------------------------------------------------------------------*/
        class GMDSCad_API GeomSurface : public GeomEntity {

        public:


            /*------------------------------------------------------------------------*/
            /** \brief  Constructor
             */
            explicit GeomSurface(const std::string& AName = "Unknown surface")
                    :GeomEntity(AName){}

            /*------------------------------------------------------------------------*/
            /** \brief  provides the dimension of the geometrical entity.
             */
            int dim() const override {return 2;}

            /** \brief Move a point AP near the surface to the closest point on the
             * 		   surface.
             *  \param AP
             */
            void project(math::Point& AP) const override;

            /*------------------------------------------------------------------------*/
            /** \brief  computes normal at the closest point to AP in 3D.
             *
             *  \param AP the point
             *  \param AV 	normal vector at the closest point of AP on this
             */
            virtual void computeNormal(	const math::Point& AP,
                                           math::Vector3d& AV) const =0;

            /*------------------------------------------------------------------------*/
            /** \brief  computes the bounding box
             *
             *  \param minXYZ The minimum coordinate of the bounding box.
             *  \param maxXYZ The maximum coordinate of the bounding box.
             */
            void computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const override =0;

            /**@brief Accessor to the adjacent points. Warning, there is no
             *  assumption about the ordering
             * @return points that are adjacent to this point
             */
            virtual std::vector<GeomPoint*>& points()=0;

            /**@brief Accessor to the adjacent curves. Warning, there is no
             *  assumption about the ordering
             * @return curves that are adjacent to this point
             */
            virtual std::vector<GeomCurve*>& curves()=0;

            /**@brief Accessor to the adjacent volumes. Warning, there is no
             *  assumption about the ordering
             * @return volumes that are adjacent to this point
             */
            virtual std::vector<GeomVolume*>& volumes()=0;


        };
/*----------------------------------------------------------------------------*/
    } // end namespace cad
/*----------------------------------------------------------------------------*/
} // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GEOM_GEOMSURFACE_H_ */

