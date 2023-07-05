/*----------------------------------------------------------------------------*/
/*
 * FACVolume.h
 *
 *  Created on: 27 juin 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_FACETEDVOLUME_H_
#define GMDS_GEOM_FACETEDVOLUME_H_
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "gmds/cad/GeomVolume.h"
#include "gmds/cadfac/FACSurface.h"
#include "GMDSCadFac_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
    namespace cad{
/*----------------------------------------------------------------------------*/
/** \class FACVolume
 *  \brief This class implements the volume services that are required by the
 *  	   mesh to the geometrical model.
 *  \param TCoord the data type used to store geometrical data
 */
/*----------------------------------------------------------------------------*/
        class GMDSCadFac_API FACVolume : public GeomVolume {

        public:
            /*------------------------------------------------------------------------*/
            /** \brief  Constructor.
             */
            FACVolume(const std::string& AName = "Unknown volume");


            /*------------------------------------------------------------------------*/
            /** \brief  Destructor.
             */
            virtual ~FACVolume();


            /*------------------------------------------------------------------------*/
            /** \brief  computes the area of the entity.
             */
            virtual TCoord computeArea() const;

            /*------------------------------------------------------------------------*/
            /** \brief  computes the bounding box
             *
             *	\param minXYZ The minimum coordinate of the bounding box.
             *	\param maxXYZ The maximum coordinate of the bounding box.
             */
            virtual void computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const;
				virtual std::tuple<TCoord,TCoord,TCoord,TCoord,TCoord,TCoord>  BBox() const;

            int id() const {return m_id;}


            /**@brief Accessor to the adjacent points. Warning, there is no
             *  assumption about the ordering
             * @return points that are adjacent to this point
             */
            virtual std::vector<GeomPoint*>& points();

            /**@brief Accessor to the adjacent curves. Warning, there is no
             *  assumption about the ordering
             * @return curves that are adjacent to this point
             */
            virtual std::vector<GeomCurve*>& curves();

            /**@brief Accessor to the adjacent surfaces. Warning, there is no
             *  assumption about the ordering
             * @return surfaces that are adjacent to this point
             */
            virtual std::vector<GeomSurface*>& surfaces();

            /**@brief Reset the global id counter to 1.
             */
             static void resetIdCounter();
        private:

            int m_id;
            static int m_next_id;

            /** adjacent geometric points that bound this volume*/
            std::vector<GeomPoint*> m_adjacent_points;
            /** adjacent geometric curves that bound this volume*/
            std::vector<GeomCurve*> m_adjacent_curves;
            /** adjacent geometric surfaces that bound this volume*/
            std::vector<GeomSurface*> m_adjacent_surfaces;
        };
/*----------------------------------------------------------------------------*/
    } // namespace cad
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GEOM_FACETEDVOLUME_H_ */
/*----------------------------------------------------------------------------*/
