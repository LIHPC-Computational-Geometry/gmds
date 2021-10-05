/*----------------------------------------------------------------------------*/
/** \file    GeomEntity.h
 *  \author  F. LEDOUX
 *  \date    09/21/2010
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_GEOMENTITY_H_
#define GMDS_GEOM_GEOMENTITY_H_
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
#include <gmds/math/Point.h>
#include "GMDSCad_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
	namespace cad{
/*----------------------------------------------------------------------------*/
/** \class Entity
 *  \brief This class provides services thar are common to all the geometric
 *  	   entities (volume, surface, curve, point).
 */
/*----------------------------------------------------------------------------*/
		class GMDSCad_API GeomEntity {
		public:

			/*------------------------------------------------------------------------*/
			/** \brief  Default constructor
             */
			GeomEntity(const std::string& AName = "Unknown entity"):
					name_(AName){;}

			/*------------------------------------------------------------------------*/
            /** \brief  sets the name of the geometrical entity.
             */
            void setName(std::string AName) {name_ = AName;}

			/*------------------------------------------------------------------------*/
			/** \brief  provides the name f the geometrical entity.
             */
			std::string getName() const {return name_;}

			/*------------------------------------------------------------------------*/
			/** \brief  provides the dimension of the geometrical entity.
             */
			virtual int getDim() const=0;

			/*------------------------------------------------------------------------*/
			/** \brief  computes the area of the entity.
             */
			virtual TCoord computeArea() const= 0;

			/*------------------------------------------------------------------------*/
			/** \brief  computes the bounding box
             *
             *	\param minXYZ The minimum coordinate of the bounding box.
             *	\param maxXYZ The maximum coordinate of the bounding box.
             */
			virtual void computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const=0;

			/*------------------------------------------------------------------------*/
			/** \brief Project the point AP unto the geometric entity.
          *
          *  \param AP the point to project
             */
			virtual void project(gmds::math::Point& AP) const =0;

			virtual int id() const=0;

		protected:

			/* name of the entity if it has one */
			std::string name_;
		};
/*----------------------------------------------------------------------------*/
	} // namespace cad
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GEOM_GEOMENTITY_H_ */

