/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_DISCRETIZATION_SCHEME_1D_H_
#define GMDS_MATH_DISCRETIZATION_SCHEME_1D_H_
/*----------------------------------------------------------------------------*/
// Gepeto File Headers
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
#include "GMDSMath_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
    namespace math {
/*----------------------------------------------------------------------------*/
/** \class Discretization1D
 *  \brief Abstract class that must be implemented by 1D discretization scheme.
 *         The only constraint is that it must have 2 end points
 *
 */
/*----------------------------------------------------------------------------*/
        class GMDSMath_API DiscretizationScheme1D {
        public:

            /** @brief  constructor.
             * @param AOrigin   origin point
             * @param ADest     destination point
             * @param ANbPoints total number of points including the origin and
             *                  destination point
             */
            DiscretizationScheme1D(const Point &AOrigin = math::Point(0, 0, 0),
                                   const Point &ADest = math::Point(1, 0, 0),
                                   const TInt ANbPoints = 10);

            void setOrigin(const math::Point &AP) { m_origin = AP; }

            void setDestination(const math::Point &AP) { m_destination = AP; }

            void setNbPoints(const int ANb) {
                if (ANb > 0)
                    m_nb_points = ANb;
            }

            TInt getNbPoints() const { return m_nb_points; }

            /** @brief Access to the @p AIndex th point of the discretization
             *         knowing that m_origin is at index 0 and ADest is the
             *         last point, i.e. at index m_nb_points-1;
             *
             *  @throw an exception if @p AIndex < 0 or >= m_nb_points
             *
             * @param AIndex the point we want to get the valule
             * @return the point location
             */
            virtual Point operator()(const int AIndex) const = 0;

        protected:
            math::Point m_origin;
            math::Point m_destination;
            TInt m_nb_points;
        };
/*----------------------------------------------------------------------------*/
/** \class DiscretizationScheme1DUniform
 *  \brief Linear discretization between 2 points with an uniform step (same
 *         distance between consecutive points)
 */
/*----------------------------------------------------------------------------*/
        class GMDSMath_API DiscretizationScheme1DUniform : public DiscretizationScheme1D {
        public:
            /** @brief  constructor.
               * @param AOrigin   origin point
               * @param ADest     destination point
               * @param ANbPoints total number of points including the origin and
               *                  destination point
               */
            DiscretizationScheme1DUniform(const Point &AOrigin,
                                          const Point &ADest,
                                          const TInt ANbPoints);

            /** @brief Access to the @p AIndex th point of the discretization
             *         knowing that m_origin is at index 0 and ADest is the
             *         last point, i.e. at index m_nb_points-1;
             *
             *  @throw an exception if @p AIndex < 0 or >= m_nb_points
             *
             * @param AIndex the point we want to get the valule
             * @return
             */
            virtual Point operator()(const int AIndex) const;
        };
/*----------------------------------------------------------------------------*/
/** \class DiscretizationScheme1DGeometric
 *  \brief 1D discretization between 2 points with a geometric step. A reason
 *          in ]0,1[ is given and the size of element is computed along it
 *
 */
/*----------------------------------------------------------------------------*/
        class GMDSMath_API DiscretizationScheme1DGeometric : public DiscretizationScheme1D {
        public:
            /** @brief  constructor.
               * @param AReason a number in ]0,1
               * @param AOrigin   origin point
               * @param ADest     destination point
               * @param ANbPoints total number of points including the origin and
               *                  destination point
               */
            DiscretizationScheme1DGeometric
                    (const double AReason,
                     const Point &AOrigin,
                     const Point &ADest,
                     const TInt ANbPoints);

            void setInverse(const bool& AInv);
            /** @brief Access to the @p AIndex th point of the discretization
             *         knowing that m_origin is at index 0 and ADest is the
             *         last point, i.e. at index m_nb_points-1;
             *
             *  @throw an exception if @p AIndex < 0 or >= m_nb_points
             *
             * @param AIndex the point we want to get the valule
             * @return
             */
            virtual Point operator()(const int AIndex) const;

        private:
            TCoord m_reason;
            bool m_inverse;
        };
/*----------------------------------------------------------------------------*/
    }
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_DISCRETIZATION_SCHEME_1D_H_ */
/*----------------------------------------------------------------------------*/
