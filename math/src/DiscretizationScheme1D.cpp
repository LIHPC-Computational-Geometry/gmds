/*-----------------------------------------------------------------*/
#include "gmds/math/DiscretizationScheme1D.h"
/*-----------------------------------------------------------------*/
#include <cmath>
/*-----------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::math;
/*-----------------------------------------------------------------*/
DiscretizationScheme1D::
DiscretizationScheme1D(const Point &AOrigin,
                       const Point &ADest,
                       const TInt ANbPoints)
        : m_origin(AOrigin),m_destination(ADest), m_nb_points(ANbPoints) {}
/*-----------------------------------------------------------------*/
DiscretizationScheme1DUniform::
DiscretizationScheme1DUniform(const Point &AOrigin, const Point &ADest,
                              const TInt ANbPoints)
        : DiscretizationScheme1D(AOrigin,ADest,ANbPoints) {}

/*-----------------------------------------------------------------*/
Point DiscretizationScheme1DUniform::operator()(const int AIndex) const
{
    if(AIndex<0 || AIndex>=m_nb_points)
        throw GMDSMathException("Range access error in Discretization1D");
    return 1.0/(m_nb_points-1)*( (m_nb_points-1-AIndex)*m_origin+(AIndex*m_destination));
}
/*-----------------------------------------------------------------*/
DiscretizationScheme1DGeometric::
DiscretizationScheme1DGeometric(const double AReason, const Point &AOrigin,
                                const Point &ADest, const TInt ANbPoints)
        : DiscretizationScheme1D(AOrigin,ADest,ANbPoints)
{
    if(AReason<=0 || AReason>=1)
        throw GMDSMathException("DiscretizationScheme1DGeometric: reason must be comprised in ]0,1[");
    m_reason = AReason;
    m_inverse=false;

}
/*-----------------------------------------------------------------*/
Point DiscretizationScheme1DGeometric::operator()(const int AIndex) const
{
    if(AIndex<0 || AIndex>m_nb_points)
        throw GMDSMathException("Range access error in Discretization1D");
    math::Vector3d v =m_destination-m_origin;
    double r = 0;
    if(AIndex==0)
        return m_origin;

    if (AIndex==m_nb_points)
        return m_destination;

    if(m_inverse){
        for(auto i = m_nb_points-1;i>=AIndex;i--) {
            r += std::pow(m_reason, m_nb_points-i);
        }
        return m_destination+(r*v).opp();
    }
    else{
        for(auto i = 1;i<=AIndex;i++) {
            r += std::pow(m_reason, i);
        }
        return m_origin+r*v;
    }
}
/*-----------------------------------------------------------------*/
void DiscretizationScheme1DGeometric::setInverse(const bool &AInv) {
    m_inverse=AInv;

}