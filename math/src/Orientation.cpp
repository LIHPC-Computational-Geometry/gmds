/*----------------------------------------------------------------------------*/
// GMDS File Headers
#include "gmds/math/Orientation.h"
/*---------------------------------------------------------------------------*/
// Predicates_psm File Headers
#include <Predicates_psm.h>
/*----------------------------------------------------------------------------*/
using  namespace gmds;
using namespace gmds::math;
/*---------------------------------------------------------------------------*/
void Orientation::initialize() {
    GEO::PCK::initialize();
}
/*---------------------------------------------------------------------------*/
void Orientation::finalize()  {
    GEO::PCK::terminate();
}
/*---------------------------------------------------------------------------*/
Orientation::Sign Orientation::orient3d(const gmds::math::Point& AP0,
                                     const gmds::math::Point& AP1,
                                     const gmds::math::Point& AP2,
                                     const gmds::math::Point& AP3)
{
    double p0[3]={AP0.X(),AP0.Y(),AP0.Z()};
    double p1[3]={AP1.X(),AP1.Y(),AP1.Z()};
    double p2[3]={AP2.X(),AP2.Y(),AP2.Z()};
    double p3[3]={AP3.X(),AP3.Y(),AP3.Z()};
    return (Orientation::Sign)(GEO::PCK::orient_3d(p0, p1, p2, p3));
}