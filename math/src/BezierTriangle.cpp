/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
 * BezierCurve.cpp
 *
 */
/*----------------------------------------------------------------------------*/
#include <gmds/math/BezierTriangle.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::math;
/*----------------------------------------------------------------------------*/
// b macros give us the flat representation of the (i,j,k) indices in control
// points


const int b_index[4][4][4] = {
        {
                {-1, -1, -1,  0}, //0,0,k
                {-1, -1,  1, -1}, //0,1,k
                {-1,  2, -1, -1}, //0,2,k
                { 3, -1, -1, -1} //0,3,k
        },
                              {{},{},{},{}},
                              {{},{},{},{}},
                              {{},{},{},{}}};
/*----------------------------------------------------------------------------*/
BezierTriangle::BezierTriangle(const Point &AP1, const Point &AP2, const Point &AP3, const Vector3d &AN1,
                               const Vector3d &AN2, const Vector3d &AN3) {}
/*----------------------------------------------------------------------------*/
Point BezierTriangle::operator()(const double AU, const double AV) const {
}
/*----------------------------------------------------------------------------*/
Vector3d BezierTriangle::normal(const double AU, const double AV) const {
}
/*----------------------------------------------------------------------------*/
void BezierTriangle::geomInfo(const double AU, const double AV, Point &AP, Vector3d &AN, Vector3d &ADU,
                              Vector3d &ADV) const
{}
/*----------------------------------------------------------------------------*/
std::vector<Point> BezierTriangle:: getDiscretization(const int ANb) const {
    if(ANb<1)
        throw GMDSException("BezierCurve discretization impossible with this parameter");


    std::vector<Point> points;

    return points;
}
/*----------------------------------------------------------------------------*/
