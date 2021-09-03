/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * The GMDS library is a computer program whose purpose is to provide a set of
 * functionnalities to represent and handle any type of meshes (2D, 3D,
 * triangles, tetrahedra, quad, hexa, polygons, polyhedra, etc.) and write
 * meshing algorithms. So it gathers many mathematical objects like points,
 * segment, quaternions, etc. and basic algorithms useful to build more evolved
 * ones.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
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
	return Point();
}
/*----------------------------------------------------------------------------*/
Vector3d BezierTriangle::normal(const double AU, const double AV) const {
	return Vector3d();
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
