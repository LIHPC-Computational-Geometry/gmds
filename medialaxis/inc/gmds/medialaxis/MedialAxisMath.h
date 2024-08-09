//
// Created by chenyt on 19/04/24.
//

#ifndef GMDS_MEDIALAXISMATH_H
#define GMDS_MEDIALAXISMATH_H

#include <gmds/math/Triangle.h>
#include <gmds/math/Point.h>
#include <gmds/math/Tetrahedron.h>
#include <gmds/math/Matrix.h>
#include <gmds/math/Vector.h>
using namespace gmds;

// Computes the maximum distance to a sphere 1 of a point of a sphere 2 which doesn't belong to the sphere 1
double maxRange(const math::Point& ACenter1, const double& ARadius1, const math::Point& ACenter2, const double& ARadius2);

#endif     // GMDS_MEDIALAXISMATH_H
