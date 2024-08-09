//
// Created by chenyt on 19/04/24.
//

#include "gmds/medialaxis/MedialAxisMath.h"
using namespace gmds;

/*----------------------------------------------------------------------------*/
// Computes the maximum distance to a sphere 1 of a point of a sphere 2
double maxRange(const math::Point& ACenter1, const double& ARadius1, const math::Point& ACenter2, const double& ARadius2)
{
	math::Point A = ACenter2 + (-1.)*ACenter1;
	double max_range = vec(A).norm() + ARadius2;
	if (max_range < ARadius1)
		return 0.;
	else
		return max_range - ARadius1;
}
