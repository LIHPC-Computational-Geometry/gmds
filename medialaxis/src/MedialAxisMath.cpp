//
// Created by chenyt on 19/04/24.
//

#include "gmds/medialaxis/MedialAxisMath.h"
using namespace gmds;

// Computes the circumcenter of a triangle
math::Point circumcenter(const math::Triangle& ATri){
	math::Point A = ATri.getPoint(0);
	math::Point B = ATri.getPoint(1);
	math::Point C = ATri.getPoint(2);
	double xA = A.X();
	double yA = A.Y();
	double xB = B.X();
	double yB = B.Y();
	double xC = C.X();
	double yC = C.Y();
	double S = ATri.area();
	double xO = ((xA*xA+yA*yA)*(yB-yC)-(xB*xB+yB*yB)*(yA-yC)+(xC*xC+yC*yC)*(yA-yB))*(1./(4.*S));
	double yO = -((xA*xA+yA*yA)*(xB-xC)-(xB*xB+yB*yB)*(xA-xC)+(xC*xC+yC*yC)*(xA-xB))*(1./(4.*S));
	math::Point Ctr(xO, yO, 0.);
	return Ctr;
}

/*----------------------------------------------------------------------------*/
// Computes the circumcenter of a tetrahedron
math::Point circumcenterTetra(const math::Tetrahedron& ATetra){
	// Get the tetra vertices
	math::Point A = ATetra.getPoint(0);
	math::Point B = ATetra.getPoint(1);
	math::Point C = ATetra.getPoint(2);
	math::Point D = ATetra.getPoint(3);
	// Get the coordinates
	double xA = A.X();
	double yA = A.Y();
	double zA = A.Z();
	double xB = B.X();
	double yB = B.Y();
	double zB = B.Z();
	double xC = C.X();
	double yC = C.Y();
	double zC = C.Z();
	double xD = D.X();
	double yD = D.Y();
	double zD = D.Z();
	math::Matrix33 M;
	M(0,0) = xB-xA;
	M(0,1) = yB-yA;
	M(0,2) = zB-zA;
	M(1,0) = xC-xB;
	M(1,1) = yC-yB;
	M(1,2) = zC-zB;
	M(2,0) = xD-xA;
	M(2,1) = yD-yA
	   ;
	M(2,2) = zD-zA;
	double xV = (xB-xA)*(1./2.)*(xA+xB)+(yB-yA)*(1./2.)*(yA+yB)+(zB-zA)*(1./2.)*(zA+zB);
	double yV = (xC-xB)*(1./2.)*(xB+xC)+(yC-yB)*(1./2.)*(yB+yC)+(zC-zB)*(1./2.)*(zB+zC);
	double zV = (xD-xA)*(1./2.)*(xA+xD)+(yD-yA)*(1./2.)*(yA+yD)+(zD-zA)*(1./2.)*(zA+zD);
	math::Vector3d V;
	V.set(0, xV);
	V.set(1, yV);
	V.set(2, zV);
	math::Vector3d X = M.solve(V);
	math::Point Ctr(X.X(), X.Y(), X.Z());
	return Ctr;
}

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
