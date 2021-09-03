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
 * QualityMeasure.cpp
 *
 *  Created on: 12 april 2017
 *      Author: ledoux f.
 */
/*----------------------------------------------------------------------------*/
#include <gmds/math/QualityMeasure.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace math {
/*----------------------------------------------------------------------------*/
    double
    QualityMeasure::sphereRatio(const Point& AP0, const Point& AP1, const Point& AP2, const Point& AP3)
    {
            //==========================================
            // Useful vectors for our computations
            //==========================================
            math::Vector3d v01(AP0, AP1);
            math::Vector3d v02(AP0, AP2);
            math::Vector3d v03(AP0, AP3);
            math::Vector3d v10(AP1, AP0);
            math::Vector3d v13(AP1, AP3);
            math::Vector3d v21(AP2, AP1);
            math::Vector3d v20(AP2, AP0);
            math::Vector3d v23(AP2, AP3);

            //==========================================
            // Volume computation
            //==========================================
            double v = v01.dot(v02.cross(v03)) / 6.0;
            //==========================================
            // area of each face
            //==========================================
            double a012 = 0.5 * (v01.cross(v02)).norm();
            double a103 = 0.5 * (v10.cross(v13)).norm();
            double a210 = 0.5 * (v21.cross(v20)).norm();
            double a023 = 0.5 * (v02.cross(v03)).norm();
            //==========================================
            // Opposite edge norm products
            //==========================================
            double p1 = v01.norm() * v23.norm();
            double p2 = v02.norm() * v13.norm();
            double p3 = v03.norm() * v21.norm();

            double num = 216 * v * v;

            double den_part1 = a012 + a103 + a210 + a023;
            double den_part2 = (p1 + p2 + p3) * (p1 + p2 - p3) * (p1 + p3 - p2) * (p2 + p3 - p1);
            double den = den_part1 * sqrt(den_part2);

            double rho = num / den;

            if (rho > 1.0)
                    rho = 1.0;

            return rho;
    }
/*----------------------------------------------------------------------------*/
    void
    QualityMeasure::extremAngles(const Triangle &AT,
                                 double &ASmallAngle, double &ALargeAngle) {
        extremAngles(AT.getPoint(0),AT.getPoint(1),AT.getPoint(2),
                     ASmallAngle, ALargeAngle);
    }
/*----------------------------------------------------------------------------*/
    void
    QualityMeasure::extremAngles(const Point &AP0, const Point &AP1, const Point &AP2,
                                 double &ASmallAngle,
                                 double &ALargeAngle) {
        double small=0, large=0;

        small = Vector3d(AP0, AP2).angle(Vector3d(AP0, AP1));
        large = small;

        double aij = Vector3d(AP1, AP0).angle(Vector3d(AP1, AP2));
        if(aij<small)
            small=aij;
        if(aij>large)
            large=aij;

        aij = Vector3d(AP2, AP1).angle(Vector3d(AP2, AP0));
        if(aij<small)
            small=aij;
        if(aij>large)
            large=aij;

        ASmallAngle = small / math::Constants::PIDIV2;
        ALargeAngle = large / math::Constants::PIDIV2;

    }
/*----------------------------------------------------------------------------*/
}  // namespace math
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/
