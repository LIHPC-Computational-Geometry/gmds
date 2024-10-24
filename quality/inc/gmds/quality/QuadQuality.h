/*----------------------------------------------------------------------------*/
#ifndef GMDS_QUAD_QUALITY_H_
#define GMDS_QUAD_QUALITY_H_
/*----------------------------------------------------------------------------*/
#include <gmds/math/Matrix.h>
#include <gmds/math/Vector.h>
#include "GMDSQuality_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace quality {
/*----------------------------------------------------------------------------*/
/** \struct QuadQuality
 *  \brief  Gathers different methods to compute quad quality
 */
struct GMDSQuality_API QuadQuality
{
	QuadQuality() = default;

	math::Vector3d p[4];

	static QuadQuality build(const math::Point& AP0, const math::Point& AP1,
									 const math::Point& AP2, const math::Point& AP3);


	double signedArea() const;
	double aspectRatio() const;
	double condition() const;
	double jacobian() const;
	double scaledJacobian() const;
	double angleDeviation() const;
	double minAngle() const;
	double maxAngle() const;
	double skew() const;

	math::Vector3d L0() const { return p[1]- p[0];}
	math::Vector3d L1() const { return p[2]- p[1];}
	math::Vector3d L2() const { return p[3]- p[2];}
	math::Vector3d L3() const { return p[0]- p[3];}

	math::Vector3d L01() const { return p[1]- p[0];}
	math::Vector3d L12() const { return p[2]- p[1];}
	math::Vector3d L23() const { return p[3]- p[2];}
	math::Vector3d L30() const { return p[0]- p[3];}

	/** Quad diagonals */
	math::Vector3d D0() const {return p[2]-p[0];}
	math::Vector3d D1() const {return p[2]-p[0];}

	/** Principal axes */
	math::Vector3d X1() const {return (p[1]-p[0])+(p[2]-p[3]);}
	math::Vector3d X2() const {return (p[2]-p[1])+(p[3]-p[0]);}
	/** Cross derivatives*/
	math::Vector3d X12() const {return (p[0]-p[1])+(p[2]-p[3]);}
	math::Vector3d X21() const {return (p[0]-p[3])+(p[2]-p[1]);}

	/** corner normals */
	math::Vector3d N0() const {return L3().cross(L0());}
	math::Vector3d N1() const {return L0().cross(L1());}
	math::Vector3d N2() const {return L1().cross(L2());}
	math::Vector3d N3() const {return L2().cross(L3());}
	/** normalized corner normals */
	math::Vector3d n0() const {return N0()/N0().norm();}
	math::Vector3d n1() const {return N1()/N1().norm();}
	math::Vector3d n2() const {return N2()/N2().norm();}
	math::Vector3d n3() const {return N3()/N3().norm();}

	double a0() const {return n0().dot(N0());}
	double a1() const {return n1().dot(N1());}
	double a2() const {return n2().dot(N2());}
	double a3() const {return n3().dot(N3());}

	/** center normal*/
	math::Vector3d NC() const {return X1().cross(X2());}
	/** normalized center normal*/
	math::Vector3d nC() const {return NC()/NC().norm();}


	double l0() const {return L0().norm();}
	double l1() const {return L1().norm();}
	double l2() const {return L2().norm();}
	double l3() const {return L3().norm();}
	double lmin() const {
		return std::min(L0().norm(), std::min( L1().norm(), std::min( L2().norm(),L3().norm() )));
	}

	double lmax() const {
		return std::max(L0().norm(), std::max( L1().norm(), std::max( L2().norm(), L3().norm() )));
	}

	double dmax() const {return std::max(D0().norm(),D1().norm());}
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_QUAD_QUALITY_H_
/*----------------------------------------------------------------------------*/
