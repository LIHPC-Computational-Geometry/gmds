/*----------------------------------------------------------------------------*/
#ifndef GMDS_HEX_QUALITY_H_
#define GMDS_HEX_QUALITY_H_
/*----------------------------------------------------------------------------*/
#include <gmds/math/Matrix.h>
#include <gmds/math/Vector.h>
#include "GMDSQuality_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace quality {
/*----------------------------------------------------------------------------*/
/** \struct HexQuality
 *  \brief  Gathers different methods to compute quad quality
 */
struct GMDSQuality_API HexQuality
{
	/** The hexahedral cells are defined with the following order, aka the
	 * "irrational" order
 	 *       P7------P6
 	 *      / |     / |
 	 *    P4------P5  |
 	 *    |   |    |  |
 	 *    |  P3----|--P2
 	 *    | /      | /
 	 *    P0------P1
 	 *
	 */
	HexQuality() = default;

	math::Vector3d p[8];

	static HexQuality build(const math::Point& AP0, const math::Point& AP1,
									const math::Point& AP2, const math::Point& AP3,
									const math::Point& AP4, const math::Point& AP5,
									const math::Point& AP6, const math::Point& AP7);


	double volume() const;
	double diagonal() const;
	double edgeRatio() const;
	double maximumEdgeRatio() const;
	double jacobian() const;
	double scaledJacobian() const;
	double shape() const;
	double maximumAspectFrobenius() const;
	double meanAspectFrobenius() const;
	double oddy() const;
	double shear() const;
	double skew() const;
	double stretch() const;
	double taper() const;

	math::Vector3d L0() const { return p[1]- p[0];}
	math::Vector3d L1() const { return p[2]- p[1];}
	math::Vector3d L2() const { return p[3]- p[2];}
	math::Vector3d L3() const { return p[3]- p[0];}
	math::Vector3d L4() const { return p[4]- p[0];}
	math::Vector3d L5() const { return p[5]- p[1];}
	math::Vector3d L6() const { return p[6]- p[2];}
	math::Vector3d L7() const { return p[7]- p[3];}
	math::Vector3d L8() const { return p[5]- p[4];}
	math::Vector3d L9() const { return p[6]- p[5];}
	math::Vector3d L10() const { return p[7]- p[6];}
	math::Vector3d L11() const { return p[7]- p[4];}

	double l0() const { return L0().norm();}
	double l1() const { return L1().norm();}
	double l2() const { return L2().norm();}
	double l3() const { return L3().norm();}
	double l4() const { return L4().norm();}
	double l5() const { return L5().norm();}
	double l6() const { return L6().norm();}
	double l7() const { return L7().norm();}
	double l8() const { return L8().norm();}
	double l9() const { return L9().norm();}
	double l10() const { return L10().norm();}
	double l11() const { return L11().norm();}

	double lmin() const;
	double lmax() const;



	/** Hex diagonals */
	math::Vector3d D0() const {return p[6]-p[0];}
	math::Vector3d D1() const {return p[7]-p[1];}
	math::Vector3d D2() const {return p[4]-p[2];}
	math::Vector3d D3() const {return p[5]-p[3];}
	double dmin() const {return std::min(D0().norm(),std::min(D1().norm(),std::min(D2().norm(),D3().norm())));}
	double dmax() const {return std::max(D0().norm(),std::max(D1().norm(),std::max(D2().norm(),D3().norm())));}
	//HERE HERE HERE

	/** Principal axes */
	math::Vector3d X1() const {return (p[1]-p[0])+(p[2]-p[3])+(p[5]-p[4])+(p[6]-p[7]);}
	math::Vector3d X2() const {return (p[3]-p[0])+(p[2]-p[1])+(p[7]-p[4])+(p[6]-p[5]);}
	math::Vector3d X3() const {return (p[4]-p[0])+(p[5]-p[1])+(p[6]-p[2])+(p[7]-p[3]);}

	/** Cross derivatives*/
	math::Vector3d X12() const {return (p[2]-p[3])-(p[1]-p[0])+(p[6]-p[7])-(p[5]-p[4]);}
	math::Vector3d X21() const {return X12();}
	math::Vector3d X13() const {return (p[5]-p[1])-(p[4]-p[0])+(p[6]-p[2])-(p[7]-p[3]);}
	math::Vector3d X31() const {return X13();}
	math::Vector3d X23() const {return (p[7]-p[4])-(p[3]-p[0])+(p[6]-p[5])-(p[2]-p[1]);}
	math::Vector3d X32() const {return X23();}

	/**  jacobian matrices */
	math::Matrix33 A0() const {return {L0(),L3(),L4()};}
	math::Matrix33 A1() const {return {L1(),-L0(),L5()};}
	math::Matrix33 A2() const {return {L2(),-L1(),L6()};}
	math::Matrix33 A3() const {return {-L3(),-L2(),L7()};}
	math::Matrix33 A4() const {return {L11(),L8(),-L4()};}
	math::Matrix33 A5() const {return {-L8(),L9(),-L5()};}
	math::Matrix33 A6() const {return {-L9(),L10(),-L6()};}
	math::Matrix33 A7() const {return {-L10(),-L11(),-L7()};}
	math::Matrix33 A8() const {return {X1(),X2(),X3()};}

	/**  normalized jacobian matrices */
	math::Matrix33 A0unit() const {return {L0().normalize(),
														L3().normalize(),
														L4().normalize()};}
	math::Matrix33 A1unit() const {return {L1().normalize(),
														-L0().normalize(),
														L5().normalize()};}
	math::Matrix33 A2unit() const {return {L2().normalize(),
														-L1().normalize(),
														L6().normalize()};}
	math::Matrix33 A3unit() const {return {-L3().normalize(),
														-L2().normalize(),
														L7().normalize()};}
	math::Matrix33 A4unit() const {return {L11().normalize(),
														L8().normalize(),
														-L4().normalize()};}
	math::Matrix33 A5unit() const {return {-L8().normalize(),
														L9().normalize(),
														-L5().normalize()};}
	math::Matrix33 A6unit() const {return {-L9().normalize(),
														L10().normalize(),
														-L6().normalize()};}
	math::Matrix33 A7unit() const {return {-L10().normalize(),
														-L11().normalize(),
														-L7().normalize()};}
	math::Matrix33 A8unit() const {return {X1().normalize(),
														X2().normalize(),
														X3().normalize()};}
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_QUAD_QUALITY_H_
/*----------------------------------------------------------------------------*/
