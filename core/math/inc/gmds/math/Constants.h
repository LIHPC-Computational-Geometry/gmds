/*----------------------------------------------------------------------------*/
/*
 * Constants.h
 *
 *  Created on: 6 févr. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_CONSTANTS_H_
#	define GMDS_MATH_CONSTANTS_H_
/*----------------------------------------------------------------------------*/
#	include <gmds/utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace math {
/*----------------------------------------------------------------------------*/
namespace Constants {
const TCoord EPSILON = static_cast<TCoord>(1.e-8);
const TCoord PI = static_cast<TCoord>(3.1415926535897931);
const TCoord PI2 = static_cast<TCoord>(6.2831853071795862);
const TCoord PIDIV2 = static_cast<TCoord>(1.5707963267948966);
const TCoord PIDIV3 = static_cast<TCoord>(1.0471975511965976);
const TCoord PIDIV4 = static_cast<TCoord>(0.78539816339744828);
const TCoord PIDIV6 = static_cast<TCoord>(0.52359877559829882);
const TCoord PIDIV8 = static_cast<TCoord>(0.39269908169872414);
const TCoord INVPIDIV180 = static_cast<TCoord>(57.295779513082323);
const TCoord PIDIV180 = static_cast<TCoord>(0.017453292519943295);

}     // namespace Constants
      /*----------------------------------------------------------------------------*/
}     // namespace math
      /*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_CONSTANTS_H_ */
/*----------------------------------------------------------------------------*/
