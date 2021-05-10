/*-----------------------------------------------------------------*/
/*
 * Cross.h
 *
 *  Created on: 03/01/2015
 */
/*-----------------------------------------------------------------*/
#ifndef GMDS_MATH_CROSS_H_
#define GMDS_MATH_CROSS_H_
/*-----------------------------------------------------------------*/
#include <gmds/math/Vector.h>
#include <gmds/math/Quaternion.h>
#include <cmath>
/*-----------------------------------------------------------------*/
#include <iostream>
/*-----------------------------------------------------------------*/
using namespace std;
namespace gmds {
  namespace math {
    /*-----------------------------------------------------------------*/
    class EXPORT_GMDS Cross
    {
    private:
      /// A cross is represented by two orthogonal vectors. 
      /* Using only one representation vector is only possible if we
       * can define a reference vector int the whole world (not done
       * in 3D, try to use transport equation along surface?)
       */
      math::Vector3d m_x;
      math::Vector3d m_y;

    public:

      /*------------------------------------------------------------*/
      /* Constructor
       *  Defines a cross with (1,0,0) and (0,1,0)
       */
      Cross();
      /*------------------------------------------------------------*/
      /* Constructor
       *  Defines a cross from two orthogonal vectors.
       */
      Cross(math::Vector3d& AV1, math::Vector3d& AV2);
      /*------------------------------------------------------------*/
      /* Constructor
       * Defines a cross from a quaternion and a normal vector to be
       * aligned with
       */
      Cross(const Quaternion& AQ, const math::Vector3d& AN);

      /*------------------------------------------------------------*/
      /* returns the right-hand chart (X,Y, XxY)
       */
      math::Chart chart() const;

      /*------------------------------------------------------------*/
      /*
       *  return the vector, which is the closest of AN between m_x,
       * m_y, -m_x, -m_y.
       */
      math::Vector3d closestVector(const math::Vector3d& AN);

      /*------------------------------------------------------------*/
      /*
       *  return the vector x
       */
      math::Vector3d X() const { return m_x; }
      /*------------------------------------------------------------*/
      /*
       *  return the vector y
       */
      math::Vector3d Y() const { return m_y; }
    };


    EXPORT_GMDS ostream & operator << (ostream & op_g, const Cross & op_d);
  }
  /*-----------------------------------------------------------------*/
}
/*-----------------------------------------------------------------*/
#endif /* GMDS_MATH_CROSS_H_ */
/*-----------------------------------------------------------------*/
