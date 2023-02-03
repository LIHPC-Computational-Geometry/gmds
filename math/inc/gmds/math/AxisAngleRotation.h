/*----------------------------------------------------------------------------*/
/*
 * AxisAngleRotation
 *
 *  Created on: December 2, 2015
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_AXIS_ANGLE_ROTATION_H_
#define GMDS_MATH_AXIS_ANGLE_ROTATION_H_
/*----------------------------------------------------------------------------*/
// gmds file headers
/*----------------------------------------------------------------------------*/
#include <gmds/math/Matrix.h>
#include <gmds/math/Quaternion.h>
#include <gmds/math/Chart.h>
#include "GMDSMath_export.h"
/*----------------------------------------------------------------------------*/
#include<cmath>
#include<cstring>
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*------------------------------------------------------------------------*/
    namespace math{
        /*--------------------------------------------------------------------*/
        /** \class Vector
         *  \brief Defines a 3D Vector
         */
        class GMDSMath_API AxisAngleRotation {
        public:
            
            /*---------------------------------------------------------------*/
            /** \brief Default Constructor with null vector and null angle
             */
            AxisAngleRotation();
            
            
            /*---------------------------------------------------------------*/
            /** \brief Constructor.
             *
             * \param[in] AV the rotation vector
             * \param[in] AA the rotation angle
             */
            AxisAngleRotation(const Vector3d& AV, double AA);
            /*---------------------------------------------------------------*/
            /** \brief Constructor.
             *
             * \param[in] AV a vector defining the rotation axis and the angle (via
             *           its norm
             */
            explicit AxisAngleRotation(const Vector3d& AV);
            
            /*---------------------------------------------------------------*/
            /** \brief Constructor from a quaternion
             *
             * \param AQ a unit quaternion corresponding to a rotation
             */
            explicit AxisAngleRotation(const Quaternion& AQ);
            
            
            
            /*---------------------------------------------------------------*/
            /** \brief Constructor from a chart \p AC. It gives the rotation
             *         bringing (OX,OY,OZ) onto the (X,Y,Z) axis of \p AC.
             *
             * \param[in] AC a chart
             */
            explicit AxisAngleRotation(const Chart& AC);
            
            
            /*---------------------------------------------------------------*/
            /**  \brief Constructs an axis-angle rotation which transforms 
             *          \p AFrom axis into \p ATo axis
             *
             * \param[in] AFrom a first vector
             * \param[in] ATo   a second vector
             */
            AxisAngleRotation(const Vector3d& AFrom, const Vector3d& ATo);
            
            /*---------------------------------------------------------------*/
            /** \brief Provides a quaternion representation of this rotation.
             */
            Quaternion quaternion() const;
            
            /*---------------------------------------------------------------*/
            /** \brief Returns the rotation vector
             */
            const Vector3d& axis() const {return m_axis;}
            Vector3d& axis() {return m_axis;}
            
            
            /*---------------------------------------------------------------*/
            /** \brief Returns the rotation angle
             */
            double angle() const {return m_axis.norm();}
            
            /*---------------------------------------------------------------*/
            /** \brief Build the inverse rotation, that is the an axis-angle
             *         rotation with opposite rotation angle */
            AxisAngleRotation inverse() const
            { return {m_axis,m_axis.norm()}; }
            
            /*---------------------------------------------------------------*/
            /** \brief Build the identity axis-angle rotation */
            static AxisAngleRotation identity()
            { return {Vector3d({1.,0.,0.}),0.}; }
            
            /*---------------------------------------------------------------*/
            /** \brief Build a corresponding 3x3 rotation matrix */
            Matrix<3,3,double> toRotationMatrix() const;
            Vector3d toRotationAxis(int AIndex) const;
            
            /*---------------------------------------------------------------*/
            /** \brief Build the chart corresponding to applying *this to
             *         the reference Chart(OX, OY,OZ)*/
            Chart toChart() const;
            /*---------------------------------------------------------------*/
            /**  \brief Gets an axis-angle rotation which transforms Z axis into
             *          \p AV axis
             *
             * \param[in]  AV the constrained axis
             * \return an axis-angle rotation
             */
            static AxisAngleRotation alignZ(const Vector3d& AV);
            /*---------------------------------------------------------------*/
            /**  \brief Gets an axis-angle rotation which transforms Y and Z 
             *          axis into \p AV1 and \p AVE respectively
             *
             * \param[in]  AY the axis Y must be aligned with
             * \param[in]  AZ the axis Z must be aligned with
             * \return an axis-angle rotation
             */

            static AxisAngleRotation alignYZ(const Vector3d& AY,
                                             const Vector3d& AZ);

        protected:
            /** rotation axis */
            Vector3d m_axis;
        };
        
        GMDSMath_API Vector3d operator*(const AxisAngleRotation& AR,
                           const Vector3d &AV);
		GMDSMath_API AxisAngleRotation operator*(const AxisAngleRotation& AR0,
                                    const AxisAngleRotation& AR1);
        /*----------------------------------------------------------------------------*/
    } // namespace math
    /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_AXIS_ANGLE_ROTATION_H_ */
/*----------------------------------------------------------------------------*/
