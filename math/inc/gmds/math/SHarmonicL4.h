/*-----------------------------------------------------------------*/
//  SHarmonicL4.h
//  Created by F. Ledoux on 02/11/2015.
/*------------------------------------------------------------------------*/
#ifndef GMDS_MATH_SHARMONICL4_H_
#define GMDS_MATH_SHARMONICL4_H_
/*------------------------------------------------------------------------*/
// GMDS Headers
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Vector.h>
#include <gmds/math/Matrix.h>
#include <gmds/math/Chart.h>
#include <gmds/math/AxisAngleRotation.h>
#include "GMDSMath_export.h"
/*------------------------------------------------------------------------*/
// STL Headers
#include <iostream>
/*------------------------------------------------------------------------*/
namespace gmds{
    /*--------------------------------------------------------------------*/
    namespace math{
        typedef Matrix<9,9,double> SHL4Matrix;
        
        inline SHL4Matrix RZ(const double ATheta){
            const double c1 = cos(    ATheta);
            const double s1 = sin(    ATheta);
            const double c2 = cos(2 * ATheta);
            const double s2 = sin(2 * ATheta);
            const double c3 = cos(3 * ATheta);
            const double s3 = sin(3 * ATheta);
            const double c4 = cos(4 * ATheta);
            const double s4 = sin(4 * ATheta);
            double val[9][9]= {
                {c4, 0, 0, 0, 0, 0, 0, 0, s4},
                {0, c3, 0, 0, 0, 0, 0, s3, 0},
                {0, 0, c2, 0, 0, 0, s2, 0, 0},
                {0, 0, 0, c1, 0, s1, 0, 0, 0},
                {0, 0, 0, 0, 1, 0, 0, 0, 0},
                {0, 0, 0, -s1, 0, c1, 0, 0, 0},
                {0, 0, -s2, 0, 0, 0, c2, 0, 0},
                {0, -s3, 0, 0, 0, 0, 0, c3, 0},
                {-s4, 0, 0, 0, 0, 0, 0, 0, c4}
            };
            
            return SHL4Matrix(val);
        }
        
        inline SHL4Matrix DRZ (const double ATheta) {
            
            const double c1 = cos(    ATheta);
            const double s1 = sin(    ATheta);
            const double c2 = cos(2 * ATheta);
            const double s2 = sin(2 * ATheta);
            const double c3 = cos(3 * ATheta);
            const double s3 = sin(3 * ATheta);
            const double c4 = cos(4 * ATheta);
            const double s4 = sin(4 * ATheta);
            double val[9][9]={
                {-4*s4, 0, 0, 0, 0, 0, 0, 0, 4*c4},
                {0, -3*s3, 0, 0, 0, 0, 0, 3*c3, 0},
                {0, 0, -2*s2, 0, 0, 0, 2*c2, 0, 0},
                {0, 0, 0, -s1, 0, c1, 0, 0, 0},
                {0, 0, 0, 0, 1, 0, 0, 0, 0},
                {0, 0, 0, -c1, 0, -s1, 0, 0, 0},
                {0, 0, -2*c2, 0, 0, 0, -2*s2, 0, 0},
                {0, -3*c3, 0, 0, 0, 0, 0, -3*s3, 0},
                {-4*c4, 0, 0, 0, 0, 0, 0, 0, -4*s4}};
            return SHL4Matrix(val);
        }
        
        inline SHL4Matrix RX90() {
            const double v_14_4 = sqrt(14.0)/4.0;
            const double v_07_4 = sqrt( 7.0)/4.0;
            const double v_05_4 = sqrt( 5.0)/4.0;
            const double v_02_4 = sqrt( 2.0)/4.0;
            const double v_35_8 = sqrt(35.0)/8.0;
            double val[9][9]={
                {0      , 0     , 0      , 0     , 0     , v_14_4, 0      , -v_02_4, 0 }     ,
                {0      , -3./4 , 0      , v_07_4, 0     , 0     , 0      , 0      , 0   }   ,
                {0      , 0     , 0      , 0     , 0     , v_02_4, 0      , v_14_4 , 0  }    ,
                {0      , v_07_4, 0      , 3./4  , 0     , 0     , 0      , 0      , 0      },
                {0      , 0     , 0      , 0     , 3./8  , 0     , v_05_4 , 0      , v_35_8 },
                {-v_14_4, 0     , -v_02_4, 0     , 0     , 0     , 0      , 0      , 0      },
                {0      , 0     , 0      , 0     , v_05_4, 0     , 1./2   , 0      , -v_07_4},
                {v_02_4 , 0     , -v_14_4, 0     , 0     , 0     , 0      , 0      , 0      },
                {0      , 0     , 0      , 0     , v_35_8, 0     , -v_07_4, 0      , 1./8   }};
            return SHL4Matrix(val);
        }
        inline SHL4Matrix RX90Inv() {
            const double v_14_4 = sqrt(14.0)/4.0;
            const double v_07_4 = sqrt( 7.0)/4.0;
            const double v_05_4 = sqrt( 5.0)/4.0;
            const double v_02_4 = sqrt( 2.0)/4.0;
            const double v_35_8 = sqrt(35.0)/8.0;
            double val[9][9]={
                {0      , 0     , 0     , 0     , 0     , -v_14_4, 0      , v_02_4 , 0      },
                {0      , -3./4 , 0     , v_07_4, 0     , 0      , 0      , 0      , 0      },
                {0      , 0     , 0     , 0     , 0     , -v_02_4, 0      , -v_14_4, 0      },
                {0      , v_07_4, 0     , 3./4  , 0     , 0      , 0      , 0      , 0      },
                {0      , 0     , 0     , 0     , 3./8  , 0      , v_05_4 , 0      , v_35_8 },
                {v_14_4 , 0     , v_02_4, 0     , 0     , 0      , 0      , 0      , 0      },
                {0      , 0     , 0     , 0     , v_05_4, 0      , 1./2   , 0      , -v_07_4},
                {-v_02_4, 0     , v_14_4, 0     , 0     , 0      , 0      , 0      , 0      },
                {0      , 0     , 0     , 0     , v_35_8, 0      , -v_07_4, 0      , 1./8   }};
            return SHL4Matrix(val);
        }
        
        inline SHL4Matrix RY (const double ATheta) {
            return RX90Inv() * (RZ(ATheta) * RX90());
        }
        inline SHL4Matrix DRY (const double ATheta) {
            return RX90Inv() * (DRZ(ATheta) * RX90());
        }
        
        inline SHL4Matrix RX (const double ATheta) {
            const double pi = std::acos(-1.0);
            SHL4Matrix RY90 = RY(pi/2);
            SHL4Matrix RY90Inv = RY90.transpose();
            return RY90 * (RZ(ATheta) * RY90Inv);
        }
        
        inline SHL4Matrix DRX (const double ATheta) {
            const double pi = std::acos(-1.0);
            SHL4Matrix RY90 = RY(pi/2);
            SHL4Matrix RY90Inv = RY90.transpose();
            return RY90 * (DRZ(ATheta) * RY90Inv);
        }
        inline SHL4Matrix EX() {
            const double s02 = sqrt(2.0);
            const double s72 = sqrt(7.0/2.0);
            const double s32 = 3.0/sqrt(2.0);
            const double s10 = sqrt(10.0);
            double val[9][9]={
                {0  , 0  , 0  , 0  , 0  , 0  , 0  ,-s02, 0  },
                {0  , 0  , 0  , 0  , 0  , 0  ,-s72, 0  ,-s02},
                {0  , 0  , 0  , 0  , 0  ,-s32, 0  ,-s72, 0  },
                {0  , 0  , 0  , 0  ,-s10, 0  ,-s32, 0  , 0  },
                {0  , 0  , 0  , s10, 0  , 0  , 0  , 0  , 0  },
                {0  , 0  , s32, 0  , 0  , 0  , 0  , 0  , 0  },
                {0  , s72, 0  , s32, 0  , 0  , 0  , 0  , 0  },
                {s02, 0  , s72, 0  , 0  , 0  , 0  , 0  , 0  },
                {0  , s02, 0  , 0  , 0  , 0  , 0  , 0  , 0  }};
            return SHL4Matrix(val);
        }
        inline SHL4Matrix EY() {
            const double s02 = sqrt(2.0);
            const double s72 = sqrt(7.0/2.0);
            const double s32 = 3.0/sqrt(2.0);
            const double s10 = sqrt(10.0);
            double val[9][9]={
                {0   , s02, 0  , 0  , 0  , 0  , 0  , 0  , 0  },
                {-s02, 0  , s72, 0  , 0  , 0  , 0  , 0  , 0  },
                {0   ,-s72, 0  , s32, 0  , 0  , 0  , 0  , 0  },
                {0   , 0  ,-s32, 0  , 0  , 0  , 0  , 0  , 0  },
                {0   , 0  , 0  , 0  , 0  ,-s10, 0  , 0  , 0  },
                {0   , 0  , 0  , 0  , s10, 0  ,-s32, 0  , 0  },
                {0   , 0  , 0  , 0  , 0  , s32, 0  ,-s72, 0  },
                {0   , 0  , 0  , 0  , 0  , 0  , s72, 0  ,-s02},
                {0   , 0  , 0  , 0  , 0  , 0  , 0  , s02, 0  }};
            return SHL4Matrix(val);
        }
        
        inline SHL4Matrix EZ() {
            double val[9][9]={
                {0, 0, 0, 0, 0, 0, 0, 0, 4},
                {0, 0, 0, 0, 0, 0, 0, 3, 0},
                {0, 0, 0, 0, 0, 0, 2, 0, 0},
                {0, 0, 0, 0, 0, 1, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0,-1, 0, 0, 0, 0, 0},
                {0, 0,-2, 0, 0, 0, 0, 0, 0},
                {0,-3, 0, 0, 0, 0, 0, 0, 0},
                {-4, 0, 0, 0, 0, 0, 0, 0, 0}};
            return SHL4Matrix(val);
        }
        /*------------------------------------------------------------*/
        /* Create the 9x9 rotation matrix from XYZ 3D angles
         */
        inline SHL4Matrix fromXYZEulerAngle(const double AX,
                                          const double AY,
                                          const double AZ) {
            return RX(AX)*(RY(AY)*RZ(AZ));
        }
        /*----------------------------------------------------------------*/
        /** \class SHarmonicL4.
         *  \brief Represents a 3D SHarmonicL4 as a function on the 3D sphere
         *         in hthe L4 spherical harmonics basis.
         */
        class GMDSMath_API SHarmonicL4: public Vector9d{
            
        public:
            /*------------------------------------------------------------*/
            /*  \brief Provides the identity SHL4 for us
             */
            static SHarmonicL4 identity();
            /*------------------------------------------------------------*/
            /*  \brief Provides the h0  SHL4
             */
            static SHarmonicL4 h0();
            
            /*------------------------------------------------------------*/
            /*  \brief Provides the h4  SHL4
             */
            static SHarmonicL4 h4();
            
            /*------------------------------------------------------------*/
            /*  \brief Provides the h8  SHL4
             */
            static SHarmonicL4 h8();
            
            /*------------------------------------------------------------*/
            /*  \brief Default constructor provides a null SHarmonicL4
             */
            SHarmonicL4();
            
            /*------------------------------------------------------------*/
            /** \brief Constructs a SHarmonicL4 from an array of 9 values
             *
             * \param[in] ATab a const pointer to the array of 9 values
             */
            SHarmonicL4(const double *ATab);
            
            
            /*------------------------------------------------------------*/
            /** \brief Constructs a SHarmonicL4 from a chart
             *
             * \param[in] AC a chart.
             */
            SHarmonicL4(const Chart& AC);
            
            /*------------------------------------------------------------*/
            /**\brief Constructs a SHarmonicL4 from a series of 9 values
             *
             * \param[in] AV1 first value
             * \param[in] AV2 second value
             * \param[in] AV3 third value
             * \param[in] AV4 fourth value
             * \param[in] AV5 fifth value
             * \param[in] AV6 sixth value
             * \param[in] AV7 seventh value
             * \param[in] AV8 eigth value
             * \param[in] AV9 ninth value
             */

            SHarmonicL4(const double AV1, const double AV2,
                        const double AV3, const double AV4,
                        const double AV5, const double AV6,
                        const double AV7, const double AV8,
                        const double AV9);
            
            /*------------------------------------------------------------*/
            /** \brief Copy constructor
             * \param[in] ASH a const reference a SHarmonicL4 to be copied
             */
            SHarmonicL4(const SHarmonicL4& ASH);
            /*------------------------------------------------------------*/
            /** \brief Copy constructor
             * \param[in] AV a const reference a 9D vector to be copied
             */
            SHarmonicL4(const Vector9d& AV);
            
            
            /*-----------------------------------------------------------*/
            /** \brief Computes the closes SHL4 of \p ATarget
             *
             * \param[in] ATarget the 9D vector we want to be the closest
             *                    as possible.
             * \param[out] ARot the axis-angle rotation used to bring XYZ
             *                  onto the chart corresponding to the 
             *                  resulting spherical harmonic
             *
             * \return The closest true spherical harmonic
             */
            static SHarmonicL4 closest(const Vector9d& ATarget,
                                       AxisAngleRotation& ARot);

            /*-----------------------------------------------------------*/
            /** \brief overloaded operator=.
             *
             * \param[in]  AV a same type generic vector
             * \return     a reference onto *this
             */
            SHarmonicL4& operator= (const SHarmonicL4& ASH) {
                memcpy(m_data, ASH.m_data, dimension * sizeof(double));
                return *this;
            }
            /*------------------------------------------------------------*/
            /** \brief Applies a \p AAlpha rotation around the X axis
             *
             * \param[in]  AAlpha the rotation angle, in radians
             */
            void rotX(const double AAlpha);
            /*------------------------------------------------------------*/
            /** \brief Applies a \p AAlpha rotation around the Y axis
             *
             * \param[in]  AAlpha the rotation angle, in radians
             */
            void rotY(const double AAlpha);
            /*------------------------------------------------------------*/
            /** \brief Applies a \p AAlpha rotation around the Z axis
             *
             * \param[in]  AAlpha the rotation angle, in radians
             */
            void rotZ(const double AAlpha);
            /*------------------------------------------------------------*/
            /** \brief Applies a series of rotation around X Y and Z axis
             *         (in this order)
             *
             * \param[in]  AX the rotation angle around X, in radians
             * \param[in]  AY the rotation angle around Y, in radians
             * \param[in]  AZ the rotation angle around Z, in radians
             */
            void rotXYZ(const double AX,
                        const double AY,
                        const double AZ);
            
            /*------------------------------------------------------------*/
            /** \brief Applies a rotation around the axis \p AAxis with
             *         angle \p AAngle. \p AAxis must be a unit vector!
             * \param[in]  AAxis  a rotation axis
             * \param[in]  AAngle a rotation angle
             */
            void rot(const math::Vector3d& AAxis,
                     const double AAngle);
            
            
            /*------------------------------------------------------------*/
            /** \brief Applies a rotation along the Axis-angle rotation
             *         \p ARot
             * \param[in]  ARot an axis-angle rotation
             */
            void rot(const math::AxisAngleRotation& ARot);
            
            /*------------------------------------------------------------*/
            /** \brief Build a chart corresponding to the rotation 
             *         of ASH
             * 
             * \return A chart
             */
            math::Chart chart() const;
            /*------------------------------------------------------------*/
            /** \brief Friend function to ease stream output
             *
             * \param[in]  AStr an output stream
             * \param[in]  ASH  a L4-basis spherical harmonic
             *
             * \return the modified output stream \p AStr
             */
            friend std::ostream & operator << (std::ostream & AStr, const SHarmonicL4 & AF);
        protected:
            /*------------------------------------------------------------*/
            /** \brief Applies a 90 degree rotation around the X axis
             */
            void rotX90();
            /*------------------------------------------------------------*/
            /** \brief Applies a -90 degree rotation around the X axis
             */
            void rotXMinus90();
            /*------------------------------------------------------------*/
            /** \brief Applies a 90 degree rotation around the Y axis
             */
            void rotY90();
            /*------------------------------------------------------------*/
            /** \brief Applies a -90 degree rotation around the Y axis
             */
            void rotYMinus90();
            
        };
        /*---------------------------------------------------------------*/
    } //namespace math
    /*-----------------------------------------------------------------*/
}//namespace gmds
/*-----------------------------------------------------------------*/
//SHFrame(const Chart& AC)
//{
//    Vector x = AC.X();
//    Vector y = AC.Y();
//    Vector z = AC.Z();
//    m_x = Eigen::Vector3d(x.X(),x.Y(),x.Z());
//    m_y = Eigen::Vector3d(y.X(),y.Y(),y.Z());
//    m_z = Eigen::Vector3d(z.X(),z.Y(),z.Z());
//    
//    Quaternion qi(AC);
//    
//    Eigen::Quaterniond eqi( qi.X(),qi.I(),qi.J(),qi.K() );
//    Eigen::AngleAxisd::Matrix3 di = eqi.matrix();
//    Eigen::Vector3d xyz = di.eulerAngles(0,1,2);
//    
//    double alpha = xyz.x();
//    double beta  = xyz.y();
//    double gamma = xyz.z();
//    
//    //Get  9x1 rotational vectors for ANode
//    math::SHMatrix m = math::SH::fromXYZEulerAngle(alpha,beta,gamma);
//    
//    m_a = m*math::SHFrame::reference().a();
//}

#endif /*GMDS_MATH_SHARMONICL4_H_*/
/*-----------------------------------------------------------------*/
