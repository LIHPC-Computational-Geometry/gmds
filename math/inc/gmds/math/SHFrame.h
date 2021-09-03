/*----------------------------------------------------------------------------*/
/*
 * SHFrame.h
 *
 *  Created on: 11/03/2015
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_FRAME_SH_H_
#define GMDS_MATH_FRAME_SH_H_
/*----------------------------------------------------------------------------*/
#include <gmds/math/Chart.h>
#include <gmds/math/Quaternion.h>
/*----------------------------------------------------------------------------*/
#include <Eigen/Sparse>
#include <Eigen/Geometry>
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*------------------------------------------------------------------------*/
    namespace math{
        typedef Eigen::Matrix<double, 9, 1> SHVector;
        typedef Eigen::Matrix<double, 9, 9> SHMatrix;
        typedef Eigen::Vector3d  ZYZVector;
        
        /*----------------------------------------------------------------*/
        class SH{
            
        public:
            static inline SHMatrix RZ (const double ATheta) {
                
                const double c1 = cos(    ATheta);
                const double s1 = sin(    ATheta);
                const double c2 = cos(2 * ATheta);
                const double s2 = sin(2 * ATheta);
                const double c3 = cos(3 * ATheta);
                const double s3 = sin(3 * ATheta);
                const double c4 = cos(4 * ATheta);
                const double s4 = sin(4 * ATheta);
                SHMatrix r;
                r <<
                c4, 0, 0, 0, 0, 0, 0, 0, s4,
                0, c3, 0, 0, 0, 0, 0, s3, 0,
                0, 0, c2, 0, 0, 0, s2, 0, 0,
                0, 0, 0, c1, 0, s1, 0, 0, 0,
                0, 0, 0, 0, 1, 0, 0, 0, 0,
                0, 0, 0, -s1, 0, c1, 0, 0, 0,
                0, 0, -s2, 0, 0, 0, c2, 0, 0,
                0, -s3, 0, 0, 0, 0, 0, c3, 0,
                -s4, 0, 0, 0, 0, 0, 0, 0, c4;
                return r;
            }
            
            static inline SHMatrix DRZ (const double ATheta) {
                
                const double c1 = cos(    ATheta);
                const double s1 = sin(    ATheta);
                const double c2 = cos(2 * ATheta);
                const double s2 = sin(2 * ATheta);
                const double c3 = cos(3 * ATheta);
                const double s3 = sin(3 * ATheta);
                const double c4 = cos(4 * ATheta);
                const double s4 = sin(4 * ATheta);
                SHMatrix r;
                r <<
                -4*s4, 0, 0, 0, 0, 0, 0, 0, 4*c4,
                0, -3*s3, 0, 0, 0, 0, 0, 3*c3, 0,
                0, 0, -2*s2, 0, 0, 0, 2*c2, 0, 0,
                0, 0, 0, -s1, 0, c1, 0, 0, 0,
                0, 0, 0, 0, 1, 0, 0, 0, 0,
                0, 0, 0, -c1, 0, -s1, 0, 0, 0,
                0, 0, -2*c2, 0, 0, 0, -2*s2, 0, 0,
                0, -3*c3, 0, 0, 0, 0, 0, -3*s3, 0,
                -4*c4, 0, 0, 0, 0, 0, 0, 0, -4*s4;
                return r;
            }
            
            static inline SHMatrix RX90() {
                const double v_14_4 = sqrt(14.0)/4.0;
                const double v_07_4 = sqrt( 7.0)/4.0;
                const double v_05_4 = sqrt( 5.0)/4.0;
                const double v_02_4 = sqrt( 2.0)/4.0;
                const double v_35_8 = sqrt(35.0)/8.0;
                SHMatrix r;
                r <<
                0      , 0     , 0      , 0     , 0     , v_14_4, 0      , -v_02_4, 0      ,
                0      , -3./4 , 0      , v_07_4, 0     , 0     , 0      , 0      , 0      ,
                0      , 0     , 0      , 0     , 0     , v_02_4, 0      , v_14_4 , 0      ,
                0      , v_07_4, 0      , 3./4  , 0     , 0     , 0      , 0      , 0      ,
                0      , 0     , 0      , 0     , 3./8  , 0     , v_05_4 , 0      , v_35_8 ,
                -v_14_4, 0     , -v_02_4, 0     , 0     , 0     , 0      , 0      , 0      ,
                0      , 0     , 0      , 0     , v_05_4, 0     , 1./2   , 0      , -v_07_4,
                v_02_4 , 0     , -v_14_4, 0     , 0     , 0     , 0      , 0      , 0      ,
                0      , 0     , 0      , 0     , v_35_8, 0     , -v_07_4, 0      , 1./8   ;
                return r;
            }
            static inline SHMatrix RX90Inv() {
                const double v_14_4 = sqrt(14.0)/4.0;
                const double v_07_4 = sqrt( 7.0)/4.0;
                const double v_05_4 = sqrt( 5.0)/4.0;
                const double v_02_4 = sqrt( 2.0)/4.0;
                const double v_35_8 = sqrt(35.0)/8.0;
                SHMatrix r;
                r <<
                0      , 0     , 0     , 0     , 0     , -v_14_4, 0      , v_02_4 , 0      ,
                0      , -3./4 , 0     , v_07_4, 0     , 0      , 0      , 0      , 0      ,
                0      , 0     , 0     , 0     , 0     , -v_02_4, 0      , -v_14_4, 0      ,
                0      , v_07_4, 0     , 3./4  , 0     , 0      , 0      , 0      , 0      ,
                0      , 0     , 0     , 0     , 3./8  , 0      , v_05_4 , 0      , v_35_8 ,
                v_14_4 , 0     , v_02_4, 0     , 0     , 0      , 0      , 0      , 0      ,
                0      , 0     , 0     , 0     , v_05_4, 0      , 1./2   , 0      , -v_07_4,
                -v_02_4, 0     , v_14_4, 0     , 0     , 0      , 0      , 0      , 0      ,
                0      , 0     , 0     , 0     , v_35_8, 0      , -v_07_4, 0      , 1./8   ;
                return r;
            }
            
            static inline SHMatrix RY (const double ATheta) {
     //            return RX90() * RZ(ATheta) * RX90Inv(); //RAY
                return RX90Inv() * RZ(ATheta) * RX90(); //HUANG
            }
            static inline SHMatrix DRY (const double ATheta) {
               return RX90Inv() * DRZ(ATheta) * RX90();//HUANG
              //  return RX90() * DRZ(ATheta) * RX90Inv();//RAY
            }
            
            static inline SHMatrix RX (const double ATheta) {
                const double pi = std::acos(-1.0);
                SHMatrix RY90 = RY(pi/2);
                SHMatrix RY90Inv = RY90.transpose();
               // return RY90Inv * RZ(ATheta) * RY90;//RAY
                return RY90 * RZ(ATheta) * RY90Inv;
            }
            
            static inline SHMatrix DRX (const double ATheta) {
                const double pi = std::acos(-1.0);
                SHMatrix RY90 = RY(pi/2);
                SHMatrix RY90Inv = RY90.transpose();
                return RY90 * DRZ(ATheta) * RY90Inv;// Huang
//                return RY90Inv * DRZ(ATheta) * RY90; //RAY
            }
            

            static inline SHMatrix EX() {
                const double s02 = sqrt(2.0);
                const double s72 = sqrt(7.0/2.0);
                const double s32 = 3.0/sqrt(2.0);
                const double s10 = sqrt(10.0);
                SHMatrix r;
                r <<
                0  , 0  , 0  , 0  , 0  , 0  , 0  ,-s02, 0  ,
                0  , 0  , 0  , 0  , 0  , 0  ,-s72, 0  ,-s02,
                0  , 0  , 0  , 0  , 0  ,-s32, 0  ,-s72, 0  ,
                0  , 0  , 0  , 0  ,-s10, 0  ,-s32, 0  , 0  ,
                0  , 0  , 0  , s10, 0  , 0  , 0  , 0  , 0  ,
                0  , 0  , s32, 0  , 0  , 0  , 0  , 0  , 0  ,
                0  , s72, 0  , s32, 0  , 0  , 0  , 0  , 0  ,
                s02, 0  , s72, 0  , 0  , 0  , 0  , 0  , 0  ,
                0  , s02, 0  , 0  , 0  , 0  , 0  , 0  , 0  ;
                return r;
            }
            static inline SHMatrix EY() {
                const double s02 = sqrt(2.0);
                const double s72 = sqrt(7.0/2.0);
                const double s32 = 3.0/sqrt(2.0);
                const double s10 = sqrt(10.0);
                SHMatrix r;
                r <<
                0   , s02, 0  , 0  , 0  , 0  , 0  , 0  , 0  ,
                -s02, 0  , s72, 0  , 0  , 0  , 0  , 0  , 0  ,
                0   ,-s72, 0  , s32, 0  , 0  , 0  , 0  , 0  ,
                0   , 0  ,-s32, 0  , 0  , 0  , 0  , 0  , 0  ,
                0   , 0  , 0  , 0  , 0  ,-s10, 0  , 0  , 0  ,
                0   , 0  , 0  , 0  , s10, 0  ,-s32, 0  , 0  ,
                0   , 0  , 0  , 0  , 0  , s32, 0  ,-s72, 0  ,
                0   , 0  , 0  , 0  , 0  , 0  , s72, 0  ,-s02,
                0   , 0  , 0  , 0  , 0  , 0  , 0  , s02, 0  ;
                return r;
            }
            
            static inline SHMatrix EZ() {
                SHMatrix r;
                r <<
                 0, 0, 0, 0, 0, 0, 0, 0, 4,
                 0, 0, 0, 0, 0, 0, 0, 3, 0,
                 0, 0, 0, 0, 0, 0, 2, 0, 0,
                 0, 0, 0, 0, 0, 1, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0,-1, 0, 0, 0, 0, 0,
                 0, 0,-2, 0, 0, 0, 0, 0, 0,
                 0,-3, 0, 0, 0, 0, 0, 0, 0,
                -4, 0, 0, 0, 0, 0, 0, 0, 0;
                return r;
            }
            
            static inline SHVector H0() {
                SHVector h;
                h << 1, 0, 0, 0, 0, 0, 0, 0, 0;
                return h;
            }
            static inline SHVector H4() {
                SHVector h;
                h << 0, 0, 0, 0, 1, 0, 0, 0, 0;
                return h;
            }
            static inline SHVector H8() {
                SHVector h;
                h << 0, 0, 0, 0, 0, 0, 0, 0, 1;
                return h;
            }

            static inline SHVector Identity() {
                SHVector h;
                h << 0, 0, 0, 0, sqrt(7.0), 0, 0, 0, sqrt(5.0);
                return h;
            }
            
//
//            /*------------------------------------------------------------*/
//            /* Create the 9x9 rotation matrix from XYZ 3D angles
//             */
//            static inline SHVector XYZToSH(const ZYZVector& AZYZ) {
//                return RZ(AZYZ[2]) * RX90Inv() * RZ(AZYZ[1]) * RX90() * RZ(AZYZ[0]) * Identity();
//            }
//
            
            
            /*------------------------------------------------------------*/
            /* Create the 9x9 rotation matrix from XYZ 3D angles
             */
            static inline SHMatrix fromXYZEulerAngle(const double AX,
                                                     const double AY,
                                                     const double AZ) {
               return RX(AX)*RY(AY)*RZ(AZ);
            }
            
            
//            inline ZYZVector project(SHVector& sh) {
//                auto zyz = sh2zyz(sh);
//                sh = zyz2sh(zyz);
//                return zyz;
//            }
        };
        
        
        /*----------------------------------------------------------------*/
        class SHFrame{
        public:
            
            /*------------------------------------------------------------*/
            /* Constructor
             *  Defines a frame equal to the reference frame
             */

            SHFrame():m_a(SH::Identity()),
            m_x(Eigen::Vector3d(1,0,0)),
            m_y(Eigen::Vector3d(0,1,0)),
            m_z(Eigen::Vector3d(0,0,1)){;}
            
            SHFrame(const SHFrame& AF){
                m_a = AF.m_a;
                m_x = AF.m_x;
                m_y = AF.m_y;
                m_z = AF.m_z;
            }
            SHFrame(const Chart& AC)
            {
                Vector x = AC.X();
                Vector y = AC.Y();
                Vector z = AC.Z();
                m_x = Eigen::Vector3d(x.X(),x.Y(),x.Z());
                m_y = Eigen::Vector3d(y.X(),y.Y(),y.Z());
                m_z = Eigen::Vector3d(z.X(),z.Y(),z.Z());
                
                Quaternion qi(AC);
                
                Eigen::Quaterniond eqi( qi.X(),qi.I(),qi.J(),qi.K() );
                Eigen::AngleAxisd::Matrix3 di = eqi.matrix();
                Eigen::Vector3d xyz = di.eulerAngles(0,1,2);
                
                double alpha = xyz.x();
                double beta  = xyz.y();
                double gamma = xyz.z();
                
                //Get  9x1 rotational vectors for ANode
                math::SHMatrix m = math::SH::fromXYZEulerAngle(alpha,beta,gamma);

                m_a = m*math::SHFrame::reference().a();
            }
            /*------------------------------------------------------------*/
            /* Provide the reference frame
             */
            SHFrame static reference(){ return SHFrame();}
            SHFrame static approx_ref(){
                SHFrame f;
                f.rotateXYZ(0.1,0.005,0.07);
                return f;
            }
            SHVector a() {return m_a;}
            
            Eigen::Vector3d X() {return m_x;}
            Eigen::Vector3d Y() {return m_y;}
            Eigen::Vector3d Z() {return m_z;}
            
            
            math::Vector VX() {return math::Vector(m_x.x(),m_x.y(),m_x.z());}
            math::Vector VY() {return math::Vector(m_y.x(),m_y.y(),m_y.z());}
            math::Vector VZ() {return math::Vector(m_z.x(),m_z.y(),m_z.z());}

            void rotateXYZ(const double& alpha,
                           const double& beta ,
                           const double& gamma)
            {
                Eigen::Matrix3d R;
                R = Eigen::AngleAxisd(gamma, Eigen::Vector3d::UnitZ())
                  * (Eigen::AngleAxisd(beta , Eigen::Vector3d::UnitY())
                  * Eigen::AngleAxisd(alpha, Eigen::Vector3d::UnitX()));
                
                m_x = R * m_x;
                m_y = R * m_y;
                m_z = R * m_z;
                
                m_a = SH::RZ(gamma)*(SH::RY(beta)* (SH::RX(alpha)*m_a));
            }
            
        //private:
            SHVector m_a;
            Eigen::Vector3d m_x;
            Eigen::Vector3d m_y;
            Eigen::Vector3d m_z;
            
        };
        
    }
    /*------------------------------------------------------------------------*/

}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_FRAME_SH_H_ */
/*----------------------------------------------------------------------------*/
