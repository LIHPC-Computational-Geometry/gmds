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
/*-----------------------------------------------------------------*/
/*
 * cross.cpp
 *
 *  Created on: Sep 05, 2014
 *      Author: franck Ledoux
 */
/*-----------------------------------------------------------------*/
#include "gmds/math/SHarmonicL4.h"
#include "gmds/math/AxisAngleRotation.h"
/*-----------------------------------------------------------------*/
// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
/*-----------------------------------------------------------------*/
namespace gmds {
  /*---------------------------------------------------------------*/
    namespace math{
        const double sh_pi = std::acos(-1.0);
        const double sh_sq2 = sqrt(2.0);
        const double sh_sq72 = sqrt(7.0/2.0);
        const double sh_3_sq2 = 3.0/sqrt(2.0);
        const double sh_sq10 = sqrt(10.0);
        const double sh_sq14_4 = sqrt(14.0)/4.0;
        const double sh_sq7_4 = sqrt( 7.0)/4.0;
        const double sh_sq5_4 = sqrt( 5.0)/4.0;
        const double sh_sq2_4 = sqrt( 2.0)/4.0;
        const double sh_sq35_8 = sqrt(35.0)/8.0;
        /*---------------------------------------------------------*/
        SHarmonicL4 SHarmonicL4::identity(){
            return SHarmonicL4(0, 0, 0, 0, sqrt(7.0/12.0),
                               0, 0, 0, sqrt(5.0/12.0));
        }
        /*---------------------------------------------------------*/
        SHarmonicL4 SHarmonicL4::h0(){
            return SHarmonicL4(1, 0, 0, 0, 0, 0, 0, 0, 0);
        }
        /*---------------------------------------------------------*/
        SHarmonicL4 SHarmonicL4::h4(){
            return SHarmonicL4(0, 0, 0, 0, 1, 0, 0, 0, 0);
        }
        /*---------------------------------------------------------*/
        SHarmonicL4 SHarmonicL4::h8(){
            return SHarmonicL4(0, 0, 0, 0, 0, 0, 0, 0, 1);
        }
        /*---------------------------------------------------------*/
        SHarmonicL4::SHarmonicL4(): Vector9d() {
            m_data[0] = sqrt(5.0/12.0);
            m_data[1] = 0;
            m_data[2] = 0;
            m_data[3] = 0;
            m_data[4] = sqrt(7.0/12.0);
            m_data[5] = 0;
            m_data[6] = 0;
            m_data[7] = 0;
            m_data[8] = sqrt(5.0/12.0);
        }
        /*---------------------------------------------------------*/
        SHarmonicL4::SHarmonicL4(const double *ATab):
        Vector9d(ATab) {;}
        /*---------------------------------------------------------*/
        SHarmonicL4::SHarmonicL4(const double AV1, const double AV2,
                     const double AV3, const double AV4,
                     const double AV5, const double AV6,
                     const double AV7, const double AV8,
                     const double AV9) {
            m_data[0] = AV1;
            m_data[1] = AV2;
            m_data[2] = AV3;
            m_data[3] = AV4;
            m_data[4] = AV5;
            m_data[5] = AV6;
            m_data[6] = AV7;
            m_data[7] = AV8;
            m_data[8] = AV9;
        }
        /*---------------------------------------------------------*/
        SHarmonicL4::SHarmonicL4(const Chart& AC):Vector9d()
        {
            m_data[0] = 0;
            m_data[1] = 0;
            m_data[2] = 0;
            m_data[3] = 0;
            m_data[4] = sqrt(7.0/12.0);
            m_data[5] = 0;
            m_data[6] = 0;
            m_data[7] = 0;
            m_data[8] = sqrt(5.0/12.0);

            Vector3d y (AC.Y()[0],AC.Y()[1], AC.Y()[2]);
            Vector3d z (AC.Z()[0],AC.Z()[1], AC.Z()[2]);
            
            AxisAngleRotation r = AxisAngleRotation::alignYZ(y,z);
            
            rot(r);
            
            //  OTHER OPTION
//            //compute Euler angles via Quaternion representation
//            double alpha, beta, gamma;
//            math::Quaternion q(AC);
//            q.toEulerAngle(alpha,beta, gamma);
//            rotXYZ(alpha,beta,gamma);
        }
        

        /*---------------------------------------------------------*/
        SHarmonicL4::SHarmonicL4(const SHarmonicL4& AF):Vector9d(AF)
        {;}
        /*---------------------------------------------------------*/
        SHarmonicL4::SHarmonicL4(const Vector9d& AF):Vector9d(AF)
        {;}
        /*---------------------------------------------------------*/
        math::Chart SHarmonicL4::chart()const{
            AxisAngleRotation r;
            SHarmonicL4::closest(*this, r);
            
            Vector3d x(1,0,0);
            Vector3d y(0,1,0);
            Vector3d z(0,0,1);
            return Chart(r*x,r*y,r*z);
        }
        /*---------------------------------------------------------*/
        void SHarmonicL4::rotXYZ(const double AX,
                                 const double AY,
                                 const double AZ)
        {
            rotX(AX);
            rotY(AY);
            rotZ(AZ);
        }
        /*---------------------------------------------------------*/
        void SHarmonicL4::rot(const math::Vector3d& AAxis,
                                 const double AAngle){
            double z_angle = AAngle;
            if (z_angle<.00001)
                return; //too small rotation, we do nothing
            
            double x = AAxis.X();
            double y = AAxis.Y();
            double z = AAxis.Z();
            
            double d = sqrt(y*y + z*z);
            
            double x_angle = (fabs(d)<.00001) ? 0 : atan2(y/d, z/d);
            double y_angle = atan2(-x/z_angle, d/z_angle);
            rotX( x_angle);
            rotY( y_angle);
            rotZ( z_angle);
            rotY(-y_angle);
            rotX(-x_angle);
        }
        /*---------------------------------------------------------*/
        void SHarmonicL4::
        rot(const math::AxisAngleRotation& ARot){
            rot(ARot.axis(), ARot.angle());
        }

        /*---------------------------------------------------------*/
        void SHarmonicL4::rotX(const double AAlpha) {
            rotYMinus90();
            rotZ(AAlpha);
            rotY90();
        }
        /*---------------------------------------------------------*/
        void SHarmonicL4::rotY(const double AAlpha) {
            rotX90();
            rotZ(AAlpha);
            rotXMinus90();
        }
        /*---------------------------------------------------------*/
        void SHarmonicL4::rotZ(const double AAlpha) {
            SHarmonicL4 old(*this);
            for (int i=0; i<4; i++) {
                m_data[i]   =  old[i]*cos((4-i)*AAlpha) + old[8-i]*sin((4-i)*AAlpha);
                m_data[8-i] = -old[i]*sin((4-i)*AAlpha) + old[8-i]*cos((4-i)*AAlpha);
            }
        }
        

        /*---------------------------------------------------------*/
        void SHarmonicL4::rotX90() {
            SHarmonicL4 old(*this);
            m_data[0] =  sh_sq14_4*old[5] -sh_sq2_4*old[7];
            m_data[1] = -0.75*old[1]+ sh_sq7_4*old[3];
            m_data[2] = sh_sq2_4*old[5]+sh_sq14_4*old[7];
            m_data[3] = sh_sq7_4*old[1]+0.75*old[3];
            m_data[4] = 0.375*old[4]+sh_sq5_4*old[6]+sh_sq35_8*old[8];
            m_data[5] = -sh_sq14_4*old[0] -sh_sq2_4*old[2];
            m_data[6] =  sh_sq5_4*old[4]+  0.5*old[6] -sh_sq7_4*old[8];
            m_data[7] = sh_sq2_4*old[0] -sh_sq14_4*old[2];
            m_data[8] =  sh_sq35_8*old[4]  -sh_sq7_4*old[6]+0.125*old[8];
        }
        
        /*---------------------------------------------------------*/
        void SHarmonicL4::rotXMinus90() {
            SHarmonicL4 old(*this);
            m_data[0] = -sh_sq14_4*old[5] + sh_sq2_4*old[7];
            m_data[1] = -0.75*old[1]+ sh_sq7_4*old[3];
            m_data[2] = -sh_sq2_4*old[5] -sh_sq14_4*old[7];
            m_data[3] = sh_sq7_4*old[1]+0.75*old[3];
            m_data[4] = 0.375*old[4]+sh_sq5_4*old[6]+sh_sq35_8*old[8];
            m_data[5] =  sh_sq14_4*old[0]+  sh_sq2_4*old[2];
            m_data[6] =  sh_sq5_4*old[4]+  0.5*old[6] -sh_sq7_4*old[8];
            m_data[7] = -sh_sq2_4*old[0]+  sh_sq14_4*old[2];
            m_data[8] =  sh_sq35_8*old[4]  -sh_sq7_4*old[6]+0.125*old[8];
        }

        /*---------------------------------------------------------*/
        void SHarmonicL4::rotY90() {
            
            const double a = 0.35355339059327373;
            const double b = 0.9354143466934854;
            const double c = 0.5590169943749475;
            const double d = 0.7395099728874521;
            const double e = 0.6614378277661477;
            SHarmonicL4 old(*this);
            m_data[0] = +a*old[1] +b*old[3];
            m_data[1] = -a*old[0] -b*old[2];
            m_data[2] = +b*old[1] -a*old[3];
            m_data[3] = -b*old[0] +a*old[2];
            m_data[4] =  0.375*old[4] -c*old[6] + d*old[8];
            m_data[5] =   0.75*old[5] -e*old[7];
            m_data[6] = -c*old[4] +0.5*old[6] + e*old[8];
            m_data[7] = -e*old[5] -0.75*old[7];
            m_data[8] =  d*old[4] +e*old[6] + 0.125*old[8];
        }

        /*---------------------------------------------------------*/
        void SHarmonicL4::rotYMinus90() {
            const double a = 0.35355339059327373;
            const double b = 0.9354143466934854;
            const double c = 0.5590169943749475;
            const double d = 0.7395099728874521;
            const double e = 0.6614378277661477;
            
            SHarmonicL4 old(*this);
            m_data[0] = -a*old[1] -b*old[3];
            m_data[1] =  a*old[0] +b*old[2];
            m_data[2] = -b*old[1] +a*old[3];
            m_data[3] =  b*old[0] -a*old[2];
            m_data[4] =  0.375*old[4] -c*old[6] + d*old[8];
            m_data[5] =  0.75*old[5] -e*old[7];
            m_data[6] = -c*old[4] +0.5*old[6] + e*old[8];
            m_data[7] = -e*old[5] -0.75*old[7];
            m_data[8] =  d*old[4] +e*old[6] + 0.125*old[8];
        }


        /*---------------------------------------------------------*/
        SHarmonicL4 SHarmonicL4::closest(const Vector9d& ATarget,
                                         AxisAngleRotation& ARot){
            double step      = 0.1;
            double grad_threshold = 1e-3;
            math::Vector9d q = ATarget;
            
            if(q.norm()==0){
                std::cout<<"Warning: null 9D vector"<<std::endl;
                return SHarmonicL4(0,0,0,0,std::sqrt(7. / 12.),
                                   0,0,0,std::sqrt(5. / 12.));
            }
            q.normalize();
            
            //================================================
            // Starting from the reference frame can lead to
            // wrong situation (due to 45 degree invariant)
            // So we enforce the descent direction into the
            // closest solution
            double init_SH[5][9] = {
                {0,0,0,0,std::sqrt(7. / 12.),0,0,0,std::sqrt(5. / 12.)},
                {0,0,0,0,-0.190941,0,-0.853913,0,0.484123},
                {0,0,0,0,-0.190941,0,0.853913,0,0.484123},
                {0,0,0,0,0.763763,0,0,0,-0.645497},
                {0,0,-0.853913,0,-0.190941,0,0,0,-0.484123}
            };
            
            Vector3d init_AARot[5] = {
                Vector3d(0,0,0),
                Vector3d(0.785398,0,0),
                Vector3d(0,0.785398,0),
                Vector3d(0,0,0.785398),
                Vector3d(0.743782,0.308085,0.743782)
            };
            
            //we look for the best starting point??
            int      i_best = 0;
            Vector9d v_test = Vector9d(init_SH[0]) - q;
            auto     d_best = v_test.norm();
            for (auto i = 1; i < 5; i++) {
                v_test = Vector9d(init_SH[i]) - q;
                auto d = v_test.norm();
                if (d<d_best) {
                    d_best = d;
                    i_best = i;
                }
            }
            
            
            ARot = AxisAngleRotation(init_AARot[i_best]);
            SHarmonicL4 sh(init_SH[i_best]);
            
            //    math::Vector3d first_gradient(qt*math::SH::EX()*a,
            //                                qt*math::SH::EY()*a,
            //                                qt*math::SH::EZ()*a);
            ////    std::cout<<"--------------------------"<<std::endl;
            //    if(first_gradient.norm()<threshold){
            //  //      std::cout<<"Change init"<<std::endl;
            //        f.rotXYZ(0,0,0.1);
            //    }
            
            auto compt=0;
            while(true){
                compt++;
                // Computation are unfolded for
                //EX_sh = math::EX()*sh;
                //EY_sh = math::EY()*sh;
                //EZ_sh = math::EZ()*sh;
                math::SHarmonicL4 EX_sh(-sh_sq2*sh[7],
                                        -sh_sq2*sh[8]-sh_sq72*sh[6],
                                        -sh_sq72*sh[7]-sh_3_sq2*sh[5],
                                        -sh_3_sq2*sh[6]-sh_sq10*sh[4],
                                        sh_sq10*sh[3],
                                        sh_3_sq2*sh[2],
                                        sh_sq72*sh[1]+sh_3_sq2*sh[3],
                                        sh_sq2*sh[0]+sh_sq72*sh[2],
                                        sh_sq2*sh[1]);
                
                math::SHarmonicL4 EY_sh(sh_sq2*sh[1],
                                        -sh_sq2*sh[0]+sh_sq72*sh[2],
                                        -sh_sq72*sh[1]+sh_3_sq2*sh[3],
                                        -sh_3_sq2*sh[2],
                                        -sh_sq10*sh[5],
                                        -sh_3_sq2*sh[6]+sh_sq10*sh[4],
                                        -sh_sq72*sh[7]+sh_3_sq2*sh[5],
                                        -sh_sq2*sh[8]+sh_sq72*sh[6],
                                        sh_sq2*sh[7]);
                
                math::SHarmonicL4 EZ_sh(4*sh[8],3*sh[7],
                                        2*sh[6],sh[5],0,
                                        -sh[3],-2*sh[2],
                                        -3*sh[1],-4*sh[0]);
                
                
                math::Vector3d gradient(q.dot(EX_sh),
                                        q.dot(EY_sh),
                                        q.dot(EZ_sh));
                
                if(gradient.norm()<grad_threshold)
                    break;
                
                if (compt>=100000) {
                    std::cout<< "Infinite loop for a projection" << std::endl;
                    return sh;
                }

                gradient = step*gradient;
                sh.rotXYZ(gradient.X(),
                          gradient.Y(),
                          gradient.Z());
                
                ARot = AxisAngleRotation(math::Vector3d(1,0,0),gradient[0])*ARot;
                ARot = AxisAngleRotation(math::Vector3d(0,1,0),gradient[1])*ARot;
                ARot = AxisAngleRotation(math::Vector3d(0,0,1),gradient[2])*ARot;
                
                
            }
            //    Eigen::Matrix<double, 1, 9>  qt=q.transpose();
            //
            //    math::Vector first_gradient(qt*math::SH::EX()*a,
            //                                qt*math::SH::EY()*a,
            //                                qt*math::SH::EZ()*a);
            //    if(first_gradient.norm()<threshold)
            //        f.rotXYZ(0,0,0.3);
            //
            //    bool not_conv = true;
            //    while(not_conv){
            //
            //        math::Vector gradient(qt*math::SH::EX()*a,
            //                              qt*math::SH::EY()*a,
            //                              qt*math::SH::EZ()*a);
            //        if(gradient.norm()<threshold)
            //            not_conv = false;
            //        else{
            //            f.rotXYZ(step*gradient.X(),
            //                        step*gradient.Y(),
            //                        step*gradient.Z());
            //
            //            a= f.a();
            //        }
            //    }
            
            //    }
            //    else{//non linear version
            //
            //        ALGLIB_SH_Target = AFrom;
            //        alglib::real_1d_array euler_angles = "[0, 0, 0]";
            //        alglib::mincgstate state;
            //        alglib::mincgreport rep;
            //        alglib::mincgcreate(euler_angles, state);
            //        alglib::mincgoptimize(state, evalFuncAndGrad);
            //        alglib::mincgresults(state, euler_angles, rep);
            //
            //        f.rotXYZ(euler_angles[0],
            //                    euler_angles[1],
            //                    euler_angles[2]);
            //
            //    }
            return sh;
        }
        
        

        /*-----------------------------------------------------------------*/
        
        std::ostream & operator << (std::ostream & AStr, const SHarmonicL4 & AF)
        {
            AStr << "[" << AF.m_data[0] << ")";
            return AStr;
        }
    }
    /*-----------------------------------------------------------------*/
}
/*-----------------------------------------------------------------*/
