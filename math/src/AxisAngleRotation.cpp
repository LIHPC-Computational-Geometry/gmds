/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
 * AxisAngleRotation.cpp
 */
/*----------------------------------------------------------------------------*/
#include <gmds/math/AxisAngleRotation.h>
#include <gmds/math/Quaternion.h>

/*-------------------------------------------------------------------------*/
namespace gmds{
    /*------------------------------------------------------------------------*/
    namespace math{
        /*--------------------------------------------------------------------*/
        AxisAngleRotation::AxisAngleRotation(){
            double val[3] = {0,0,0};
            m_axis=Vector3d(val);
        }
        
        /*--------------------------------------------------------------------*/
        AxisAngleRotation::
        AxisAngleRotation(const Vector3d& AV, const double AA)
        : m_axis(AA*AV[0],AA*AV[1],AA*AV[2]){
        }
        
        /*--------------------------------------------------------------------*/
        AxisAngleRotation::
        AxisAngleRotation(const Vector3d& AV)
        : m_axis(AV){
        }
        /*--------------------------------------------------------------------*/
        Quaternion AxisAngleRotation::quaternion() const
        {
            double angle = m_axis.norm();
            if (angle<1e-5) {
                return Quaternion(1.0, 0.0, 0.0, 0.0);
            }
            Vector3d n = m_axis;
            n.normalize();
            Vector3d v = sin(angle/2.0)*n;
            return Quaternion(cos(angle/2.0),v);
        }
        
        /*--------------------------------------------------------------------*/
        Matrix<3,3,double> AxisAngleRotation::toRotationMatrix() const
        {
            // Rodiguez formulae
            // R = I + (sin a)K + (1 − cos a)K*K
            // with K the cross product matrix of the rotation angle
            
            Vector3d a = m_axis;
            a.normalize();
            double K_val[3][3]= {
                { 0.0  ,-a.Z(), a.Y()},
                { a.Z(), 0.0  ,-a.X()},
                {-a.Y(),a.X() , 0.0   }};
            
            Matrix<3,3,double> K(K_val);
            Matrix<3,3,double> K2 = K*K;
            Matrix<3,3,double> I = Matrix<3,3,double>::identity();
            
            double angle = m_axis.norm();
            return I+ sin(angle)*K + (1-cos(angle))*K2;
        }
        
        
        /*--------------------------------------------------------------------*/
        Chart AxisAngleRotation::toChart() const
        {
            return Chart((*this)*math::Vector3d(1,0,0),
                         (*this)*math::Vector3d(0,1,0),
                         (*this)*math::Vector3d(0,0,1));
        }
        
        /*--------------------------------------------------------------------*/
        AxisAngleRotation::AxisAngleRotation(const Quaternion& AQ)
        {
            double real_part = AQ.realPart();
            double bounded_real_part = (real_part>-1.0)?real_part:-1.0;
            double v = (bounded_real_part<1.0)?bounded_real_part:1.0;
            double alpha = 2.0*acos(v);
            
            if (AQ.imaginaryPart().norm()<1.0e-5) {
                m_axis=Vector3d(0., 0., 0.);
            }
            else {
                Vector3d ax =AQ.imaginaryPart();
                ax.normalize();
                m_axis = alpha*ax;
            }
        }
        
        /*--------------------------------------------------------------------*/
        AxisAngleRotation::AxisAngleRotation(const Chart& AC)
        {
            Vector3d y (AC.Y()[0],AC.Y()[1], AC.Y()[2]);
            Vector3d z (AC.Z()[0],AC.Z()[1], AC.Z()[2]);
            AxisAngleRotation r = AxisAngleRotation::alignYZ(y,z);
            m_axis = r.m_axis;
        }
        /*--------------------------------------------------------------------*/
        AxisAngleRotation::AxisAngleRotation(const Vector3d& AFrom,
                                             const Vector3d& ATo)
        {
            Vector3d rotation_vector(0.0, 0.0, 0.0);
            
            if ((AFrom+ATo).norm() < 1.0e-2) {
                // vectors are opposite
                rotation_vector = Vector3d(0.0,M_PI,0.0);
            }
            else {
                rotation_vector = AFrom.cross(ATo);
                double min      = rotation_vector.norm();
                
                if(min>1.)
                    min =1.;
                
                double as = ::asin(min);
                rotation_vector.normalize();

                if (ATo.dot(AFrom)>0)
                    rotation_vector = as * rotation_vector;
                else
                    rotation_vector= (M_PI-as) * rotation_vector;
                
            }
            m_axis=rotation_vector;
            
        }
        
        /*--------------------------------------------------------------------*/
        AxisAngleRotation AxisAngleRotation::alignZ(const Vector3d& AV)
        {
            return AxisAngleRotation(Vector3d(0., 0., 1.), AV);
        }
        
        Vector3d AxisAngleRotation::toRotationAxis(const int AIndex) const{
            Quaternion q = quaternion();
            Vector3d v = q.imaginaryPart();
            double s = q.realPart();
            double t, xs, ys, zs, wx, wy, wz, xx, xy, xz, yy, yz, zz;
            t = 2.0 / (v.dot( v) + (s* s));
            
            xs = v.X() * t;
            ys = v.Y() * t;
            zs = v.Z() * t;
            
            wx = s * xs;
            wy = s * ys;
            wz = s * zs;
            
            xx = v.X() * xs;
            xy = v.X() * ys;
            xz = v.X() * zs;
            
            yy = v.Y() * ys;
            yz = v.Y() * zs;
            zz = v.Z() * zs;
            
            switch(AIndex) {
                case 0:
                    return Vector3d(1.0 - (yy + zz), xy + wz, xz - wy);
                case 1:
                    return Vector3d(xy - wz, 1.0 - (xx + zz), yz + wx);
                case 2:
                    return Vector3d(xz + wy, yz - wx, 1.0 - (xx + yy));
            }
            return Vector3d(0.0, 0.0, 0.0);
        }
        /*--------------------------------------------------------------------*/
        AxisAngleRotation AxisAngleRotation::alignYZ(const Vector3d& AY,
                                                     const Vector3d& AZ){
            AxisAngleRotation res =AxisAngleRotation::alignZ(AZ);
            Vector3d rotation_y = res * Vector3d(0.0, 1.0, 0.0);

            double angle = atan2(AZ.dot(AY.cross(rotation_y)),
                                 AY.dot(rotation_y)          );
            if(angle!=0)
                res = res*AxisAngleRotation(-angle*Vector3d(0.0, 0.0, 1.0));
            return res;
        }
        /*--------------------------------------------------------------------*/
        Vector3d operator*(const AxisAngleRotation& AR,
                           const Vector3d &AV) {
            Matrix<4,4, double> M = AR.quaternion().toMatrix();
            return Vector3d(
                        M(0,0)*AV[0] + M(0,1)*AV[1] + M(0,2)*AV[2],
                        M(1,0)*AV[0] + M(1,1)*AV[1] + M(1,2)*AV[2],
                        M(2,0)*AV[0] + M(2,1)*AV[1] + M(2,2)*AV[2]
                        );
        }
        /*--------------------------------------------------------------------*/
        AxisAngleRotation operator*(const AxisAngleRotation& AR0,
                                    const AxisAngleRotation& AR1){
            return AxisAngleRotation(AR0.quaternion()*AR1.quaternion());
        }

    } // namespace math
    /*------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
