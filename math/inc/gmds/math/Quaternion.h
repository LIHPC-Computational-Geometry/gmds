/*----------------------------------------------------------------------------*/
/*
 * Quaternion.h
 *
 *  Created on: 03/01/2015
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_QUATERNION_H_
#define GMDS_MATH_QUATERNION_H_
/*----------------------------------------------------------------------------*/
#include <cmath>
#include <iostream>
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Vector.h>
#include <gmds/math/Vector.h>
#include <gmds/math/Chart.h>
#include <gmds/math/Matrix.h>
/*----------------------------------------------------------------------------*/
using namespace std;
/*----------------------------------------------------------------------------*/
namespace gmds {
    /*----------------------------------------------------------------------------*/
    namespace math {
        /*----------------------------------------------------------------------------*/
        /** \class Quaternion
         *  \brief Define a 4D vector, i.e. a quaternion
         */
        class EXPORT_GMDS Quaternion {
            
            
        public:
            
            /*------------------------------------------------------------------------*/
            /** \brief Default Constructor.
             */
            Quaternion():m_r(1.), m_i(0.,0.,0.) {};
            /*------------------------------------------------------------------------*/
            /** \brief Constructor.
             *
             * \param[in] AX real part
             * \param[in] AI first component of the imaginary part
             * \param[in] AJ second component of the imaginary part
             * \param[in] AK third component of the imaginary part
             */
            Quaternion(const TCoord AX, const TCoord AI,
                       const TCoord AJ, const TCoord AK)
            :m_r(AX),m_i(AI,AJ,AK)
            {
                this->normalize();
            };
            
            /*------------------------------------------------------------------------*/
            /** \brief Constructor from a chart
             */
            explicit Quaternion(const Chart & t);
            
            /*------------------------------------------------------------------------*/
            /** \brief Copy constructor
             * \param[in] AQ A quaternion to be copied
             */
            Quaternion(const Quaternion& AQ) : m_r(AQ.m_r), m_i(AQ.m_i) {}
            
            /*------------------------------------------------------------------------*/
            /** \brief Constructs Quaternion from a vector and a scalar
             *
             * \param[in] AR a "real" scalar
             * \param[in] AI an "imaginary" vector
             */
            Quaternion(const TCoord AR, const Vector3d& AI) : m_r(AR), m_i(AI){}
            
            /*------------------------------------------------------------------------*/
            /** \brief Overloaded operator=
             *
             * \param[in] AQ the Quaternion to be copied
             * \return    a reference to (*this)
             */
            Quaternion& operator = ( const Quaternion &AQ ) {
                m_r = AQ.m_r;
                m_i = AQ.m_i;
                return *this ;
            }

            /*------------------------------------------------------------------------*/
            /** \brief Accessors
             */
            TCoord X()const{ return m_r;   }
            TCoord I()const{ return m_i[0]; }
            TCoord J()const{ return m_i[1]; }
            TCoord K()const{ return m_i[2]; }
            
            vector<TCoord> getVal()const;
            /*------------------------------------------------------------------------*/
            /** \brief Gets the imaginary vector component
             *
             * \return the vector component
             */
            const Vector3d& imaginaryPart() const {
                return m_i;
            }
            /*------------------------------------------------------------------------*/
            /** \brief Gets the real scalar component
             *
             * \return the scalar component
             */
            double realPart() const {
                return m_r;
            }

            /*------------------------------------------------------------------------*/
            /** \brief Set the rotation angle.
             *
             * \param[in] AA the rotation angle
             */
            void  setAngle(const TCoord& AA) {
                m_r = cos(AA) / 2.0 ;
                m_i = sin(AA / 2.0)*axis();
            }
            /*------------------------------------------------------------------------*/
            /** \brief Get the rotation angle
             *
             * \return the rotation angle
             */
            double angle() const {
                return 2.0 * acos(m_r) ;
            }
            /*------------------------------------------------------------------------*/
            /**
             * \brief Get the rotation axis.
             * \return The rotation  axis.
             */
            Vector3d axis() const;

            /*------------------------------------------------------------------------*/
            /** \brief Multiplication between (*this) and \p AQ
             *
             * \param[in] AQ a quaternion to multiply with (*this)
             *
             * \return *this *\p AQ
             */
            inline Quaternion operator*(const Quaternion & AQ)const{
                math::Vector3d v1 = m_i, v2 = AQ.m_i;
                double         s1 = m_r, s2 = AQ.m_r;
                
                return Quaternion(s1*s2 - v1.dot(v2), s1*v2 + s2*v1 + v1.cross(v2));
            }
            
//            /*------------------------------------------------------------------------*/
//            /** \brief Difference
//             */
//            TCoord operator -(const Quaternion & AQ)const
//            {
//                return pow((x - AQ.x), 2) + pow((i - AQ.i), 2) + pow((j - AQ.j), 2) +
//                pow((k - AQ.k), 2);
//            }
            /*------------------------------------------------------------------------*/
            /** \brief Sum between (*this) and \p AQ
             *
             * \param[in] AQ a quaternion to add to (*this)
             *
             * \return *this +\p AQ
             */
            inline Quaternion operator +(const Quaternion & AQ)const
            {
                return Quaternion(m_r + AQ.m_r, m_i + AQ.m_i);
            }
            
            /*------------------------------------------------------------------------*/
            /** \brief Difference between (*this) and \p AQ
             *
             * \param[in] AQ a quaternion to substract from (*this)
             *
             * \return *this -\p AQ
             */
            inline  Quaternion operator -(const Quaternion & AQ)const
            {
                return Quaternion(m_r - AQ.m_r, m_i - AQ.m_i);
            }
            
            /// \brief Comparison operator "different form"
            bool operator!=(const Quaternion & AQ)const;
            
            
            inline TCoord dot(const Quaternion & AQ)const
            {
                return (m_r*AQ.m_r + m_i[0]*AQ.m_i[0] +
                        m_i[1]*AQ.m_i[1] + m_i[2]*AQ.m_i[2]);
            };
            
            Quaternion opposite()const {
                return Quaternion( -m_r, -m_i);
            }
            
            Quaternion conjugate()const
            {
                return Quaternion(m_r, -m_i);
            }
            
            TCoord angle(const Quaternion& AQ)const
            {
                double d = dot(AQ.closestImg(*this));
                return  acos(d);
            }
            
            //this method rotate *this to be aligned with AVec. The minimal rotation is
            //performed
            Quaternion alignWith(const math::Vector3d& AVec);
            
            
            /// \brief Recuperation des angles d'Euler correspondant à *this
            void toEulerAngle(TCoord& angleX, TCoord& angleY, TCoord& angleZ) const;
            void setFromEulerAngle(TCoord& angleX, TCoord& angleY, TCoord& angleZ);
            
            void reset(const TCoord X, const TCoord I, const TCoord J, const TCoord K)
            {
                m_r    = X;
                m_i[0] = I;
                m_i[1] = J;
                m_i[2] = K;
                normalize();
            }
            
            Quaternion closestImg(const Quaternion & item)const;
            
            
            /*------------------------------------------------------------------------*/
            /** \brief Computes the spherical interpolation between two Quaternions (SLERP)
             *
             * \param[in] AFrom the first quaternion we come from
             * \param[in] ATO   the second quaternion we go to
             * \param[in] AParam the interpoloation parameter  in [0.0,1.0]
             * \return a Quaternion equals to (\p AParam-1)*\p AFrom +\p AParam*\p ATo
             */
            static Quaternion SLERP(const Quaternion& AFrom,
                                    const Quaternion& ATo,
                                    const double APAram);
            /*------------------------------------------------------------------------*/
            /** \brief A simple pondered mean
             *
             * This methods implement a one-pass direct mean pondered by the
             * weights. Please note that the elements of AQuats are not used directly
             * to compute the mean. Instead it is their image closest to ref which is
             * used.
             *
             * \param AQuats   quaternions we want to compute the mean
             * \param Aweights weights associated to each cross in the mean computation
             * \param ARef     a reference used to compute the closest image
             */
            static EXPORT_GMDS Quaternion simpleMean(const vector<Quaternion> & AQuats,
                                                     const vector<TCoord> & AWeights,
                                                     const Quaternion & ARef);
            
            /*------------------------------------------------------------------------*/
            /** \brief A multi-pass mean
             *
             * This method implement a multi-pass mean. It is intended for first
             * definition of quaternions, when no reference is available. It uses
             * internally simpleMeanQuat multiples times to converge towards a stable
             * mean.
             *
             * \param AQuats   quaternions we want to compute the mean
             * \param Aweights weights associated to each cross in the mean computation
             */
            static EXPORT_GMDS Quaternion mean(const vector<Quaternion> & AQuats,
                                               const vector<TCoord> & AWeights);
            static EXPORT_GMDS Quaternion mean(const Quaternion& AQ1, 
                                               const TCoord& AC1, 
                                               const Quaternion& AQ2, 
                                               const TCoord& AC2);
            
            static EXPORT_GMDS Quaternion getAleat();
            
            
            /*------------------------------------------------------------------------*/
            /** \brief Check consistency between 4 Quaternions
             */
            static EXPORT_GMDS int testSingularity(Quaternion& q1,
                                                   Quaternion& q2,
                                                   Quaternion& q3,
                                                   Quaternion& q4);
            
            /*------------------------------------------------------------------------*/
            /** \brief Check consistency between 3 Quaternions
             */
            static EXPORT_GMDS int testSingularity(Quaternion& q1,
                                                   Quaternion& q2,
                                                   Quaternion& q3);
            
            
            
            /**
             * \brief Converts this Quat into a matrix
             * \return a matrix that represents this Quat
             */
            Matrix<4,4,double> toMatrix() const {
                double t, xs, ys, zs, wx, wy, wz, xx, xy, xz, yy, yz, zz;
                
                t= 2.0/(m_r*m_r+m_i.dot(m_i));
                xs = m_i[0] * t ;
                ys = m_i[1] * t ;
                zs = m_i[2] * t ;
                
                wx = m_r * xs ;
                wy = m_r * ys ;
                wz = m_r * zs ;
                
                xx = m_i[0] * xs ;
                xy = m_i[0] * ys ;
                xz = m_i[0] * zs ;
                
                yy = m_i[1] * ys ;
                yz = m_i[1] * zs ;
                zz = m_i[2] * zs ;
                
                Matrix<4,4,double> matrix ;
                
                matrix(0,0) = 1.0 - (yy+zz) ;
                matrix(1,0) = xy + wz ;
                matrix(2,0) = xz - wy ;
                matrix(0,1) = xy - wz ;
                matrix(1,1) = 1.0 - (xx+zz) ;
                matrix(2,1) = yz+wx ;
                matrix(0,2) = xz + wy ;
                matrix(1,2) = yz - wx ;
                matrix(2,2) = 1.0 - (xx+yy) ;
                
                return matrix;
            }

        private:
            /** real part or scalar part */
            TCoord   m_r;
            /** Pure imaginary part or vectorial part*/
            Vector3d m_i;

            
            inline void normalize()
            {
                const double n = sqrt(m_r*m_r+m_i.dot(m_i));
                if (n != 0.0){
                    m_r /= n;
                    m_i[0]/= n;
                    m_i[1]/= n;
                    m_i[2]/= n;
                }
            }
            static const Quaternion img[24];
        }; 
        /*--------------------------------------------------------------------*/
        EXPORT_GMDS ostream & operator << (ostream & AStr,
                                           const Quaternion & AQ);
        /*--------------------------------------------------------------------*/
  }
    /*------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_QUATERNION_H_ */
/*----------------------------------------------------------------------------*/
