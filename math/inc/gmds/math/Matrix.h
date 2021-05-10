/*----------------------------------------------------------------------------*/
/** \file    Matrix.h
 *  \author  F. LEDOUX
 *  \date    14/06/2011
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_MATRIX_H_
#define GMDS_MATH_MATRIX_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <cmath>
#include <cassert>
#include <string.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*------------------------------------------------------------------------*/
    namespace math{
        /*--------------------------------------------------------------------*/
#include <gmds/math/Matrix_def.h>
        /*--------------------------------------------------------------------*/
        /** \struct DetTraits
         *  \brief  Template structure providing a tailored process for 
         *          computing the determinant of an TLinxTCol matrix
         *
         *  \tparam TLin	number of lines
         *  \tparam TCol	number of columns
         *  \tparam TType	value type
         */
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType> struct DetTraits {
            static TType perform(const TType t[TLin][TCol])
            {throw GMDSException();}
        };
        /*-------------------------------------------------------------------*/
        template<typename TType> struct DetTraits<2, 2,TType>  {
            static TType perform(const TType t[2][2])
            {
                return  t[0][0] * t[1][1] - t[1][0] * t[0][1];
            }
        };
        /*-------------------------------------------------------------------*/
        template<typename TType> struct DetTraits<3, 3, TType> {
            static TType perform(const TType t[3][3]) {
                return (t[0][0] * (t[1][1] * t[2][2] - t[1][2] * t[2][1]) -
                        t[0][1] * (t[1][0] * t[2][2] - t[1][2] * t[2][0]) +
                        t[0][2] * (t[1][0] * t[2][1] - t[1][1] * t[2][0]));
            }
        };
        /*-------------------------------------------------------------------*/
        template<typename TType> struct DetTraits<4, 4, TType> {
            static TType perform(const TType t[4][4]) {
                
                return (
                        t[0][0] * (t[1][1] * (t[2][2] * t[3][3] - t[3][2] * t[2][3]) -
                                   t[1][2] * (t[2][1] * t[3][3] - t[3][1] * t[2][3]) +
                                   t[1][3] * (t[2][1] * t[3][2] - t[3][1] * t[2][2]))
                        -
                        t[0][1] * (t[1][0] * (t[2][2] * t[3][3] - t[3][2] * t[2][3]) -
                                   t[1][2] * (t[2][0] * t[3][3] - t[3][0] * t[2][3]) +
                                   t[1][3] * (t[2][0] * t[3][2] - t[3][0] * t[2][2]))
                        +
                        t[0][2] * (t[1][0] * (t[2][1] * t[3][3] - t[3][1] * t[2][3]) -
                                   t[1][1] * (t[2][0] * t[3][3] - t[3][0] * t[2][3]) +
                                   t[1][3] * (t[2][0] * t[3][1] - t[3][0] * t[2][1]))
                        -
                        t[0][3] * (t[1][0] * (t[2][1] * t[3][2] - t[3][1] * t[2][2]) -
                                   t[1][1] * (t[2][0] * t[3][2] - t[3][0] * t[2][2]) +
                                   t[1][2] * (t[2][0] * t[3][1] - t[3][0] * t[2][1])));
            }
        };
        /*-------------------------------------------------------------------*/
        /** \struct InvTraits
         *  \brief  Template structure providing a tailored process for 
         *          computing the inverse matrix of an TLinxTCol matrix
         *
         *  \tparam TLin	number of lines
         *  \tparam TCol	number of columns
         *  \tparam TType   value type
         */
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType> struct InvTraits {
            static bool perform(const TType t[TLin][TCol],
                                TType(&i)[TLin][TCol]) {
                throw GMDSException("Not yet implemented");
            }
        };
        /*-------------------------------------------------------------------*/
        template<typename TType> struct InvTraits<3, 3, TType>
        {
            static bool perform(const TType t[3][3], TType(&i)[3][3])
            {
                double det = DetTraits<3, 3, TType>::perform(t);
                if (det){
                    double inv_det = 1. / det;
                    i[0][0] =  (t[1][1] * t[2][2] - t[1][2] * t[2][1]) * inv_det;
                    i[1][0] = -(t[1][0] * t[2][2] - t[1][2] * t[2][0]) * inv_det;
                    i[2][0] =  (t[1][0] * t[2][1] - t[1][1] * t[2][0]) * inv_det;
                    
                    i[0][1] = -(t[0][1] * t[2][2] - t[0][2] * t[2][1]) * inv_det;
                    i[1][1] =  (t[0][0] * t[2][2] - t[0][2] * t[2][0]) * inv_det;
                    i[2][1] = -(t[0][0] * t[2][1] - t[0][1] * t[2][0]) * inv_det;
                    
                    i[0][2] =  (t[0][1] * t[1][2] - t[0][2] * t[1][1]) * inv_det;
                    i[1][2] = -(t[0][0] * t[1][2] - t[0][2] * t[1][0]) * inv_det;
                    i[2][2] =  (t[0][0] * t[1][1] - t[0][1] * t[1][0]) * inv_det;
                    
                    return true;
                }
                
                return false;
            }
        };
        /*-------------------------------------------------------------------*/
        template<typename TType> struct InvTraits<4, 4, TType>
        {
            static bool perform(const TType t[4][4], TType(&i)[4][4])
            {
                i[0][0] =
                t[1][2] * t[2][3] * t[3][1] -
                t[1][3] * t[2][2] * t[3][1] +
                t[1][3] * t[2][1] * t[3][2] -
                t[1][1] * t[2][3] * t[3][2] -
                t[1][2] * t[2][1] * t[3][3] +
                t[1][1] * t[2][2] * t[3][3];
                
                i[0][1] =
                t[0][3] * t[2][2] * t[3][1] -
                t[0][2] * t[2][3] * t[3][1] -
                t[0][3] * t[2][1] * t[3][2] +
                t[0][1] * t[2][3] * t[3][2] +
                t[0][2] * t[2][1] * t[3][3] -
                t[0][1] * t[2][2] * t[3][3];
                
                i[0][2] =
                t[0][2] * t[1][3] * t[3][1] -
                t[0][3] * t[1][2] * t[3][1] +
                t[0][3] * t[1][1] * t[3][2] -
                t[0][1] * t[1][3] * t[3][2] -
                t[0][2] * t[1][1] * t[3][3] +
                t[0][1] * t[1][2] * t[3][3];
                
                i[0][3] =
                t[0][3] * t[1][2] * t[2][1] -
                t[0][2] * t[1][3] * t[2][1] -
                t[0][3] * t[1][1] * t[2][2] +
                t[0][1] * t[1][3] * t[2][2] +
                t[0][2] * t[1][1] * t[2][3] -
                t[0][1] * t[1][2] * t[2][3];
                
                i[1][0] =
                t[1][3] * t[2][2] * t[3][0] -
                t[1][2] * t[2][3] * t[3][0] -
                t[1][3] * t[2][0] * t[3][2] +
                t[1][0] * t[2][3] * t[3][2] +
                t[1][2] * t[2][0] * t[3][3] -
                t[1][0] * t[2][2] * t[3][3];
                
                i[1][1] =
                t[0][2] * t[2][3] * t[3][0] -
                t[0][3] * t[2][2] * t[3][0] +
                t[0][3] * t[2][0] * t[3][2] -
                t[0][0] * t[2][3] * t[3][2] -
                t[0][2] * t[2][0] * t[3][3] +
                t[0][0] * t[2][2] * t[3][3];
                
                i[1][2] =
                t[0][3] * t[1][2] * t[3][0] -
                t[0][2] * t[1][3] * t[3][0] -
                t[0][3] * t[1][0] * t[3][2] +
                t[0][0] * t[1][3] * t[3][2] +
                t[0][2] * t[1][0] * t[3][3] -
                t[0][0] * t[1][2] * t[3][3];
                
                i[1][3] =
                t[0][2] * t[1][3] * t[2][0] -
                t[0][3] * t[1][2] * t[2][0] +
                t[0][3] * t[1][0] * t[2][2] -
                t[0][0] * t[1][3] * t[2][2] -
                t[0][2] * t[1][0] * t[2][3] +
                t[0][0] * t[1][2] * t[2][3];
                
                i[2][0] =
                t[1][1] * t[2][3] * t[3][0] -
                t[1][3] * t[2][1] * t[3][0] +
                t[1][3] * t[2][0] * t[3][1] -
                t[1][0] * t[2][3] * t[3][1] -
                t[1][1] * t[2][0] * t[3][3] +
                t[1][0] * t[2][1] * t[3][3];
                
                i[2][1] =
                t[0][3] * t[2][1] * t[3][0] -
                t[0][1] * t[2][3] * t[3][0] -
                t[0][3] * t[2][0] * t[3][1] +
                t[0][0] * t[2][3] * t[3][1] +
                t[0][1] * t[2][0] * t[3][3] -
                t[0][0] * t[2][1] * t[3][3];
                
                i[2][2] =
                t[0][1] * t[1][3] * t[3][0] -
                t[0][3] * t[1][1] * t[3][0] +
                t[0][3] * t[1][0] * t[3][1] -
                t[0][0] * t[1][3] * t[3][1] -
                t[0][1] * t[1][0] * t[3][3] +
                t[0][0] * t[1][1] * t[3][3];
                i[2][3] =
                t[0][3] * t[1][1] * t[2][0] -
                t[0][1] * t[1][3] * t[2][0] -
                t[0][3] * t[1][0] * t[2][1] +
                t[0][0] * t[1][3] * t[2][1] +
                t[0][1] * t[1][0] * t[2][3] -
                t[0][0] * t[1][1] * t[2][3];
                i[3][0] =
                t[1][2] * t[2][1] * t[3][0] -
                t[1][1] * t[2][2] * t[3][0] -
                t[1][2] * t[2][0] * t[3][1] +
                t[1][0] * t[2][2] * t[3][1] +
                t[1][1] * t[2][0] * t[3][2] -
                t[1][0] * t[2][1] * t[3][2];
                i[3][1] =
                t[0][1] * t[2][2] * t[3][0] -
                t[0][2] * t[2][1] * t[3][0] +
                t[0][2] * t[2][0] * t[3][1] -
                t[0][0] * t[2][2] * t[3][1] -
                t[0][1] * t[2][0] * t[3][2] +
                t[0][0] * t[2][1] * t[3][2];
                i[3][2] =
                t[0][2] * t[1][1] * t[3][0] -
                t[0][1] * t[1][2] * t[3][0] -
                t[0][2] * t[1][0] * t[3][1] +
                t[0][0] * t[1][2] * t[3][1] +
                t[0][1] * t[1][0] * t[3][2] -
                t[0][0] * t[1][1] * t[3][2];
                i[3][3] =
                t[0][1] * t[1][2] * t[2][0] -
                t[0][2] * t[1][1] * t[2][0] +
                t[0][2] * t[1][0] * t[2][1] -
                t[0][0] * t[1][2] * t[2][1] -
                t[0][1] * t[1][0] * t[2][2] +
                t[0][0] * t[1][1] * t[2][2];
                
                TType det = DetTraits<4, 4, TType>::perform(t);
                
                if (det == 0)
                    return false;
                
                det = 1.0 / det;
                
                for (int j = 0; j < 4; j++)
                    for (int k = 0; k < 4; k++)
                        i[j][k] = i[j][k] * det;
                
                return true;
                
            }
        };
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType>
        Matrix<TLin, TCol, TType>::Matrix()
        {
            for (unsigned int i = 0; i < TLin; i++)
            {
                for (unsigned int j = 0; j < TCol; j++)
                    m_data[i][j] = 0;
            }
        }
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType>
        Matrix<TLin, TCol, TType>::Matrix(TType m[TLin][TCol])
        {
            for (unsigned int i = 0; i < TLin; i++)
                memcpy(&m_data[i][0], &m[i][0], TCol*sizeof(TType));
        }
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType>
        Matrix<TLin, TCol, TType>::Matrix(const Matrix<TLin, TCol, TType>& mat)
        {
            for (unsigned int i = 0; i < TLin; i++)
                memcpy(&m_data[i][0], &mat.m_data[i][0], TCol*sizeof(TType));
        }
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType>
        Matrix<TLin, TCol, TType>& Matrix<TLin, TCol, TType>::operator=(const Matrix<TLin, TCol, TType>& mat)
        {
            if (mat == *this)
                return *this;
            
            for (unsigned int i = 0; i < TLin; i++)
                memcpy(&m_data[i][0], &mat.m_data[i][0], TCol*sizeof(TType));
            
            return *this;
        }
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType>
        bool Matrix<TLin, TCol, TType>::operator==(const Matrix<TLin, TCol, TType>& mat) const
        {
            bool result = true;
            for (unsigned int i = 0; i < TLin; i++)
                for (unsigned int j = 0; j < TCol; j++)
                    if (m_data[i][j] != mat.m_data[i][j])
                        return false;
            
            return result;
        }
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType> Matrix<TLin, TCol, TType>::~Matrix()
        {}
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType> TType Matrix<TLin, TCol, TType>::det() const
        {
            return DetTraits<TLin, TCol, TType>::perform(m_data);
        }
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType>
        TType Matrix<TLin, TCol, TType>::get(const int& i, const int& j) const
        {
            assert(i>=0 && i<TLin && j>=0 && j<TCol);
            return m_data[i][j];
        }
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType>
        void Matrix<TLin, TCol, TType>::set(const int& i, const int& j, const TType& AVal)
        {
            assert(i>=0 && i<TLin && j>=0 && j<TCol);
            m_data[i][j] = AVal;
        }
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType>
        void Matrix<TLin, TCol, TType>::set(const TType param[TLin][TCol])
        {
            for (unsigned int i = 0; i < TLin; i++)
                memcpy(&m_data[i][0], &param[i][0], TCol*sizeof(TType));
        }
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType>
        Matrix<TCol, TLin, TType> Matrix<TLin, TCol, TType>::transpose() const
        {
            Matrix<TCol, TLin, TType> T;
            for (int i = 0; i < TCol; i++)
                for (int j = 0; j < TLin; j++)
                    T.set(i, j, (*this).get(j, i));
            
            return T;
        }
	/*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType>
        TType Matrix<TLin, TCol, TType>::
        frobeniusNorm2() const {

		TType frob = 0;
		for(int i=0; i<TLin; i++) {
	                for(int j=0; j<TCol; j++) {
                   	 	frob += m_data[i][j] * m_data[i][j];
			}
		}
		return frob;
	}
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType>
        Matrix<TLin, TCol, TType> Matrix<TLin, TCol, TType>::zero()
        {
            TType null_vect[TLin][TCol];
            for (int i = 0; i < TLin; i++)
                for (int j = 0; j < TCol; j++)
                    null_vect[i][j] = 0.0;
            
            return Matrix<TLin, TCol, TType>(null_vect);
        }
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType>
        Matrix<TLin, TCol, TType> Matrix<TLin, TCol, TType>::identity()
        {
            if (TLin != TCol)
                throw GMDSException("Only NxN matrices have an identity!!");
            
            Matrix<TLin, TCol, TType> T = zero();
            for (int i = 0; i < TLin; i++)
                T.set(i, i, 1.0);
            
            return T;
        }
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType>
        Matrix<TLin, TCol, TType> Matrix<TLin, TCol, TType>::inverse() const
        {
            TType inv_tab[TLin][TCol];
            bool done = InvTraits<TLin, TCol, TType>::perform(m_data, inv_tab);
            if (!done) {
                throw GMDSException("Matrix::inverse - Null Determinant");
            }
            return Matrix<TLin, TCol, TType>(inv_tab);
        }
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType>
        VectorND<TLin, TType>
        Matrix<TLin, TCol, TType>::solve(const VectorND<TLin, TType>& AB) const
        {
            Matrix<TLin, TCol, TType> Ainv = inverse();
            
            VectorND<TLin, TType> x;
            
            for(int i=0; i<TLin; i++) {
                x[i] = 0.0;
                for(int j=0; j<TLin; j++) {
                    x[i] += Ainv(i,j)*AB[j];
                }
            }
            return x;
        }
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType>
        VectorND<3,TType> Matrix<TLin, TCol, TType>::
        eulerAngles(const int& AR0, const int& AR1, const int& AR2) const {
            if(TLin!=3 || TCol!=3)
                throw GMDSException("Only work for 3x3 rotation matrices");

            VectorND<3,TType> coeff;
            
            const int odd = ((AR0+1)%3 == AR1) ? 0 : 1;
            const int i = AR0;
            const int j = (AR0 + 1 + odd)%3;
            const int k = (AR0 + 2 - odd)%3;
            
            if (AR0==AR2)
            {
                coeff[0] = atan2(m_data[j][i], m_data[k][i]);
                TType s2 = Vector2d(m_data[j][i],m_data[k][i]).norm();
                
                if(( odd   && coeff[0]<TType(0)) ||
                   ((!odd) && coeff[0]>TType(0))) {
                    coeff[0] = (coeff[0] > TType(0)) ? coeff[0]-TType(M_PI) : coeff[0]+TType(M_PI);
                    coeff[1] = -atan2(s2, m_data[i][i]);
                }
                else  {
                    coeff[1] = atan2(s2, m_data[i][i]);
                }
                
                TType s1 = sin(coeff[0]);
                TType c1 = cos(coeff[0]);
                coeff[2] = atan2(c1*m_data[j][k] - s1*m_data[k][k],
                                 c1*m_data[j][j] - s1*m_data[k][j]);
            }
            else {
                coeff[0] = atan2(m_data[j][k], m_data[k][k]);
                double c2_x = m_data[i][i];
                double c2_y = m_data[i][j];
                TType c2bis = Vector3d(c2_x, c2_y,0).norm();
                TType c2 = Vector2d(c2_x, c2_y).norm();

                if(  (odd  && coeff[0]<TType(0)) ||
                   ((!odd) && coeff[0]>TType(0))) {
                    coeff[0] = (coeff[0] > TType(0)) ? coeff[0]-TType(M_PI) : coeff[0]+TType(M_PI);
                    coeff[1] = atan2(-m_data[i][k], -c2);
                }
                else
                    coeff[1] = atan2(-m_data[i][k], c2);
                TType s1 = sin(coeff[0]);
                TType c1 = cos(coeff[0]);
                coeff[2] = atan2(s1*m_data[k][i] - c1*m_data[j][i],
                                 c1*m_data[j][j] - s1*m_data[k][j]);
            }
            if (!odd)
                coeff = -coeff;
            
            return coeff;
        }
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType>
        std::ostream& operator<<(std::ostream& stream, const Matrix<TLin, TCol, TType>& mat)
        {
            for (int i = 0; i < TLin; i++)
            {
                stream << "[ ";
                for (int j = 0; j < TCol; j++)
                    stream << mat(i,j) << " ";
                stream << "] ";
            }
            
            return stream;
        }
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType> Matrix<TLin, TCol, TType>
        operator-(const Matrix<TLin, TCol, TType>& a, const Matrix<TLin, TCol, TType>& b)
        {
            Matrix<TLin, TCol, TType> m;
            for (int i = 0; i < TLin; i++)
                for (int j = 0; j < TCol; j++)
                    m(i,j) = a(i,j) - b(i,j);
            return m;
        }
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType> Matrix<TLin, TCol, TType>
        operator+(const Matrix<TLin, TCol, TType>& a, const Matrix<TLin, TCol, TType>& b)
        {
            Matrix<TLin, TCol, TType> m;
            for (int i = 0; i < TLin; i++)
                for (int j = 0; j < TCol; j++)
                    m(i,j) = a(i,j) + b(i,j);
            return m;
        }
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCommon, int TCol, typename TType>
        Matrix<TLin, TCol, TType>
        operator*(const Matrix<TLin, TCommon, TType>& a,
                  const Matrix<TCommon, TCol, TType>& b)
        {
            Matrix<TLin, TCol, TType> m;
            for (int i = 0; i < TLin; i++)
                for (int j = 0; j < TCol; j++)
                {
                    TType val = 0.0;
                    for (int k = 0; k < TCommon; k++){
                        val = val + a(i, k)*b(k, j);
                    }
                    m(i, j)= val;
                    
                }
            return m;
        }
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType>
        VectorND<TLin, TType>
        operator*(const Matrix<TLin, TCol, TType>& a,
                  const VectorND<TCol, TType>& v)
        {
            VectorND<TLin, TType> r;
            for (int i = 0; i < TLin; i++)
                for (int j = 0; j < TCol; j++) {
                    r[i] += a(i, j)*v[j];
                }
            return r;
        }
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType>
        Matrix<TLin, TCol, TType>
        operator*(const Matrix<TLin, TCol, TType>& a, const TType& k)
        {
            Matrix<TLin, TCol, TType> m;
            for (int i = 0; i < TLin; i++)
                for (int j = 0; j < TCol; j++)
                    m(i,j) = k*a(i,j);
            return m;
        }
        
        /*-------------------------------------------------------------------*/
        template<int TLin, int TCol, typename TType>
        Matrix<TLin, TCol, TType>
        operator/(const Matrix<TLin, TCol, TType>& a, const TType& k)
        {
            Matrix<TLin, TCol, TType> m;
            for (int i = 0; i < TLin; i++)
                for (int j = 0; j < TCol; j++)
                    m(i,j) = a(i,j)/k;
            return m;
        }
        /*-------------------------------------------------------------------*/
        template<int I, int J, typename T> Matrix<I, J, T>
        operator*(const T& k, const Matrix<I, J, T>& a)
        {
            Matrix<I, J, T> m;
            for (int i = 0; i < I; i++)
                for (int j = 0; j < J; j++)
                    m(i,j) = k*a(i,j);
            return m;
        }

        /*-------------------------------------------------------------------*/
    } // namespace math
    /*------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_MATRIX_H_ */
/*----------------------------------------------------------------------------*/
