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
#include <cstring>
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

		  /*----------------------------------------------------------------------------*/
		  /** \class Matrix Nx2
	  *  \brief template class implementing a full mathematical matrix
	  *
	  *  \tparam TRow   number of rows
	  *  \tparam TCol   number of columns
	  *  \tparam TType type of elements stored in the matrix
			*/
		  /*----------------------------------------------------------------------------*/
		  template <int TRow, int TCol, typename TType>
		  struct Matrix {

			  /** type of matrix components */
			  typedef TType value_type;


	        /*------------------------------------------------------------------------*/
			  /** \brief  Accessor to the row \param AI.
				*/

	        VectorND<TCol, TType>& operator[] (const int AI){
		        assert(AI>=0 && AI<TRow);
		        return static_cast<VectorND<TCol, double> &>(m_rows[AI]);
	        }
				const VectorND<TCol, TType>& operator[] (const int AI) const {
		         assert(AI>=0 && AI<TRow);
		         return m_rows[AI];
	         }
	         
	         /*------------------------------------------------------------------------*/
	         /** \brief  Accessor to the column \param AI.
	          */
	         VectorND<TRow, TType> col(const int AI) const {
		         assert(AI>=0 && AI<TCol);
		         VectorND<TRow, TType> ret;
		         for (int i=TRow; i--; ret[i]=m_rows[i][AI]);
		         return ret;
	         }
	         
	         void set_col(const int AI, const VectorND<TRow, TType>  &AV) {
		         assert(AI>=0 && AI<TCol);
		         for (int i=TRow; i--; m_rows[i][AI]=AV[i]);
	         }
			  /*------------------------------------------------------------------------*/
			  /** \brief  Creation of the null matrix
				*/
			  static Matrix<TRow, TCol,TType> zero();

			  /*------------------------------------------------------------------------*/
			  /** \brief  Creation of the identity matrix
				*/
			  static Matrix<TRow, TCol,TType> identity();


	        Matrix<TRow-1, TCol-1,TType> get_minor(const int ARow, const int ACol) const {
		        Matrix<TRow-1, TCol-1,TType> res;
		        for (int i=TRow-1; i--; )
					  for (int j=TCol-1;j--; res[i][j]=m_rows[i<ARow?i:i+1][j<ACol?j:j+1]);
				  return res;
			  }

			  double cofactor(const int ARow, const int ACol) const {
				  return get_minor(ARow, ACol).det()*((ARow+ACol)%2 ? -1 : 1);
			  }

	        Matrix<TRow, TCol,TType>  adjugate() const {
		        Matrix<TRow, TCol,TType> res;
				  for (int i=TRow; i--; )
					  for (int j=TCol; j--; res[i][j]=cofactor(i,j));
				  return res;
			  }

	        Matrix<TRow, TCol,TType> inverse_transpose() const {
		        Matrix<TRow, TCol,TType> res = adjugate();
				  return res/(res[0].dot(m_rows[0]));
			  }


	        /*------------------------------------------------------------------------*/
	        /** \brief  Creation of the inverse matrix
				*/
	        Matrix<TRow, TCol,TType> inverse() const {
				  return inverse_transpose().transpose();
			  }

	        /*------------------------------------------------------------------------*/
	        /** \brief  Computes the transpose matrix
	         */
	        Matrix<TRow, TCol,TType> transpose() const {
		        Matrix<TRow, TCol,TType> res;
				  for (int i=TCol; i--; res[i]=this->col(i));
				  return res;
			  }

			  /*------------------------------------------------------------------------*/
			  /** \brief Solve (*this)* X = \p AB with X the vector solution returned by
			*         this method. This method computes the inverse of *this. It is
			*         only implemented for 3x3 and 4x4 matrices.
			*
			* \param[in] AB the right-hand side vector, which must have \p TRow
			*            components
			*
			* \return The TRow-dimension solution vector
				*/
			  math::VectorND<TRow, TType> solve(const math::VectorND<TRow, TType>& AB) const;

			  /*------------------------------------------------------------------------*/
			  /** \brief  Read-only Access to element (i,j)
		  *
		  *	\param i line number
		  *	\param j column number
		  *
		  *  \return the (i,j) th element of *this
		  *
		  *  \assert if i<0 or i>TRow or j<0 or j>TCol
				*/
			  value_type get(const int& i, const int& j) const;


			  /*------------------------------------------------------------------------*/
			  /** \brief  Read-only Access to element (i,j)
			*
			*	\param i line number
			*	\param j column number
			*
			*  \return the (i,j) th element of *this
			*
			*  \assert if i<0 or i>TRow or j<0 or j>TCol
				*/
			  value_type operator()(const int& i, const int& j) const {
				  return  m_rows[i][j];
			  }

			  /*------------------------------------------------------------------------*/
			  /** \brief  Read and Write Access to element (i,j)
			*
			*	\param i line number
			*	\param j column number
			*
			*  \return the (i,j) th element of *this
			*
			*  \assert if i<0 or i>TRow or j<0 or j>TCol
				*/
			  value_type& operator()(const int& i, const int& j) {
				  return  (m_rows[i])[j];
			  }

			  /*------------------------------------------------------------------------*/
			  /** \brief  Return a tabular view of the matrix
				*/
			  void getTab(value_type (&t)[TRow][TCol]) const
			  {
				  for (unsigned int i = 0; i < TRow; i++)
				  {
					  for (unsigned int j = 0; j < TCol; j++)
						  t[i][j] = m_rows[i][j];
				  }
			  }


			  /*------------------------------------------------------------------------*/
			  /** \brief  Set the elemment (i,j)
		  *
		  *	\param i line number
		  *	\param j column number
		  *	\param AVal the new value for element (i,j)
		  *
		  *  \exception if i<0 or i>T1 or j<0 or j>T2
				*/
			  void set(const int& i, const int& j, const value_type& AVal);

			  /*------------------------------------------------------------------------*/
			  /** \brief  Computes the matrix determinant
		  *
		  *  \return the determinant value
				*/
			  value_type det() const;


			  /*------------------------------------------------------------------------*/
			  /** \brief  Computes the squared forbenius norm of the matrix
				*/
			  value_type frobeniusNorm2() const;


			  /*------------------------------------------------------------------------*/
			  /** \brief Returns the Euler-angles of the 3x3 rotation matrix \c *this
			*         from the parameters (\p AR0,\p AR1,\p AR2) where each parameter
			*         defines a rotation axis as an integer in {0,1,2}. 0 means X,
			*         1 means Y and 2 means Z. So writing
			*         \code M.eulerAngles(0, 2, 1); \endcode
			*         means to perform the X rotation then the Z rotation and finally
			*         the Y rotation. Warning this function must only be used for 3x3
			*         rotation matrix in practice. The implementation comes from
			*         Graphics Gems IV.
			*
			* \param[in] AR0 First  rotation axis
			* \param[in] AR1 Second rotation axis
			* \param[in] AR2 Third  rotation axis
			*
			* \return A vector of corresponding rotation angles. Component i
			*         corresponds to the i^th parameter
				*/
			  math::VectorND<3,value_type> eulerAngles(const int& AR0, const int& AR1,
																	  const int& AR2) const;


			  VectorND<TCol,TType> m_rows[TRow]={};
		  };
		  /*----------------------------------------------------------------------------*/

        /*--------------------------------------------------------------------*/
        /** \struct DetTraits
         *  \brief  Template structure providing a tailored process for 
         *          computing the determinant of an TRowxTCol matrix
         *
         *  \tparam TRow	number of lines
         *  \tparam TCol	number of columns
         *  \tparam TType	value type
         */
        /*-------------------------------------------------------------------*/
        template<int TRow, int TCol, typename TType> struct DetTraits {
            static TType perform(const VectorND<TCol,TType> (&t)[TRow])
            {throw GMDSException();}
        };
        /*-------------------------------------------------------------------*/
        template<typename TType> struct DetTraits<2, 2,TType>  {
            static TType perform(const VectorND<2,TType> (&t)[2])
            {
                return  t[0][0] * t[1][1] - t[1][0] * t[0][1];
            }
        };
        /*-------------------------------------------------------------------*/
        template<typename TType> struct DetTraits<3, 3, TType> {
            static TType perform(const VectorND<3,TType> (&t)[3]) {
                return (t[0][0] * (t[1][1] * t[2][2] - t[1][2] * t[2][1]) -
                        t[0][1] * (t[1][0] * t[2][2] - t[1][2] * t[2][0]) +
                        t[0][2] * (t[1][0] * t[2][1] - t[1][1] * t[2][0]));
            }
        };
        /*-------------------------------------------------------------------*/
        template<typename TType> struct DetTraits<4, 4, TType> {
            static TType perform(const VectorND<4,TType> (&t)[4]) {
                
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
         *          computing the inverse matrix of an TRowxTCol matrix
         *
         *  \tparam TRow	number of lines
         *  \tparam TCol	number of columns
         *  \tparam TType   value type
         */
        /*-------------------------------------------------------------------*/
        template<int TRow, int TCol, typename TType> struct InvTraits {
            static bool perform(const TType t[TRow][TCol],
                                TType(&i)[TRow][TCol]) {
                throw GMDSException("Not yet implemented");
            }
        };
        /*-------------------------------------------------------------------*/
        template<typename TType> struct InvTraits<3, 3, TType>
        {
            static bool perform(const TType t[3][3], TType(&i)[3][3])
            {
                double det = DetTraits<3, 3, TType>::perform(t);
                if (det!=0){
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
                
                for (auto & j : i)
                    for (int k = 0; k < 4; k++)
                        j[k] = j[k] * det;
                
                return true;
                
            }
        };
        /*-------------------------------------------------------------------*/
        template<int TRow, int TCol, typename TType> TType Matrix<TRow, TCol, TType>::det() const
        {
            return DetTraits<TRow, TCol, TType>::perform(m_rows);
        }
        /*-------------------------------------------------------------------*/
        template<int TRow, int TCol, typename TType>
        TType Matrix<TRow, TCol, TType>::get(const int& i, const int& j) const
        {
            assert(i>=0 && i<TRow && j>=0 && j<TCol);
            return m_rows[i][j];
        }
        /*-------------------------------------------------------------------*/
        template<int TRow, int TCol, typename TType>
        void Matrix<TRow, TCol, TType>::set(const int& i, const int& j, const TType& AVal)
        {
            assert(i>=0 && i<TRow && j>=0 && j<TCol);
            m_rows[i][j] = AVal;
        }
        
        /*-------------------------------------------------------------------*/
        template<int TRow, int TCol, typename TType>
        TType Matrix<TRow, TCol, TType>::
        frobeniusNorm2() const {

			TType frob = 0;
			for(int i=0; i<TRow; i++) {
	                for(int j=0; j<TCol; j++) {
                   	 	frob += m_rows[i][j] * m_rows[i][j];
			}
		}
		return frob;
	}
        /*-------------------------------------------------------------------*/
        template<int TRow, int TCol, typename TType>
        Matrix<TRow, TCol, TType> Matrix<TRow, TCol, TType>::zero()
        {
            TType null_vect[TRow][TCol];
            for (int i = 0; i < TRow; i++)
                for (int j = 0; j < TCol; j++)
                    null_vect[i][j] = 0.0;
            
            return Matrix<TRow, TCol, TType>(null_vect);
        }
        /*-------------------------------------------------------------------*/
        template<int TRow, int TCol, typename TType>
        Matrix<TRow, TCol, TType> Matrix<TRow, TCol, TType>::identity()
        {
	        Matrix<TRow, TCol, TType> res;
	        for (int i=TRow; i--; )
		        for (int j=TCol;j--; res[i][j]=(i==j));
	        return res;
        } 
        /*-------------------------------------------------------------------*/
        template<int TRow, int TCol, typename TType>
        VectorND<TRow, TType>
        Matrix<TRow, TCol, TType>::solve(const VectorND<TRow, TType>& AB) const
        {
            Matrix<TRow, TCol, TType> Ainv = inverse();

            VectorND<TRow, TType> x;
            
            for(int i=0; i<TRow; i++) {
                x[i] = 0.0;
                for(int j=0; j<TRow; j++) {
                    x[i] += Ainv(i,j)*AB[j];
                }
            }
            return x;
        }
        /*-------------------------------------------------------------------*/
        template<int TRow, int TCol, typename TType>
        VectorND<3,TType> Matrix<TRow, TCol, TType>::
        eulerAngles(const int& AR0, const int& AR1, const int& AR2) const {
            if(TRow!=3 || TCol!=3)
                throw GMDSException("Only work for 3x3 rotation matrices");

            VectorND<3,TType> coeff;
            
            const int odd = ((AR0+1)%3 == AR1) ? 0 : 1;
            const int i = AR0;
            const int j = (AR0 + 1 + odd)%3;
            const int k = (AR0 + 2 - odd)%3;
            
            if (AR0==AR2)
            {
                coeff[0] = atan2(m_rows[j][i], m_rows[k][i]);
                TType s2 = Vector2d({m_rows[j][i], m_rows[k][i]}).norm();
                
                if(( odd   && coeff[0]<TType(0)) ||
                   ((!odd) && coeff[0]>TType(0))) {
                    coeff[0] = (coeff[0] > TType(0)) ? coeff[0]-TType(M_PI) : coeff[0]+TType(M_PI);
                    coeff[1] = -atan2(s2, m_rows[i][i]);
                }
                else  {
                    coeff[1] = atan2(s2, m_rows[i][i]);
                }
                
                TType s1 = sin(coeff[0]);
                TType c1 = cos(coeff[0]);
                coeff[2] = atan2(c1*m_rows[j][k] - s1*m_rows[k][k],
                                 c1*m_rows[j][j] - s1*m_rows[k][j]);
            }
            else {
                coeff[0] = atan2(m_rows[j][k], m_rows[k][k]);
                double c2_x = m_rows[i][i];
                double c2_y = m_rows[i][j];
                TType c2bis = Vector3d({c2_x, c2_y, 0}).norm();
                TType c2 = Vector2d({c2_x, c2_y}).norm();

                if(  (odd  && coeff[0]<TType(0)) ||
                   ((!odd) && coeff[0]>TType(0))) {
                    coeff[0] = (coeff[0] > TType(0)) ? coeff[0]-TType(M_PI) : coeff[0]+TType(M_PI);
                    coeff[1] = atan2(-m_rows[i][k], -c2);
                }
                else
                    coeff[1] = atan2(-m_rows[i][k], c2);
                TType s1 = sin(coeff[0]);
                TType c1 = cos(coeff[0]);
                coeff[2] = atan2(s1*m_rows[k][i] - c1*m_rows[j][i],
                                 c1*m_rows[j][j] - s1*m_rows[k][j]);
            }
            if (!odd)
                coeff = -coeff;
            
            return coeff;
        }
        /*-------------------------------------------------------------------*/
        template<int TRow, int TCol, typename TType>
        std::ostream& operator<<(std::ostream& stream, const Matrix<TRow, TCol, TType>& mat)
        {
            for (int i = 0; i < TRow; i++)
            {
                stream << "[ ";
                for (int j = 0; j < TCol; j++)
                    stream << mat(i,j) << " ";
                stream << "] ";
            }
            
            return stream;
        }
        /*-------------------------------------------------------------------*/
        template<int TRow, int TCol, typename TType> Matrix<TRow, TCol, TType>
        operator-(const Matrix<TRow, TCol, TType>& a, const Matrix<TRow, TCol, TType>& b)
        {
            Matrix<TRow, TCol, TType> m;
            for (int i = 0; i < TRow; i++)
                for (int j = 0; j < TCol; j++)
                    m(i,j) = a(i,j) - b(i,j);
            return m;
        }
        /*-------------------------------------------------------------------*/
        template<int TRow, int TCol, typename TType> Matrix<TRow, TCol, TType>
        operator+(const Matrix<TRow, TCol, TType>& a, const Matrix<TRow, TCol, TType>& b)
        {
            Matrix<TRow, TCol, TType> m;
            for (int i = 0; i < TRow; i++)
                for (int j = 0; j < TCol; j++)
                    m(i,j) = a(i,j) + b(i,j);
            return m;
        }
        /*-------------------------------------------------------------------*/
        template<int TRow, int TCommon, int TCol, typename TType>
        Matrix<TRow, TCol, TType>
        operator*(const Matrix<TRow, TCommon, TType>& a,
                  const Matrix<TCommon, TCol, TType>& b)
        {
            Matrix<TRow, TCol, TType> m;
            for (int i = 0; i < TRow; i++)
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
        template<int TRow, int TCol, typename TType>
        VectorND<TRow, TType>
        operator*(const Matrix<TRow, TCol, TType>& a,
                  const VectorND<TCol, TType>& v)
        {
            VectorND<TRow, TType> r;
            for (int i = 0; i < TRow; i++)
                for (int j = 0; j < TCol; j++) {
                    r[i] += a(i, j)*v[j];
                }
            return r;
        }
        /*-------------------------------------------------------------------*/
        template<int TRow, int TCol, typename TType>
        Matrix<TRow, TCol, TType>
        operator*(const Matrix<TRow, TCol, TType>& a, const TType& k)
        {
            Matrix<TRow, TCol, TType> m;
            for (int i = 0; i < TRow; i++)
                for (int j = 0; j < TCol; j++)
                    m(i,j) = k*a(i,j);
            return m;
        }
        
        /*-------------------------------------------------------------------*/
        template<int TRow, int TCol, typename TType>
        Matrix<TRow, TCol, TType>
        operator/(const Matrix<TRow, TCol, TType>& a, const TType& k)
        {
            Matrix<TRow, TCol, TType> m;
            for (int i = 0; i < TRow; i++)
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

        typedef Matrix<4,4,double> Matrix44;
        typedef Matrix<3,3,double> Matrix33;
        typedef Matrix<2,2,double> Matrix22;
        /*-------------------------------------------------------------------*/

    } // namespace math
    /*------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_MATRIX_H_ */
/*----------------------------------------------------------------------------*/
