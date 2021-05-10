/*----------------------------------------------------------------------------*/
/*
 * Matrix_def.h
 *
 *  Created on: 14 juin 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
/** \class Matrix 2x2
 *  \brief template class implementing a full mathematical matrix
 *
 *  \tparam TLin   number of lines
 *  \tparam TCol   number of columns
 *  \tparam TType type of elements stored in the matrix
 */
/*----------------------------------------------------------------------------*/
template <int TLin, int TCol, typename TType>
class EXPORT_GMDS Matrix {
public:

    /** type of matrix components */
    typedef TType value_type;

	/*------------------------------------------------------------------------*/
	/** \brief  Default constructor.
	 */
	Matrix();

	/*------------------------------------------------------------------------*/
	/** \brief ND constructor.
	 */
	Matrix(value_type param[TLin][TCol]);

	/*------------------------------------------------------------------------*/
	/** \brief  Copy constructor.
	 */
	Matrix(const Matrix<TLin, TCol,TType>&);

	/*------------------------------------------------------------------------*/
	/** \brief  Overloaded operator=
	 */
	virtual Matrix<TLin, TCol,TType>& operator= (const Matrix<TLin, TCol,TType>&);

	/*------------------------------------------------------------------------*/
	/** \brief  Overloaded operator==
	 */
	bool operator== (const Matrix<TLin, TCol,TType>&) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Destructor.
	 */
	virtual ~Matrix();

	/*------------------------------------------------------------------------*/
	/** \brief  Creation of the null matrix
	 */
	static Matrix<TLin, TCol,TType> zero();

	/*------------------------------------------------------------------------*/
	/** \brief  Creation of the identity matrix
	 */
	static Matrix<TLin, TCol,TType> identity();


    /*------------------------------------------------------------------------*/
    /** \brief  Creation of the inverse matrix
     */
    Matrix<TLin, TCol,TType> inverse() const;

    /*------------------------------------------------------------------------*/
    /** \brief Solve (*this)* X = \p AB with X the vector solution returned by
     *         this method. This method computes the inverse of *this. It is
     *         only implemented for 3x3 and 4x4 matrices.
     *
     * \param[in] AB the right-hand side vector, which must have \p TLin
     *            components
     *
     * \return The TLin-dimension solution vector
     */
    math::VectorND<TLin, TType> solve(const math::VectorND<TLin, TType>& AB) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Read-only Access to element (i,j)
	 *
	 *	\param i line number
	 *	\param j column number
	 *
	 *  \return the (i,j) th element of *this
	 *
	 *  \assert if i<0 or i>TLin or j<0 or j>TCol
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
     *  \assert if i<0 or i>TLin or j<0 or j>TCol
     */
	value_type operator()(const int& i, const int& j) const {
	  return  m_data[i][j];
	}

    /*------------------------------------------------------------------------*/
    /** \brief  Read and Write Access to element (i,j)
     *
     *	\param i line number
     *	\param j column number
     *
     *  \return the (i,j) th element of *this
     *
     *  \assert if i<0 or i>TLin or j<0 or j>TCol
     */
	value_type& operator()(const int& i, const int& j) {
	  return  m_data[i][j];
	}

	/*------------------------------------------------------------------------*/
	/** \brief  Return a tabular view of the matrix
	*/
	void getTab(value_type (&t)[TLin][TCol]) const
	{
		for (unsigned int i = 0; i < TLin; i++)
		{
			for (unsigned int j = 0; j < TCol; j++)
				t[i][j] = m_data[i][j];
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
	/** \brief  Computes the transpose matrix
	 */
	Matrix<TCol, TLin,TType> transpose() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Computes the squared forbenius norm of the matrix
         */
	value_type frobeniusNorm2() const;


	/*------------------------------------------------------------------------*/
	/** \brief  set vector
	 *
	 *  \param a tabular of TDim TBase-type element
	 */
	void set(const value_type param[TLin][TCol]);

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


protected:
	value_type m_data[TLin][TCol];
};
/*----------------------------------------------------------------------------*/
