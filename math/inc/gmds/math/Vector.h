/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
//  Vector.h
//
//  Created by F. Ledoux on 02/11/2015. Inspirated by the generic
//  vectors of Eigen and and Geogram (INRIA)
//
/*------------------------------------------------------------------------*/
#ifndef GMDS_MATH_VECTOR_H_
#	define GMDS_MATH_VECTOR_H_
/*------------------------------------------------------------------------*/
// STL Headers
#	include <cassert>
#	include <string.h>
/*------------------------------------------------------------------------*/
// GMDS Headers
#	include "GMDSMath_export.h"
#	include <gmds/math/Constants.h>
#	include <gmds/math/Point.h>
/*------------------------------------------------------------------------*/
namespace gmds {
/*--------------------------------------------------------------------*/
namespace math {

/*--------------------------------------------------------------------*/
/** \struct CrossNDPolicy
 *  \brief  Template structure providing a tailored process for
 *          computing the cross product between two N-dim vectors
 *
 *  \tparam TDim	vector dimension
 *  \tparam TType	value type
 */
/*-------------------------------------------------------------------*/
template<int TDim, typename TType> struct CrossNDPolicy
{
	static TType X(const TType t1[TDim], const TType t2[TDim])
	{
		throw GMDSException("Not yet implemented");
	}
	static TType Y(const TType t1[TDim], const TType t2[TDim])
	{
		throw GMDSException("Not yet implemented");
	}
	static TType Z(const TType t1[TDim], const TType t2[TDim])
	{
		throw GMDSException("Not yet implemented");
	}
};
/*-------------------------------------------------------------------*/
template<typename TType> struct CrossNDPolicy<3, TType>
{
	static TType X(const TType t1[3], const TType t2[3])
	{
		return t1[1] * t2[2] - t1[2] * t2[1];
	}
	static TType Y(const TType t1[3], const TType t2[3])
	{
		return -t1[0] * t2[2] + t1[2] * t2[0];
	}
	static TType Z(const TType t1[3], const TType t2[3])
	{
		return t1[0] * t2[1] - t1[1] * t2[0];
	}

};

/*----------------------------------------------------------------*/
/**
 * \brief Mathematical vector of dimension \p TDim and where each
 *        component is of type \p TType. Usual mathematical
 *        operations are provided to handle them in a simple manner.
 *
 * \tparam TDim  vector dimension
 * \tparam TType component type
 */
template<int TDim, class TType> struct GMDSMath_API VectorND
{

	/** alias on this vector type */
	typedef VectorND<TDim, TType> vector_type;
	/** alias on the type of components */
	typedef TType value_type;
	/** alias on the vector dimension */
	static const int dimension = TDim;

		/*-----------------------------------------------------------*/
		/** \brief Default constructor (all component are zero)
		 */
		VectorND()=default;

	/*-----------------------------------------------------------*/
	/** \brief Gets read-write access to the vector data
	 *
	 * \return a pointer onto the first item of the vector
	 */
	value_type *data()
	{
		return m_data;
	}

	/*-----------------------------------------------------------*/
	/** \brief Gets read-only access to the vector data
	 *
	 * \return a pointer onto the first item of the vector
	 */
	const value_type *data() const
	{
		return m_data;
	}
	/*-----------------------------------------------------------*/
	/** \brief Gets read-write access to a vector component
	 *
	 * \param[in]  AI component index
	 * \return     a reference onto the \p i th component
	 */
	value_type &operator[](const int AI)
	{
		assert(AI >= 0 && AI < dimension);
		return m_data[AI];
	}

	/*-----------------------------------------------------------*/
	/** \brief Gets read-only access to a vector component
	 *
	 * \param[in]  AI component index
	 * \return     a reference onto the \p i th component
	 */
	const value_type &operator[](const int AI) const
	{
		assert(AI >= 0 && AI < dimension);
		return m_data[AI];
	}

	/*-----------------------------------------------------------*/
	/** \brief Gets a vector component
	 *
	 * \param[in]  AI component index
	 * \return     the value of  the \p i th component
	 */
	value_type get(const int AI) const
	{
		return m_data[AI];
	}

	/*-----------------------------------------------------------*/
	/** \brief Sets a vector component
	 *
	 * \param[in]  AI component index
	 * \param[in]  the value to set
	 */
	void set(const int AI, const value_type AVal)
	{
		m_data[AI] = AVal;
	}

	/*-----------------------------------------------------------*/
	/** \brief Gets the vector dimension
	 * \return the value of \p TDIM
	 */
	int dim() const
	{
		return dimension;
	}

	/*-----------------------------------------------------------*/
	/** \brief Gets the L1 norm of *this
	 */
	value_type normL1() const
	{
		value_type r = value_type(0.);
		for (int i = 0; i < dimension; i++) {
			r += std::fabs(m_data[i]);
		}
		return r;
	}

	/*-----------------------------------------------------------*/
	/** \brief Gets the squared L2 norm of *this
	 */
	value_type norm2() const
	{
		value_type r = value_type(0.);
		for (int i = 0; i < dimension; i++) {
			r += m_data[i] * m_data[i];
		}
		return r;
	}

	/*-----------------------------------------------------------*/
	/** \brief Gets the  L2 norm of *this
	 */
	value_type norm() const
	{
		return std::sqrt(norm2());
	}

	/*-----------------------------------------------------------*/
	/** \brief Normalizes this vector if its norm is not zero
	 */
	vector_type &normalize()
	{
		value_type n = norm();
		if (n > 1e-24) {
			for (int i = 0; i < dimension; i++) {
				m_data[i] /= n;
			}
		}
		return *this;
	}
	vector_type getNormalize() const
	{
		vector_type v(*this);
		v.normalize();
		return v;
	}
	/*-----------------------------------------------------------*/
	/** \brief Adds vector \p AV to *this
	 *
	 * \param[in] AV a same type vector
	 * \return    a reference to *this
	 */
	vector_type &operator+=(const vector_type &AV)
	{
		for (int i = 0; i < dimension; i++) {
			m_data[i] += AV.m_data[i];
		}
		return *this;
	}
	/*-----------------------------------------------------------*/
	/** \brief Substracts vector \p AV to *this
	 *
	 * \param[in] AV a same type vector
	 * \return    a reference to *this
	 */
	vector_type &operator-=(const vector_type &AV)
	{
		for (int i = 0; i < dimension; i++) {
			m_data[i] -= AV.m_data[i];
		}
		return *this;
	}

	/*-----------------------------------------------------------*/
	/** \brief Multiplies all component of *this by the scalar
	 *         \p AS
	 *
	 * \tparam    T  scalar type compatible with AType
	 * \param[in] AS a \p T-type calar
	 * \return    a reference to *this
	 */
	template<class T> vector_type &operator*=(const T AS)
	{
		for (int i = 0; i < dimension; i++) {
			m_data[i] *= value_type(AS);
		}
		return *this;
	}

	/*-----------------------------------------------------------*/
	/** \brief Divides all component of *this by the scalar \p AS
	 *
	 * \tparam    T  scalar type compatible with AType
	 * \param[in] AS a \p T-type calar
	 * \return    a reference to *this
	 */
	template<class T> vector_type &operator/=(const T AS)
	{
		for (int i = 0; i < dimension; i++) {
			m_data[i] /= value_type(AS);
		}
		return *this;
	}

	/*-----------------------------------------------------------*/
	/** \brief Addition of 2 vectors
	 *
	 * \param[in]  AV another vector
	 * \return     Vector \p *this + \p AV
	 */
	vector_type operator+(const vector_type &AV) const
	{
		vector_type r(*this);
		for (int i = 0; i < dimension; i++) {
			r.m_data[i] += AV.m_data[i];
		}
		return r;
	}
	/*-----------------------------------------------------------*/
	/** \brief Difference of 2 vectors
	 *
	 * \param[in]  AV another vector
	 * \return     Vector \p *this - \p AV
	 */
	vector_type operator-(const vector_type &AV) const
	{
		vector_type r(*this);
		for (int i = 0; i < dimension; i++) {
			r.m_data[i] -= AV.m_data[i];
		}
		return r;
	}

	/*-----------------------------------------------------------*/
	/** \brief Get the opposite of this vector
	 *
	 * \return  opposite vector
	 */
	vector_type opp() const
	{
		vector_type r;
		for (int i = 0; i < dimension; i++) {
			r.m_data[i] = -m_data[i];
		}
		return r;
	}

	/*-----------------------------------------------------------*/
	/** \brief Multiplies a vector \p AV by a scalar \p AS (of a
	 *         compatible type \p T)
	 *
	 * \tparam    T scalar type
	 * \param[in] AS a scalar value
	 * \return    Vector \p *this * \p AS
	 */
	template<class T> vector_type operator*(const T AS) const
	{
		vector_type r(*this);
		for (int i = 0; i < dimension; i++) {
			r.m_data[i] *= value_type(AS);
		}
		return r;
	}
	/*-----------------------------------------------------------*/
	/** \brief Divides a vector \p AV by a scalar \p AS (of a
	 *         compatible type \p T)
	 *
	 * \tparam    T scalar type
	 * \param[in] AS a scalar value
	 * \return    Vector \p *this / \p AS
	 */
	template<class T> vector_type operator/(const T AS) const
	{
		vector_type r(*this);
		for (int i = 0; i < dimension; i++) {
			r.m_data[i] /= value_type(AS);
		}
		return r;
	}

	/*-----------------------------------------------------------*/
	/** \brief Provides the opposite vector
	 *
	 * \return Vector -\p *this
	 */
	vector_type operator-() const
	{
		vector_type r;
		for (int i = 0; i < dimension; i++) {
			r.m_data[i] = -m_data[i];
		}
		return r;
	}

	/*-----------------------------------------------------------*/
	/** \brief Strict equality comparison
	 * \return true if equals, false otherwise
	 */
	bool operator==(const vector_type &AV) const
	{
		for (int i = 0; i < dimension; i++) {
			if (AV.m_data[i] != m_data[i]) return false;
		}
		return true;
	}
	bool operator!=(const vector_type &AV) const
	{
		return !(this->operator==(AV));
	}

	/*------------------------------------------------------------------------*/
	/** \brief Return if the vector is zero
	 */
	bool isZero() const
	{
		for (int i = 0; i < dimension; i++) {
			if (m_data[i] != 0) return false;
		}
		return true;
	}
	/**
	 *
	 * @return One vector orthogonal to this
	 */
	vector_type getOneOrtho() const
	{
		vector_type v_ortho;
		if (TDim == 2) {
			v_ortho[0] = m_data[1];
			v_ortho[1] = -m_data[0];
		}
		else if (TDim == 3) {
			vector_type v_ortho;
			v_ortho[0] = m_data[2] - m_data[1];
			v_ortho[1] = m_data[0] - m_data[2];
			v_ortho[2] = m_data[1] - m_data[0];
		}
		else {
			auto first_non_zero_field = -1;

			for (auto i = 0; i < dimension && first_non_zero_field == -1; i++) {
				if (m_data[i] != 0) first_non_zero_field = i;
			}

			vector_type v_ortho;
			auto j = first_non_zero_field;
			auto k = (first_non_zero_field + 1) % dimension;
			v_ortho[j] = m_data[k] / m_data[j];

			v_ortho[k] = -1;
		}

		return v_ortho;
	}
	/*-----------------------------------------------------------*/
	/** \brief Gets the dot product of this vector with \p AV
	 *
	 * \param[in] AV a vector
	 * \return the dot product (\p *this . \p AV)
	 */
	value_type dot(const vector_type &AV) const
	{
		value_type r = 0;
		for (int i = 0; i < dimension; i++) {
			r += m_data[i] * AV.m_data[i];
		}
		return r;
	}

	/*-----------------------------------------------------------*/
	/** \brief Gets the cross product of this vector with \p AV
	 *
	 * \param[in] AV a vector
	 * \return the cross product (\p *this x \p AV)
	 */
	vector_type cross(const vector_type &AV) const
	{
		return VectorND<dimension, value_type>({CrossNDPolicy<dimension, value_type>::X(m_data, AV.m_data),
		CrossNDPolicy<dimension, value_type>::Y(m_data, AV.m_data),
		CrossNDPolicy<dimension, value_type>::Z(m_data, AV.m_data)});
	}
	/*-----------------------------------------------------------*/
	/** \brief Compute the angle between 0 and 360 degrees
	 *
	 * \param[in] AV a vector
	 * \return the angle between \p *this and \p AV
	 */
	value_type angle(const vector_type &AV) const
	{

		vector_type v1 = *this;
		vector_type v2 = AV;
		v1.normalize();
		v2.normalize();

		if (v1.isZero() || v2.isZero()) {
			throw GMDSException("VectorND::angle one of the vectors is zero.");
		}
		value_type dotProduct = v2.dot(v1);

		if (dotProduct > 1.) {
			return 0.;
		}
		else if (dotProduct < -1.) {
			return gmds::math::Constants::PI;
		}
		return acos(dotProduct);
	}
	/*-----------------------------------------------------------*/
	/** \brief Indicate if all the component of *this are 0.0
	 * \return true if zero, false otherwise
	 */
	bool isNull() const
	{
		for (int i = 0; i < TDim; i++) {
			if (m_data[i] != 0.0) return false;
		}
		return true;
	}
	/*------------------------------------------------------------------------*/
	/** \brief operator power
	 */
	inline vector_type operator^(const int AN)
	{
		vector_type res;
		for (int i = 0; i < TDim; i++) {
			res.m_data[i] = pow(m_data[i], AN);
		}
		return res;
	}

	/*-----------------------------------------------------------*/
	/** \brief Setter on X, Y, Z
	 */
	void setX(const TCoord &AX)
	{
		m_data[0] = AX;
	}

	void setY(const TCoord &AY)
	{
		m_data[1] = AY;
	}

	void setZ(const TCoord &AZ)
	{
		m_data[2] = AZ;
	}

	/*------------------------------------------------------------------------*/
	/** \brief update all XYZ components
	 */
	inline void setXYZ(const TCoord AValX, const TCoord AValY, const TCoord AValZ)
	{
		m_data[0] = AValX;
		m_data[1] = AValY;
		m_data[2] = AValZ;
	}

	/*-----------------------------------------------------------*/
	/** \brief Gets read-only access to the X component
	 *
	 * \return X Component
	 */
	const double &X() const
	{
		return m_data[0];
	}

	/*-----------------------------------------------------------*/
	/** \brief Gets read-only access to the Y component
	 *
	 * \return Y Component
	 */
	const double &Y() const
	{
		return m_data[1];
	}

	/*-----------------------------------------------------------*/
	/** \brief Gets read-only access to the Z component
	 *
	 * \return Z Component
	 */
	const double &Z() const
	{
			return m_data[2];
	}

	/*------------------------------------------------------------------------*/
	/** \brief compute the oriented angle from this to AV. This operation
	 *         requires to give a reference vector, that must be orthogonal to
	 *         the plane containing (*this) and AV. By default, we put it
	 *         equals to Oz(0,0,1)
	 *
	 *  \param AV a vector
	 */
	inline TCoord orientedAngle(const VectorND &AV, const VectorND &AOrtho = VectorND({0, 0, 1})) const
	{
		VectorND ref = cross(AOrtho);
		double a = (AV.dot(ref) > 0) ? -angle(AV) : angle(AV);
		return a;
	}

	/*------------------------------------------------------------------------*/
	/** \brief compute the oriented angle from this to AV. This operation
	 *         requires to give a reference vector, that must be orthogonal to
	 *         the plane containing (*this) and AV. By default, we put it
	 *         equals to Oz(0,0,1)
	 *
	 *  \param AV a vector
	 */
	inline TCoord angleIn02PI(const VectorND &AV, const VectorND &AOrtho = VectorND({0, 0, 1})) const
	{
		VectorND ref = cross(AOrtho);
		double oriented_angle = orientedAngle(AV, AOrtho);
		double a = 0;
		if (oriented_angle > 0.0)
			a = oriented_angle;
		else
			a = Constants::PI2 + oriented_angle;

		return a;
	}

	/*------------------------------------------------------------------------*/
	/** \brief provides the max absolute value component index
	 */
	TInt getMaxAbsComponentIndex() const
	{
		TInt index = 0;
		for (int i = 1; i < 3; ++i) {
			if (fabs(m_data[i]) > fabs(m_data[index])) {
				index = i;
			}
		}
		return index;
	}

	/*------------------------------------------------------------------------*/
	/** \brief sum of all the components of the vector
	 */
	TCoord sumComponents() const
	{
		return m_data[0] + m_data[1] + m_data[2];
	}

	/*------------------------------------------------------------------------*/
	/** \brief point formed by the components of the vector. Only 3D.
	 */
	gmds::math::Point getPoint() const
	{
		return gmds::math::Point(m_data[0], m_data[1], m_data[2]);
	}

	/*------------------------------------------------------------------------*/
	/** \brief  indicates if two vector are colinear
	 *
	 *  \param AV a second vector
	 *
	 *  \return true if vector are colinear, no otherwise
	 */
	bool isColinear(const VectorND &AV) const
	{

		VectorND v(*this);
		VectorND v2 = AV;
		v.normalize();
		v2.normalize();

		TCoord result = fabs(v.dot(v2));
		return (result == 1e+0);
	}

	/*------------------------------------------------------------------------*/
	/** \brief  indicates if two vector are orthogonal
	 *
	 *  \param AV a second vector
	 *
	 *  \return GEOM_YES if vector are orthogonal, GEOM_NO if they are
	 *  		not, GEOM_UNDEF otherwise
	 */
	bool isOrthogonal(const VectorND &AV) const
	{
		return (dot(AV) == 0.0);
	}

	TType m_data[TDim] = {0};
};     // class VectorND

/*--------------------------------------------------------------*/
/** \brief Multiplies a scalar by a vector like done in the class
 *         VectorND
 *
 * \tparam TScalar scalar type
 * \tparam TDim    vector dimension
 * \tparam TVec    vector component type
 *
 * \param[in] AS a scalar value of type \p TScalar
 * \param[in] AV a vector of type dimension \p TDim and component
 *                 type \p TVec
 *
 * \return Vector \p AS * \p AV
 */
template<typename TScalar, int TDim, typename TVec> VectorND<TDim, TVec>
operator*(const TScalar AS, const VectorND<TDim, TVec> &AV)
{
	VectorND<TDim, TVec> r;
	for (int i = 0; i < TDim; i++) {
		r[i] = TVec(AS) * AV[i];
	}
	return r;
}
template<int TDim, typename TVec> Point
operator+(const Point &AP, const VectorND<TDim, TVec> &AV)
{
	Point r;
	for (int i = 0; i < TDim; i++) {
		r[i] = AP[i] + AV[i];
	}
	return r;
}
template<int TDim, typename TVec> Point
operator+(const VectorND<TDim, TVec> &AV,const Point &AP)
{
	Point r;
	for (int i = 0; i < TDim; i++) {
		r[i] = AP[i] + AV[i];
	}
	return r;
}

template<int TDim, typename TType> TType
operator*(const VectorND<TDim, TType> &AV1, const VectorND<TDim, TType> &AV2)
{
	TType r = 0;
	for (int i = TDim; i--; r += AV1[i] * AV2[i])
		;
	return r;
}


typedef VectorND<2, TCoord> Vector2d;
typedef VectorND<3, TCoord> Vector3d;
typedef VectorND<4, TCoord> Vector4d;
typedef VectorND<9, TCoord> Vector9d;
typedef VectorND<3, int64_t> Vector3i;

typedef VectorND<3, TCoord> Vector;

Vector3d GMDSMath_API operator-(const Point& AP1, const Point& AP2);
/** Convert a point into a 3D vector
 * @param AP a 3D point
 * @return the corresponding 3D vector
 */
Vector3d GMDSMath_API vec(const Point&AP);
/*---------------------------------------------------------------------*/
/** \brief Writes a vector to a stream
 *
 * \param[in] AStream the output stream
 * \param[in] AV      the vector to write
 * \return a reference to the output stream \p AStream
 */
template<int TDim, class TType> std::ostream &
operator<<(std::ostream &AStream, const VectorND<TDim, TType> &AV)
{
	AStream << "(";
	for (int i = 0; i < TDim - 1; i++) {
		AStream << AV[i] << ", ";
	}
	AStream << AV[TDim - 1] << ")";
	return AStream;
}

}     // namespace math
/*-----------------------------------------------------------------*/
}     // namespace gmds
/*---------------------------------------------------------------------*/
#endif /*GMDS_MATH_VECTOR_ND_H_*/
/*---------------------------------------------------------------------*/
