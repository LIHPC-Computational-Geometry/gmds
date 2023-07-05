/*----------------------------------------------------------------------------*/
/*
 * Plane.h
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_PLANE_H_
#define GMDS_MATH_PLANE_H_
/*----------------------------------------------------------------------------*/
// GMDS File Headers
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
#include "GMDSMath_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
class Segment;
class Triangle;
/*----------------------------------------------------------------------------*/
/** \class Plane
 *  \brief template class implementing a geometrical plane in space
 */
/*----------------------------------------------------------------------------*/
	class GMDSMath_API Plane
{

public:

    enum IntersectionType{
        NO_INTERSECTION=0,
        SEGMENT_FIRST_PNT,
        SEGMENT_SECOND_PNT,
        SEGMENT_MIDDLE
    } ;
	/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 */
	Plane();

	/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 * \param AP point of the plane
	 * \param AN normal vector
	 */
	Plane(const Point&AP, const Vector3d& AN);

	/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 * \param AP1 a point of the plane
	 * \param AP2 a point of the plane
	 * \param AP3 a point of the plane
	 */
	Plane(const Point&AP1, const Point &AP2, const Point &AP3);



        /*------------------------------------------------------------------------*/
        /** \brief  constructor.
         *
         * \param AT A 3D Triangle
         */
        explicit Plane(const Triangle& AT);

	/*------------------------------------------------------------------------*/
	/** \brief  Overloaded operator=
	 */
	 Plane& operator= (const Plane&);

	/*------------------------------------------------------------------------*/
	/** \brief  Overloaded operator==
	 */
	bool operator==(const Plane&) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Overloaded operator!=
	 */
	bool operator!=(const Plane&) const;

	/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 * \param AP point of the plane
	 * \param AN normal vector
	 */
	void set(const Point &AP, const Vector3d & AN) {
		m_pnt = AP;
		m_normal = AN;
        m_normal.normalize();
    }

	/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 * \param AP1 a point of the plane
	 * \param AP2 a point of the plane
	 * \param AP3 a point of the plane
	 */
	void set(const Point &AP1, const Point &AP2, const Point &AP3) {
		m_pnt = AP1;
		Vector3d v1=AP2-AP1;
		Vector3d v2=AP3-AP1;
		m_normal = v1.cross(v2);
        m_normal.normalize();
    }

	/*------------------------------------------------------------------------*/
	/** \brief  provides the coeff AA, AB, AC and AD for the parametric
	 * 			representation of the plane where AAx,+ABy+ACz=AD is the plane
	 * 			equation
	 *
	 * \param AA first component
	 * \param AB first component
	 * \param AC first component
	 * \param AD first component
	 */
	void getEquationCoeffs(TCoord& AA, TCoord& AB, TCoord& AC, TCoord& AD) const
	{
		AA = m_normal.X();
		AB = m_normal.Y();
		AC = m_normal.Z();
		Vector3d v({m_pnt.X(), m_pnt.Y(), m_pnt.Z()});
		AD = m_normal.dot(v);
	}

	/*------------------------------------------------------------------------*/
	/** \brief  predicate indicating if a point is on a plane
	 *
	 * \param AP  a point
	 */
	bool isIn(const Point& AP) const;

	/*------------------------------------------------------------------------*/
	/** \brief  predicate indicating if this point is on the left of the plane.
	 *
	 * \param AP a point
	 */
	bool isStrictlyOnLeft(const Point& AP) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Compute the distance of a point AP to a plane APlane
	 *			Warning, only meaningful in 3D.
	 *
	 * \param AP 	 a point
	 *
	 * \return the distance between AP and AR
	 */
	TCoord distance(const Point& AP) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Orthogonal projection of a point AP onto a plane.
	 *
	 * \param AP 	 a point
	 *
	 * \return the orthogonal projection of AP onto APlane
	 */
	Point project(const Point& AP) const;

	/*------------------------------------------------------------------------*/
	/** \brief  predicate indicating if a segment and a plane intersect each
	 * 			other.
	 * \param AS a segment
	 */
	bool intersect(const Segment& AS, bool AProper = false) const;

	/**
	 * @brief checf it the current plan intersect segment @AS
	 * @param[in] AS one segment
	 * @param[out] PI the intersection point
	 * @param[out] AW0  weight for barycentric coord of PI along AS for
	 * 					the first point of AS
	 * @param[out] AW1  weight for barycentric coord of PI along AS for
	 * 					the second point of AS
	 * @return an IntersectionType info
	 */
    IntersectionType intersect(const Segment& AS, Point &PI,
    		double& AW0,
    		double& AW1,
    		bool AProper = false) const;

	/*------------------------------------------------------------------------*/
	/** \brief  predicate indicating if two planes intersect each other.
	 *
	 * \param AP a plane
	 * \param AProper indicates if limit cases must be true or false
	 *
	 */
	bool intersect(const Plane& AP, bool AProper = false) const;



	/*------------------------------------------------------------------------*/
	/** \brief  Getter for the plane point
	 *
	 * \return a point
	 */
	const Point& getPoint() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Getter for the plane normal
	 *
	 * \return a vector
	 */
	const Vector3d& getNormal() const;
	private:

		/* point */
		Point m_pnt;

		/* vector normal to the plane */
		Vector3d m_normal;


	};
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_PLANE_H_ */
/*----------------------------------------------------------------------------*/
