/*----------------------------------------------------------------------------*/
/*
 * Segment.h
 *
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_SEGMENT_H_
#define GMDS_MATH_SEGMENT_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <vector>
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
#include "GMDSMath_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
class Plane;
/*----------------------------------------------------------------------------*/
/** \class Segment
 *  \brief class implementing a geometrical segment
 */
/*----------------------------------------------------------------------------*/
class GMDSMath_API Segment
{
	private:

	/* points */
	Point pnts_[2];

	/*------------------------------------------------------------------------*/
	/* derived attributes */
	mutable bool isVunit_;
	mutable Vector3d vunit_;
	// TODO dans le cas du segment... vaudrait mieux pas !!??


	/*------------------------------------------------------------------------*/
	/** \brief  reset the "evaluators" stocked are false as center,
	 * sphere including, plane including
	 *
	 */
	void _reset();

	public:

		/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 */
	Segment();

	/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 * \param AP1 the first point of the segment
	 * \param AP2 the second point of the segment
	 */
	Segment(const Point&AP1, const Point&AP2);

	/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 * \param AP1 the first point of the segment
	 * \param AP2 the second point of the segment
	 */
	void set(const Point &AP1, const Point  &AP2);

	/*------------------------------------------------------------------------*/
	/** \brief  Overloaded operator=
	 */
	virtual Segment& operator=(const Segment&);


	/*------------------------------------------------------------------------*/
	/** \brief  Getter for the segment points
	 *
	 * \param AIndex an integer, either 0 ot 1
	 *
	 * \return a point
	 */
	Point getPoint(const int AIndex) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Getter for having the vector corresponding to the segment
	 *
	 * \return a vector
	 */
	Vector3d getDir() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Return the unit vector defining the segment direction
	 *
	 * \return  the vector
	 */
	Vector3d& getUnitVector() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Return the center of the segment.
	 *
	 * \return  the center
	 */
	Point computeCenter() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Return the length of the segment.
	 *
	 * \return  the length value
	 */
	TCoord computeLength2() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Return the length of the segment.
	 *
	 * \return  the length value
	 */
	TCoord computeLength() const;

    /*------------------------------------------------------------------------*/
    /** \brief  Computes the bounding box of the segment.
     *
     * \param AMinXYZ the lower front left coordinates
     * \param AMaxXYZ the upper back right coordinates
     */
    void computeBoundingBox(TCoord AMinXYZ[3], TCoord AMaxXYZ[3]) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Overloaded operator==
	 */
	bool operator==(const Segment&) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Overloaded operator!=
	 */
	bool operator!=(const Segment&) const;

	/*------------------------------------------------------------------------*/
	/** \brief  predicate indicating if a point is on a segment
	 *
	 * \param AP a point
	 */
	bool isIn(const Point& AP, const bool AProper = false) const;
	
	/*------------------------------------------------------------------------*/
	/** \brief  predicate indicating if a point is on a segment (check with +/- EPSILON)
	 *
	 * \param AP a point
	 */
	bool isIn2ndMethod(const Point& AP, const bool AProper = false) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Compute the distance of a point AP to a segment
	 *
	 * \param AP a point
	 *
	 * \return the distance between AP and AS
	 */
	TCoord distance(const Point& AP) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Compute the distance of a segment ASegment to a segment
	 *
	 * \param ASegmùent a segment
	 *
	 * \return the distance between ASegment and this segment
	 */
	TCoord distanceInf(const Segment& ASegment) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Orthogonal projection of a point AP onto a segment
	 *
	 * \param AP a point
	 *
	 * \return the orthogonal projection of AP
	 */
	Point project(const Point& AP) const;


	/*------------------------------------------------------------------------*/
	/** \brief  predicate indicating if a segment and a plane intersect each
	 * 			other.
	 *			With AProper=true, if the segment lies in the plane we get
	 * 			GEOM_UNDEF
	 *
	 * 			Warning this predicate is only meaningful in 3D
	 *
	 * \param AP a plane
	 * \param AProper indicates if limit cases must be GEOM_YES or GEOM_UNDEF
	 *
	 * \return GEOM_YES if they intersect, GEOM_NO if they don't and GEOM_UNDEF
	 * 		   if we cannot conclude
	 */
	bool intersect(const Plane& AP, const bool AProper = false) const;
	bool intersect(const Plane& AP, Point &API, const bool AProper = false) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Orthogonal projection of a point AP onto a segment
         *
         * \param AP a point
         *
         * \return the orthogonal projection of AP
         */
        bool intersect2D(const Segment& ASeg, const bool AProper = false) const;

       
	/*------------------------------------------------------------------------*/
        /** \brief 2D Intersection with a segment
         *
         * \param[IN]  AS the segment we want to get the intersection with
	 * \param[OUT] AP the intersection point if it exists
         *
         * \return true if (*this) intersects AS
         */
        bool intersect2D(const Segment& AS, Point& AP) const;

	   bool SecondMetIntersect2D(const Segment& AS, Point& AP, double& AParam, double& tempEpsilon) const;
	    
    	bool intersect3D(const Segment& AS, Point& AP, double& AParamSeg, double& AParamThis) const;

		bool intersect3D(const Segment &AS, Point &AP, double &AParamSeg, double &AParamThis, double &tempEpsilon) const;
	/*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator<< for output
         */
        friend GMDSMath_API std::ostream& operator<<(std::ostream&, const Segment&);
};
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_SEGMENT_H_ */
/*----------------------------------------------------------------------------*/
