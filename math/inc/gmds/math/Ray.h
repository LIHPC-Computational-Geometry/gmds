/*----------------------------------------------------------------------------*/
/*
 * Ray.h
 *
 *  Created on: 22 oct. 2014
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_RAY_H_
#define GMDS_MATH_RAY_H_
/*----------------------------------------------------------------------------*/
// Gepeto File Headers
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
class Segment;
class Plane;
class Triangle;
/*----------------------------------------------------------------------------*/
/** \class Ray
 *  \brief template class implementing a geometrical ray in space
 */
/*----------------------------------------------------------------------------*/
	class EXPORT_GMDS Ray
{
public:

	/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 * \param AP origin point of the ray
	 * \param ADir vector
	 */
	Ray(const Point& AP, const Vector3d& ADir)
	:m_pnt(AP),m_dir(ADir) {
	  m_isDirUnit = true;
	  m_dirUnit = m_dir.getNormalize();
	}

	/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 * \param AP1 a point of the plane
	 * \param AP2 a point of the plane
	 */
	Ray(const Point& AP1, const Point& AP2)
	:m_pnt(AP1),m_dir(Vector3d(AP1,AP2)) {
	  m_isDirUnit = true;
	  m_dirUnit = m_dir.getNormalize();
	}

        /*------------------------------------------------------------------------*/
        /** \brief  constructor.
         *
         * \param ARay a ray
         */
        Ray(const Ray& ARay)
        :m_pnt(ARay.m_pnt),m_dir(ARay.m_dir),
	  m_isDirUnit(ARay.m_isDirUnit),m_dirUnit(ARay.m_dirUnit) {
        }
	

	/*------------------------------------------------------------------------*/
	/** \brief  Overloaded operator=
	 */
	virtual Ray& operator= (const Ray&);

	/*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator<< for output
         */
        friend EXPORT_GMDS std::ostream& operator<<(std::ostream&, const Ray&);

	/*------------------------------------------------------------------------*/
	/** \brief  Getter for the plane point
	 *
	 * \return a point
	 */
	const Point& getPoint() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Getter for the directional vector
	 *
	 * \return a vector
	 */
	const Vector3d& getDir() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Getter for the directional unit vector
	 *
	 * \return a vector
	 */
	const Vector3d& getDirUnit() const;
    
    /*------------------------------------------------------------------------*/
    /** \brief 2D Intersection with a segment
     *
     * \param[in]  AS the segment we want to get the intersection with
     * \param[out] AP the intersection point if it exists
     * \param[out] AParam the a parameter such that AP = a AS[0]+ (1-a) AS[1]
     *
     * \return true if (*this) intersects AS
     */
    bool intersect2D(const Segment& AS, Point& AP, double& AParam) const;
    
    bool SecondMetIntersect2D(const Segment& AS, Point& AP, double& AParam, double& tempEpsilon) const;
    
    /*------------------------------------------------------------------------*/
    /** \brief 3D Intersection with a segment
     *
     * \param[in]  AS the segment we want to get the intersection with
     * \param[out] AP the intersection point if it exists
     * \param[out] AParamSeg the a parameter such that AP = a AS[0]+ (1-a) AS[1]
     * \param[out] AParamRay the a parameter along the ray
     *
     * \return true if (*this) intersects AS
     */
    bool intersect3D(const Segment& AS, Point& AP, double& AParamSeg,
                     double& AParamRay) const;

	/*------------------------------------------------------------------------*/
        /** \brief 2D Intersection with another ray
         *
         * \param[IN]  AR the ray we want to get the intersection with
	 * \param[OUT] AP the intersection point if it exists
         *
         * \return true if (*this) intersects AS
         */
        bool intersect2D(const Ray& AR, Point& AP) const;

	/*------------------------------------------------------------------------*/
        /** \brief 3D Intersection with a plane
         *
         * \param[IN]  APlane the plane we want to get the intersection with
         * \param[OUT] AP the intersection point if it exists
         *
         * \return true if (*this) intersects APlane
         */
        bool intersect3D(const Plane& APlane, Point& AP) const;

	/*------------------------------------------------------------------------*/
        /** \brief 3D Intersection with a triangle
         *
         * \param[IN]  ATri the triangle we want to get the intersection with
         * \param[OUT] AP the intersection point if it exists
         *
         * \return true if (*this) intersects ATri
         */
        bool intersect3D(const Triangle& ATri, Point& AP) const;

	/*------------------------------------------------------------------------*/
        /** \brief Compute the orthogonal projection of AP on (*this)
         *
         * \param[IN]  AP the point we want to project onto (*this)
         *
         * \return the projected point
         */
        Point project(const Point& AP) const;

private:

	/* point */
	Point m_pnt;

	/* directional vector  */
	Vector3d m_dir;

	/* unit directional vector */
	mutable bool m_isDirUnit;
	mutable Vector3d m_dirUnit;

};
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_RAY_H_ */
/*----------------------------------------------------------------------------*/
