/*----------------------------------------------------------------------------*/
/*
 * Point.h
 *
 *  Created on: 6 févr. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_POINT_H_
#define GMDS_MATH_POINT_H_
/*----------------------------------------------------------------------------*/
// gmds file headers
#include <gmds/utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*----------------------------------------------------------------------------*/
    namespace math{
        /*----------------------------------------------------------------------------*/
        /** \class Point
         *  \brief Defines a 3D point
         */
        /*----------------------------------------------------------------------------*/
        class EXPORT_GMDS Point
        {
            
        public:

            /**
             * \brief Constructor
             * \param AX x coordinnate
             */
            Point(const TCoord& AX=0.0, const TCoord& AY=0.0, const TCoord& AZ=0.0);
            
            virtual ~Point();
            
            inline TCoord X() const {return m_coord[0];}
            inline TCoord Y() const {return m_coord[1];}
            inline TCoord Z() const {return m_coord[2];}
            inline TCoord operator[](const int& AI) const{
                return m_coord[AI];
            }
            inline TCoord& X() {return m_coord[0];}
            inline TCoord& Y() {return m_coord[1];}
            inline TCoord& Z() {return m_coord[2];}
            
            inline TCoord& operator[](const int& AI){
                return m_coord[AI];
            }
            
            inline void setX(const TCoord AVal){m_coord[0] = AVal;}
            inline void setY(const TCoord AVal){m_coord[1] = AVal;}
            inline void setZ(const TCoord AVal){m_coord[2] = AVal;}
            
            inline void setXYZ(const TCoord& AX, const TCoord& AY, const TCoord& AZ)
            {m_coord[0]=AX; m_coord[1]=AY; m_coord[2]=AZ; }
            
            /*------------------------------------------------------------------------*/
            /** \brief  Compute the distance of a point AP to a point
             *
             * \param AP a point
             *
             * \return the distance between AP and this
             */
            TCoord distance(const Point& AP) const;

            /*------------------------------------------------------------------------*/
            /** \brief  Compute the square distance of a point AP to a point
             *
             * \param AP a point
             *
             * \return the distance between AP and this
             */
            TCoord distance2(const Point& AP) const;

            /*------------------------------------------------------------------------*/
            /** \brief  predicate indicating if three points are colinear
             *
             * \param AP2 a second point
             * \param AP3 a third point
             */
            bool areColinear(const Point &AP2, const Point& AP3) const;
		  
		  bool areColinear2ndMethod(const Point& AP2, const Point& AP3) const;

            /*------------------------------------------------------------------------*/
            /** \brief  predicate indicating if four points are coplanar
             *
             * 			Warning this predicate is only meaningful in 3D
             *
             *  \param AP2 a second point
             *  \param AP3 a third point
             *  \param AP4 a fourth point
             */
            bool areCoplanar(const Point& AP2, const Point& AP3, const Point& AP4) const;

            /*------------------------------------------------------------------------*/
            /** \brief  predicate indicating if a point is on the left of the line formed
             *          by the two other points
             *
             *  \param AP1 a point
             *  \param AP2 a second point
             */
            bool isStrictlyOnLeft2D(const Point& AP1, const Point& AP2) const;

            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator==
             */
            bool operator==(const Point&) const;

            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator<
             */
            bool operator<(const Point&) const;

            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator<=
             */
            bool operator<=(const Point&) const;

            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator!=
             */
            bool operator!=(const Point&) const;
            
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator+ to create a new point from 2 points
             */
            friend EXPORT_GMDS Point operator+(const Point&, const Point&);
            
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator- to create a new point from 2 points
             */
            friend EXPORT_GMDS Point operator-(const Point&, const Point&);
            
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator* to create a new point
             */
            friend EXPORT_GMDS Point operator*(const TCoord&, const Point&);
            friend EXPORT_GMDS Point operator*(const Point&, const TCoord&);
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator<< for output
             */
            friend EXPORT_GMDS std::ostream& operator<<(std::ostream&, const Point&);
            
            static void computeBarycentric(const math::Point& AT1,
                                           const math::Point& AT2,
                                           const math::Point& AT3,
                                           const math::Point& AT4,
                                           const math::Point& AP,
                                           TCoord& A1,
                                           TCoord& A2,
                                           TCoord& A3,
                                           TCoord& A4);

            static void computeBarycentric(const math::Point& AT1,
                                           const math::Point& AT2,
                                           const math::Point& AT3,
                                           const math::Point& AP,
                                           TCoord& AX, TCoord& AY, TCoord& AZ);
		  
		  static void computeBarycentric2ndMethod(const math::Point& AT1,
                                           		  const math::Point& AT2,
                                           		  const math::Point& AT3,
                                           		  const math::Point& AP,
                                           		  TCoord& AX, TCoord& AY, TCoord& AZ);


            static void computeBarycentric(const std::vector< math::Point>& AT,
                                           const math::Point& AP,
                                           std::vector<TCoord>& ACoeff);

            static Point massCenter(const std::vector<Point>& AT);

        protected:
            static void computeBarycentric2D(const math::Point& AT1, const math::Point& AT2,
                                             const math::Point& AT3,	const math::Point& AP,
                                             TCoord& AX, TCoord& AY, TCoord& AZ);
		  static void computeBarycentric2D2ndMethod(const math::Point& AT1, const math::Point& AT2,
                                             const math::Point& AT3,	const math::Point& AP,
                                             TCoord& AX, TCoord& AY, TCoord& AZ);
        protected:
            TCoord m_coord[3];
        };
        /*--------------------------------------------------------------------*/
    } // namespace math
    /*------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_POINT_H_ */
/*----------------------------------------------------------------------------*/
