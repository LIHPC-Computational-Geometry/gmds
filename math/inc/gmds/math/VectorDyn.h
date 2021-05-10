/*----------------------------------------------------------------------------*//*
 * VectorDyn.h
 *
 *  Created on: 15 juin 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_VECTORDYN_H_
#define GMDS_MATH_VECTORDYN_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <iostream>
#include <vector>
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
  /*----------------------------------------------------------------------------*/
  /** \class VectorDyn
   *  \brief class implementing a mathematical vector whose size can be
   *  		dynamically changed
   *
   */
  /*----------------------------------------------------------------------------*/
  class  EXPORT_GMDS VectorDyn {
  public:
    /*------------------------------------------------------------------------*/
    /** \brief  Default constructor. In this case, the vector has no
     * 			components.
     */
    VectorDyn();

    /*------------------------------------------------------------------------*/
    /** \brief 2D constructor.
     *
     *  \param a1 first coordinate
     */
    VectorDyn(const TCoord& a1, const TCoord& a2);

    /*------------------------------------------------------------------------*/
    /** \brief 3D constructor.
     */
    VectorDyn(const TCoord& a1, const TCoord& a2, const TCoord& a3);

    /*------------------------------------------------------------------------*/
    /** \brief nD constructor, full of zeros.
     */
    VectorDyn(const TInt& n);

    /*------------------------------------------------------------------------*/
    /** \brief  Copy constructor.
     */
    VectorDyn(const VectorDyn& );

    /*------------------------------------------------------------------------*/
    /** \brief  Constructor from specif static numeric vectors
     */
    VectorDyn(const Vector& );
    /*------------------------------------------------------------------------*/
    /** \brief  Overloaded operator =
     */
    VectorDyn& operator=(const VectorDyn& vec);
    /*------------------------------------------------------------------------*/
    /** \brief Destructor
     */
    virtual ~VectorDyn();

    /*------------------------------------------------------------------------*/
    /** \brief give the size of the vector
     */
    TInt size() const;

    /*------------------------------------------------------------------------*/
    /** \brief  Access to i th element of the vector
     *       \return the i th element of *this
     *       \exception if i<0 or i>size()
     */
    /*------------------------------------------------------------------------*/
    inline TCoord const& operator[](const TInt& i) const{
      return tab_[i];
    }
    /*------------------------------------------------------------------------*/
    inline TCoord& operator[](const TInt& i) {
      return tab_[i];
    }

    /*------------------------------------------------------------------------*/
    /** \brief  set a vector from a C tab
     *
     *  \param ATab a tabular of ANb TCoord-type element
     *  \param ANb	the number of elements in ATab
     */
    void set(const TCoord* ATab, const TInt ANb);
    /*------------------------------------------------------------------------*/
    /** \brief normalize
     */
    void normalize();

    /*------------------------------------------------------------------------*/
    /** \brief compute the square norm of the vector
     */
    TCoord norm2() const;

    /*------------------------------------------------------------------------*/
    /** \brief compute the L2 norm of the vector
     */
    TCoord norm() const;
    TCoord normL2() const;

    /*------------------------------------------------------------------------*/
    /** \brief compute the L1 norm of the vector
     */
    TCoord normL1() const;

    /*------------------------------------------------------------------------*/
    /** \brief compute the Linf norm of the vector
     */
    TCoord normLinf() const;

    /*------------------------------------------------------------------------*/
    /** \brief compute the Lp norm of the vector
     *
     *  \param AP the value of p in Lp
     */
    TCoord normLp(const TInt& AP) const;

    /*------------------------------------------------------------------------*/
    /** \brief  Dot product
     *
     *  \param AV a second vector
     *
     *  \return this.AV
     */
    TCoord dot(const VectorDyn& AV) const;

    /*------------------------------------------------------------------------*/
    /** \brief Cross product
     */
    VectorDyn cross(const VectorDyn& v) const;

    /*------------------------------------------------------------------------*/
    /** \brief  indicates if two vector are colinear
     *
     *  \param AV a second vector
     *
     *  \return GEOM_YES if vector are colinear, GEOM_NO if they are
     *  		not, GEOM_UNDEF otherwise
     */
    bool isColinear(const VectorDyn&) const;

    /*------------------------------------------------------------------------*/
    /** \brief  indicates if two vector are orthogonal
     *
     *  \param AV a second vector
     *
     *  \return GEOM_YES if vector are orthogonal, GEOM_NO if they are
     *  		not, GEOM_UNDEF otherwise
     */
    bool isOrthogonal(const VectorDyn&) const;

    /*------------------------------------------------------------------------*/
    /** \brief sum of all the components of the vector
     */
    TCoord sumComponents() const;

    /*------------------------------------------------------------------------*/
    /** \brief operator equal
     */
    bool operator==(const VectorDyn& vec) const;


    /*------------------------------------------------------------------------*/
    /** \brief  Overloaded operator- to get the difference between two vectors
     */
    friend VectorDyn operator-(const VectorDyn&, const VectorDyn&);

    /*------------------------------------------------------------------------*/
    /** \brief  Overloaded operator+ to get the sum oftwo vectors
     */
    friend VectorDyn operator+(const VectorDyn&, const VectorDyn&);

    /*------------------------------------------------------------------------*/
    /** \brief  Overloaded operator* the product of a vector by a scalar value
     */
    friend VectorDyn operator*(const TCoord&, const VectorDyn&);

    /*------------------------------------------------------------------------*/
    /** \brief  Overloaded operator* the product of a vector by a scalar value
     */
    friend VectorDyn operator*(const VectorDyn&, const TCoord&);

    /*------------------------------------------------------------------------*/
    /** \brief  Overloaded operator/ of a vector by a scalar value
     */
    friend VectorDyn operator/(const VectorDyn&, const TCoord&);

    /*------------------------------------------------------------------------*/
    /** \brief operator power
     */
    VectorDyn operator^(const TInt n);

    /*------------------------------------------------------------------------*/
    /** \brief  Overloaded operator<< for output
     */
    friend std::ostream& operator<<(std::ostream&, const VectorDyn&);



    /*------------------------------------------------------------------------*/
    /** \brief operator += for vectors
     */
    VectorDyn& operator+=(const VectorDyn& v);


  private:
    TInt size_;
    std::vector<TCoord> tab_;
  };
  /*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH__VECTORDYN_H_ */
/*----------------------------------------------------------------------------*/

