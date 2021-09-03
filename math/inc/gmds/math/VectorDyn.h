/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * The GMDS library is a computer program whose purpose is to provide a set of
 * functionnalities to represent and handle any type of meshes (2D, 3D,
 * triangles, tetrahedra, quad, hexa, polygons, polyhedra, etc.) and write
 * meshing algorithms. So it gathers many mathematical objects like points,
 * segment, quaternions, etc. and basic algorithms useful to build more evolved
 * ones.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability. 
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or 
 * data to be ensured and,  more generally, to use and operate it in the 
 * same conditions as regards security. 
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
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
#include "GMDSMath_export.h"
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
  class GMDSMath_API VectorDyn {
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

