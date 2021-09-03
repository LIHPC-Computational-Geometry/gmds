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
/*----------------------------------------------------------------------------*/
/** \file    Exception.h
 *  \author  F. LEDOUX
 *  \date    03/10/2017
 */
/*----------------------------------------------------------------------------*/
#ifndef KMDS_EXCEPTION_H_
#define KMDS_EXCEPTION_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <exception>
#include <string>
/*----------------------------------------------------------------------------*/
namespace kmds {
/*----------------------------------------------------------------------------*/
/** \class KException
 *  \brief Class which handles exceptions.
 */
/*----------------------------------------------------------------------------*/
class KException : public std::exception
{
 public:
        /*------------------------------------------------------------------------*/
        /** \brief Default constructor.
         *
         *  \param AWhat The exception message
    */
        KException(const std::string& AWhat = "")
         : std::exception()
         , what_(AWhat)
        {
        }

        /*------------------------------------------------------------------------*/
        /** \brief Copy constructor.
         *
         *  \param AExc the original exception
    */
        KException(const KException& AExc)
         : std::exception(AExc)
         , what_(AExc.what_)
        {
        }

        /*------------------------------------------------------------------------*/
        /** \brief Overload of the = operator.
         *
         *  \param AExc the exception to be equal
    */
        KException&
        operator=(const KException& AExc)
        {
                if (&AExc != this)
                        what_ = AExc.what_;

                return *this;
        }

        /*------------------------------------------------------------------------*/
        /** \brief Destructor.
    */
        virtual ~KException() throw()
        {
        }

        /*------------------------------------------------------------------------*/
        /** \brief Returns a C-style character string describing the general cause
             *				   of the current error.
         *
         *  \return the exception message
         */
        virtual const char*
        what() const throw()
        {
                return what_.c_str();
        }

 protected:
        /** exception message */
        std::string what_;
};

/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* KMDS_EXCEPTION_H_ */
/*----------------------------------------------------------------------------*/
