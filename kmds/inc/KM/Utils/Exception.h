/*----------------------------------------------------------------------------*/
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
