/*----------------------------------------------------------------------------*/
/** \file    Exception.h
 *  \author  F. LEDOUX
 *  \date    09/07/2007
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_EXCEPTION_H_
#define GMDS_EXCEPTION_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <string>
#include <exception>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** \class GMDSException
 *  \brief Class which handles exceptions.
 */
/*----------------------------------------------------------------------------*/
class GMDSException: public std::exception
{

public :

    /*------------------------------------------------------------------------*/
    /** \brief Default constructor.
     *
     *  \param AWhat The exception message
*/
    GMDSException (const std::string& AWhat = "")
    : std::exception(), what_ (AWhat)
    {}

    /*------------------------------------------------------------------------*/
    /** \brief Copy constructor.
     *
     *  \param AExc the original exception
*/
    GMDSException (const GMDSException& AExc)
    : std::exception(AExc), what_ (AExc.what_)
    {}

    /*------------------------------------------------------------------------*/
    /** \brief Overload of the = operator.
     *
     *  \param AExc the exception to be equal
*/
    GMDSException& operator= (const GMDSException& AExc)
    {
        if (&AExc != this)
            what_  = AExc.what_;

        return *this;
    }

    /*------------------------------------------------------------------------*/
    /** \brief Destructor.
*/
    virtual ~GMDSException() throw()
    {}

    /*------------------------------------------------------------------------*/
    /** \brief Returns a C-style character string describing the general cause
	 *				   of the current error.
     *
     *  \return the exception message
     */
    virtual const char* what() const throw() { return what_.c_str(); }

protected :

    /** exception message */
    std::string   what_;
};
/*----------------------------------------------------------------------------*/
/** \class GMDSMathException
 *  \brief Class which handles exceptions due to mathematical troubles.
 */
/*----------------------------------------------------------------------------*/
class GMDSMathException: public GMDSException
{
public :

	GMDSMathException(const std::string& what = ""):GMDSException(what){}
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif  /* GMDS_EXCEPTION_H_ */
/*----------------------------------------------------------------------------*/
