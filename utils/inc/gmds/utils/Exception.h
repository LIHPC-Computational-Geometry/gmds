/*----------------------------------------------------------------------------*/
/** \file    Exception.h
 *  \author  F. LEDOUX
 *  \date    09/07/2007
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_EXCEPTION_H_
#	define GMDS_EXCEPTION_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <string>
#include <stdexcept>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
/** \class GMDSException
 *  \brief Class which handles exceptions.
 */
/*----------------------------------------------------------------------------*/
class GMDSException : public std::runtime_error
{
 public:
	/*------------------------------------------------------------------------*/
	/** \brief Default constructor.
     *
     *  \param AWhat The exception message
	 */
	explicit GMDSException(const std::string &AWhat = "") : std::runtime_error(AWhat) {}

	/*------------------------------------------------------------------------*/
	/** \brief Copy constructor.
     *
     *  \param AExc the original exception
	 */
	GMDSException(const GMDSException &AExc) : std::runtime_error(AExc.what()) {}
};
/*----------------------------------------------------------------------------*/
/** \class GMDSMathException
 *  \brief Class which handles exceptions due to mathematical troubles.
 */
/*----------------------------------------------------------------------------*/
class GMDSMathException : public GMDSException
{
 public:
	explicit GMDSMathException(const std::string &what = "") : GMDSException(what) {}
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_EXCEPTION_H_ */
/*----------------------------------------------------------------------------*/
