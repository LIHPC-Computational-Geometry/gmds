//
// Created by rochec on 15/04/2022.
//

#ifndef GMDS_AEROEXCEPTION_H
#define GMDS_AEROEXCEPTION_H

/*----------------------------------------------------------------------------*/
// STL File Headers
#include <string>
#include <exception>
#include <gmds/utils/Exception.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** \class AeroException
 *  \brief Class which handles exceptions.
 */
/*----------------------------------------------------------------------------*/
class AeroException: public GMDSException
{
 public :

	AeroException(const std::string& what = "");
};
/*----------------------------------------------------------------------------*/
}

#endif     // GMDS_AEROEXCEPTION_H
