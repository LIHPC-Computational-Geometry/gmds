//
// Created by rochec on 15/04/2022.
//

#include <gmds/aero/AeroException.h>

using namespace gmds;

AeroException::AeroException(const std::string& what):
  GMDSException(what)
{

}