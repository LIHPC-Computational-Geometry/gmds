/*----------------------------------------------------------------------------*/
//
// Created by totoro on 2019-01-12.
//
/*----------------------------------------------------------------------------*/
#include <gmds/io/IMeshIOService.h>
#include <gmds/math/Vector.h>
#include <gmds/math/Cross.h>
#include <gmds/math/Cross2D.h>
#include <gmds/math/Quaternion.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
IMeshIOService::EVariableType
IMeshIOService::getType(const gmds::VariableItf *AVar) {

    IMeshIOService::EVariableType t;
    if (dynamic_cast<const Variable<int>* >(AVar) != nullptr)
        t = var_int;
    else  if (dynamic_cast<const Variable<double>* >(AVar) != nullptr)
        t = var_double;
    else  if (dynamic_cast<const Variable<math::Vector2d>* >(AVar) != nullptr)
        t = var_double_vec;
    else  if (dynamic_cast<const Variable<math::Vector3d>* >(AVar) != nullptr)
        t = var_double_vec;
    else  if (dynamic_cast<const Variable<math::Cross>* >(AVar) != nullptr)
        t = var_cross;
    else  if (dynamic_cast<const Variable<math::Cross2D>* >(AVar) != nullptr)
        t = var_cross_2D;
    else  if (dynamic_cast<const Variable<math::Quaternion>* >(AVar) != nullptr)
        t = var_quaternion;
    else
        t = var_unknown;
    return t;
}
/*----------------------------------------------------------------------------*/
