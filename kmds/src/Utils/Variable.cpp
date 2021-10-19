/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
 * VariableManager.cpp
 *
 *  Created on: 26 juil. 2010
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include "KM/Utils/Variable.h"
/*----------------------------------------------------------------------------*/
// Kokkos headers
/*----------------------------------------------------------------------------*/
// KMDS headers
#include <gmds/math/Vector.h>
//#include <GMDS/Math/VectorND.h>
/*----------------------------------------------------------------------------*/
// STL headers
/*----------------------------------------------------------------------------*/
namespace kmds {
/*----------------------------------------------------------------------------*/

    /*----------------------------------------------------------------------------*/
    VariableItf::export_type VariableItf::getType()
    {
            export_type t;
            if (dynamic_cast<Variable<int>* >(this) != 0)
                    t = var_int;
            else  if (dynamic_cast<Variable<double>* >(this) != 0)
                    t = var_double;
//            else  if (dynamic_cast<Variable<gmds::math::Vector>* >(this) != 0)
//                    t = var_double_vec;
            else  if (dynamic_cast<Variable<gmds::math::Vector3d>* >(this) != 0)
                    t = var_double_vec;
            else
                    t = var_unknown;
            return t;
    }
    /*----------------------------------------------------------------------------*/
//template <>
//void
//Variable<TCellID>::initialize(const TCellID AI)
//{
//        m_data(AI) = NullID;
//}
///*----------------------------------------------------------------------------*/
//template <>
//void
//Variable<TInt16>::initialize(const TCellID AI)
//{
//        m_data(AI) = 0;
//}
///*----------------------------------------------------------------------------*/
//template <>
//void
//Variable<TInt32>::initialize(const TCellID AI)
//{
//        m_data(AI) = 0;
//}
///*----------------------------------------------------------------------------*/
//template <>
//void
//Variable<TInt64>::initialize(const TCellID AI)
//{
//        m_data(AI) = 0;
//}
///*----------------------------------------------------------------------------*/
//template <>
//void
//Variable<TSize>::initialize(const TCellID AI)
//{
//        m_data(AI) = 0;
//}
///*----------------------------------------------------------------------------*/
//template <>
//void
//Variable<TCoord>::initialize(const TCellID AI)
//{
//        m_data(AI) = 0;
//}
///*----------------------------------------------------------------------------*/
//template <>
//void
//Variable<bool>::initialize(const TCellID AI)
//{
//        m_data(AI) = false;
//}
///*----------------------------------------------------------------------------*/
//template <>
//void
//Variable<gmds::math::Point>::initialize(const TCellID AI)
//{
//        m_data(AI) = gmds::math::Point(0.,0.,0.);
//}
///*----------------------------------------------------------------------------*/
//template <>
//void
//Variable<gmds::math::Vector>::initialize(const TCellID AI)
//{
//        m_data(AI) = gmds::math::Vector(0.,0.,0.);
//}
/*----------------------------------------------------------------------------*/
}  // namespace kmds
/*----------------------------------------------------------------------------*/
