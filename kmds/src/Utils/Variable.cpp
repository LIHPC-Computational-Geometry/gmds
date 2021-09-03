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
