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
/* Variable.t.h
 *
 *  Created on: 3 aout 2010
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef KMDS_VARIABLE_H_
#define KMDS_VARIABLE_H_
/*----------------------------------------------------------------------------*/
#include "KM/Utils/Variable_def.h"
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// KMDS headers
#include "KM/DS/CellHandle.h"
#include "KM/Utils/Exception.h"
#include "KM/Utils/KTypes.h"
/*----------------------------------------------------------------------------*/
// STL headers
#include <iostream>
#include <vector>
/*----------------------------------------------------------------------------*/
namespace kmds {
/*----------------------------------------------------------------------------*/
template <typename T>
void
Variable<T>::serialize(std::ostream& stream)
{
        // // necessary to use iterators
        // m_data.update();

        // const size_t name_size = m_name.size();
        // const size_t nb_values = m_data.top();

        // const size_t total_size = sizeof(int)                /* total size */
        //                           + sizeof(int)              /* name size*/
        //                           + name_size * sizeof(char) /* name */
        //                           + sizeof(int)              /*nb values*/
        //                           + nb_values * sizeof(T);
        // //							(sizeof(id)	+		/* id per value */
        // //							 sizeof(T));		/* one value*/

        // stream.write((char*) &total_size, sizeof(int));
        // stream.write((char*) &name_size, sizeof(int));
        // stream.write(m_name.c_str(), sizeof(char) * m_name.size());

        // m_data.serialize(stream);
}
/*----------------------------------------------------------------------------*/
template <typename T>
void
Variable<T>::unserialize(std::istream& stream)
{
        //   int total_size = 0;
        // int name_size = 0;

        // stream.read((char*) &total_size, sizeof(int));
        // stream.read((char*) &name_size, sizeof(int));
        // char* n = new char[name_size];
        // stream.read(n, sizeof(char) * name_size);
        // m_name.clear();
        // m_name.assign(n, n + name_size);
        // delete[] n;
        // std::cout << "VAR UNSERIALIZATION: " << name_size << " - " << m_name << std::endl;
        // m_data.unserialize(stream);
}
/*----------------------------------------------------------------------------*/
template <typename T>
void
Variable<T>::compact()
{
        //      m_data.compact();
}

/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* KMDS_VARIABLE_H_ */
/*----------------------------------------------------------------------------*/
