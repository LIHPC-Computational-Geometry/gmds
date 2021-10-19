/*----------------------------------------------------------------------------*/
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
