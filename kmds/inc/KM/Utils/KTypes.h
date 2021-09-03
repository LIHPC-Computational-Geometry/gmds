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
 * CommonTypes.h
 *
 *  Created on: 6 janv. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef KMDS_TYPES_H_
#define KMDS_TYPES_H_
/*----------------------------------------------------------------------------*/
// STL File headers
#include <stdlib.h>
#include <climits>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>
/*----------------------------------------------------------------------------*/
#include "KFlags.h"
/*----------------------------------------------------------------------------*/
namespace kmds {
/*----------------------------------------------------------------------------*/
typedef std::int16_t TInt16;
typedef std::int32_t TInt32;
typedef std::int64_t TInt64;
typedef std::size_t TSize;
/*----------------------------------------------------------------------------*/
/** \type type used to define cell ids locally to a part */
/*----------------------------------------------------------------------------*/
typedef std::uint32_t TCellID;

const TCellID NullID = UINT_MAX;
const TCellID InfinityID = NullID - 1;
const TCellID MaxID = InfinityID - 1;
const TSize NullTSize = UINT64_MAX;
/*----------------------------------------------------------------------------*/
// Default size for containers
const TSize DefaultContainerSize = 2048;
/*----------------------------------------------------------------------------*/
/* \type TMark used to mark the different map entities. */
/*----------------------------------------------------------------------------*/
typedef std::int32_t TMark;
/*----------------------------------------------------------------------------*/
/* \type TCoord used for the coordinate of a point in the space. */
/*----------------------------------------------------------------------------*/
typedef double TFloat;
typedef TFloat TCoord;
/*----------------------------------------------------------------------------*/
const TCoord TCoord_Epsilon = 1e-4;  // 1e-11;

/*----------------------------------------------------------------------------*/
/** Types of the different usual cells. This type allows us to template element
 * building and editing */
/*----------------------------------------------------------------------------*/
typedef enum {
        KMDS_UNDEF,
        KMDS_NODE, /* 0D*/
        KMDS_EDGE, /* 1D*/
        KMDS_FACE,
        KMDS_QUAD,
        KMDS_TRIANGLE,
        KMDS_PENTAGON,
        KMDS_POLYGON, /* 2D*/
        KMDS_REGION,
        KMDS_HEX,
        KMDS_TETRA,
        KMDS_PRISM3,
        GMDS_PRISM5, /* 3D*/
        KMDS_PRISM6,
        KMDS_PYRAMID,
        KMDS_POLYHEDRA /* 3D*/
} ECellType;
/*----------------------------------------------------------------------------*/
typedef enum {
        N2N = 1,
        N2E = 1 << 1,
        N2F = 1 << 2,
        N2R = 1 << 3,
        E2N = 1 << 4,
        E2E = 1 << 5,
        E2F = 1 << 6,
        E2R = 1 << 7,
        F2N = 1 << 8,
        F2E = 1 << 9,
        F2F = 1 << 10,
        F2R = 1 << 11,
        R2N = 1 << 12,
        R2E = 1 << 13,
        R2F = 1 << 14,
        R2R = 1 << 15,
        DIM0 = 1 << 16,
        DIM1 = 1 << 17,
        DIM2 = 1 << 18,
        DIM3 = 1 << 19,
        N = 1 << 20,
        E = 1 << 21,
        F = 1 << 22,
        R = 1 << 23,
        F2F_byN = 1 << 24,
        R2R_byN = 1 << 25

} EMeshDefinition;
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* KMDS_TYPES_H_ */
/*----------------------------------------------------------------------------*/
