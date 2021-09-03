/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * This software is a computer program whose purpose is to provide a set of
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
 * data to be ensured and, more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    InitData.h
 *  \author  legoff
 *  \date    22/11/2017
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_INITDATA_H_
#define ELG3D_INITDATA_H_
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>

#include <KM/DS/Mesh.h>
/*----------------------------------------------------------------------------*/
// elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/DATACMPNT/FracPres.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/
    void initData_3x3_2D(kmds::Mesh* AMesh, elg3d::FracPres* Afp);
/*----------------------------------------------------------------------------*/
    void initData_2x2_2D_nonmanifold(kmds::Mesh* AMesh, elg3d::FracPres* Afp);
/*----------------------------------------------------------------------------*/
    void initData_3x3x3_3D(kmds::Mesh* AMesh, elg3d::FracPres* Afp);
/*----------------------------------------------------------------------------*/
    void initData_2x1x2_3D_nonmanifold(kmds::Mesh* AMesh, elg3d::FracPres* Afp);
/*----------------------------------------------------------------------------*/
    void initData_I_2D_nonmanifold(kmds::Mesh* AMesh, elg3d::FracPres* Afp, int nb_i);
/*----------------------------------------------------------------------------*/
    void initData_IxJxI_3D_nonmanifold(kmds::Mesh* AMesh, elg3d::FracPres* Afp, int nb_i, int nb_j);
/*----------------------------------------------------------------------------*/
    void initData_rainbow_3D_nonmanifold(kmds::Mesh* AMesh, elg3d::FracPres* Afp);
/*----------------------------------------------------------------------------*/
    void initData_internalFormat(kmds::Mesh* AMesh, elg3d::FracPres* Afp, std::string AFilename);
/*----------------------------------------------------------------------------*/
    void initData_unstructFormat(kmds::Mesh* AMesh, elg3d::FracPres* Afp, std::string AFilename);
/*----------------------------------------------------------------------------*/
    void initData_2x2x2_3D(kmds::Mesh* AMesh, elg3d::FracPres* Afp);
/*----------------------------------------------------------------------------*/
    void initData_2x2x2_badpillow_3D(kmds::Mesh* AMesh, elg3d::FracPres* Afp);
/*----------------------------------------------------------------------------*/
    void initData_exodusReader_2D(gmds::Mesh* AMesh, std::string AFileName);
/*----------------------------------------------------------------------------*/
    void initData_exodusReader_3D(gmds::Mesh* AMesh, std::string AFileName);
/*----------------------------------------------------------------------------*/
    void initData_fromExodus_2D(kmds::Mesh* AMesh, elg3d::FracPres* Afp, std::string AFileName, const double Axyz_min[3], const double Axyz_max[3], int ANi, int ANj, bool ADeduceVoid);
/*----------------------------------------------------------------------------*/
    void initData_fromExodus_3D(kmds::Mesh* AMesh, elg3d::FracPres* Afp, std::string AFileName, const double Axyz_min[3], const double Axyz_max[3], int ANi, int ANj, int ANk, bool ADeduceVoid);
    /*----------------------------------------------------------------------------*/
    void initData_fromExodus_vf_3D(int AVoidID, kmds::Mesh* AMesh, elg3d::FracPres* Afp, std::string AFileName);
/*----------------------------------------------------------------------------*/
    void initData_TEClike_2D(kmds::Mesh* AMesh, elg3d::FracPres* Afp, std::string AFilename);
    void initData_TEClike_3D(kmds::Mesh* AMesh, elg3d::FracPres* Afp, std::string AFilename);
    /*----------------------------------------------------------------------------*/
    void initData_fromTXT_3D(kmds::Mesh* AMesh, elg3d::FracPres* Afp, std::string AFileName, const double Axyz_min[3], const double Axyz_max[3], int ANi, int ANj, int ANk, int ANbMat);
/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_INITDATA_H_ */
