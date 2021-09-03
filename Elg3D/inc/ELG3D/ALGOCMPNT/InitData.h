/*----------------------------------------------------------------------------*/

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
