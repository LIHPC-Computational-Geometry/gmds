/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    FracPresEnforcement.h
 *  \author  legoff
 *  \date    01/14/2020
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_FRACPRESENFORCEMENT_H_
#define ELG3D_FRACPRESENFORCEMENT_H_
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <KM/DS/Mesh.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/DATACMPNT/FracPres.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {

/*------------------------------------------------------------------------*/
    /** \brief  Increases the specified materials frac pres so that they become the most present ones.
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Afp_ref the fracpres carried by AMesh
     *  \param[out] Afp_new the newly computed fracpres
     *  \param[in]  AMatNames the names of the materials to maintain, ordered by decreasing priority
     */
    double FracPresEnforcement_maintain_xD(const kmds::GrowingView<kmds::TCellID>* cellIDs,
                                           const elg3d::FracPres* Afp_ref,
                                           elg3d::FracPres* Afp_new,
                                           const std::vector<std::string> AMatNames);

    /*------------------------------------------------------------------------*/
    /** \brief  Fuses the specified materials into one new material.
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Afp_ref the fracpres carried by AMesh
     *  \param[out] Afp_new the newly computed fracpres
     *  \param[in]  AMatNames the names of the materials to fuse
     *  \param[in]  AMatName the name of the resulting new material
     *  \param[in]  AKeepNumbering whether or not we keep or discard the fused material (useful for outputs in order to
     *              preserve the color scheme)
     */
    double FracPresEnforcement_fuse_xD(const kmds::GrowingView<kmds::TCellID>* cellIDs,
                                       const elg3d::FracPres* Afp_ref,
                                       elg3d::FracPres* Afp_new,
                                       const std::vector<std::string> AMatNames,
                                       const std::string AMatName,
                                       const bool AKeepNumbering);

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_FRACPRESENFORCEMENT_H_ */