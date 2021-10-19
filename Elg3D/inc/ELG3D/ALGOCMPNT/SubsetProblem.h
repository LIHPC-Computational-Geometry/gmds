/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    SubsetProblem.h
 *  \author  legoff
 *  \date    01/23/2020
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_SUBSETPROBLEM_H_
#define ELG3D_SUBSETPROBLEM_H_
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
    /** \brief  Writes the mesh
     *
     *  \param[in]  ACellIDs the ids of the cells of the subset we want to extract
     *  \param[in]  AMesh the reference mesh
     *  \param[in]  Afp_ref the fracpres carried by AMesh
     *  \param[out]  AMesh the extracted mesh subset
     *  \param[out] Afp_new the fracpres of the extracted mesh subset
     *  \param[in]  AKeepMaterials we keep all the materials, even those no longer appearing in the extracted subset
     */
    void SubsetProblem_extract_3D(const kmds::GrowingView<kmds::TCellID>* ACellIDs,
                                  kmds::Mesh* AMesh_ref,
                                  const elg3d::FracPres* Afp_ref,
                                  kmds::Mesh* AMesh_new,
                                  elg3d::FracPres* Afp_new,
                                  const bool AKeepMaterials);


/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_SUBSETPROBLEM_H_ */