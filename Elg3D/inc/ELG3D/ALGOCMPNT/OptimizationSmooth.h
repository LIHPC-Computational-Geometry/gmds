/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    OptimizationSmooth.h
 *  \author  legoff
 *  \date    10/02/2019
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_OPTIMIZATIONSMOOTH_H_
#define ELG3D_OPTIMIZATIONSMOOTH_H_
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
#include "ELG3D/DATACMPNT/FacetedCurveGeomServices.h"
#include "ELG3D/DATACMPNT/FacetedSurfaceGeomServices.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/

    const int optimizationSmooth_NBMAXCOLORS = 20;


    /*------------------------------------------------------------------------*/
    /** \brief  Smooth the mesh
     *
     *  \param[in]  AMesh the mesh
     *
     *  // TODO add the "smart" part
     */
    void optimizationSmooth_2D(const int nbIter,
                               kmds::Mesh* AMesh,
                               const kmds::Connectivity* c_N2F,
                               const kmds::Connectivity* c_N2N,
                               const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                               const kmds::GrowingView<kmds::TCellID>* AFixedNodes);

    /*------------------------------------------------------------------------*/
    /** \brief  Smooth the mesh
     *
     *  \param[in]  AMesh the mesh
     *
     *  // TODO add the "smart" part
     */
    void optimizationSmooth_3D(const int nbIter,
                               kmds::Mesh* AMesh,
                               const kmds::Connectivity* c_N2R,
                               const kmds::Connectivity* c_N2N,
                               const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                               const kmds::GrowingView<kmds::TCellID>* AFixedNodes);

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_OPTIMIZATIONSMOOTH_H_ */

