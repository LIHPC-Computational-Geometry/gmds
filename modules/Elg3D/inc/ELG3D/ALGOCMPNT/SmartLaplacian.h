/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    SmartLaplacian.h
 *  \author  legoff
 *  \date    05/11/2018
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_SMARTLAPLACIAN_H_
#define ELG3D_SMARTLAPLACIAN_H_
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

    const int smartLaplacian_NBMAXCOLORS = 20;


    /*------------------------------------------------------------------------*/
    /** \brief  Smooth the mesh
     *
     *  \param[in]  AMesh the mesh
     *
     *  // TODO add the "smart" part
     */
    void smartLaplacian_2D(const int nbIter,
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
    void smartLaplacian_3D(const int nbIter,
                           kmds::Mesh* AMesh,
                           const kmds::Connectivity* c_N2R,
                           const kmds::Connectivity* c_N2N,
                           const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                           const kmds::GrowingView<kmds::TCellID>* AFixedNodes);


    void smartLaplacian_interface_3D(const int nbIter,
                                     kmds::Mesh* AMesh,
                                     const kmds::Connectivity* c_N2R,
                                     const kmds::Connectivity* c_N2N,
                                     const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                                     const kmds::Variable<std::uintptr_t>* AVarNodeGeomInterface,
                                     const kmds::GrowingView<kmds::TCellID>* AFixedNodes);

    void smartLaplacian_interface_fromF_3D(const int nbIter,
                                           kmds::Mesh* AMesh,
                                           const kmds::Connectivity* c_N2F,
                                           const kmds::Connectivity* c_N2N,
                                           const gmds::cad::FACSurface * ASurf,
                                           const elg3d::FacetedSurfaceGeomServices* AGeomServices);


    void smartLaplacian_interface_fromF_doubleGeom_3D(const int nbIter,
                                                      kmds::Mesh* AMesh,
                                                      const kmds::Connectivity* c_N2F,
                                                      const kmds::Connectivity* c_N2N,
                                                      const kmds::Variable <std::uintptr_t> *AGeomassoc_interface_pixels,
                                                      const kmds::Variable <std::uintptr_t> *AGeomassoc_interface_boundingbox,
                                                      const elg3d::FacetedSurfaceGeomServices* AGeomServicesSurfaces,
                                                      const elg3d::FacetedCurveGeomServices* AGeomServicesCurves);
/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_SMARTLAPLACIAN_H_ */

