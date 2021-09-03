/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    QualityLoop.h
 *  \author  legoff
 *  \date    04/18/2019
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_QUALITYLOOP_H_
#define ELG3D_QUALITYLOOP_H_
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
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/



//    /*------------------------------------------------------------------------*/
//    /** \brief  Pillowing
//     *
//     *  \param[in]  AMesh the mesh
//     *
//     */
//    void pillow_2D(const kmds::GrowingView<kmds::TCellID>* AInterfaceEdges,
//                   const kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
//                   kmds::Mesh* AMesh,
//                   const kmds::Connectivity* c_F2C,
//                   const kmds::Connectivity* c_N2C,
//                   elg3d::MaterialAssignment* Ama,
//                   kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation);
//
//    /*------------------------------------------------------------------------*/
//    /** \brief  Pillowing
//     *
//     *  \param[in]  AMesh the mesh
//     *
//     */
//    void pillow_3D(const kmds::GrowingView<kmds::TCellID>* AInterfaceFaces,
//                   const kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
//                   kmds::Mesh* AMesh,
//                   const kmds::Connectivity* c_F2C,
//                   const kmds::Connectivity* c_N2C,
//                   elg3d::MaterialAssignment* Ama,
//                   kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation);
//
//
//
////    void pillow_markNodes_basedondist_xD(const kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
////                                         const kmds::Variable<double>* AVarNodeDist,
////                                         kmds::Variable<bool>* AVarIsMarkedNode);
////
////    void pillow_markCells_firstround_2D(kmds::Mesh* AMesh,
////                                        const kmds::Variable<bool>* AVarIsMarkedNode,
////                                        const int AImat,
////                                        const elg3d::MaterialAssignment* Ama,
////                                        kmds::Variable<bool>* AVarIsMarkedCell);
//
//    void pillow_execute_2D(
//            const kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
//            kmds::Mesh* AMesh,
//            const kmds::Connectivity* c_A2C,
////                           const kmds::Connectivity* c_N2C,
//            const kmds::Connectivity* c_N2A,
//            elg3d::MaterialAssignment* Ama,
//            kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
//
//            kmds::Variable<double>* AVarNodeDist);
//
//    void pillow_execute_3D();

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_QUALITYLOOP_H_ */

