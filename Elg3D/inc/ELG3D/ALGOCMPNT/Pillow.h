/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    Pillow.h
 *  \author  legoff
 *  \date    05/13/2018
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_PILLOW_H_
#define ELG3D_PILLOW_H_
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

    const int pillow_MAXNBMATPERNODE = 8;


    /*------------------------------------------------------------------------*/
    /** \brief  Pillowing
     *
     *  \param[in]  AMesh the mesh
     *
     */
    void pillow_2D(const kmds::GrowingView<kmds::TCellID>* AInterfaceEdges,
                   const kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
                   kmds::Mesh* AMesh,
                   const kmds::Connectivity* c_F2C,
                   const kmds::Connectivity* c_N2C,
                   elg3d::MaterialAssignment* Ama,
                   kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation);

    /*------------------------------------------------------------------------*/
    /** \brief  Pillowing
     *
     *  \param[in]  AMesh the mesh
     *
     */
    void pillow_3D(const kmds::GrowingView<kmds::TCellID>* AInterfaceFaces,
                   const kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
                   kmds::Mesh* AMesh,
                   const kmds::Connectivity* c_F2C,
                   const kmds::Connectivity* c_N2C,
                   elg3d::MaterialAssignment* Ama,
                   kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation);



//    void pillow_markNodes_basedondist_xD(const kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
//                                         const kmds::Variable<double>* AVarNodeDist,
//                                         kmds::Variable<bool>* AVarIsMarkedNode);
//
//    void pillow_markCells_firstround_2D(kmds::Mesh* AMesh,
//                                        const kmds::Variable<bool>* AVarIsMarkedNode,
//                                        const int AImat,
//                                        const elg3d::MaterialAssignment* Ama,
//                                        kmds::Variable<bool>* AVarIsMarkedCell);

    void pillow_execute_2D(
                           const kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
                           kmds::Mesh* AMesh,
                           const kmds::Connectivity* c_A2C,
                           const kmds::Connectivity* c_N2A,
                           const kmds::Connectivity* c_N2C,
                           elg3d::MaterialAssignment* Ama,
                           kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,

                           kmds::Variable<double>* AVarNodeDist);

    void pillow_execute_propagate_2D(const kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
                                     kmds::Mesh* AMesh,
                                     const kmds::Connectivity* c_A2C,
                                     const kmds::Connectivity* c_N2A,
                                     const kmds::Connectivity* c_N2C,
                                     elg3d::MaterialAssignment* Ama,
                                     kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                                     kmds::Variable<double>* AVarNodeDist);


    void pillow_execute_propagate_3D(const kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
                                     kmds::Mesh* AMesh,
                                     const kmds::Connectivity* c_A2C,
                                     const kmds::Connectivity* c_N2A,
                                     const kmds::Connectivity* c_N2C,
                                     elg3d::MaterialAssignment* Ama,
                                     kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                                     kmds::Variable<double>* AVarNodeDist);

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_PILLOW_H_ */

