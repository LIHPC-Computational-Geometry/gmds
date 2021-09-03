/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    Cavity.h
 *  \author  legoff
 *  \date    05/01/2019
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_CAVITY_H_
#define ELG3D_CAVITY_H_
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
#include <ELG3D/DATACMPNT/MaterialAssignment.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/


    void cavity_pillow_create_2D(const kmds::GrowingView<kmds::TCellID>* ANodeIDs,
                                 kmds::Mesh* AMesh_origin,
                                 kmds::Mesh* AMesh_cavity,
                                 const kmds::Connectivity* c_N2C,
                                 const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                                 kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation_cavity,
                                 const kmds::Variable<bool>* AVarMarkedNodes,
                                 const kmds::Variable<bool>* AVarMarkedEdges,
                                 const kmds::Variable<bool>* AVarMarkedEdges_withBoundary,
                                 const kmds::Variable<bool>* AVarMarkedCells,
                                 kmds::Variable<bool>* AVarMarkedNodes_cavity,
                                 kmds::Variable<bool>* AVarMarkedEdges_cavity,
                                 kmds::Variable<bool>* AVarIsMarkedEdge_withBoundary_cavity,
                                 kmds::Variable<bool>* AVarMarkedCells_cavity,
                                 kmds::Variable<kmds::TCellID>* AVarCavity2originN,
                                 kmds::Variable<kmds::TCellID>* AVarOrigin2cavityN,
                                 kmds::Variable<kmds::TCellID>* AVarCavity2originC,
                                 kmds::Variable<kmds::TCellID>* AVarOrigin2cavityC,
                                 kmds::GrowingView<kmds::TCellID>* ANodes_boundary_cavity);


    void cavity_pillow_insert_2D(kmds::Mesh* AMesh_origin,
                                 kmds::Mesh* AMesh_cavity,
                                 kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                                 const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation_cavity,
                                 kmds::Variable<kmds::TCellID>* AVarCavity2originN,
                                 const kmds::Variable<kmds::TCellID>* AVarOrigin2cavityN,
                                 const kmds::Variable<kmds::TCellID>* AVarCavity2originC,
                                 const kmds::Variable<kmds::TCellID>* AVarOrigin2cavityC,
                                 const int AImat,
                                 elg3d::MaterialAssignment* Ama);

    int cavity_computeNbNewNodes_xD(kmds::Mesh* AMesh_cavity,
                                    const kmds::Variable<kmds::TCellID>* AVarCavity2originN);

    int cavity_computeNbNewCells_2D(kmds::Mesh* AMesh_cavity,
                                    const kmds::Variable<kmds::TCellID>* AVarCavity2originC);

    int cavity_computeNbNewCells_3D(kmds::Mesh* AMesh_cavity,
                                    const kmds::Variable<kmds::TCellID>* AVarCavity2originC);

    void cavity_pillow_create_3D(const kmds::GrowingView<kmds::TCellID>* ANodeIDs,
                                 kmds::Mesh* AMesh_origin,
                                 kmds::Mesh* AMesh_cavity,
                                 const kmds::Connectivity* c_N2C,
                                 const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                                 kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation_cavity,
                                 const kmds::Variable<bool>* AVarMarkedNodes,
                                 const kmds::Variable<bool>* AVarMarkedEdges,
                                 const kmds::Variable<bool>* AVarMarkedEdges_withBoundary,
                                 const kmds::Variable<bool>* AVarMarkedCells,
                                 kmds::Variable<bool>* AVarMarkedNodes_cavity,
                                 kmds::Variable<bool>* AVarMarkedEdges_cavity,
                                 kmds::Variable<bool>* AVarIsMarkedEdge_withBoundary_cavity,
                                 kmds::Variable<bool>* AVarMarkedCells_cavity,
                                 kmds::Variable<kmds::TCellID>* AVarCavity2originN,
                                 kmds::Variable<kmds::TCellID>* AVarOrigin2cavityN,
                                 kmds::Variable<kmds::TCellID>* AVarCavity2originC,
                                 kmds::Variable<kmds::TCellID>* AVarOrigin2cavityC,
                                 kmds::GrowingView<kmds::TCellID>* ANodes_boundary_cavity);

    void cavity_pillow_insert_3D(kmds::Mesh* AMesh_origin,
                                 kmds::Mesh* AMesh_cavity,
                                 kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                                 const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation_cavity,
                                 kmds::Variable<kmds::TCellID>* AVarCavity2originN,
                                 const kmds::Variable<kmds::TCellID>* AVarOrigin2cavityN,
                                 const kmds::Variable<kmds::TCellID>* AVarCavity2originC,
                                 const kmds::Variable<kmds::TCellID>* AVarOrigin2cavityC,
                                 const int AImat,
                                 elg3d::MaterialAssignment* Ama);

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_CAVITY_H_ */

