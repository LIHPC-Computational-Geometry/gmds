/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    MaterialInterfaces.h
 *  \author  legoff
 *  \date    02/22/2018
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_MATERIALINTERFACES_H_
#define ELG3D_MATERIALINTERFACES_H_
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

    /*------------------------------------------------------------------------*/
    /** \brief  Detects the nodes on interfaces between materials.
     *          Does not include nodes at the boundary of only one material (domain boundary)
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  c_N2C the node to cell connectivity
     *  \param[in]  Ama the material assignment for the cells of AMesh
     *  \param[out]  ASelection a container of the nodes on the interfaces
     */
    void MaterialInterfaces_getNodeOnInterfaces(kmds::Mesh* AMesh,
                                                const kmds::Connectivity* c_N2C,
                                                const elg3d::MaterialAssignment* Ama,
                                                kmds::GrowingView<kmds::TCellID>* ASelection);

    /*------------------------------------------------------------------------*/
    /** \brief  Detects the nodes on the interface of one material
     *          Does not include nodes at the boundary of only one material (domain boundary)
     *
     *  \param[in]  AMat the material
     *  \param[in]  AMesh the mesh
     *  \param[in]  c_N2C the node to cell connectivity
     *  \param[in]  Ama the material assignment for the cells of AMesh
     *  \param[out]  ASelection a container of the nodes on the interfaces
     */
    void MaterialInterfaces_getNodeOnInterfaces(int AMat,
                                                kmds::Mesh* AMesh,
                                                const kmds::Connectivity* c_N2C,
                                                const elg3d::MaterialAssignment* Ama,
                                                kmds::GrowingView<kmds::TCellID>* ASelection);

    /*------------------------------------------------------------------------*/
    /** \brief  Detects the edges on interfaces between materials.
     *          Does not include edges at the boundary of only one material (domain boundary)
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  c_N2F the node to face connectivity
     *  \param[in]  Ama the material assignment for the cells of AMesh
     *  \param[out]  ASelection a container of the nodes on the interfaces
     */
    void MaterialInterfaces_getEdgeOnInterfaces(kmds::Mesh* AMesh,
                                                const kmds::Connectivity* c_A2C,
                                                const elg3d::MaterialAssignment* Ama,
                                                kmds::GrowingView<kmds::TCellID>* ASelection);

    /*------------------------------------------------------------------------*/
    /** \brief  Detects the edges on interfaces between materials.
     *          Includes edges at the boundary of only one material (domain boundary)
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  c_N2F the node to face connectivity
     *  \param[in]  Ama the material assignment for the cells of AMesh
     *  \param[out]  ASelection a container of the nodes on the interfaces
     */
    void MaterialInterfaces_getEdgeOnInterfacesWithBoundary(kmds::Mesh* AMesh,
                                                            const kmds::Connectivity* c_A2C,
                                                            const elg3d::MaterialAssignment* Ama,
                                                            kmds::GrowingView<kmds::TCellID>* ASelection);

    /*------------------------------------------------------------------------*/
    /** \brief  Detects the faces on interfaces between materials.
     *          Does not include faces at the boundary of only one material (domain boundary)
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  c_N2F the node to face connectivity
     *  \param[in]  Ama the material assignment for the cells of AMesh
     *  \param[out]  ASelection a container of the nodes on the interfaces
     */
    void MaterialInterfaces_getFaceOnInterfaces(kmds::Mesh* AMesh,
                                                const kmds::Connectivity* c_F2C,
                                                const elg3d::MaterialAssignment* Ama,
                                                kmds::GrowingView<kmds::TCellID>* ASelection);

    /*------------------------------------------------------------------------*/
    /** \brief  Detects the faces on interfaces between materials.
     *          Includes faces at the boundary of only one material (domain boundary)
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  c_N2F the node to face connectivity
     *  \param[in]  Ama the material assignment for the cells of AMesh
     *  \param[out]  ASelection a container of the nodes on the interfaces
     */
    void MaterialInterfaces_getFaceOnInterfacesWithBoundary(kmds::Mesh* AMesh,
                                                            const kmds::Connectivity* c_F2C,
                                                            const elg3d::MaterialAssignment* Ama,
                                                            kmds::GrowingView<kmds::TCellID>* ASelection);

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_MATERIALINTERFACES_H_ */
