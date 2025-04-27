/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    BadPillowDetection.h
 *  \author  legoff
 *  \date    02/01/2019
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_BADPILLOWDETECTION_H_
#define ELG3D_BADPILLOWDETECTION_H_
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

#include <gmds/math/Triangle.h>
#include <gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
#include "ELG3D/DATACMPNT/Parameters.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/


    /*------------------------------------------------------------------------*/
    /** \brief  Detects whether the set of cells at a node has a bad configuration for pillowing
     *             Currently determined as having two interface faces with opposite normals
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  ANodeIDs the id of the node
     *  \param[in]  ACellIDs the ids of the set of cells, locally to the node
     *  \param[in]  Affs the fakefaces of the cells
     *  \param[in]  AffIsOutward whether the fakefaces are outward oriented compared to the cell
     *  \param[in]  AffNormals the normals, computed at the node, of the fakefaces
     *
     *  \return  true if the cells form a bad configuration, false otherwise
     */

    bool BadPillowDetection_IsNodeBad_3D(kmds::Mesh* AMesh,
                                         const kmds::TCellID ANodeID,
                                         const std::vector<kmds::TCellID> ACellIDs,
                                         const std::vector<std::vector<kmds::FakeFace> >& Affs,
                                         const std::vector<std::vector<bool> >& AffIsOutward,
                                         const std::vector<std::vector<gmds::math::Vector> >& AffNormals);

    bool BadPillowDetection_IsNodeBad_glpk_3D(kmds::Mesh* AMesh,
                                              const kmds::TCellID ANodeID,
                                              const std::vector<kmds::TCellID> ACellIDs,
                                              const std::vector<std::vector<kmds::FakeFace> >& Affs,
                                              const std::vector<std::vector<bool> >& AffIsOutward,
                                              const std::vector<std::vector<gmds::math::Vector> >& AffNormals);

    bool BadPillowDetection_IsNodeBad_r3d_3D(kmds::Mesh* AMesh,
                                             const kmds::TCellID ANodeID,
                                             const std::vector<kmds::TCellID> ACellIDs,
                                             const std::vector<std::vector<kmds::FakeFace> >& Affs,
                                             const std::vector<std::vector<bool> >& AffIsOutward,
                                             const std::vector<std::vector<gmds::math::Vector> >& AffNormals,
                                             const std::vector<std::vector<gmds::math::Triangle> >& AffTriangles);

    bool BadPillowDetection_IsNodeBad_ortho_3D(kmds::Mesh* AMesh,
                                               const kmds::TCellID ANodeID,
                                               const std::vector<kmds::TCellID> ACellIDs,
                                               const std::vector<std::vector<kmds::FakeFace> >& Affs,
                                               const std::vector<std::vector<bool> >& AffIsOutward,
                                               const std::vector<std::vector<gmds::math::Vector> >& AffNormals,
                                               const std::vector<std::vector<gmds::math::Triangle> >& AffTriangles,
                                               const std::vector<std::vector<bool> >& ffBoundary);

    /*------------------------------------------------------------------------*/
    /** \brief  Pillowing
     *
     *  \param[in]  AMesh the mesh
     *
     */
    void BadPillowDetection_detect_3D(const kmds::GrowingView<kmds::TCellID>* AInterfaceFaces,
                                      const kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
                                      kmds::Mesh* AMesh,
                                      const kmds::Connectivity* c_N2F,
                                      const kmds::Connectivity* c_F2C,
                                      const elg3d::MaterialAssignment* Ama,
                                      kmds::GrowingView<kmds::TCellID>* ASelectionBadNodes);

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_BADPILLOWDETECTION_H_ */

