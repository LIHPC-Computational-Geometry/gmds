/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * This software is a computer program whose purpose is to provide a set of
 * functionnalities to represent and handle any type of meshes (2D, 3D,
 * triangles, tetrahedra, quad, hexa, polygons, polyhedra, etc.) and write
 * meshing algorithms. So it gathers many mathematical objects like points,
 * segment, quaternions, etc. and basic algorithms useful to build more evolved
 * ones.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and, more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
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

