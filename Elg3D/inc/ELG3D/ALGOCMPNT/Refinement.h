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
/** \file    Refinement.h
 *  \author  legoff
 *  \date    03/23/2020
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_REFINEMENT_H_
#define ELG3D_REFINEMENT_H_
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
//#include <Kokkos_UnorderedMap.hpp>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <KM/Utils/Graph.h>
#include <KM/DS/Mesh.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
//#include "ELG3D/DATACMPNT/FracPres.h"
//#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/

//    struct ManifoldDetection_options
//    {
//        bool badpillow_active;
//
//        ManifoldDetection_options(){
//            badpillow_active = false;
//        }
//    };
//
//    const ManifoldDetection_options manifolddetection_options_default;


//    // WARNING hard-coded
//    // number of conficts solutions stored per nodes
//    const static int MANIFOLDDETECTION_NB_STORED_SOLUTIONS = 10;
//
//    // WARNING hard-coded
//    // Maximum number of potential solutions
//    // higher would be too much and we exit
//    //const static int MANIFOLDDETECTION_MAX_NB_SOLUTIONS = 100000;
//    const static int MANIFOLDDETECTION_MAX_NB_SOLUTIONS = std::pow(5,8)+1;
//
//
//    struct ManifoldDetection_configStorage
//    {
//        int m_conf[MANIFOLDDETECTION_NB_STORED_SOLUTIONS];
//    };

    /*------------------------------------------------------------------------*/
    /** \brief  Checks that the refinement can be applied (currently, only checks that the mesh is
     * a full quad/hex-mesh)
     *
     *  \param[in]  AMesh the mesh
     */
    bool Refinement_validityCheck_2D(const kmds::Mesh* AMesh);

    bool Refinement_validityCheck_3D(const kmds::Mesh* AMesh);

    /*------------------------------------------------------------------------*/
    /** \brief
     *
     *  \param[in]  AInterfaceNodes the nodes to test
     *  \param[in]  AMesh the mesh
     *  \param[in]  c_N2F the node to face connectivity
     *  \param[in]  Ama the material assignment for the cells of AMesh
     *  \param[out]  ASelection a container of the nonmanifold nodes among those on the interface
     */
    kmds::TCellID  Refinement_getInvalidCells_2D(const Kokkos::View<bool *> markedNodes,
                                                 const kmds::Mesh* AMesh,
                                                 kmds::GrowingView<kmds::TCellID>* ASelection);

    kmds::TCellID  Refinement_getInvalidCells_3D(const Kokkos::View<bool *> markedNodes,
                                                 const kmds::Mesh* AMesh,
                                                 kmds::GrowingView<kmds::TCellID>* ASelection);

    /*------------------------------------------------------------------------*/
    /** \brief
     *
     *  \param[in]  AInterfaceNodes the nodes to test
     *  \param[in]  AMesh the mesh
     *  \param[in]  c_N2F the node to face connectivity
     *  \param[in]  Ama the material assignment for the cells of AMesh
     *  \param[out]  ASelection a container of the nonmanifold nodes among those on the interface
     */
    kmds::TCellID Refinement_avoidBadConfig_2D(Kokkos::View<bool *> markedCells,
                                      Kokkos::View<bool *> markedNodes,
                                      const kmds::Mesh* AMesh,
                                      const kmds::Connectivity* c_C2C_byN);

    kmds::TCellID Refinement_avoidBadConfig_3D(Kokkos::View<bool *> markedCells,
                                               Kokkos::View<bool *> markedNodes,
                                               const kmds::Mesh* AMesh,
                                               const kmds::Connectivity* c_C2C_byN);

    /*------------------------------------------------------------------------*/
    /** \brief
     *
     *  \param[in]  AInterfaceNodes the nodes to test
     *  \param[in]  AMesh the mesh
     *  \param[in]  c_N2F the node to face connectivity
     *  \param[in]  Ama the material assignment for the cells of AMesh
     *  \param[out]  ASelection a container of the nonmanifold nodes among those on the interface
     */
    kmds::TCellID Refinement_createNodesEdges_2D(const Kokkos::View<bool *> markedEdges,
                                                 const Kokkos::View<bool *> markedNodes,
                                                 kmds::Mesh* AMesh,
                                                 const kmds::Connectivity* c_E2C,
                                                 Kokkos::View<kmds::TSize *> ACell2NodesExt,
                                                 Kokkos::View<kmds::TCellID *[8]> ANodesExt);

    kmds::TCellID Refinement_createNodesEdges_3D(const Kokkos::View<bool *> markedEdges,
                                                 const Kokkos::View<bool *> markedNodes,
                                                 kmds::Mesh* AMesh,
                                                 const kmds::Connectivity* c_E2C,
                                                 Kokkos::View<kmds::TSize *> ACell2NodesExt,
                                                 Kokkos::View<kmds::TCellID *[48]> ANodesExt);

    kmds::TCellID Refinement_createNodesFaces_3D(const Kokkos::View<bool *> markedFaces,
                                                 const Kokkos::View<bool *> markedNodes,
                                                 kmds::Mesh* AMesh,
                                                 const kmds::Connectivity* c_F2C,
                                                 Kokkos::View<kmds::TSize *> ACell2NodesExt,
                                                 Kokkos::View<kmds::TCellID *[48]> ANodesExt);

    /*------------------------------------------------------------------------*/
    /** \brief  Detects the non manifold nodes on material interfaces
     *
     *  \param[in]  AInterfaceNodes the nodes to test
     *  \param[in]  AMesh the mesh
     *  \param[in]  c_N2F the node to face connectivity
     *  \param[in]  Ama the material assignment for the cells of AMesh
     *  \param[out]  ASelection a container of the nonmanifold nodes among those on the interface
     */
    void Refinement_refine_2D(const kmds::GrowingView<kmds::TCellID>* ACells2Refine,
                              const bool AIsCellList,
                              kmds::Mesh* AMesh,
                              const kmds::Connectivity* c_C2C_byN,
                              const kmds::Connectivity* c_E2C);

    void Refinement_refine_3D(const kmds::GrowingView<kmds::TCellID>* ACells2Refine,
                              const bool AIsCellList,
                              kmds::Mesh* AMesh,
                              const kmds::Connectivity* c_C2C_byN,
                              const kmds::Connectivity* c_E2C,
                              const kmds::Connectivity* c_F2C);



/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_REFINEMENT_H_ */
