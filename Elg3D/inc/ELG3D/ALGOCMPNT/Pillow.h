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

