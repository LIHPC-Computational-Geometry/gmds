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

