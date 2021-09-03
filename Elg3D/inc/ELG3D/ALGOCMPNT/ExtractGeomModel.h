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
/** \file    ExtractGeomModel.h
 *  \author  legoff
 *  \date    01/18/2019
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_EXTRACTGEOMMODEL_H_
#define ELG3D_EXTRACTGEOMMODEL_H_
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

#include <gmds/cad/FACManager.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/


    /*------------------------------------------------------------------------*/
    /** \brief  Pillowing
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]
     *
     */
    void extractGeomModel_extract_2D(
            kmds::Mesh* AMesh,
            const elg3d::MaterialAssignment* Ama,
            const kmds::Connectivity* c_E2C,
            gmds::cad::FACManager* AGeomModel);

    void extractGeomModel_extract_3D(
            kmds::Mesh* AMesh,
            const elg3d::MaterialAssignment* Ama,
            const kmds::Connectivity* c_A2C,
            gmds::cad::FACManager* AGeomModel,
            kmds::Variable <std::uintptr_t> *AGeomassoc);

    void extractGeomModel_extract_AandN_2D(
            kmds::Mesh* AMesh_source,
            const elg3d::MaterialAssignment* Ama,
            const kmds::Connectivity* c_A2C,
            kmds::GrowingView<kmds::TCellID>* ASelectionN,
            kmds::GrowingView<kmds::TCellID>* ASelectionA,
            kmds::Mesh* AMesh_interface,
            kmds::Variable<std::pair<int, int> >* AVarMatAssign);

    void extractGeomModel_extract_AandN_3D(
            kmds::Mesh* AMesh_source,
            const elg3d::MaterialAssignment* Ama,
            const kmds::Connectivity* c_A2C,
            kmds::GrowingView<kmds::TCellID>* ASelectionN,
            kmds::GrowingView<kmds::TCellID>* ASelectionA,
            kmds::Mesh* AMesh_interface,
            kmds::Variable<std::pair<int, int> >* AVarMatAssign,
            kmds::Variable<kmds::TCellID > *AVarNew2Old_N);

    void extractGeomModel_keep_EandN_3D(
            kmds::Mesh* AMesh_interface,
            const kmds::Connectivity* c_E2F,
            const kmds::Variable<std::pair<int, int> >* AVarMatAssign,
            kmds::Variable<std::array<char, 100> >* AVarCurveMat,
            kmds::GrowingView<kmds::TCellID>* ASelectionN_vertices);

    void extractGeomModel_vertexMat_N_3D(
            kmds::Mesh* AMesh_interface,
            const kmds::Connectivity* c_N2E,
            const kmds::Variable<std::array<char, 100> >* AVarCurveMat,
            const kmds::GrowingView<kmds::TCellID>* ASelectionN_vertices,
            kmds::Variable<std::array<char, 100> >* AVarVertexMat);

    void extractGeomModel_buildGeomModel_2D(
            kmds::Mesh* AMesh,
            const kmds::Variable<std::pair<int, int> >* AVarMatAssign,
            const kmds::Connectivity* c_N2E,
            gmds::cad::FACManager* AGeomModel);

    void extractGeomModel_buildGeomModel_3D(
            kmds::Mesh* AMesh,
            const kmds::Variable<std::pair<int, int> >* AVarMatAssign,
            const kmds::Variable<std::array<char, 100> >* AVarCurveMat,
            const kmds::Variable<std::array<char, 100> >* AVarVertexMat,
            const kmds::GrowingView<kmds::TCellID>* ASelectionN_vertices,
            gmds::cad::FACManager* AGeomModel,
            kmds::Variable <std::uintptr_t> *AGeomassoc);

    /*------------------------------------------------------------------------*/
    /** \brief  Pillowing
     *
     *  \param[in]  AMesh the mesh
     *
     */
//    void extractGeomModel_extract_3D(
//            kmds::Mesh* AMesh,
//            const elg3d::MaterialAssignment* Ama,
//            const kmds::Connectivity* c_A2C,
//            gmds::cad::FACManager* AGeomModel,
//            kmds::Variable <std::uintptr_t> *AGeomassoc);

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_EXTRACTGEOMMODEL_H_ */

