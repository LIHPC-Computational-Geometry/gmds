/*----------------------------------------------------------------------------*/

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

