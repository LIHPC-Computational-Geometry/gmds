/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    AssignCells.h
 *  \author  legoff
 *  \date    22/11/2017
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_ASSIGNCELLS_H_
#define ELG3D_ASSIGNCELLS_H_
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
#include "ELG3D/DATACMPNT/FracPres.h"
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
#include "ELG3D/ALGOCMPNT/ManifoldDetection.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/

    /*------------------------------------------------------------------------*/
    /** \brief  Assigns each cell to a material
     *
     *  \param[in]  AMesh the mesh
     *  \param[in,out]  Afp the fracpres carried by AMesh
     *  \param[out]  Ama the material assignment for the cells of AMesh
     */
// TODO distributed memory, treat only owned cells
    void assignCells_3D(kmds::Mesh* AMesh, elg3d::FracPres* Afp, elg3d::MaterialAssignment* Ama);

    /*------------------------------------------------------------------------*/
    /** \brief  First assignment of cells to a material, based on the highest fracpres
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Afp the fracpres carried by AMesh
     *  \param[out]  Ama the material assignment for the cells of AMesh
     */
    void assignCellsMajorityCriteria_2D(kmds::Mesh* AMesh, elg3d::FracPres* Afp, elg3d::MaterialAssignment* Ama);

    /*------------------------------------------------------------------------*/
    /** \brief  First assignment of cells to a material, based on the highest fracpres
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Afp the fracpres carried by AMesh
     *  \param[out]  Ama the material assignment for the cells of AMesh
     */
    void assignCellsMajorityCriteria_3D(kmds::Mesh* AMesh, elg3d::FracPres* Afp, elg3d::MaterialAssignment* Ama);

    /*------------------------------------------------------------------------*/
    /** \brief  First assignment of cells to a material, based on the highest fracpres
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Afp the fracpres carried by AMesh
     *  \param[out]  Ama the material assignment for the cells of AMesh
     *  \param[in]  AMatNames the names of the materials to maintain, ordered by decreasing priority
     */
    void assignCellsMajorityCriteria_maintain_2D(const kmds::Mesh* AMesh,
                                                 const elg3d::FracPres* Afp,
                                                 elg3d::MaterialAssignment* Ama,
                                                 const std::vector<std::string> AMatNames);

    /*------------------------------------------------------------------------*/
    /** \brief  First assignment of cells to a material, based on the highest fracpres
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Afp the fracpres carried by AMesh
     *  \param[out]  Ama the material assignment for the cells of AMesh
     *  \param[in]  AMatNames the names of the materials to maintain, ordered by decreasing priority
     */
    void assignCellsMajorityCriteria_maintain_3D(const kmds::Mesh* AMesh,
                                                 const elg3d::FracPres* Afp,
                                                 elg3d::MaterialAssignment* Ama,
                                                 const std::vector<std::string> AMatNames);

    /*------------------------------------------------------------------------*/
    /** \brief  Defeaturing
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  AC_C2C cell to cell connectivity
     *  \param[in,out]  Ama the material assignment for the cells of AMesh
     */
    void assignCells_defeaturing_3D(const kmds::Mesh* AMesh,
                                    const kmds::Connectivity* AC_C2C_byN,
                                    elg3d::MaterialAssignment* Ama);

    /*------------------------------------------------------------------------*/
    /** \brief  Refeaturing
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  AC_C2C cell to cell connectivity
     *  \param[in,out]  Ama the material assignment for the cells of AMesh
     */
    void assignCells_refeaturing_XD(const kmds::GrowingView<kmds::TCellID>* ACellIDs,
                                    const elg3d::FracPres* Afp,
                                    const elg3d::MaterialAssignment* Ama,
                                    const kmds::Connectivity* AC_C2C_byN,
                                    kmds::GrowingView<kmds::TCellID>* ACellsList);

    /*------------------------------------------------------------------------*/
    /** \brief  Corrects the assignment of cells based on manifold, defeaturing, ...
     *
     *  \param[in]  AMesh the mesh
     *  \param[in,out]  Afp the fracpres carried by AMesh
     *  \param[in,out]  Ama the material assignment for the cells of AMesh
     */
    void assignCellsCorrection_3D(kmds::Mesh* AMesh, const kmds::Connectivity* c_N2R, elg3d::FracPres* Afp, elg3d::MaterialAssignment* Ama);

    /*------------------------------------------------------------------------*/
    /** \brief  Corrects the assignment of cells based on manifold, defeaturing, ...
     *
     *  \param[in]  AMesh the mesh
     *  \param[in,out]  Afp the fracpres carried by AMesh
     *  \param[in,out]  Ama the material assignment for the cells of AMesh
     */
    void assignCellsCorrection_2D(kmds::Mesh* AMesh, const kmds::Connectivity* c_N2F, elg3d::FracPres* Afp, elg3d::MaterialAssignment* Ama);

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_ASSIGNCELLS_H_ */
