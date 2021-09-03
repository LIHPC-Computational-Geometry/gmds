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
