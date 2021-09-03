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
/** \file    MaterialGradientComputation.h
 *  \author  legoff
 *  \date    05/01/2018
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_MATERIALGRADIENTCOMPUTATION_H_
#define ELG3D_MATERIALGRADIENTCOMPUTATION_H_
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
#include <KM/DS/Mesh.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/DATACMPNT/FracPres.h"
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/

    const int MaterialGradientComputation_MAXNBMATPERCELLBALL = 10;

    struct MaterialGradientComputation_midcellGradients
    {
        int m_matindex[MaterialGradientComputation_MAXNBMATPERCELLBALL];
        gmds::math::Vector m_grad[MaterialGradientComputation_MAXNBMATPERCELLBALL];

        MaterialGradientComputation_midcellGradients& operator=(const MaterialGradientComputation_midcellGradients& Amatgrad) {

            for(int i=0; i<MaterialGradientComputation_MAXNBMATPERCELLBALL; i++) {
                this->m_matindex[i] = Amatgrad.m_matindex[i];
                this->m_grad[i] = Amatgrad.m_grad[i];
            }

            return *this;
        }

        MaterialGradientComputation_midcellGradients (const MaterialGradientComputation_midcellGradients& Amatgrad) {
            for(int i=0; i<MaterialGradientComputation_MAXNBMATPERCELLBALL; i++) {
                this->m_matindex[i] = Amatgrad.m_matindex[i];
                this->m_grad[i] = Amatgrad.m_grad[i];
            }
        }

        MaterialGradientComputation_midcellGradients () {
            for(int i=0; i<MaterialGradientComputation_MAXNBMATPERCELLBALL; i++) {
                this->m_matindex[i] = -1;
                this->m_grad[i] = gmds::math::Vector (-HUGE_VALF, -HUGE_VALF, -HUGE_VALF);
            }
        }
    };

        const MaterialGradientComputation_midcellGradients MaterialGradientComputation_MIDCELLGRADIENTS_NULL;


    /*------------------------------------------------------------------------*/
    /** \brief  Computes the gradient at one cell for one mat using a least-square approx
     *
     *  \param[in]  Acid the id of the cell
     *  \param[in]  AMat the material id
     *  \param[in]  Afp the fracpres
     *  \param[in]  c_C2C the cell to cell (by nodes) connectivity
     *  \param[in]  AVarCenter a variable carried by the cells storing the centroids
     *  \param[out]  AGrad the computed position
     *
     *  \return true if the system has a solution, false otherwise
     */
    bool MaterialGradientComputation_computeGrad_leastsquare_onecell_onemat_2D(const kmds::TCellID Acid,
                                                                               const int AMat,
                                                                               const elg3d::FracPres* Afp,
                                                                               const kmds::Connectivity* c_C2C,
                                                                               const kmds::Variable<gmds::math::Point>* AVarCenter,
                                                                               gmds::math::Vector& AGrad);

    /*------------------------------------------------------------------------*/
    /** \brief  Computes the gradient from the fracpres, restricted for each cell to
     *          all the materials assigned to said cell and its neighbors
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Afp the fracpres carried by AMesh
     *  \param[in]  Ama the material assignment for the cells of AMesh
     */
    void MaterialGradientComputation_leastsquare_2D(const kmds::GrowingView<kmds::TCellID>* ACellIDsAccessor,
                                                    const kmds::Connectivity* c_C2C_byN,
                                                    const elg3d::FracPres* Afp,
                                                    const elg3d::MaterialAssignment *Ama,
                                                    const kmds::Variable<gmds::math::Point>* AVarMidpoints,
                                                    kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads);

    /*------------------------------------------------------------------------*/
    /** \brief  Computes the gradient from the fracpres, for all the materials
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Afp the fracpres carried by AMesh
     */
    void MaterialGradientComputation_leastsquare_allmat_2D(const kmds::GrowingView<kmds::TCellID>* ACellIDsAccessor,
                                                           const kmds::Connectivity* c_C2C_byN,
                                                           const elg3d::FracPres* Afp,
                                                           const kmds::Variable<gmds::math::Point>* AVarMidpoints,
                                                           kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads);

    /*------------------------------------------------------------------------*/
    /** \brief  Computes the gradient at one cell for one mat using a least-square approx
     *
     *  \param[in]  Acid the id of the cell
     *  \param[in]  AMat the material id
     *  \param[in]  Afp the fracpres
     *  \param[in]  c_C2C the cell to cell (by nodes) connectivity
     *  \param[in]  AVarCenter a variable carried by the cells storing the centroids
     *  \param[out]  AGrad the computed position
     *
     *  \return true if the system has a solution, false otherwise
     */
    bool MaterialGradientComputation_computeGrad_leastsquare_onecell_onemat_3D(const kmds::TCellID Acid,
                                                                     const int AMat,
                                                                     const elg3d::FracPres* Afp,
                                                                     const kmds::Connectivity* c_C2C,
                                                                     const kmds::Variable<gmds::math::Point>* AVarCenter,
                                                                     gmds::math::Vector& AGrad);

    /*------------------------------------------------------------------------*/
    /** \brief  Corrects the assignment of cells based on manifold, defeaturing, ...
     *
     *  \param[in]  AMesh the mesh
     *  \param[in,out]  Afp the fracpres carried by AMesh
     *  \param[in,out]  Ama the material assignment for the cells of AMesh
     */
    void MaterialGradientComputation_leastsquare_3D(const kmds::GrowingView<kmds::TCellID>* ACellIDsAccessor,
                                                    const kmds::Connectivity* c_C2C_byN,
                                                    const elg3d::FracPres* Afp,
                                                    const elg3d::MaterialAssignment *Ama,
                                                    const kmds::Variable<gmds::math::Point>* AVarMidpoints,
                                                    kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads);

    /*------------------------------------------------------------------------*/
    /** \brief  Computes the gradient from the fracpres, for all the materials
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Afp the fracpres carried by AMesh
     */
    void MaterialGradientComputation_leastsquare_allmat_3D(const kmds::GrowingView<kmds::TCellID>* ACellIDsAccessor,
                                                           const kmds::Connectivity* c_C2C_byN,
                                                           const elg3d::FracPres* Afp,
                                                           const kmds::Variable<gmds::math::Point>* AVarMidpoints,
                                                           kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads);


/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_MATERIALGRADIENTCOMPUTATION_H_ */
