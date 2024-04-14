/*----------------------------------------------------------------------------*/

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
                this->m_grad[i] = gmds::math::Vector3d ({-HUGE_VALF, -HUGE_VALF, -HUGE_VALF});
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
