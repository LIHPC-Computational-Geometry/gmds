

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    Tools.h
 *  \author  legoff
 *  \date    05/01/2018
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_TOOLS_H_
#define ELG3D_TOOLS_H_
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>

#include <KM/DS/Mesh.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/MaterialGradientComputation.h"
#include "ELG3D/DATACMPNT/FracPres.h"
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {

    /*----------------------------------------------------------------------------*/
    void
    Tools_computeMidPointVar_2D(const kmds::GrowingView<kmds::TCellID>* ASelectionCells,
                                const kmds::Mesh* AMesh,
                                kmds::Variable<gmds::math::Point>* AVarMidpoints);

    /*----------------------------------------------------------------------------*/
    void
    Tools_computeMidPointVar_3D(const kmds::GrowingView<kmds::TCellID>* ASelectionCells,
                                const kmds::Mesh* AMesh,
                                kmds::Variable<gmds::math::Point>* AVarMidpoints);

    /*----------------------------------------------------------------------------*/
    void
    Tools_storeNodePos_xD(const kmds::Mesh* AMesh,
                          kmds::Variable<gmds::math::Point>* AVarNodePos);

    /*----------------------------------------------------------------------------*/
    void
    Tools_loadNodePos_xD(kmds::Mesh* AMesh,
                         const kmds::Variable<gmds::math::Point>* AVarNodePos);

    /*----------------------------------------------------------------------------*/
    kmds::TCoord
    Tools_computeScaledJacobian_2D(const kmds::Mesh* AMesh);

    /*----------------------------------------------------------------------------*/
    kmds::TCoord
    Tools_computeScaledJacobian_3D(const kmds::Mesh* AMesh);

    /*----------------------------------------------------------------------------*/
    kmds::TCoord
    Tools_computeScaledJacobian_2D(kmds::Mesh* AMesh,
                                   kmds::Variable<kmds::TCoord>* AVarQuality);

    /*----------------------------------------------------------------------------*/
    kmds::TCoord
    Tools_computeScaledJacobian_3D(kmds::Mesh* AMesh,
                                   kmds::Variable<kmds::TCoord>* AVarQuality);

    /*----------------------------------------------------------------------------*/
    // returns the maximum distance
    kmds::TCoord
    Tools_computeDistance_xD(const kmds::GrowingView<kmds::TCellID>* ASelectionNodes,
                             kmds::Mesh* AMesh,
                             const kmds::Variable<gmds::math::Point>* AVarTargetPosition,
                             kmds::Variable<kmds::TCoord>* AVarDist);

    /*----------------------------------------------------------------------------*/
    void
    Tools_callGETME_2D(kmds::Mesh* AMesh,
                       const kmds::GrowingView<kmds::TCellID>* ASelectionFixedNodes,
                       const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation);

    /*----------------------------------------------------------------------------*/
    void
    Tools_callGETME_3D(kmds::Mesh* AMesh,
                       const kmds::GrowingView<kmds::TCellID>* ASelectionFixedNodes,
                       const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation);

/*------------------------------------------------------------------------*/
    /** \brief  Writes the mesh
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Afp the fracpres carried by AMesh
     *  \param[in]  Ama the material assignment for the cells of AMesh
     */
    void Tools_write_2D(kmds::Mesh* AMesh,
                        const elg3d::FracPres* Afp,
                        const elg3d::MaterialAssignment *Ama,
                        const kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads,
                        const std::string AFileName);

    /*------------------------------------------------------------------------*/
    /** \brief  Writes the mesh
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Afp the fracpres carried by AMesh
     *  \param[in]  Ama the material assignment for the cells of AMesh
     */
    void Tools_write_2D(kmds::Mesh* AMesh,
                        const elg3d::FracPres* Afp,
                        const kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads,
                        const std::string AFileName);

    /*------------------------------------------------------------------------*/
    /** \brief  Writes the mesh
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Afp the fracpres carried by AMesh
     *  \param[in]  Ama the material assignment for the cells of AMesh
     */
    void Tools_write_3D(kmds::Mesh* AMesh,
                        const elg3d::FracPres* Afp,
                        const elg3d::MaterialAssignment *Ama,
                        const kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads,
                        const std::string AFileName);

    /*------------------------------------------------------------------------*/
    /** \brief  Writes the mesh
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Afp the fracpres carried by AMesh
     *  \param[in]  Ama the material assignment for the cells of AMesh
     */
    void Tools_write_3D(kmds::Mesh* AMesh,
                        const elg3d::FracPres* Afp,
                        const kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads,
                        const std::string AFileName);

    /*------------------------------------------------------------------------*/
    /** \brief  Writes the mesh
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Afp the fracpres carried by AMesh
     *  \param[in]  Ama the material assignment for the cells of AMesh
     */
    void Tools_write_2D(kmds::Mesh* AMesh,
                        const elg3d::FracPres* Afp,
                        const elg3d::MaterialAssignment *Ama,
                        const std::string AFileName);

    /*------------------------------------------------------------------------*/
    /** \brief  Writes the mesh, with only one variable for material assignement
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Ama the material assignment for the cells of AMesh
     */
    void Tools_write_lite_2D(kmds::Mesh* AMesh,
                             const elg3d::MaterialAssignment *Ama,
                             const std::string AFileName);

    /*------------------------------------------------------------------------*/
    /** \brief  Writes the mesh
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Afp the fracpres carried by AMesh
     *  \param[in]  Ama the material assignment for the cells of AMesh
     */
    void Tools_write_3D(kmds::Mesh* AMesh,
                        const elg3d::FracPres* Afp,
                        const elg3d::MaterialAssignment *Ama,
                        const std::string AFileName);

    /*------------------------------------------------------------------------*/
    /** \brief  Writes the mesh, with only one variable for material assignement
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Ama the material assignment for the cells of AMesh
     */
    void Tools_write_lite_3D(kmds::Mesh* AMesh,
                             const elg3d::MaterialAssignment *Ama,
                             const std::string AFileName);

    /*------------------------------------------------------------------------*/
    /** \brief  Writes the interfaces of the mesh
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Ama the material assignment for the cells of AMesh
     */
    void Tools_write_interfaces_3D(const kmds::GrowingView<kmds::TCellID>* AInterfaceFaces,
                                   const kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
                                   kmds::Mesh* AMesh,
                                   const kmds::Connectivity* c_F2C,
                                   const elg3d::MaterialAssignment *Ama,
                                   const std::string AFileName);


    /*------------------------------------------------------------------------*/
    /** \brief  Reads the mesh
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Afp the fracpres carried by AMesh
     *  \param[in]  Ama the material assignment for the cells of AMesh
     */
    void Tools_read_fracpres_vtk_3D(kmds::Mesh *AMesh,
                                    elg3d::FracPres *Afp,
                                    const std::string AFileName,
                                    bool ADeduceVoid);

    /*------------------------------------------------------------------------*/
    /** \brief  Reads the mesh
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Ama the material assignment carried by AMesh
     */
    void Tools_read_ma_vtk_3D(kmds::Mesh *AMesh,
                              elg3d::MaterialAssignment *Ama,
                              const std::string AFileName);

    /*------------------------------------------------------------------------*/
    /** \brief  Reads the mesh
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Afp the fracpres carried by AMesh
     *  \param[in]  Ama the material assignment for the cells of AMesh
     */
    void Tools_read_fracpres_txt_3D(kmds::Mesh *AMesh,
                                    elg3d::FracPres *Afp,
                                    const std::string AFileName,
                                    bool ADeduceVoid);

    /*------------------------------------------------------------------------*/
    /** \brief  Reads the mesh
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Afp the fracpres carried by AMesh
     *  \param[in]  Ama the material assignment for the cells of AMesh
     */
    void Tools_read_fracpres_vtk_2D(kmds::Mesh *AMesh,
                                    elg3d::FracPres *Afp,
                                    const std::string AFileName,
                                    bool ADeduceVoid);

    /*------------------------------------------------------------------------*/
    /** \brief  Reads the mesh
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  Afp the fracpres carried by AMesh
     *  \param[in]  Ama the material assignment for the cells of AMesh
     */
    void Tools_read_fracpres_txt_2D(kmds::Mesh *AMesh,
                                    elg3d::FracPres *Afp,
                                    const std::string AFileName,
                                    bool ADeduceVoid);

    /*------------------------------------------------------------------------*/
    /** \brief  Compute volume fractions from a GMDS mesh with surfaces
     *
     *  \param[in]  AMesh the kmds mesh
     *  \param[in]  AgmdsMesh the gmds mesh
     *  \param[in]  ADeduceVoid
     *  \param[out]  Afp the fracpres carried by AMesh
     */
    void Tools_compute_fracpres_2D(kmds::Mesh* AMesh,
                                   gmds::Mesh* AgmdsMesh,
                                   bool ADeduceVoid,
                                   elg3d::FracPres* Afp);

    /*------------------------------------------------------------------------*/
    /** \brief  Compute volume fractions from a GMDS mesh with volumes
     *
     *  \param[in]  AMesh the kmds mesh
     *  \param[in]  AgmdsMesh the gmds mesh
     *  \param[in]  ADeduceVoid
     *  \param[out]  Afp the fracpres carried by AMesh
     */
    void Tools_compute_fracpres_3D(kmds::Mesh* AMesh,
                                   gmds::Mesh* AgmdsMesh,
                                   bool ADeduceVoid,
                                   elg3d::FracPres* Afp);


    /*------------------------------------------------------------------------*/
    /** \brief  Compute the discrepancy
     *
     */
    void Tools_copy_mesh_2D(kmds::Mesh* AMesh_ref,
                            kmds::Mesh* AMesh_new);


    void Tools_compute_fracpres_2D(kmds::Mesh* AMesh_ref,
                                   kmds::Mesh* AMesh_new,
                                   const elg3d::MaterialAssignment* Ama_new,
                                   elg3d::FracPres* Afp_ref);

    void Tools_compute_fracpres_3D(kmds::Mesh* AMesh_ref,
                                   kmds::Mesh* AMesh_new,
                                   const elg3d::MaterialAssignment* Ama_new,
                                   elg3d::FracPres* Afp_ref);

    double Tools_compute_discrepancy_2D(kmds::Mesh* AMesh,
                                        const elg3d::FracPres* Afp_ref,
                                        const elg3d::FracPres* Afp_new,
                                        kmds::Variable<double>* AVarDiscrepancy);

    double Tools_compute_discrepancy_3D(kmds::Mesh* AMesh,
                                        const elg3d::FracPres* Afp_ref,
                                        const elg3d::FracPres* Afp_new,
                                        kmds::Variable<double>* AVarDiscrepancy);

    void Tools_compute_fracpres_source_and_submesh_2D(const kmds::Mesh* AMesh_source,
                                                      const elg3d::FracPres* Afp_source,
                                                      const kmds::Variable<kmds::TCellID>* AVarOldCells2firstSubCells,
                                                      const kmds::Mesh* AMesh_submesh,
                                                      const elg3d::MaterialAssignment* Ama_submesh,
                                                      kmds::Mesh* AMesh_target,
                                                      elg3d::FracPres* Afp_target);

    void Tools_compute_fracpres_source_and_submesh_3D(const kmds::Mesh* AMesh_source,
                                                      const elg3d::FracPres* Afp_source,
                                                      const kmds::Variable<kmds::TCellID>* AVarOldCells2firstSubCells,
                                                      const kmds::Mesh* AMesh_submesh,
                                                      const elg3d::MaterialAssignment* Ama_submesh,
                                                      kmds::Mesh* AMesh_target,
                                                      elg3d::FracPres* Afp_target);
/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_TOOLS_H_ */