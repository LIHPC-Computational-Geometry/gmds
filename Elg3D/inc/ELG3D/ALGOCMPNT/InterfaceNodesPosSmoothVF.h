/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    InterfaceNodesPosSmoothVF.h
 *  \author  legoff
 *  \date    10/29/2018
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_INTERFACENODESPOSSMOOTHVF_H_
#define ELG3D_INTERFACENODESPOSSMOOTHVF_H_
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
#include <KM/Utils/Graph.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/DATACMPNT/FracPres.h"
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/

    const int InterfaceNodesPosSmoothVF_NBSUB = 3;

    const int InterfaceNodesPosSmoothVF_NBSUB2 = InterfaceNodesPosSmoothVF_NBSUB * InterfaceNodesPosSmoothVF_NBSUB;

    const int InterfaceNodesPosSmoothVF_NBSUB3 = InterfaceNodesPosSmoothVF_NBSUB * InterfaceNodesPosSmoothVF_NBSUB * InterfaceNodesPosSmoothVF_NBSUB;

    const double InterfaceNodesPosSmoothVF_SIMULATEDANNEALING_TEMPERATURE = 0.25;
//    const int InterfaceNodesPosSmoothVF_SIMULATEDANNEALING_NBITER = 10000000;
    const int InterfaceNodesPosSmoothVF_SIMULATEDANNEALING_NBITER = 0;

/*----------------------------------------------------------------------------*/

    /*------------------------------------------------------------------------*/
    /** \brief  Computes the local indices, in an IJ fashion, of a node inside a face.
     *
     *  \param[in]  AF a face
     *  \param[in]  AID the id of a node of the face
     *  \param[out]  AI the i index of AID in AF
     *  \param[out]  AJ the j index of AID in AF
     *
     */
    void
    InterfaceNodesPosSmoothVF_getNodeSubIndices_2D(const kmds::Face AF,
                                                   const kmds::TCellID AID,
                                                   int *AI,
                                                   int *AJ);

    /*------------------------------------------------------------------------*/
    /** \brief  Computes the local indices, in an IJK fashion, of a node inside a region.
     *
     *  \param[in]  AR a region
     *  \param[in]  AID the id of a node of the region
     *  \param[out]  AI the i index of AID in AR
     *  \param[out]  AJ the j index of AID in AR
     *  \param[out]  AK the k index of AID in AR
     *
     */
    void
    InterfaceNodesPosSmoothVF_getNodeSubIndices_3D(const kmds::Region AR,
                                                   const kmds::TCellID AID,
                                                   int *AI,
                                                   int *AJ,
                                                   int *AK);

    /*------------------------------------------------------------------------*/
    /** \brief  Computes the local starting and ending indices, in an (i,j) fashion,
     *          of an edge inside a face. We exclude the end-nodes.
     *
     *  \param[in]  AF a face
     *  \param[in]  AID1 the id of the first node of an edge of AF
     *  \param[in]  AID2 the id of the second node of an edge of AF
     *  \param[out]  AImin the starting i index of the edge in AF
     *  \param[out]  AImax the starting i index of the edge in AF
     *  \param[out]  AJmin the starting j index of the edge in AF
     *  \param[out]  AJmax the starting j index of the edge in AF
     *  \param[out]  ADirect whether AID1->AID2 goes in increasing order for i or j
     *
     */
    void
    InterfaceNodesPosSmoothVF_getEdgeSubIndices_2D(const kmds::Face AF,
                                                   const kmds::TCellID AID1,
                                                   const kmds::TCellID AID2,
                                                   int *AImin,
                                                   int *AImax,
                                                   int *AJmin,
                                                   int *AJmax,
                                                   bool *ADirect);

    /*------------------------------------------------------------------------*/
    /** \brief  Computes the local starting and ending indices, in an (i,j,k) fashion,
     *          of an edge inside a region. We exclude the end-nodes.
     *
     *  \param[in]  AR a region
     *  \param[in]  AID1 the id of the first node of an edge of AR
     *  \param[in]  AID2 the id of the second node of an edge of AR
     *  \param[out]  AImin the starting i index of the edge in AR
     *  \param[out]  AImax the starting i index of the edge in AR
     *  \param[out]  AJmin the starting j index of the edge in AR
     *  \param[out]  AJmax the starting j index of the edge in AR
     *  \param[out]  AKmin the starting k index of the edge in AR
     *  \param[out]  AKmax the starting k index of the edge in AR
     *  \param[out]  ADirect whether AID1->AID2 goes in increasing order for i, j or k
     *
     */
    void
    InterfaceNodesPosSmoothVF_getEdgeSubIndices_3D(const kmds::Region AR,
                                                   const kmds::TCellID AID1,
                                                   const kmds::TCellID AID2,
                                                   int *AImin,
                                                   int *AImax,
                                                   int *AJmin,
                                                   int *AJmax,
                                                   int *AKmin,
                                                   int *AKmax,
                                                   bool *ADirect);

    /*------------------------------------------------------------------------*/
    /** \brief  Computes the local starting and ending indices, in an (i,j,k) fashion,
     *          of an face inside a region. We exclude the end-nodes.
     *
     *  \param[in]  AF a face
     *  \param[in]  AID1 the id of the first node of an edge of AF
     *  \param[in]  AID2 the id of the second node of an edge of AF
     *  \param[out]  AImin the starting i index of the edge in AF
     *  \param[out]  AImax the starting i index of the edge in AF
     *  \param[out]  AJmin the starting j index of the edge in AF
     *  \param[out]  AJmax the starting j index of the edge in AF
     *  \param[out]  AKmin the starting k index of the edge in AF
     *  \param[out]  AKmax the starting k index of the edge in AF
     *  \param[out]  ADirect whether AID1->AID2 goes in increasing order for i, j or k
     *
     */
    void
    InterfaceNodesPosSmoothVF_getFaceSubIndices_3D(const kmds::Region AR,
                                                   const kmds::TCellID AID1,
                                                   const kmds::TCellID AID2,
                                                   const kmds::TCellID AID3,
                                                   const kmds::TCellID AID4,
                                                   int *AImin,
                                                   int *AImax,
                                                   int *AJmin,
                                                   int *AJmax,
                                                   int *AKmin,
                                                   int *AKmax,
                                                   int *AFaceRel);

    /*------------------------------------------------------------------------*/
    /** \brief
     *
     */
    void
    InterfaceNodesPosSmoothVF_buildSubMesh_2D(kmds::Mesh* AMesh_source,
                                              kmds::Mesh* AMesh_target,
                                              const kmds::Connectivity* c_N2F,
                                              const kmds::Connectivity* c_E2F,
                                              const FracPres* Afp_source,
                                              MaterialAssignment* Ama_target,
                                              kmds::Variable<kmds::TCellID>* AVarSubCells2oldCells,
                                              kmds::Variable<kmds::TCellID>* AVarOldCells2firstSubCells,
                                              kmds::Variable<bool> *AVarMixedCells_source,
                                              kmds::Variable<bool> *AVarMixedCells_target,
                                              kmds::Variable<double> *AVarSurfvol_source,
                                              kmds::Variable<double> *AVarSurfvol_target);

    void
    InterfaceNodesPosSmoothVF_buildSubMesh_3D(kmds::Mesh* AMesh_source,
                                              kmds::Mesh* AMesh_target,
                                              const kmds::Connectivity* c_N2C,
                                              const kmds::Connectivity* c_E2C,
                                              const kmds::Connectivity* c_F2C,
                                              const FracPres* Afp_source,
                                              MaterialAssignment* Ama_target,
                                              kmds::Variable<kmds::TCellID>* AVarSubCells2oldCells,
                                              kmds::Variable<kmds::TCellID>* AVarOldCells2firstSubCells,
                                              kmds::Variable<bool> *AVarMixedCells_source,
                                              kmds::Variable<bool> *AVarMixedCells_target,
                                              kmds::Variable<double> *AVarSurfvol_source,
                                              kmds::Variable<double> *AVarSurfvol_target);


    void
    InterfaceNodesPosSmoothVF_selectCells_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                                             const kmds::Connectivity *c_C2C_byN,
                                             const kmds::Variable<bool> *AVarMixedCells_target,
                                             kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                             kmds::GrowingView<kmds::TCellID> *ASelectionCellsIDs_target);

    void
    InterfaceNodesPosSmoothVF_buildGraph_xD(const kmds::GrowingView<kmds::TCellID> *ASelectionCellsIDs_target,
                                            const kmds::Connectivity *c_F2F_byN,
                                            const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                            kmds::Graph* AGraph);

    void
    InterfaceNodesPosSmoothVF_assignHeuristic_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
                                                 const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                                                 const FracPres *Afp_source,
                                                 MaterialAssignment *Ama_target,
                                                 const kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                                                 const kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                                                 const kmds::Variable<bool> *AVarMixedCells_source,
                                                 const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                                 const kmds::Variable<double> *AVarSurfvol_source,
                                                 const kmds::Variable<double> *AVarSurfvol_target,
                                                 const kmds::Graph *AGraph,
                                                 const int ANbSubPixels);

    void
    InterfaceNodesPosSmoothVF_assignHeuristic_withVolumeCheck_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
                                                                 const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                                                                 const FracPres *Afp_source,
                                                                 MaterialAssignment *Ama_target,
                                                                 const kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                                                                 const kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                                                                 const kmds::Variable<bool> *AVarMixedCells_source,
                                                                 const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                                                 const kmds::Variable<double> *AVarSurfvol_source,
                                                                 const kmds::Variable<double> *AVarSurfvol_target,
                                                                 const kmds::Graph *AGraph,
                                                                 const int ANbSubPixels);

    void
    InterfaceNodesPosSmoothVF_assignHeuristic_withVolumeCheck_horizontal_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
                                                                            const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                                                                            const FracPres *Afp_source,
                                                                            MaterialAssignment *Ama_target,
                                                                            const kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                                                                            const kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                                                                            const kmds::Variable<bool> *AVarMixedCells_source,
                                                                            const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                                                            const kmds::Variable<double> *AVarSurfvol_source,
                                                                            const kmds::Variable<double> *AVarSurfvol_target,
                                                                            const kmds::Graph *AGraph,
                                                                            const kmds::Variable<gmds::math::Point> *AVarMidPoints_target,
                                                                            const int ANbSubPixels);

    void
    InterfaceNodesPosSmoothVF_assignHeuristic_print_xD(kmds::Mesh* AMesh_target,
                                                       const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
                                                       const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                                                       const FracPres *Afp_source,
                                                       MaterialAssignment *Ama_target,
                                                       const kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                                                       const kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                                                       const kmds::Variable<bool> *AVarMixedCells_source,
                                                       const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                                       const kmds::Variable<double> *AVarSurfvol_source,
                                                       const kmds::Variable<double> *AVarSurfvol_target,
                                                       const kmds::Graph *AGraph,
                                                       const int ANbSubPixels);

    void
    InterfaceNodesPosSmoothVF_assignHeuristic_withVolumeCheck_print_xD(kmds::Mesh* AMesh_target,
                                                                       const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
                                                                       const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                                                                       const FracPres *Afp_source,
                                                                       MaterialAssignment *Ama_target,
                                                                       const kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                                                                       const kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                                                                       const kmds::Variable<bool> *AVarMixedCells_source,
                                                                       const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                                                       const kmds::Variable<double> *AVarSurfvol_source,
                                                                       const kmds::Variable<double> *AVarSurfvol_target,
                                                                       const kmds::Graph *AGraph,
                                                                       const int ANbSubPixels);

    void
    InterfaceNodesPosSmoothVF_assignMIP_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
                                           const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                                           const FracPres *Afp_source,
                                           MaterialAssignment *Ama_target,
                                           const kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                                           const kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                                           const kmds::Variable<bool> *AVarMixedCells_source,
                                           const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                           const kmds::Graph *AGraph,
                                           const int ANbSubPixels);


    void
    InterfaceNodesPosSmoothVF_assignLP_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
                                          const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                                          const FracPres *Afp_source,
                                          FracPres *Afp_target,
                                          MaterialAssignment *Ama_target,
                                          const kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                                          const kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                                          const kmds::Variable<bool> *AVarMixedCells_source,
                                          const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                          const kmds::Graph *AGraph,
                                          const int ANbSubPixels);


    void
    InterfaceNodesPosSmoothVF_assignGraphCut_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
                                                const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                                                const FracPres *Afp_source,
                                                MaterialAssignment *Ama_target,
                                                const kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                                                const kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                                                const kmds::Variable<bool> *AVarMixedCells_source,
                                                const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                                const kmds::Graph *AGraph,
                                                const kmds::Variable<gmds::math::Point> *AVarMidPoints_target,
                                                const int ANbSubPixels);

    void
    InterfaceNodesPosSmoothVF_assignSimulatedAnnealing_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
                                                          const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                                                          const FracPres *Afp_source,
                                                          MaterialAssignment *Ama_target,
                                                          const kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                                                          const kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                                                          const kmds::Variable<bool> *AVarMixedCells_source,
                                                          const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                                          const kmds::Graph *AGraph,
                                                          const kmds::Variable<gmds::math::Point> *AVarMidPoints_target,
                                                          const int ANbSubPixels);

    void
    InterfaceNodesPosSmoothVF_assignSimulatedAnnealing_horizontal_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
                                                                     const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                                                                     const FracPres *Afp_source,
                                                                     MaterialAssignment *Ama_target,
                                                                     const kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                                                                     const kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                                                                     const kmds::Variable<bool> *AVarMixedCells_source,
                                                                     const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                                                     const kmds::Graph *AGraph,
                                                                     const kmds::Variable<gmds::math::Point> *AVarMidPoints_target,
                                                                     const int ANbSubPixels);

    void
    InterfaceNodesPosSmoothVF_summary_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
                                         const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                                         kmds::Mesh* AMesh_source,
                                         kmds::Mesh* AMesh_target,
                                         const FracPres *Afp_source,
                                         MaterialAssignment *Ama_target,
                                         const kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                                         const kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                                         const kmds::Variable<bool> *AVarMixedCells_source,
                                         const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                         const kmds::Graph *AGraph,
                                         const kmds::Variable<double> *AVarSurfVol_source,
                                         const kmds::Variable<double> *AVarSurfVol_target,
                                         const int ANbSubPixels);

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_INTERFACENODESPOSSMOOTHVF_H_ */
