/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    ManifoldDetection.h
 *  \author  legoff
 *  \date    02/22/2018
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_MANIFOLDDETECTION_H_
#define ELG3D_MANIFOLDDETECTION_H_
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <Kokkos_UnorderedMap.hpp>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <KM/Utils/Graph.h>
#include <KM/DS/Mesh.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/DATACMPNT/FracPres.h"
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/

//    struct ManifoldDetection_options
//    {
//        bool badpillow_active;
//
//        ManifoldDetection_options(){
//            badpillow_active = false;
//        }
//    };
//
//    const ManifoldDetection_options manifolddetection_options_default;


    // WARNING hard-coded
    // number of conficts solutions stored per nodes
    const static int MANIFOLDDETECTION_NB_STORED_SOLUTIONS = 10;

    // WARNING hard-coded
    // Maximum number of potential solutions
    // higher would be too much and we exit
    //const static int MANIFOLDDETECTION_MAX_NB_SOLUTIONS = 100000;
    const static int MANIFOLDDETECTION_MAX_NB_SOLUTIONS = std::pow(5,8)+1;


    struct ManifoldDetection_configStorage
    {
        int m_conf[MANIFOLDDETECTION_NB_STORED_SOLUTIONS];
    };

    /*------------------------------------------------------------------------*/
    /** \brief  Detects the non manifold nodes on material interfaces
     *
     *  \param[in]  AInterfaceNodes the nodes to test
     *  \param[in]  AMesh the mesh
     *  \param[in]  c_N2F the node to face connectivity
     *  \param[in]  Ama the material assignment for the cells of AMesh
     *  \param[out]  ASelection a container of the nonmanifold nodes among those on the interface
     */
    void ManifoldDetection_getNonManifoldNodes_2D(kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
                                                  kmds::Mesh* AMesh,
                                                  const kmds::Connectivity* c_N2F,
                                                  elg3d::MaterialAssignment* Ama,
                                                  kmds::GrowingView<kmds::TCellID>* ASelection);

    /*------------------------------------------------------------------------*/
    /** \brief  Detects the non manifold nodes on material interfaces
     *
     *  \param[in]  AInterfaceNodes the nodes to test
     *  \param[in]  AMesh the mesh
     *  \param[in]  c_N2R the node to region connectivity
     *  \param[in]  Ama the material assignment for the cells of AMesh
     *  \param[out]  ASelection a container of the nonmanifold nodes among those on the interface
     */
    void ManifoldDetection_getNonManifoldNodes_3D(kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
                                                  kmds::Mesh* AMesh,
                                                  const kmds::Connectivity* c_N2R,
                                                  elg3d::MaterialAssignment* Ama,
                                                  kmds::GrowingView<kmds::TCellID>* ASelection);

    /*------------------------------------------------------------------------*/
    /** \brief  Detects the non manifold nodes on material interfaces
     *
     *  \param[in]  AInterfaceNodes the nodes to test
     *  \param[in]  AMesh the mesh
     *  \param[in]  c_N2R the node to region connectivity
     *  \param[in]  Ama the material assignment for the cells of AMesh
     *  \param[out]  ASelection a container of the nodes among those on the interface
     */
    void ManifoldDetection_getNonManifoldNodes_3D(kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
                                                  int AMat,
                                                  kmds::Mesh* AMesh,
                                                  const kmds::Connectivity* c_N2R,
                                                  elg3d::MaterialAssignment* Ama,
                                                  kmds::GrowingView<kmds::TCellID>* ASelection);


    /*------------------------------------------------------------------------*/
    /** \brief  Detects whether the set of cells is non-manifold
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  ACellIDs the ids of the cells
     *
     *  \return  true if the cells form a non-manifold, false if not
     */
    bool ManifoldDetection_IsNonManifoldOneMaterial_3D(kmds::Mesh* AMesh,
                                                       std::vector<kmds::TCellID> ACellIDs);

    bool ManifoldDetection_IsNonManifoldOneMaterial_3D(kmds::Mesh* AMesh,
                                                       std::vector<kmds::TCellID> ACellIDs,
                                                       std::vector<std::vector<kmds::FakeFace> >& Affs);

    bool ManifoldDetection_IsNonManifoldOneMaterial_2D(kmds::Mesh* AMesh,
                                                       std::vector<kmds::TCellID> ACellIDs,
                                                       std::vector<std::vector<kmds::FakeEdge> >& Afes);

    /*------------------------------------------------------------------------*/
    /** \brief  Resolves a non-manifold in a ball of cells
     *
     *  \param[in]  AMesh the mesh
     *  \param[in]  ACellIDs the ids of the cells
     *
     *  \return  true if the cells form a non-manifold, false if not
     */
    void ManifoldDetection_SolveNonManifoldAtOneNode_3D(kmds::Mesh* AMesh,
                                                        const kmds::TCellID AId,
                                                        const std::vector<kmds::TCellID> ACellIDs,
                                                        const FracPres* Afp,
                                                        const kmds::Variable<kmds::TCoord>* AVarVolumes,
                                                        std::vector<int>& cells2mat,
                                                        Kokkos::UnorderedMap<kmds::TCellID, ManifoldDetection_configStorage>* AStoredSolutions);

    kmds::TCoord ManifoldDetection_ComputeCost(std::vector<kmds::TCellID> ACellIDs,
                                               std::vector<int> cells2mat,
                                               const FracPres* Afp,
                                               const kmds::Variable<kmds::TCoord>* AVarVolumes);

    void ManifoldDetection_buildGraph_N_N2F2N(const kmds::GrowingView<kmds::TCellID>* ASelection,
                                              kmds::Mesh* AMesh,
                                              const kmds::Connectivity* c_N2C,
                                              kmds::Graph* AGraph,
                                              Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID>* kmap);

    void ManifoldDetection_buildGraph_N_N2R2N(const kmds::GrowingView<kmds::TCellID>* ASelection,
                                              kmds::Mesh* AMesh,
                                              const kmds::Connectivity* c_N2C,
                                              kmds::Graph* AGraph,
                                              Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID>* kmap);

    void ManifoldDetection_solveIndset_2D(const kmds::GrowingView<kmds::TCellID>* ASelection,
                                          kmds::Mesh* AMesh,
                                          const kmds::Connectivity* c_N2C,
                                          const FracPres* Afp,
                                          MaterialAssignment* Ama,
                                          const kmds::Variable<kmds::TCoord>* AVarVolumes,
                                          Kokkos::UnorderedMap<kmds::TCellID, ManifoldDetection_configStorage>* AStoredSolutions);

    void ManifoldDetection_solveIndset_3D(const kmds::GrowingView<kmds::TCellID>* ASelection,
                                          kmds::Mesh* AMesh,
                                          const kmds::Connectivity* c_N2C,
                                          const FracPres* Afp,
                                          MaterialAssignment* Ama,
                                          const kmds::Variable<kmds::TCoord>* AVarVolumes,
                                          Kokkos::UnorderedMap<kmds::TCellID, ManifoldDetection_configStorage>* AStoredSolutions);

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_MANIFOLDDETECTION_H_ */
