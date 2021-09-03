/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    SubsetProblem.cpp
 *  \author  legoff
 *  \date    01/23/2020
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/SubsetProblem.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <numeric>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <Kokkos_UnorderedMap.hpp>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/
    /*----------------------------------------------------------------------------*/
//    struct SubsetProblem_extractNodes_xD
//    {
//
//        const kmds::GrowingView<kmds::TCellID>* nodesIDs;
//        const kmds::Variable<bool>* var;
//
//        kmds::GrowingView<kmds::TCellID>* selection;
//
//        SmartLaplacian_nonFixedNodesSelection(
//
//                const kmds::GrowingView<kmds::TCellID>* nodesIDs_,
//
//                const kmds::Variable<bool>* var_,
//                kmds::GrowingView<kmds::TCellID>* selection_
//        )
//                : nodesIDs(nodesIDs_)
//                , var(var_)
//                , selection(selection_)
//        {
//        }
//
//        KOKKOS_INLINE_FUNCTION
//        void
//        operator()(int i) const {
//            int nid = nodesIDs->get(i);
//
//            if(!(*var)[nid]) {
//                selection->push_back(nid);
//            }
//        }
//    };
/*----------------------------------------------------------------------------*/
    void SubsetProblem_extractNodes_xD(const Kokkos::View<bool*, Kokkos::MemoryTraits<Kokkos::Atomic> >* AMarkNodes2keep,
                                       kmds::GrowingView<kmds::TCellID>* ANodes2keepIDs,
                                       kmds::Mesh* AMesh_ref,
                                       kmds::Mesh* AMesh_new,
                                       Kokkos::View<kmds::TCellID*>* ANodesOld2New)
    {
        const kmds::TSize nbNodes_ref_maxid = AMesh_ref->getNodeCapacity();

        kmds::TSize nbNodes_kept = 0;
        Kokkos::parallel_reduce(nbNodes_ref_maxid,
                                KOKKOS_LAMBDA(const kmds::TSize i, kmds::TSize& ls) {
                                    if((*AMarkNodes2keep)(i)) {
                                        ls = ls + 1;
                                        kmds::TCellID nid = ANodes2keepIDs->push_back(i);
                                        (*ANodesOld2New)(i) = nid;
                                    }
                                },
                                nbNodes_kept);

        AMesh_new->updateNodeCapacity(nbNodes_kept);

        AMesh_new->addNodes(nbNodes_kept);
        Kokkos::parallel_for(nbNodes_kept,
                             KOKKOS_LAMBDA(const kmds::TSize i) {
                                 gmds::math::Point pt = AMesh_ref->getNodeLocation(ANodes2keepIDs->get(i));
                                 AMesh_new->setNodeLocation(i, pt);
                             });

    }
/*----------------------------------------------------------------------------*/
    void SubsetProblem_extractFracPres_xD(const kmds::GrowingView<kmds::TCellID>* ACellIDs,
                                          const elg3d::FracPres* Afp_ref,
                                          elg3d::FracPres* Afp_new,
                                          const Kokkos::View<kmds::TCellID*>* ACellsNew2old,
                                          const std::vector<int>* AMatOld2new)
    {
        const kmds::TSize nbCells = ACellIDs->getNbElems();
        const int nbMat_ref = Afp_ref->getNbMaterials();



        for(int imat=0; imat<nbMat_ref; imat++) {
            if((*AMatOld2new)[imat] != -1) {
                Kokkos::View<double*> fpmat_value("fpmat_value", nbCells);
                Kokkos::View<bool*> fpmat_check("fpmat_check", nbCells);

                Kokkos::parallel_for(nbCells,
                                     KOKKOS_LAMBDA(const kmds::TSize i) {
                                         const kmds::TCellID cid_ref = ACellIDs->get(i);

                                         double fp = Afp_ref->getFracPres(imat, cid_ref);

                                         if(fp > 0.) {
                                             fpmat_check[i] = true;
                                             fpmat_value[i] = fp;
                                         }
                                     });

                for(kmds::TSize i_c=0; i_c<nbCells; i_c++) {
                    if(fpmat_check[i_c]) {
                        Afp_new->setFracPres((*AMatOld2new)[imat], i_c, fpmat_value[i_c]);
                    }
                }
            }
        }

    }
/*----------------------------------------------------------------------------*/
    void SubsetProblem_extract_3D(const kmds::GrowingView<kmds::TCellID>* ACellIDs,
                                  kmds::Mesh* AMesh_ref,
                                  const elg3d::FracPres* Afp_ref,
                                  kmds::Mesh* AMesh_new,
                                  elg3d::FracPres* Afp_new,
                                  const bool AKeepMaterials)
    {
        Afp_new->clear();

        // mark the nodes to keep
        const kmds::TSize nbNodes_ref_maxid = AMesh_ref->getNodeCapacity();
        Kokkos::View<bool*, Kokkos::MemoryTraits<Kokkos::Atomic> > markNodes2keep("markNodes2keep", nbNodes_ref_maxid);

        const kmds::TSize nbCells = ACellIDs->getNbElems();

        Kokkos::parallel_for(nbCells,
                             KOKKOS_LAMBDA(const kmds::TCellID i) {
                                 const kmds::TCellID cid = ACellIDs->get(i);

                                 const kmds::Region c = AMesh_ref->getRegion(cid);

                                 Kokkos::View<kmds::TCellID *> nids;
                                 c.nodeIds(nids);

                                 for(int i_n=0; i_n<nids.size(); i_n++) {
                                     // WARNING: this assignment must be atomic
                                     markNodes2keep(nids(i_n)) = true;
                                 }

                             });

        kmds::GrowingView<kmds::TCellID> nodes2keepIDs("nodes2keepIDs", nbNodes_ref_maxid);
        Kokkos::View<kmds::TCellID*> nodesOld2New("nodesOld2New", nbNodes_ref_maxid);

        SubsetProblem_extractNodes_xD(&markNodes2keep,
                                      &nodes2keepIDs,
                                      AMesh_ref,
                                      AMesh_new,
                                      &nodesOld2New);

        // extract the cells
        // WARNING: HEX valid only for full-hex mesh

        AMesh_new->updateRegionCapacity(nbCells);

        Kokkos::View<kmds::TCellID*> cellsNew2Old("cellsNewd2Old", nbCells);
        AMesh_new->addHexahedra(nbCells);
        Kokkos::parallel_for(nbCells,
                             KOKKOS_LAMBDA(const kmds::TSize i) {
                                 const kmds::TCellID cid_ref = ACellIDs->get(i);

                                 kmds::Region c_ref = AMesh_ref->getRegion(cid_ref);
                                 Kokkos::View<kmds::TCellID *> nids_ref;
                                 c_ref.nodeIds(nids_ref);

                                 kmds::TCellID nids_new[nids_ref.size()];
                                 for(int i_n=0; i_n<nids_ref.size(); i_n++) {
                                     nids_new[i_n] = nodesOld2New(nids_ref(i_n));
                                 }
                                 kmds::Region c_new = AMesh_new->getRegion(i);
                                 c_new.setNodes(nids_new, nids_ref.size());

                                 cellsNew2Old(i) = cid_ref;
                             });

        // select the frac pres to keep
        const std::map<int, std::string> matList_ref = Afp_ref->getMaterialList();
        const int nbMat = Afp_ref->getNbMaterials();

        std::vector<int> matOld2new(nbMat);

        if(AKeepMaterials) {
            Afp_new->setMaterialList(matList_ref);
            std::iota(matOld2new.begin(), matOld2new.end(), 0);
        } else {
            Kokkos::View<bool*, Kokkos::MemoryTraits<Kokkos::Atomic> > markMaterials("markMaterials", nbMat);

            Kokkos::parallel_for(nbCells,
                                 KOKKOS_LAMBDA(const kmds::TSize i) {
                                     const kmds::TCellID cid_ref = ACellIDs->get(i);

                                     for(int imat=0; imat<nbMat; imat++) {
                                         if(Afp_ref->getFracPres(imat, cid_ref) > 0.) {
                                             markMaterials(imat) = true;
                                         }
                                     }
                                 });


            std::fill(matOld2new.begin(), matOld2new.end(), -1);
            for(int imat=0; imat<nbMat; imat++) {
                if(markMaterials(imat)) {
                    matOld2new[imat] = Afp_new->createMaterial(matList_ref.at(imat));
                }
            }
        }

        // extract the frac pres
        SubsetProblem_extractFracPres_xD(ACellIDs,
                                         Afp_ref,
                                         Afp_new,
                                         &cellsNew2Old,
                                         &matOld2new);
    }
/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
