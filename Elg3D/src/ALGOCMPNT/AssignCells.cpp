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
/** \file    AssignCells.cpp
 *  \author  legoff
 *  \date    22/11/2017
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/AssignCells.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <Kokkos_Core.hpp>
#include <Kokkos_UnorderedMap.hpp>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "KM/Utils/Graph.h"
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/ManifoldDetection.h"
#include "ELG3D/ALGOCMPNT/MaterialInterfaces.h"
#include "ELG3D/ALGOCMPNT/Tools.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {

struct AssignCells_GetMajorityFracPres
{
        const kmds::GrowingView<kmds::TCellID>* cellIDs;

        const elg3d::FracPres* fp;
        elg3d::MaterialAssignment* ma;

        AssignCells_GetMajorityFracPres(const kmds::GrowingView<kmds::TCellID>* r_, const elg3d::FracPres* fp_, elg3d::MaterialAssignment* ma_)
         :
           cellIDs(r_)
         , fp(fp_)
         , ma(ma_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const
        {
            const int id = cellIDs->get(i);
            const int imat = fp->getMaxMatFracPresIndex(id);
            ma->setMaterial(imat, id);
        }
};

/*----------------------------------------------------------------------------*/
    void
    assignCells_3D(kmds::Mesh* AMesh, elg3d::FracPres* Afp, elg3d::MaterialAssignment* Ama)
    {
//        assignCellsMajorityCriteria_3D(AMesh, Afp, Ama);
//        assignCellsCorrection_3D(AMesh, Afp, Ama);
    }
    /*----------------------------------------------------------------------------*/
    void
    assignCellsMajorityCriteria_2D(kmds::Mesh* AMesh, elg3d::FracPres* Afp, elg3d::MaterialAssignment* Ama)
    {
        Ama->createMaterials(Afp->getMaterialList());

        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", AMesh->getNbFaces());
        AMesh->getFaceIDs_dummy(&cellIDs);

        Kokkos::parallel_for(cellIDs.getNbElems(), AssignCells_GetMajorityFracPres(&cellIDs, Afp, Ama));
    }
    /*----------------------------------------------------------------------------*/
    void
    assignCellsMajorityCriteria_3D(kmds::Mesh* AMesh, elg3d::FracPres* Afp, elg3d::MaterialAssignment* Ama)
    {
        Ama->createMaterials(Afp->getMaterialList());

        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", AMesh->getNbRegions());
        AMesh->getRegionIDs_dummy(&cellIDs);

        Kokkos::parallel_for(cellIDs.getNbElems(), AssignCells_GetMajorityFracPres(&cellIDs, Afp, Ama));
    }
    /*----------------------------------------------------------------------------*/
    void
    assignCellsMajorityCriteria_maintain_2D(const kmds::Mesh* AMesh,
                                            const elg3d::FracPres* Afp,
                                            elg3d::MaterialAssignment* Ama,
                                            const std::vector<std::string> AMatNames)
    {
        Ama->createMaterials(Afp->getMaterialList());

        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", AMesh->getNbFaces());
        AMesh->getFaceIDs_dummy(&cellIDs);

        std::vector<int> matId2maintain_vec;
        for(auto name: AMatNames) {
            matId2maintain_vec.push_back(Afp->getMaterialID(name));
        }
        const int nbMat2maintain = matId2maintain_vec.size();

        const kmds::TCellID nbCells = cellIDs.getNbElems();

        Kokkos::parallel_for(nbCells,
                             KOKKOS_LAMBDA(const kmds::TCellID i) {
                                 const kmds::TCellID cid = cellIDs.get(i);

                                 bool found = false;

                                 for(int imat=0; imat<nbMat2maintain; imat++) {
                                     double fp = Afp->getFracPres(matId2maintain_vec[imat], cid);

                                     if(fp > 0.) {
                                         Ama->setMaterial(matId2maintain_vec[imat], cid);
                                         found = true;
                                         break;
                                     }
                                 }
                                 if(!found) {
                                     const int imat = Afp->getMaxMatFracPresIndex(cid);
                                     Ama->setMaterial(imat, cid);
                                 }
                             });
    }
    /*----------------------------------------------------------------------------*/
    void
    assignCellsMajorityCriteria_maintain_3D(const kmds::Mesh* AMesh,
                                            const elg3d::FracPres* Afp,
                                            elg3d::MaterialAssignment* Ama,
                                            const std::vector<std::string> AMatNames)
    {
        Ama->createMaterials(Afp->getMaterialList());

        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", AMesh->getNbRegions());
        AMesh->getRegionIDs_dummy(&cellIDs);

        std::vector<int> matId2maintain_vec;
        for(auto name: AMatNames) {
            matId2maintain_vec.push_back(Afp->getMaterialID(name));
        }
        const int nbMat2maintain = matId2maintain_vec.size();

        const kmds::TCellID nbCells = cellIDs.getNbElems();

        Kokkos::parallel_for(nbCells,
                             KOKKOS_LAMBDA(const kmds::TCellID i) {
                                 const kmds::TCellID cid = cellIDs.get(i);

                                 bool found = false;

                                 for(int imat=0; imat<nbMat2maintain; imat++) {
                                     double fp = Afp->getFracPres(matId2maintain_vec[imat], cid);

                                     if(fp > 0.) {
                                         Ama->setMaterial(matId2maintain_vec[imat], cid);
                                         found = true;
                                         break;
                                     }
                                 }
                                 if(!found) {
                                     const int imat = Afp->getMaxMatFracPresIndex(cid);
                                     Ama->setMaterial(imat, cid);
                                 }
                             });
    }
    /*----------------------------------------------------------------------------*/
    void
    assignCells_defeaturing_3D(const kmds::Mesh* AMesh,
                               const kmds::Connectivity* AC_C2C_byN,
                               elg3d::MaterialAssignment* Ama)
    {
        struct assignCells_defeaturing_checkcell
        {
            const kmds::GrowingView<kmds::TCellID>* selection;
            const kmds::Connectivity* c_C2C_byN;
            const elg3d::MaterialAssignment* Ama;
            kmds::GrowingView<kmds::TCellID>* cells2Change;


            assignCells_defeaturing_checkcell(
                    const kmds::GrowingView<kmds::TCellID>* Selection_,
                    const kmds::Connectivity* c_C2C_byN_,
                    const elg3d::MaterialAssignment* Ama_,
                    kmds::GrowingView<kmds::TCellID>* cells2Change_
            )
                    : selection(Selection_)
                    , c_C2C_byN(c_C2C_byN_)
                    , Ama(Ama_)
                    , cells2Change(cells2Change_)
            {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(size_t i) const {
                const kmds::TCellID cid = selection->get(i);
                const int mat = Ama->getMaterial(cid);

                int nbCellsOtherMat = 0;

                Kokkos::View<kmds::TCellID *> cids;
                c_C2C_byN->get(cid, cids);

                for (auto i_c = 0; i_c < cids.size(); i_c++) {
                    if (mat != Ama->getMaterial(cids[i_c])) {
                        nbCellsOtherMat++;
                    }
                }

                if ((double) nbCellsOtherMat / (double) cids.size() > 0.70) {
                    cells2Change->push_back(cid);
                }
            }
        };

        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", AMesh->getNbRegions());
        AMesh->getRegionIDs_dummy(&cellIDs);

        const kmds::TCellID nbCells = cellIDs.getNbElems();

        const int nbMat = Ama->getNbMaterials();

        kmds::GrowingView<kmds::TCellID> cellIDs2Change("CELLS", AMesh->getNbRegions());


        // temporary storage for the new assignment
        Kokkos::View<int *> cells2Mat("cells2Mat", AMesh->getRegionSupID()+1);

        bool checkAssign = true;

        while(checkAssign) {

            cellIDs2Change.clear();

//            Kokkos::parallel_for(nbCells,
//                                 KOKKOS_LAMBDA(const kmds::TCellID i) {
//                                     const kmds::TCellID cid = cellIDs.get(i);
//                                     const int mat = Ama->getMaterial(cid);
//
//                                     int nbCellsOtherMat = 0;
//
//                                     Kokkos::View<kmds::TCellID *> cids;
//                                     AC_C2C_byN->get(cid, cids);
//
//                                     for (auto i_c = 0; i_c < cids.size(); i_c++) {
//                                         if (mat != Ama->getMaterial(cids[i_c])) {
//                                             nbCellsOtherMat++;
//                                         }
//                                     }
//
//                                     if (nbCellsOtherMat >= 20) {
//                                         cellIDs2Change.push_back(cid);
//                                     }
//                                 });
            Kokkos::parallel_for(nbCells,
                                 assignCells_defeaturing_checkcell(&cellIDs,
                                                                   AC_C2C_byN,
                                                                   Ama,
                                                                   &cellIDs2Change));

            const kmds::TCellID nbCellIDs2sChange = cellIDs2Change.getNbElems();

            std::cout<<"nbCellIDs2sChange "<<nbCellIDs2sChange<<std::endl;

            if(nbCellIDs2sChange == 0) {
                checkAssign = false;
            } else {

                Kokkos::parallel_for(nbCellIDs2sChange,
                                     KOKKOS_LAMBDA(const kmds::TCellID i) {
                                         const kmds::TCellID cid = cellIDs2Change.get(i);
                                         const int mat = Ama->getMaterial(cid);

                                         int matsOccurence[nbMat];
                                         for(int imat=0; imat<nbMat; imat++) {
                                             matsOccurence[imat] = 0;
                                         }

                                         Kokkos::View<kmds::TCellID *> cids;
                                         AC_C2C_byN->get(cid, cids);

                                         for (auto i_c = 0; i_c < cids.size(); i_c++) {

                                             const int imat = Ama->getMaterial(cids[i_c]);

                                             // we exclude the previously assigned mat
                                             if (mat != imat) {
                                                 matsOccurence[imat]++;
                                             }
                                         }

                                         int maxMatIndex = -1;
                                         int maxMat = -1;
                                         for(int imat=0; imat<nbMat; imat++) {
                                             if(matsOccurence[imat] > maxMat) {
                                                 maxMat = matsOccurence[imat];
                                                 maxMatIndex = imat;
                                             }
                                         }

                                         cells2Mat(cid) = maxMatIndex;
                                     });
            }

            Kokkos::parallel_for(nbCellIDs2sChange,
                                 KOKKOS_LAMBDA(const kmds::TCellID i) {
                                     const kmds::TCellID cid = cellIDs2Change.get(i);

                                     Ama->setMaterial(cells2Mat(cid), cid);
                                 });
        }  // while(checkAssign)
    }

    /*----------------------------------------------------------------------------*/
    void
    assignCells_refeaturing_XD(const kmds::GrowingView<kmds::TCellID>* ACellIDs,
                               const elg3d::FracPres* Afp,
                               const elg3d::MaterialAssignment* Ama,
                               const kmds::Connectivity* AC_C2C_byN,
                               kmds::GrowingView<kmds::TCellID>* ACellsList)
    {
        struct assignCells_refeaturing_isolated
        {
            const kmds::GrowingView<kmds::TCellID>* selection;
            const kmds::Connectivity* c_C2C_byN;
            const elg3d::MaterialAssignment* Ama;
            kmds::GrowingView<kmds::TCellID>* cells2Change;


            assignCells_refeaturing_isolated(
                    const kmds::GrowingView<kmds::TCellID>* Selection_,
                    const kmds::Connectivity* c_C2C_byN_,
                    const elg3d::MaterialAssignment* Ama_,
                    kmds::GrowingView<kmds::TCellID>* cells2Change_
            )
                    : selection(Selection_)
                    , c_C2C_byN(c_C2C_byN_)
                    , Ama(Ama_)
                    , cells2Change(cells2Change_)
            {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(size_t i) const {
                const kmds::TCellID cid = selection->get(i);
                const int mat = Ama->getMaterial(cid);

                int nbCellsOtherMat = 0;

                Kokkos::View<kmds::TCellID *> cids;
                c_C2C_byN->get(cid, cids);

                for (auto i_c = 0; i_c < cids.size(); i_c++) {
                    if (mat != Ama->getMaterial(cids[i_c])) {
                        nbCellsOtherMat++;
                    }
                }

                if ((double) nbCellsOtherMat / (double) cids.size() > 0.70) {
                    cells2Change->push_back(cid);
                }
            }
        };

        struct assignCells_refeaturing_noassign
        {
            const kmds::GrowingView<kmds::TCellID>* selection;
            const kmds::Connectivity* c_C2C_byN;
            const elg3d::FracPres* Afp;
            const elg3d::MaterialAssignment* Ama;
            const int nbMat;
            kmds::GrowingView<kmds::TCellID>* cells2Change;


            assignCells_refeaturing_noassign(
                    const kmds::GrowingView<kmds::TCellID>* Selection_,
                    const kmds::Connectivity* c_C2C_byN_,
                    const elg3d::FracPres* Afp_,
                    const elg3d::MaterialAssignment* Ama_,
                    const int nbMat_,
                    kmds::GrowingView<kmds::TCellID>* cells2Change_
            )
                    : selection(Selection_)
                    , c_C2C_byN(c_C2C_byN_)
                    , Afp(Afp_)
                    , Ama(Ama_)
                    , nbMat(nbMat_)
                    , cells2Change(cells2Change_)
            {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(size_t i) const {
                const kmds::TCellID cid = selection->get(i);
                const int mat = Ama->getMaterial(cid);

                Kokkos::View<kmds::TCellID *> cids;
                c_C2C_byN->get(cid, cids);

                for(int imat=0; imat<nbMat; imat++) {
                    if(Afp->getFracPres(imat, cid) > 0.) {
                        bool found = false;
                        if(mat == imat) {
                            found = true;
                        } else {
                            for (auto i_c = 0; i_c < cids.size(); i_c++) {
                                if (imat == Ama->getMaterial(cids[i_c])) {
                                    found = true;
                                }
                            }
                        }
                        if(!found) {
                            cells2Change->push_back(cid);
                            return;
                        }
                    }
                }

            }
        };

        const kmds::TCellID nbCells = ACellIDs->getNbElems();

        kmds::GrowingView<kmds::TCellID> cellIDs2Change("CELLS", nbCells);

        cellIDs2Change.clear();

        Kokkos::parallel_for(nbCells,
                             assignCells_refeaturing_isolated(ACellIDs,
                                                              AC_C2C_byN,
                                                              Ama,
                                                              &cellIDs2Change));

        kmds::TCellID nbCellIDs2sChange = cellIDs2Change.getNbElems();
        std::cout<<"nbCellIDs2sChange isolated_components "<<nbCellIDs2sChange<<std::endl;

        for(int icell=0; icell<nbCellIDs2sChange; icell++) {
            std::cout<<"cid "<<cellIDs2Change.get(icell)<<std::endl;
        }

        const int nbMat = Afp->getNbMaterials();
        Kokkos::parallel_for(nbCells,
                             assignCells_refeaturing_noassign(ACellIDs,
                                                              AC_C2C_byN,
                                                              Afp,
                                                              Ama,
                                                              nbMat,
                                                              &cellIDs2Change));

        nbCellIDs2sChange = cellIDs2Change.getNbElems();
        std::cout<<"nbCellIDs2sChange noassign "<<cellIDs2Change.getNbElems()<<std::endl;

        for(int icell=0; icell<nbCellIDs2sChange; icell++) {
            std::cout<<"cid "<<cellIDs2Change.get(icell)<<std::endl;
        }

        Kokkos::UnorderedMap<kmds::TCellID, void> kmap(nbCells);

        Kokkos::parallel_for(nbCellIDs2sChange, KOKKOS_LAMBDA(const int i) {
            const kmds::TCellID cid = cellIDs2Change.get(i);
            Kokkos::UnorderedMapInsertResult res = kmap.insert(cid);

            if (res.success()) {
                ACellsList->push_back(cid);
            }

            Kokkos::View<kmds::TCellID *> cids;
            AC_C2C_byN->get(cid, cids);

            for (auto i_c = 0; i_c < cids.size(); i_c++) {
                Kokkos::UnorderedMapInsertResult res = kmap.insert(cids[i_c]);
                if (res.success()) {
                    ACellsList->push_back(cids[i_c]);
                }
            }

        });

        std::cout<<"ACellsList "<<ACellsList->getNbElems()<<std::endl;
    }

    /*----------------------------------------------------------------------------*/
    void
    assignCellsCorrection_2D(kmds::Mesh* AMesh, const kmds::Connectivity* c_N2F, elg3d::FracPres* Afp, elg3d::MaterialAssignment* Ama)
    {
        // necessary init

        Kokkos::Timer timer;
        timer.reset();


        // compute and store cells volume
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", AMesh->getNbFaces());
        AMesh->getFaceIDs_dummy(&cellIDs);
        kmds::Variable<kmds::TCoord>* v = AMesh->createVariable<kmds::TCoord>(-1., kmds::KMDS_FACE, "cellVolume");
        Kokkos::parallel_for(cellIDs.getNbElems(), KOKKOS_LAMBDA(const int i) { (*v)[i] = AMesh->getFace(cellIDs.get(i)).surfvol(); });

        // non-manifold
        int nbNonManifoldNodes = 0;

        // storage of already evaluated conflicts
        // WARNING we hard-coded the number of stored solutions per node
        Kokkos::UnorderedMap<kmds::TCellID, ManifoldDetection_configStorage> storedConflictSolutions(AMesh->getNbNodes());


        std::cout<<"nonmanifold preptime "<<timer.seconds()<<std::endl;
        timer.reset();

        int nbIter = 0;
        do{

//            std::string filename_prefix("nonmanifold_resolution_step_");
//            std::string filename = filename_prefix + std::to_string(nbIter);
//            elg3d::Tools_write_2D(AMesh, Afp, Ama, filename);

            // filter nodes on interfaces between materials
            kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", AMesh->getNbNodes());
            elg3d::MaterialInterfaces_getNodeOnInterfaces(AMesh, c_N2F, Ama, &nodesInterfaces);

            std::cout<<"nonmanifold getNodeOnInterfaces "<<timer.seconds()<<" "<<nodesInterfaces.getNbElems()<<std::endl;
            timer.reset();

            // filter non-manifolds nodes
            kmds::GrowingView<kmds::TCellID> nodesNonManifold("NODES_NONMANIFOLD", nodesInterfaces.getNbElems());
            ManifoldDetection_getNonManifoldNodes_2D(&nodesInterfaces, AMesh, c_N2F, Ama, &nodesNonManifold);

            nbNonManifoldNodes = nodesNonManifold.getNbElems();

            std::cout<<"nonmanifold getNonManifoldNodes "<<timer.seconds()<<" "<<nbNonManifoldNodes<<std::endl;
            timer.reset();

            if(nbNonManifoldNodes > 0) {

                //TODO necessary if we reassign on-the-fly

                // build graph
                // WARNING we hard-coded 20 as the max number of neighbors
                kmds::Graph graph("NODES_NONMANIFOLD_GRAPH", nbNonManifoldNodes, 40);
                Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> kmapNode2Vertex(nbNonManifoldNodes);
                ManifoldDetection_buildGraph_N_N2F2N(&nodesNonManifold, AMesh, c_N2F, &graph, &kmapNode2Vertex);

                std::cout<<"nonmanifold buildGraph_N_N2F2N "<<timer.seconds()<<std::endl;
                timer.reset();

                // extract independent set
                kmds::GrowingView<kmds::TCellID> nodesNonManifold_indepset("NODES_NONMANIFOLD_INDEPENDENTSET", nodesNonManifold.getNbElems());
                graph.getIndependentSet(&nodesNonManifold_indepset);

                std::cout<<"nonmanifold getIndependentSet "<<timer.seconds()<<std::endl;
                timer.reset();

                std::cout<<"POYOP "<<nbIter++<<" "<<nbNonManifoldNodes<<" "<<nodesNonManifold_indepset.getNbElems()<<" "<<nodesNonManifold.get(nodesNonManifold_indepset.get(0))<<std::endl;


                // convert back from graph vertex to nodes IDs

                struct AssignCells_vertex2Node
                {
                    kmds::GrowingView<kmds::TCellID>* selection;
                    kmds::GrowingView<kmds::TCellID>* selection_map;

                    AssignCells_vertex2Node(kmds::GrowingView<kmds::TCellID>* selection_, kmds::GrowingView<kmds::TCellID>* selection_map_)
                            : selection(selection_)
                            , selection_map(selection_map_)
                    {
                    }

                    KOKKOS_INLINE_FUNCTION
                    void
                    operator()(int i) const {
                        int nid = selection->get(i);

                        selection->set(i, selection_map->get(nid));
                    }
                };

                Kokkos::parallel_for(nodesNonManifold_indepset.getNbElems(), AssignCells_vertex2Node(&nodesNonManifold_indepset, &nodesNonManifold));

                //TODO Steve's algorithm IMR 2010
                if(0) {
//                    std::map<int, std::string> materials = Afp->getMaterialList();
//                    for (auto mat: materials) {
//
//                        // filter nodes on interface of this material only
//                        kmds::GrowingView<kmds::TCellID> nodesInterfaces_mat("NODES_ON_INTERFACES_MAT",
//                                                                             AMesh->getNbNodes());
//                        elg3d::MaterialInterfaces_getNodeOnInterfaces(mat.first, AMesh, c_N2F, Ama,
//                                                                      &nodesInterfaces_mat);
//
//                        // filter non-manifolds nodes for this material only
//                        kmds::GrowingView<kmds::TCellID> nodesNonManifoldMat("NODES_NONMANIFOLD_MAT",
//                                                                             nodesNonManifold.getNbElems());
//                        ManifoldDetection_getNonManifoldNodes_2D(&nodesInterfaces_mat, mat.first, AMesh, c_N2F, Ama,
//                                                                 &nodesNonManifoldMat);
//
//                        // fill L_add and L_sub
//
//                        // modify frac pres
//
//                    }
//                    // reassign with majority criteria


                } else {

                    //TODO solving all materials in one go
                    ManifoldDetection_solveIndset_2D(&nodesNonManifold_indepset, AMesh, c_N2F, Afp, Ama, v, &storedConflictSolutions);

                    std::cout<<"nonmanifold solveIndset "<<timer.seconds()<<" "<<nodesNonManifold_indepset.getNbElems()<<std::endl;
                    timer.reset();

                }

            }

        } while (nbNonManifoldNodes > 0);

        // clean-up temporary data
        AMesh->deleteVariable(kmds::KMDS_FACE, "cellVolume");

    }
    /*----------------------------------------------------------------------------*/
    void
    assignCellsCorrection_3D(kmds::Mesh* AMesh, const kmds::Connectivity* c_N2R, elg3d::FracPres* Afp, elg3d::MaterialAssignment* Ama)
    {
        // necessary init

        Kokkos::Timer timer;
        timer.reset();


        // compute and store cells volume
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", AMesh->getNbRegions());
        AMesh->getRegionIDs_dummy(&cellIDs);
        kmds::Variable<kmds::TCoord>* v = AMesh->createVariable<kmds::TCoord>(-1., kmds::KMDS_REGION, "cellVolume");
        Kokkos::parallel_for(cellIDs.getNbElems(), KOKKOS_LAMBDA(const int i) { (*v)[i] = AMesh->getRegion(
                cellIDs.get(i)).surfvol(); });

        // non-manifold
        int nbNonManifoldNodes = 0;

        // storage of already evaluated conflicts
        // WARNING we hard-coded the number of stored solutions per node
        Kokkos::UnorderedMap<kmds::TCellID, ManifoldDetection_configStorage> storedConflictSolutions(AMesh->getNbNodes());


        std::cout<<"nonmanifold preptime "<<timer.seconds()<<std::endl;
        timer.reset();


        int nbIter = 0;
        do{

            std::string filename_prefix("nonmanifold_resolution_step_");
            std::string filename = filename_prefix + std::to_string(nbIter);
            elg3d::Tools_write_3D(AMesh, Afp, Ama, filename);


            // filter nodes on interfaces between materials
            kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", AMesh->getNbNodes());
            elg3d::MaterialInterfaces_getNodeOnInterfaces(AMesh, c_N2R, Ama, &nodesInterfaces);

            std::cout<<"nonmanifold getNodeOnInterfaces "<<timer.seconds()<<" "<<nodesInterfaces.getNbElems()<<std::endl;
            timer.reset();

            // filter non-manifolds nodes
            kmds::GrowingView<kmds::TCellID> nodesNonManifold("NODES_NONMANIFOLD", nodesInterfaces.getNbElems());
            ManifoldDetection_getNonManifoldNodes_3D(&nodesInterfaces, AMesh, c_N2R, Ama, &nodesNonManifold);

            nbNonManifoldNodes = nodesNonManifold.getNbElems();


            std::cout<<"nonmanifold getNonManifoldNodes "<<timer.seconds()<<" "<<nbNonManifoldNodes<<std::endl;
            timer.reset();


            if(nbNonManifoldNodes > 0) {

                //TODO necessary if we reassign on-the-fly

                // build graph
                // WARNING we hard-coded 20 as the max number of neighbors
                kmds::Graph graph("NODES_NONMANIFOLD_GRAPH", nbNonManifoldNodes, 40);
                Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> kmapNode2Vertex(nbNonManifoldNodes);
                ManifoldDetection_buildGraph_N_N2R2N(&nodesNonManifold, AMesh, c_N2R, &graph, &kmapNode2Vertex);

                std::cout<<"nonmanifold buildGraph_N_N2R2N "<<timer.seconds()<<std::endl;
                timer.reset();

                // extract independent set
                kmds::GrowingView<kmds::TCellID> nodesNonManifold_indepset("NODES_NONMANIFOLD_INDEPENDENTSET", nodesNonManifold.getNbElems());
                graph.getIndependentSet(&nodesNonManifold_indepset);

                std::cout<<"nonmanifold getIndependentSet "<<timer.seconds()<<std::endl;
                timer.reset();

                std::cout<<"POYOP "<<nbIter++<<" "<<nbNonManifoldNodes<<" "<<nodesNonManifold_indepset.getNbElems()<<" "<<nodesNonManifold.get(nodesNonManifold_indepset.get(0))<<std::endl;

                // convert back from graph vertex to nodes IDs

                struct AssignCells_vertex2Node
                {
                    kmds::GrowingView<kmds::TCellID>* selection;
                    kmds::GrowingView<kmds::TCellID>* selection_map;

                    AssignCells_vertex2Node(kmds::GrowingView<kmds::TCellID>* selection_, kmds::GrowingView<kmds::TCellID>* selection_map_)
                            : selection(selection_)
                            , selection_map(selection_map_)
                    {
                    }

                    KOKKOS_INLINE_FUNCTION
                    void
                    operator()(int i) const {
                        int nid = selection->get(i);

                        selection->set(i, selection_map->get(nid));
                    }
                };

                Kokkos::parallel_for(nodesNonManifold_indepset.getNbElems(), AssignCells_vertex2Node(&nodesNonManifold_indepset, &nodesNonManifold));

                //TODO Steve's algorithm IMR 2010
                if(0) {
                    std::map<int, std::string> materials = Afp->getMaterialList();
                    for (auto mat: materials) {

                        // filter nodes on interface of this material only
                        kmds::GrowingView<kmds::TCellID> nodesInterfaces_mat("NODES_ON_INTERFACES_MAT",
                                                                             AMesh->getNbNodes());
                        elg3d::MaterialInterfaces_getNodeOnInterfaces(mat.first, AMesh, c_N2R, Ama,
                                                                      &nodesInterfaces_mat);

                        // filter non-manifolds nodes for this material only
                        kmds::GrowingView<kmds::TCellID> nodesNonManifoldMat("NODES_NONMANIFOLD_MAT",
                                                                             nodesNonManifold.getNbElems());
                        ManifoldDetection_getNonManifoldNodes_3D(&nodesInterfaces_mat, mat.first, AMesh, c_N2R, Ama,
                                                                 &nodesNonManifoldMat);

                        // fill L_add and L_sub

                        // modify frac pres

                    }
                    // reassign with majority criteria


                } else {

                    //TODO solving all materials in one go
                    ManifoldDetection_solveIndset_3D(&nodesNonManifold_indepset, AMesh, c_N2R, Afp, Ama, v, &storedConflictSolutions);

                    std::cout<<"nonmanifold solveIndset "<<timer.seconds()<<" "<<nodesNonManifold_indepset.getNbElems()<<std::endl;
                    timer.reset();

                }

            }

        } while (nbNonManifoldNodes > 0);

        // clean-up temporary data
        AMesh->deleteVariable(kmds::KMDS_REGION, "cellVolume");

    }
/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
