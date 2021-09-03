/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    ManifoldDetection.cpp
 *  \author  legoff
 *  \date    02/22/2018
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/ManifoldDetection.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <bitset>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/math/Triangle.h>
#include <gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/BadPillowDetection.h"
#include "ELG3D/DATACMPNT/Parameters.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {


    struct ManifoldDetection_NodeIsNonManifold_2D
    {
        kmds::GrowingView<kmds::TCellID>* interfaceNodes;

        kmds::Mesh* mesh;
        const kmds::Connectivity* c_N2C;
        elg3d::MaterialAssignment* ma;
        kmds::GrowingView<kmds::TCellID>* selection;


        ManifoldDetection_NodeIsNonManifold_2D(kmds::GrowingView<kmds::TCellID>* AInterfaceNodes_,
                                               kmds::Mesh* AMesh_,
                                               const kmds::Connectivity* c_N2C_,
                                               elg3d::MaterialAssignment* ma_,
                                               kmds::GrowingView<kmds::TCellID>* Selection_)
                : interfaceNodes(AInterfaceNodes_)
                , mesh(AMesh_)
                , c_N2C(c_N2C_)
                , ma(ma_)
                , selection(Selection_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const
        {
            int nid = interfaceNodes->get(i);

            Kokkos::View<kmds::TCellID *> cellIDs;
            c_N2C->get(nid, cellIDs);

            const int nbCells = cellIDs.size();

            // preparing data
            std::vector<std::vector<kmds::FakeEdge> > fes(nbCells);
            for(int icell=0; icell<nbCells; icell++) {

                kmds::Face f = mesh->getFace(cellIDs[icell]);

                std::vector<kmds::FakeEdge> fakeEdges = f.getFakeEdges();

                fes[icell] = fakeEdges;
            }

            // we build the subsets of the cells assigned to the materials
            std::map<int, std::vector<kmds::TCellID> > mat2cells;


            for(int icell=0; icell<nbCells; icell++) {
                //mat2cells[ma->getMaterial(cellIDs[icell])].push_back(cellIDs[icell]);
                mat2cells[ma->getMaterial(cellIDs[icell])].push_back(icell);
            }

            // we study the manifoldness for each material
            for(auto mat: mat2cells) {

                if(ManifoldDetection_IsNonManifoldOneMaterial_2D(mesh, mat.second, fes)) {
                    selection->push_back(nid);

                    // node was detected as nonmanifold for one material, we can stop now
                    break;
                }

            }

        }
    };

    struct ManifoldDetection_NodeIsNonManifold_3D
    {
        kmds::GrowingView<kmds::TCellID>* interfaceNodes;

        kmds::Mesh* mesh;
        const kmds::Connectivity* c_N2C;
        elg3d::MaterialAssignment* ma;
        kmds::GrowingView<kmds::TCellID>* selection;


        ManifoldDetection_NodeIsNonManifold_3D(kmds::GrowingView<kmds::TCellID>* AInterfaceNodes_,
                                            kmds::Mesh* AMesh_,
                                             const kmds::Connectivity* c_N2C_,
                                             elg3d::MaterialAssignment* ma_,
                                             kmds::GrowingView<kmds::TCellID>* Selection_)
                : interfaceNodes(AInterfaceNodes_)
                , mesh(AMesh_)
                , c_N2C(c_N2C_)
                , ma(ma_)
                , selection(Selection_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const
        {
            int nid = interfaceNodes->get(i);

            Kokkos::View<kmds::TCellID *> cellIDs;
            c_N2C->get(nid, cellIDs);

            const int nbCells = cellIDs.size();

            // preparing data
            std::vector<std::vector<kmds::FakeFace> > ffs(nbCells);

            // TODO : BadPillowing
            std::vector<std::vector<bool> > ffIsOutward(nbCells);
            std::vector<std::vector<gmds::math::Vector> > ffNormals(nbCells);
            std::vector<std::vector<gmds::math::Triangle> > ffTriangles(nbCells);


            std::vector<std::vector<kmds::FakeFace> > ffs_internal(nbCells);
            std::map<kmds::FakeFace, std::pair<int, int> > internalFF;
            std::vector<std::vector<bool> > ffBoundary(nbCells);

            for(int icell=0; icell<nbCells; icell++) {

                kmds::Region r = mesh->getRegion(cellIDs[icell]);

                std::vector<kmds::FakeFace> fakeFaces = r.getFakeFaces();

                ffs[icell] = fakeFaces;

                ffBoundary[icell] = std::vector<bool> (fakeFaces.size(), true);

                // TODO : BadPillowing
                if(Parameters::badpillowing_detection_ON) {
                    ffIsOutward[icell] = r.getFakeFacesOrientation();

                    ffNormals[icell].resize(ffs[icell].size());
                    ffTriangles[icell].resize(ffs[icell].size());
                    for(int i_f=0; i_f<ffs[icell].size(); i_f++) {

                        if(internalFF.find(ffs[icell][i_f]) == internalFF.end()) {
                            internalFF.emplace(ffs[icell][i_f], std::pair<int, int> (icell, i_f));
                        } else {
                            ffBoundary[icell][i_f] = false;
                            ffBoundary[internalFF[ffs[icell][i_f]].first][internalFF[ffs[icell][i_f]].second] = false;
                        }

                        if (ffs[icell][i_f].hasNode(nid)) {

                            std::vector<kmds::TCellID> tri = ffs[icell][i_f].getTriangle(nid);

                            kmds::TCoord xyz_0[3];
                            kmds::TCoord xyz_1[3];
                            kmds::TCoord xyz_2[3];
                            mesh->getNodeLocation(tri[0], xyz_0[0], xyz_0[1], xyz_0[2]);
                            mesh->getNodeLocation(tri[1], xyz_1[0], xyz_1[1], xyz_1[2]);
                            mesh->getNodeLocation(tri[2], xyz_2[0], xyz_2[1], xyz_2[2]);

                            gmds::math::Point pt_0(xyz_0[0], xyz_0[1], xyz_0[2]);
                            gmds::math::Point pt_1(xyz_1[0], xyz_1[1], xyz_1[2]);
                            gmds::math::Point pt_2(xyz_2[0], xyz_2[1], xyz_2[2]);

                            gmds::math::Triangle t(pt_0, pt_1, pt_2);

                            gmds::math::Vector n = t.getNormal();
                            n.normalize();

                            ffNormals[icell][i_f] = n;
                            ffTriangles[icell][i_f] = t;

                        } else {
                            ffNormals[icell][i_f] = gmds::math::Vector(0., 0., 0.);
                        }
                    }
                }

            }

            if(nid == 13) {
                std::cout<<"nid==13 "<<"nbCells "<<nbCells<<std::endl;
            }

            // we build the subsets of the cells assigned to the materials
            std::map<int, std::vector<kmds::TCellID> > mat2cells;

            for(int icell=0; icell<nbCells; icell++) {
//                mat2cells[ma->getMaterial(cellIDs[icell])].push_back(cellIDs[icell]);
                mat2cells[ma->getMaterial(cellIDs[icell])].push_back(icell);
            }

            // we study the manifoldness for each material
            for(auto mat: mat2cells) {

//                if(ManifoldDetection_IsNonManifoldOneMaterial_3D(mesh, mat.second)) {
                if(ManifoldDetection_IsNonManifoldOneMaterial_3D(mesh, mat.second, ffs)) {
                    selection->push_back(nid);

                    // node was detected as nonmanifold for one material, we can stop now
                    break;
                }

                // TODO : BadPillowing
                bool isBadPillow = false;

                if(Parameters::badpillowing_detection_ON) {

                    switch (Parameters::badpillowing_detection_method) {
                        case Parameter_badpillowing_method::glpk :
                            if (BadPillowDetection_IsNodeBad_glpk_3D(mesh, nid, mat.second, ffs, ffIsOutward,
                                                                     ffNormals)) {
                                isBadPillow = true;
                            }
                            break;
                        case Parameter_badpillowing_method::r3d :
                            if (BadPillowDetection_IsNodeBad_r3d_3D(mesh, nid, mat.second, ffs, ffIsOutward, ffNormals,
                                                                    ffTriangles)) {
                                isBadPillow = true;
                            }
                            break;
                        case Parameter_badpillowing_method::ortho :
                            std::cout << "mat " << mat.first << std::endl;
                            if (BadPillowDetection_IsNodeBad_ortho_3D(mesh, nid, mat.second, ffs, ffIsOutward,
                                                                      ffNormals, ffTriangles, ffBoundary)) {
                                isBadPillow = true;
                            }
                            break;
                        default:
                            std::cerr
                                    << "ManifoldDetection_SolveNonManifoldAtOneNode_3D : badpillowing_detection_method incorrect."
                                    << std::endl;
                            exit(-1);
                            break;
                    }
                }

                if (isBadPillow) {
                    selection->push_back(nid);
                    break;
                }

            }

        }
    };

    struct ManifoldDetection_NodeIsNonManifoldOneMaterial_3D
    {
        kmds::GrowingView<kmds::TCellID>* interfaceNodes;

        int matIndex;
        kmds::Mesh* mesh;
        const kmds::Connectivity* c_N2C;
        elg3d::MaterialAssignment* ma;
        kmds::GrowingView<kmds::TCellID>* selection;


        ManifoldDetection_NodeIsNonManifoldOneMaterial_3D(kmds::GrowingView<kmds::TCellID>* AInterfaceNodes_,
                                               int AMat,
                                               kmds::Mesh* AMesh_,
                                               const kmds::Connectivity* c_N2C_,
                                               elg3d::MaterialAssignment* ma_,
                                               kmds::GrowingView<kmds::TCellID>* Selection_)
                : interfaceNodes(AInterfaceNodes_)
                , matIndex(AMat)
                , mesh(AMesh_)
                , c_N2C(c_N2C_)
                , ma(ma_)
                , selection(Selection_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const
        {
            int nid = interfaceNodes->get(i);

            Kokkos::View<kmds::TCellID *> cellIDs;
            c_N2C->get(nid, cellIDs);

            // we build the subset of the cells assigned to matIndex
            std::vector<kmds::TCellID> cellsofmat;

            int nbCells = cellIDs.size();
            for(int icell=0; icell<nbCells; icell++) {
                if(ma->getMaterial(cellIDs(icell)) == matIndex) {
                    cellsofmat.push_back(cellIDs(icell));
                }
            }

            // we study the manifoldness for the matIndex material
            if(ManifoldDetection_IsNonManifoldOneMaterial_3D(mesh, cellsofmat)) {
                selection->push_back(nid);
            }
        }
    };
/*----------------------------------------------------------------------------*/
    void
    ManifoldDetection_getNonManifoldNodes_2D(kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
                                             kmds::Mesh* AMesh,
                                             const kmds::Connectivity* Ac_N2C,
                                             elg3d::MaterialAssignment* Ama,
                                             kmds::GrowingView<kmds::TCellID>* ASelection)
    {
        Kokkos::parallel_for(AInterfaceNodes->getNbElems(), ManifoldDetection_NodeIsNonManifold_2D(AInterfaceNodes, AMesh, Ac_N2C, Ama, ASelection));
    }
/*----------------------------------------------------------------------------*/
    void
    ManifoldDetection_getNonManifoldNodes_3D(kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
                                             kmds::Mesh* AMesh,
                                             const kmds::Connectivity* Ac_N2C,
                                             elg3d::MaterialAssignment* Ama,
                                             kmds::GrowingView<kmds::TCellID>* ASelection)
    {
        Kokkos::parallel_for(AInterfaceNodes->getNbElems(), ManifoldDetection_NodeIsNonManifold_3D(AInterfaceNodes, AMesh, Ac_N2C, Ama, ASelection));
    }
/*----------------------------------------------------------------------------*/
    void
    ManifoldDetection_getNonManifoldNodes_3D(kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
                                             int AMat,
                                             kmds::Mesh* AMesh,
                                             const kmds::Connectivity* Ac_N2C,
                                             elg3d::MaterialAssignment* Ama,
                                             kmds::GrowingView<kmds::TCellID>* ASelection)
    {
        Kokkos::parallel_for(AInterfaceNodes->getNbElems(), ManifoldDetection_NodeIsNonManifoldOneMaterial_3D(AInterfaceNodes, AMat, AMesh, Ac_N2C, Ama, ASelection));
    }
/*----------------------------------------------------------------------------*/
    bool
    ManifoldDetection_IsNonManifoldOneMaterial_3D(kmds::Mesh* AMesh, std::vector<kmds::TCellID> ACellIDs)
    {
        // if there is only one cell, it is always a manifold
        if(ACellIDs.size() == 1) {
            return false;
        }

        // we build this set of regions FakeFaces and we compute how many there are considering unicity
        std::map<kmds::FakeFace, int> nbFakeFaces;

        for(auto cellid: ACellIDs) {

            kmds::Region r = AMesh->getRegion(cellid);

            std::vector<kmds::FakeFace> fakeFaces = r.getFakeFaces();
            for (auto ff: fakeFaces) {
                nbFakeFaces[ff]++;
            }
        }

        // we compute the Euler number :
        // we keep the faces that appeared only once (non-internal faces)
        // and we get their edges and nodes
        int nbInterfaceFaces = 0;
        std::set<kmds::FakeEdge> interfaceEdges;
        std::set<kmds::TCellID> interfaceNodes;
        for (auto ff: nbFakeFaces) {
            if (ff.second < 2) {
                nbInterfaceFaces++;

                std::vector<kmds::FakeEdge> fakeEdges = ff.first.getFakeEdges();
                for (auto fe: fakeEdges) {
                    interfaceEdges.insert(fe);
                }

                std::vector<kmds::TCellID> nids = ff.first.node_ids();
                for (auto n: nids) {
                    interfaceNodes.insert(n);
                }
            }
        }

        int nbInterfaceEdges = interfaceEdges.size();
        int nbInterfaceNodes = interfaceNodes.size();

        if (nbInterfaceNodes - nbInterfaceEdges + nbInterfaceFaces != 2) {
            return true;
        }

        return false;
    }
/*----------------------------------------------------------------------------*/
    bool
    ManifoldDetection_IsNonManifoldOneMaterial_3D(kmds::Mesh* AMesh,
                                                  std::vector<kmds::TCellID> ACellIDs,
                                                  std::vector<std::vector<kmds::FakeFace> >& Affs)
    {
        // if there is only one cell, it is always a manifold
        if(ACellIDs.size() == 1) {
            return false;
        }

        // we build this set of regions FakeFaces and we compute how many there are considering unicity
        std::map<kmds::FakeFace, int> nbFakeFaces;

        for(auto cellid: ACellIDs) {

            for (auto ff: Affs[cellid]) {
                nbFakeFaces[ff]++;
            }
        }

        // we compute the Euler number :
        // we keep the faces that appeared only once (non-internal faces)
        // and we get their edges and nodes
        int nbInterfaceFaces = 0;
        std::set<kmds::FakeEdge> interfaceEdges;
        std::set<kmds::TCellID> interfaceNodes;
        for (auto ff: nbFakeFaces) {
            if (ff.second < 2) {
                nbInterfaceFaces++;

                std::vector<kmds::FakeEdge> fakeEdges = ff.first.getFakeEdges();
                for (auto fe: fakeEdges) {
                    interfaceEdges.insert(fe);
                }

                std::vector<kmds::TCellID> nids = ff.first.node_ids();
                for (auto n: nids) {
                    interfaceNodes.insert(n);
                }
            }
        }

        int nbInterfaceEdges = interfaceEdges.size();
        int nbInterfaceNodes = interfaceNodes.size();

        if (nbInterfaceNodes - nbInterfaceEdges + nbInterfaceFaces != 2) {
            return true;
        }

        return false;
    }
/*----------------------------------------------------------------------------*/
    bool
    ManifoldDetection_IsNonManifoldOneMaterial_2D(kmds::Mesh* AMesh,
                                                  std::vector<kmds::TCellID> ACellIDs,
                                                  std::vector<std::vector<kmds::FakeEdge> >& Afes)
    {
        // if there is only one cell, it is always a manifold
        if(ACellIDs.size() == 1) {
            return false;
        }

        // we build this set of regions FakeFaces and we compute how many there are considering unicity
        std::map<kmds::FakeEdge, int> nbFakeEdges;

        for(auto cellid: ACellIDs) {

            for (auto fe: Afes[cellid]) {
                nbFakeEdges[fe]++;
            }
        }

        // we compute the Euler number :
        // we keep the faces that appeared only once (non-internal faces)
        // and we get their edges and nodes
        int nbInterfaceEdges = 0;
        std::set<kmds::TCellID> interfaceNodes;
        for (auto fe: nbFakeEdges) {
            if (fe.second < 2) {
                nbInterfaceEdges++;

                std::vector<kmds::TCellID> nids = fe.first.node_ids();
                for (auto n: nids) {
                    interfaceNodes.insert(n);
                }
            }
        }

        int nbInterfaceNodes = interfaceNodes.size();

        if (nbInterfaceNodes - nbInterfaceEdges != 0) {
            return true;
        }

        return false;
    }
/*----------------------------------------------------------------------------*/
    void
    ManifoldDetection_SolveNonManifoldAtOneNode_3D(kmds::Mesh* AMesh,
                                                   const kmds::TCellID AId,
                                                   const std::vector<kmds::TCellID> ACellIDs,
                                                   const FracPres* Afp,
                                                   const kmds::Variable<kmds::TCoord>* AVarVolumes,
                                                   std::vector<int>& cells2mat,
                                                   Kokkos::UnorderedMap<kmds::TCellID, ManifoldDetection_configStorage>* AStoredSolutions)
    {
        int nbCells = ACellIDs.size();
        cells2mat.resize(nbCells);

        // get the possible set of materials
        std::set<int> materialSet;

        // WARNING the set of possible materials is formed from the fracpres because we suppose that it does not change
        // WARNING between two calls to this function, otherwise storing the best solutions would make no sense
//        for(auto mat: cells2mat) {
//            materialSet.insert(mat);
//        }
        for (auto c: ACellIDs) {
            for (int mat = 0; mat < Afp->getNbMaterials(); mat++) {
                if (Afp->getFracPres(mat, c) > 0.) {
                    materialSet.insert(mat);
                }
            }
        }

        int nbMaterials = materialSet.size();

        // order the materials
        std::vector<int> matIDs;
        for (auto mat: materialSet) {
            matIDs.push_back(mat);
        }


        // first we check whether this node was already treated
        // if so, we have its best solutions stored
        int index = AStoredSolutions->find(AId);
        if(AStoredSolutions->valid_at(index)) {

            ManifoldDetection_configStorage storage = AStoredSolutions->value_at(index);

            int iconf;
            bool found = false;
            for(int istore=0; istore<MANIFOLDDETECTION_NB_STORED_SOLUTIONS; istore++) {

                // looking for an available solution
                if(storage.m_conf[istore] != -1) {
                    iconf = storage.m_conf[istore];
                    found = true;
                    storage.m_conf[istore] = -1; // we disable this solution
                    AStoredSolutions->value_at(index) = storage;
                    break;
                }
            }

            // build cell assignment for this config

            // conversion into base nbMaterials
            int remain = iconf;

            for (int icell = nbCells - 1; icell >= 0; icell--) {
                int div = std::floor((double) remain / (double) std::pow(nbMaterials, icell));
                cells2mat[icell] = matIDs[div];
                remain -= div * std::pow(nbMaterials, icell);
            }

            if(!found) {
                std::cerr<<"ManifoldDetection_SolveNonManifoldAtOneNode_3D no remaining suitable solution"<<std::endl;
                exit(-1);
            }

        } else {
            // we compute and store the best solutions

            // initialization
            kmds::TCoord min_cost_current = HUGE_VALF;

            if (std::pow(nbMaterials, nbCells) > MANIFOLDDETECTION_MAX_NB_SOLUTIONS) {
                std::cout << "The number of potential solutions is too high" << std::endl;
                exit(-1);
            }

            std::vector<kmds::TCoord> storedCosts(MANIFOLDDETECTION_NB_STORED_SOLUTIONS, HUGE_VALF);
            ManifoldDetection_configStorage storage;
            for(int istore=0; istore<MANIFOLDDETECTION_NB_STORED_SOLUTIONS; istore++) {
                storage.m_conf[istore] = -1;
            }

            // preparing data
            std::vector<std::vector<kmds::FakeFace> > ffs(nbCells);

            // TODO : BadPillowing
            std::vector<std::vector<bool> > ffIsOutward(nbCells);
            std::vector<std::vector<gmds::math::Vector> > ffNormals(nbCells);
            std::vector<std::vector<gmds::math::Triangle> > ffTriangles(nbCells);


            std::vector<std::vector<kmds::FakeFace> > ffs_internal(nbCells);
            std::map<kmds::FakeFace, std::pair<int, int> > internalFF;
            std::vector<std::vector<bool> > ffBoundary(nbCells);


            for(int icell=0; icell<nbCells; icell++) {

                kmds::Region r = AMesh->getRegion(ACellIDs[icell]);

                std::vector<kmds::FakeFace> fakeFaces = r.getFakeFaces();

                ffs[icell] = fakeFaces;

                ffBoundary[icell] = std::vector<bool> (fakeFaces.size(), true);

                // TODO : BadPillowing
                if(Parameters::badpillowing_detection_ON) {

                    ffIsOutward[icell] = r.getFakeFacesOrientation();

                    ffNormals[icell].resize(ffs[icell].size());
                    ffTriangles[icell].resize(ffs[icell].size());
                    for (int i_f = 0; i_f < ffs[icell].size(); i_f++) {


                        if(internalFF.find(ffs[icell][i_f]) == internalFF.end()) {
                            internalFF.emplace(ffs[icell][i_f], std::pair<int, int> (icell, i_f));
                        } else {
                            ffBoundary[icell][i_f] = false;
                            ffBoundary[internalFF[ffs[icell][i_f]].first][internalFF[ffs[icell][i_f]].second] = false;
                        }


                        std::vector<kmds::TCellID> tri = ffs[icell][i_f].getTriangle(AId);

                        kmds::TCoord xyz_0[3];
                        kmds::TCoord xyz_1[3];
                        kmds::TCoord xyz_2[3];
                        AMesh->getNodeLocation(tri[0], xyz_0[0], xyz_0[1], xyz_0[2]);
                        AMesh->getNodeLocation(tri[1], xyz_1[0], xyz_1[1], xyz_1[2]);
                        AMesh->getNodeLocation(tri[2], xyz_2[0], xyz_2[1], xyz_2[2]);

                        gmds::math::Point pt_0(xyz_0[0], xyz_0[1], xyz_0[2]);
                        gmds::math::Point pt_1(xyz_1[0], xyz_1[1], xyz_1[2]);
                        gmds::math::Point pt_2(xyz_2[0], xyz_2[1], xyz_2[2]);

                        gmds::math::Triangle t(pt_0, pt_1, pt_2);

                        gmds::math::Vector n = t.getNormal();
                        n.normalize();

                        ffNormals[icell][i_f] = n;
                        ffTriangles[icell][i_f] = t;
                    }
                }
            }

            int nbConfigs = std::pow(nbMaterials, nbCells);

//            std::cout<<"nbConfigs "<<nbConfigs<<" "<<nbMaterials<<" "<<nbCells<<std::endl;

            // brute force
            for (int iconf = 0; iconf < nbConfigs; iconf++) {

                // build cell assignment for this config
                std::vector<int> cells2mat_tmp(nbCells, 0);

                // conversion into base nbMaterials
                int remain = iconf;

                for (int icell = nbCells - 1; icell >= 0; icell--) {
                    int div = std::floor((double) remain / (double) std::pow(nbMaterials, icell));
                    cells2mat_tmp[icell] = matIDs[div];
                    remain -= div * std::pow(nbMaterials, icell);
                }

            //std::cout<<"iconfig "<<iconf<<std::endl;

//                for (int icell = 0; icell < nbCells; icell++) {
//                    std::cout << cells2mat_tmp[icell] << " ";
//                }
//                std::cout << std::endl;

                // check non-manifold
                // we build the subsets of the cells assigned to the materials
                std::map<int, std::vector<kmds::TCellID> > mat2cells;

                for (int icell = 0; icell < nbCells; icell++) {

//                    mat2cells[cells2mat_tmp[icell]].push_back(ACellIDs[icell]);
                    mat2cells[cells2mat_tmp[icell]].push_back(icell);
                }

                // compute cost
                kmds::TCoord cost_current = ManifoldDetection_ComputeCost(ACellIDs, cells2mat_tmp, Afp, AVarVolumes);

                // study config only if it is not worst than the worst stored one
                if (cost_current <= storedCosts[MANIFOLDDETECTION_NB_STORED_SOLUTIONS - 1]) {


                    // we study the manifoldness for each material
                    bool isNonManifold = false;
                    for (auto mat: mat2cells) {
                        if (ManifoldDetection_IsNonManifoldOneMaterial_3D(AMesh, mat.second, ffs)) {

                            isNonManifold = true;
                            break;
                        }

                        // TODO : BadPillowing
                        if(Parameters::badpillowing_detection_ON) {

                            switch (Parameters::badpillowing_detection_method) {
                                case Parameter_badpillowing_method::glpk :
                                    if (BadPillowDetection_IsNodeBad_glpk_3D(AMesh, AId, mat.second, ffs, ffIsOutward, ffNormals)) {
                                        isNonManifold = true;
                                    }
                                    break;
                                case Parameter_badpillowing_method::r3d :
                                    if (BadPillowDetection_IsNodeBad_r3d_3D(AMesh, AId, mat.second, ffs, ffIsOutward, ffNormals, ffTriangles)) {
                                        isNonManifold = true;
                                    }
                                    break;
                                case Parameter_badpillowing_method::ortho :
                                    if (BadPillowDetection_IsNodeBad_ortho_3D(AMesh, AId, mat.second, ffs, ffIsOutward,
                                                                              ffNormals, ffTriangles, ffBoundary)) {
                                        isNonManifold = true;

                                    }
                                    break;
                                default:
                                    std::cerr<<"ManifoldDetection_SolveNonManifoldAtOneNode_3D : badpillowing_detection_method incorrect."<<std::endl;
                                    exit(-1);
                                    break;
                            }
                            if (isNonManifold) {
                                break;
                            }
                        }
                    }

                    // when non-manifold, we ignore this configuration
                    if (isNonManifold) {
                        continue;
                    }


//            std::cout<<"cost "<<cost_current<<std::endl;
//
//            if(min_cost_current > cost_current) {
//                min_cost_current = cost_current;
//            }

                    // find the position this solution should be put at
                    for (int istore = 0; istore < MANIFOLDDETECTION_NB_STORED_SOLUTIONS; istore++) {
                        if (storedCosts[istore] > cost_current) {

                            for (int ileft = MANIFOLDDETECTION_NB_STORED_SOLUTIONS - 1; ileft > istore; ileft--) {
                                storedCosts[ileft] = storedCosts[ileft - 1];
                                storage.m_conf[ileft] = storage.m_conf[ileft - 1];
                            }

                            storedCosts[istore] = cost_current;
                            storage.m_conf[istore] = iconf;
                            break;
                        }
                    }

                }

            }

            // this is the best solution
            int remain = storage.m_conf[0];

            // it is taken, so remove it from the solutions
            storage.m_conf[0] = -1;

            // we store the best solutions
            AStoredSolutions->insert(AId, storage);

            // conversion into base nbMaterials
            // build cell assignment for this config
            for (int icell = nbCells - 1; icell >= 0; icell--) {
                int div = std::floor((double) remain / (double) std::pow(nbMaterials, icell));
                cells2mat[icell] = matIDs[div];
                remain -= div * std::pow(nbMaterials, icell);
            }

        }

//        std::cout<<"min_cost_current "<<min_cost_current<<std::endl;

    }
/*----------------------------------------------------------------------------*/
    void
    ManifoldDetection_SolveNonManifoldAtOneNode_2D(kmds::Mesh* AMesh,
                                                   const kmds::TCellID AId,
                                                   const std::vector<kmds::TCellID> ACellIDs,
                                                   const FracPres* Afp,
                                                   const kmds::Variable<kmds::TCoord>* AVarVolumes,
                                                   std::vector<int>& cells2mat,
                                                   Kokkos::UnorderedMap<kmds::TCellID, ManifoldDetection_configStorage>* AStoredSolutions)
    {
        int nbCells = ACellIDs.size();
        cells2mat.resize(nbCells);

        // get the possible set of materials
        std::set<int> materialSet;

        // WARNING the set of possible materials is formed from the fracpres because we suppose that it does not change
        // WARNING between two calls to this function, otherwise storing the best solutions would make no sense
//        for(auto mat: cells2mat) {
//            materialSet.insert(mat);
//        }
        for (auto c: ACellIDs) {
            for (int mat = 0; mat < Afp->getNbMaterials(); mat++) {
                if (Afp->getFracPres(mat, c) > 0.) {
                    materialSet.insert(mat);
                }
            }
        }

        int nbMaterials = materialSet.size();

        // order the materials
        std::vector<int> matIDs;
        for (auto mat: materialSet) {
            matIDs.push_back(mat);
        }


        // first we check whether this node was already treated
        // if so, we have its best solutions stored
        int index = AStoredSolutions->find(AId);
        if(AStoredSolutions->valid_at(index)) {

            ManifoldDetection_configStorage storage = AStoredSolutions->value_at(index);

            int iconf;
            bool found = false;
            for(int istore=0; istore<MANIFOLDDETECTION_NB_STORED_SOLUTIONS; istore++) {

                // looking for an available solution
                if(storage.m_conf[istore] != -1) {
                    iconf = storage.m_conf[istore];
                    found = true;
                    storage.m_conf[istore] = -1; // we disable this solution
                    AStoredSolutions->value_at(index) = storage;
                    break;
                }
            }

            // build cell assignment for this config

            // conversion into base nbMaterials
            int remain = iconf;

            for (int icell = nbCells - 1; icell >= 0; icell--) {
                int div = std::floor((double) remain / (double) std::pow(nbMaterials, icell));
                cells2mat[icell] = matIDs[div];
                remain -= div * std::pow(nbMaterials, icell);
            }

            if(!found) {
                std::cerr<<"ManifoldDetection_SolveNonManifoldAtOneNode_2D no remaining suitable solution"<<std::endl;
                exit(-1);
            }

        } else {
            // we compute and store the best solutions

            // initialization
            kmds::TCoord min_cost_current = HUGE_VALF;

            if (std::pow(nbMaterials, nbCells) > MANIFOLDDETECTION_MAX_NB_SOLUTIONS) {
                std::cout << "The number of potential solutions is too high" << std::endl;
                exit(-1);
            }

            std::vector<kmds::TCoord> storedCosts(MANIFOLDDETECTION_NB_STORED_SOLUTIONS, HUGE_VALF);
            ManifoldDetection_configStorage storage;
            for(int istore=0; istore<MANIFOLDDETECTION_NB_STORED_SOLUTIONS; istore++) {
                storage.m_conf[istore] = -1;
            }

            // preparing data
            std::vector<std::vector<kmds::FakeEdge> > fes(nbCells);
            for(int icell=0; icell<nbCells; icell++) {

                kmds::Face f = AMesh->getFace(ACellIDs[icell]);

                std::vector<kmds::FakeEdge> fakeEdges = f.getFakeEdges();

                fes[icell] = fakeEdges;
            }

            int nbConfigs = std::pow(nbMaterials, nbCells);

//            std::cout<<"nbConfigs "<<nbConfigs<<" "<<nbMaterials<<" "<<nbCells<<std::endl;

            // brute force
            for (int iconf = 0; iconf < nbConfigs; iconf++) {

                // build cell assignment for this config
                std::vector<int> cells2mat_tmp(nbCells, 0);

                // conversion into base nbMaterials
                int remain = iconf;

                for (int icell = nbCells - 1; icell >= 0; icell--) {
                    int div = std::floor((double) remain / (double) std::pow(nbMaterials, icell));
                    cells2mat_tmp[icell] = matIDs[div];
                    remain -= div * std::pow(nbMaterials, icell);
                }

                //std::cout<<"iconfig "<<iconf<<std::endl;

//                for (int icell = 0; icell < nbCells; icell++) {
//                    std::cout << cells2mat_tmp[icell] << " ";
//                }
//                std::cout << std::endl;

                // check non-manifold
                // we build the subsets of the cells assigned to the materials
                std::map<int, std::vector<kmds::TCellID> > mat2cells;

                for (int icell = 0; icell < nbCells; icell++) {

//                    mat2cells[cells2mat_tmp[icell]].push_back(ACellIDs[icell]);
                    mat2cells[cells2mat_tmp[icell]].push_back(icell);
                }

                // compute cost
                kmds::TCoord cost_current = ManifoldDetection_ComputeCost(ACellIDs, cells2mat_tmp, Afp, AVarVolumes);

                // study config only if it is not worst than the worst stored one
                if (cost_current <= storedCosts[MANIFOLDDETECTION_NB_STORED_SOLUTIONS - 1]) {


                    // we study the manifoldness for each material
                    bool isNonManifold = false;
                    for (auto mat: mat2cells) {
                        if (ManifoldDetection_IsNonManifoldOneMaterial_2D(AMesh, mat.second, fes)) {
//                    if (ManifoldDetection_IsNonManifoldOneMaterial_3D(AMesh, mat.second)) {
                            isNonManifold = true;
                            break;
                        }
                    }

                    // when non-manifold, we ignore this configuration
                    if (isNonManifold) {
                        continue;
                    }


//            std::cout<<"cost "<<cost_current<<std::endl;
//
//            if(min_cost_current > cost_current) {
//                min_cost_current = cost_current;
//            }

                    // find the position this solution should be put at
                    for (int istore = 0; istore < MANIFOLDDETECTION_NB_STORED_SOLUTIONS; istore++) {
                        if (storedCosts[istore] > cost_current) {

                            for (int ileft = MANIFOLDDETECTION_NB_STORED_SOLUTIONS - 1; ileft > istore; ileft--) {
                                storedCosts[ileft] = storedCosts[ileft - 1];
                                storage.m_conf[ileft] = storage.m_conf[ileft - 1];
                            }

                            storedCosts[istore] = cost_current;
                            storage.m_conf[istore] = iconf;
                            break;
                        }
                    }

                }

            }

            // this is the best solution
            int remain = storage.m_conf[0];

            // it is taken, so remove it from the solutions
            storage.m_conf[0] = -1;

            // we store the best solutions
            AStoredSolutions->insert(AId, storage);

            // conversion into base nbMaterials
            // build cell assignment for this config
            for (int icell = nbCells - 1; icell >= 0; icell--) {
                int div = std::floor((double) remain / (double) std::pow(nbMaterials, icell));
                cells2mat[icell] = matIDs[div];
                remain -= div * std::pow(nbMaterials, icell);
            }

        }

//        std::cout<<"min_cost_current "<<min_cost_current<<std::endl;

    }
/*----------------------------------------------------------------------------*/
    kmds::TCoord
    ManifoldDetection_ComputeCost(std::vector<kmds::TCellID> ACellIDs,
                                  std::vector<int> cells2mat,
                                  const FracPres* Afp,
                                  const kmds::Variable<kmds::TCoord>* AVarVolumes)
    {
        kmds::TCoord cost_tot = 0.;

        for(int icell=0; icell<ACellIDs.size(); icell++) {

            kmds::TCoord diffFracPres = 1. - Afp->getFracPres(cells2mat[icell], ACellIDs[icell]);;
            kmds::TCoord cost = diffFracPres * (*AVarVolumes)[icell];

            cost_tot += cost;
        }

        return cost_tot;
    }
/*----------------------------------------------------------------------------*/
    void
    ManifoldDetection_buildGraph_N_N2F2N(const kmds::GrowingView<kmds::TCellID>* ASelection,
                                         kmds::Mesh* AMesh,
                                         const kmds::Connectivity* c_N2C,
                                         kmds::Graph* AGraph,
                                         Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID>* kmap)
    {
        struct ManifoldDetection_buildNodeNeighbors
        {
            const kmds::GrowingView<kmds::TCellID>* selection;

            const Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID>* kmap;
            kmds::Mesh* mesh;
            const kmds::Connectivity* c_N2C;

            kmds::Graph* graph;


            ManifoldDetection_buildNodeNeighbors(
                    const kmds::GrowingView<kmds::TCellID>* Selection_,
                    const Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID>* kmap_,
                    kmds::Mesh* AMesh_,
                    const kmds::Connectivity* c_N2C_,
                    kmds::Graph* graph_
            )
                    : selection(Selection_)
                    , kmap(kmap_)
                    , mesh(AMesh_)
                    , c_N2C(c_N2C_)
                    , graph(graph_)
            {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(int i) const {
                int nid = selection->get(i);

                Kokkos::View<kmds::TCellID*> cids;
                c_N2C->get(nid, cids);

                std::set<kmds::TCellID> dist_N2C2N_nodes;

                for(int ic=0; ic<cids.size(); ic++) {

                    kmds::Face c = mesh->getFace(cids[ic]);

                    Kokkos::View<kmds::TCellID*> nids_tmp;
                    c.nodeIds(nids_tmp);

                    for(int in=0; in<nids_tmp.size(); in++) {
                        if(nids_tmp[in] != nid) {

                            // we keep as neighbors the nodes in the selection only
                            kmds::TCellID index = kmap->find(nids_tmp[in]);
                            if(kmap->valid_at(index)) {
                                dist_N2C2N_nodes.insert(kmap->value_at(index));
                            }
                        }
                    }
                }

//                graph->setNeighbors(nid, dist_N2C2N_nodes);
                graph->setNeighbors(i, dist_N2C2N_nodes);
            }
        };

        struct ManifoldDetection_buildNodeIndexMappinqg
        {
            const kmds::GrowingView<kmds::TCellID>* selection;
            Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID>* kmap;

            ManifoldDetection_buildNodeIndexMappinqg(
                    const kmds::GrowingView<kmds::TCellID>* selection_,
                    Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID>* kmap_
            )
                    : selection(selection_)
                    , kmap(kmap_)
            {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(int i) const {
                int nid = selection->get(i);

                kmap->insert(nid,i);
            }
        };


        kmds::TCellID nbVec = ASelection->getNbElems();

        std::cout<<"nbVec "<<nbVec<<std::endl;

        // build mapping between nodes ids and graph vertices indices
        Kokkos::parallel_for(nbVec, ManifoldDetection_buildNodeIndexMappinqg(ASelection, kmap));

        Kokkos::parallel_for(nbVec, ManifoldDetection_buildNodeNeighbors(ASelection, kmap, AMesh, c_N2C, AGraph));
    }
/*----------------------------------------------------------------------------*/
    void
    ManifoldDetection_buildGraph_N_N2R2N(const kmds::GrowingView<kmds::TCellID>* ASelection,
                                         kmds::Mesh* AMesh,
                                         const kmds::Connectivity* c_N2C,
                                         kmds::Graph* AGraph,
                                         Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID>* kmap)
    {
        struct ManifoldDetection_buildNodeNeighbors
        {
            const kmds::GrowingView<kmds::TCellID>* selection;

            const Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID>* kmap;
            kmds::Mesh* mesh;
            const kmds::Connectivity* c_N2C;

            kmds::Graph* graph;


            ManifoldDetection_buildNodeNeighbors(
                    const kmds::GrowingView<kmds::TCellID>* Selection_,
                    const Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID>* kmap_,
                    kmds::Mesh* AMesh_,
                    const kmds::Connectivity* c_N2C_,
                    kmds::Graph* graph_
            )
                    : selection(Selection_)
                    , kmap(kmap_)
                    , mesh(AMesh_)
                    , c_N2C(c_N2C_)
                    , graph(graph_)
            {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(int i) const {
                int nid = selection->get(i);

                Kokkos::View<kmds::TCellID*> cids;
                c_N2C->get(nid, cids);

                std::set<kmds::TCellID> dist_N2C2N_nodes;

                for(int ic=0; ic<cids.size(); ic++) {

                    kmds::Region c = mesh->getRegion(cids[ic]);

                    Kokkos::View<kmds::TCellID*> nids_tmp;
                    c.nodeIds(nids_tmp);

                    for(int in=0; in<nids_tmp.size(); in++) {
                        if(nids_tmp[in] != nid) {

                            // we keep as neighbors the nodes in the selection only
                            kmds::TCellID index = kmap->find(nids_tmp[in]);
                            if(kmap->valid_at(index)) {
                                dist_N2C2N_nodes.insert(kmap->value_at(index));
                            }
                        }
                    }
                }

//                graph->setNeighbors(nid, dist_N2C2N_nodes);
                graph->setNeighbors(i, dist_N2C2N_nodes);
            }
        };

        struct ManifoldDetection_buildNodeIndexMappinqg
        {
            const kmds::GrowingView<kmds::TCellID>* selection;
            Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID>* kmap;

            ManifoldDetection_buildNodeIndexMappinqg(
                    const kmds::GrowingView<kmds::TCellID>* selection_,
                    Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID>* kmap_
            )
                    : selection(selection_)
                    , kmap(kmap_)
            {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(int i) const {
                int nid = selection->get(i);

                kmap->insert(nid,i);
            }
        };


        kmds::TCellID nbVec = ASelection->getNbElems();

        std::cout<<"nbVec "<<nbVec<<std::endl;

        // build mapping between nodes ids and graph vertices indices
        Kokkos::parallel_for(nbVec, ManifoldDetection_buildNodeIndexMappinqg(ASelection, kmap));

        Kokkos::parallel_for(nbVec, ManifoldDetection_buildNodeNeighbors(ASelection, kmap, AMesh, c_N2C, AGraph));
    }
/*----------------------------------------------------------------------------*/
    void
    ManifoldDetection_solveIndset_3D(const kmds::GrowingView<kmds::TCellID>* ASelection,
                                     kmds::Mesh* AMesh,
                                     const kmds::Connectivity* c_N2C,
                                     const FracPres* Afp,
                                     MaterialAssignment* Ama,
                                     const kmds::Variable<kmds::TCoord>* AVarVolumes,
                                     Kokkos::UnorderedMap<kmds::TCellID, ManifoldDetection_configStorage>* AStoredSolutions)
    {
        struct  ManifoldDetection_callForOneNode_3D
        {
            const kmds::GrowingView<kmds::TCellID>* selection;
            kmds::Mesh* mesh;
            const kmds::Connectivity* c_N2C;
            const FracPres* fp;
            MaterialAssignment* ma;
            const kmds::Variable<kmds::TCoord>* varVolumes;
            Kokkos::UnorderedMap<kmds::TCellID, ManifoldDetection_configStorage>* storedSolutions;

            ManifoldDetection_callForOneNode_3D(const kmds::GrowingView<kmds::TCellID>* selection_,
                                                kmds::Mesh* mesh_,
                                                const kmds::Connectivity* c_N2C_,
                                                const FracPres* fp_,
                                                MaterialAssignment* ma_,
                                                const kmds::Variable<kmds::TCoord>* varVolumes_,
                                                Kokkos::UnorderedMap<kmds::TCellID, ManifoldDetection_configStorage>* storedSolutions_
            )
                    : selection(selection_)
                    , mesh(mesh_)
                    , c_N2C(c_N2C_)
                    , fp(fp_)
                    , ma(ma_)
                    , varVolumes(varVolumes_)
                    , storedSolutions(storedSolutions_)
            {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(int i) const
            {
                int nid = selection->get(i);

                // form the ball of cells adjacent to nid
                Kokkos::View<kmds::TCellID *> cellIDs_tmp;
                c_N2C->get(nid, cellIDs_tmp);

                std::vector<kmds::TCellID> cellIDs(cellIDs_tmp.size());
                for(int icell=0; icell<cellIDs_tmp.size(); icell++) {
                    cellIDs[icell] = cellIDs_tmp[icell];
                }

                std::vector<int> cells2mat(cellIDs_tmp.size());

                // get a solution to the non-manifold
                ManifoldDetection_SolveNonManifoldAtOneNode_3D(mesh,
                                                               nid,
                                                               cellIDs,
                                                               fp,
                                                               varVolumes,
                                                               cells2mat,
                                                               storedSolutions);

                // reassign the cells
                for(int icell=0; icell<cellIDs_tmp.size(); icell++) {
                    ma->setMaterial(cells2mat[icell], cellIDs[icell]);
                }

            }

        };

        Kokkos::parallel_for(ASelection->getNbElems(), ManifoldDetection_callForOneNode_3D(ASelection, AMesh, c_N2C, Afp, Ama, AVarVolumes, AStoredSolutions));
    }
/*----------------------------------------------------------------------------*/
    void
    ManifoldDetection_solveIndset_2D(const kmds::GrowingView<kmds::TCellID>* ASelection,
                                     kmds::Mesh* AMesh,
                                     const kmds::Connectivity* c_N2C,
                                     const FracPres* Afp,
                                     MaterialAssignment* Ama,
                                     const kmds::Variable<kmds::TCoord>* AVarVolumes,
                                     Kokkos::UnorderedMap<kmds::TCellID, ManifoldDetection_configStorage>* AStoredSolutions)
    {
        struct  ManifoldDetection_callForOneNode_2D
        {
            const kmds::GrowingView<kmds::TCellID>* selection;
            kmds::Mesh* mesh;
            const kmds::Connectivity* c_N2C;
            const FracPres* fp;
            MaterialAssignment* ma;
            const kmds::Variable<kmds::TCoord>* varVolumes;
            Kokkos::UnorderedMap<kmds::TCellID, ManifoldDetection_configStorage>* storedSolutions;

            ManifoldDetection_callForOneNode_2D(const kmds::GrowingView<kmds::TCellID>* selection_,
                                                kmds::Mesh* mesh_,
                                                const kmds::Connectivity* c_N2C_,
                                                const FracPres* fp_,
                                                MaterialAssignment* ma_,
                                                const kmds::Variable<kmds::TCoord>* varVolumes_,
                                                Kokkos::UnorderedMap<kmds::TCellID, ManifoldDetection_configStorage>* storedSolutions_
            )
                    : selection(selection_)
                    , mesh(mesh_)
                    , c_N2C(c_N2C_)
                    , fp(fp_)
                    , ma(ma_)
                    , varVolumes(varVolumes_)
                    , storedSolutions(storedSolutions_)
            {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(int i) const
            {
                int nid = selection->get(i);

                // form the ball of cells adjacent to nid
                Kokkos::View<kmds::TCellID *> cellIDs_tmp;
                c_N2C->get(nid, cellIDs_tmp);

                std::vector<kmds::TCellID> cellIDs(cellIDs_tmp.size());
                for(int icell=0; icell<cellIDs_tmp.size(); icell++) {
                    cellIDs[icell] = cellIDs_tmp[icell];
                }

                std::vector<int> cells2mat(cellIDs_tmp.size());

                // get a solution to the non-manifold
                ManifoldDetection_SolveNonManifoldAtOneNode_2D(mesh,
                                                               nid,
                                                               cellIDs,
                                                               fp,
                                                               varVolumes,
                                                               cells2mat,
                                                               storedSolutions);

                // reassign the cells
                for(int icell=0; icell<cellIDs_tmp.size(); icell++) {
                    ma->setMaterial(cells2mat[icell], cellIDs[icell]);
                }

            }

        };

        Kokkos::parallel_for(ASelection->getNbElems(), ManifoldDetection_callForOneNode_2D(ASelection, AMesh, c_N2C, Afp, Ama, AVarVolumes, AStoredSolutions));
    }
/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
