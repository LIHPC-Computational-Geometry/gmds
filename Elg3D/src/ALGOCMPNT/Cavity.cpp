/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    Cavity.cpp
 *  \author  legoff
 *  \date    05/01/2019
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/Cavity.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <set>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <Kokkos_Core.hpp>
#include <Kokkos_UnorderedMap.hpp>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <KM/IO/VTKWriter.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include <ELG3D/DATACMPNT/MaterialAssignment.h>
/*----------------------------------------------------------------------------*/
namespace elg3d {

/*----------------------------------------------------------------------------*/
//    struct InterfaceNodesPos_computePos_onenode_2D {
//
//        const kmds::GrowingView<kmds::TCellID> *nodeIDsAccessor;
//
//        kmds::Mesh *mesh;
//        const kmds::Connectivity *c_N2F;
//        const kmds::Variable<gmds::math::Point> *varPCrossPt;
//        const kmds::Variable<gmds::math::Vector> *varPCrossVec;
//
//        kmds::Variable<gmds::math::Point> *varNewPos;
//
//
//        InterfaceNodesPos_computePos_onenode_2D(const kmds::GrowingView<kmds::TCellID> *ANodeIDsAccessor_,
//                                                kmds::Mesh *AMesh_,
//                                                const kmds::Connectivity *c_N2F_,
//                                                const kmds::Variable<gmds::math::Point> *AVarPCrossPt_,
//                                                const kmds::Variable<gmds::math::Vector> *AVarPCrossVec_,
//                                                kmds::Variable<gmds::math::Point> *AVarNewPos
//        )
//                : nodeIDsAccessor(ANodeIDsAccessor_)
//                , mesh(AMesh_)
//                , c_N2F(c_N2F_)
//                , varPCrossPt(AVarPCrossPt_)
//                , varPCrossVec(AVarPCrossVec_)
//                , varNewPos(AVarNewPos) {
//        }
//
//        KOKKOS_INLINE_FUNCTION
//        void
//        operator()(int i) const {
//            int nid = nodeIDsAccessor->get(i);
//
//            kmds::TCoord coords[3];
//            mesh->getNodeLocation(nid, coords[0], coords[1], coords[2]);
//            gmds::math::Point oldPos(coords[0], coords[1], coords[2]);
//
//            Kokkos::View<kmds::TCellID*> faces;
//            c_N2F->get(nid, faces);
//
//
//            gmds::math::Point newPos(0.,0.,0.);
//            int nbplanes = 0;
//            for(int i_f=0; i_f<faces.size(); i_f++) {
//                gmds::math::Point pcrosspt = (*varPCrossPt)[faces(i_f)];
//                gmds::math::Vector pcrossvec = (*varPCrossVec)[faces(i_f)];
//
//                if(pcrosspt != gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF)) {
//
//                    gmds::math::Vector pcrossvec_orth(-pcrossvec.Y(), pcrossvec.X(), 0.);
//
//                    gmds::math::Line pl(pcrosspt, pcrossvec_orth);
//                    gmds::math::Point proj = pl.project(oldPos);
//                    newPos = newPos + proj;
//                    nbplanes++;
//
//                }
//            }
//
//            newPos = (1./nbplanes) * newPos;
//
//            (*varNewPos)[nid] = newPos;
//
////            exit(-1);
//        }
//    };

/*----------------------------------------------------------------------------*/
    void
    cavity_pillow_create_2D(const kmds::GrowingView<kmds::TCellID>* ANodeIDs,
                            kmds::Mesh* AMesh_origin,
                            kmds::Mesh* AMesh_cavity,
                            const kmds::Connectivity* c_N2C,
                            const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                            kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation_cavity,
                            const kmds::Variable<bool>* AVarIsMarkedNode,
                            const kmds::Variable<bool>* AVarIsMarkedEdge,
                            const kmds::Variable<bool>* AVarIsMarkedEdge_withBoundary,
                            const kmds::Variable<bool>* AVarIsMarkedCell,
                            kmds::Variable<bool>* AVarIsMarkedNode_cavity,
                            kmds::Variable<bool>* AVarIsMarkedEdge_cavity,
                            kmds::Variable<bool>* AVarIsMarkedEdge_withBoundary_cavity,
                            kmds::Variable<bool>* AVarIsMarkedCell_cavity,
                            kmds::Variable<kmds::TCellID>* AVarCavity2originN,
                            kmds::Variable<kmds::TCellID>* AVarOrigin2cavityN,
                            kmds::Variable<kmds::TCellID>* AVarCavity2originC,
                            kmds::Variable<kmds::TCellID>* AVarOrigin2cavityC,
                            kmds::GrowingView<kmds::TCellID>* ANodes_boundary_cavity)
    {
        std::cout<<"cavity_pillow_create_2D "<<std::endl;

        // copy the nodes
        Kokkos::parallel_for(ANodeIDs->getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID nid = ANodeIDs->get(i);

                                 Kokkos::View<kmds::TCellID *> cids;
                                 c_N2C->get(nid, cids);

                                 for (int i_c = 0; i_c < cids.size(); i_c++) {

                                     if ((*AVarIsMarkedCell)[cids[i_c]]) {
                                         kmds::TCellID nid_cavity = AMesh_cavity->newNode(AMesh_origin->getNodeLocation(nid));
                                         (*AVarCavity2originN)[nid_cavity] = nid;
                                         (*AVarOrigin2cavityN)[nid] = nid_cavity;

                                         (*AVarNodeGeomAssociation_cavity)[nid_cavity] = (*AVarNodeGeomAssociation)[nid];
                                         (*AVarIsMarkedNode_cavity)[nid_cavity] = (*AVarIsMarkedNode)[nid];
                                         break;
                                     }
                                 }
                             });

        // copy the cells
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", AMesh_origin->getNbFaces());
        AMesh_origin->getFaceIDs(&cellIDs);

        Kokkos::parallel_for(cellIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID cid = cellIDs.get(i);

                                 if((*AVarIsMarkedCell)[cid]) {

                                     kmds::Face c = AMesh_origin->getFace(cid);

                                     Kokkos::View<kmds::TCellID *> nids;
                                     c.nodeIds(nids);


                                     kmds::TCellID cid_cavity = AMesh_cavity->newQuad((*AVarOrigin2cavityN)[nids[0]],
                                                                                      (*AVarOrigin2cavityN)[nids[1]],
                                                                                      (*AVarOrigin2cavityN)[nids[2]],
                                                                                      (*AVarOrigin2cavityN)[nids[3]]);

                                     (*AVarCavity2originC)[cid_cavity] = cid;
                                     (*AVarOrigin2cavityC)[cid] = cid_cavity;
                                     (*AVarIsMarkedCell_cavity)[cid_cavity] = true;
                                 }
                             });

        kmds::VTKWriter<kmds::Mesh> w(*AMesh_cavity);
        w.write("cavity_NandC", kmds::F);

        // create the edges and connectivities
        kmds::ConnectivityHelper ch(AMesh_cavity);

        kmds::Connectivity* c_E2F = AMesh_cavity->createConnectivity(kmds::E2F);
        ch.buildEandE2F_2D_variant_0();
        kmds::Connectivity* c_N2E = AMesh_cavity->createConnectivity(kmds::N2E);
        ch.buildN2E();


        kmds::GrowingView<kmds::TCellID> edgeIDs("EDGES", AMesh_origin->getNbEdges());
        AMesh_origin->getEdgeIDs(&edgeIDs);

        Kokkos::UnorderedMap<uint64_t, kmds::TCellID > e2n_single(AMesh_cavity->getNbEdges());

        const int nbNodes_cavity = AMesh_cavity->getNbNodes();


        Kokkos::parallel_for(edgeIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID eid = edgeIDs.get(i);

                                 kmds::Edge e = AMesh_origin->getEdge(eid);

                                 Kokkos::View<kmds::TCellID *> nids;
                                 e.nodeIds(nids);

                                 // keep only edges which nodes are in the cavity
                                 if(((*AVarOrigin2cavityN)[nids[0]] != kmds::NullID) && ((*AVarOrigin2cavityN)[nids[1]] != kmds::NullID)) {

                                     if ((nids[0] == 22 && nids[1] == 28)
                                         || (nids[1] == 22 && nids[0] == 28)) {
                                         std::cout << "edge_origin " << eid << " marked " << (*AVarIsMarkedEdge)[eid]
                                                   << " boundary " << (*AVarIsMarkedEdge_withBoundary)[eid]
                                                   << std::endl;
                                     }

                                     // unique id for the edge, defined as its expression in base nbNodes_cavity
                                     uint64_t eid_unique;
                                     if ((*AVarOrigin2cavityN)[nids[0]] > (*AVarOrigin2cavityN)[nids[1]]) {
                                         eid_unique = (*AVarOrigin2cavityN)[nids[0]] +
                                                      (*AVarOrigin2cavityN)[nids[1]] * nbNodes_cavity;
                                     } else {
                                         eid_unique = (*AVarOrigin2cavityN)[nids[1]] +
                                                      (*AVarOrigin2cavityN)[nids[0]] * nbNodes_cavity;
                                     }

                                     if (eid == 43) {
                                         std::cout << "edge_43 " << nids[0] << " " << nids[1] << " cavity "
                                                   << (*AVarOrigin2cavityN)[nids[0]] << " "
                                                   << (*AVarOrigin2cavityN)[nids[1]] << " " << eid_unique << std::endl;
                                     }

                                     //Kokkos::UnorderedMapInsertResult res = e2fid_single.insert(eid_unique, fid);
                                     Kokkos::UnorderedMapInsertResult res = e2n_single.insert(eid_unique, eid);

                                 }
                             });


        kmds::GrowingView<kmds::TCellID> edgeIDs_cavity("EDGES_CAVITY", AMesh_cavity->getNbEdges());
        AMesh_cavity->getEdgeIDs(&edgeIDs_cavity);

        Kokkos::parallel_for(edgeIDs_cavity.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID eid_cavity = edgeIDs_cavity.get(i);

                                 kmds::Edge e = AMesh_cavity->getEdge(eid_cavity);

                                 Kokkos::View<kmds::TCellID *> nids;
                                 e.nodeIds(nids);

//                                 std::cout<<"edge_cavity "<<nids[0]<<" "<<nids[1]<<std::endl;


//                                 if((nids[0] == 16 && nids[1] == 18)
//                                 || (nids[1] == 16 && nids[0] == 18)) {
//                                     std::cout<<"edge_cavity "<<eid_cavity<<std::endl;
//                                 }

                                 // unique id for the edge, defined as its expression in base nb_nodes
                                 uint64_t eid_unique;
                                 if (nids[0] > nids[1]) {
                                     eid_unique = nids[0] + nids[1] * nbNodes_cavity;
                                 } else {
                                     eid_unique = nids[1] + nids[0] * nbNodes_cavity;
                                 }

                                 int index = e2n_single.find(eid_unique);
                                 if(e2n_single.valid_at(index)) {
                                     kmds::TCellID eid_origin = e2n_single.value_at(index);

                                     (*AVarIsMarkedEdge_cavity)[eid_cavity] = (*AVarIsMarkedEdge)[eid_origin];
                                     (*AVarIsMarkedEdge_withBoundary_cavity)[eid_cavity] = (*AVarIsMarkedEdge_withBoundary)[eid_origin];

//                                     if(eid_cavity == 26) {
//                                         std::cout<<"reverse_edge "<<eid_cavity<<" "<<eid_origin<<std::endl;
//                                     }

                                 } else {
                                     std::cerr<<"cavity_pillow_create_2D : pb with an edge "<<std::endl;
                                 }

                             });


        // identify and fill the container of boundary nodes
        ANodes_boundary_cavity->reserve(AMesh_cavity->getNbNodes());
        Kokkos::UnorderedMap<kmds::TCellID, void> nodes_boundary_set(AMesh_cavity->getNbNodes());

        Kokkos::parallel_for(edgeIDs_cavity.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID eid_cavity = edgeIDs_cavity.get(i);

                                 kmds::Edge e = AMesh_cavity->getEdge(eid_cavity);

                                 Kokkos::View<kmds::TCellID *> cids;
                                 c_E2F->get(eid_cavity, cids);

                                 if(cids.size() == 1) {

                                     Kokkos::View<kmds::TCellID *> nids;
                                     e.nodeIds(nids);

                                     Kokkos::UnorderedMapInsertResult res0 = nodes_boundary_set.insert(nids[0]);
                                     bool success0 = res0.success();

                                     if (success0) {
                                         ANodes_boundary_cavity->push_back(nids[0]);
                                     }

                                     Kokkos::UnorderedMapInsertResult res1 = nodes_boundary_set.insert(nids[1]);
                                     bool success1 = res1.success();

                                     if (success1) {
                                         ANodes_boundary_cavity->push_back(nids[1]);
                                     }
                                 }
                             });


        // clean-up
//        AMesh->deleteVariable(kmds::KMDS_FACE, "tmp_pcrosspt");
//        AMesh->deleteVariable(kmds::KMDS_FACE, "tmp_pcrossvec");
    }

    /*----------------------------------------------------------------------------*/
    void
    cavity_pillow_insert_2D(kmds::Mesh* AMesh_origin,
                            kmds::Mesh* AMesh_cavity,
                            kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                            const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation_cavity,
                            kmds::Variable<kmds::TCellID>* AVarCavity2originN,
                            const kmds::Variable<kmds::TCellID>* AVarOrigin2cavityN,
                            const kmds::Variable<kmds::TCellID>* AVarCavity2originC,
                            const kmds::Variable<kmds::TCellID>* AVarOrigin2cavityC,
                            const int AImat,
                            elg3d::MaterialAssignment* Ama)
    {

        // insert the nodes

        kmds::GrowingView<kmds::TCellID> nodeIDs("NODES_CAVITY", AMesh_cavity->getNbNodes());
        AMesh_cavity->getNodeIDs(&nodeIDs);

        Kokkos::parallel_for(nodeIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 kmds::TCellID nid_cavity = nodeIDs.get(i);

                                 if((*AVarCavity2originN)[nid_cavity] == kmds::NullID) {

                                     kmds::TCellID nid_origin = AMesh_origin->newNode(AMesh_cavity->getNodeLocation(nid_cavity));
                                     (*AVarCavity2originN)[nid_cavity] = nid_origin;

                                     (*AVarNodeGeomAssociation)[nid_origin] = (*AVarNodeGeomAssociation_cavity)[nid_cavity];

                                 } else {
                                     AMesh_origin->setNodeLocation((*AVarCavity2originN)[nid_cavity], AMesh_cavity->getNodeLocation(nid_cavity));
                                 }
                             });

        // insert the cells

        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS_CAVITY", AMesh_cavity->getNbFaces());
        AMesh_cavity->getFaceIDs(&cellIDs);

        Kokkos::parallel_for(cellIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 kmds::TCellID cid_cavity = cellIDs.get(i);

                                 kmds::Face c = AMesh_cavity->getFace(cid_cavity);

                                 Kokkos::View<kmds::TCellID *> nids;
                                 c.nodeIds(nids);

                                 for(int i_n=0; i_n<nids.size(); i_n++) {
                                     nids[i_n] = (*AVarCavity2originN)[nids[i_n]];
                                 }

                                 if((*AVarCavity2originC)[cid_cavity] == kmds::NullID) {

                                     kmds::TCellID cid_origin = AMesh_origin->newQuad(nids[0], nids[1], nids[2], nids[3]);
//                                     (*AVarCavity2originC)[cid_cavity] = cid_origin;
                                     Ama->setMaterial(AImat, cid_origin);

                                 } else {
                                     kmds::Face c_origin = AMesh_origin->getFace((*AVarCavity2originC)[cid_cavity]);
                                     c_origin.setNodes(nids);
                                 }
                             });
    }

    /*----------------------------------------------------------------------------*/
    int
    cavity_computeNbNewNodes_xD(kmds::Mesh* AMesh_cavity,
                                const kmds::Variable<kmds::TCellID>* AVarCavity2originN)

    {
        kmds::GrowingView<kmds::TCellID> nodeIDs("NODES_CAVITY", AMesh_cavity->getNbNodes());
        AMesh_cavity->getNodeIDs(&nodeIDs);

        int nbNewNodes = 0;

        Kokkos::parallel_reduce(nodeIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i, int& sum) {

                                 kmds::TCellID nid_cavity = nodeIDs.get(i);

                                 if ((*AVarCavity2originN)[nid_cavity] == kmds::NullID) {
                                    sum++;
                                 }
                             },
                             nbNewNodes);

        return nbNewNodes;
    }

    /*----------------------------------------------------------------------------*/
    int
    cavity_computeNbNewCells_2D(kmds::Mesh* AMesh_cavity,
                                const kmds::Variable<kmds::TCellID>* AVarCavity2originC)

    {
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS_CAVITY", AMesh_cavity->getNbFaces());
        AMesh_cavity->getFaceIDs(&cellIDs);

        int nbNewCells = 0;

        Kokkos::parallel_reduce(cellIDs.getNbElems(),
                                KOKKOS_LAMBDA(const int i, int& sum) {

                                    kmds::TCellID cid_cavity = cellIDs.get(i);

                                    if ((*AVarCavity2originC)[cid_cavity] == kmds::NullID) {
                                        sum++;
                                    }
                                },
                                nbNewCells);

        return nbNewCells;
    }

    /*----------------------------------------------------------------------------*/
    int
    cavity_computeNbNewCells_3D(kmds::Mesh* AMesh_cavity,
                                const kmds::Variable<kmds::TCellID>* AVarCavity2originC)

    {
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS_CAVITY", AMesh_cavity->getNbRegions());
        AMesh_cavity->getRegionIDs(&cellIDs);

        int nbNewCells = 0;

        Kokkos::parallel_reduce(cellIDs.getNbElems(),
                                KOKKOS_LAMBDA(const int i, int& sum) {

                                    kmds::TCellID cid_cavity = cellIDs.get(i);

                                    if ((*AVarCavity2originC)[cid_cavity] == kmds::NullID) {
                                        sum++;
                                    }
                                },
                                nbNewCells);

        return nbNewCells;
    }

    /*----------------------------------------------------------------------------*/
    void
    cavity_pillow_create_3D(const kmds::GrowingView<kmds::TCellID>* ANodeIDs,
                            kmds::Mesh* AMesh_origin,
                            kmds::Mesh* AMesh_cavity,
                            const kmds::Connectivity* c_N2C,
                            const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                            kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation_cavity,
                            const kmds::Variable<bool>* AVarIsMarkedNode,
                            const kmds::Variable<bool>* AVarIsMarkedA,
                            const kmds::Variable<bool>* AVarIsMarkedA_withBoundary,
                            const kmds::Variable<bool>* AVarIsMarkedCell,
                            kmds::Variable<bool>* AVarIsMarkedNode_cavity,
                            kmds::Variable<bool>* AVarIsMarkedA_cavity,
                            kmds::Variable<bool>* AVarIsMarkedA_withBoundary_cavity,
                            kmds::Variable<bool>* AVarIsMarkedCell_cavity,
                            kmds::Variable<kmds::TCellID>* AVarCavity2originN,
                            kmds::Variable<kmds::TCellID>* AVarOrigin2cavityN,
                            kmds::Variable<kmds::TCellID>* AVarCavity2originC,
                            kmds::Variable<kmds::TCellID>* AVarOrigin2cavityC,
                            kmds::GrowingView<kmds::TCellID>* ANodes_boundary_cavity)
    {
        std::cout<<"cavity_pillow_create_3D "<<std::endl;

        // copy the nodes
        Kokkos::parallel_for(ANodeIDs->getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID nid = ANodeIDs->get(i);

                                 Kokkos::View<kmds::TCellID *> cids;
                                 c_N2C->get(nid, cids);

                                 for (int i_c = 0; i_c < cids.size(); i_c++) {

                                     if ((*AVarIsMarkedCell)[cids[i_c]]) {
                                         kmds::TCellID nid_cavity = AMesh_cavity->newNode(AMesh_origin->getNodeLocation(nid));
                                         (*AVarCavity2originN)[nid_cavity] = nid;
                                         (*AVarOrigin2cavityN)[nid] = nid_cavity;

                                         (*AVarNodeGeomAssociation_cavity)[nid_cavity] = (*AVarNodeGeomAssociation)[nid];
                                         (*AVarIsMarkedNode_cavity)[nid_cavity] = (*AVarIsMarkedNode)[nid];
                                         break;
                                     }
                                 }
                             });

        // copy the cells
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", AMesh_origin->getNbRegions());
        AMesh_origin->getRegionIDs(&cellIDs);

        Kokkos::parallel_for(cellIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID cid = cellIDs.get(i);

                                 if((*AVarIsMarkedCell)[cid]) {

                                     kmds::Region c = AMesh_origin->getRegion(cid);

                                     Kokkos::View<kmds::TCellID *> nids;
                                     c.nodeIds(nids);


                                     kmds::TCellID cid_cavity = AMesh_cavity->newHexahedron((*AVarOrigin2cavityN)[nids[0]],
                                                                                            (*AVarOrigin2cavityN)[nids[1]],
                                                                                            (*AVarOrigin2cavityN)[nids[2]],
                                                                                            (*AVarOrigin2cavityN)[nids[3]],
                                                                                            (*AVarOrigin2cavityN)[nids[4]],
                                                                                            (*AVarOrigin2cavityN)[nids[5]],
                                                                                            (*AVarOrigin2cavityN)[nids[6]],
                                                                                            (*AVarOrigin2cavityN)[nids[7]]);

                                     (*AVarCavity2originC)[cid_cavity] = cid;
                                     (*AVarOrigin2cavityC)[cid] = cid_cavity;
                                     (*AVarIsMarkedCell_cavity)[cid_cavity] = true;
                                 }
                             });

        kmds::VTKWriter<kmds::Mesh> w(*AMesh_cavity);
        w.write("cavity_NandC", kmds::R);

        // create the edges and connectivities
        kmds::ConnectivityHelper ch(AMesh_cavity);

        kmds::Connectivity* c_F2R = AMesh_cavity->createConnectivity(kmds::F2R);
        ch.buildFandF2R_variant_0();
        kmds::Connectivity* c_N2F = AMesh_cavity->createConnectivity(kmds::N2F);
        ch.buildN2F();


        kmds::GrowingView<kmds::TCellID> faceIDs("FACES_cavity", AMesh_origin->getNbFaces());
        AMesh_origin->getFaceIDs(&faceIDs);

        Kokkos::UnorderedMap<uint64_t, kmds::TCellID > a2n_single(AMesh_cavity->getNbFaces());

        const int nbNodes_cavity = AMesh_cavity->getNbNodes();


        Kokkos::parallel_for(faceIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID aid = faceIDs.get(i);

                                 kmds::Face a = AMesh_origin->getFace(aid);

                                 Kokkos::View<kmds::TCellID *> nids;
                                 a.nodeIds(nids);

                                 // keep only edges which nodes are in the cavity
                                 if(((*AVarOrigin2cavityN)[nids[0]] != kmds::NullID)
                                    && ((*AVarOrigin2cavityN)[nids[1]] != kmds::NullID)
                                    && ((*AVarOrigin2cavityN)[nids[2]] != kmds::NullID)
                                    && ((*AVarOrigin2cavityN)[nids[3]] != kmds::NullID)) {

//                                     if ((nids[0] == 22 && nids[1] == 28)
//                                         || (nids[1] == 22 && nids[0] == 28)) {
//                                         std::cout << "edge_origin " << eid << " marked " << (*AVarIsMarkedEdge)[eid]
//                                                   << " boundary " << (*AVarIsMarkedEdge_withBoundary)[eid]
//                                                   << std::endl;
//                                     }

                                     // unique id for the face, defined as its expression in base nbNodes_cavity
                                     // WARNING : we consider only the first three nodes of the face

                                     uint64_t aid_unique;

                                     kmds::TCellID order[3];
                                     order[0] = (*AVarOrigin2cavityN)[nids[0]];
                                     order[1] = (*AVarOrigin2cavityN)[nids[1]];
                                     order[2] = (*AVarOrigin2cavityN)[nids[2]];

                                     if(order[0] < order[1]) {
                                         std::swap(order[0], order[1]);
                                     }
                                     if(order[1] < order[2]) {
                                         std::swap(order[1], order[2]);
                                     }
                                     if(order[0] < order[1]) {
                                         std::swap(order[0], order[1]);
                                     }

                                     aid_unique = order[0] +
                                                  order[1] * nbNodes_cavity +
                                                  order[2] * std::pow(nbNodes_cavity, 2);


//                                     if (eid == 43) {
//                                         std::cout << "edge_43 " << nids[0] << " " << nids[1] << " cavity "
//                                                   << (*AVarOrigin2cavityN)[nids[0]] << " "
//                                                   << (*AVarOrigin2cavityN)[nids[1]] << " " << eid_unique << std::endl;
//                                     }

                                     //Kokkos::UnorderedMapInsertResult res = e2fid_single.insert(eid_unique, fid);
                                     Kokkos::UnorderedMapInsertResult res = a2n_single.insert(aid_unique, aid);

                                 }
                             });


        kmds::GrowingView<kmds::TCellID> faceIDs_cavity("FACES_CAVITY", AMesh_cavity->getNbFaces());
        AMesh_cavity->getFaceIDs(&faceIDs_cavity);

        Kokkos::parallel_for(faceIDs_cavity.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID aid_cavity = faceIDs_cavity.get(i);

                                 kmds::Face a = AMesh_cavity->getFace(aid_cavity);

                                 Kokkos::View<kmds::TCellID *> nids;
                                 a.nodeIds(nids);

//                                 std::cout<<"edge_cavity "<<nids[0]<<" "<<nids[1]<<std::endl;


//                                 if((nids[0] == 16 && nids[1] == 18)
//                                 || (nids[1] == 16 && nids[0] == 18)) {
//                                     std::cout<<"edge_cavity "<<eid_cavity<<std::endl;
//                                 }

                                 // unique id for the face, defined as its expression in base nb_nodes
                                 uint64_t aid_unique;

                                 kmds::TCellID order[3];
                                 order[0] = nids[0];
                                 order[1] = nids[1];
                                 order[2] = nids[2];

                                 if(order[0] < order[1]) {
                                     std::swap(order[0], order[1]);
                                 }
                                 if(order[1] < order[2]) {
                                     std::swap(order[1], order[2]);
                                 }
                                 if(order[0] < order[1]) {
                                     std::swap(order[0], order[1]);
                                 }

                                 aid_unique = order[0] +
                                              order[1] * nbNodes_cavity +
                                              order[2] * std::pow(nbNodes_cavity, 2);


                                 int index = a2n_single.find(aid_unique);
                                 if(a2n_single.valid_at(index)) {
                                     kmds::TCellID aid_origin = a2n_single.value_at(index);

                                     (*AVarIsMarkedA_cavity)[aid_cavity] = (*AVarIsMarkedA)[aid_origin];
                                     (*AVarIsMarkedA_withBoundary_cavity)[aid_cavity] = (*AVarIsMarkedA_withBoundary)[aid_origin];

//                                     if(eid_cavity == 26) {
//                                         std::cout<<"reverse_edge "<<eid_cavity<<" "<<eid_origin<<std::endl;
//                                     }

                                 } else {
                                     std::cerr<<"cavity_pillow_create_3D : pb with a face "<<std::endl;
                                 }

                             });


        // identify and fill the container of boundary nodes
        ANodes_boundary_cavity->reserve(AMesh_cavity->getNbNodes());
        Kokkos::UnorderedMap<kmds::TCellID, void> nodes_boundary_set(AMesh_cavity->getNbNodes());

        Kokkos::parallel_for(faceIDs_cavity.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID aid_cavity = faceIDs_cavity.get(i);

                                 kmds::Face a = AMesh_cavity->getFace(aid_cavity);

                                 Kokkos::View<kmds::TCellID *> cids;
                                 c_F2R->get(aid_cavity, cids);

                                 if(cids.size() == 1) {

                                     Kokkos::View<kmds::TCellID *> nids;
                                     a.nodeIds(nids);

                                     Kokkos::UnorderedMapInsertResult res0 = nodes_boundary_set.insert(nids[0]);
                                     bool success0 = res0.success();

                                     if (success0) {
                                         ANodes_boundary_cavity->push_back(nids[0]);
                                     }

                                     Kokkos::UnorderedMapInsertResult res1 = nodes_boundary_set.insert(nids[1]);
                                     bool success1 = res1.success();

                                     if (success1) {
                                         ANodes_boundary_cavity->push_back(nids[1]);
                                     }

                                     Kokkos::UnorderedMapInsertResult res2 = nodes_boundary_set.insert(nids[2]);
                                     bool success2 = res2.success();

                                     if (success2) {
                                         ANodes_boundary_cavity->push_back(nids[2]);
                                     }

                                     Kokkos::UnorderedMapInsertResult res3 = nodes_boundary_set.insert(nids[3]);
                                     bool success3 = res3.success();

                                     if (success3) {
                                         ANodes_boundary_cavity->push_back(nids[3]);
                                     }
                                 }
                             });


        // clean-up
//        AMesh->deleteVariable(kmds::KMDS_FACE, "tmp_pcrosspt");
//        AMesh->deleteVariable(kmds::KMDS_FACE, "tmp_pcrossvec");
    }

    /*----------------------------------------------------------------------------*/
    void
    cavity_pillow_insert_3D(kmds::Mesh* AMesh_origin,
                            kmds::Mesh* AMesh_cavity,
                            kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                            const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation_cavity,
                            kmds::Variable<kmds::TCellID>* AVarCavity2originN,
                            const kmds::Variable<kmds::TCellID>* AVarOrigin2cavityN,
                            const kmds::Variable<kmds::TCellID>* AVarCavity2originC,
                            const kmds::Variable<kmds::TCellID>* AVarOrigin2cavityC,
                            const int AImat,
                            elg3d::MaterialAssignment* Ama)
    {

        // insert the nodes

        kmds::GrowingView<kmds::TCellID> nodeIDs("NODES_CAVITY", AMesh_cavity->getNbNodes());
        AMesh_cavity->getNodeIDs(&nodeIDs);

        Kokkos::parallel_for(nodeIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 kmds::TCellID nid_cavity = nodeIDs.get(i);

                                 if((*AVarCavity2originN)[nid_cavity] == kmds::NullID) {

                                     kmds::TCellID nid_origin = AMesh_origin->newNode(AMesh_cavity->getNodeLocation(nid_cavity));
                                     (*AVarCavity2originN)[nid_cavity] = nid_origin;

                                     // TODO : not valid in unstructured cases
                                     (*AVarNodeGeomAssociation)[nid_origin] = (*AVarNodeGeomAssociation_cavity)[nid_cavity];

                                 } else {
                                     AMesh_origin->setNodeLocation((*AVarCavity2originN)[nid_cavity], AMesh_cavity->getNodeLocation(nid_cavity));
                                 }
                             });

        // insert the cells

        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS_CAVITY", AMesh_cavity->getNbRegions());
        AMesh_cavity->getRegionIDs(&cellIDs);

        Kokkos::parallel_for(cellIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 kmds::TCellID cid_cavity = cellIDs.get(i);

                                 kmds::Region c = AMesh_cavity->getRegion(cid_cavity);

                                 Kokkos::View<kmds::TCellID *> nids;
                                 c.nodeIds(nids);

                                 for(int i_n=0; i_n<nids.size(); i_n++) {
                                     nids[i_n] = (*AVarCavity2originN)[nids[i_n]];
                                 }

                                 if((*AVarCavity2originC)[cid_cavity] == kmds::NullID) {

                                     kmds::TCellID cid_origin = AMesh_origin->newHexahedron(nids[0], nids[1], nids[2], nids[3],
                                                                                            nids[4], nids[5], nids[6], nids[7]);
//                                     (*AVarCavity2originC)[cid_cavity] = cid_origin;
                                     Ama->setMaterial(AImat, cid_origin);

                                 } else {
                                     kmds::Region c_origin = AMesh_origin->getRegion((*AVarCavity2originC)[cid_cavity]);
                                     c_origin.setNodes(nids);
                                 }
                             });
    }

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
