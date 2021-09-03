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
/** \file    ExtractGeomModel.cpp
 *  \author  legoff
 *  \date    01/18/2019
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/ExtractGeomModel.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <list>
#include <set>
#include <utility>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <Kokkos_UnorderedMap.hpp>
#include <gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
//#include "GMDS/CAD/GeomEntity.h"

#include <KM/IO/VTKWriter.h>

/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/MaterialInterfaces.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {




//    /*----------------------------------------------------------------------------*/
//    struct Pillow_createNodes {
//        const kmds::GrowingView<kmds::TCellID> *selection;
//        kmds::Mesh *mesh;
//        const kmds::Connectivity* c_N2C;
//        const elg3d::MaterialAssignment* ma;
//        kmds::Variable<std::uintptr_t> *varNodeGeomAssociation;
//        std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> >* node2newNodes;
//
//
//        Pillow_createNodes(
//
//                const kmds::GrowingView<kmds::TCellID> *selection_,
//                kmds::Mesh *mesh_,
//                const kmds::Connectivity* c_N2C_,
//                const elg3d::MaterialAssignment* ma_,
//                kmds::Variable<std::uintptr_t> *AVarNodeGeomAssociation_,
//                std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> >* node2newNodes_
//        )
//                : selection(selection_)
//                , mesh(mesh_)
//                , c_N2C(c_N2C_)
//                , ma(ma_)
//                , varNodeGeomAssociation(AVarNodeGeomAssociation_)
//                , node2newNodes(node2newNodes_)
//        {
//        }
//
//        KOKKOS_INLINE_FUNCTION
//        void
//        operator()(int i) const {
//            int nid = selection->get(i);
//
//            Kokkos::View<kmds::TCellID*> cells;
//            c_N2C->get(nid, cells);
//
//            std::set<int> materials;
//            for(int i_c=0; i_c<cells.size(); i_c++) {
//                materials.insert(ma->getMaterial(cells(i_c)));
//            }
//
//
//            kmds::TCellID id0 = mesh->addNodes(materials.size());
//            kmds::TCoord xyz[3];
//            mesh->getNodeLocation(nid, xyz[0], xyz[1], xyz[2]);
//
//            int offset = 0;
//            for(auto mat: materials) {
//
//                // compute the new node position
//                gmds::math::Point pt_old(xyz[0], xyz[1], xyz[2]);
//                gmds::math::Point pt_new(0., 0., 0.);
//                int nbContrib = 0;
//                for(int i_c=0; i_c<cells.size(); i_c++) {
//                    if(ma->getMaterial(cells(i_c)) == mat) {
//                        gmds::math::Point pt_midpoint = mesh->getFace(cells(i_c)).midpoint();
//                        gmds::math::Point pt_tmp(pt_old + (1/5.) * gmds::math::Vector(pt_old, pt_midpoint));
//                        pt_new = pt_new + pt_tmp;
//                        nbContrib++;
//                    }
//                }
//
//                pt_new = (1./(double) nbContrib) * pt_new;
//                if((*varNodeGeomAssociation)[nid] != reinterpret_cast<std::uintptr_t>(nullptr)) {
//                    reinterpret_cast<gmds::cad::GeomEntity *>((*varNodeGeomAssociation)[nid])->project(pt_new);
//                }
//
//
//                mesh->setNodeLocation(id0 + offset, pt_new.X(), pt_new.Y(), pt_new.Z());
//
////                mesh->setNodeLocation(id0 + offset, xyz[0], xyz[1], xyz[2]);
//
//                (*node2newNodes)[mat].insert(nid, id0 + offset);
//                (*varNodeGeomAssociation)[id0 + offset] = (*varNodeGeomAssociation)[nid];
//                offset++;
//            }
//
//        }
//
//    };
//
//    /*----------------------------------------------------------------------------*/
//    struct Pillow_createTwoCells_2D {
//        const kmds::GrowingView<kmds::TCellID> *selection;
//        kmds::Mesh *mesh;
//        const kmds::Connectivity* c_F2C;
//        elg3d::MaterialAssignment* ma;
//        const std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> >* node2newNodes;
//
//
//        Pillow_createTwoCells_2D(
//                const kmds::GrowingView<kmds::TCellID> *selection_,
//                kmds::Mesh *mesh_,
//                const kmds::Connectivity* c_F2C_,
//                elg3d::MaterialAssignment* ma_,
//                std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> >* node2newNodes_
//        )
//                : selection(selection_)
//                , mesh(mesh_)
//                , c_F2C(c_F2C_)
//                , ma(ma_)
//                , node2newNodes(node2newNodes_)
//        {
//        }
//
//        KOKKOS_INLINE_FUNCTION
//        void
//        operator()(int i) const {
//            int fid = selection->get(i);
//
//            Kokkos::View<kmds::TCellID*> cells;
//            c_F2C->get(fid, cells);
//
//            kmds::Edge f = mesh->getEdge(fid);
//            Kokkos::View<kmds::TCellID*> nodes;
//            f.nodeIds(nodes);
//
//
//            // create two cells, which are always quads in 2D
//            kmds::TCellID cid0_new = mesh->addQuads(2);
//
//
//            // determine which new nodes will be used for which cell
//            int mat0 = ma->getMaterial(cells(0));
//            int mat1 = ma->getMaterial(cells(1));
//
//            kmds::Face c0 = mesh->getFace(cells(0));
//
//            bool isOutward0 = c0.isEdgeOrientedOutward(nodes);
//
//            kmds::TCellID newIDs_0[4];
//            kmds::TCellID newIDs_1[4];
//
//            if(isOutward0) {
//                newIDs_0[0] = nodes(0);
//                newIDs_0[1] = nodes(1);
//                newIDs_0[2] = (*node2newNodes)[mat0].value_at((*node2newNodes)[mat0].find(nodes(1)));
//                newIDs_0[3] = (*node2newNodes)[mat0].value_at((*node2newNodes)[mat0].find(nodes(0)));
//
//                newIDs_1[0] = nodes(1);
//                newIDs_1[1] = nodes(0);
//                newIDs_1[2] = (*node2newNodes)[mat1].value_at((*node2newNodes)[mat1].find(nodes(0)));
//                newIDs_1[3] = (*node2newNodes)[mat1].value_at((*node2newNodes)[mat1].find(nodes(1)));
//
//            } else {
//                newIDs_0[0] = nodes(0);
//                newIDs_0[1] = (*node2newNodes)[mat0].value_at((*node2newNodes)[mat0].find(nodes(0)));
//                newIDs_0[2] = (*node2newNodes)[mat0].value_at((*node2newNodes)[mat0].find(nodes(1)));
//                newIDs_0[3] = nodes(1);
//
//                newIDs_1[0] = nodes(0);
//                newIDs_1[1] = nodes(1);
//                newIDs_1[2] = (*node2newNodes)[mat1].value_at((*node2newNodes)[mat1].find(nodes(1)));
//                newIDs_1[3] = (*node2newNodes)[mat1].value_at((*node2newNodes)[mat1].find(nodes(0)));
//
//            }
//
//            kmds::Face c0_new = mesh->getFace(cid0_new);
//            kmds::Face c1_new = mesh->getFace(cid0_new + 1);
//            c0_new.setNodes(newIDs_0, 4);
//            c1_new.setNodes(newIDs_1, 4);
//
//            ma->setMaterial(mat0, cid0_new);
//            ma->setMaterial(mat1, cid0_new + 1);
//        }
//
//    };
//
//    /*----------------------------------------------------------------------------*/
//    struct Pillow_createTwoCells_3D {
//        const kmds::GrowingView<kmds::TCellID> *selection;
//        kmds::Mesh *mesh;
//        const kmds::Connectivity* c_F2C;
//        elg3d::MaterialAssignment* ma;
//        const std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> >* node2newNodes;
//
//
//        Pillow_createTwoCells_3D(
//                const kmds::GrowingView<kmds::TCellID> *selection_,
//                kmds::Mesh *mesh_,
//                const kmds::Connectivity* c_F2C_,
//                elg3d::MaterialAssignment* ma_,
//                std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> >* node2newNodes_
//        )
//                : selection(selection_)
//                , mesh(mesh_)
//                , c_F2C(c_F2C_)
//                , ma(ma_)
//                , node2newNodes(node2newNodes_)
//        {
//        }
//
//        KOKKOS_INLINE_FUNCTION
//        void
//        operator()(int i) const {
//            int fid = selection->get(i);
//
//            Kokkos::View<kmds::TCellID*> cells;
//            c_F2C->get(fid, cells);
//
//            kmds::Face f = mesh->getFace(fid);
//            Kokkos::View<kmds::TCellID*> nodes;
//            f.nodeIds(nodes);
//
//
//            // create two cells
//            kmds::TCellID cid0_new;
//            if(nodes.size() == 4) {
//                cid0_new = mesh->addHexahedra(2);
//            } else {
//                cid0_new = mesh->addPrism3s(2);
//            }
//
//
//            // determine which new nodes will be used for which cell
//            int mat0 = ma->getMaterial(cells(0));
//            int mat1 = ma->getMaterial(cells(1));
//
//            kmds::Region c0 = mesh->getRegion(cells(0));
//
//            bool isOutward0 = c0.isFaceOrientedOutward(nodes);
//
//            kmds::TCellID newIDs_0[8];
//            kmds::TCellID newIDs_1[8];
//
//            if(isOutward0) {
//                for(int i_n=0; i_n<nodes.size(); i_n++) {
//                    newIDs_0[i_n] = (*node2newNodes)[mat0].value_at((*node2newNodes)[mat0].find(nodes(i_n)));
//                    newIDs_0[i_n + nodes.size()] = nodes(i_n);
//                    newIDs_1[i_n] = nodes(i_n);
//                    newIDs_1[i_n + nodes.size()] = (*node2newNodes)[mat1].value_at((*node2newNodes)[mat1].find(nodes(i_n)));
//                }
//            } else {
//                for(int i_n=0; i_n<nodes.size(); i_n++) {
//                    newIDs_0[i_n] = nodes(i_n);
//                    newIDs_0[i_n + nodes.size()] = (*node2newNodes)[mat0].value_at((*node2newNodes)[mat0].find(nodes(i_n)));
//                    newIDs_1[i_n] = (*node2newNodes)[mat1].value_at((*node2newNodes)[mat1].find(nodes(i_n)));
//                    newIDs_1[i_n + nodes.size()] = nodes(i_n);
//                }
//            }
//
//            kmds::Region c0_new = mesh->getRegion(cid0_new);
//            kmds::Region c1_new = mesh->getRegion(cid0_new + 1);
//            c0_new.setNodes(newIDs_0, nodes.size() * 2);
//            c1_new.setNodes(newIDs_1, nodes.size() * 2);
//
//            ma->setMaterial(mat0, cid0_new);
//            ma->setMaterial(mat1, cid0_new + 1);
//        }
//
//    };
//
//    /*----------------------------------------------------------------------------*/
//    struct Pillow_replaceInCells_2D {
//        const kmds::GrowingView<kmds::TCellID> *selection;
//        kmds::Mesh *mesh;
//        const elg3d::MaterialAssignment* ma;
//        const kmds::Variable<bool>* varIsInterfaceNode;
//        const std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> >* node2newNodes;
//
//
//        Pillow_replaceInCells_2D(
//                const kmds::GrowingView<kmds::TCellID> *selection_,
//                kmds::Mesh *mesh_,
//                const elg3d::MaterialAssignment* ma_,
//                const kmds::Variable<bool>* varIsInterfaceNode_,
//                const std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> >* node2newNodes_
//        )
//                : selection(selection_)
//                , mesh(mesh_)
//                , ma(ma_)
//                , varIsInterfaceNode(varIsInterfaceNode_)
//                , node2newNodes(node2newNodes_)
//        {
//        }
//
//        KOKKOS_INLINE_FUNCTION
//        void
//        operator()(int i) const {
//            int cid = selection->get(i);
//
//            kmds::Face c = mesh->getFace(cid);
//            Kokkos::View<kmds::TCellID *> nodes;
//            c.nodeIds(nodes);
//
//            for (int i_n = 0; i_n < nodes.size(); i_n++) {
//
//                if((*varIsInterfaceNode)[nodes(i_n)]) {
//                    int mat = ma->getMaterial(cid);
//                    nodes(i_n) = (*node2newNodes)[mat].value_at((*node2newNodes)[mat].find(nodes(i_n)));
//                }
//            }
//
//            c.setNodes(nodes);
//        }
//
//    };
//
//    /*----------------------------------------------------------------------------*/
//    struct Pillow_replaceInCells_3D {
//        const kmds::GrowingView<kmds::TCellID> *selection;
//        kmds::Mesh *mesh;
//        const elg3d::MaterialAssignment* ma;
//        const kmds::Variable<bool>* varIsInterfaceNode;
//        const std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> >* node2newNodes;
//
//
//        Pillow_replaceInCells_3D(
//                const kmds::GrowingView<kmds::TCellID> *selection_,
//                kmds::Mesh *mesh_,
//                const elg3d::MaterialAssignment* ma_,
//                const kmds::Variable<bool>* varIsInterfaceNode_,
//                const std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> >* node2newNodes_
//        )
//                : selection(selection_)
//                , mesh(mesh_)
//                , ma(ma_)
//                , varIsInterfaceNode(varIsInterfaceNode_)
//                , node2newNodes(node2newNodes_)
//        {
//        }
//
//        KOKKOS_INLINE_FUNCTION
//        void
//        operator()(int i) const {
//            int rid = selection->get(i);
//
//            kmds::Region c = mesh->getRegion(rid);
//            Kokkos::View<kmds::TCellID *> nodes;
//            c.nodeIds(nodes);
//
//            for (int i_n = 0; i_n < nodes.size(); i_n++) {
//
//                if((*varIsInterfaceNode)[nodes(i_n)]) {
//                    int mat = ma->getMaterial(rid);
//                    nodes(i_n) = (*node2newNodes)[mat].value_at((*node2newNodes)[mat].find(nodes(i_n)));
//                }
//            }
//
//            c.setNodes(nodes);
//        }
//
//    };

    /*----------------------------------------------------------------------------*/
    void
    extractGeomModel_extract_2D(kmds::Mesh *AMesh,
                                const elg3d::MaterialAssignment *Ama,
                                const kmds::Connectivity *c_A2C,
                                gmds::cad::FACManager *AGeomModel)
    {

        // mesh that will be used as the mesh representation for the geom model
        kmds::Mesh mesh_tmp;

        kmds::Variable<std::pair<int, int> > *varMatAssign = mesh_tmp.createVariable<std::pair<int, int> >(
                std::pair<int, int>(-1, -1), kmds::KMDS_EDGE,
                "tmp_varMatAssign");

        kmds::GrowingView<kmds::TCellID> selectionN("NODES", AMesh->getNbNodes());
        kmds::GrowingView<kmds::TCellID> selectionA("FACETS", AMesh->getNbEdges());
        extractGeomModel_extract_AandN_2D(AMesh, Ama, c_A2C, &selectionN, &selectionA, &mesh_tmp, varMatAssign);


        kmds::Connectivity* c_N2E = mesh_tmp.createConnectivity(kmds::N2E);
        kmds::ConnectivityHelper ch(&mesh_tmp);
        ch.buildN2E();

        // build the geom model from the extracted mesh and its varMatAssign
        extractGeomModel_buildGeomModel_2D(&mesh_tmp, varMatAssign, c_N2E, AGeomModel);


        // clean-up
        mesh_tmp.deleteVariable(kmds::KMDS_EDGE, "tmp_varMatAssign");

        mesh_tmp.deleteConnectivity(kmds::N2E);
    }

    /*----------------------------------------------------------------------------*/
    void
    extractGeomModel_extract_3D(kmds::Mesh *AMesh,
                                const elg3d::MaterialAssignment *Ama,
                                const kmds::Connectivity *c_A2C,
                                gmds::cad::FACManager *AGeomModel,
                                kmds::Variable <std::uintptr_t> *AGeomassoc)
    {

        // mesh that will be used as the mesh representation for the geom model
        kmds::Mesh mesh_tmp;

        kmds::Variable<std::pair<int, int> > *varMatAssign = mesh_tmp.createVariable<std::pair<int, int> >(
                std::pair<int, int>(-1, -1), kmds::KMDS_FACE,
                "tmp_varMatAssign");

        kmds::GrowingView<kmds::TCellID> selectionN("NODES", AMesh->getNbNodes());
        kmds::GrowingView<kmds::TCellID> selectionA("FACETS", AMesh->getNbFaces());

        kmds::Variable<kmds::TCellID > *varInterface2Ref_N = mesh_tmp.createVariable<kmds::TCellID >(
                kmds::NullID, kmds::KMDS_NODE, "tmp_varInterface2Ref_N");

        extractGeomModel_extract_AandN_3D(AMesh, Ama, c_A2C, &selectionN, &selectionA, &mesh_tmp, varMatAssign, varInterface2Ref_N);

        std::cout<<"selectionN "<<selectionN.getNbElems()<<std::endl;
        std::cout<<"selectionA "<<selectionA.getNbElems()<<std::endl;

        kmds::Connectivity* c_E2F = mesh_tmp.createConnectivity(kmds::E2F);
        kmds::ConnectivityHelper ch(&mesh_tmp);
        ch.buildEandE2F_3D();

        kmds::Variable<std::array<char, 100> > *varCurveMat = mesh_tmp.createVariable<std::array<char, 100> >(
                std::array<char, 100>(), kmds::KMDS_EDGE,
                "tmp_varCurveMat");

        kmds::GrowingView<kmds::TCellID> selectionN_vertices("NODES", AMesh->getNbNodes());

        std::cout<<"mesh_tmp.getNbEdges() before keep "<<mesh_tmp.getNbEdges()<<std::endl;
        extractGeomModel_keep_EandN_3D(&mesh_tmp, c_E2F ,varMatAssign, varCurveMat, &selectionN_vertices);
        std::cout<<"mesh_tmp.getNbEdges() after keep "<<mesh_tmp.getNbEdges()<<std::endl;

        std::cout<<"nbNodes after keep "<<selectionN_vertices.getNbElems()<<std::endl;


        kmds::Connectivity* c_N2E = mesh_tmp.createConnectivity(kmds::N2E);
        ch.buildN2E();

        kmds::Variable<std::array<char, 100> > *varVertexMat = mesh_tmp.createVariable<std::array<char, 100> >(
                std::array<char, 100>(), kmds::KMDS_NODE,
                "tmp_varVertexMat");

        extractGeomModel_vertexMat_N_3D(&mesh_tmp, c_N2E, varCurveMat, &selectionN_vertices, varVertexMat);
        std::cout<<"selectionN_vertices.getNbElems() "<<selectionN_vertices.getNbElems()<<std::endl;

        kmds::VTKWriter<kmds::Mesh> w(mesh_tmp);
        w.write("debug_model", kmds::F);

        // build the geom model from the extracted mesh and its varMatAssign
        kmds::Variable <std::uintptr_t> *varGeomassoc_tmp = mesh_tmp.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr),
                                                                                                    kmds::KMDS_NODE, "tmp_varGeomassoc_tmp");

        extractGeomModel_buildGeomModel_3D(&mesh_tmp, varMatAssign, varCurveMat, varVertexMat, &selectionN_vertices, AGeomModel, varGeomassoc_tmp);


        // get the geomassoc back to the mesh
        const kmds::TCellID nbNodes_interface = mesh_tmp.getNbNodes();
        kmds::GrowingView<kmds::TCellID> nodesIDs_interface("NODES", nbNodes_interface);
        mesh_tmp.getNodeIDs(&nodesIDs_interface);

        Kokkos::parallel_for(nbNodes_interface,
                             KOKKOS_LAMBDA(const int i)
                             {
                                 const kmds::TCellID nid_interface = nodesIDs_interface.get(i);
                                 const kmds::TCellID nid = (*varInterface2Ref_N)[nid_interface];

                                 (*AGeomassoc)[nid] = (*varGeomassoc_tmp)[nid_interface];
                             }
        );


        // clean-up
        mesh_tmp.deleteVariable(kmds::KMDS_NODE, "tmp_varInterface2Ref_N");

        mesh_tmp.deleteVariable(kmds::KMDS_FACE, "tmp_varMatAssign");

        mesh_tmp.deleteVariable(kmds::KMDS_EDGE, "tmp_varCurveMat");
        mesh_tmp.deleteVariable(kmds::KMDS_NODE, "tmp_varVertexMat");

        mesh_tmp.deleteConnectivity(kmds::E2F);
        mesh_tmp.deleteConnectivity(kmds::N2E);
    }

    /*----------------------------------------------------------------------------*/
    void extractGeomModel_extract_AandN_2D(kmds::Mesh* AMesh,
                                           const elg3d::MaterialAssignment* Ama,
                                           const kmds::Connectivity* c_A2C,
                                           kmds::GrowingView<kmds::TCellID>* ASelectionN,
                                           kmds::GrowingView<kmds::TCellID>* ASelectionA,
                                           kmds::Mesh* AMesh_interface,
                                           kmds::Variable<std::pair<int, int> >* AVarMatAssign)

    {

        // get the facets that will represent the interfaces
        MaterialInterfaces_getEdgeOnInterfacesWithBoundary(AMesh,
                                                           c_A2C,
                                                           Ama,
                                                           ASelectionA);

        int nbFacets_tokeep = ASelectionA->getNbElems();
        std::cout<<"nbFacets_tokeep "<<nbFacets_tokeep<<std::endl;

        // identify the nodes of those facets
        Kokkos::View<bool*, Kokkos::MemoryTraits<Kokkos::Atomic> > markNodesToKeep("NODES_TO_KEEP", AMesh->getNbNodes());

        Kokkos::parallel_for(nbFacets_tokeep,
                             KOKKOS_LAMBDA(const int i)
                             {
                                 const kmds::TCellID aid = ASelectionA->get(i);
                                 kmds::Edge a = AMesh->getEdge(aid);

                                 Kokkos::View<kmds::TCellID*> nids;
                                 a.nodeIds(nids);

                                 markNodesToKeep(nids(0)) = true;
                                 markNodesToKeep(nids(1)) = true;
                             }
        );

        Kokkos::parallel_for(AMesh->getNbNodes(),
                             KOKKOS_LAMBDA(const int i)
                             {
                                 if(markNodesToKeep(i)) {
                                     ASelectionN->push_back(i);
                                 }
                             }
        );

        int nbNodes_tokeep = ASelectionN->getNbElems();
        std::cout<<"nbNodes_tokeep "<<nbNodes_tokeep<<std::endl;


        // store the relevent nodes and facets in the interface mesh
        AMesh_interface->updateNodeCapacity(nbNodes_tokeep);
        AMesh_interface->updateNodeCapacity(nbFacets_tokeep);

        // create the nodes
        kmds::Variable<kmds::TCellID > *varNew2Old_N = AMesh->createVariable<kmds::TCellID >(
                kmds::NullID, kmds::KMDS_NODE, "tmp_varNew2Old_N");

        const kmds::TCellID nid_interface = AMesh_interface->addNodes(nbNodes_tokeep);

        Kokkos::parallel_for(nbNodes_tokeep,
                             KOKKOS_LAMBDA(const int i)
                             {
                                 const kmds::TCellID nid = ASelectionN->get(i);

                                 kmds::TCoord x, y, z;
                                 AMesh->getNodeLocation(nid, x, y, z);

                                 AMesh_interface->setNodeLocation(nid_interface + i, x, y, z);

                                 (*varNew2Old_N)[nid] = nid_interface + i;
                             }
        );

        // create the facets
        const kmds::TCellID aid_interface = AMesh_interface->addEdges(nbFacets_tokeep);
        if(aid_interface != 0) {
            throw kmds::KException("extractGeomModel_extract_AandN_2D aid_interface should be equal to zero");
        }

        Kokkos::parallel_for(nbFacets_tokeep,
                             KOKKOS_LAMBDA(const int i)
                             {
                                 const kmds::TCellID aid = ASelectionA->get(i);
                                 kmds::Edge a = AMesh->getEdge(aid);

                                 Kokkos::View<kmds::TCellID*> nids;
                                 a.nodeIds(nids);

                                 kmds::Edge a_interface = AMesh_interface->getEdge(aid_interface + i);

                                 a_interface.setNodes((*varNew2Old_N)[nids(0)], (*varNew2Old_N)[nids(1)]);
                             }
        );

        // store the material assignement of the kept facets
        kmds::GrowingView<kmds::TCellID> facetsIDs_source("FACETS", AMesh->getNbEdges());
        AMesh->getEdgeIDs(&facetsIDs_source);

        Kokkos::parallel_for(nbFacets_tokeep,
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID aid = ASelectionA->get(i);

                                 Kokkos::View<kmds::TCellID *> cids;
                                 c_A2C->get(aid, cids);

                                 int first = -1;
                                 int second = Ama->getMaterial(cids(0));

                                 if(cids.size() == 2) {
                                     if(Ama->getMaterial(cids(1)) > second) {
                                        std::swap(first, second);
                                     }

                                     (*AVarMatAssign)[aid_interface + i] = std::pair<int, int>(first, second);
                                 }

                             });


        // clean-up
        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varNew2Old_N");
    }

    /*----------------------------------------------------------------------------*/
    void extractGeomModel_extract_AandN_3D(kmds::Mesh* AMesh,
                                           const elg3d::MaterialAssignment* Ama,
                                           const kmds::Connectivity* c_A2C,
                                           kmds::GrowingView<kmds::TCellID>* ASelectionN,
                                           kmds::GrowingView<kmds::TCellID>* ASelectionA,
                                           kmds::Mesh* AMesh_interface,
                                           kmds::Variable<std::pair<int, int> >* AVarMatAssign,
                                           kmds::Variable<kmds::TCellID >* AVarInterface2Ref_N)

    {

        // get the facets that will represent the interfaces
        MaterialInterfaces_getFaceOnInterfacesWithBoundary(AMesh,
                                                           c_A2C,
                                                           Ama,
                                                           ASelectionA);

        const int nbFacets_tokeep = ASelectionA->getNbElems();

        // identify the nodes of those facets
        Kokkos::View<bool*, Kokkos::MemoryTraits<Kokkos::Atomic> > markNodesToKeep("NODES_TO_KEEP", AMesh->getNbNodes());

        Kokkos::parallel_for(nbFacets_tokeep,
                             KOKKOS_LAMBDA(const int i)
                             {
                                 const kmds::TCellID aid = ASelectionA->get(i);
                                 kmds::Face a = AMesh->getFace(aid);

                                 Kokkos::View<kmds::TCellID*> nids;
                                 a.nodeIds(nids);

                                 for(int n=0; n<nids.extent(0); n++) {
                                     markNodesToKeep(nids(n)) = true;
                                 }
                             }
        );

        const int nbNodes = AMesh->getNbNodes();
        Kokkos::parallel_for(nbNodes,
                             KOKKOS_LAMBDA(const int i)
                             {
                                 if(markNodesToKeep(i)) {
                                     ASelectionN->push_back(i);
                                 }
                             }
        );

        const int nbNodes_tokeep = ASelectionN->getNbElems();


        // store the relevent nodes and facets in the interface mesh
        AMesh_interface->updateNodeCapacity(nbNodes_tokeep);
        AMesh_interface->updateFaceCapacity(nbFacets_tokeep);

        // create the nodes
        kmds::Variable<kmds::TCellID > *varNew2Old_N = AMesh->createVariable<kmds::TCellID >(
                kmds::NullID, kmds::KMDS_NODE, "tmp_varNew2Old_N");

        const kmds::TCellID nid_interface = AMesh_interface->addNodes(nbNodes_tokeep);

        Kokkos::parallel_for(nbNodes_tokeep,
                             KOKKOS_LAMBDA(const int i)
                             {
                                 const kmds::TCellID nid = ASelectionN->get(i);

                                 kmds::TCoord x, y, z;
                                 AMesh->getNodeLocation(nid, x, y, z);

                                 AMesh_interface->setNodeLocation(nid_interface + i, x, y, z);

                                 (*varNew2Old_N)[nid] = nid_interface + i;

                                 (*AVarInterface2Ref_N)[nid_interface + i] = nid;
                             }
        );

        // create the facets
        kmds::Variable<kmds::TCellID > *varNew2Old = AMesh_interface->createVariable<kmds::TCellID >(
                kmds::NullID, kmds::KMDS_FACE, "tmp_varNew2Old");

        Kokkos::parallel_for(nbFacets_tokeep,
                             KOKKOS_LAMBDA(const int i)
                             {
                                 const kmds::TCellID aid = ASelectionA->get(i);
                                 kmds::Face a = AMesh->getFace(aid);

                                 Kokkos::View<kmds::TCellID*> nids;
                                 a.nodeIds(nids);

                                 kmds::TCellID nids_new[nids.extent(0)];
                                 for(int n=0; n<nids.extent(0); n++) {
                                     nids_new[n] = (*varNew2Old_N)[nids(n)];
                                 }

                                 (*varNew2Old)[AMesh_interface->newFace(nids_new, nids.extent(0))] = aid;
                             }
        );

        // store the material assignement of the kept facets
        kmds::GrowingView<kmds::TCellID> facetsIDs_interface("FACETS", AMesh_interface->getNbFaces());
        AMesh_interface->getFaceIDs(&facetsIDs_interface);

        Kokkos::parallel_for(nbFacets_tokeep,
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID aid = facetsIDs_interface.get(i);

                                 kmds::TCellID old_aid = (*varNew2Old)[aid];

                                 Kokkos::View<kmds::TCellID *> cids;
                                 c_A2C->get(old_aid, cids);

                                 int first = -1;
                                 int second = Ama->getMaterial(cids(0));

                                 if(cids.size() == 2) {

                                     first = Ama->getMaterial(cids(1));

                                     if (first > second) {
                                         std::swap(first, second);
                                     }
                                 }

                                 (*AVarMatAssign)[aid] = std::pair<int, int>(first, second);
                             });


        // clean-up
        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varNew2Old_N");
        AMesh_interface->deleteVariable(kmds::KMDS_FACE, "tmp_varNew2Old");
    }

    /*----------------------------------------------------------------------------*/
    void extractGeomModel_keep_EandN_3D(kmds::Mesh* AMesh_interface,
                                    const kmds::Connectivity* c_E2F,
                                    const kmds::Variable<std::pair<int, int> >* AVarMatAssign,
                                    kmds::Variable<std::array<char, 100> >* AVarCurveMat,
                                    kmds::GrowingView<kmds::TCellID>* ASelectionN_vertices)
    {

        // keep count of the number of curve edges adjacent to the nodes
        Kokkos::View<int*, Kokkos::MemoryTraits<Kokkos::Atomic> > nodesToEdges("NODES_TO_EDGES", AMesh_interface->getNbNodes());


        const kmds::TCellID nbEdges_tot = AMesh_interface->getNbEdges();

        kmds::GrowingView<kmds::TCellID> edgesIDs("EDGES", nbEdges_tot);
        AMesh_interface->getEdgeIDs(&edgesIDs);


        Kokkos::parallel_for(nbEdges_tot,
                             KOKKOS_LAMBDA(const int i)
                             {
                                 const kmds::TCellID eid = edgesIDs.get(i);

                                 Kokkos::View<kmds::TCellID *> fids;
                                 c_E2F->get(eid, fids);

                                 bool tokeep = false;

                                 // is this a boundary edge (size == 1) or curve edge (size > 2) ?
                                 if(fids.size() != 2) {
                                     tokeep = true;
                                 }

                                 if(!tokeep) {
                                     // remove the edge
                                     AMesh_interface->removeEdge(eid);
                                 } else {
                                     // get the ordered set of materials
                                     std::set<int> matSet;

                                     for(int i_f=0; i_f<fids.size(); i_f++) {
                                         matSet.insert(((*AVarMatAssign)[fids[i_f]]).first);
                                         matSet.insert(((*AVarMatAssign)[fids[i_f]]).second);
                                     }

                                     std::string curve_name("curve_");
                                     for(auto mat: matSet) {
                                         curve_name += std::to_string(mat) + std::string("_");
                                     }

                                     for(int l=0; l<curve_name.length(); l++) {
                                         (*AVarCurveMat)[eid][l] = curve_name[l];
                                     }

//                                     std::cout<<"curve_name "<<curve_name<<std::endl;


                                    // update the count for the nodes
                                     kmds::Edge e = AMesh_interface->getEdge(eid);

                                     Kokkos::View<kmds::TCellID *> nids;
                                     e.nodeIds(nids);

                                     nodesToEdges[nids[0]]++;
                                     nodesToEdges[nids[1]]++;
                                 }
                             }
        );


        const kmds::TCellID nbNodes = AMesh_interface->getNbNodes();

        kmds::GrowingView<kmds::TCellID> nodesIDs("NODES", nbNodes);
        AMesh_interface->getNodeIDs(&nodesIDs);


        Kokkos::parallel_for(nbNodes,
                             KOKKOS_LAMBDA(const int i) {
                                 const kmds::TCellID nid = nodesIDs.get(i);

                                 if(nodesToEdges[nid] > 2) {
                                     ASelectionN_vertices->push_back(nid);

                                     std::cout<<"nid "<<nid<<" nbmarked "<<nodesToEdges[nid]<<std::endl;
                                 }
                             }
        );

        // clean-up

    }

    /*----------------------------------------------------------------------------*/
    void extractGeomModel_vertexMat_N_3D(kmds::Mesh* AMesh_interface,
                                         const kmds::Connectivity* c_N2E,
                                         const kmds::Variable<std::array<char, 100> >* AVarCurveMat,
                                         const kmds::GrowingView<kmds::TCellID>* ASelectionN_vertices,
                                         kmds::Variable<std::array<char, 100> >* AVarVertexMat)
    {
        Kokkos::parallel_for(ASelectionN_vertices->getNbElems(),
                             KOKKOS_LAMBDA(const int i)
                             {
                                 const kmds::TCellID nid = ASelectionN_vertices->get(i);

                                 Kokkos::View<kmds::TCellID *> eids;
                                 c_N2E->get(nid, eids);

                                 std::cout<<"nid "<<nid<<" eids.size() "<<eids.size()<<std::endl;

                                 std::set<int> matSet;


                                 for(int i_e=0; i_e<eids.size(); i_e++) {
                                     std::string str((*AVarCurveMat)[eids[i_e]].begin(), (*AVarCurveMat)[eids[i_e]].end());

                                     std::cout<<"str "<<str<<std::endl;

                                     std::string::size_type current_pos = str.find_first_of('_');
                                     std::string::size_type next_pos = str.find('_', current_pos+1);

                                     std::cout<<"str "<<str<<" "<<current_pos<<" "<<next_pos<<std::endl;

                                     while(std::string::npos != next_pos) {

                                         const int mat = atoi((*AVarCurveMat)[eids[i_e]].data()+current_pos+1);
                                         matSet.insert(mat);

                                         current_pos = next_pos;
                                         next_pos = str.find('_', current_pos+1);
                                     }
                                 }

                                 std::cout<<"matSet.size() "<<matSet.size()<<std::endl;

                                 std::string vertex_name("vertex_");
                                 for(auto mat: matSet) {
                                     vertex_name += std::to_string(mat) + std::string("_");
                                 }

                                 for(int l=0; l<vertex_name.length(); l++) {
                                     (*AVarVertexMat)[nid][l] = vertex_name[l];
                                 }

                                 std::cout<<"vertex_name "<<vertex_name<<std::endl;
                             }
        );


        // clean-up

    }

    /*----------------------------------------------------------------------------*/
    void
    extractGeomModel_buildGeomModel_2D(kmds::Mesh *AMesh,
                                       const kmds::Variable<std::pair<int, int> >* AVarMatAssign,
                                       const kmds::Connectivity* c_N2E,
                                       gmds::cad::FACManager *AGeomModel)

    {

        // identify the nodes that correspond to geometric points
        // based on |adjacentE(n)| >= 3
        kmds::GrowingView<kmds::TCellID> nodesIDs("NODES", AMesh->getNbNodes());
        AMesh->getNodeIDs(&nodesIDs);

        const kmds::TCellID nbNodes = nodesIDs.getNbElems();

        std::set<kmds::TCellID> nodeAsGeomPoints;

        for(int i=0; i<nbNodes; i++) {

            const kmds::TCellID nid = nodesIDs.get(i);

            Kokkos::View<kmds::TCellID *> eids;
            c_N2E->get(nid, eids);

            if(eids.extent(0) >= 3) {

                // create a GeomPoint
                kmds::TCoord xyz[3];
                AMesh->getNodeLocation(nid, xyz[0], xyz[1], xyz[2]);
                std::cout<<"point created "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<std::endl;

                nodeAsGeomPoints.insert(nid);
            }
        }


        // identify the curves
        kmds::GrowingView<kmds::TCellID> edgesIDs("EDGES", AMesh->getNbEdges());
        AMesh->getEdgeIDs(&edgesIDs);

        const kmds::TCellID nbEdges = edgesIDs.getNbElems();

        kmds::Variable<bool> *varFreeEdges = AMesh->createVariable<bool>(true, kmds::KMDS_EDGE, "tmp_varFreeEdges");

        std::set<kmds::TCellID> freeEdges;


        for(int i=0; i<nbEdges; i++) {

            const kmds::TCellID eid = edgesIDs.get(i);
            freeEdges.insert(eid);
        }


        kmds::TCellID nbEdges_free = nbEdges;

//        while(nbEdges_free > 0) {
        while(!freeEdges.empty()) {

            // select an edge
            std::set<unsigned int>::iterator it = freeEdges.begin();
            kmds::TCellID eid_first = *it;

            freeEdges.erase(it);
            (*varFreeEdges)[eid_first] = false;
            std::list<kmds::TCellID> edgesCurve {eid_first};

            kmds::Edge e = AMesh->getEdge(eid_first);

            Kokkos::View<kmds::TCellID *> nids;
            e.nodeIds(nids);

            kmds::TCellID n_first = nids(0);
            kmds::TCellID n_last = nids(1);

            std::list<kmds::TCellID> nodesCurve {n_first, n_last};

            // propagate until the next node is
            // - a geom point
            // - the first node in case of a loop
            kmds::TCellID n_prev = n_first;
            kmds::TCellID n_next = n_last;

            while((n_next != n_first) && (nodeAsGeomPoints.find(n_next) == nodeAsGeomPoints.end()) ) {

                Kokkos::View<kmds::TCellID *> eids;
                c_N2E->get(n_next, eids);

                // |eids| should be equal to two
                kmds::TCellID eid_next;
                if((*varFreeEdges)[eids(0)]) {
                    eid_next = eids(0);
                } else {
                    eid_next = eids(1);
                }

                freeEdges.erase(eid_next);
                (*varFreeEdges)[eid_next] = false;
                edgesCurve.push_back(eid_next);

                kmds::Edge e_next = AMesh->getEdge(eid_next);

                Kokkos::View<kmds::TCellID *> nids;
                e_next.nodeIds(nids);

                if(n_next == nids(0)) {
                    n_next = nids(1);
                } else {
                    n_next = nids(0);
                }

                nodesCurve.push_back(n_next);

            }

            n_last = n_next;

            // case of a loop
            if(n_last == n_first) {

                // create a geom point
                //TODO : is it necessary ?
                nodeAsGeomPoints.insert(n_first);


                continue;
            }

            // propagate until the prev node is
            // - a geom point

            while(nodeAsGeomPoints.find(n_prev) == nodeAsGeomPoints.end()) {

                Kokkos::View<kmds::TCellID *> eids;
                c_N2E->get(n_prev, eids);

                // |eids| should be equal to two
                kmds::TCellID eid_prev;
                if((*varFreeEdges)[eids(0)]) {
                    eid_prev = eids(0);
                } else {
                    eid_prev = eids(1);
                }

                freeEdges.erase(eid_prev);
                (*varFreeEdges)[eid_prev] = false;
                edgesCurve.push_front(eid_prev);

                kmds::Edge e_prev = AMesh->getEdge(eid_prev);

                Kokkos::View<kmds::TCellID *> nids;
                e_prev.nodeIds(nids);

                if(n_prev == nids(0)) {
                    n_prev = nids(1);
                } else {
                    n_prev = nids(0);
                }

                nodesCurve.push_front(n_prev);
            }

            n_first = n_prev;


            // create the curve
            std::cout<<"curve "<<std::endl;
            for(auto n: nodesCurve) {
                double xyz[3];
                AMesh->getNodeLocation(n, xyz[0], xyz[1], xyz[2]);
                std::cout<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<std::endl;
            }


        }



    }

    /*----------------------------------------------------------------------------*/
    void
    extractGeomModel_buildGeomModel_3D(kmds::Mesh *AMesh,
                                       const kmds::Variable<std::pair<int, int> >* AVarMatAssign,
                                       const kmds::Variable<std::array<char, 100> >* AVarCurveMat,
                                       const kmds::Variable<std::array<char, 100> >* AVarVertexMat,
                                       const kmds::GrowingView<kmds::TCellID>* ASelectionN_vertices,
                                       gmds::cad::FACManager *AGeomModel,
                                       kmds::Variable <std::uintptr_t> *AVarGeomassoc)

    {
        // get the geometric model underlying gmds mesh representation
        gmds::Mesh& gmdsmesh = AGeomModel->getMeshView();


        // get the nodes and recreate them in gmds
        kmds::GrowingView<kmds::TCellID> nodesIDs("NODES", AMesh->getNbNodes());
        AMesh->getNodeIDs_dummy(&nodesIDs);

        kmds::Variable<kmds::TCellID > *varInterface2gmds_N = AMesh->createVariable<kmds::TCellID >(
                kmds::NullID, kmds::KMDS_NODE, "tmp_varInterface2gmds_N");


        for(int i=0; i<nodesIDs.getNbElems(); i++) {
            kmds::TCellID nid = nodesIDs.get(i);
            kmds::TCoord xyz[3];
            AMesh->getNodeLocation(nid, xyz[0], xyz[1], xyz[2]);
            gmds::Node n = gmdsmesh.newNode(xyz[0], xyz[1], xyz[2]);

            (*varInterface2gmds_N)[nid] = n.id();
        }


        // get the faces and recreate them in gmds using triangles
        kmds::GrowingView<kmds::TCellID> facesIDs("FACES", AMesh->getNbFaces());
        AMesh->getFaceIDs_dummy(&facesIDs);

        kmds::Variable<kmds::TCellID > *varOld2new_F = AMesh->createVariable<kmds::TCellID >(
                kmds::NullID, kmds::KMDS_FACE, "tmp_varOld2new_F");
        kmds::Variable<int > *varOld2new_F_nbTri = AMesh->createVariable<int >(
                0, kmds::KMDS_FACE, "tmp_varOld2new_F_nbTri");

        for(int i=0; i<facesIDs.getNbElems(); i++) {
            kmds::TCellID fid = facesIDs.get(i);

            kmds::Face f = AMesh->getFace(fid);

            Kokkos::View<kmds::TCellID *> nids;
            f.nodeIds(nids);

            if(kmds::KMDS_TRIANGLE == f.computeType()) {
                gmds::Face f_first = gmdsmesh.newTriangle((*varInterface2gmds_N)[nids(0)], (*varInterface2gmds_N)[nids(1)], (*varInterface2gmds_N)[nids(2)]);
                (*varOld2new_F)[fid] = f_first.id();
                (*varOld2new_F_nbTri)[fid] = 1;

            } else if(kmds::KMDS_QUAD == f.computeType()) {
                gmds::math::Point pt = f.midpoint();
                gmds::Node n = gmdsmesh.newNode(pt);

                gmds::Face f_first = gmdsmesh.newTriangle((*varInterface2gmds_N)[nids(0)], (*varInterface2gmds_N)[nids(1)], n.id());
                gmdsmesh.newTriangle((*varInterface2gmds_N)[nids(1)], (*varInterface2gmds_N)[nids(2)], n.id());
                gmdsmesh.newTriangle((*varInterface2gmds_N)[nids(2)], (*varInterface2gmds_N)[nids(3)], n.id());
                gmdsmesh.newTriangle((*varInterface2gmds_N)[nids(3)], (*varInterface2gmds_N)[nids(0)], n.id());

                (*varOld2new_F)[fid] = f_first.id();
                (*varOld2new_F_nbTri)[fid] = 4;
            }
        }


        // mark for the already treated faces
        kmds::Variable<bool > *varTreatedFace = AMesh->createVariable<bool >(
                false, kmds::KMDS_FACE, "tmp_varTreatedFace");

        kmds::TCellID nbTreatedFaces = 0;
        const kmds::TCellID nbFaces = facesIDs.getNbElems();


        // dispatch the faces (their triangles) into several surfaces
        while(nbFaces != nbTreatedFaces) {

            std::vector<gmds::Face> faces_to_add;
            std::pair<int, int> reference_pair;

            gmds::cad::FACSurface *surf = dynamic_cast<gmds::cad::FACSurface *> (AGeomModel->newSurface());

            std::set<kmds::TCellID> nodes_of_surf;
            std::set<kmds::TCellID> faces_of_surf;

            bool found_pair = false;
            for(int i=0; i<nbFaces; i++) {
                kmds::TCellID fid = facesIDs.get(i);

                // check whether face was already treated
                if(!(*varTreatedFace)[fid]) {

                    if(!found_pair) {

                        // we are working on a new interface
                        reference_pair = (*AVarMatAssign)[fid];

                        for(int itri=0; itri<(*varOld2new_F_nbTri)[fid]; itri++) {

                            faces_to_add.push_back(gmdsmesh.get<gmds::Face>((*varOld2new_F)[fid] + itri));
                        }

                        (*varTreatedFace)[fid] = true;
                        found_pair = true;
                        nbTreatedFaces++;

                        // reinitialize the set of kmds interface nodes
                        nodes_of_surf.clear();
                        faces_of_surf.clear();
                        faces_of_surf.insert(fid);

                        kmds::Face f = AMesh->getFace(fid);
                        Kokkos::View<kmds::TCellID *> nids;
                        f.nodeIds(nids);

                        for(int i_n=0; i_n<nids.size(); i_n++) {
                            nodes_of_surf.insert(nids[i_n]);
                        }

                    } else {
                        // we are identifying the faces of that interface
                        if(reference_pair == (*AVarMatAssign)[fid]) {

                            for(int itri=0; itri<(*varOld2new_F_nbTri)[fid]; itri++) {

                                faces_to_add.push_back(gmdsmesh.get<gmds::Face>((*varOld2new_F)[fid] + itri));
                            }

                            (*varTreatedFace)[fid] = true;
                            nbTreatedFaces++;

                            kmds::Face f = AMesh->getFace(fid);
                            Kokkos::View<kmds::TCellID *> nids;
                            f.nodeIds(nids);

                            for(int i_n=0; i_n<nids.size(); i_n++) {
                                nodes_of_surf.insert(nids[i_n]);
                            }
                            faces_of_surf.insert(fid);
                        }

                    }
                }
            }

            if(!found_pair) {
                throw kmds::KException("extractGeomModel_buildGeomModel_3D no reference_pair found.");
            }

            std::string surfname = std::string("interface_") + std::to_string(reference_pair.first) + std::string("_") + std::to_string(reference_pair.second);
            surf->setName(surfname);
            surf->setMeshFaces(faces_to_add);

            std::cout<<"surfname "<<surfname<<" faces_to_add.size() "<<faces_to_add.size()<<std::endl;

//            nbTreatedFaces += faces_to_add.size();

            // associate the nodes of the kmds interface
            for(auto nid: nodes_of_surf) {
                (*AVarGeomassoc)[nid] = reinterpret_cast<std::uintptr_t>(surf);

                std::cout<<"nodes_of_surf "<<nid<<std::endl;
            }
            std::cout<<"nodes_of_surf.size() "<<nodes_of_surf.size()<<std::endl;
            std::cout<<"faces_of_surf.size() "<<faces_of_surf.size()<<std::endl;
        }


        // get the edges
        kmds::GrowingView<kmds::TCellID> edgesIDs("EDGES", AMesh->getNbEdges());
        AMesh->getEdgeIDs(&edgesIDs);

        // mark for the already treated edges
        kmds::Variable<bool > *varTreatedEdge = AMesh->createVariable<bool >(
                false, kmds::KMDS_EDGE, "tmp_varTreatedEdge");

        kmds::TCellID nbTreatedEdges = 0;
        const kmds::TCellID nbEdges = edgesIDs.getNbElems();


        // dispatch the edges into several curves
        while(nbEdges != nbTreatedEdges) {

            std::vector<gmds::Edge> edges_to_add;
            std::array<char, 100> reference_string{};

            gmds::cad::FACCurve *curv = dynamic_cast<gmds::cad::FACCurve *> (AGeomModel->newCurve());

            std::set<kmds::TCellID> nodes_of_curv;

            bool found_curve = false;
            for(int i=0; i<nbEdges; i++) {
                kmds::TCellID eid = edgesIDs.get(i);

                // check whether edge was already treated
                if(!(*varTreatedEdge)[eid]) {

                    if(!found_curve) {

                        // we are working on a new curve
                        reference_string = (*AVarCurveMat)[eid];

                        // create the edge
                        kmds::Edge e = AMesh->getEdge(eid);

                        Kokkos::View<kmds::TCellID *> nids;
                        e.nodeIds(nids);

                        gmds::Edge e_model = gmdsmesh.newEdge((*varInterface2gmds_N)[nids[0]], (*varInterface2gmds_N)[nids[1]]);

                        edges_to_add.push_back(e_model);

                        (*varTreatedEdge)[eid] = true;
                        found_curve = true;
                        nbTreatedEdges++;

                        nodes_of_curv.insert(nids[0]);
                        nodes_of_curv.insert(nids[1]);

                    } else {
                        // we are identifying the edges of that curve
                        if(reference_string == (*AVarCurveMat)[eid]) {

                            // create the edge
                            kmds::Edge e = AMesh->getEdge(eid);

                            Kokkos::View<kmds::TCellID *> nids;
                            e.nodeIds(nids);

                            gmds::Edge e_model = gmdsmesh.newEdge((*varInterface2gmds_N)[nids[0]], (*varInterface2gmds_N)[nids[1]]);

                            edges_to_add.push_back(e_model);

                            (*varTreatedEdge)[eid] = true;
                            nbTreatedEdges++;

                            nodes_of_curv.insert(nids[0]);
                            nodes_of_curv.insert(nids[1]);
                        }

                    }
                }
            }

            if(!found_curve) {
                throw kmds::KException("extractGeomModel_buildGeomModel_3D no reference_string found for curve.");
            }

            std::string curvname(reference_string.begin(), reference_string.end());
            curv->setName(curvname);
            curv->setMeshEdges(edges_to_add);

            std::cout<<"curvname "<<curvname<<" edges_to_add.size() "<<edges_to_add.size()<<std::endl;

            // associate the nodes of the kmds interface
            for(auto nid: nodes_of_curv) {
                (*AVarGeomassoc)[nid] = reinterpret_cast<std::uintptr_t>(curv);
            }

            std::cout<<"nodes_of_curv.size() "<<nodes_of_curv.size()<<std::endl;
        }



        // mark for the already treated nodes
        kmds::Variable<bool > *varTreatedNode = AMesh->createVariable<bool >(
                false, kmds::KMDS_NODE, "tmp_varTreatedNode");

        kmds::TCellID nbTreatedNodes = 0;
        const kmds::TCellID nbNodes_vertices = ASelectionN_vertices->getNbElems();

        std::cout<<"========================================="<<std::endl;
        for(int i=0; i<nbNodes_vertices; i++) {
            kmds::TCellID nid = ASelectionN_vertices->get(i);
            std::string poyop((*AVarVertexMat)[nid].begin(), (*AVarVertexMat)[nid].end());
            std::cout<<poyop<<std::endl;

            std::cout<<"position "<<nid<<" "<<AMesh->getNodeLocation(nid)<<std::endl;
        }
        std::cout<<"========================================="<<std::endl;

        // dispatch the nodes into several vertices
        while(nbNodes_vertices != nbTreatedNodes) {

            std::vector<gmds::Node> nodes_to_add;
            std::array<char, 100> reference_string{};

            gmds::cad::FACPoint *vertex = dynamic_cast<gmds::cad::FACPoint *> (AGeomModel->newPoint());

            std::set<kmds::TCellID> nodes_of_vert;

            bool found_vertex = false;
            for(int i=0; i<nbNodes_vertices; i++) {
                kmds::TCellID nid = ASelectionN_vertices->get(i);

                // check whether the node was already treated
                if(!(*varTreatedNode)[nid]) {

                    if(!found_vertex) {

                        // we are working on a new vertex
                        reference_string = (*AVarVertexMat)[nid];

                        nodes_to_add.push_back(gmdsmesh.get<gmds::Node>((*varInterface2gmds_N)[nid]));

                        (*varTreatedNode)[nid] = true;
                        found_vertex = true;
                        nbTreatedNodes++;

                        nodes_of_vert.insert(nid);

                    } else {
                        // we are identifying the nodes of that vertex
                        if(reference_string == (*AVarVertexMat)[nid]) {

                            nodes_to_add.push_back(gmdsmesh.get<gmds::Node>((*varInterface2gmds_N)[nid]));

                            (*varTreatedNode)[nid] = true;
                            nbTreatedNodes++;

                            nodes_of_vert.insert(nid);
                        }

                    }
                }
            }

            if(!found_vertex) {
                throw kmds::KException("extractGeomModel_buildGeomModel_3D no reference_string found for vertex.");
            }

            std::string vertexname(reference_string.begin(), reference_string.end());
            vertex->setName(vertexname);
            vertex->activateMultiNodes();
            vertex->setMeshNodes(nodes_to_add);

            std::cout<<"vertexname "<<vertexname<<" nodes_to_add.size() "<<nodes_to_add.size()<<std::endl;

            // associate the nodes of the kmds interface
            for(auto nid: nodes_of_vert) {
                (*AVarGeomassoc)[nid] = reinterpret_cast<std::uintptr_t>(vertex);
            }

            std::cout<<"nodes_of_vert.size() "<<nodes_of_vert.size()<<std::endl;
        }


        // clean-up
        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varInterface2gmds_N");
        AMesh->deleteVariable(kmds::KMDS_FACE, "tmp_varOld2new_F");

        AMesh->deleteVariable(kmds::KMDS_FACE, "tmp_varTreatedFace");
        AMesh->deleteVariable(kmds::KMDS_EDGE, "tmp_varTreatedEdge");
        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varTreatedNode");
    }
/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
