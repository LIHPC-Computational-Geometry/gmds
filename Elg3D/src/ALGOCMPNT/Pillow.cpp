/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    Pillow.cpp
 *  \author  legoff
 *  \date    05/13/2018
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/Pillow.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <utility>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <Kokkos_UnorderedMap.hpp>
#include <gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/cad/GeomEntity.h>
#include <gmds/math/Point.h>
#include <gmds/math/Triangle.h>

#include <KM/IO/VTKWriter.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include <ELG3D/ALGOCMPNT/Cavity.h>
#include <ELG3D/ALGOCMPNT/Tools.h>
#include <ELG3D/DATACMPNT/Parameters.h>
#include <ELG3D/ALGOCMPNT/SmartLaplacian.h>
#include <ELG3D/ALGOCMPNT/OptimizationSmooth.h>
/*----------------------------------------------------------------------------*/
namespace elg3d {

    struct duplicate_of_node {
        kmds::TCellID nid;
        kmds::TCellID eid0;
        kmds::TCellID eid1;
    };

    struct duplicate_of_node_3D {
        kmds::TCellID nid;
    };

    /*----------------------------------------------------------------------------*/
    struct Pillow_createNodes_xD {
        const kmds::GrowingView<kmds::TCellID> *selection;
        kmds::Mesh *mesh;
        const kmds::Connectivity *c_N2C;
        const elg3d::MaterialAssignment *ma;
        kmds::Variable<std::uintptr_t> *varNodeGeomAssociation;
        std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> > *node2newNodes;


        Pillow_createNodes_xD(

                const kmds::GrowingView<kmds::TCellID> *selection_,
                kmds::Mesh *mesh_,
                const kmds::Connectivity *c_N2C_,
                const elg3d::MaterialAssignment *ma_,
                kmds::Variable<std::uintptr_t> *AVarNodeGeomAssociation_,
                std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> > *node2newNodes_
        )
                : selection(selection_), mesh(mesh_), c_N2C(c_N2C_), ma(ma_),
                  varNodeGeomAssociation(AVarNodeGeomAssociation_), node2newNodes(node2newNodes_) {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int nid = selection->get(i);

            Kokkos::View<kmds::TCellID *> cells;
            c_N2C->get(nid, cells);

            std::set<int> materials;
            for (int i_c = 0; i_c < cells.size(); i_c++) {
                materials.insert(ma->getMaterial(cells(i_c)));
            }


            kmds::TCellID id0 = mesh->addNodes(materials.size());
            kmds::TCoord xyz[3];
            mesh->getNodeLocation(nid, xyz[0], xyz[1], xyz[2]);

            int offset = 0;
            for (auto mat: materials) {

                // compute the new node position
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

                mesh->setNodeLocation(id0 + offset, xyz[0], xyz[1], xyz[2]);

                (*node2newNodes)[mat].insert(nid, id0 + offset);
                (*varNodeGeomAssociation)[id0 + offset] = (*varNodeGeomAssociation)[nid];
                offset++;
            }

        }

    };

    /*----------------------------------------------------------------------------*/
    struct Pillow_createNodes_2D {
        const kmds::GrowingView<kmds::TCellID> *selection;
        kmds::Mesh *mesh;
        const kmds::Connectivity *c_N2C;
        const elg3d::MaterialAssignment *ma;
        kmds::Variable<std::uintptr_t> *varNodeGeomAssociation;
        std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> > *node2newNodes;


        Pillow_createNodes_2D(

                const kmds::GrowingView<kmds::TCellID> *selection_,
                kmds::Mesh *mesh_,
                const kmds::Connectivity *c_N2C_,
                const elg3d::MaterialAssignment *ma_,
                kmds::Variable<std::uintptr_t> *AVarNodeGeomAssociation_,
                std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> > *node2newNodes_
        )
                : selection(selection_), mesh(mesh_), c_N2C(c_N2C_), ma(ma_),
                  varNodeGeomAssociation(AVarNodeGeomAssociation_), node2newNodes(node2newNodes_) {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int nid = selection->get(i);

            Kokkos::View<kmds::TCellID *> cells;
            c_N2C->get(nid, cells);

            std::set<int> materials;
            for (int i_c = 0; i_c < cells.size(); i_c++) {
                materials.insert(ma->getMaterial(cells(i_c)));
            }


            kmds::TCellID id0 = mesh->addNodes(materials.size());
            kmds::TCoord xyz[3];
            mesh->getNodeLocation(nid, xyz[0], xyz[1], xyz[2]);

            int offset = 0;
            for (auto mat: materials) {

                // compute the new node position
                gmds::math::Point pt_old(xyz[0], xyz[1], xyz[2]);
                gmds::math::Point pt_new(0., 0., 0.);
                int nbContrib = 0;
                for (int i_c = 0; i_c < cells.size(); i_c++) {
                    if (ma->getMaterial(cells(i_c)) == mat) {
                        gmds::math::Point pt_midpoint = mesh->getFace(cells(i_c)).midpoint();
                        gmds::math::Point pt_tmp(pt_old + (1 / 5.) * gmds::math::Vector3d(pt_midpoint - pt_old));
                        pt_new = pt_new + pt_tmp;
                        nbContrib++;
                    }
                }

                pt_new = (1. / (double) nbContrib) * pt_new;
                if ((*varNodeGeomAssociation)[nid] != reinterpret_cast<std::uintptr_t>(nullptr)) {
                    reinterpret_cast<gmds::cad::GeomEntity *>((*varNodeGeomAssociation)[nid])->project(pt_new);
                }


                mesh->setNodeLocation(id0 + offset, pt_new.X(), pt_new.Y(), pt_new.Z());

//                mesh->setNodeLocation(id0 + offset, xyz[0], xyz[1], xyz[2]);

                (*node2newNodes)[mat].insert(nid, id0 + offset);
                (*varNodeGeomAssociation)[id0 + offset] = (*varNodeGeomAssociation)[nid];
                offset++;
            }

        }

    };

    /*----------------------------------------------------------------------------*/
    struct Pillow_createNodes_3D {
        const kmds::GrowingView<kmds::TCellID> *selection;
        kmds::Mesh *mesh;
        const kmds::Connectivity *c_N2C;
        const elg3d::MaterialAssignment *ma;
        kmds::Variable<std::uintptr_t> *varNodeGeomAssociation;
        std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> > *node2newNodes;


        Pillow_createNodes_3D(

                const kmds::GrowingView<kmds::TCellID> *selection_,
                kmds::Mesh *mesh_,
                const kmds::Connectivity *c_N2C_,
                const elg3d::MaterialAssignment *ma_,
                kmds::Variable<std::uintptr_t> *AVarNodeGeomAssociation_,
                std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> > *node2newNodes_
        )
                : selection(selection_), mesh(mesh_), c_N2C(c_N2C_), ma(ma_),
                  varNodeGeomAssociation(AVarNodeGeomAssociation_), node2newNodes(node2newNodes_) {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int nid = selection->get(i);

            Kokkos::View<kmds::TCellID *> cells;
            c_N2C->get(nid, cells);

            const int nbCells = cells.size();


            // preparing data
            std::vector<std::vector<kmds::FakeFace> > ffs(nbCells);

            std::vector<std::vector<bool> > ffIsOutward(nbCells);
            std::vector<std::vector<gmds::math::Vector> > ffNormals(nbCells);
            std::vector<std::vector<gmds::math::Triangle> > ffTriangles(nbCells);

            std::vector<std::vector<kmds::FakeFace> > ffs_internal(nbCells);
            std::map<kmds::FakeFace, std::pair<int, int> > internalFF;
            std::vector<std::vector<bool> > ffBoundary(nbCells);

            for(int icell=0; icell<nbCells; icell++) {

                kmds::Region r = mesh->getRegion(cells[icell]);

                std::vector<kmds::FakeFace> fakeFaces = r.getFakeFaces();

                ffs[icell] = fakeFaces;

                ffBoundary[icell] = std::vector<bool>(fakeFaces.size(), true);


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
                }
            }



            std::set<int> materials;
            for (int i_c = 0; i_c < cells.size(); i_c++) {
                materials.insert(ma->getMaterial(cells(i_c)));
            }


            kmds::TCellID id0 = mesh->addNodes(materials.size());
            gmds::math::Point pt_old = mesh->getNodeLocation(nid);

            int offset = 0;
            for (auto mat: materials) {

                double length = HUGE_VALF;

                // we build this set of regions FakeFaces and we compute how many there are considering unicity
                std::map<kmds::FakeFace, int> nbFakeFaces;
                std::map<kmds::FakeFace, std::pair<bool, gmds::math::Vector> > pair_orientation_normals;
                std::map<kmds::FakeFace, gmds::math::Triangle> triangles;

                for(int ic=0; ic<nbCells; ic++) {

                    kmds::TCellID cellid = cells[ic];

                    if(ma->getMaterial(cellid) == mat) {

                        gmds::math::Point pt_midpoint = mesh->getRegion(cellid).midpoint();
                        length = std::min(length, pt_old.distance(pt_midpoint));


                        for (int iff = 0; iff < ffs[ic].size(); iff++) {

                            if (ffs[ic][iff].hasNode(nid)) {

                                if (!ffBoundary[ic][iff]) {

                                    nbFakeFaces[ffs[ic][iff]]++;
                                    pair_orientation_normals[ffs[ic][iff]] = std::pair<bool, gmds::math::Vector>(
                                            ffIsOutward[ic][iff], ffNormals[ic][iff]);
                                    triangles[ffs[ic][iff]] = ffTriangles[ic][iff];
                                }
                            }
                        }
                    }
                }


                // we keep the faces that appeared only once (non-internal faces)
                std::vector<gmds::math::Vector> normals_kept;
                std::vector<gmds::math::Triangle> triangles_kept;
                for (auto ff: nbFakeFaces) {
                    if (ff.second < 2) {

                        //WARNING : remember that the normals are computed at ANodeID of the fakeface
                        if(pair_orientation_normals[ff.first].first) {
                            normals_kept.push_back(pair_orientation_normals[ff.first].second);
                        } else {
                            normals_kept.push_back((-1) * pair_orientation_normals[ff.first].second);
                        }
                        triangles_kept.push_back(triangles[ff.first]);
                    }
                }

                gmds::math::Vector v;
                double area = 0.;

                for(int i=0; i<normals_kept.size(); i++) {
//            double weight = triangles_kept[i].area();
                    double weight = triangles_kept[i].angle();

                    v = v + weight * normals_kept[i];
                    area += weight;
                }

                v.normalize();

                v = -1. * (length / 5.) * v;

                // compute the new node position
                gmds::math::Point pt_new(pt_old + v);

                if((*varNodeGeomAssociation)[nid] != reinterpret_cast<std::uintptr_t>(nullptr)) {
                    reinterpret_cast<gmds::cad::GeomEntity *>((*varNodeGeomAssociation)[nid])->project(pt_new);
                }

                mesh->setNodeLocation(id0 + offset, pt_new.X(), pt_new.Y(), pt_new.Z());

                (*node2newNodes)[mat].insert(nid, id0 + offset);
                (*varNodeGeomAssociation)[id0 + offset] = (*varNodeGeomAssociation)[nid];
                offset++;
            }

        }

    };

    /*----------------------------------------------------------------------------*/
    struct Pillow_createTwoCells_2D {
        const kmds::GrowingView<kmds::TCellID> *selection;
        kmds::Mesh *mesh;
        const kmds::Connectivity *c_F2C;
        elg3d::MaterialAssignment *ma;
        const std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> > *node2newNodes;


        Pillow_createTwoCells_2D(
                const kmds::GrowingView<kmds::TCellID> *selection_,
                kmds::Mesh *mesh_,
                const kmds::Connectivity *c_F2C_,
                elg3d::MaterialAssignment *ma_,
                std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> > *node2newNodes_
        )
                : selection(selection_), mesh(mesh_), c_F2C(c_F2C_), ma(ma_), node2newNodes(node2newNodes_) {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int fid = selection->get(i);

            Kokkos::View<kmds::TCellID *> cells;
            c_F2C->get(fid, cells);

            kmds::Edge f = mesh->getEdge(fid);
            Kokkos::View<kmds::TCellID *> nodes;
            f.nodeIds(nodes);


            // create two cells, which are always quads in 2D
            kmds::TCellID cid0_new = mesh->addQuads(2);


            // determine which new nodes will be used for which cell
            int mat0 = ma->getMaterial(cells(0));
            int mat1 = ma->getMaterial(cells(1));

            kmds::Face c0 = mesh->getFace(cells(0));

            bool isOutward0 = c0.isEdgeOrientedOutward(nodes);

            kmds::TCellID newIDs_0[4];
            kmds::TCellID newIDs_1[4];

            if (isOutward0) {
                newIDs_0[0] = nodes(0);
                newIDs_0[1] = nodes(1);
                newIDs_0[2] = (*node2newNodes)[mat0].value_at((*node2newNodes)[mat0].find(nodes(1)));
                newIDs_0[3] = (*node2newNodes)[mat0].value_at((*node2newNodes)[mat0].find(nodes(0)));

                newIDs_1[0] = nodes(1);
                newIDs_1[1] = nodes(0);
                newIDs_1[2] = (*node2newNodes)[mat1].value_at((*node2newNodes)[mat1].find(nodes(0)));
                newIDs_1[3] = (*node2newNodes)[mat1].value_at((*node2newNodes)[mat1].find(nodes(1)));

            } else {
                newIDs_0[0] = nodes(0);
                newIDs_0[1] = (*node2newNodes)[mat0].value_at((*node2newNodes)[mat0].find(nodes(0)));
                newIDs_0[2] = (*node2newNodes)[mat0].value_at((*node2newNodes)[mat0].find(nodes(1)));
                newIDs_0[3] = nodes(1);

                newIDs_1[0] = nodes(0);
                newIDs_1[1] = nodes(1);
                newIDs_1[2] = (*node2newNodes)[mat1].value_at((*node2newNodes)[mat1].find(nodes(1)));
                newIDs_1[3] = (*node2newNodes)[mat1].value_at((*node2newNodes)[mat1].find(nodes(0)));

            }

            kmds::Face c0_new = mesh->getFace(cid0_new);
            kmds::Face c1_new = mesh->getFace(cid0_new + 1);
            c0_new.setNodes(newIDs_0, 4);
            c1_new.setNodes(newIDs_1, 4);

            ma->setMaterial(mat0, cid0_new);
            ma->setMaterial(mat1, cid0_new + 1);
        }

    };

    /*----------------------------------------------------------------------------*/
    struct Pillow_createTwoCells_3D {
        const kmds::GrowingView<kmds::TCellID> *selection;
        kmds::Mesh *mesh;
        const kmds::Connectivity *c_F2C;
        elg3d::MaterialAssignment *ma;
        const std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> > *node2newNodes;


        Pillow_createTwoCells_3D(
                const kmds::GrowingView<kmds::TCellID> *selection_,
                kmds::Mesh *mesh_,
                const kmds::Connectivity *c_F2C_,
                elg3d::MaterialAssignment *ma_,
                std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> > *node2newNodes_
        )
                : selection(selection_), mesh(mesh_), c_F2C(c_F2C_), ma(ma_), node2newNodes(node2newNodes_) {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int fid = selection->get(i);

            Kokkos::View<kmds::TCellID *> cells;
            c_F2C->get(fid, cells);

            kmds::Face f = mesh->getFace(fid);
            Kokkos::View<kmds::TCellID *> nodes;
            f.nodeIds(nodes);


            // create two cells
            kmds::TCellID cid0_new;
            if (nodes.size() == 4) {
                cid0_new = mesh->addHexahedra(2);
            } else {
                cid0_new = mesh->addPrism3s(2);
            }


            // determine which new nodes will be used for which cell
            int mat0 = ma->getMaterial(cells(0));
            int mat1 = ma->getMaterial(cells(1));

            kmds::Region c0 = mesh->getRegion(cells(0));

            bool isOutward0 = c0.isFaceOrientedOutward(nodes);

            kmds::TCellID newIDs_0[8];
            kmds::TCellID newIDs_1[8];

            if (isOutward0) {
                for (int i_n = 0; i_n < nodes.size(); i_n++) {
                    newIDs_0[i_n] = (*node2newNodes)[mat0].value_at((*node2newNodes)[mat0].find(nodes(i_n)));
                    newIDs_0[i_n + nodes.size()] = nodes(i_n);
                    newIDs_1[i_n] = nodes(i_n);
                    newIDs_1[i_n + nodes.size()] = (*node2newNodes)[mat1].value_at(
                            (*node2newNodes)[mat1].find(nodes(i_n)));
                }
            } else {
                for (int i_n = 0; i_n < nodes.size(); i_n++) {
                    newIDs_0[i_n] = nodes(i_n);
                    newIDs_0[i_n + nodes.size()] = (*node2newNodes)[mat0].value_at(
                            (*node2newNodes)[mat0].find(nodes(i_n)));
                    newIDs_1[i_n] = (*node2newNodes)[mat1].value_at((*node2newNodes)[mat1].find(nodes(i_n)));
                    newIDs_1[i_n + nodes.size()] = nodes(i_n);
                }
            }

            kmds::Region c0_new = mesh->getRegion(cid0_new);
            kmds::Region c1_new = mesh->getRegion(cid0_new + 1);
            c0_new.setNodes(newIDs_0, nodes.size() * 2);
            c1_new.setNodes(newIDs_1, nodes.size() * 2);

            ma->setMaterial(mat0, cid0_new);
            ma->setMaterial(mat1, cid0_new + 1);
        }

    };

    /*----------------------------------------------------------------------------*/
    struct Pillow_replaceInCells_2D {
        const kmds::GrowingView<kmds::TCellID> *selection;
        kmds::Mesh *mesh;
        const elg3d::MaterialAssignment *ma;
        const kmds::Variable<bool> *varIsInterfaceNode;
        const std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> > *node2newNodes;


        Pillow_replaceInCells_2D(
                const kmds::GrowingView<kmds::TCellID> *selection_,
                kmds::Mesh *mesh_,
                const elg3d::MaterialAssignment *ma_,
                const kmds::Variable<bool> *varIsInterfaceNode_,
                const std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> > *node2newNodes_
        )
                : selection(selection_), mesh(mesh_), ma(ma_), varIsInterfaceNode(varIsInterfaceNode_),
                  node2newNodes(node2newNodes_) {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int cid = selection->get(i);

            kmds::Face c = mesh->getFace(cid);
            Kokkos::View<kmds::TCellID *> nodes;
            c.nodeIds(nodes);

            for (int i_n = 0; i_n < nodes.size(); i_n++) {

                if ((*varIsInterfaceNode)[nodes(i_n)]) {
                    int mat = ma->getMaterial(cid);
                    nodes(i_n) = (*node2newNodes)[mat].value_at((*node2newNodes)[mat].find(nodes(i_n)));
                }
            }

            c.setNodes(nodes);
        }

    };

    /*----------------------------------------------------------------------------*/
    struct Pillow_replaceInCells_3D {
        const kmds::GrowingView<kmds::TCellID> *selection;
        kmds::Mesh *mesh;
        const elg3d::MaterialAssignment *ma;
        const kmds::Variable<bool> *varIsInterfaceNode;
        const std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> > *node2newNodes;


        Pillow_replaceInCells_3D(
                const kmds::GrowingView<kmds::TCellID> *selection_,
                kmds::Mesh *mesh_,
                const elg3d::MaterialAssignment *ma_,
                const kmds::Variable<bool> *varIsInterfaceNode_,
                const std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> > *node2newNodes_
        )
                : selection(selection_), mesh(mesh_), ma(ma_), varIsInterfaceNode(varIsInterfaceNode_),
                  node2newNodes(node2newNodes_) {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int rid = selection->get(i);

            kmds::Region c = mesh->getRegion(rid);
            Kokkos::View<kmds::TCellID *> nodes;
            c.nodeIds(nodes);

            for (int i_n = 0; i_n < nodes.size(); i_n++) {

                if ((*varIsInterfaceNode)[nodes(i_n)]) {
                    int mat = ma->getMaterial(rid);
                    nodes(i_n) = (*node2newNodes)[mat].value_at((*node2newNodes)[mat].find(nodes(i_n)));
                }
            }

            c.setNodes(nodes);
        }

    };

    /*----------------------------------------------------------------------------*/
    void
    pillow_2D(const kmds::GrowingView<kmds::TCellID> *AInterfaceEdges,
              const kmds::GrowingView<kmds::TCellID> *AInterfaceNodes,
              kmds::Mesh *AMesh,
              const kmds::Connectivity *c_F2C,
              const kmds::Connectivity *c_N2C,
              elg3d::MaterialAssignment *Ama,
              kmds::Variable<std::uintptr_t> *AVarNodeGeomAssociation) {

        // necessary init
        kmds::Variable<bool> *varIsInterfaceNode = AMesh->createVariable<bool>(false, kmds::KMDS_NODE,
                                                                               "pillow_varIsInterfaceNode");
        Kokkos::parallel_for(AInterfaceNodes->getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 (*varIsInterfaceNode)[AInterfaceNodes->get(i)] = true;
                             });

        AMesh->updateNodeCapacity(AMesh->getNbNodes() + pillow_MAXNBMATPERNODE * AInterfaceNodes->getNbElems());
        AMesh->updateFaceCapacity(AMesh->getNbFaces() + 2 * AInterfaceEdges->getNbElems());

        Ama->updateCapacity(AMesh->getNbFaces() + 2 * AInterfaceEdges->getNbElems());

        std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> > node2newNodes;
        for (int imat = 0; imat < Ama->getNbMaterials(); imat++) {
            node2newNodes.push_back(Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID>(AInterfaceNodes->getNbElems()));
        }


        // for each interface node create one for each mat (valid because each material is manifold)
        Kokkos::parallel_for(AInterfaceNodes->getNbElems(), Pillow_createNodes_2D(AInterfaceNodes,
                                                                                  AMesh,
                                                                                  c_N2C,
                                                                                  Ama,
                                                                                  AVarNodeGeomAssociation,
                                                                                  &node2newNodes));


        // WARNING we need to get this cell container before creating the new cells
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", AMesh->getNbFaces());
        AMesh->getFaceIDs_dummy(&cellIDs);


        //
        Kokkos::parallel_for(AInterfaceEdges->getNbElems(), Pillow_createTwoCells_2D(AInterfaceEdges,
                                                                                     AMesh,
                                                                                     c_F2C,
                                                                                     Ama,
                                                                                     &node2newNodes));


        //
// TODO get regions container, limit it to regions adjacent to interface nodes
        Kokkos::parallel_for(cellIDs.getNbElems(), Pillow_replaceInCells_2D(&cellIDs,
                                                                            AMesh,
                                                                            Ama,
                                                                            varIsInterfaceNode,
                                                                            &node2newNodes));


        // clean-up temporary data
        AMesh->deleteVariable(kmds::KMDS_NODE, "pillow_varIsInterfaceNode");

    }

    /*----------------------------------------------------------------------------*/
    void
    pillow_3D(const kmds::GrowingView<kmds::TCellID> *AInterfaceFaces,
              const kmds::GrowingView<kmds::TCellID> *AInterfaceNodes,
              kmds::Mesh *AMesh,
              const kmds::Connectivity *c_F2C,
              const kmds::Connectivity *c_N2C,
              elg3d::MaterialAssignment *Ama,
              kmds::Variable<std::uintptr_t> *AVarNodeGeomAssociation) {

        // necessary init
        kmds::Variable<bool> *varIsInterfaceNode = AMesh->createVariable<bool>(false, kmds::KMDS_NODE,
                                                                               "pillow_varIsInterfaceNode");
        Kokkos::parallel_for(AInterfaceNodes->getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 (*varIsInterfaceNode)[AInterfaceNodes->get(i)] = true;
                             });

        AMesh->updateNodeCapacity(AMesh->getNbNodes() + pillow_MAXNBMATPERNODE * AInterfaceNodes->getNbElems());
        AMesh->updateRegionCapacity(AMesh->getNbRegions() + 2 * AInterfaceFaces->getNbElems());

        Ama->updateCapacity(AMesh->getNbRegions() + 2 * AInterfaceFaces->getNbElems());

        std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> > node2newNodes;
        for (int imat = 0; imat < Ama->getNbMaterials(); imat++) {
            node2newNodes.push_back(Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID>(AInterfaceNodes->getNbElems()));
        }


        // for each interface node create one for each mat (valid because each material is manifold)
        Kokkos::parallel_for(AInterfaceNodes->getNbElems(), Pillow_createNodes_3D(AInterfaceNodes,
                                                                                  AMesh,
                                                                                  c_N2C,
                                                                                  Ama,
                                                                                  AVarNodeGeomAssociation,
                                                                                  &node2newNodes));


        // we need to get this cell container before creating the new cells
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", AMesh->getNbRegions());
        AMesh->getRegionIDs_dummy(&cellIDs);


        //
        Kokkos::parallel_for(AInterfaceFaces->getNbElems(), Pillow_createTwoCells_3D(AInterfaceFaces,
                                                                                     AMesh,
                                                                                     c_F2C,
                                                                                     Ama,
                                                                                     &node2newNodes));


        //
// TODO get regions container, limit it to regions adjacent to interface nodes
        Kokkos::parallel_for(cellIDs.getNbElems(), Pillow_replaceInCells_3D(&cellIDs,
                                                                            AMesh,
                                                                            Ama,
                                                                            varIsInterfaceNode,
                                                                            &node2newNodes));


        // clean-up temporary data
        AMesh->deleteVariable(kmds::KMDS_NODE, "pillow_varIsInterfaceNode");

    }

    /*----------------------------------------------------------------------------*/
    void
    pillow_markNodes_basedondist_xD(const kmds::GrowingView<kmds::TCellID> *AInterfaceNodes,
                                    const kmds::Variable<double> *AVarNodeDist,
                                    kmds::Variable<bool> *AVarIsMarkedNode) {


        Kokkos::parallel_for(AInterfaceNodes->getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID nid = AInterfaceNodes->get(i);

                                 if ((*AVarNodeDist)[nid] > 0.) {
                                     (*AVarIsMarkedNode)[nid] = true;
                                 }
                             });
    }

    /*----------------------------------------------------------------------------*/
    void
    pillow_markFacets_fromcells_xD(const kmds::GrowingView<kmds::TCellID>* AFacetIDs,
                                   kmds::Mesh *AMesh,
                                   const kmds::Connectivity *c_A2C,
                                   const kmds::Variable<bool> *AVarIsMarkedCell,
                                   kmds::Variable<bool> *AVarIsMarkedA,
                                   kmds::Variable<bool> *AVarIsMarkedA_withBoundary)
    {

        Kokkos::parallel_for(AFacetIDs->getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID aid = AFacetIDs->get(i);

                                 Kokkos::View<kmds::TCellID *> cids;
                                 c_A2C->get(aid, cids);

//                                 if(cids.size() == 1) {
//                                     std::cout << "eid " << eid << " " << cids[0] << std::endl;
//                                 } else {
//                                     std::cout << "eid " << eid << " " << cids[0] <<" "<<cids[1]<< std::endl;
//                                 }


                                 if (cids.size() == 1) {
                                     // in case the edge is on the boundary we do not want to pillow along it but we
                                     // will mark it anyway
                                     if ((*AVarIsMarkedCell)[cids[0]]) {
                                         (*AVarIsMarkedA_withBoundary)[aid] = true;
//                                         std::cout << "eid " << eid << " " << cids[0] << std::endl;
                                     }
                                 } else {
                                     if ((*AVarIsMarkedCell)[cids[0]] != (*AVarIsMarkedCell)[cids[1]]) {
                                         (*AVarIsMarkedA)[aid] = true;
                                         (*AVarIsMarkedA_withBoundary)[aid] = true;
//                                         std::cout << "eid " << eid << " " << cids[0] <<" "<<cids[1]<< std::endl;
                                     }
                                 }
                             });

    }

    /*----------------------------------------------------------------------------*/
    void
    pillow_markNodes_completefromtheedges_2D(kmds::Mesh *AMesh,
                                             const kmds::Connectivity *c_N2A,
                                             const kmds::Variable<bool> *AVarIsMarkedEdge,
                                             const kmds::Variable<bool> *AVarIsDuplicateNode,
                                             kmds::Variable<bool> *AVarIsMarkedNode)
    {
        kmds::GrowingView<kmds::TCellID> nodeIDs("NODES", AMesh->getNbNodes());
        AMesh->getNodeIDs(&nodeIDs);

        Kokkos::parallel_for(nodeIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID nid = nodeIDs.get(i);

                                 if (!(*AVarIsDuplicateNode)[nid]) {

                                     Kokkos::View<kmds::TCellID *> eids;
                                     c_N2A->get(nid, eids);

                                     for (int i_e = 0; i_e < eids.size(); i_e++) {
                                         if ((*AVarIsMarkedEdge)[eids[i_e]]) {
                                             (*AVarIsMarkedNode)[nid] = true;

//                                             std::cout<<"markednode "<<nid<<std::endl;

                                             break;
                                         }
                                     }
                                 }
                             });

    }

    /*----------------------------------------------------------------------------*/
    void
    pillow_markNodes_completefromthefacets_propagate_xD(kmds::GrowingView<kmds::TCellID>* ANodeIDs_base,
                                                       kmds::Mesh *AMesh,
                                                       const kmds::Connectivity *c_N2A,
                                                       const kmds::Variable<bool> *AVarIsMarkedA,
                                                       const kmds::Variable<bool> *AVarIsDuplicateNode,
                                                       kmds::Variable<bool> *AVarIsMarkedNode)
    {

        Kokkos::parallel_for(ANodeIDs_base->getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID nid = ANodeIDs_base->get(i);

                                 if (!(*AVarIsDuplicateNode)[nid]) {

                                     Kokkos::View<kmds::TCellID *> aids;
                                     c_N2A->get(nid, aids);

                                     for (int i_a = 0; i_a < aids.size(); i_a++) {
                                         if ((*AVarIsMarkedA)[aids[i_a]]) {
                                             (*AVarIsMarkedNode)[nid] = true;

//                                             std::cout<<"markednode "<<nid<<std::endl;

                                             break;
                                         }
                                     }
                                 }
                             });

    }

    /*----------------------------------------------------------------------------*/
    void
    pillow_markCells_fromnodes_2D(kmds::Mesh *AMesh,
                                  const kmds::Variable<bool> *AVarIsMarkedNode,
                                  const int AImat,
                                  const elg3d::MaterialAssignment *Ama,
                                  kmds::Variable<bool> *AVarIsMarkedCell) {
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", AMesh->getNbFaces());
        AMesh->getFaceIDs(&cellIDs);

        Kokkos::parallel_for(cellIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID cid = cellIDs.get(i);

                                 if (AImat == Ama->getMaterial(cid)) {

                                     kmds::Face c = AMesh->getFace(cid);
                                     Kokkos::View<kmds::TCellID *> nodes;
                                     c.nodeIds(nodes);

                                     for (int i_n = 0; i_n < nodes.size(); i_n++) {
                                         if ((*AVarIsMarkedNode)[nodes(i_n)]) {
                                             (*AVarIsMarkedCell)[cid] = true;
//                                             std::cout<<"CCCCCCCCCCCCCCCC "<<cid<<std::endl;
                                             break;
                                         }
                                     }
                                 }

                             });
    }

    /*----------------------------------------------------------------------------*/
    void
    pillow_markCells_fromnodes_propagate_2D(const kmds::GrowingView<kmds::TCellID>* ACellIDs,
                                            kmds::Mesh *AMesh,
                                            const kmds::Variable<bool> *AVarIsMarkedNode,
                                            const int AImat,
                                            const elg3d::MaterialAssignment *Ama,
                                            kmds::Variable<bool> *AVarIsMarkedCell,
                                            const kmds::Variable<std::pair<kmds::TCellID, int> > *AVarSeedRangeN,
                                            kmds::Variable<std::pair<kmds::TCellID, int> > *AVarSeedRangeC
    ) {
//        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", AMesh->getNbFaces());
//        AMesh->getFaceIDs(&cellIDs);

        Kokkos::parallel_for(ACellIDs->getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID cid = ACellIDs->get(i);

                                 if (!(*AVarIsMarkedCell)[cid]) {

                                     if (AImat == Ama->getMaterial(cid)) {

                                         kmds::Face c = AMesh->getFace(cid);
                                         Kokkos::View<kmds::TCellID *> nodes;
                                         c.nodeIds(nodes);


                                         for (int i_n = 0; i_n < nodes.size(); i_n++) {

                                             kmds::TCellID nid = nodes(i_n);

                                             if ((*AVarIsMarkedNode)[nid]) {

                                                 if ((*AVarSeedRangeN)[nid].second >= 1) {
                                                     (*AVarIsMarkedCell)[cid] = true;

                                                     (*AVarSeedRangeC)[cid] = std::make_pair(
                                                             (*AVarSeedRangeN)[nid].first,
                                                             (*AVarSeedRangeN)[nid].second);

                                                     break;
                                                 }

                                             }
                                         }

                                     }
                                 }

                             });
    }

    /*----------------------------------------------------------------------------*/
    void
    pillow_markCells_fromnodes_propagate_3D(const kmds::GrowingView<kmds::TCellID>* ACellIDs,
                                            kmds::Mesh *AMesh,
                                            const kmds::Variable<bool> *AVarIsMarkedNode,
                                            const int AImat,
                                            const elg3d::MaterialAssignment *Ama,
                                            kmds::Variable<bool> *AVarIsMarkedCell,
                                            const kmds::Variable<std::pair<kmds::TCellID, int> > *AVarSeedRangeN,
                                            kmds::Variable<std::pair<kmds::TCellID, int> > *AVarSeedRangeC)
    {

        Kokkos::parallel_for(ACellIDs->getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID cid = ACellIDs->get(i);

                                 if (!(*AVarIsMarkedCell)[cid]) {

                                     if (AImat == Ama->getMaterial(cid)) {

                                         kmds::Region c = AMesh->getRegion(cid);
                                         Kokkos::View<kmds::TCellID *> nodes;
                                         c.nodeIds(nodes);


                                         for (int i_n = 0; i_n < nodes.size(); i_n++) {

                                             kmds::TCellID nid = nodes(i_n);

                                             if ((*AVarIsMarkedNode)[nid]) {

                                                 if ((*AVarSeedRangeN)[nid].second >= 1) {
                                                     (*AVarIsMarkedCell)[cid] = true;

                                                     (*AVarSeedRangeC)[cid] = std::make_pair(
                                                             (*AVarSeedRangeN)[nid].first,
                                                             (*AVarSeedRangeN)[nid].second);

                                                     break;
                                                 }

                                             }
                                         }

                                     }
                                 }

                             });
    }

    /*----------------------------------------------------------------------------*/
    void
    pillow_markNodes_fromCells_2D(kmds::Mesh *AMesh,
                                  const kmds::Connectivity *c_N2C,
                                  const int AImat,
                                  const elg3d::MaterialAssignment *Ama,
                                  const kmds::Variable<bool> *AVarIsMarkedCell,
                                  const kmds::Variable<bool> *AVarIsDuplicateNode,
                                  kmds::Variable<bool> *AVarIsMarkedNode) {
        kmds::GrowingView<kmds::TCellID> nodeIDs("NODES", AMesh->getNbNodes());
        AMesh->getNodeIDs(&nodeIDs);

        Kokkos::parallel_for(nodeIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID nid = nodeIDs.get(i);

                                 if (!(*AVarIsDuplicateNode)[nid]) {

                                     Kokkos::View<kmds::TCellID *> cids;
                                     c_N2C->get(nid, cids);

                                     for (int i_c = 0; i_c < cids.size(); i_c++) {
                                         if ((*AVarIsMarkedCell)[cids[i_c]]) {
                                             if (AImat == Ama->getMaterial(cids[i_c])) {
                                                 (*AVarIsMarkedNode)[nid] = true;

//                                         std::cout<<"markednode "<<nid<<std::endl;

                                                 break;
                                             }
                                         }
                                     }
                                 }
                             });

    }

    /*----------------------------------------------------------------------------*/
    void
    pillow_markNodes_fromCells_propagate_xD(const kmds::GrowingView<kmds::TCellID>* ANodeIDs_base,
                                            kmds::Mesh *AMesh,
                                            const kmds::Connectivity *c_N2C,
                                            const int AImat,
                                            const elg3d::MaterialAssignment *Ama,
                                            const kmds::Variable<bool> *AVarIsMarkedCell,
                                            const kmds::Variable<bool> *AVarIsDuplicateNode,
                                            kmds::Variable<bool> *AVarIsMarkedNode,
                                            kmds::Variable<std::pair<kmds::TCellID, int> > *AVarSeedRangeN,
                                            const kmds::Variable<std::pair<kmds::TCellID, int> > *AVarSeedRangeC)
    {

        Kokkos::parallel_for(ANodeIDs_base->getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID nid = ANodeIDs_base->get(i);

                                 if (!(*AVarIsDuplicateNode)[nid]) {

                                     if (!(*AVarIsMarkedNode)[nid]) {

                                         Kokkos::View<kmds::TCellID *> cids;
                                         c_N2C->get(nid, cids);

                                         for (int i_c = 0; i_c < cids.size(); i_c++) {

                                             kmds::TCellID cid = cids[i_c];

                                             if ((*AVarIsMarkedCell)[cid]) {
                                                 if (AImat == Ama->getMaterial(cid)) {
                                                     if ((*AVarSeedRangeC)[cid].second >= 1) {
                                                         (*AVarIsMarkedNode)[nid] = true;

                                                         (*AVarSeedRangeN)[nid] = std::make_pair(
                                                                 (*AVarSeedRangeC)[cid].first,
                                                                 (*AVarSeedRangeC)[cid].second - 1);

                                                         break;
                                                     }
                                                 }
                                             }
                                         }
                                     }
                                 }
                             });

    }

    /*----------------------------------------------------------------------------*/
    void
    pillow_createNodes_withnonmanifold_2D(kmds::Mesh *AMesh,
                                          const kmds::Connectivity *c_N2A,
                                          const kmds::Connectivity *c_A2C,
//                                          const int AImat,
//                                          const elg3d::MaterialAssignment* Ama,
                                          const kmds::Variable<bool> *AVarIsMarkedNode,
                                          const kmds::Variable<bool> *AVarIsMarkedEdge,
                                          const kmds::Variable<bool> *AVarIsMarkedEdge_withBoundary,
                                          const kmds::Variable<bool> *AVarIsMarkedCell,
                                          kmds::GrowingView<duplicate_of_node> *ANodes_duplicates,
                                          kmds::Variable<int> *AVarNbDuplicates,
                                          kmds::Variable<int> *AVarDuplicatesFirstIndex,
                                          kmds::Variable<std::uintptr_t> *AVarNodeGeomAssociation,
                                          kmds::Variable<bool> *AVarIsDuplicateNode) {
        kmds::GrowingView<kmds::TCellID> nodeIDs("NODES", AMesh->getNbNodes());
        AMesh->getNodeIDs(&nodeIDs);


        Kokkos::parallel_for(nodeIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
//for(int i=0; i<nodeIDs.getNbElems(); i++) {
                                 const kmds::TCellID nid = nodeIDs.get(i);

                                 if ((*AVarIsMarkedNode)[nid]) {

                                     gmds::math::Point pt_n = AMesh->getNodeLocation(nid);


                                     Kokkos::View<kmds::TCellID *> eids;
                                     c_N2A->get(nid, eids);

                                     // get number of adjacent marked edges, including those on the domain boundary
                                     int nbMarkedEdges = 0;

                                     for (int i_e = 0; i_e < eids.size(); i_e++) {
                                         if ((*AVarIsMarkedEdge_withBoundary)[eids[i_e]]) {
                                             nbMarkedEdges++;
                                         }
                                     }


                                     // the number of areas, ie the number of duplicate nodes this node will spawn
                                     const int nbAreas = nbMarkedEdges / 2;

                                     (*AVarNbDuplicates)[nid] = nbAreas;

                                     const int first_duplicate = ANodes_duplicates->addElems(nbAreas);

                                     (*AVarDuplicatesFirstIndex)[nid] = first_duplicate;

                                     const kmds::TCellID firstID = AMesh->addNodes(nbAreas);


                                     // compute the duplicate node(s) position
                                     std::vector<bool> used_marked_edges(eids.size(), false);

                                     // prevent the future use of the edges not adjacent to at
                                     // least one marked cell
                                     for (int i_e = 0; i_e < eids.size(); i_e++) {
                                         Kokkos::View<kmds::TCellID *> cids;
                                         c_A2C->get(eids[i_e], cids);

                                         if (!(*AVarIsMarkedCell)[cids[0]]) {
                                             if (cids.size() == 1) {
                                                 used_marked_edges[i_e] = true;
                                             } else {
                                                 if (!(*AVarIsMarkedCell)[cids[1]]) {
                                                     used_marked_edges[i_e] = true;
                                                 }
                                             }
                                         }

                                     }

//                                     std::cout << "nid " << nid << " nbMarkedEdges " << nbMarkedEdges << " nbAreas "
//                                               << nbAreas << " eids.size " << eids.size() << std::endl;


                                     for (int i_n = 0; i_n < nbAreas; i_n++) {

                                         // first select an unused marked edge
                                         kmds::TCellID eid0 = kmds::NullID;
                                         for (int i_e = 0; i_e < eids.size(); i_e++) {
                                             if (!used_marked_edges[i_e] &&
                                                 (*AVarIsMarkedEdge_withBoundary)[eids[i_e]]) {
                                                 eid0 = eids[i_e];
                                                 used_marked_edges[i_e] = true;
                                                 break;
                                             }
                                         }


                                         Kokkos::View<kmds::TCellID *> cids0;
                                         c_A2C->get(eid0, cids0);

                                         kmds::TCellID cid0 = kmds::NullID;
                                         if (cids0.size() == 1) {
                                             cid0 = cids0[0];
                                         } else {
                                             if ((*AVarIsMarkedCell)[cids0[0]]) {
                                                 cid0 = cids0[0];
                                             } else {
                                                 cid0 = cids0[1];
                                             }

                                         }

                                         kmds::TCellID eid_current = eid0;
                                         kmds::TCellID eid_next = kmds::NullID;

                                         kmds::TCellID eid1 = kmds::NullID;
                                         kmds::TCellID cid1 = kmds::NullID;

                                         kmds::TCellID cid_current = cid0;

                                         bool ending_edge_found = false;

                                         while (!ending_edge_found) {

                                             // find the next edge
                                             for (int i_e = 0; i_e < eids.size(); i_e++) {

                                                 kmds::TCellID eid = eids[i_e];

                                                 if (!used_marked_edges[i_e]) {

                                                     Kokkos::View<kmds::TCellID *> cids;
                                                     c_A2C->get(eid, cids);

                                                     {
                                                         kmds::Edge e_current = AMesh->getEdge(eid);

                                                         Kokkos::View<kmds::TCellID *> nids_current;
                                                         e_current.nodeIds(nids_current);

//                                                         std::cout << "e_current " << nids_current[0] << " "
//                                                                   << nids_current[1] << std::endl;
                                                     }


                                                     if (cids.size() == 1) {

                                                         if (cids[0] == cid_current) {

                                                             used_marked_edges[i_e] = true;
                                                             ending_edge_found = true;
                                                             eid1 = eid;
                                                             cid1 = cid_current;
                                                             break;

                                                         } else {
                                                             continue;
                                                         }
                                                     } else {

                                                         if (cids[0] == cid_current) {

                                                             if ((*AVarIsMarkedEdge)[eid]) {
                                                                 used_marked_edges[i_e] = true;
                                                                 ending_edge_found = true;
                                                                 eid1 = eid;
                                                                 cid1 = cid_current;
                                                                 break;
                                                             } else {
                                                                 used_marked_edges[i_e] = true;
//                                                                 eid_current = eid;
                                                                 cid_current = cids[1];
                                                                 break;
                                                             }

                                                         } else {
                                                             if (cids[1] == cid_current) {

                                                                 if ((*AVarIsMarkedEdge)[eid]) {
                                                                     used_marked_edges[i_e] = true;
                                                                     ending_edge_found = true;
                                                                     eid1 = eid;
                                                                     cid1 = cid_current;
                                                                     break;
                                                                 } else {
                                                                     used_marked_edges[i_e] = true;
//                                                                 eid_current = eid;
                                                                     cid_current = cids[0];
                                                                     break;
                                                                 }

                                                             } else {
                                                                 continue;
                                                             }

                                                         }

                                                     }

                                                 }
                                             }

                                         }


                                         duplicate_of_node dn = {firstID + i_n, eid0, eid1};
                                         ANodes_duplicates->set(first_duplicate + i_n, dn);

                                         // now to compute the position of the duplicate
                                         kmds::Edge e0 = AMesh->getEdge(eid0);
                                         kmds::Edge e1 = AMesh->getEdge(eid1);

                                         kmds::Face c0 = AMesh->getFace(cid0);
                                         kmds::Face c1 = AMesh->getFace(cid1);

                                         Kokkos::View<kmds::TCellID *> nids0;
                                         e0.nodeIds(nids0);
                                         Kokkos::View<kmds::TCellID *> nids1;
                                         e1.nodeIds(nids1);


                                         gmds::math::Point pt0_0 = AMesh->getNodeLocation(nids0[0]);
                                         gmds::math::Point pt0_1 = AMesh->getNodeLocation(nids0[1]);
                                         gmds::math::Point pt1_0 = AMesh->getNodeLocation(nids1[0]);
                                         gmds::math::Point pt1_1 = AMesh->getNodeLocation(nids1[1]);
                                         gmds::math::Vector3d n0({pt0_0.Y() - pt0_1.Y(), pt0_1.X() - pt0_0.X(), 0.});
                                         gmds::math::Vector3d n1({pt1_0.Y() - pt1_1.Y(), pt1_1.X() - pt1_0.X(), 0.});

                                         n0.normalize();
                                         n1.normalize();


                                         Kokkos::View<kmds::TCellID *> cids0bis;
                                         c_A2C->get(eid0, cids0bis);
//                                         std::cout << "nid " << nid << " EDGE0 " << nids0[0] << " " << nids0[1]
//                                                   << " c0 " << cid0 << " cids0[0] " << cids0bis[0] << std::endl;

                                         if (!c0.isEdgeOrientedOutward(nids0)) {
                                             n0 = (-1) * n0;
                                         }

                                         Kokkos::View<kmds::TCellID *> cids1;
                                         c_A2C->get(eid1, cids1);
//                                         std::cout << "nid " << nid << " EDGE1 " << nids1[0] << " " << nids1[1]
//                                                   << " c1 " << cid1 << " cids1[1] " << cids1[0] << std::endl;


                                         if (!c1.isEdgeOrientedOutward(nids1)) {
                                             n1 = (-1) * n1;
                                         }


                                         // WARNING : hard-coded one fifth ot the length
                                         kmds::TCoord length = 1. / 5. * std::min(e0.surfvol(), e1.surfvol());

                                         // we sum the normals, ignoring the one from a boundary edge, if any
                                         gmds::math::Vector normal;
                                         if ((*AVarIsMarkedEdge)[eid0]) {
                                             normal = normal + n0;
                                         }
                                         if ((*AVarIsMarkedEdge)[eid1]) {
                                             normal = normal + n1;
                                         }
                                         normal.normalize();

                                         gmds::math::Point pt_new = pt_n + length * normal;


                                         // do not forget to project on the associated geometric entity, if any
                                         // and to associate the duplicate to the same entity
                                         if ((*AVarNodeGeomAssociation)[nid] !=
                                             reinterpret_cast<std::uintptr_t>(nullptr)) {
                                             reinterpret_cast<gmds::cad::GeomEntity *>((*AVarNodeGeomAssociation)[nid])->project(
                                                     pt_new);
                                             (*AVarNodeGeomAssociation)[firstID +
                                                                        i_n] = (*AVarNodeGeomAssociation)[nid];
                                         }


                                         AMesh->setNodeLocation(firstID + i_n, pt_new);

                                         (*AVarIsDuplicateNode)[firstID + i_n] = true;


                                     }  // for(int i_n=0; i_n<nbAreas; i_n++)

                                 }
                             });

    }

    /*----------------------------------------------------------------------------*/
    void
    pillow_createNodes_withnonmanifold_3D(kmds::Mesh *AMesh,
                                          const kmds::Connectivity *c_N2A,
                                          const kmds::Connectivity *c_A2C,
//                                          const int AImat,
//                                          const elg3d::MaterialAssignment* Ama,
                                          const kmds::Variable<bool> *AVarIsMarkedNode,
                                          const kmds::Variable<bool> *AVarIsMarkedA,
                                          const kmds::Variable<bool> *AVarIsMarkedA_withBoundary,
                                          const kmds::Variable<bool> *AVarIsMarkedCell,
                                          kmds::GrowingView<duplicate_of_node_3D> *ANodes_duplicates,
                                          kmds::Variable<int> *AVarNbDuplicates,
                                          kmds::Variable<int> *AVarDuplicatesFirstIndex,
                                          kmds::Variable<std::uintptr_t> *AVarNodeGeomAssociation,
                                          kmds::Variable<bool> *AVarIsDuplicateNode)
    {
        kmds::GrowingView<kmds::TCellID> nodeIDs("NODES", AMesh->getNbNodes());
        AMesh->getNodeIDs(&nodeIDs);


        Kokkos::parallel_for(nodeIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID nid = nodeIDs.get(i);

                                 if ((*AVarIsMarkedNode)[nid]) {

                                     gmds::math::Point pt_n = AMesh->getNodeLocation(nid);


                                     Kokkos::View<kmds::TCellID *> aids;
                                     c_N2A->get(nid, aids);


                                     std::vector<bool> used_marked_faces(aids.size(), false);

                                     // get number of adjacent marked faces, including those on the domain boundary
                                     int nbMarkedFaces = 0;

                                     for (int i_a = 0; i_a < aids.size(); i_a++) {
//                                         if ((*AVarIsMarkedA_withBoundary)[aids[i_a]]) {
                                         if ((*AVarIsMarkedA)[aids[i_a]]) {
                                             nbMarkedFaces++;
                                         } else {

                                             // prevent selecting a non-marked face
                                             used_marked_faces[i_a] = true;
                                         }
                                     }

                                     std::vector<std::vector<kmds::TCellID> > areasF;


                                     std::set<kmds::TCellID > used_faces;


                                     int nbTreatedFaces = 0;
                                     while(nbTreatedFaces < nbMarkedFaces) {

                                         // select an unused marked face
                                         int iFace;

                                         for (int i_a = 0; i_a < aids.size(); i_a++) {

                                             if(!used_marked_faces[i_a]) {
                                                 used_marked_faces[i_a] = true;
                                                 iFace = i_a;
                                                 break;
                                             }
                                         }

                                         const kmds::TCellID aid = aids[iFace];

                                         nbTreatedFaces++;
                                         areasF.push_back(std::vector<kmds::TCellID > (1, aid));

                                         const int area_index = areasF.size() - 1;

                                         Kokkos::View<kmds::TCellID *> cids;
                                         c_A2C->get(aid, cids);


                                         used_faces.insert(aid);

                                         std::set<kmds::TCellID > areaR;

                                         if((*AVarIsMarkedCell)[cids[0]]) {
                                             areaR.insert(cids[0]);
                                         }
                                         if(cids.size() == 2 && (*AVarIsMarkedCell)[cids[1]]) {
                                             areaR.insert(cids[1]);
                                         }

                                         bool addedForR = false;
                                         do{
                                             addedForR = false;

                                             for (int i_a = 0; i_a < aids.size(); i_a++) {

                                                 const kmds::TCellID aid = aids[i_a];

                                                 if(used_faces.find(aid) == used_faces.end()) {

                                                     Kokkos::View<kmds::TCellID *> cids;
                                                     c_A2C->get(aids[i_a], cids);

                                                     if (cids.size() == 1) {
                                                         if ((*AVarIsMarkedCell)[cids[0]]) {
                                                             if (areaR.find(cids[0]) != areaR.end()) {
                                                                 if((*AVarIsMarkedA)[aid]) {
                                                                     used_marked_faces[i_a] = true;
                                                                     used_faces.insert(aid);
                                                                     areasF[area_index].push_back(aid);
                                                                     nbTreatedFaces++;
                                                                     addedForR = true;
                                                                 } else {
                                                                     used_faces.insert(aid);
                                                                 }
                                                             }
                                                         }
                                                     } else {
                                                         if ((*AVarIsMarkedCell)[cids[0]]) {
                                                             if (areaR.find(cids[0]) != areaR.end()) {

                                                                 if ((*AVarIsMarkedCell)[cids[1]]) {
                                                                     if (areaR.find(cids[1]) == areaR.end()) {
                                                                         areaR.insert(cids[1]);
                                                                         used_faces.insert(aid);
                                                                         addedForR = true;
                                                                     }
                                                                 } else {
                                                                     if((*AVarIsMarkedA)[aid]) {
                                                                         used_marked_faces[i_a] = true;
                                                                         used_faces.insert(aid);
                                                                         areasF[area_index].push_back(aid);
                                                                         nbTreatedFaces++;
                                                                         addedForR = true;
                                                                     }
                                                                 }
                                                             } else {
                                                                 if ((*AVarIsMarkedCell)[cids[1]]) {
                                                                     if (areaR.find(cids[1]) != areaR.end()) {
                                                                         areaR.insert(cids[0]);
                                                                         used_faces.insert(aid);
                                                                         addedForR = true;
                                                                     }
                                                                 }
                                                             }


                                                         } else {
                                                             if ((*AVarIsMarkedCell)[cids[1]]) {
                                                                 if (areaR.find(cids[1]) != areaR.end()) {

                                                                     if((*AVarIsMarkedA)[aid]) {
                                                                         used_marked_faces[i_a] = true;
                                                                         used_faces.insert(aid);
                                                                         areasF[area_index].push_back(aid);
                                                                         nbTreatedFaces++;
                                                                         addedForR = true;
                                                                     }
                                                                 }
                                                             }
                                                         }
                                                     }


                                                 }

                                             }


                                         } while (addedForR);

                                     }


                                     // the number of areas, ie the number of duplicate nodes this node will spawn
                                     const int nbAreas = areasF.size();


                                     (*AVarNbDuplicates)[nid] = nbAreas;

                                     const int first_duplicate = ANodes_duplicates->addElems(nbAreas);

                                     (*AVarDuplicatesFirstIndex)[nid] = first_duplicate;

                                     const kmds::TCellID firstID = AMesh->addNodes(nbAreas);


                                     // compute the duplicate node(s) position
                                     for (int i_n = 0; i_n < nbAreas; i_n++) {


                                         duplicate_of_node_3D dn = {firstID + i_n};
                                         ANodes_duplicates->set(first_duplicate + i_n, dn);

                                         // now to compute the position of the duplicate
                                         double surf_min = HUGE_VALF;

                                         gmds::math::Vector3d v({0., 0., 0.});


                                         for(int i_f=0; i_f<areasF[i_n].size(); i_f++) {

                                             const kmds::TCellID aid = areasF[i_n][i_f];

                                             const kmds::Face a = AMesh->getFace(aid);

                                             const double surf = a.surfvol();
                                             if(surf < surf_min) {
                                                 surf_min = surf;
                                             }

                                             kmds::TCellID cid;

                                             Kokkos::View<kmds::TCellID *> cids;
                                             c_A2C->get(aid, cids);

                                             if((*AVarIsMarkedCell)[cids[0]]) {
                                                 cid = cids[0];
                                             } else {
                                                 cid = cids[1];
                                             }


                                             const gmds::math::Triangle tri = AMesh->getFace(aid).getTriangle(nid);
                                             double weight = tri.angle();
                                             gmds::math::Vector n = tri.getNormal();


                                             n.normalize();

                                             const kmds::Region c = AMesh->getRegion(cid);


                                             Kokkos::View<kmds::TCellID *> nids;
                                             a.nodeIds(nids);

                                             if(c.isFaceOrientedOutward(nids)) {
                                                 n = (-1) * n;
                                             }

                                             v = v + weight * n;
                                         }

                                         v.normalize();

                                         // WARNING: hard-coded 1/5 length
                                         double length = (1./5.) * std::sqrt(surf_min);

                                         gmds::math::Point pt_new = pt_n + length * v;

//                                         if(std::isnan(pt_new.X())) {
//                                             std::cout<<"NAN "<<std::endl;
//                                             std::cout<<"pt_n   "<<pt_n<<std::endl;
//                                             std::cout<<"pt_new "<<pt_new<<std::endl;
//                                             std::cout<<"v "<<v<<std::endl;
//                                             std::cout<<"length "<<length<<" surf_min "<<surf_min<<std::endl;
//                                             std::cout<<"nid "<<nid<<" "<<" i_n "<<i_n<<" nbAreas "<<nbAreas<<std::endl;
//
//                                             exit(-1);
//                                         }


                                         // do not forget to project on the associated geometric entity, if any
                                         // and to associate the duplicate to the same entity
                                         if ((*AVarNodeGeomAssociation)[nid] !=
                                             reinterpret_cast<std::uintptr_t>(nullptr)) {
                                             reinterpret_cast<gmds::cad::GeomEntity *>((*AVarNodeGeomAssociation)[nid])->project(
                                                     pt_new);
                                             (*AVarNodeGeomAssociation)[firstID +
                                                                        i_n] = (*AVarNodeGeomAssociation)[nid];
                                         }


                                         AMesh->setNodeLocation(firstID + i_n, pt_new);

                                         (*AVarIsDuplicateNode)[firstID + i_n] = true;


                                     }  // for(int i_n=0; i_n<nbAreas; i_n++)

                                 }
                             });

    }

    /*----------------------------------------------------------------------------*/
    void
    pillow_createCells_withnonmanifold_2D(kmds::Mesh *AMesh,
                                          const kmds::Connectivity *c_A2C,
                                          const kmds::Variable<bool> *AVarIsMarkedEdge,
                                          const kmds::Variable<bool> *AVarIsMarkedCell,
                                          const kmds::GrowingView<duplicate_of_node> *ANodes_duplicates,
                                          const kmds::Variable<int> *AVarNbDuplicates,
                                          const kmds::Variable<int> *AVarDuplicatesFirstIndex,
                                          const int AImat,
                                          elg3d::MaterialAssignment *Ama) {
        kmds::GrowingView<kmds::TCellID> edgeIDs("EDGES", AMesh->getNbEdges());
        AMesh->getEdgeIDs(&edgeIDs);


        Kokkos::parallel_for(edgeIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID eid = edgeIDs.get(i);

                                 if ((*AVarIsMarkedEdge)[eid]) {

                                     // get the cell of the material
                                     Kokkos::View<kmds::TCellID *> cids;
                                     c_A2C->get(eid, cids);

                                     kmds::TCellID cid = kmds::NullID;

                                     if ((*AVarIsMarkedCell)[cids[0]]) {
                                         cid = cids[0];
                                     } else {
                                         cid = cids[1];
                                     }

                                     kmds::Face c = AMesh->getFace(cid);


                                     kmds::Edge e = AMesh->getEdge(eid);


                                     Kokkos::View<kmds::TCellID *> nids;
                                     e.nodeIds(nids);

                                     // find the two duplicate nodes
                                     const kmds::TCellID nid0 = nids[0];
                                     const kmds::TCellID nid1 = nids[1];

                                     kmds::TCellID nid0_copy = kmds::NullID;
                                     kmds::TCellID nid1_copy = kmds::NullID;

                                     for (int i_d = 0; i_d < (*AVarNbDuplicates)[nid0]; i_d++) {

                                         int index = (*AVarDuplicatesFirstIndex)[nid0] + i_d;
                                         if ((ANodes_duplicates->get(index)).eid0 == eid ||
                                             (ANodes_duplicates->get(index)).eid1 == eid) {
                                             nid0_copy = (ANodes_duplicates->get(index)).nid;
                                             break;
                                         }

                                     }

                                     for (int i_d = 0; i_d < (*AVarNbDuplicates)[nid1]; i_d++) {

                                         int index = (*AVarDuplicatesFirstIndex)[nid1] + i_d;
                                         if ((ANodes_duplicates->get(index)).eid0 == eid ||
                                             (ANodes_duplicates->get(index)).eid1 == eid) {
                                             nid1_copy = (ANodes_duplicates->get(index)).nid;
                                             break;
                                         }

                                     }

                                     bool oriented_ok = c.isEdgeOrientedOutward(nids);
                                     kmds::TCellID cid_new = kmds::NullID;
                                     if (oriented_ok) {
                                         cid_new = AMesh->newQuad(nid0, nid1, nid1_copy, nid0_copy);
                                     } else {
                                         cid_new = AMesh->newQuad(nid1, nid0, nid0_copy, nid1_copy);
                                     }

                                     // assign the newly created cell to the material
                                     Ama->setMaterial(AImat, cid_new);
                                 }

                             });
    }

    /*----------------------------------------------------------------------------*/
    void
    pillow_createCells_withnonmanifold_propagate_2D(kmds::Mesh *AMesh,
                                                    const kmds::Connectivity *c_A2C,
                                                    const kmds::Variable<bool> *AVarIsMarkedEdge,
                                                    const kmds::Variable<bool> *AVarIsMarkedCell,
                                                    const kmds::GrowingView<duplicate_of_node> *ANodes_duplicates,
                                                    const kmds::Variable<int> *AVarNbDuplicates,
                                                    const kmds::Variable<int> *AVarDuplicatesFirstIndex)
    {
        kmds::GrowingView<kmds::TCellID> edgeIDs("EDGES", AMesh->getNbEdges());
        AMesh->getEdgeIDs(&edgeIDs);


        Kokkos::parallel_for(edgeIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID eid = edgeIDs.get(i);

                                 if ((*AVarIsMarkedEdge)[eid]) {

                                     // get the cell of the material
                                     Kokkos::View<kmds::TCellID *> cids;
                                     c_A2C->get(eid, cids);

                                     kmds::TCellID cid = kmds::NullID;

                                     if ((*AVarIsMarkedCell)[cids[0]]) {
                                         cid = cids[0];
                                     } else {
                                         cid = cids[1];
                                     }

                                     kmds::Face c = AMesh->getFace(cid);


                                     kmds::Edge e = AMesh->getEdge(eid);


                                     Kokkos::View<kmds::TCellID *> nids;
                                     e.nodeIds(nids);

                                     // find the two duplicate nodes
                                     const kmds::TCellID nid0 = nids[0];
                                     const kmds::TCellID nid1 = nids[1];

                                     kmds::TCellID nid0_copy = kmds::NullID;
                                     kmds::TCellID nid1_copy = kmds::NullID;

                                     for (int i_d = 0; i_d < (*AVarNbDuplicates)[nid0]; i_d++) {

                                         int index = (*AVarDuplicatesFirstIndex)[nid0] + i_d;
                                         if ((ANodes_duplicates->get(index)).eid0 == eid ||
                                             (ANodes_duplicates->get(index)).eid1 == eid) {
                                             nid0_copy = (ANodes_duplicates->get(index)).nid;
                                             break;
                                         }

                                     }

                                     for (int i_d = 0; i_d < (*AVarNbDuplicates)[nid1]; i_d++) {

                                         int index = (*AVarDuplicatesFirstIndex)[nid1] + i_d;
                                         if ((ANodes_duplicates->get(index)).eid0 == eid ||
                                             (ANodes_duplicates->get(index)).eid1 == eid) {
                                             nid1_copy = (ANodes_duplicates->get(index)).nid;
                                             break;
                                         }

                                     }

                                     bool oriented_ok = c.isEdgeOrientedOutward(nids);
                                     kmds::TCellID cid_new = kmds::NullID;
                                     if (oriented_ok) {
                                         cid_new = AMesh->newQuad(nid0, nid1, nid1_copy, nid0_copy);
                                     } else {
                                         cid_new = AMesh->newQuad(nid1, nid0, nid0_copy, nid1_copy);
                                     }

                                 }

                             });
    }

    /*----------------------------------------------------------------------------*/
    void
    pillow_createCells_withnonmanifold_propagate_3D(kmds::Mesh *AMesh,
                                                    const kmds::Connectivity *c_A2C,
                                                    const kmds::Variable<bool> *AVarIsMarkedA,
                                                    const kmds::Variable<bool> *AVarIsMarkedCell,
                                                    const kmds::GrowingView<duplicate_of_node_3D> *ANodes_duplicates,
                                                    const kmds::Variable<int> *AVarNbDuplicates,
                                                    const kmds::Variable<int> *AVarDuplicatesFirstIndex)
    {
        kmds::GrowingView<kmds::TCellID> faceIDs("FACES", AMesh->getNbFaces());
        AMesh->getFaceIDs(&faceIDs);


        Kokkos::parallel_for(faceIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

            const kmds::TCellID aid = faceIDs.get(i);

            if ((*AVarIsMarkedA)[aid]) {

                // get the cell of the material
                Kokkos::View<kmds::TCellID *> cids;
                c_A2C->get(aid, cids);

                kmds::TCellID cid = kmds::NullID;

                if ((*AVarIsMarkedCell)[cids[0]]) {
                    cid = cids[0];
                } else {
                    cid = cids[1];
                }

                const kmds::Region c = AMesh->getRegion(cid);
                const kmds::Face a = AMesh->getFace(aid);

                const gmds::math::Point midpoint = c.midpoint();


                Kokkos::View<kmds::TCellID *> nids;
                a.nodeIds(nids);

                // find the four duplicate nodes
                const kmds::TCellID nid0 = nids[0];
                const kmds::TCellID nid1 = nids[1];
                const kmds::TCellID nid2 = nids[2];
                const kmds::TCellID nid3 = nids[3];

                kmds::TCellID nid0_copy = kmds::NullID;
                kmds::TCellID nid1_copy = kmds::NullID;
                kmds::TCellID nid2_copy = kmds::NullID;
                kmds::TCellID nid3_copy = kmds::NullID;

                double dist0 = HUGE_VALF;
                double dist1 = HUGE_VALF;
                double dist2 = HUGE_VALF;
                double dist3 = HUGE_VALF;

                for (int i_d = 0; i_d < (*AVarNbDuplicates)[nid0]; i_d++) {

                    int index = (*AVarDuplicatesFirstIndex)[nid0] + i_d;

                    if(midpoint.distance(AMesh->getNodeLocation((ANodes_duplicates->get(index)).nid)) < dist0) {
                        nid0_copy = (ANodes_duplicates->get(index)).nid;
                    }
                }
                for (int i_d = 0; i_d < (*AVarNbDuplicates)[nid1]; i_d++) {

                    int index = (*AVarDuplicatesFirstIndex)[nid1] + i_d;

                    if(midpoint.distance(AMesh->getNodeLocation((ANodes_duplicates->get(index)).nid)) < dist1) {
                        nid1_copy = (ANodes_duplicates->get(index)).nid;
                    }
                }
                for (int i_d = 0; i_d < (*AVarNbDuplicates)[nid2]; i_d++) {

                    int index = (*AVarDuplicatesFirstIndex)[nid2] + i_d;

                    if(midpoint.distance(AMesh->getNodeLocation((ANodes_duplicates->get(index)).nid)) < dist2) {
                        nid2_copy = (ANodes_duplicates->get(index)).nid;
                    }
                }
                for (int i_d = 0; i_d < (*AVarNbDuplicates)[nid3]; i_d++) {

                    int index = (*AVarDuplicatesFirstIndex)[nid3] + i_d;

                    if(midpoint.distance(AMesh->getNodeLocation((ANodes_duplicates->get(index)).nid)) < dist3) {
                        nid3_copy = (ANodes_duplicates->get(index)).nid;
                    }
                }


                bool oriented_ok = c.isFaceOrientedOutward(nids);

//                std::cout<<"newface "<<nid0<<" "<<nid1<<" "<<nid2<<" "<<nid3<<std::endl;
//                std::cout<<"copy "<<nid0_copy<<" "<<nid1_copy<<" "<<nid2_copy<<" "<<nid3_copy<<std::endl;
//                std::cout<<"orient "<<oriented_ok<<std::endl;

                kmds::TCellID cid_new = kmds::NullID;
                if (oriented_ok) {
                    cid_new = AMesh->newHexahedron(nid0_copy, nid1_copy, nid2_copy, nid3_copy, nid0, nid1, nid2, nid3);
                } else {
                    cid_new = AMesh->newHexahedron(nid0, nid1, nid2, nid3, nid0_copy, nid1_copy, nid2_copy, nid3_copy);
                }



            }

        });
    }

    /*----------------------------------------------------------------------------*/
    void
    pillow_replaceinCells_withnonmanifold_2D(kmds::Mesh *AMesh,
                                             const kmds::Variable<bool> *AVarIsMarkedCell,
                                             const kmds::GrowingView<duplicate_of_node> *ANodes_duplicates,
                                             const kmds::Variable<int> *AVarNbDuplicates,
                                             const kmds::Variable<int> *AVarDuplicatesFirstIndex) {
        kmds::GrowingView<kmds::TCellID> faceIDs("FACES", AMesh->getNbFaces());
        AMesh->getFaceIDs(&faceIDs);


        Kokkos::parallel_for(faceIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID cid = faceIDs.get(i);

                                 if ((*AVarIsMarkedCell)[cid]) {

                                     kmds::Face c = AMesh->getFace(cid);

                                     Kokkos::View<kmds::TCellID *> nids;
                                     c.nodeIds(nids);


                                     for (int i_n = 0; i_n < nids.size(); i_n++) {

                                         kmds::TCellID nid = nids[i_n];

                                         if ((*AVarDuplicatesFirstIndex)[nid] != kmds::NullID) {

                                             // we have to find the good duplicate
                                             // WARNING : for the time being, based on the geometric position
                                             // in relation with the midpoint of the cell. Should be done
                                             // in another way

                                             gmds::math::Point midpoint = c.midpoint();
                                             kmds::TCoord mindist = HUGE_VALF;


                                             for (int i_d = 0; i_d < (*AVarNbDuplicates)[nid]; i_d++) {

                                                 kmds::TCellID nid_new = ANodes_duplicates->get(
                                                         (*AVarDuplicatesFirstIndex)[nid] + i_d).nid;

                                                 gmds::math::Point pt = AMesh->getNodeLocation(nid_new);

                                                 double dist = midpoint.distance2(pt);

                                                 if (dist < mindist) {
                                                     nids(i_n) = nid_new;
                                                     mindist = dist;
                                                 }
                                             }

                                         }

                                     }


                                     c.setNodes(nids);

                                 }

                             });

    }

    /*----------------------------------------------------------------------------*/
    void
    pillow_replaceinCells_withnonmanifold_3D(kmds::Mesh *AMesh,
                                             const kmds::Variable<bool> *AVarIsMarkedCell,
                                             const kmds::GrowingView<duplicate_of_node_3D> *ANodes_duplicates,
                                             const kmds::Variable<int> *AVarNbDuplicates,
                                             const kmds::Variable<int> *AVarDuplicatesFirstIndex)
    {
        kmds::GrowingView<kmds::TCellID> cellIDs("REGIONS", AMesh->getNbRegions());
        AMesh->getRegionIDs(&cellIDs);


        Kokkos::parallel_for(cellIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

            const kmds::TCellID cid = cellIDs.get(i);

            if ((*AVarIsMarkedCell)[cid]) {

                kmds::Region c = AMesh->getRegion(cid);

                Kokkos::View<kmds::TCellID *> nids;
                c.nodeIds(nids);


                for (int i_n = 0; i_n < nids.size(); i_n++) {

                    kmds::TCellID nid = nids[i_n];

                    if ((*AVarDuplicatesFirstIndex)[nid] != kmds::NullID) {

                        // we have to find the good duplicate
                        // WARNING : for the time being, based on the geometric position
                        // in relation with the midpoint of the cell. Should be done
                        // in another way

                        gmds::math::Point midpoint = c.midpoint();
                        kmds::TCoord mindist = HUGE_VALF;


                        for (int i_d = 0; i_d < (*AVarNbDuplicates)[nid]; i_d++) {

                            kmds::TCellID nid_new = ANodes_duplicates->get(
                                    (*AVarDuplicatesFirstIndex)[nid] + i_d).nid;

                            gmds::math::Point pt = AMesh->getNodeLocation(nid_new);

                            double dist = midpoint.distance2(pt);

                            if (dist < mindist) {
                                nids(i_n) = nid_new;
                                mindist = dist;
                            }
                        }

                    }

                }


                c.setNodes(nids);

            }

        });

    }

    /*----------------------------------------------------------------------------*/
    int
    pillow_seed_feedback_2D(const kmds::GrowingView<kmds::TCellID>* ANodeIDs,
                            kmds::Mesh* AMesh_cavity,
                            const kmds::Variable<double>* AVarQuality_cavity,
                            const kmds::Variable<kmds::TCellID>* AVarCavity2originN,
                            const kmds::Variable<kmds::TCellID>* AVarCavity2originC,
                            kmds::Variable<std::pair<kmds::TCellID, int> >* AVarSeedRangeN,
                            kmds::Variable<std::pair<kmds::TCellID, int> >* AVarSeedRangeC)
    {

        // We decrease the range for the nodes that are seeds to bad quality cells
        kmds::GrowingView<kmds::TCellID> cellIDs_cavity("CELLS_CAVITY", AMesh_cavity->getNbFaces());
        AMesh_cavity->getFaceIDs(&cellIDs_cavity);


        Kokkos::UnorderedMap<kmds::TCellID, void> kmap_seedNodesToDecrease(cellIDs_cavity.getNbElems());


        Kokkos::parallel_for(cellIDs_cavity.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID cid = cellIDs_cavity.get(i);

                                 if((*AVarQuality_cavity)[cid] < Parameters::quality_threshold) {

                                     const kmds::TCellID cid_origin = (*AVarCavity2originC)[cid];

                                     if(cid_origin == kmds::NullID) {

                                         // cell created from an edge
                                         kmds::Face c = AMesh_cavity->getFace(cid);

                                         Kokkos::View<kmds::TCellID *> nids;
                                         c.nodeIds(nids);

                                         for(int i_n = 0; i_n<nids.size(); i_n++) {
                                             const kmds::TCellID nid_origin = (*AVarCavity2originN)[nids[i_n]];

                                             if(nid_origin != kmds::NullID) {
                                                 Kokkos::UnorderedMapInsertResult res = kmap_seedNodesToDecrease.insert(
                                                         (*AVarSeedRangeN)[nid_origin].first);
                                             }
                                         }

                                     } else {
                                         Kokkos::UnorderedMapInsertResult res = kmap_seedNodesToDecrease.insert(
                                                 (*AVarSeedRangeC)[cid_origin].first);
                                     }
                                 }

                             });

        AVarSeedRangeC->initialize();


        Kokkos::parallel_for(kmap_seedNodesToDecrease.capacity(),
                             KOKKOS_LAMBDA(
                                     const int i) {
                                 if(kmap_seedNodesToDecrease.valid_at(i)) {
                                     kmds::TCellID nid = kmap_seedNodesToDecrease.key_at(i);

                                    // it should not get lower than zero because in that case the cell would
                                    // not have been marked with this seed in the first place
                                     if(((*AVarSeedRangeN)[nid].second - 1) == 0) {
                                         AVarSeedRangeN->initialize(nid);
                                     } else {
                                         (*AVarSeedRangeN)[nid].second = (*AVarSeedRangeN)[nid].second - 1;
                                     }
                                 }
                             });



        // remove the seed range for all the nodes except the seeds
        int sumSeedRange = 0;
        Kokkos::parallel_reduce(ANodeIDs->getNbElems(),
                                KOKKOS_LAMBDA(const int i,int &sum) {
                                    kmds::TCellID nid = ANodeIDs->get(i);

                                    if ((*AVarSeedRangeN)[nid].first != nid) {
                                        AVarSeedRangeN->initialize(nid);
                                    } else {
                                        sum += (*AVarSeedRangeN)[nid].second;
                                    }

                                },
                                sumSeedRange);

        return sumSeedRange;
    }

    /*----------------------------------------------------------------------------*/
    int
    pillow_seed_feedback_3D(const kmds::GrowingView<kmds::TCellID>* ANodeIDs,
                            kmds::Mesh* AMesh_cavity,
                            const kmds::Variable<double>* AVarQuality_cavity,
                            const kmds::Variable<kmds::TCellID>* AVarCavity2originN,
                            const kmds::Variable<kmds::TCellID>* AVarCavity2originC,
                            kmds::Variable<std::pair<kmds::TCellID, int> >* AVarSeedRangeN,
                            kmds::Variable<std::pair<kmds::TCellID, int> >* AVarSeedRangeC)
    {

        // We decrease the range for the nodes that are seeds to bad quality cells
        kmds::GrowingView<kmds::TCellID> cellIDs_cavity("CELLS_CAVITY", AMesh_cavity->getNbRegions());
        AMesh_cavity->getRegionIDs(&cellIDs_cavity);


        Kokkos::UnorderedMap<kmds::TCellID, void> kmap_seedNodesToDecrease(cellIDs_cavity.getNbElems());


        Kokkos::parallel_for(cellIDs_cavity.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 const kmds::TCellID cid = cellIDs_cavity.get(i);

                                 if((*AVarQuality_cavity)[cid] < Parameters::quality_threshold) {

                                     const kmds::TCellID cid_origin = (*AVarCavity2originC)[cid];

                                     if(cid_origin == kmds::NullID) {

                                         // cell created from an edge
                                         kmds::Region c = AMesh_cavity->getRegion(cid);

                                         Kokkos::View<kmds::TCellID *> nids;
                                         c.nodeIds(nids);

                                         for(int i_n = 0; i_n<nids.size(); i_n++) {
                                             const kmds::TCellID nid_origin = (*AVarCavity2originN)[nids[i_n]];

                                             if(nid_origin != kmds::NullID) {
                                                 Kokkos::UnorderedMapInsertResult res = kmap_seedNodesToDecrease.insert(
                                                         (*AVarSeedRangeN)[nid_origin].first);
                                             }
                                         }

                                     } else {
                                         Kokkos::UnorderedMapInsertResult res = kmap_seedNodesToDecrease.insert(
                                                 (*AVarSeedRangeC)[cid_origin].first);
                                     }
                                 }

                             });

        AVarSeedRangeC->initialize();


        Kokkos::parallel_for(kmap_seedNodesToDecrease.capacity(),
                             KOKKOS_LAMBDA(
                                     const int i) {
                                 if(kmap_seedNodesToDecrease.valid_at(i)) {
                                     kmds::TCellID nid = kmap_seedNodesToDecrease.key_at(i);

                                     // it should not get lower than zero because in that case the cell would
                                     // not have been marked with this seed in the first place
                                     if(((*AVarSeedRangeN)[nid].second - 1) == 0) {
                                         AVarSeedRangeN->initialize(nid);
                                     } else {
                                         (*AVarSeedRangeN)[nid].second = (*AVarSeedRangeN)[nid].second - 1;
                                     }
                                 }
                             });



        // remove the seed range for all the nodes except the seeds
        int sumSeedRange = 0;
        Kokkos::parallel_reduce(ANodeIDs->getNbElems(),
                                KOKKOS_LAMBDA(const int i,int &sum) {
                                    kmds::TCellID nid = ANodeIDs->get(i);

                                    if ((*AVarSeedRangeN)[nid].first != nid) {
                                        AVarSeedRangeN->initialize(nid);
                                    } else {
                                        sum += (*AVarSeedRangeN)[nid].second;
                                    }

                                },
                                sumSeedRange);

        return sumSeedRange;
    }

    /*----------------------------------------------------------------------------*/
    void
    pillow_execute_2D(
                      const kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,

                      kmds::Mesh* AMesh,
                      const kmds::Connectivity* c_A2C,
                      const kmds::Connectivity* c_N2A,
                      const kmds::Connectivity* c_N2C,
                      elg3d::MaterialAssignment* Ama,
                      kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,

                      kmds::Variable<double>* AVarNodeDist
                      )
    {
        // necessary init
        kmds::Variable<bool>* varIsMarkedNode = AMesh->createVariable<bool> (false, kmds::KMDS_NODE, "tmp_varIsMarkedNode");
        kmds::Variable<bool>* varIsMarkedEdge = AMesh->createVariable<bool> (false, kmds::KMDS_EDGE, "tmp_varIsMarkedEdge");
        kmds::Variable<bool>* varIsMarkedEdge_withBoundary = AMesh->createVariable<bool> (false, kmds::KMDS_EDGE, "tmp_varIsMarkedEdge_withBoundary");
        kmds::Variable<bool>* varIsMarkedCell = AMesh->createVariable<bool> (false, kmds::KMDS_FACE, "tmp_varIsMarkedCell");

        kmds::Variable<bool>* varIsDuplicateNode = AMesh->createVariable<bool> (false, kmds::KMDS_NODE, "tmp_varIsDuplicateNode");

        // WARNING : hard-coded value of 4 used here
        kmds::Variable<int>* varNbDuplicates = AMesh->createVariable<int> (0, kmds::KMDS_NODE, "tmp_varNbDuplicates");
        kmds::Variable<int>* varDuplicatesFirstIndex = AMesh->createVariable<int> (kmds::NullID, kmds::KMDS_NODE, "tmp_varDuplicatesFirstIndex");



        for(int imat=0; imat<Ama->getNbMaterials(); imat++) {
//        for(int imat=1; imat<2; imat++) {

            varIsMarkedNode->setValuesTo(false);
            pillow_markNodes_basedondist_xD(AInterfaceNodes,
//                                            c_N2C,
//                                            imat,
//                                            Ama,
                                            AVarNodeDist,
                                            varIsMarkedNode);


            varIsMarkedCell->setValuesTo(false);

            // loop that marks the cells at distance nb iterations
//            parameters* params = new parameters;
            std::cout<<"pillow_neighbors_dist_mark "<<Parameters::pillow_node_seed_range<<std::endl;
            for(int ineighborhood=0; ineighborhood<Parameters::pillow_node_seed_range; ineighborhood++) {

                // mark the cells from the nodes
                pillow_markCells_fromnodes_2D(AMesh,
                                              varIsMarkedNode,
                                              imat,
                                              Ama,
                                              varIsMarkedCell
                );

                // unmark cells around forbidden nodes
                // at a topological distance of



                // mark the nodes from the cells
                pillow_markNodes_fromCells_2D(AMesh,
                                              c_N2C,
                                              imat,
                                              Ama,
                                              varIsMarkedCell,
                                              varIsDuplicateNode,
                                              varIsMarkedNode

                );

            }

            varIsMarkedEdge->setValuesTo(false);
            varIsMarkedEdge_withBoundary->setValuesTo(false);

            kmds::GrowingView<kmds::TCellID> edgeIDs("EDGES", AMesh->getNbEdges());
            AMesh->getEdgeIDs(&edgeIDs);

            pillow_markFacets_fromcells_xD(&edgeIDs,
                                           AMesh,
                                           c_A2C,
                                           varIsMarkedCell,
                                           varIsMarkedEdge,
                                           varIsMarkedEdge_withBoundary
            );

            varIsMarkedNode->setValuesTo(false);
            pillow_markNodes_completefromtheedges_2D(AMesh,
                                                     c_N2A,
                                                     varIsMarkedEdge,
                                                     varIsDuplicateNode,
                                                     varIsMarkedNode
            );

            kmds::GrowingView<kmds::TCellID> nodeIDs("NODES", AMesh->getNbNodes());
            AMesh->getNodeIDs(&nodeIDs);

            int nbMarkedNodes = 0;
            Kokkos::parallel_reduce(nodeIDs.getNbElems(),
                                    KOKKOS_LAMBDA(const int i, int& sum) {
                                        if((*varIsMarkedNode)[nodeIDs.get(i)]) {
                                            sum += 1;
                                        }

                                    },
                                    nbMarkedNodes);

            kmds::GrowingView<duplicate_of_node> nodes_duplicates("NODES", nbMarkedNodes * 4);

            // WARNING hard-coded value of 4
            AMesh->updateNodeCapacity(AMesh->getNbNodes() + nbMarkedNodes * 4);
//            AMesh->updateNodeCapacity(AMesh->getNbNodes() * 2);


            pillow_createNodes_withnonmanifold_2D(AMesh,
                                                  c_N2A,
                                                  c_A2C,
                                                  varIsMarkedNode,
                                                  varIsMarkedEdge,
                                                  varIsMarkedEdge_withBoundary,
                                                  varIsMarkedCell,
                                                  &nodes_duplicates,
                                                  varNbDuplicates,
                                                  varDuplicatesFirstIndex,
                                                  AVarNodeGeomAssociation,
                                                  varIsDuplicateNode
                                                  );

            std::cout<<"ZZZZ "<<AMesh->getNbNodes()<<std::endl;


            int nbMarkedEdges = 0;
            Kokkos::parallel_reduce(edgeIDs.getNbElems(),
                                    KOKKOS_LAMBDA(const int i, int& sum) {
                                        if((*varIsMarkedEdge)[edgeIDs.get(i)]) {
                                            sum += 1;
                                        }

                                    },
                                    nbMarkedEdges);

            // WARNING hard-coded value of 4
            AMesh->updateFaceCapacity(AMesh->getNbFaces() + nbMarkedEdges * 4);
            Ama->updateCapacity(AMesh->getNbFaces() + nbMarkedEdges * 4);



            pillow_createCells_withnonmanifold_2D(AMesh,
                                                  c_A2C,
                                                  varIsMarkedEdge,
                                                  varIsMarkedCell,
                                                  &nodes_duplicates,
                                                  varNbDuplicates,
                                                  varDuplicatesFirstIndex,
                                                  imat,
                                                  Ama);

            pillow_replaceinCells_withnonmanifold_2D(AMesh,
                                                     varIsMarkedCell,
                                                     &nodes_duplicates,
                                                     varNbDuplicates,
                                                     varDuplicatesFirstIndex);



//            c_N2A->setCapacity(c_N2A->getCapacity() + nbMarkedNodes * 4);
//            c_N2A->setAbsoluteCapacity(c_N2A->getAbsoluteCapacity() + nbMarkedNodes);
//            c_N2C->setCapacity(c_N2C->getCapacity() + nbMarkedNodes * 4);
//            c_N2C->setAbsoluteCapacity(c_N2C->getAbsoluteCapacity() + nbMarkedNodes);

        }


        // clean-up temporary data
        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varIsMarkedNode");
        AMesh->deleteVariable(kmds::KMDS_EDGE, "tmp_varIsMarkedEdge");
        AMesh->deleteVariable(kmds::KMDS_EDGE, "tmp_varIsMarkedEdge_withBoundary");
        AMesh->deleteVariable(kmds::KMDS_FACE, "tmp_varIsMarkedCell");

        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varIsDuplicateNode");

        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varNbDuplicates");
        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varDuplicatesFirstIndex");

    }

    /*----------------------------------------------------------------------------*/
    void
    pillow_execute_propagate_2D(const kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
                                kmds::Mesh* AMesh,
                                const kmds::Connectivity* c_A2C,
                                const kmds::Connectivity* c_N2A,
                                const kmds::Connectivity* c_N2C,
                                elg3d::MaterialAssignment* Ama,
                                kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                                kmds::Variable<double>* AVarNodeDist
    )
    {
        // necessary init
        kmds::Variable<bool>* varIsMarkedNode = AMesh->createVariable<bool> (false, kmds::KMDS_NODE, "tmp_varIsMarkedNode");
        kmds::Variable<bool>* varIsMarkedEdge = AMesh->createVariable<bool> (false, kmds::KMDS_EDGE, "tmp_varIsMarkedEdge");
        kmds::Variable<bool>* varIsMarkedEdge_withBoundary = AMesh->createVariable<bool> (false, kmds::KMDS_EDGE, "tmp_varIsMarkedEdge_withBoundary");
        kmds::Variable<bool>* varIsMarkedCell = AMesh->createVariable<bool> (false, kmds::KMDS_FACE, "tmp_varIsMarkedCell");

        kmds::Variable<bool>* varIsDuplicateNode = AMesh->createVariable<bool> (false, kmds::KMDS_NODE, "tmp_varIsDuplicateNode");

        // WARNING : hard-coded value of 4 used here
//        kmds::Variable<int>* varNbDuplicates = AMesh->createVariable<int> (0, kmds::KMDS_NODE, "tmp_varNbDuplicates");
//        kmds::Variable<int>* varDuplicatesFirstIndex = AMesh->createVariable<int> (kmds::NullID, kmds::KMDS_NODE, "tmp_varDuplicatesFirstIndex");


        kmds::Variable<std::pair<kmds::TCellID, int> >* varSeedRangeN =
                AMesh->createVariable<std::pair<kmds::TCellID, int> > (std::make_pair(kmds::NullID, -1), kmds::KMDS_NODE, "tmp_varSeedRangeN");
        kmds::Variable<std::pair<kmds::TCellID, int> >* varSeedRangeC =
                AMesh->createVariable<std::pair<kmds::TCellID, int> > (std::make_pair(kmds::NullID, -1), kmds::KMDS_FACE, "tmp_varSeedRangeC");



        kmds::GrowingView<kmds::TCellID> nodeIDs_base("NODES", AMesh->getNbNodes());
        AMesh->getNodeIDs(&nodeIDs_base);


        for(int imat=0; imat<Ama->getNbMaterials(); imat++) {
//        for(int imat=0; imat<1; imat++) {

            varIsMarkedNode->setValuesTo(false);
            pillow_markNodes_basedondist_xD(AInterfaceNodes,
//                                            c_N2C,
//                                            imat,
//                                            Ama,
                                            AVarNodeDist,
                                            varIsMarkedNode);


            Kokkos::parallel_for(AInterfaceNodes->getNbElems(),
                                 KOKKOS_LAMBDA(const int i) {

                                     kmds::TCellID nid = AInterfaceNodes->get(i);

                                     if((*varIsMarkedNode)[nid]) {
                                         (*varSeedRangeN)[nid] = std::pair<kmds::TCellID, int> (nid, Parameters::pillow_node_seed_range);
                                     }

                                 });




            bool pillowStop = false;
            int numIter = -1;

            while (!pillowStop) {

                numIter++;

                std::string suffix("_mat_" + std::to_string(imat) + "_iter_" + std::to_string(numIter));
                std::cout<<"mat "<<imat<<" iter "<<numIter<<std::endl;


                // mark the nodes
                varIsMarkedNode->setValuesTo(false);
                Kokkos::parallel_for(AInterfaceNodes->getNbElems(),
                                     KOKKOS_LAMBDA(const int i) {

                                         kmds::TCellID nid = AInterfaceNodes->get(i);

//                                         std::cout<<"seedRangeN "<<nid<<" (*varSeedRangeN)[nid].first "<<(*varSeedRangeN)[nid].first<<" (*varSeedRangeN)[nid].second "<<(*varSeedRangeN)[nid].second<<std::endl;

                                         if((*varSeedRangeN)[nid].second > 0) {
                                             (*varIsMarkedNode)[nid] = true;
                                         }

                                     });


                varIsMarkedCell->setValuesTo(false);

                // loop that marks the cells at distance nb iterations
                std::cout << "pillow_neighbors_dist_mark " << Parameters::pillow_node_seed_range << std::endl;
                for (int ineighborhood = 0; ineighborhood < Parameters::pillow_node_seed_range; ineighborhood++) {


                    kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", AMesh->getNbFaces());
                    AMesh->getFaceIDs(&cellIDs);


                    // mark the cells from the nodes
                    pillow_markCells_fromnodes_propagate_2D(&cellIDs,
                                                            AMesh,
                                                            varIsMarkedNode,
                                                            imat,
                                                            Ama,
                                                            varIsMarkedCell,
                                                            varSeedRangeN,
                                                            varSeedRangeC
                    );


                    int nbMarkedCells = 0;
                    Kokkos::parallel_reduce(cellIDs.getNbElems(),
                                            KOKKOS_LAMBDA(const int i, int &sum) {
                                                if ((*varIsMarkedCell)[cellIDs.get(i)]) {
                                                    sum += 1;
                                                }

                                            },
                                            nbMarkedCells);


                    kmds::VTKWriter<kmds::Mesh> w_origin_mark(*AMesh);
                    w_origin_mark.write("origin_cell_mark_" + std::to_string(ineighborhood), kmds::F);


                    // TODO : unmark cells around forbidden nodes
                    // at a topological distance of



                    // mark the nodes from the cells
                    pillow_markNodes_fromCells_propagate_xD(&nodeIDs_base,
                                                            AMesh,
                                                            c_N2C,
                                                            imat,
                                                            Ama,
                                                            varIsMarkedCell,
                                                            varIsDuplicateNode,
                                                            varIsMarkedNode,
                                                            varSeedRangeN,
                                                            varSeedRangeC
                    );

                }

                varIsMarkedEdge->setValuesTo(false);
                varIsMarkedEdge_withBoundary->setValuesTo(false);

                kmds::GrowingView<kmds::TCellID> edgeIDs("EDGES", AMesh->getNbEdges());
                AMesh->getEdgeIDs(&edgeIDs);

                pillow_markFacets_fromcells_xD(&edgeIDs,
                                               AMesh,
                                               c_A2C,
                                               varIsMarkedCell,
                                               varIsMarkedEdge,
                                               varIsMarkedEdge_withBoundary
                );

                varIsMarkedNode->setValuesTo(false);
                pillow_markNodes_completefromthefacets_propagate_xD(&nodeIDs_base,
                                                                   AMesh,
                                                                   c_N2A,
                                                                   varIsMarkedEdge,
                                                                   varIsDuplicateNode,
                                                                   varIsMarkedNode
                );


                int nbMarkedNodes = 0;
                Kokkos::parallel_reduce(nodeIDs_base.getNbElems(),
                                        KOKKOS_LAMBDA(const int i, int &sum) {
                                            if ((*varIsMarkedNode)[nodeIDs_base.get(i)]) {
                                                sum += 1;
//                                                std::cout<<"marked node "<<nodeIDs.get(i)<<std::endl;
                                            }

                                        },
                                        nbMarkedNodes);


                // we only update the mesh capacity on the first try
//                if(firstTry) {
//
//                    // WARNING hard-coded value of 4
//                    AMesh->updateNodeCapacity(AMesh->getNbNodes() + nbMarkedNodes * 4);
//


                    int nbMarkedEdges = 0;
                    Kokkos::parallel_reduce(edgeIDs.getNbElems(),
                                            KOKKOS_LAMBDA(const int i, int &sum) {
                                                if ((*varIsMarkedEdge)[edgeIDs.get(i)]) {
                                                    sum += 1;
                                                }

                                            },
                                            nbMarkedEdges);
//
//                    // WARNING hard-coded value of 4
//                    AMesh->updateFaceCapacity(AMesh->getNbFaces() + nbMarkedEdges * 4);
//                    Ama->updateCapacity(AMesh->getNbFaces() + nbMarkedEdges * 4);
//
//                    firstTry = false;
//                }

                kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", AMesh->getNbFaces());
                AMesh->getFaceIDs(&cellIDs);

                int nbMarkedCells = 0;
                Kokkos::parallel_reduce(cellIDs.getNbElems(),
                                        KOKKOS_LAMBDA(const int i, int &sum) {
                                            if ((*varIsMarkedCell)[cellIDs.get(i)]) {
                                                sum += 1;
//                                                std::cout<<"marked cell "<<cellIDs.get(i)<<std::endl;
                                            }

                                        },
                                        nbMarkedCells);


                kmds::Mesh mesh_cavity;

                kmds::Variable<double>* varQuality_cavity =
                        mesh_cavity.createVariable<double> (-HUGE_VALF, kmds::KMDS_FACE, "tmp_cavity_varQuality");

                // fill the cavity mesh

                // WARNING : hard-coded value of 4
                mesh_cavity.updateNodeCapacity(nbMarkedCells*4 + nbMarkedNodes*4);
                mesh_cavity.updateFaceCapacity(nbMarkedCells + nbMarkedEdges);

                kmds::Variable<kmds::TCellID >* varCavity2originN =
                        mesh_cavity.createVariable<kmds::TCellID > (kmds::NullID, kmds::KMDS_NODE, "tmp_cavity_varCavity2originN");
                kmds::Variable<kmds::TCellID >* varOrigin2cavityN =
                        AMesh->createVariable<kmds::TCellID > (kmds::NullID, kmds::KMDS_NODE, "tmp_varOrigin2cavityN");
                kmds::Variable<kmds::TCellID >* varCavity2originC =
                        mesh_cavity.createVariable<kmds::TCellID > (kmds::NullID, kmds::KMDS_FACE, "tmp_cavity_varCavity2originC");
                kmds::Variable<kmds::TCellID >* varOrigin2cavityC =
                        AMesh->createVariable<kmds::TCellID > (kmds::NullID, kmds::KMDS_FACE, "tmp_cavity_varOrigin2cavityC");

                kmds::Variable<std::uintptr_t>* varNodeGeomAssociation_cavity = mesh_cavity.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr), kmds::KMDS_NODE, "tmp_geomEntity_cavity");

                kmds::Variable<bool>* varIsMarkedNode_cavity =
                        mesh_cavity.createVariable<bool> (false, kmds::KMDS_NODE, "tmp_varIsMarkedNode_cavity");
                kmds::Variable<bool>* varIsMarkedEdge_cavity =
                        mesh_cavity.createVariable<bool> (false, kmds::KMDS_EDGE, "tmp_varIsMarkedEdge_cavity");
                kmds::Variable<bool>* varIsMarkedEdge_withBoundary_cavity =
                        mesh_cavity.createVariable<bool> (false, kmds::KMDS_EDGE, "tmp_varIsMarkedEdge_withBoundary_cavity");
                kmds::Variable<bool>* varIsMarkedCell_cavity =
                        mesh_cavity.createVariable<bool> (false, kmds::KMDS_FACE, "tmp_varIsMarkedCell_cavity");

                kmds::GrowingView<kmds::TCellID> nodes_boundary_cavity("NODES_BOUNDARY_CAVITY", 1);

                cavity_pillow_create_2D(&nodeIDs_base,
                                        AMesh,
                                        &mesh_cavity,
                                        c_N2C,
                                        AVarNodeGeomAssociation,
                                        varNodeGeomAssociation_cavity,
                                        varIsMarkedNode,
                                        varIsMarkedEdge,
                                        varIsMarkedEdge_withBoundary,
                                        varIsMarkedCell,
                                        varIsMarkedNode_cavity,
                                        varIsMarkedEdge_cavity,
                                        varIsMarkedEdge_withBoundary_cavity,
                                        varIsMarkedCell_cavity,
                                        varCavity2originN,
                                        varOrigin2cavityN,
                                        varCavity2originC,
                                        varOrigin2cavityC,
                                        &nodes_boundary_cavity
                );

                kmds::VTKWriter<kmds::Mesh> w_origin_mark(*AMesh);
                w_origin_mark.write("origin_cavity_" + std::to_string(numIter), kmds::F);


                kmds::VTKWriter<kmds::Mesh> w(mesh_cavity);
                w.write("cavity_create" + suffix, kmds::F);


                // create the pillow
                kmds::ConnectivityHelper ch(&mesh_cavity);

                kmds::Connectivity* c_E2F_cavity = mesh_cavity.getConnectivity(kmds::E2F);
                kmds::Connectivity* c_N2E_cavity = mesh_cavity.getConnectivity(kmds::N2E);



                kmds::GrowingView<duplicate_of_node> nodes_duplicates("NODES", nbMarkedNodes * 4);


                kmds::Variable<bool>* varIsDuplicateNode_cavity = mesh_cavity.createVariable<bool> (false, kmds::KMDS_NODE, "tmp_varIsDuplicateNode_cavity");

                // WARNING : hard-coded value of 4 used here
                kmds::Variable<int>* varNbDuplicates_cavity = mesh_cavity.createVariable<int> (0, kmds::KMDS_NODE, "tmp_varNbDuplicates_cavity");
                kmds::Variable<int>* varDuplicatesFirstIndex_cavity = mesh_cavity.createVariable<int> (kmds::NullID, kmds::KMDS_NODE, "tmp_varDuplicatesFirstIndex_cavity");


                pillow_createNodes_withnonmanifold_2D(&mesh_cavity,
                                                      c_N2E_cavity,
                                                      c_E2F_cavity,
                                                      varIsMarkedNode_cavity,
                                                      varIsMarkedEdge_cavity,
                                                      varIsMarkedEdge_withBoundary_cavity,
                                                      varIsMarkedCell_cavity,
                                                      &nodes_duplicates,
                                                      varNbDuplicates_cavity,
                                                      varDuplicatesFirstIndex_cavity,
                                                      varNodeGeomAssociation_cavity,
                                                      varIsDuplicateNode_cavity
                );

                w.write("cavity_duplicate" + suffix, kmds::F);

                pillow_createCells_withnonmanifold_propagate_2D(&mesh_cavity,
                                                                c_E2F_cavity,
                                                                varIsMarkedEdge_cavity,
                                                                varIsMarkedCell_cavity,
                                                                &nodes_duplicates,
                                                                varNbDuplicates_cavity,
                                                                varDuplicatesFirstIndex_cavity);

                w.write("cavity_cells" + suffix, kmds::F);

                pillow_replaceinCells_withnonmanifold_2D(&mesh_cavity,
                                                         varIsMarkedCell_cavity,
                                                         &nodes_duplicates,
                                                         varNbDuplicates_cavity,
                                                         varDuplicatesFirstIndex_cavity);

                w.write("cavity_replace" + suffix, kmds::F);


                // smooth the cavity
//                mesh_cavity.deleteConnectivity(kmds::N2F);
                kmds::Connectivity* c_N2F_cavity = mesh_cavity.createConnectivity(kmds::N2F);
                ch.buildN2F();

                ch.buildN2N(kmds::F);
                kmds::Connectivity* c_N2N_cavity = mesh_cavity.getConnectivity(kmds::N2N);


                elg3d::smartLaplacian_2D(10,
                                         &mesh_cavity,
                                         c_N2F_cavity,
                                         c_N2N_cavity,
                                         varNodeGeomAssociation_cavity,
                                         &nodes_boundary_cavity);

                double minQual = Tools_computeScaledJacobian_2D(&mesh_cavity, varQuality_cavity);
                w.write("cavity_smooth" + suffix, kmds::F);


                // evaluate whether we have to continue

                if(minQual >= Parameters::quality_threshold) {

                    // the mesh of the cavity is ok, now replace it into the origin mesh
                    pillowStop = true;


                    // update origin mesh capacity
                    const int nbNewNodes = cavity_computeNbNewNodes_xD(&mesh_cavity, varCavity2originN);
                    const int nbNewFaces = cavity_computeNbNewCells_2D(&mesh_cavity, varCavity2originC);
                    AMesh->updateNodeCapacity(AMesh->getNbNodes() + nbNewNodes);
                    AMesh->updateFaceCapacity(AMesh->getNbFaces() + nbNewFaces);

                    Ama->updateCapacity(AMesh->getNbFaces() + nbNewFaces);

                    cavity_pillow_insert_2D(AMesh,
                                            &mesh_cavity,
                                            AVarNodeGeomAssociation,
                                            varNodeGeomAssociation_cavity,
                                            varCavity2originN,
                                            varOrigin2cavityN,
                                            varCavity2originC,
                                            varOrigin2cavityC,
                                            imat,
                                            Ama
                    );

//                    kmds::VTKWriter<kmds::Mesh> w_origin(*AMesh);
//                    w_origin.write("origin_insert" + suffix, kmds::F);
                    elg3d::Tools_write_lite_2D(AMesh, Ama, "origin_insert" + suffix);

                } else {

                    // the mesh of the cavity is not ok, now change the values for the seed nodes and retry.
                    const int sumSeedRange = pillow_seed_feedback_2D(&nodeIDs_base,
                                                                     &mesh_cavity,
                                                                     varQuality_cavity,
                                                                     varCavity2originN,
                                                                     varCavity2originC,
                                                                     varSeedRangeN,
                                                                     varSeedRangeC);

                    std::cout<<"mat "<<imat<<" numIter "<<numIter<<" sumSeedRange "<<sumSeedRange<<std::endl;


                    if(sumSeedRange == 0) {

                        // the next cavity will be empty, so we stop
                        pillowStop = true;

                    }
                }

                mesh_cavity.deleteVariable(kmds::KMDS_FACE, "tmp_cavity_varQuality");

                mesh_cavity.deleteVariable(kmds::KMDS_NODE, "tmp_cavity_varCavity2originN");
                AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varOrigin2cavityN");
                mesh_cavity.deleteVariable(kmds::KMDS_FACE, "tmp_cavity_varCavity2originC");
                AMesh->deleteVariable(kmds::KMDS_FACE, "tmp_cavity_varOrigin2cavityC");

                mesh_cavity.deleteVariable(kmds::KMDS_NODE, "tmp_varIsMarkedNode_cavity");
                mesh_cavity.deleteVariable(kmds::KMDS_EDGE, "tmp_varIsMarkedEdge_cavity");
                mesh_cavity.deleteVariable(kmds::KMDS_EDGE, "tmp_varIsMarkedEdge_withBoundary_cavity");
                mesh_cavity.deleteVariable(kmds::KMDS_FACE, "tmp_varIsMarkedCell_cavity");

                mesh_cavity.deleteVariable(kmds::KMDS_NODE, "tmp_geomEntity_cavity");

                mesh_cavity.deleteConnectivity(kmds::N2N);
                mesh_cavity.deleteConnectivity(kmds::E2F);
                mesh_cavity.deleteConnectivity(kmds::N2E);


            }  // while (!pillowStop)

        }  // for(int imat=0; imat<Ama->getNbMaterials(); imat++)


        // clean-up temporary data
        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varIsMarkedNode");
        AMesh->deleteVariable(kmds::KMDS_EDGE, "tmp_varIsMarkedEdge");
        AMesh->deleteVariable(kmds::KMDS_EDGE, "tmp_varIsMarkedEdge_withBoundary");
        AMesh->deleteVariable(kmds::KMDS_FACE, "tmp_varIsMarkedCell");

        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varIsDuplicateNode");

//        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varNbDuplicates");
//        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varDuplicatesFirstIndex");

    }

    /*----------------------------------------------------------------------------*/
    void
    pillow_execute_propagate_3D(const kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
                                kmds::Mesh* AMesh,
                                const kmds::Connectivity* c_A2C,
                                const kmds::Connectivity* c_N2A,
                                const kmds::Connectivity* c_N2C,
                                elg3d::MaterialAssignment* Ama,
                                kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                                kmds::Variable<double>* AVarNodeDist
    )
    {
        // necessary init
        kmds::Variable<bool>* varIsMarkedNode = AMesh->createVariable<bool> (false, kmds::KMDS_NODE, "tmp_varIsMarkedNode");
        kmds::Variable<bool>* varIsMarkedFace = AMesh->createVariable<bool> (false, kmds::KMDS_FACE, "tmp_varIsMarkedFace");
        kmds::Variable<bool>* varIsMarkedFace_withBoundary = AMesh->createVariable<bool> (false, kmds::KMDS_FACE, "tmp_varIsMarkedFace_withBoundary");
        kmds::Variable<bool>* varIsMarkedCell = AMesh->createVariable<bool> (false, kmds::KMDS_REGION, "tmp_varIsMarkedCell");

        kmds::Variable<bool>* varIsDuplicateNode = AMesh->createVariable<bool> (false, kmds::KMDS_NODE, "tmp_varIsDuplicateNode");


        kmds::Variable<std::pair<kmds::TCellID, int> >* varSeedRangeN =
                AMesh->createVariable<std::pair<kmds::TCellID, int> > (std::make_pair(kmds::NullID, -1), kmds::KMDS_NODE, "tmp_varSeedRangeN");
        kmds::Variable<std::pair<kmds::TCellID, int> >* varSeedRangeC =
                AMesh->createVariable<std::pair<kmds::TCellID, int> > (std::make_pair(kmds::NullID, -1), kmds::KMDS_REGION, "tmp_varSeedRangeC");



        kmds::GrowingView<kmds::TCellID> nodeIDs_base("NODES", AMesh->getNbNodes());
        AMesh->getNodeIDs(&nodeIDs_base);


        for(int imat=0; imat<Ama->getNbMaterials(); imat++) {
//        for(int imat=0; imat<1; imat++) {

            varIsMarkedNode->setValuesTo(false);
            pillow_markNodes_basedondist_xD(AInterfaceNodes,
//                                            c_N2C,
//                                            imat,
//                                            Ama,
                                            AVarNodeDist,
                                            varIsMarkedNode);


            Kokkos::parallel_for(AInterfaceNodes->getNbElems(),
                                 KOKKOS_LAMBDA(const int i) {

                                     kmds::TCellID nid = AInterfaceNodes->get(i);

                                     if((*varIsMarkedNode)[nid]) {
                                         (*varSeedRangeN)[nid] = std::pair<kmds::TCellID, int> (nid, Parameters::pillow_node_seed_range);
                                     }

                                 });




            bool pillowStop = false;
            int numIter = -1;

            while (!pillowStop) {

                numIter++;

                std::string suffix("_mat_" + std::to_string(imat) + "_iter_" + std::to_string(numIter));
                std::cout<<"mat "<<imat<<" iter "<<numIter<<std::endl;


                // mark the nodes
                varIsMarkedNode->setValuesTo(false);
                Kokkos::parallel_for(AInterfaceNodes->getNbElems(),
                                     KOKKOS_LAMBDA(const int i) {

                                         kmds::TCellID nid = AInterfaceNodes->get(i);

//                                         std::cout<<"seedRangeN "<<nid<<" (*varSeedRangeN)[nid].first "<<(*varSeedRangeN)[nid].first<<" (*varSeedRangeN)[nid].second "<<(*varSeedRangeN)[nid].second<<std::endl;

                                         if((*varSeedRangeN)[nid].second > 0) {
                                             (*varIsMarkedNode)[nid] = true;
                                         }

                                     });


                varIsMarkedCell->setValuesTo(false);

                // loop that marks the cells at distance nb iterations
                std::cout << "pillow_neighbors_dist_mark " << Parameters::pillow_node_seed_range << std::endl;
                for (int ineighborhood = 0; ineighborhood < Parameters::pillow_node_seed_range; ineighborhood++) {


                    kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", AMesh->getNbRegions());
                    AMesh->getRegionIDs(&cellIDs);


                    // mark the cells from the nodes
                    pillow_markCells_fromnodes_propagate_3D(&cellIDs,
                                                            AMesh,
                                                            varIsMarkedNode,
                                                            imat,
                                                            Ama,
                                                            varIsMarkedCell,
                                                            varSeedRangeN,
                                                            varSeedRangeC
                    );


                    int nbMarkedCells = 0;
                    Kokkos::parallel_reduce(cellIDs.getNbElems(),
                                            KOKKOS_LAMBDA(const int i, int &sum) {
                                                if ((*varIsMarkedCell)[cellIDs.get(i)]) {
                                                    sum += 1;
//                                                    std::cout<<"markedcell "<<cellIDs.get(i)<<" ineighborhood "<<ineighborhood<<std::endl;
                                                }

                                            },
                                            nbMarkedCells);




                    kmds::VTKWriter<kmds::Mesh> w_origin_mark(*AMesh);
                    w_origin_mark.write("origin_cell_mark_" + std::to_string(ineighborhood), kmds::R);


                    // TODO : unmark cells around forbidden nodes
                    // at a topological distance of



                    // mark the nodes from the cells
                    pillow_markNodes_fromCells_propagate_xD(&nodeIDs_base,
                                                            AMesh,
                                                            c_N2C,
                                                            imat,
                                                            Ama,
                                                            varIsMarkedCell,
                                                            varIsDuplicateNode,
                                                            varIsMarkedNode,
                                                            varSeedRangeN,
                                                            varSeedRangeC
                    );

                }

                varIsMarkedFace->setValuesTo(false);
                varIsMarkedFace_withBoundary->setValuesTo(false);

                kmds::GrowingView<kmds::TCellID> faceIDs("FACES", AMesh->getNbFaces());
                AMesh->getFaceIDs(&faceIDs);

                pillow_markFacets_fromcells_xD(&faceIDs,
                                               AMesh,
                                               c_A2C,
                                               varIsMarkedCell,
                                               varIsMarkedFace,
                                               varIsMarkedFace_withBoundary
                );

                varIsMarkedNode->setValuesTo(false);
                pillow_markNodes_completefromthefacets_propagate_xD(&nodeIDs_base,
                                                                   AMesh,
                                                                   c_N2A,
                                                                   varIsMarkedFace,
                                                                   varIsDuplicateNode,
                                                                   varIsMarkedNode
                );


                int nbMarkedNodes = 0;
                Kokkos::parallel_reduce(nodeIDs_base.getNbElems(),
                                        KOKKOS_LAMBDA(const int i, int &sum) {
                                            if ((*varIsMarkedNode)[nodeIDs_base.get(i)]) {
                                                sum += 1;
                                            }

                                        },
                                        nbMarkedNodes);


//                // we only update the mesh capacity on the first try
////                if(firstTry) {
////
////                    // WARNING hard-coded value of 4
////                    AMesh->updateNodeCapacity(AMesh->getNbNodes() + nbMarkedNodes * 4);
////
//                kmds::GrowingView<kmds::TCellID> edgeIDs("EDGES", AMesh->getNbEdges());
//                AMesh->getEdgeIDs(&edgeIDs);
//
                int nbMarkedFaces = 0;
                Kokkos::parallel_reduce(faceIDs.getNbElems(),
                                        KOKKOS_LAMBDA(const int i, int &sum) {
                                            if ((*varIsMarkedFace)[faceIDs.get(i)]) {
                                                sum += 1;
                                            }

                                        },
                                        nbMarkedFaces);
////
////                    // WARNING hard-coded value of 4
////                    AMesh->updateFaceCapacity(AMesh->getNbFaces() + nbMarkedEdges * 4);
////                    Ama->updateCapacity(AMesh->getNbFaces() + nbMarkedEdges * 4);
////
////                    firstTry = false;
////                }
//
                kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", AMesh->getNbRegions());
                AMesh->getRegionIDs(&cellIDs);

                int nbMarkedCells = 0;
                Kokkos::parallel_reduce(cellIDs.getNbElems(),
                                        KOKKOS_LAMBDA(const int i, int &sum) {
                                            if ((*varIsMarkedCell)[cellIDs.get(i)]) {
                                                sum += 1;
//                                                std::cout<<"marked cell "<<cellIDs.get(i)<<std::endl;
                                            }

                                        },
                                        nbMarkedCells);


                kmds::Mesh mesh_cavity;

                kmds::Variable<double>* varQuality_cavity =
                        mesh_cavity.createVariable<double> (-HUGE_VALF, kmds::KMDS_REGION, "tmp_cavity_varQuality");

                // fill the cavity mesh

                // WARNING : hard-coded value of 8
                mesh_cavity.updateNodeCapacity(nbMarkedCells*8 + nbMarkedNodes*8);
                mesh_cavity.updateRegionCapacity(nbMarkedCells + nbMarkedFaces);

                kmds::Variable<kmds::TCellID >* varCavity2originN =
                        mesh_cavity.createVariable<kmds::TCellID > (kmds::NullID, kmds::KMDS_NODE, "tmp_cavity_varCavity2originN");
                kmds::Variable<kmds::TCellID >* varOrigin2cavityN =
                        AMesh->createVariable<kmds::TCellID > (kmds::NullID, kmds::KMDS_NODE, "tmp_varOrigin2cavityN");
                kmds::Variable<kmds::TCellID >* varCavity2originC =
                        mesh_cavity.createVariable<kmds::TCellID > (kmds::NullID, kmds::KMDS_REGION, "tmp_cavity_varCavity2originC");
                kmds::Variable<kmds::TCellID >* varOrigin2cavityC =
                        AMesh->createVariable<kmds::TCellID > (kmds::NullID, kmds::KMDS_REGION, "tmp_cavity_varOrigin2cavityC");

                kmds::Variable<std::uintptr_t>* varNodeGeomAssociation_cavity = mesh_cavity.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr), kmds::KMDS_NODE, "tmp_geomEntity_cavity");

                kmds::Variable<bool>* varIsMarkedNode_cavity =
                        mesh_cavity.createVariable<bool> (false, kmds::KMDS_NODE, "tmp_varIsMarkedNode_cavity");
                kmds::Variable<bool>* varIsMarkedFace_cavity =
                        mesh_cavity.createVariable<bool> (false, kmds::KMDS_FACE, "tmp_varIsMarkedFace_cavity");
                kmds::Variable<bool>* varIsMarkedFace_withBoundary_cavity =
                        mesh_cavity.createVariable<bool> (false, kmds::KMDS_FACE, "tmp_varIsMarkedFace_withBoundary_cavity");
                kmds::Variable<bool>* varIsMarkedCell_cavity =
                        mesh_cavity.createVariable<bool> (false, kmds::KMDS_REGION, "tmp_varIsMarkedCell_cavity");

                kmds::GrowingView<kmds::TCellID> nodes_boundary_cavity("NODES_BOUNDARY_CAVITY", 1);

                cavity_pillow_create_3D(&nodeIDs_base,
                                        AMesh,
                                        &mesh_cavity,
                                        c_N2C,
                                        AVarNodeGeomAssociation,
                                        varNodeGeomAssociation_cavity,
                                        varIsMarkedNode,
                                        varIsMarkedFace,
                                        varIsMarkedFace_withBoundary,
                                        varIsMarkedCell,
                                        varIsMarkedNode_cavity,
                                        varIsMarkedFace_cavity,
                                        varIsMarkedFace_withBoundary_cavity,
                                        varIsMarkedCell_cavity,
                                        varCavity2originN,
                                        varOrigin2cavityN,
                                        varCavity2originC,
                                        varOrigin2cavityC,
                                        &nodes_boundary_cavity
                );

                kmds::VTKWriter<kmds::Mesh> w_origin_mark(*AMesh);
                w_origin_mark.write("origin_cavity_" + std::to_string(numIter), kmds::R);


                kmds::VTKWriter<kmds::Mesh> w(mesh_cavity);
                w.write("cavity_create" + suffix, kmds::R);


                // create the pillow
                kmds::ConnectivityHelper ch(&mesh_cavity);

                kmds::Connectivity* c_F2R_cavity = mesh_cavity.getConnectivity(kmds::F2R);
                kmds::Connectivity* c_N2F_cavity = mesh_cavity.getConnectivity(kmds::N2F);


                // WARNING : hard-coded value of 4 used here
                kmds::GrowingView<duplicate_of_node_3D> nodes_duplicates("NODES", nbMarkedNodes * 8);


                kmds::Variable<bool>* varIsDuplicateNode_cavity = mesh_cavity.createVariable<bool> (false, kmds::KMDS_NODE, "tmp_varIsDuplicateNode_cavity");


                kmds::Variable<int>* varNbDuplicates_cavity = mesh_cavity.createVariable<int> (0, kmds::KMDS_NODE, "tmp_varNbDuplicates_cavity");
                kmds::Variable<int>* varDuplicatesFirstIndex_cavity = mesh_cavity.createVariable<int> (kmds::NullID, kmds::KMDS_NODE, "tmp_varDuplicatesFirstIndex_cavity");


                pillow_createNodes_withnonmanifold_3D(&mesh_cavity,
                                                      c_N2F_cavity,
                                                      c_F2R_cavity,
                                                      varIsMarkedNode_cavity,
                                                      varIsMarkedFace_cavity,
                                                      varIsMarkedFace_withBoundary_cavity,
                                                      varIsMarkedCell_cavity,
                                                      &nodes_duplicates,
                                                      varNbDuplicates_cavity,
                                                      varDuplicatesFirstIndex_cavity,
                                                      varNodeGeomAssociation_cavity,
                                                      varIsDuplicateNode_cavity
                );

                w.write("cavity_duplicate" + suffix, kmds::R); //exit(-1);

                pillow_createCells_withnonmanifold_propagate_3D(&mesh_cavity,
                                                                c_F2R_cavity,
                                                                varIsMarkedFace_cavity,
                                                                varIsMarkedCell_cavity,
                                                                &nodes_duplicates,
                                                                varNbDuplicates_cavity,
                                                                varDuplicatesFirstIndex_cavity);

                w.write("cavity_cells" + suffix, kmds::R);

                pillow_replaceinCells_withnonmanifold_3D(&mesh_cavity,
                                                         varIsMarkedCell_cavity,
                                                         &nodes_duplicates,
                                                         varNbDuplicates_cavity,
                                                         varDuplicatesFirstIndex_cavity);

                w.write("cavity_replace" + suffix, kmds::R);


                // smooth the cavity
//                mesh_cavity.deleteConnectivity(kmds::N2F);
                kmds::Connectivity* c_N2R_cavity = mesh_cavity.createConnectivity(kmds::N2R);
                ch.buildN2R();

                ch.buildN2N(kmds::R);
                kmds::Connectivity* c_N2N_cavity = mesh_cavity.getConnectivity(kmds::N2N);


//                elg3d::smartLaplacian_3D(10,
//                                         &mesh_cavity,
//                                         c_N2R_cavity,
//                                         c_N2N_cavity,
//                                         varNodeGeomAssociation_cavity,
//                                         &nodes_boundary_cavity);
                elg3d::optimizationSmooth_3D(10,
                                         &mesh_cavity,
                                         c_N2R_cavity,
                                         c_N2N_cavity,
                                         varNodeGeomAssociation_cavity,
                                         &nodes_boundary_cavity);

                double minQual = Tools_computeScaledJacobian_3D(&mesh_cavity, varQuality_cavity);
                w.write("cavity_smooth" + suffix, kmds::R);


                // evaluate whether we have to continue

                if(minQual >= Parameters::quality_threshold) {

                    // the mesh of the cavity is ok, now replace it into the origin mesh
                    pillowStop = true;


                    // update origin mesh capacity
                    const int nbNewNodes = cavity_computeNbNewNodes_xD(&mesh_cavity, varCavity2originN);
                    const int nbNewRegions = cavity_computeNbNewCells_3D(&mesh_cavity, varCavity2originC);
                    AMesh->updateNodeCapacity(AMesh->getNbNodes() + nbNewNodes);
                    AMesh->updateRegionCapacity(AMesh->getNbRegions() + nbNewRegions);

                    Ama->updateCapacity(AMesh->getNbRegions() + nbNewRegions);

                    cavity_pillow_insert_3D(AMesh,
                                            &mesh_cavity,
                                            AVarNodeGeomAssociation,
                                            varNodeGeomAssociation_cavity,
                                            varCavity2originN,
                                            varOrigin2cavityN,
                                            varCavity2originC,
                                            varOrigin2cavityC,
                                            imat,
                                            Ama
                    );

                    kmds::VTKWriter<kmds::Mesh> w_origin(*AMesh);
                    w_origin.write("origin_insert" + suffix, kmds::R);

                } else {

                    // the mesh of the cavity is not ok, now change the values for the seed nodes and retry.
                    const int sumSeedRange = pillow_seed_feedback_3D(&nodeIDs_base,
                                                                     &mesh_cavity,
                                                                     varQuality_cavity,
                                                                     varCavity2originN,
                                                                     varCavity2originC,
                                                                     varSeedRangeN,
                                                                     varSeedRangeC);

                    std::cout<<"mat "<<imat<<" numIter "<<numIter<<" sumSeedRange "<<sumSeedRange<<std::endl;


                    if(sumSeedRange == 0) {

                        // the next cavity will be empty, so we stop
                        pillowStop = true;

                    }
                }

                mesh_cavity.deleteVariable(kmds::KMDS_REGION, "tmp_cavity_varQuality");

                mesh_cavity.deleteVariable(kmds::KMDS_NODE, "tmp_cavity_varCavity2originN");
                AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varOrigin2cavityN");
                mesh_cavity.deleteVariable(kmds::KMDS_REGION, "tmp_cavity_varCavity2originC");
                AMesh->deleteVariable(kmds::KMDS_REGION, "tmp_cavity_varOrigin2cavityC");

                mesh_cavity.deleteVariable(kmds::KMDS_NODE, "tmp_varIsMarkedNode_cavity");
                mesh_cavity.deleteVariable(kmds::KMDS_FACE, "tmp_varIsMarkedFace_cavity");
                mesh_cavity.deleteVariable(kmds::KMDS_FACE, "tmp_varIsMarkedFace_withBoundary_cavity");
                mesh_cavity.deleteVariable(kmds::KMDS_REGION, "tmp_varIsMarkedCell_cavity");

                mesh_cavity.deleteVariable(kmds::KMDS_NODE, "tmp_geomEntity_cavity");

                mesh_cavity.deleteConnectivity(kmds::N2N);
                mesh_cavity.deleteConnectivity(kmds::F2R);
                mesh_cavity.deleteConnectivity(kmds::N2F);


            }  // while (!pillowStop)

        }  // for(int imat=0; imat<Ama->getNbMaterials(); imat++)
//
//
//        // clean-up temporary data
//        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varIsMarkedNode");
//        AMesh->deleteVariable(kmds::KMDS_EDGE, "tmp_varIsMarkedEdge");
//        AMesh->deleteVariable(kmds::KMDS_EDGE, "tmp_varIsMarkedEdge_withBoundary");
//        AMesh->deleteVariable(kmds::KMDS_FACE, "tmp_varIsMarkedCell");
//
//        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varIsDuplicateNode");
//
////        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varNbDuplicates");
////        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varDuplicatesFirstIndex");

    }

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
