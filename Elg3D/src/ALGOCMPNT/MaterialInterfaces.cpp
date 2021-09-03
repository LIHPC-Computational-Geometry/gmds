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
/** \file    MaterialInterfaces.cpp
 *  \author  legoff
 *  \date    22/11/2017
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/MaterialInterfaces.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
namespace elg3d {

    struct MaterialInterfaces_NodeIsOnInterface
    {
        kmds::GrowingView<kmds::TCellID>* nodeIDsAccessor;

        const kmds::Connectivity* c_N2C;
        const elg3d::MaterialAssignment* ma;
        kmds::GrowingView<kmds::TCellID>* selection;


        MaterialInterfaces_NodeIsOnInterface(kmds::GrowingView<kmds::TCellID>* nodeIDsAccessor_,
                                        const kmds::Connectivity* c_N2C_,
                                        const elg3d::MaterialAssignment* ma_,
                                        kmds::GrowingView<kmds::TCellID>* Selection_)
                : nodeIDsAccessor(nodeIDsAccessor_)
                , c_N2C(c_N2C_)
                , ma(ma_)
                , selection(Selection_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const
        {
            const int nid = nodeIDsAccessor->get(i);

            Kokkos::View<kmds::TCellID *> cellIDs;
            c_N2C->get(nid, cellIDs);

            // get the material assignement for one cell
            const int matRef = ma->getMaterial(cellIDs(0));

            // check whether some cells have another material assignement
            const int nbCells = cellIDs.size();
            for(int icell=1; icell<nbCells; icell++) {
                if(ma->getMaterial(cellIDs(icell)) != matRef) {
                    selection->push_back(nid);
                    break;
                }
            }
        }
    };

    struct MaterialInterfaces_NodeIsOnInterfaceOneMaterial
    {
        kmds::GrowingView<kmds::TCellID>* nodeIDsAccessor;

        int matIndex;
        const kmds::Connectivity* c_N2C;
        const elg3d::MaterialAssignment* ma;
        kmds::GrowingView<kmds::TCellID>* selection;


        MaterialInterfaces_NodeIsOnInterfaceOneMaterial(kmds::GrowingView<kmds::TCellID>* nodeIDsAccessor_,
                                             int matIndex_,
                                             const kmds::Connectivity* c_N2C_,
                                             const elg3d::MaterialAssignment* ma_,
                                             kmds::GrowingView<kmds::TCellID>* Selection_)
                : nodeIDsAccessor(nodeIDsAccessor_)
                , matIndex(matIndex_)
                , c_N2C(c_N2C_)
                , ma(ma_)
                , selection(Selection_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const
        {
            int nid = nodeIDsAccessor->get(i);

            Kokkos::View<kmds::TCellID *> cellIDs;
            c_N2C->get(nid, cellIDs);

            int nbCells = cellIDs.size();

            // check whether the node is part of material matIndex
            bool nodeIsOfMat = false;
            for(int icell=0; icell<nbCells; icell++) {

                if(matIndex == ma->getMaterial(cellIDs(icell))) {
                    nodeIsOfMat = true;
                    break;
                }
            }

            if(!nodeIsOfMat) {
                return;
            }

            // get the material assignement for one cell
            int matRef = ma->getMaterial(cellIDs(0));

            // check whether some cells have another material assignement
            for(int icell=1; icell<nbCells; icell++) {

                if(ma->getMaterial(cellIDs(icell)) != matRef) {
                    selection->push_back(nid);
                    break;
                }
            }
        }
    };
/*----------------------------------------------------------------------------*/
    void
    MaterialInterfaces_getNodeOnInterfaces(kmds::Mesh* AMesh, const kmds::Connectivity* Ac_N2C, const elg3d::MaterialAssignment* Ama, kmds::GrowingView<kmds::TCellID>* ASelection)
    {
        kmds::GrowingView<kmds::TCellID> nodeIDsAccessor("NODES", AMesh->getNbNodes());
        AMesh->getNodeIDs_dummy(&nodeIDsAccessor);

        Kokkos::parallel_for(nodeIDsAccessor.getNbElems(), MaterialInterfaces_NodeIsOnInterface(&nodeIDsAccessor, Ac_N2C, Ama, ASelection));
    }
/*----------------------------------------------------------------------------*/
    void
    MaterialInterfaces_getNodeOnInterfaces(int AMat, kmds::Mesh* AMesh, const kmds::Connectivity* Ac_N2C, const elg3d::MaterialAssignment* Ama, kmds::GrowingView<kmds::TCellID>* ASelection)
    {
        kmds::GrowingView<kmds::TCellID> nodeIDsAccessor("NODES", AMesh->getNbNodes());
        AMesh->getNodeIDs_dummy(&nodeIDsAccessor);

        Kokkos::parallel_for(nodeIDsAccessor.getNbElems(), MaterialInterfaces_NodeIsOnInterfaceOneMaterial(&nodeIDsAccessor, AMat, Ac_N2C, Ama, ASelection));
    }
/*----------------------------------------------------------------------------*/
    void
    MaterialInterfaces_getEdgeOnInterfaces(kmds::Mesh* AMesh, const kmds::Connectivity* Ac_A2C, const elg3d::MaterialAssignment* Ama, kmds::GrowingView<kmds::TCellID>* ASelection)
    {
        kmds::GrowingView<kmds::TCellID> edgeIDsAccessor("EDGES", AMesh->getNbEdges());
        AMesh->getEdgeIDs_dummy(&edgeIDsAccessor);

        Kokkos::parallel_for(edgeIDsAccessor.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID eid = edgeIDsAccessor.get(i);

                                 Kokkos::View<kmds::TCellID *> cids;
                                 Ac_A2C->get(eid, cids);

                                 // check if it is an interface face
                                 // we ignore boundary faces
                                 if(cids.size() == 2) {

                                     if(Ama->getMaterial(cids(0)) != Ama->getMaterial(cids(1))) {
                                         ASelection->push_back(eid);
                                     }
                                 }
                             });
    }

    /*----------------------------------------------------------------------------*/
    void
    MaterialInterfaces_getEdgeOnInterfacesWithBoundary(kmds::Mesh* AMesh,
                                                       const kmds::Connectivity* Ac_A2C,
                                                       const elg3d::MaterialAssignment* Ama,
                                                       kmds::GrowingView<kmds::TCellID>* ASelection)
    {
        kmds::GrowingView<kmds::TCellID> edgeIDsAccessor("EDGES", AMesh->getNbEdges());
        AMesh->getEdgeIDs_dummy(&edgeIDsAccessor);

        Kokkos::parallel_for(edgeIDsAccessor.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID eid = edgeIDsAccessor.get(i);

                                 Kokkos::View<kmds::TCellID *> cids;
                                 Ac_A2C->get(eid, cids);

                                 // check if it is an interface face
                                 if(cids.size() == 2) {

                                     if(Ama->getMaterial(cids(0)) != Ama->getMaterial(cids(1))) {
                                         ASelection->push_back(eid);
                                     }
                                 } else {
                                     ASelection->push_back(eid);
                                 }
                             });
    }

/*----------------------------------------------------------------------------*/
    void
    MaterialInterfaces_getFaceOnInterfaces(kmds::Mesh* AMesh, const kmds::Connectivity* Ac_A2C, const elg3d::MaterialAssignment* Ama, kmds::GrowingView<kmds::TCellID>* ASelection)
    {
        kmds::GrowingView<kmds::TCellID> faceIDsAccessor("FACES", AMesh->getNbFaces());
        AMesh->getFaceIDs_dummy(&faceIDsAccessor);

        Kokkos::parallel_for(faceIDsAccessor.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID fid = faceIDsAccessor.get(i);

                                 Kokkos::View<kmds::TCellID *> cids;
                                 Ac_A2C->get(fid, cids);

                                 // check if it is an interface face
                                 // we ignore boundary faces
                                 if(cids.size() == 2) {

                                     if(Ama->getMaterial(cids(0)) != Ama->getMaterial(cids(1))) {
                                         ASelection->push_back(fid);
                                     }
                                 };
                             });
    }

    /*----------------------------------------------------------------------------*/
    void
    MaterialInterfaces_getFaceOnInterfacesWithBoundary(kmds::Mesh* AMesh,
                                                       const kmds::Connectivity* Ac_A2C,
                                                       const elg3d::MaterialAssignment* Ama,
                                                       kmds::GrowingView<kmds::TCellID>* ASelection)
    {
        kmds::GrowingView<kmds::TCellID> faceIDsAccessor("FACES", AMesh->getNbFaces());
        AMesh->getFaceIDs_dummy(&faceIDsAccessor);

        Kokkos::parallel_for(faceIDsAccessor.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID fid = faceIDsAccessor.get(i);

                                 Kokkos::View<kmds::TCellID *> cids;
                                 Ac_A2C->get(fid, cids);

                                 // check if it is an interface face
                                 if(cids.size() == 2) {

                                     if(Ama->getMaterial(cids(0)) != Ama->getMaterial(cids(1))) {
                                         ASelection->push_back(fid);
                                     }
                                 } else {
//                                     ASelection->push_back(fid);
                                 }
                             });
    }
/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
