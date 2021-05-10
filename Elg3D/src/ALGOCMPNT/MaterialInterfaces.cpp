/*----------------------------------------------------------------------------*/

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
