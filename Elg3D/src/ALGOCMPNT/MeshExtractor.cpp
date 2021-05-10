/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    MeshExtractor.cpp
 *  \author  legoff
 *  \date    05/22/2019
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/MeshExtractor.h"
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
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/MaterialInterfaces.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {





    /*----------------------------------------------------------------------------*/
    void MeshExtractor_extract_NandF(kmds::Mesh* AMesh_source,
                                     kmds::Mesh* AMesh_target,
                                     const kmds::GrowingView<kmds::TCellID>* ASelectionN,
                                     const kmds::GrowingView<kmds::TCellID>* ASelectionF,
                                     kmds::GrowingView<kmds::TCellID>* ANid_old2new,
                                     kmds::GrowingView<kmds::TCellID>* AFid_old2new)

    {

        const int nbNodes_target = ASelectionN->getNbElems();
        const int nbFaces_target = ASelectionF->getNbElems();

        AMesh_target->updateNodeCapacity(nbNodes_target);
        AMesh_target->updateFaceCapacity(nbFaces_target);

        const kmds::TCellID firstNID = AMesh_target->addNodes(nbNodes_target);

//        ANid_old2new->resize(nbNodes_target);

        kmds::Variable<kmds::TCellID> *tmp_varOld2NewNid =
                AMesh_source->createVariable<kmds::TCellID>(kmds::NullID,
                                                       kmds::KMDS_NODE, "tmp_varOld2NewNid");

        // WARNING : currently only considers quads
        const kmds::TCellID firstFID = AMesh_target->addQuads(nbFaces_target);


        Kokkos::parallel_for(nbNodes_target,
                             KOKKOS_LAMBDA(const int i)
                             {
                                 const kmds::TCellID nid = ASelectionN->get(i);
                                 kmds::Node n = AMesh_source->getNode(nid);

                                 AMesh_target->setNodeLocation(firstNID + i, AMesh_source->getNodeLocation(nid));

                                 ANid_old2new->set(i, firstNID + i);
                                 (*tmp_varOld2NewNid)[nid] = firstNID + i;
                             }
        );


        Kokkos::parallel_for(nbFaces_target,
                             KOKKOS_LAMBDA(const int i)
                             {
                                 const kmds::TCellID fid = ASelectionF->get(i);
                                 kmds::Face f_source = AMesh_source->getFace(fid);

                                 Kokkos::View<kmds::TCellID*> nids;
                                 f_source.nodeIds(nids);

                                 kmds::TCellID nids_target[nids.size()];
                                 for(int i_n=0; i_n<nids.size(); i_n++) {
                                     nids_target[i_n] = (*tmp_varOld2NewNid)[nids[i_n]];
                                 }

                                 AFid_old2new->set(i, firstFID + i);
                                 kmds::Face f_target = AMesh_target->getFace(firstFID + i);
                                 f_target.setNodes(nids_target, nids.size());
                             }
        );


        AMesh_source->deleteVariable(kmds::KMDS_NODE, "tmp_varOld2NewNid");
    }

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
