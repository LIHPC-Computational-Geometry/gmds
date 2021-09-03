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
