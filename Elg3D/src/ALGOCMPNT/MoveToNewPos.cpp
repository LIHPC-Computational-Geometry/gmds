

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
/** \file    MoveToNewPos.cpp
 *  \author  legoff
 *  \date    22/11/2017
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/MoveToNewPos.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <cmath>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <Kokkos_UnorderedMap.hpp>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <KM/Utils/Graph.h>
#include <gmds/cad/GeomEntity.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/ManifoldDetection.h"
#include "ELG3D/ALGOCMPNT/MaterialInterfaces.h"
#include "ELG3D/ALGOCMPNT/Tools.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {

    /*----------------------------------------------------------------------------*/
    struct moveToNewPos_noBadMove_oneNode_2D
    {
        const kmds::GrowingView<kmds::TCellID>* selection;
        const kmds::TCoord ratioIter;
        const kmds::TCoord qualThreshold;
        kmds::Mesh* mesh;
        const kmds::Connectivity* c_N2C;
        const kmds::Variable<gmds::math::Point>* varNodeOrigPos;
        const kmds::Variable<gmds::math::Point>* varNodeFinalPos;
        const kmds::Variable<std::uintptr_t>* varNodeGeomAssociation;
        kmds::Variable<int>* AVarNodeNbMoves;


        moveToNewPos_noBadMove_oneNode_2D(

                const kmds::GrowingView<kmds::TCellID>* selection_,
                const kmds::TCoord ratioIter_,
                const kmds::TCoord qualThreshold_,
                kmds::Mesh* mesh_,
                const kmds::Connectivity* c_N2C_,
                const kmds::Variable<gmds::math::Point>* varNodeOrigPos_,
                const kmds::Variable<gmds::math::Point>* varNodeFinalPos_,
                const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation_,
                kmds::Variable<int>* AVarNodeNbMoves_
        )
                : selection(selection_)
                , ratioIter(ratioIter_)
                , qualThreshold(qualThreshold_)
                , mesh(mesh_)
                , c_N2C(c_N2C_)
                , varNodeOrigPos(varNodeOrigPos_)
                , varNodeFinalPos(varNodeFinalPos_)
                , varNodeGeomAssociation(AVarNodeGeomAssociation_)
                , AVarNodeNbMoves(AVarNodeNbMoves_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int nid = selection->get(i);

            gmds::math::Point pt_current = mesh->getNodeLocation(nid);

            // compute this iteration destination
            gmds::math::Point pt_next = (1. - ratioIter) * (*varNodeOrigPos)[nid] + ratioIter * (*varNodeFinalPos)[nid];

            if((*varNodeGeomAssociation)[nid] != reinterpret_cast<std::uintptr_t>(nullptr)) {
                reinterpret_cast<gmds::cad::GeomEntity *>((*varNodeGeomAssociation)[nid])->project(pt_next);
            }

            mesh->setNodeLocation(nid, pt_next);

            Kokkos::View<kmds::TCellID*> cells;
            c_N2C->get(nid, cells);

            bool reverted = false;

            for(int i_c=0; i_c<cells.size(); i_c++) {
                kmds::TCellID cid = cells[i_c];

                kmds::Face c = mesh->getFace(cid);

                if(c.scaledJacobian() < qualThreshold) {
                    // revert location
                    mesh->setNodeLocation(nid, pt_current);
                    reverted = true;
                    break;
                }
            }

            if(!reverted) {
                // we count this as an effective move
                (*AVarNodeNbMoves)[nid]++;
            }
        }

    };

    /*----------------------------------------------------------------------------*/
    struct moveToNewPos_noBadMove_oneNode_3D
    {
        const kmds::GrowingView<kmds::TCellID>* selection;
        const kmds::TCoord ratioIter;
        const kmds::TCoord qualThreshold;
        kmds::Mesh* mesh;
        const kmds::Connectivity* c_N2C;
        const kmds::Variable<gmds::math::Point>* varNodeOrigPos;
        const kmds::Variable<gmds::math::Point>* varNodeFinalPos;
        const kmds::Variable<std::uintptr_t>* varNodeGeomAssociation;
        kmds::Variable<int>* AVarNodeNbMoves;


        moveToNewPos_noBadMove_oneNode_3D(

                const kmds::GrowingView<kmds::TCellID>* selection_,
                const kmds::TCoord ratioIter_,
                const kmds::TCoord qualThreshold_,
                kmds::Mesh* mesh_,
                const kmds::Connectivity* c_N2C_,
                const kmds::Variable<gmds::math::Point>* varNodeOrigPos_,
                const kmds::Variable<gmds::math::Point>* varNodeFinalPos_,
                const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation_,
                kmds::Variable<int>* AVarNodeNbMoves_
        )
                : selection(selection_)
                , ratioIter(ratioIter_)
                , qualThreshold(qualThreshold_)
                , mesh(mesh_)
                , c_N2C(c_N2C_)
                , varNodeOrigPos(varNodeOrigPos_)
                , varNodeFinalPos(varNodeFinalPos_)
                , varNodeGeomAssociation(AVarNodeGeomAssociation_)
                , AVarNodeNbMoves(AVarNodeNbMoves_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int nid = selection->get(i);

            gmds::math::Point pt_current = mesh->getNodeLocation(nid);

            // compute this iteration destination
            gmds::math::Point pt_next = (1. - ratioIter) * (*varNodeOrigPos)[nid] + ratioIter * (*varNodeFinalPos)[nid];

            if((*varNodeGeomAssociation)[nid] != reinterpret_cast<std::uintptr_t>(nullptr)) {
                reinterpret_cast<gmds::cad::GeomEntity *>((*varNodeGeomAssociation)[nid])->project(pt_next);
            }

            mesh->setNodeLocation(nid, pt_next);

            Kokkos::View<kmds::TCellID*> cells;
            c_N2C->get(nid, cells);

            bool reverted = false;

            for(int i_c=0; i_c<cells.size(); i_c++) {
                kmds::TCellID cid = cells[i_c];

                kmds::Region c = mesh->getRegion(cid);

                if(c.scaledJacobian() < qualThreshold) {
                    // revert location
                    mesh->setNodeLocation(nid, pt_current);
                    reverted = true;
                    break;
                }
            }

            if(!reverted) {
                // we count this as an effective move
                (*AVarNodeNbMoves)[nid]++;
            }
        }

    };

    /*----------------------------------------------------------------------------*/
    void
    moveToNewPos_basicMove(const kmds::GrowingView<kmds::TCellID>* ASelectionInterfaceNodes,
                           kmds::Mesh* AMesh,
                           const kmds::Variable<gmds::math::Point>* AVarNewPos,
                           const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation)
    {
        Kokkos::parallel_for(ASelectionInterfaceNodes->getNbElems(),
                             KOKKOS_LAMBDA(const int i)
                             {
                                 const kmds::TCellID nid = ASelectionInterfaceNodes->get(i);

                                 kmds::TCoord xyz[3];
                                 xyz[0] = (*AVarNewPos)[nid].X();
                                 xyz[1] = (*AVarNewPos)[nid].Y();
                                 xyz[2] = (*AVarNewPos)[nid].Z();

                                 if(std::isnan(xyz[0]) ||
                                    std::isnan(xyz[1]) ||
                                    std::isnan(xyz[2]) ||
                                    std::isinf(xyz[0]) ||
                                    std::isinf(xyz[1]) ||
                                    std::isinf(xyz[2])) {


                                     std::cerr<<"moveToNewPos_basicMove : AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAa"<<std::endl;
                                     exit(-1);

                                 } else {

                                     if ((*AVarNodeGeomAssociation)[nid] != reinterpret_cast<std::uintptr_t>(nullptr)) {
                                         gmds::math::Point pt(xyz[0], xyz[1], xyz[2]);
                                         reinterpret_cast<gmds::cad::GeomEntity *>((*AVarNodeGeomAssociation)[nid])->project(
                                                 pt);
                                         xyz[0] = pt.X();
                                         xyz[1] = pt.Y();
                                         xyz[2] = pt.Z();
                                     }

                                     AMesh->setNodeLocation(nid, xyz[0], xyz[1], xyz[2]);
                                 }
                             });

    }

    /*----------------------------------------------------------------------------*/
    void
    moveToNewPos_noBadMove_2D(const int nbIter,
                              const kmds::TCoord qualThreshold,
                              const kmds::GrowingView<kmds::TCellID>* ASelectionNodes,
                              kmds::Mesh* AMesh,
                              const kmds::Connectivity* c_N2C,
                              const kmds::Variable<gmds::math::Point>* AVarNodeDestination,
                              const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                              kmds::Variable<int>* AVarNodeNbMoves)
    {
        Kokkos::Timer timer;
        timer.reset();

        // necessary init
        const int nbNodes = AMesh->getNbNodes();
        kmds::GrowingView<kmds::TCellID> nodesIDs("NODES", nbNodes);
        AMesh->getNodeIDs(&nodesIDs);

        const int nbCells = AMesh->getNbFaces();
        kmds::GrowingView<kmds::TCellID> cellsIDs("CELLS", nbCells);
        AMesh->getFaceIDs(&cellsIDs);


        // we store the current node positions
        kmds::Variable <gmds::math::Point> *varRefPos =
                AMesh->createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF),
                                                         kmds::KMDS_NODE, "tmp_varRefPos");

        elg3d::Tools_storeNodePos_xD(AMesh, varRefPos);

        // set the starting value of the number of effective moves to 0
        AVarNodeNbMoves->setValuesTo(0);


        kmds::Variable<double>* varquality = AMesh->createVariable<double>(-HUGE_VALF, kmds::KMDS_FACE, "tmp_varcellquality");
        kmds::TCoord minScaledJacobian_before = Tools_computeScaledJacobian_2D(AMesh, varquality);

        std::cout<<"noBadMove_2D minScaledJacobian_before "<<minScaledJacobian_before<<std::endl;

        if(minScaledJacobian_before < qualThreshold) {
            throw kmds::KException("moveToNewPos_noBadMove_2D : cannot be called when current quality is below the threshold.");
        }


        std::cout<<"noBadMove_2D preptime "<<timer.seconds()<<std::endl;
        timer.reset();


        const int nbNonFixedNodes = ASelectionNodes->getNbElems();

        // build graph
        // WARNING we hard-coded 20 as the max number of neighbors
        kmds::Graph graph("NODES_GRAPH", nbNonFixedNodes, 40);

        Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> kmapNode2Vertex(nbNonFixedNodes);

        ManifoldDetection_buildGraph_N_N2F2N(ASelectionNodes, AMesh, c_N2C, &graph, &kmapNode2Vertex);

        std::cout<<"noBadMove_2D buildGraph_N_N2C2N "<<timer.seconds()<<std::endl;
        timer.reset();

        // extract independent set
        graph.buildColoring(MoveToNewPos_NBMAXCOLORS);

        std::cout<<"noBadMove_2D buildColoring "<<timer.seconds()<<std::endl;
        timer.reset();

        // convert back from graph vertex to nodes IDs

        struct MoveToNewPos_vertex2Node
        {
            kmds::GrowingView<kmds::TCellID>* selection;
            const kmds::GrowingView<kmds::TCellID>* selection_map;

            MoveToNewPos_vertex2Node(kmds::GrowingView<kmds::TCellID>* selection_,
                                     const kmds::GrowingView<kmds::TCellID>* selection_map_)
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

        int colorIndexMax = -1;

        for(int icolor=0; icolor<MoveToNewPos_NBMAXCOLORS; icolor++) {

            kmds::GrowingView<kmds::TCellID> *nodesNonFixed_indepset = graph.getColoring(icolor);

            Kokkos::parallel_for(nodesNonFixed_indepset->getNbElems(),
                                 MoveToNewPos_vertex2Node(nodesNonFixed_indepset, ASelectionNodes));

            if(nodesNonFixed_indepset->getNbElems() > 0 ) {
                colorIndexMax = icolor;
            }
            std::cout<<"noBadMove_2D colors "<<icolor<<" "<<nodesNonFixed_indepset->getNbElems()<<std::endl;
        }

        std::cout<<"noBadMove_2D vertex2nodes "<<timer.seconds()<<std::endl;
        timer.reset();


        int iter = 0;
        while(iter < nbIter) {


//            for(int icolor=0; icolor<smartLaplacian_NBMAXCOLORS; icolor++) {
            for(int icolor=0; icolor<=colorIndexMax; icolor++) {

                kmds::GrowingView<kmds::TCellID>* nodesNonFixed_indepset = graph.getColoring(icolor);


                Kokkos::parallel_for(nodesNonFixed_indepset->getNbElems(),
                                     moveToNewPos_noBadMove_oneNode_2D(nodesNonFixed_indepset,
                                                                       (kmds::TCoord) (iter+1) / (kmds::TCoord) nbIter,
                                                                       qualThreshold,
                                                                       AMesh,
                                                                       c_N2C,
                                                                       varRefPos,
                                                                       AVarNodeDestination,
                                                                       AVarNodeGeomAssociation,
                                                                       AVarNodeNbMoves
                                     ));

//                std::cout<<"noBadMove_2D colors "<<icolor<<" "<<nodesNonFixed_indepset->getNbElems()<<" "<<timer.seconds()<<std::endl;
//                timer.reset();
            }

            iter++;
        }

        std::cout<<"noBadMove_2D afteriter "<<timer.seconds()<<std::endl;

        kmds::TCoord minScaledJacobian_after = Tools_computeScaledJacobian_2D(AMesh, varquality);
        std::cout<<"noBadMove_2D minScaledJacobian "<<minScaledJacobian_before<<" "<<minScaledJacobian_after<<std::endl;


        // clean-up temporary data
        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varRefPos");
        AMesh->deleteVariable(kmds::KMDS_FACE, "tmp_varcellquality");

    }

    /*----------------------------------------------------------------------------*/
    void
    moveToNewPos_noBadMove_3D(const int nbIter,
                              const kmds::TCoord qualThreshold,
                              const kmds::GrowingView<kmds::TCellID>* ASelectionNodes,
                              kmds::Mesh* AMesh,
                              const kmds::Connectivity* c_N2C,
                              const kmds::Variable<gmds::math::Point>* AVarNodeDestination,
                              const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                              kmds::Variable<int>* AVarNodeNbMoves)
    {
        Kokkos::Timer timer;
        timer.reset();

        // necessary init
        const int nbNodes = AMesh->getNbNodes();
        kmds::GrowingView<kmds::TCellID> nodesIDs("NODES", nbNodes);
        AMesh->getNodeIDs(&nodesIDs);

        const int nbCells = AMesh->getNbRegions();
        kmds::GrowingView<kmds::TCellID> cellsIDs("CELLS", nbCells);
        AMesh->getRegionIDs(&cellsIDs);


        // we store the current node positions
        kmds::Variable <gmds::math::Point> *varRefPos =
                AMesh->createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF),
                                                       kmds::KMDS_NODE, "tmp_varRefPos");

        elg3d::Tools_storeNodePos_xD(AMesh, varRefPos);

        // set the starting value of the number of effective moves to 0
        AVarNodeNbMoves->setValuesTo(0);


        kmds::Variable<double>* varquality = AMesh->createVariable<double>(-HUGE_VALF, kmds::KMDS_REGION, "tmp_varcellquality");
        kmds::TCoord minScaledJacobian_before = Tools_computeScaledJacobian_3D(AMesh, varquality);

        std::cout<<"noBadMove_3D minScaledJacobian_before "<<minScaledJacobian_before<<std::endl;

        if(minScaledJacobian_before < qualThreshold) {
            throw kmds::KException("moveToNewPos_noBadMove_3D : cannot be called when current quality is below the threshold.");
        }


        std::cout<<"noBadMove_3D preptime "<<timer.seconds()<<std::endl;
        timer.reset();


        const int nbNonFixedNodes = ASelectionNodes->getNbElems();

        // build graph
        // WARNING we hard-coded 20 as the max number of neighbors
        kmds::Graph graph("NODES_GRAPH", nbNonFixedNodes, 40);

        Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> kmapNode2Vertex(nbNonFixedNodes);

        ManifoldDetection_buildGraph_N_N2R2N(ASelectionNodes, AMesh, c_N2C, &graph, &kmapNode2Vertex);

        std::cout<<"noBadMove_3D buildGraph_N_N2C2N "<<timer.seconds()<<std::endl;
        timer.reset();

        // extract independent set
        graph.buildColoring(MoveToNewPos_NBMAXCOLORS);

        std::cout<<"noBadMove_3D buildColoring "<<timer.seconds()<<std::endl;
        timer.reset();

        // convert back from graph vertex to nodes IDs

        struct MoveToNewPos_vertex2Node
        {
            kmds::GrowingView<kmds::TCellID>* selection;
            const kmds::GrowingView<kmds::TCellID>* selection_map;

            MoveToNewPos_vertex2Node(kmds::GrowingView<kmds::TCellID>* selection_,
                                     const kmds::GrowingView<kmds::TCellID>* selection_map_)
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

        int colorIndexMax = -1;

        for(int icolor=0; icolor<MoveToNewPos_NBMAXCOLORS; icolor++) {

            kmds::GrowingView<kmds::TCellID> *nodesNonFixed_indepset = graph.getColoring(icolor);

            Kokkos::parallel_for(nodesNonFixed_indepset->getNbElems(),
                                 MoveToNewPos_vertex2Node(nodesNonFixed_indepset, ASelectionNodes));

            if(nodesNonFixed_indepset->getNbElems() > 0 ) {
                colorIndexMax = icolor;
            }
            std::cout<<"noBadMove_3D colors "<<icolor<<" "<<nodesNonFixed_indepset->getNbElems()<<std::endl;
        }

        std::cout<<"noBadMove_3D vertex2nodes "<<timer.seconds()<<std::endl;
        timer.reset();


        int iter = 0;
        while(iter < nbIter) {


//            for(int icolor=0; icolor<smartLaplacian_NBMAXCOLORS; icolor++) {
            for(int icolor=0; icolor<=colorIndexMax; icolor++) {

                kmds::GrowingView<kmds::TCellID>* nodesNonFixed_indepset = graph.getColoring(icolor);


                Kokkos::parallel_for(nodesNonFixed_indepset->getNbElems(),
                                     moveToNewPos_noBadMove_oneNode_3D(nodesNonFixed_indepset,
                                                                       (kmds::TCoord) (iter+1) / (kmds::TCoord) nbIter,
                                                                       qualThreshold,
                                                                       AMesh,
                                                                       c_N2C,
                                                                       varRefPos,
                                                                       AVarNodeDestination,
                                                                       AVarNodeGeomAssociation,
                                                                       AVarNodeNbMoves
                                     ));

//                std::cout<<"noBadMove_3D colors "<<icolor<<" "<<nodesNonFixed_indepset->getNbElems()<<" "<<timer.seconds()<<std::endl;
//                timer.reset();
            }

            iter++;
        }

        std::cout<<"noBadMove_3D afteriter "<<timer.seconds()<<std::endl;

        kmds::TCoord minScaledJacobian_after = Tools_computeScaledJacobian_3D(AMesh, varquality);
        std::cout<<"noBadMove_3D minScaledJacobian "<<minScaledJacobian_before<<" "<<minScaledJacobian_after<<std::endl;


        // clean-up temporary data
        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varRefPos");
        AMesh->deleteVariable(kmds::KMDS_REGION, "tmp_varcellquality");

    }

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
