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
/** \file    Refinement.cpp
 *  \author  legoff
 *  \date    03/23/2020
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/Refinement.h"
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
//#include <gmds/math/Triangle.h>
//#include <gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/DATACMPNT/Parameters.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {
    /*----------------------------------------------------------------------------*/


    /*----------------------------------------------------------------------------*/
    struct refinement_checkcell_2D
    {
        const kmds::GrowingView<kmds::TCellID>* selection;
        const kmds::Mesh* mesh;
        Kokkos::View<bool *> markedNodes;


        refinement_checkcell_2D(
                const kmds::GrowingView<kmds::TCellID>* Selection_,
                const kmds::Mesh* mesh_,
                Kokkos::View<bool *> markedNodes_
        )
                : selection(Selection_)
                , mesh(mesh_)
                , markedNodes(markedNodes_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(size_t i) const {
            const kmds::TCellID cid = selection->get(i);

            const kmds::Face c = mesh->getFace(cid);

            Kokkos::View<kmds::TCellID*> nodes;
            c.nodeIds(nodes);
            const int nbNodes = nodes.size();

            std::array<bool, 4> tableBool{ {false, false, false, false} };

//            bool isSplit = false;
            for(unsigned int i_n=0; i_n<4; i_n++) {
                tableBool[i_n] = markedNodes[nodes[i_n]];

//                if(this->mesh_.isMarked(nodes[iNode],AMarkNodeToSplit)) {
//                    isSplit = true;
//                }
            }




            bool toMark = false;
            if(toMark) {
//                cells2Mark->push_back(cid);
            }
        }
    };

    /*----------------------------------------------------------------------------*/
    void
    Refinement_buildGraph_C_C2C(const kmds::GrowingView<kmds::TCellID>* ASelection,
                                const kmds::Mesh* AMesh,
                                const kmds::Connectivity* c_C2C_byN,
                                kmds::Graph* AGraph,
                                Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID>* kmap)
    {
        struct Refinement_buildCellNeighbors
        {
            const kmds::GrowingView<kmds::TCellID>* selection;

            const Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID>* kmap;
            const kmds::Mesh* mesh;
            const kmds::Connectivity* c_C2C_byN;

            kmds::Graph* graph;


            Refinement_buildCellNeighbors(
                    const kmds::GrowingView<kmds::TCellID>* Selection_,
                    const Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID>* kmap_,
                    const kmds::Mesh* AMesh_,
                    const kmds::Connectivity* c_C2C_byN_,
                    kmds::Graph* graph_
            )
                    : selection(Selection_)
                    , kmap(kmap_)
                    , mesh(AMesh_)
                    , c_C2C_byN(c_C2C_byN_)
                    , graph(graph_)
            {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(int i) const {
                kmds::TCellID cid = selection->get(i);

                Kokkos::View<kmds::TCellID*> cids;
                c_C2C_byN->get(cid, cids);

                // we keep as neighbors the cells in the selection only
                std::set<kmds::TCellID> cellsKept;

                const kmds::TSize nbCells = cids.size();

                for(kmds::TSize i_c=0; i_c<nbCells; i_c++) {
                    kmds::TCellID index = kmap->find(cids[i_c]);
                    if(kmap->valid_at(index)) {
                        cellsKept.insert(kmap->value_at(index));
                    }
                }

                graph->setNeighbors(i, cellsKept);
            }
        };

        struct Refinement_buildCellIndexMappinqg
        {
            const kmds::GrowingView<kmds::TCellID>* selection;
            Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID>* kmap;

            Refinement_buildCellIndexMappinqg(
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
                kmds::TCellID cid = selection->get(i);

                kmap->insert(cid,i);
            }
        };


        kmds::TCellID nbVec = ASelection->getNbElems();

        std::cout<<"nbVec "<<nbVec<<std::endl;

        // build mapping between nodes ids and graph vertices indices
        Kokkos::parallel_for(nbVec, Refinement_buildCellIndexMappinqg(ASelection, kmap));

        Kokkos::parallel_for(nbVec, Refinement_buildCellNeighbors(ASelection, kmap, AMesh, c_C2C_byN, AGraph));
    }

/*----------------------------------------------------------------------------*/
bool
    Refinement_validityCheck_2D(const kmds::Mesh* AMesh)
{
    // check that the mesh is a full quad-mesh
    if(AMesh->getNbQuads() == AMesh->getNbFaces()) {
        return true;
    } else {
        return false;
    }
}

/*----------------------------------------------------------------------------*/
bool
Refinement_validityCheck_3D(const kmds::Mesh* AMesh)
{
    // check that the mesh is a full hex-mesh
    if(AMesh->getNbHexahedra() == AMesh->getNbRegions()) {
        return true;
    } else {
        return false;
    }
}

/*----------------------------------------------------------------------------*/
kmds::TCellID
Refinement_getInvalidCells_2D(const Kokkos::View<bool *> markedNodes,
                              const kmds::Mesh* AMesh,
                              kmds::GrowingView<kmds::TCellID>* ASelection)
{
    const kmds::TSize nbCells = AMesh->getNbFaces();
    kmds::GrowingView <kmds::TCellID> cellIDs("CELLS", nbCells);
    AMesh->getFaceIDs(&cellIDs);

    kmds::TCellID nbCellInvalid = 0;

    Kokkos::parallel_reduce(nbCells,
                            KOKKOS_LAMBDA(const kmds::TCellID i, kmds::TCellID& sum) {
                                const kmds::TCellID cid = cellIDs.get(i);

                                const kmds::Face c = AMesh->getFace(cid);

                                Kokkos::View<kmds::TCellID*> nodes;
                                c.nodeIds(nodes);
                                const int nbNodes = nodes.size();

                                std::array<bool, 4> tableBool{ {false, false, false, false} };

                                for(unsigned int i_n=0; i_n<4; i_n++) {
                                    tableBool[i_n] = markedNodes[nodes[i_n]];
                                }

                                elg3d::tableMarkedNodes<4> tableMark(tableBool[0],tableBool[1],tableBool[2],tableBool[3]);

                                bool isInvalid = false;

                                // not relevant in this case in 2D : all marked node configurations are valid
//                             if(Refinement_tables::lookupNodesValid_3_2D.find(tableMark) == Refinement_tables::lookupNodesValid_3_2D.end()) {
//                                 isInvalid = true;
//                             }

                                if(isInvalid) {
                                    ASelection->push_back(cid);
                                    sum++;
                                }
                            },
                            nbCellInvalid);

    return nbCellInvalid;
}

/*----------------------------------------------------------------------------*/
kmds::TCellID
Refinement_getInvalidCells_3D(const Kokkos::View<bool *> markedNodes,
                              const kmds::Mesh* AMesh,
                              kmds::GrowingView<kmds::TCellID>* ASelection)
{
    const kmds::TSize nbCells = AMesh->getNbRegions();
    kmds::GrowingView <kmds::TCellID> cellIDs("CELLS", nbCells);
    AMesh->getRegionIDs(&cellIDs);

    kmds::TCellID nbCellInvalid = 0;

    Kokkos::parallel_reduce(nbCells,
                            KOKKOS_LAMBDA(const kmds::TCellID i, kmds::TCellID& sum) {
                                const kmds::TCellID cid = cellIDs.get(i);

                                const kmds::Region c = AMesh->getRegion(cid);

                                Kokkos::View<kmds::TCellID*> nodes;
                                c.nodeIds(nodes);
                                const int nbNodes = nodes.size();

                                std::array<bool, 8> tableBool{ {false, false, false, false, false, false, false, false} };

                                for(unsigned int i_n=0; i_n<8; i_n++) {
                                    tableBool[i_n] = markedNodes[nodes[i_n]];
                                }

                                elg3d::tableMarkedNodes<8> tableMark(tableBool[0],tableBool[1],tableBool[2],tableBool[3],
                                                                     tableBool[4],tableBool[5],tableBool[6],tableBool[7]);

                                bool isInvalid = false;

                                if(Refinement_tables::lookupNodesValid_3_3D.find(tableMark) == Refinement_tables::lookupNodesValid_3_3D.end()) {
                                    isInvalid = true;
                                }

                                if(isInvalid) {
                                    ASelection->push_back(cid);
                                    sum++;
                                }
                            },
                            nbCellInvalid);

    return nbCellInvalid;
}

/*----------------------------------------------------------------------------*/
    void
    Refinement_solveInvalidCells_2D(const kmds::GrowingView<kmds::TCellID>* ASelection,
                                    const kmds::Mesh* AMesh,
                                    Kokkos::View<bool *> markedNodes)
    {
        const kmds::TSize nbCells = ASelection->getNbElems();

        Kokkos::parallel_for(nbCells,
                             KOKKOS_LAMBDA(const kmds::TCellID i) {
                                 const kmds::TCellID cid = ASelection->get(i);

                                 const kmds::Face c = AMesh->getFace(cid);

                                 Kokkos::View<kmds::TCellID*> nodes;
                                 c.nodeIds(nodes);
                                 const int nbNodes = nodes.size();

                                 std::array<bool, 4> tableBool{ {false, false, false, false} };

                                 for(unsigned int i_n=0; i_n<4; i_n++) {
                                     tableBool[i_n] = markedNodes[nodes[i_n]];
                                 }

                                 elg3d::tableMarkedNodes<4> tableMark(tableBool[0],tableBool[1],tableBool[2],tableBool[3]);

                                 // not relevant in this case in 2D : all marked node configurations are valid
                                 tableMarkedNodes<4> tableMarkTarget;
                                 if(Refinement_tables::lookupNodesSolve_3_2D.find(tableMark) != Refinement_tables::lookupNodesSolve_3_2D.end()) {
                                     tableMarkTarget = Refinement_tables::lookupNodesSolve_3_2D[tableMark];
                                 } else {
                                     tableMarkTarget = tableMarkedNodes<4> {1,1,1,1};
                                 }

                                 //
                                 if(tableMarkTarget != tableMark) {

                                     for(unsigned int i_n=0; i_n<4; i_n++) {

                                         markedNodes[nodes[i_n]] = tableMarkTarget.get(i_n);
                                     }

                                 }
                             });

    }

    /*----------------------------------------------------------------------------*/
    void
    Refinement_solveInvalidCells_3D(const kmds::GrowingView<kmds::TCellID>* ASelection,
                                    const kmds::Mesh* AMesh,
                                    Kokkos::View<bool *> markedNodes)
    {
        const kmds::TSize nbCells = ASelection->getNbElems();

        Kokkos::parallel_for(nbCells,
                             KOKKOS_LAMBDA(const kmds::TCellID i) {
                                 const kmds::TCellID cid = ASelection->get(i);

                                 const kmds::Region c = AMesh->getRegion(cid);

                                 Kokkos::View<kmds::TCellID*> nodes;
                                 c.nodeIds(nodes);
                                 const int nbNodes = nodes.size();

                                 std::array<bool, 8> tableBool{ {false, false, false, false, false, false, false, false} };

                                 for(unsigned int i_n=0; i_n<8; i_n++) {
                                     tableBool[i_n] = markedNodes[nodes[i_n]];
                                 }

                                 elg3d::tableMarkedNodes<8> tableMark(tableBool[0],tableBool[1],tableBool[2],tableBool[3],
                                                                      tableBool[4],tableBool[5],tableBool[6],tableBool[7]);

                                 // not relevant in this case in 2D : all marked node configurations are valid
                                 tableMarkedNodes<8> tableMarkTarget;
                                 if(Refinement_tables::lookupNodesSolve_3_3D.find(tableMark) != Refinement_tables::lookupNodesSolve_3_3D.end()) {
                                     tableMarkTarget = Refinement_tables::lookupNodesSolve_3_3D[tableMark];
                                 } else {
                                     tableMarkTarget = tableMarkedNodes<8> {1,1,1,1,1,1,1,1};
                                 }

                                 //
                                 if(tableMarkTarget != tableMark) {

                                     for(unsigned int i_n=0; i_n<8; i_n++) {

                                         markedNodes[nodes[i_n]] = tableMarkTarget.get(i_n);
                                     }
                                 }
                             });

    }

/*----------------------------------------------------------------------------*/
    kmds::TCellID Refinement_avoidBadConfig_2D(Kokkos::View<bool *> markedCells,
                                  Kokkos::View<bool *> markedNodes,
                                  const kmds::Mesh* AMesh,
                                  const kmds::Connectivity* c_C2C_byN)
{
    const kmds::TSize nbCells = AMesh->getNbFaces();
    kmds::GrowingView <kmds::TCellID> cellIDs("CELLS", nbCells);
    AMesh->getFaceIDs(&cellIDs);

    int nbCellsInvalid = nbCells;

    do{

        // get cells that we will investigate
        kmds::GrowingView <kmds::TCellID> cellsInvalid("cellSuspicious", nbCells);

        Refinement_getInvalidCells_2D(markedNodes,
                                      AMesh,
                                      &cellsInvalid);

        nbCellsInvalid = cellsInvalid.getNbElems();


        // build graph
        // WARNING we hard-coded 20 as the max number of neighbors
        kmds::Graph graph("NODES_NONMANIFOLD_GRAPH", nbCellsInvalid, 20);
        Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> kmapCell2Vertex(nbCellsInvalid);
        Refinement_buildGraph_C_C2C(&cellsInvalid, AMesh, c_C2C_byN, &graph, &kmapCell2Vertex);

        // extract independent set
        kmds::GrowingView<kmds::TCellID> cellsInvalid_indepset("cellsInvalid_indepset", nbCellsInvalid);
        graph.getIndependentSet(&cellsInvalid_indepset);

        // convert back from graph vertex to cell IDs
        struct AssignCells_vertex2Cell
        {
            kmds::GrowingView<kmds::TCellID>* selection;
            kmds::GrowingView<kmds::TCellID>* selection_map;

            AssignCells_vertex2Cell(kmds::GrowingView<kmds::TCellID>* selection_, kmds::GrowingView<kmds::TCellID>* selection_map_)
                    : selection(selection_)
                    , selection_map(selection_map_)
            {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(int i) const {
                const kmds::TCellID cid = selection->get(i);

                selection->set(i, selection_map->get(cid));
            }
        };

        Kokkos::parallel_for(cellsInvalid_indepset.getNbElems(), AssignCells_vertex2Cell(&cellsInvalid_indepset, &cellsInvalid));

        // solve this indset by marking the necessary nodes avoid bad configurations
        Refinement_solveInvalidCells_2D(&cellsInvalid_indepset,
                                        AMesh,
                                        markedNodes);

    } while (nbCellsInvalid > 0);

}

/*----------------------------------------------------------------------------*/
    kmds::TCellID Refinement_avoidBadConfig_3D(Kokkos::View<bool *> markedCells,
                                  Kokkos::View<bool *> markedNodes,
                                  const kmds::Mesh* AMesh,
                                  const kmds::Connectivity* c_C2C_byN)
{
    const kmds::TSize nbCells = AMesh->getNbRegions();
    kmds::GrowingView <kmds::TCellID> cellIDs("CELLS", nbCells);
    AMesh->getRegionIDs(&cellIDs);

    int nbCellsInvalid = nbCells;

    do{

        // get cells that we will investigate
        kmds::GrowingView <kmds::TCellID> cellsInvalid("cellSuspicious", nbCells);

        Refinement_getInvalidCells_3D(markedNodes,
                                      AMesh,
                                      &cellsInvalid);

        nbCellsInvalid = cellsInvalid.getNbElems();

        std::cout<<"nbCellsInvalid "<<nbCellsInvalid<<std::endl;


        // build graph
        // WARNING we hard-coded 40 as the max number of neighbors
        kmds::Graph graph("NODES_NONMANIFOLD_GRAPH", nbCellsInvalid, 40);
        Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> kmapCell2Vertex(nbCellsInvalid);
        Refinement_buildGraph_C_C2C(&cellsInvalid, AMesh, c_C2C_byN, &graph, &kmapCell2Vertex);

        // extract independent set
        kmds::GrowingView<kmds::TCellID> cellsInvalid_indepset("cellsInvalid_indepset", nbCellsInvalid);
        graph.getIndependentSet(&cellsInvalid_indepset);

        // convert back from graph vertex to cell IDs
        struct AssignCells_vertex2Cell
        {
            kmds::GrowingView<kmds::TCellID>* selection;
            kmds::GrowingView<kmds::TCellID>* selection_map;

            AssignCells_vertex2Cell(kmds::GrowingView<kmds::TCellID>* selection_, kmds::GrowingView<kmds::TCellID>* selection_map_)
                    : selection(selection_)
                    , selection_map(selection_map_)
            {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(int i) const {
                const kmds::TCellID cid = selection->get(i);

                selection->set(i, selection_map->get(cid));
            }
        };

        Kokkos::parallel_for(cellsInvalid_indepset.getNbElems(), AssignCells_vertex2Cell(&cellsInvalid_indepset, &cellsInvalid));

        // solve this indset by marking the necessary nodes avoid bad configurations
        Refinement_solveInvalidCells_3D(&cellsInvalid_indepset,
                                        AMesh,
                                        markedNodes);

    } while (nbCellsInvalid > 0);

}

/*----------------------------------------------------------------------------*/
    kmds::TCellID
    Refinement_createNodesEdges_2D(const Kokkos::View<bool *> markedEdges,
                                   const Kokkos::View<bool *> markedNodes,
                                   kmds::Mesh* AMesh,
                                   const kmds::Connectivity* c_E2C,
                                   Kokkos::View<kmds::TSize *> ACell2NodesExt,
                                   Kokkos::View<kmds::TCellID *[8]> ANodesExt)
    {
        kmds::TSize nbEdges = AMesh->getNbEdges();
        kmds::GrowingView <kmds::TCellID> edgeIDs("EDGES", nbEdges);
        AMesh->getEdgeIDs(&edgeIDs);


        // first update the capacity of the mesh for nodes
        kmds::TSize nbNewNodes = 0;
        Kokkos::parallel_reduce(nbEdges,
                                KOKKOS_LAMBDA(const kmds::TCellID i, kmds::TSize& sum) {
                                    const kmds::TCellID eid = edgeIDs.get(i);

                                    const kmds::Edge e = AMesh->getEdge(eid);

                                    Kokkos::View<kmds::TCellID *> nodes;
                                    e.nodeIds(nodes);

                                    if (markedNodes[nodes[0]] && markedNodes[nodes[1]]) {
                                        sum += 2;
                                    } else {
                                        if (markedNodes[nodes[0]] || markedNodes[nodes[1]]) {
                                            sum++;
                                        }
                                    }
                                },
                                nbNewNodes);

        AMesh->updateNodeCapacity(AMesh->getNbNodes() + nbNewNodes);


        Kokkos::parallel_for(nbEdges,
                             KOKKOS_LAMBDA(const kmds::TCellID i) {
                                 const kmds::TCellID eid = edgeIDs.get(i);

                                 const kmds::Edge e = AMesh->getEdge(eid);

                                 Kokkos::View<kmds::TCellID *> nodes;
                                 e.nodeIds(nodes);

                                 if(markedNodes[nodes[0]] || markedNodes[nodes[1]]) {

                                     gmds::math::Point pt0 = AMesh->getNodeLocation(nodes[0]);
                                     gmds::math::Point pt1 = AMesh->getNodeLocation(nodes[1]);

                                     kmds::TCellID firstID  = kmds::NullID;
                                     kmds::TCellID secondID = kmds::NullID;

                                     if(markedNodes(nodes[0]) && markedNodes(nodes[1])) {
                                         firstID = AMesh->addNodes(2);
                                         secondID = firstID + 1;
                                         AMesh->setNodeLocation(firstID, pt0 + (1./3.)*(pt1-pt0));
                                         AMesh->setNodeLocation(secondID, pt0 + (2./3.)*(pt1-pt0));
                                     } else {
                                         if(markedNodes[nodes[0]]) {
                                             firstID = AMesh->newNode(pt0 + (1./3.)*(pt1-pt0));
                                         } else {
                                             secondID = AMesh->newNode(pt0 + (2./3.)*(pt1-pt0));
                                         }

                                     }

                                     // put the newly created nodes in their cells storage
                                     Kokkos::View<kmds::TCellID*> cids;
                                     c_E2C->get(eid, cids);

                                     const kmds::TSize nbCells = cids.size();

                                     for(kmds::TSize i_c=0; i_c<nbCells; i_c++) {

                                         const kmds::Face c = AMesh->getFace(cids[i_c]);
                                         bool isIJK;
                                         int numEdge = -1;
                                         c.getEdgeInfo(nodes, numEdge, isIJK);

                                         if(isIJK) {
                                             ANodesExt(ACell2NodesExt[cids[i_c]],numEdge*2)   = firstID;
                                             ANodesExt(ACell2NodesExt[cids[i_c]],numEdge*2+1) = secondID;
                                         } else {
                                             ANodesExt(ACell2NodesExt[cids[i_c]],numEdge*2)   = secondID;
                                             ANodesExt(ACell2NodesExt[cids[i_c]],numEdge*2+1) = firstID;
                                         }

                                     }
                                 }
                             });

    }

/*----------------------------------------------------------------------------*/
kmds::TCellID
Refinement_createNodesEdges_3D(const Kokkos::View<bool *> markedEdges,
                               const Kokkos::View<bool *> markedNodes,
                               kmds::Mesh* AMesh,
                               const kmds::Connectivity* c_E2C,
                               Kokkos::View<kmds::TSize *> ACell2NodesExt,
                               Kokkos::View<kmds::TCellID *[48]> ANodesExt)
{
    kmds::TSize nbEdges = AMesh->getNbEdges();
    kmds::GrowingView <kmds::TCellID> edgeIDs("EDGES", nbEdges);
    AMesh->getEdgeIDs(&edgeIDs);

    // first update the capacity of the mesh for nodes
    kmds::TSize nbNewNodes = 0;
    Kokkos::parallel_reduce(nbEdges,
                            KOKKOS_LAMBDA(const kmds::TCellID i, kmds::TSize& sum) {
                                const kmds::TCellID eid = edgeIDs.get(i);

                                const kmds::Edge e = AMesh->getEdge(eid);

                                Kokkos::View<kmds::TCellID *> nodes;
                                e.nodeIds(nodes);

                                if (markedNodes[nodes[0]] && markedNodes[nodes[1]]) {
                                    sum += 2;
                                } else {
                                    if (markedNodes[nodes[0]] || markedNodes[nodes[1]]) {
                                        sum++;
                                    }
                                }
                            },
                            nbNewNodes);

    AMesh->updateNodeCapacity(AMesh->getNbNodes() + nbNewNodes);


    Kokkos::parallel_for(nbEdges,
                         KOKKOS_LAMBDA(const kmds::TCellID i) {
                             const kmds::TCellID eid = edgeIDs.get(i);

                             const kmds::Edge e = AMesh->getEdge(eid);

                             Kokkos::View<kmds::TCellID *> nodes;
                             e.nodeIds(nodes);

                             if(markedNodes[nodes[0]] || markedNodes[nodes[1]]) {

                                 gmds::math::Point pt0 = AMesh->getNodeLocation(nodes[0]);
                                 gmds::math::Point pt1 = AMesh->getNodeLocation(nodes[1]);

                                 kmds::TCellID firstID  = kmds::NullID;
                                 kmds::TCellID secondID = kmds::NullID;

                                 if(markedNodes(nodes[0]) && markedNodes(nodes[1])) {
                                     firstID = AMesh->addNodes(2);
                                     secondID = firstID + 1;
                                     AMesh->setNodeLocation(firstID, pt0 + (1./3.)*(pt1-pt0));
                                     AMesh->setNodeLocation(secondID, pt0 + (2./3.)*(pt1-pt0));
                                 } else {
                                     if(markedNodes[nodes[0]]) {
                                         firstID = AMesh->newNode(pt0 + (1./3.)*(pt1-pt0));
                                     } else {
                                         secondID = AMesh->newNode(pt0 + (2./3.)*(pt1-pt0));
                                     }

                                 }

                                 // put the newly created nodes in their cells storage
                                 Kokkos::View<kmds::TCellID*> cids;
                                 c_E2C->get(eid, cids);

                                 const kmds::TSize nbCells = cids.size();

                                 for(kmds::TSize i_c=0; i_c<nbCells; i_c++) {

                                     const kmds::Region c = AMesh->getRegion(cids[i_c]);
                                     bool isIJK;
                                     int numEdge = -1;
                                     c.getEdgeInfo(nodes, numEdge, isIJK);

                                     if(isIJK) {
                                         ANodesExt(ACell2NodesExt[cids[i_c]],numEdge*2)   = firstID;
                                         ANodesExt(ACell2NodesExt[cids[i_c]],numEdge*2+1) = secondID;
                                     } else {
                                         ANodesExt(ACell2NodesExt[cids[i_c]],numEdge*2)   = secondID;
                                         ANodesExt(ACell2NodesExt[cids[i_c]],numEdge*2+1) = firstID;
                                     }

                                 }
                             }
                         });

}

/*----------------------------------------------------------------------------*/
kmds::TCellID
Refinement_createNodesFaces_3D(const Kokkos::View<bool *> markedFaces,
                               const Kokkos::View<bool *> markedNodes,
                               kmds::Mesh* AMesh,
                               const kmds::Connectivity* c_F2C,
                               Kokkos::View<kmds::TSize *> ACell2NodesExt,
                               Kokkos::View<kmds::TCellID *[48]> ANodesExt)
{
    kmds::TSize nbFaces = AMesh->getNbFaces();
    kmds::GrowingView<kmds::TCellID> faceIDs("FACES", nbFaces);
    AMesh->getFaceIDs(&faceIDs);

    // first update the capacity of the mesh for nodes
    kmds::TSize nbNewNodes = 0;
    Kokkos::parallel_reduce(nbFaces,
                            KOKKOS_LAMBDA(const kmds::TCellID i, kmds::TSize& sum) {
                                const kmds::TCellID fid = faceIDs.get(i);

                                const kmds::Face f = AMesh->getFace(fid);

                                Kokkos::View<kmds::TCellID *> nodes;
                                f.nodeIds(nodes);

                                elg3d::tableMarkedNodes<4> tableMark(markedNodes[nodes[0]], markedNodes[nodes[1]],
                                                                     markedNodes[nodes[2]], markedNodes[nodes[3]]);

                                sum += Refinement_tables::nbNodesOnFaceToBuild_3b_2D[tableMark];
                            },
                            nbNewNodes);

    AMesh->updateNodeCapacity(AMesh->getNbNodes() + nbNewNodes);


    Kokkos::parallel_for(nbFaces,
                         KOKKOS_LAMBDA(const kmds::TCellID i) {
                             const kmds::TCellID fid = faceIDs.get(i);

                             const kmds::Face f = AMesh->getFace(fid);

                             Kokkos::View<kmds::TCellID *> nodes;
                             f.nodeIds(nodes);

                             elg3d::tableMarkedNodes<4> tableMark(markedNodes[nodes[0]], markedNodes[nodes[1]],
                                                                  markedNodes[nodes[2]], markedNodes[nodes[3]]);

                             if(tableMark != elg3d::tableMarkedNodes<4> {0,0,0,0}) {
                                 gmds::math::Point pt0 = AMesh->getNodeLocation(nodes[0]);
                                 gmds::math::Point pt1 = AMesh->getNodeLocation(nodes[1]);
                                 gmds::math::Point pt2 = AMesh->getNodeLocation(nodes[2]);
                                 gmds::math::Point pt3 = AMesh->getNodeLocation(nodes[3]);


                                 const int nbNewNodes = Refinement_tables::nbNodesOnFaceToBuild_3b_2D[tableMark];
                                 kmds::TCellID newNodeIDs[4];
                                 newNodeIDs[0] = AMesh->addNodes(nbNewNodes);
                                 newNodeIDs[1] = kmds::NullID;
                                 newNodeIDs[2] = kmds::NullID;
                                 newNodeIDs[3] = kmds::NullID;

                                 for(int i_n=0; i_n<nbNewNodes; i_n++) {
                                     const double ix = Refinement_tables::nodesOnFaceToBuild_3b_2D[tableMark][i_n][0];
                                     const double iy = Refinement_tables::nodesOnFaceToBuild_3b_2D[tableMark][i_n][1];
                                     gmds::math::Point pt =
                                             (6 - ix) / 6. * (6 - iy) / 6. * AMesh->getNodeLocation(nodes[0])
                                             + ix / 6. * (6 - iy) / 6. * AMesh->getNodeLocation(nodes[1])
                                             + ix / 6. * iy / 6. * AMesh->getNodeLocation(nodes[2])
                                             + (6 - ix) / 6. * iy / 6. * AMesh->getNodeLocation(nodes[3]);

                                     newNodeIDs[i_n] = newNodeIDs[0] + i_n;
                                     AMesh->setNodeLocation(newNodeIDs[0] + i_n, pt);
                                 }

//                                 if(markedNodes(nodes[0]) && markedNodes(nodes[1])) {
//                                     firstID = AMesh->addNodes(2);
//                                     secondID = firstID + 1;
//                                     AMesh->setNodeLocation(firstID, pt0 + (1./3.)*(pt1-pt0));
//                                     AMesh->setNodeLocation(secondID, pt0 + (2./3.)*(pt1-pt0));
//                                 } else {
//                                     if(markedNodes[nodes[0]]) {
//                                         firstID = AMesh->newNode(pt0 + (1./3.)*(pt1-pt0));
//                                     } else {
//                                         secondID = AMesh->newNode(pt0 + (2./3.)*(pt1-pt0));
//                                     }
//
//                                 }

                                 // put the newly created nodes in their cells storage
                                 Kokkos::View<kmds::TCellID*> cids;
                                 c_F2C->get(fid, cids);

                                 const kmds::TSize nbCells = cids.size();

                                 for(kmds::TSize i_c=0; i_c<nbCells; i_c++) {

                                     const kmds::Region c = AMesh->getRegion(cids[i_c]);
                                     bool isIJK;
                                     int numFace = -1;
                                     int offset = 0;
                                     c.getFaceInfo(nodes[0],
                                                   nodes[1],
                                                   nodes[2],
                                                   nodes[3],
                                                   numFace,
                                                   offset,
                                                   isIJK);

                                     if(isIJK) {
                                         switch(offset) {
                                             case 0:
                                                 ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 0) = newNodeIDs[0];
                                                 ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 1) = newNodeIDs[1];
                                                 ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 2) = newNodeIDs[2];
                                                 ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 3) = newNodeIDs[3];
                                                 break;
                                             case 1:
                                                 ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 0) = newNodeIDs[2];
                                                 ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 1) = newNodeIDs[0];
                                                 ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 2) = newNodeIDs[3];
                                                 ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 3) = newNodeIDs[1];
                                                 break;
                                             case 2:
                                                 ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 0) = newNodeIDs[3];
                                                 ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 1) = newNodeIDs[2];
                                                 ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 2) = newNodeIDs[1];
                                                 ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 3) = newNodeIDs[0];
                                                 break;
                                             case 3:
                                                 ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 0) = newNodeIDs[1];
                                                 ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 1) = newNodeIDs[3];
                                                 ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 2) = newNodeIDs[0];
                                                 ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 3) = newNodeIDs[2];
                                                 break;
                                             default :
                                                 exit(-1);
                                                 break;
                                         }
                                     } else {
                                        switch(offset) {
                                            case 0:
                                                ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 0) = newNodeIDs[0];
                                                ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 1) = newNodeIDs[2];
                                                ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 2) = newNodeIDs[1];
                                                ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 3) = newNodeIDs[3];
                                                break;
                                            case 1:
                                                ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 0) = newNodeIDs[1];
                                                ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 1) = newNodeIDs[0];
                                                ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 2) = newNodeIDs[3];
                                                ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 3) = newNodeIDs[2];
                                                break;
                                            case 2:
                                                ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 0) = newNodeIDs[3];
                                                ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 1) = newNodeIDs[1];
                                                ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 2) = newNodeIDs[2];
                                                ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 3) = newNodeIDs[0];
                                                break;
                                            case 3:
                                                ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 0) = newNodeIDs[2];
                                                ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 1) = newNodeIDs[3];
                                                ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 2) = newNodeIDs[0];
                                                ANodesExt(ACell2NodesExt[cids[i_c]], 24 + numFace * 4 + 3) = newNodeIDs[1];
                                                break;
                                            default :
                                                exit(-1);
                                                break;
                                        }
                                     }
                                 }

                             }

                         });
}

/*----------------------------------------------------------------------------*/
void
Refinement_refine_2D(const kmds::GrowingView<kmds::TCellID>* ACells2Refine,
                     const bool AIsCellList,
                     kmds::Mesh* AMesh,
                     const kmds::Connectivity* c_C2C_byN,
                     const kmds::Connectivity* c_E2C)
{
    // check the validity domain
    bool isValid = Refinement_validityCheck_2D(AMesh);
    if(!isValid) {
        std::cerr<<"Refinement_refine_2D invalid mesh."<<std::endl;
        exit(-1);
    }

    const kmds::TSize nbCells = AMesh->getNbFaces();
    kmds::GrowingView <kmds::TCellID> cellIDs("CELLS", nbCells);
    AMesh->getFaceIDs(&cellIDs);

    // mark for the cells that must be refined
    Kokkos::View<bool *> markedCells("markedCells", AMesh->getFaceSupID()+1);

    // mark for the edges that must be refined
    Kokkos::View<bool *> markedEdges("markedEdges", AMesh->getEdgeSupID()+1);

    // mark for the nodes that must be refined
    Kokkos::View<bool *, Kokkos::MemoryTraits<Kokkos::Atomic> > markedNodes("markedNodes", AMesh->getNodeSupID()+1);


    // initial marking for the cells and nodes that we explicitly want the refinement of
    if(AIsCellList) {
        Kokkos::parallel_for(ACells2Refine->getNbElems(),
                             KOKKOS_LAMBDA(const kmds::TCellID i) {
                                 const kmds::TCellID cid = ACells2Refine->get(i);

                                 markedCells(cid) = true;

                                 // mark their nodes
                                 const kmds::Face c = AMesh->getFace(cid);

                                 Kokkos::View<kmds::TCellID *> nodes;
                                 c.nodeIds(nodes);
                                 const int nbNodes = nodes.size();

                                 for (auto i_n = 0; i_n < nbNodes; i_n++) {
                                     markedNodes(nodes[i_n]) = true;
                                 }
                             });
    } else {
        Kokkos::parallel_for(ACells2Refine->getNbElems(),
                             KOKKOS_LAMBDA(const kmds::TCellID i) {
                                 const kmds::TCellID nid = ACells2Refine->get(i);
                                 markedNodes(nid) = true;
                             });

        Kokkos::parallel_for(nbCells,
                             KOKKOS_LAMBDA(const kmds::TCellID i) {
                                 const kmds::TCellID cid = cellIDs.get(i);

                                 const kmds::Face c = AMesh->getFace(cid);

                                 Kokkos::View<kmds::TCellID*> nodes;
                                 c.nodeIds(nodes);
                                 const int nbNodes = nodes.size();

                                 for(auto i_n=0; i_n<nbNodes; i_n++) {
                                     if(markedNodes(nodes[i_n])) {
                                         markedCells(cid) = true;
                                         break;
                                     }
                                 }
                             });
    }

    Refinement_avoidBadConfig_2D(markedCells, markedNodes, AMesh, c_C2C_byN);

    // mark the cells that have marked nodes
    Kokkos::parallel_for(nbCells,
                         KOKKOS_LAMBDA(const kmds::TCellID i) {
                             const kmds::TCellID cid = cellIDs.get(i);

                             const kmds::Face c = AMesh->getFace(cid);

                             Kokkos::View<kmds::TCellID*> nodes;
                             c.nodeIds(nodes);
                             const int nbNodes = nodes.size();

                             for(auto i_n=0; i_n<nbNodes; i_n++) {
                                 if(markedNodes(nodes[i_n])) {
                                     markedCells(cid) = true;
                                     break;
                                 }
                             }

                             // permut the nodes
                         });

// update the capacity of the mesh


//    kmds::GrowingView <kmds::TCellID> cells2Mark("cells2Refine", AMesh->getNbFaces());
//    kmds::TSize nbCells2Mark = ACells2Refine->getNbElems();
//    for(auto i=0; i<nbCells2Mark; i++) {
//        cells2Mark.push_back(ACells2Refine->get(i));
//    }

        // storage of the new nodes created by edges and faces for the cells that must be refined
        const kmds::TCellID maxID_c = AMesh->getFaceSupID();
        Kokkos::View<kmds::TSize *> cell2nodesExt("markedCells", maxID_c+1);

        kmds::TCellID nbMarkedCells = 0;
        for(kmds::TSize i_c=0; i_c<=maxID_c; i_c++) {
            if(markedCells[i_c]) {
                cell2nodesExt(i_c) = nbMarkedCells;
                nbMarkedCells++;
            }
        }

        // storage initialized to NullID
        Kokkos::View<kmds::TCellID *[8]> nodesExt("markedCells", nbMarkedCells);
        Kokkos::parallel_for(nbMarkedCells,
                             KOKKOS_LAMBDA(const kmds::TCellID i) {
                                 nodesExt(i,0) = kmds::NullID;
                                 nodesExt(i,1) = kmds::NullID;
                                 nodesExt(i,2) = kmds::NullID;
                                 nodesExt(i,3) = kmds::NullID;
                                 nodesExt(i,4) = kmds::NullID;
                                 nodesExt(i,5) = kmds::NullID;
                                 nodesExt(i,6) = kmds::NullID;
                                 nodesExt(i,7) = kmds::NullID;
                             });

        // create nodes for edges
        Refinement_createNodesEdges_2D(markedEdges,
                                       markedNodes,
                                       AMesh,
                                       c_E2C,
                                       cell2nodesExt,
                                       nodesExt);


        // update the mesh capacity for nodes and cells
        kmds::TSize nbNewNodes = 0;
        kmds::TSize nbNewCells = 0;
        Kokkos::parallel_reduce(maxID_c+1,
                             KOKKOS_LAMBDA(const kmds::TCellID i, kmds::TSize& sum) {

                                 const kmds::TCellID cid = i;

                                 if (markedCells[cid]) {
                                     const kmds::Face c = AMesh->getFace(cid);

                                     Kokkos::View<kmds::TCellID *> nodes;
                                     c.nodeIds(nodes);
                                     const int nbNodes = nodes.size();

                                     elg3d::tableMarkedNodes<4> tableMark(markedNodes[nodes[0]], markedNodes[nodes[1]],
                                                                          markedNodes[nodes[2]], markedNodes[nodes[3]]);

                                     sum += Refinement_tables::nbNodesOnFaceToBuild_3b_2D[tableMark];
                                 }
                             },
                             nbNewNodes);

        Kokkos::parallel_reduce(maxID_c+1,
                                KOKKOS_LAMBDA(const kmds::TCellID i, kmds::TSize& sum) {

                                    const kmds::TCellID cid = i;

                                    if (markedCells[cid]) {
                                        const kmds::Face c = AMesh->getFace(cid);

                                        Kokkos::View<kmds::TCellID *> nodes;
                                        c.nodeIds(nodes);
                                        const int nbNodes = nodes.size();

                                        elg3d::tableMarkedNodes<4> tableMark(markedNodes[nodes[0]], markedNodes[nodes[1]],
                                                                             markedNodes[nodes[2]], markedNodes[nodes[3]]);

                                        sum += Refinement_tables::nbFacesOnFaceToBuild_3b_2D[tableMark];
                                    }
                                },
                                nbNewCells);

        AMesh->updateNodeCapacity(AMesh->getNbNodes() + nbNewNodes);
        AMesh->updateFaceCapacity(AMesh->getNbFaces() + nbNewCells);


        // nodes for cells and the new cells
        Kokkos::parallel_for(maxID_c+1,
                             KOKKOS_LAMBDA(const kmds::TCellID i) {

                                const kmds::TCellID cid = i;

                                if(markedCells[cid]) {
                                    const kmds::Face c = AMesh->getFace(cid);

                                    Kokkos::View<kmds::TCellID*> nodes;
                                    c.nodeIds(nodes);
                                    const int nbNodes = nodes.size();

                                    elg3d::tableMarkedNodes<4> tableMark(markedNodes[nodes[0]], markedNodes[nodes[1]],
                                                                         markedNodes[nodes[2]], markedNodes[nodes[3]]);

                                    // regroup the nodes :
                                    kmds::TCellID subNodeIDs[7*7];

                                    // nodes of the region
                                    subNodeIDs[0*7 + 0] = nodes[0];
                                    subNodeIDs[6*7 + 0] = nodes[1];
                                    subNodeIDs[6*7 + 6] = nodes[2];
                                    subNodeIDs[0*7 + 6] = nodes[3];

                                    // nodes created by edges
                                    for(int i_n=0; i_n<8; i_n++) {
                                        const int index = Refinement_tables::nodesOnFaceFromExternal_3ab_2D[i_n];
                                        subNodeIDs[index] = nodesExt(cell2nodesExt[cid], i_n);
                                    }

                                    // new nodes built exclusively for this region
                                    const int nbNewNodes = Refinement_tables::nbNodesOnFaceToBuild_3b_2D[tableMark];

                                    kmds::TCellID firstIDN = AMesh->addNodes(nbNewNodes);

                                    for(int i_n=0; i_n<nbNewNodes; i_n++) {
                                        const double ix = Refinement_tables::nodesOnFaceToBuild_3b_2D[tableMark][i_n][0];
                                        const double iy = Refinement_tables::nodesOnFaceToBuild_3b_2D[tableMark][i_n][1];
                                        gmds::math::Point pt = (6 - ix)/6. * (6 - iy)/6. * AMesh->getNodeLocation(nodes[0])
                                                               + ix/6. * (6 - iy)/6. * AMesh->getNodeLocation(nodes[1])
                                                               + ix/6. * iy/6. * AMesh->getNodeLocation(nodes[2])
                                                               + (6 - ix)/6. * iy/6. * AMesh->getNodeLocation(nodes[3]);


                                        const int index = ix*7 + iy;
                                        subNodeIDs[index] = firstIDN + i_n;

                                        AMesh->setNodeLocation(firstIDN + i_n, pt);
                                    }


                                    const int nbNewCells = Refinement_tables::nbFacesOnFaceToBuild_3b_2D[tableMark];
                                    std::vector<std::array<int, 4> > cells2Build = Refinement_tables::facesOnFaceToBuild_3b_2D[tableMark];

                                    // the first one will take the place of the current region
                                    nodes[0] = subNodeIDs[cells2Build[0][0]];
                                    nodes[1] = subNodeIDs[cells2Build[0][1]];
                                    nodes[2] = subNodeIDs[cells2Build[0][2]];
                                    nodes[3] = subNodeIDs[cells2Build[0][3]];

//                             c.setNodes(nodes);

                                    kmds::TCellID firstIDC = AMesh->addQuads(nbNewCells-1);

                                    for(int i_c=1; i_c<nbNewCells; i_c++) {
                                        kmds::Face c_bis = AMesh->getFace(firstIDC + i_c - 1);

                                        Kokkos::View<kmds::TCellID*> nodes_bis;
                                        c_bis.nodeIds(nodes_bis);

                                        nodes_bis[0] = subNodeIDs[cells2Build[i_c][0]];
                                        nodes_bis[1] = subNodeIDs[cells2Build[i_c][1]];
                                        nodes_bis[2] = subNodeIDs[cells2Build[i_c][2]];
                                        nodes_bis[3] = subNodeIDs[cells2Build[i_c][3]];

//                                 c_bis.setNodes(nodes_bis);
                                    }

                                }  // if(markedCells[cid])

                             });



}

/*----------------------------------------------------------------------------*/
void
Refinement_refine_3D(const kmds::GrowingView<kmds::TCellID>* ACells2Refine,
                     const bool AIsCellList,
                     kmds::Mesh* AMesh,
                     const kmds::Connectivity* c_C2C_byN,
                     const kmds::Connectivity* c_E2C,
                     const kmds::Connectivity* c_F2C)
{
       // check the validity domain
       bool isValid = Refinement_validityCheck_3D(AMesh);
       if (!isValid) {
           std::cerr << "Refinement_refine_3D invalid mesh." << std::endl;
           exit(-1);
       }

    const kmds::TSize nbCells = AMesh->getNbRegions();
    kmds::GrowingView <kmds::TCellID> cellIDs("CELLS", nbCells);
    AMesh->getRegionIDs(&cellIDs);

    // mark for the cells that must be refined
    Kokkos::View<bool *> markedCells("markedCells", AMesh->getRegionSupID()+1);

    // mark for the faces that must be refined
    Kokkos::View<bool *> markedFaces("markedFaces", AMesh->getFaceSupID()+1);

    // mark for the edges that must be refined
    Kokkos::View<bool *> markedEdges("markedEdges", AMesh->getEdgeSupID()+1);

    // mark for the nodes that must be refined
    Kokkos::View<bool *, Kokkos::MemoryTraits<Kokkos::Atomic> > markedNodes("markedNodes", AMesh->getNodeSupID()+1);


    // initial marking for the cells and nodes that we explicitly want the refinement of
    if(AIsCellList) {
        Kokkos::parallel_for(ACells2Refine->getNbElems(),
                             KOKKOS_LAMBDA(const kmds::TCellID i) {
                                 const kmds::TCellID cid = ACells2Refine->get(i);

                                 markedCells(cid) = true;

                                 // mark their nodes
                                 const kmds::Region c = AMesh->getRegion(cid);

                                 Kokkos::View<kmds::TCellID *> nodes;
                                 c.nodeIds(nodes);
                                 const int nbNodes = nodes.size();

                                 for (auto i_n = 0; i_n < nbNodes; i_n++) {
                                     markedNodes(nodes[i_n]) = true;
                                 }
                             });
    } else {
        Kokkos::parallel_for(ACells2Refine->getNbElems(),
                             KOKKOS_LAMBDA(const kmds::TCellID i) {
                                 const kmds::TCellID nid = ACells2Refine->get(i);
                                 markedNodes(nid) = true;
                             });

        Kokkos::parallel_for(nbCells,
                             KOKKOS_LAMBDA(const kmds::TCellID i) {
                                 const kmds::TCellID cid = cellIDs.get(i);

                                 const kmds::Region c = AMesh->getRegion(cid);

                                 Kokkos::View<kmds::TCellID*> nodes;
                                 c.nodeIds(nodes);
                                 const int nbNodes = nodes.size();

                                 bool isMarked = false;

                                 for(auto i_n=0; i_n<nbNodes; i_n++) {
                                     if(markedNodes(nodes[i_n])) {
                                         markedCells(cid) = true;
                                         isMarked = true;
                                         break;
                                     }
                                 }
                             });
    }

    Refinement_avoidBadConfig_3D(markedCells, markedNodes, AMesh, c_C2C_byN);

    // mark the cells that have marked nodes and permute theirs nodes
    Kokkos::parallel_for(nbCells,
                         KOKKOS_LAMBDA(const kmds::TCellID i) {
                             const kmds::TCellID cid = cellIDs.get(i);

                             const kmds::Region c = AMesh->getRegion(cid);

                             Kokkos::View<kmds::TCellID*> nodes;
                             c.nodeIds(nodes);
                             const int nbNodes = nodes.size();

                             bool isMarked = false;

                             for(auto i_n=0; i_n<nbNodes; i_n++) {
                                 if(markedNodes(nodes[i_n])) {
                                     markedCells(cid) = true;
                                     isMarked = true;
                                     break;
                                 }
                             }

                             if(isMarked) {

                                 // permute the nodes of the cell
                                 elg3d::tableMarkedNodes<8> tableMark(markedNodes[nodes[0]], markedNodes[nodes[1]],
                                                                      markedNodes[nodes[2]], markedNodes[nodes[3]],
                                                                      markedNodes[nodes[4]], markedNodes[nodes[5]],
                                                                      markedNodes[nodes[6]], markedNodes[nodes[7]]);

                                 std::array<int, 8> permutIDs = Refinement_tables::permutNodesOfRegion_3b_3D[tableMark];

                                 for (unsigned int i_n = 0; i_n < 8; i_n++) {
                                     permutIDs[i_n] = nodes[permutIDs[i_n]];
                                 }

                                 for (unsigned int i_n = 0; i_n < 8; i_n++) {
                                     nodes[i_n] = permutIDs[i_n];
                                 }
                             }
                         });


    // marked entities selection
//    kmds::GrowingView <kmds::TCellID> markedCellIDs("CELLS", nbCells);
//    const kmds::TSize nbFaces = AMesh->getNbFaces();
//    kmds::GrowingView <kmds::TCellID> markedFaceIDs("CELLS", nbFaces);
//    const kmds::TSize nbEdges = AMesh->getNbEdges();
//    kmds::GrowingView <kmds::TCellID> markedEdgeIDs("CELLS", nbEdges);

//    const kmds::TCellID nbMarkedCells = markedCellIDs.getNbElems();
//    const kmds::TCellID nbMarkedFaces = markedFaceIDs.getNbElems();
//    const kmds::TCellID nbMarkedEdges = markedEdgeIDs.getNbElems();

    // storage of the new nodes created by edges and faces for the cells that must be refined
    const kmds::TCellID maxID_c = AMesh->getRegionSupID();
    Kokkos::View<kmds::TSize *> cell2nodesExt("markedCells", maxID_c+1);

    kmds::TCellID nbMarkedCells = 0;
    for(kmds::TSize i_c=0; i_c<=maxID_c; i_c++) {
        if(markedCells[i_c]) {
            cell2nodesExt(i_c) = nbMarkedCells;
            nbMarkedCells++;
        }
    }

        std::cout<<"nbMarkedCells "<<nbMarkedCells<<std::endl;

        // storage initialized to NullID
        Kokkos::View<kmds::TCellID *[48]> nodesExt("markedCells", nbMarkedCells);
        Kokkos::parallel_for(nbMarkedCells,
                             KOKKOS_LAMBDA(const kmds::TCellID i) {
                                 for(int i_n=0; i_n<48; i_n++) {
                                     nodesExt(i,i_n) = kmds::NullID;
                                 }
                             });

        // create nodes for edges
        Refinement_createNodesEdges_3D(markedEdges,
                                       markedNodes,
                                       AMesh,
                                       c_E2C,
                                       cell2nodesExt,
                                       nodesExt);


        // create nodes for faces
        Refinement_createNodesFaces_3D(markedFaces,
                                       markedNodes,
                                       AMesh,
                                       c_F2C,
                                       cell2nodesExt,
                                       nodesExt);


        // update the mesh capacity for nodes and cells
        kmds::TSize nbNewNodes = 0;
        kmds::TSize nbNewCells = 0;
        Kokkos::parallel_reduce(maxID_c+1,
                                KOKKOS_LAMBDA(const kmds::TCellID i, kmds::TSize& sum) {

                                    const kmds::TCellID cid = i;

                                    if (markedCells[cid]) {
                                        const kmds::Region c = AMesh->getRegion(cid);

                                        Kokkos::View<kmds::TCellID *> nodes;
                                        c.nodeIds(nodes);

                                        elg3d::tableMarkedNodes<8> tableMark(markedNodes[nodes[0]], markedNodes[nodes[1]],
                                                                             markedNodes[nodes[2]], markedNodes[nodes[3]],
                                                                             markedNodes[nodes[4]], markedNodes[nodes[5]],
                                                                             markedNodes[nodes[6]], markedNodes[nodes[7]]);

                                        sum += Refinement_tables::nbNodesOnRegionToBuild_3b_3D[tableMark];
                                    }
                                },
                                nbNewNodes);

        Kokkos::parallel_reduce(maxID_c+1,
                                KOKKOS_LAMBDA(const kmds::TCellID i, kmds::TSize& sum) {

                                    const kmds::TCellID cid = i;

                                    if (markedCells[cid]) {
                                        const kmds::Region c = AMesh->getRegion(cid);

                                        Kokkos::View<kmds::TCellID *> nodes;
                                        c.nodeIds(nodes);

                                        elg3d::tableMarkedNodes<8> tableMark(markedNodes[nodes[0]], markedNodes[nodes[1]],
                                                                             markedNodes[nodes[2]], markedNodes[nodes[3]],
                                                                             markedNodes[nodes[4]], markedNodes[nodes[5]],
                                                                             markedNodes[nodes[6]], markedNodes[nodes[7]]);

                                        sum += Refinement_tables::nbRegionsOnRegionToBuild_3b_3D[tableMark];
                                    }
                                },
                                nbNewCells);

        AMesh->updateNodeCapacity(AMesh->getNbNodes() + nbNewNodes);
        AMesh->updateRegionCapacity(AMesh->getNbRegions() + nbNewCells);


    // nodes for cells and the new cells
    Kokkos::parallel_for(nbCells,
                         KOKKOS_LAMBDA(const kmds::TCellID i) {
                             const kmds::TCellID cid = cellIDs.get(i);

                             if(markedCells[cid]) {

                                 const kmds::Region c = AMesh->getRegion(cid);

                                 Kokkos::View<kmds::TCellID *> nodes;
                                 c.nodeIds(nodes);
                                 const int nbNodes = nodes.size();

                                 elg3d::tableMarkedNodes<8> tableMark(markedNodes[nodes[0]], markedNodes[nodes[1]],
                                                                      markedNodes[nodes[2]], markedNodes[nodes[3]],
                                                                      markedNodes[nodes[4]], markedNodes[nodes[5]],
                                                                      markedNodes[nodes[6]], markedNodes[nodes[7]]);

                                 // regroup the nodes :
                                 kmds::TCellID subNodeIDs[7 * 7 * 7];

                                 // nodes of the region
                                 subNodeIDs[0 * 7 * 7 + 0 * 7 + 0] = nodes[0];
                                 subNodeIDs[6 * 7 * 7 + 0 * 7 + 0] = nodes[1];
                                 subNodeIDs[6 * 7 * 7 + 6 * 7 + 0] = nodes[2];
                                 subNodeIDs[0 * 7 * 7 + 6 * 7 + 0] = nodes[3];
                                 subNodeIDs[0 * 7 * 7 + 0 * 7 + 6] = nodes[4];
                                 subNodeIDs[6 * 7 * 7 + 0 * 7 + 6] = nodes[5];
                                 subNodeIDs[6 * 7 * 7 + 6 * 7 + 6] = nodes[6];
                                 subNodeIDs[0 * 7 * 7 + 6 * 7 + 6] = nodes[7];

                                 // nodes created by edges and faces
                                 for (int i_n = 0; i_n < 48; i_n++) {
                                     const int index = Refinement_tables::nodesOnRegionFromExternal_3b_3D[i_n];
                                     subNodeIDs[index] = nodesExt(cell2nodesExt[cid], i_n);
                                 }

                                 // new nodes built exclusively for this region
                                 const int nbNewNodes = Refinement_tables::nbNodesOnRegionToBuild_3b_3D[tableMark];

                                 kmds::TCellID firstIDN = AMesh->addNodes(nbNewNodes);

                                 for (int i_n = 0; i_n < nbNewNodes; i_n++) {
                                     const double ix = Refinement_tables::nodesOnRegionToBuild_3b_3D[tableMark][i_n][0];
                                     const double iy = Refinement_tables::nodesOnRegionToBuild_3b_3D[tableMark][i_n][1];
                                     const double iz = Refinement_tables::nodesOnRegionToBuild_3b_3D[tableMark][i_n][2];
                                     gmds::math::Point pt = (6 - ix) / 6. * (6 - iy) / 6. * (6 - iz) / 6. *
                                                            AMesh->getNodeLocation(nodes[0])
                                                            + (ix) / 6. * (6 - iy) / 6. * (6 - iz) / 6. *
                                                              AMesh->getNodeLocation(nodes[1])
                                                            + (ix) / 6. * (iy) / 6. * (6 - iz) / 6. *
                                                              AMesh->getNodeLocation(nodes[2])
                                                            + (6 - ix) / 6. * (iy) / 6. * (6 - iz) / 6. *
                                                              AMesh->getNodeLocation(nodes[3])
                                                            + (6 - ix) / 6. * (6 - iy) / 6. * (iz) / 6. *
                                                              AMesh->getNodeLocation(nodes[4])
                                                            + (ix) / 6. * (6 - iy) / 6. * (iz) / 6. *
                                                              AMesh->getNodeLocation(nodes[5])
                                                            + (ix) / 6. * (iy) / 6. * (iz) / 6. *
                                                              AMesh->getNodeLocation(nodes[6])
                                                            + (6 - ix) / 6. * (iy) / 6. * (iz) / 6. *
                                                              AMesh->getNodeLocation(nodes[7]);

                                     const int index = ix * 7 * 7 + iy * 7 + iz;
                                     subNodeIDs[index] = firstIDN + i_n;

                                     AMesh->setNodeLocation(firstIDN + i_n, pt);
                                 }

                                 const int nbNewRegions = Refinement_tables::nbRegionsOnRegionToBuild_3b_3D[tableMark];
                                 std::vector<std::array<int, 8> > regions2Build = Refinement_tables::regionsOnRegionToBuild_3b_3D[tableMark];


                                 // the first one will take the place of the current region
                                 nodes[0] = subNodeIDs[regions2Build[0][0]];
                                 nodes[1] = subNodeIDs[regions2Build[0][1]];
                                 nodes[2] = subNodeIDs[regions2Build[0][2]];
                                 nodes[3] = subNodeIDs[regions2Build[0][3]];
                                 nodes[4] = subNodeIDs[regions2Build[0][4]];
                                 nodes[5] = subNodeIDs[regions2Build[0][5]];
                                 nodes[6] = subNodeIDs[regions2Build[0][6]];
                                 nodes[7] = subNodeIDs[regions2Build[0][7]];

//                             c.setNodes(nodes);

                                 kmds::TCellID firstIDC = AMesh->addHexahedra(nbNewRegions - 1);

                                 for (int i_c = 1; i_c < nbNewRegions; i_c++) {
                                     kmds::Region c_bis = AMesh->getRegion(firstIDC + i_c - 1);

                                     Kokkos::View<kmds::TCellID *> nodes_bis;
                                     c_bis.nodeIds(nodes_bis);

                                     nodes_bis[0] = subNodeIDs[regions2Build[i_c][0]];
                                     nodes_bis[1] = subNodeIDs[regions2Build[i_c][1]];
                                     nodes_bis[2] = subNodeIDs[regions2Build[i_c][2]];
                                     nodes_bis[3] = subNodeIDs[regions2Build[i_c][3]];
                                     nodes_bis[4] = subNodeIDs[regions2Build[i_c][4]];
                                     nodes_bis[5] = subNodeIDs[regions2Build[i_c][5]];
                                     nodes_bis[6] = subNodeIDs[regions2Build[i_c][6]];
                                     nodes_bis[7] = subNodeIDs[regions2Build[i_c][7]];

//                                 c_bis.setNodes(nodes_bis);
                                 }
                             }  // if(markedCells[cid])
                         });

}

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
