/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    InterfaceNodesPosSmoothVF.cpp
 *  \author  legoff
 *  \date    09/30/2018
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/InterfaceNodesPosSmoothVF.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <set>
#include <random>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <Kokkos_Core.hpp>
#include <Kokkos_UnorderedMap.hpp>

#include <glpk.h>

#include <GCoptimization.h>
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/math/Point.h>
#include <gmds/math/VectorDyn.h>

#include <KM/Utils/Graph.h>
#include <KM/DS/Mesh.h>
#include <KM/IO/VTKWriter.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include <ELG3D/DATACMPNT/FracPres.h>
#include <ELG3D/DATACMPNT/MaterialAssignment.h>
#include <ELG3D/DATACMPNT/Parameters.h>
#include <ELG3D/ALGOCMPNT/Tools.h>
/*----------------------------------------------------------------------------*/
using namespace elg3d;
namespace milp {

    /*----------------------------------------------------------------------------*/
    void
    InterfaceNodesPosSmoothVF_getNodeSubIndices_2D(const kmds::Face AF,
                                                   const kmds::TCellID AID,
                                                   int *AI,
                                                   int *AJ) {
        Kokkos::View<kmds::TCellID *> nids;
        AF.nodeIds(nids);

        int localIndex = -1;
        for (int i = 0; i < nids.size(); i++) {
            if (nids[i] == AID) {
                localIndex = i;
            }
        }

        switch (localIndex) {
            case 0 :
                *AI = 0;
                *AJ = 0;
                break;
            case 1 :
                *AI = InterfaceNodesPosSmoothVF_NBSUB;
                *AJ = 0;
                break;
            case 2 :
                *AI = InterfaceNodesPosSmoothVF_NBSUB;
                *AJ = InterfaceNodesPosSmoothVF_NBSUB;
                break;
            case 3 :
                *AI = 0;
                *AJ = InterfaceNodesPosSmoothVF_NBSUB;
                break;

            default:
                break;
        }
    }

    /*----------------------------------------------------------------------------*/
    void
    InterfaceNodesPosSmoothVF_getNodeSubIndices_3D(const kmds::Region AR,
                                                   const kmds::TCellID AID,
                                                   int *AI,
                                                   int *AJ,
                                                   int *AK) {
        Kokkos::View<kmds::TCellID *> nids;
        AR.nodeIds(nids);

        int localIndex = -1;
        for (int i = 0; i < nids.size(); i++) {
            if (nids[i] == AID) {
                localIndex = i;
                break;
            }
        }

        switch (localIndex) {
            case 0 :
                *AI = 0;
                *AJ = 0;
                *AK = 0;
                break;
            case 1 :
                *AI = InterfaceNodesPosSmoothVF_NBSUB;
                *AJ = 0;
                *AK = 0;
                break;
            case 2 :
                *AI = InterfaceNodesPosSmoothVF_NBSUB;
                *AJ = InterfaceNodesPosSmoothVF_NBSUB;
                *AK = 0;
                break;
            case 3 :
                *AI = 0;
                *AJ = InterfaceNodesPosSmoothVF_NBSUB;
                *AK = 0;
                break;
            case 4 :
                *AI = 0;
                *AJ = 0;
                *AK = InterfaceNodesPosSmoothVF_NBSUB;
                break;
            case 5 :
                *AI = InterfaceNodesPosSmoothVF_NBSUB;
                *AJ = 0;
                *AK = InterfaceNodesPosSmoothVF_NBSUB;
                break;
            case 6 :
                *AI = InterfaceNodesPosSmoothVF_NBSUB;
                *AJ = InterfaceNodesPosSmoothVF_NBSUB;
                *AK = InterfaceNodesPosSmoothVF_NBSUB;
                break;
            case 7 :
                *AI = 0;
                *AJ = InterfaceNodesPosSmoothVF_NBSUB;
                *AK = InterfaceNodesPosSmoothVF_NBSUB;
                break;
            default:
                break;
        }
    }

    /*----------------------------------------------------------------------------*/
    void
    InterfaceNodesPosSmoothVF_getEdgeSubIndices_2D(const kmds::Face AF,
                                                   const kmds::TCellID AID1,
                                                   const kmds::TCellID AID2,
                                                   int *AImin,
                                                   int *AImax,
                                                   int *AJmin,
                                                   int *AJmax,
                                                   bool *direct) {
        Kokkos::View<kmds::TCellID *> nids;
        AF.nodeIds(nids);

        int edgeIndex = -1;

        // IJ_direct
        if ((nids[0] == AID1) && (nids[1] == AID2)) {
            edgeIndex = 0;
            *direct = true;
        }
        if ((nids[1] == AID1) && (nids[2] == AID2)) {
            edgeIndex = 1;
            *direct = true;
        }
        if ((nids[3] == AID1) && (nids[2] == AID2)) {
            edgeIndex = 2;
            *direct = true;
        }
        if ((nids[0] == AID1) && (nids[3] == AID2)) {
            edgeIndex = 3;
            *direct = true;
        }

        // IJ_indirect
        if ((nids[0] == AID2) && (nids[1] == AID1)) {
            edgeIndex = 0;
            *direct = false;
        }
        if ((nids[1] == AID2) && (nids[2] == AID1)) {
            edgeIndex = 1;
            *direct = false;
        }
        if ((nids[3] == AID2) && (nids[2] == AID1)) {
            edgeIndex = 2;
            *direct = false;
        }
        if ((nids[0] == AID2) && (nids[3] == AID1)) {
            edgeIndex = 3;
            *direct = false;
        }

        *AImin = -1;
        *AImax = -1;
        *AJmin = -1;
        *AJmax = -1;

        switch (edgeIndex) {
            case 0 :
                *AImin = 1;
                *AImax = InterfaceNodesPosSmoothVF_NBSUB - 1;
                *AJmin = 0;
                *AJmax = 0;
                break;
            case 1 :
                *AImin = InterfaceNodesPosSmoothVF_NBSUB;
                *AImax = InterfaceNodesPosSmoothVF_NBSUB;
                *AJmin = 1;
                *AJmax = InterfaceNodesPosSmoothVF_NBSUB - 1;
                break;
            case 2 :
                *AImin = 1;
                *AImax = InterfaceNodesPosSmoothVF_NBSUB - 1;
                *AJmin = InterfaceNodesPosSmoothVF_NBSUB;
                *AJmax = InterfaceNodesPosSmoothVF_NBSUB;
                break;
            case 3 :
                *AImin = 0;
                *AImax = 0;
                *AJmin = 1;
                *AJmax = InterfaceNodesPosSmoothVF_NBSUB - 1;
                break;
            default:
                break;
        }
    }

    /*----------------------------------------------------------------------------*/
    void
    InterfaceNodesPosSmoothVF_getEdgeSubIndices_3D(const kmds::Region AR,
                                                   const kmds::TCellID AID1,
                                                   const kmds::TCellID AID2,
                                                   int *AImin,
                                                   int *AImax,
                                                   int *AJmin,
                                                   int *AJmax,
                                                   int *AKmin,
                                                   int *AKmax,
                                                   bool *direct) {
        Kokkos::View<kmds::TCellID *> nids;
        AR.nodeIds(nids);

        int edgeIndex = -1;

        // IJK_direct
        // bottom
        if ((nids[0] == AID1) && (nids[1] == AID2)) {
            edgeIndex = 0;
            *direct = true;
        }
        if ((nids[1] == AID1) && (nids[2] == AID2)) {
            edgeIndex = 1;
            *direct = true;
        }
        if ((nids[3] == AID1) && (nids[2] == AID2)) {
            edgeIndex = 2;
            *direct = true;
        }
        if ((nids[0] == AID1) && (nids[3] == AID2)) {
            edgeIndex = 3;
            *direct = true;
        }

        // top
        if ((nids[4] == AID1) && (nids[5] == AID2)) {
            edgeIndex = 4;
            *direct = true;
        }
        if ((nids[5] == AID1) && (nids[6] == AID2)) {
            edgeIndex = 5;
            *direct = true;
        }
        if ((nids[7] == AID1) && (nids[6] == AID2)) {
            edgeIndex = 6;
            *direct = true;
        }
        if ((nids[4] == AID1) && (nids[7] == AID2)) {
            edgeIndex = 7;
            *direct = true;
        }

        // side
        if ((nids[0] == AID1) && (nids[4] == AID2)) {
            edgeIndex = 8;
            *direct = true;
        }
        if ((nids[1] == AID1) && (nids[5] == AID2)) {
            edgeIndex = 9;
            *direct = true;
        }
        if ((nids[2] == AID1) && (nids[6] == AID2)) {
            edgeIndex = 10;
            *direct = true;
        }
        if ((nids[3] == AID1) && (nids[7] == AID2)) {
            edgeIndex = 11;
            *direct = true;
        }

        // IJ_indirect
        // bottom
        if ((nids[0] == AID2) && (nids[1] == AID1)) {
            edgeIndex = 0;
            *direct = false;
        }
        if ((nids[1] == AID2) && (nids[2] == AID1)) {
            edgeIndex = 1;
            *direct = false;
        }
        if ((nids[3] == AID2) && (nids[2] == AID1)) {
            edgeIndex = 2;
            *direct = false;
        }
        if ((nids[0] == AID2) && (nids[3] == AID1)) {
            edgeIndex = 3;
            *direct = false;
        }

        // top
        if ((nids[4] == AID2) && (nids[5] == AID1)) {
            edgeIndex = 4;
            *direct = false;
        }
        if ((nids[5] == AID2) && (nids[6] == AID1)) {
            edgeIndex = 5;
            *direct = false;
        }
        if ((nids[7] == AID2) && (nids[6] == AID1)) {
            edgeIndex = 6;
            *direct = false;
        }
        if ((nids[4] == AID2) && (nids[7] == AID1)) {
            edgeIndex = 7;
            *direct = false;
        }

        // side
        if ((nids[0] == AID2) && (nids[4] == AID1)) {
            edgeIndex = 8;
            *direct = false;
        }
        if ((nids[1] == AID2) && (nids[5] == AID1)) {
            edgeIndex = 9;
            *direct = false;
        }
        if ((nids[2] == AID2) && (nids[6] == AID1)) {
            edgeIndex = 10;
            *direct = false;
        }
        if ((nids[3] == AID2) && (nids[7] == AID1)) {
            edgeIndex = 11;
            *direct = false;
        }

        *AImin = -1;
        *AImax = -1;
        *AJmin = -1;
        *AJmax = -1;
        *AKmin = -1;
        *AKmax = -1;

        switch (edgeIndex) {
            case 0 :
                *AImin = 1;
                *AImax = InterfaceNodesPosSmoothVF_NBSUB - 1;
                *AJmin = 0;
                *AJmax = 0;
                *AKmin = 0;
                *AKmax = 0;
                break;
            case 1 :
                *AImin = InterfaceNodesPosSmoothVF_NBSUB;
                *AImax = InterfaceNodesPosSmoothVF_NBSUB;
                *AJmin = 1;
                *AJmax = InterfaceNodesPosSmoothVF_NBSUB - 1;
                *AKmin = 0;
                *AKmax = 0;
                break;
            case 2 :
                *AImin = 1;
                *AImax = InterfaceNodesPosSmoothVF_NBSUB - 1;
                *AJmin = InterfaceNodesPosSmoothVF_NBSUB;
                *AJmax = InterfaceNodesPosSmoothVF_NBSUB;
                *AKmin = 0;
                *AKmax = 0;
                break;
            case 3 :
                *AImin = 0;
                *AImax = 0;
                *AJmin = 1;
                *AJmax = InterfaceNodesPosSmoothVF_NBSUB - 1;
                *AKmin = 0;
                *AKmax = 0;
                break;
            case 4 :
                *AImin = 1;
                *AImax = InterfaceNodesPosSmoothVF_NBSUB - 1;
                *AJmin = 0;
                *AJmax = 0;
                *AKmin = InterfaceNodesPosSmoothVF_NBSUB;
                *AKmax = InterfaceNodesPosSmoothVF_NBSUB;
                break;
            case 5 :
                *AImin = InterfaceNodesPosSmoothVF_NBSUB;
                *AImax = InterfaceNodesPosSmoothVF_NBSUB;
                *AJmin = 1;
                *AJmax = InterfaceNodesPosSmoothVF_NBSUB - 1;
                *AKmin = InterfaceNodesPosSmoothVF_NBSUB;
                *AKmax = InterfaceNodesPosSmoothVF_NBSUB;
                break;
            case 6 :
                *AImin = 1;
                *AImax = InterfaceNodesPosSmoothVF_NBSUB - 1;
                *AJmin = InterfaceNodesPosSmoothVF_NBSUB;
                *AJmax = InterfaceNodesPosSmoothVF_NBSUB;
                *AKmin = InterfaceNodesPosSmoothVF_NBSUB;
                *AKmax = InterfaceNodesPosSmoothVF_NBSUB;
                break;
            case 7 :
                *AImin = 0;
                *AImax = 0;
                *AJmin = 1;
                *AJmax = InterfaceNodesPosSmoothVF_NBSUB - 1;
                *AKmin = InterfaceNodesPosSmoothVF_NBSUB;
                *AKmax = InterfaceNodesPosSmoothVF_NBSUB;
                break;
            case 8 :
                *AImin = 0;
                *AImax = 0;
                *AJmin = 0;
                *AJmax = 0;
                *AKmin = 1;
                *AKmax = InterfaceNodesPosSmoothVF_NBSUB - 1;
                break;
            case 9 :
                *AImin = InterfaceNodesPosSmoothVF_NBSUB;
                *AImax = InterfaceNodesPosSmoothVF_NBSUB;
                *AJmin = 0;
                *AJmax = 0;
                *AKmin = 1;
                *AKmax = InterfaceNodesPosSmoothVF_NBSUB - 1;
                break;
            case 10 :
                *AImin = InterfaceNodesPosSmoothVF_NBSUB;
                *AImax = InterfaceNodesPosSmoothVF_NBSUB;
                *AJmin = InterfaceNodesPosSmoothVF_NBSUB;
                *AJmax = InterfaceNodesPosSmoothVF_NBSUB;
                *AKmin = 1;
                *AKmax = InterfaceNodesPosSmoothVF_NBSUB - 1;
                break;
            case 11 :
                *AImin = 0;
                *AImax = 0;
                *AJmin = InterfaceNodesPosSmoothVF_NBSUB;
                *AJmax = InterfaceNodesPosSmoothVF_NBSUB;
                *AKmin = 1;
                *AKmax = InterfaceNodesPosSmoothVF_NBSUB - 1;
                break;
            default:
                break;
        }
    }

    /*----------------------------------------------------------------------------*/
    void
    InterfaceNodesPosSmoothVF_getFaceSubIndices_3D(const kmds::Region AR,
                                                   const kmds::TCellID AID1,
                                                   const kmds::TCellID AID2,
                                                   const kmds::TCellID AID3,
                                                   const kmds::TCellID AID4,
                                                   int *AImin,
                                                   int *AImax,
                                                   int *AJmin,
                                                   int *AJmax,
                                                   int *AKmin,
                                                   int *AKmax,
                                                   int *AFaceRel) {

        Kokkos::View<kmds::TCellID *> nids;
        AR.nodeIds(nids);

        // find which facet it is
        int numFacet = -1;

        if ((AID1 == nids[0] || AID1 == nids[1] || AID1 == nids[2] || AID1 == nids[3])
            && (AID2 == nids[0] || AID2 == nids[1] || AID2 == nids[2] || AID2 == nids[3])
            && (AID3 == nids[0] || AID3 == nids[1] || AID3 == nids[2] || AID3 == nids[3])
            && (AID4 == nids[0] || AID4 == nids[1] || AID4 == nids[2] || AID4 == nids[3])) {
            // bottom
            numFacet = 0;

            *AImin = 1;
            *AImax = InterfaceNodesPosSmoothVF_NBSUB - 1;
            *AJmin = 1;
            *AJmax = InterfaceNodesPosSmoothVF_NBSUB - 1;
            *AKmin = 0;
            *AKmax = 0;

        } else if ((AID1 == nids[4] || AID1 == nids[5] || AID1 == nids[6] || AID1 == nids[7])
                   && (AID2 == nids[4] || AID2 == nids[5] || AID2 == nids[6] || AID2 == nids[7])
                   && (AID3 == nids[4] || AID3 == nids[5] || AID3 == nids[6] || AID3 == nids[7])
                   && (AID4 == nids[4] || AID4 == nids[5] || AID4 == nids[6] || AID4 == nids[7])) {
            // top
            numFacet = 1;

            *AImin = 1;
            *AImax = InterfaceNodesPosSmoothVF_NBSUB - 1;
            *AJmin = 1;
            *AJmax = InterfaceNodesPosSmoothVF_NBSUB - 1;
            *AKmin = InterfaceNodesPosSmoothVF_NBSUB;
            *AKmax = InterfaceNodesPosSmoothVF_NBSUB;

        } else if ((AID1 == nids[0] || AID1 == nids[3] || AID1 == nids[7] || AID1 == nids[4])
                   && (AID2 == nids[0] || AID2 == nids[3] || AID2 == nids[7] || AID2 == nids[4])
                   && (AID3 == nids[0] || AID3 == nids[3] || AID3 == nids[7] || AID3 == nids[4])
                   && (AID4 == nids[0] || AID4 == nids[3] || AID4 == nids[7] || AID4 == nids[4])) {
            // left
            numFacet = 2;

            *AImin = 0;
            *AImax = 0;
            *AJmin = 1;
            *AJmax = InterfaceNodesPosSmoothVF_NBSUB - 1;
            *AKmin = 1;
            *AKmax = InterfaceNodesPosSmoothVF_NBSUB - 1;

        } else if ((AID1 == nids[1] || AID1 == nids[2] || AID1 == nids[6] || AID1 == nids[5])
                   && (AID2 == nids[1] || AID2 == nids[2] || AID2 == nids[6] || AID2 == nids[5])
                   && (AID3 == nids[1] || AID3 == nids[2] || AID3 == nids[6] || AID3 == nids[5])
                   && (AID4 == nids[1] || AID4 == nids[2] || AID4 == nids[6] || AID4 == nids[5])) {
            // right
            numFacet = 3;

            *AImin = InterfaceNodesPosSmoothVF_NBSUB;
            *AImax = InterfaceNodesPosSmoothVF_NBSUB;
            *AJmin = 1;
            *AJmax = InterfaceNodesPosSmoothVF_NBSUB - 1;
            *AKmin = 1;
            *AKmax = InterfaceNodesPosSmoothVF_NBSUB - 1;

        } else if ((AID1 == nids[0] || AID1 == nids[1] || AID1 == nids[5] || AID1 == nids[4])
                   && (AID2 == nids[0] || AID2 == nids[1] || AID2 == nids[5] || AID2 == nids[4])
                   && (AID3 == nids[0] || AID3 == nids[1] || AID3 == nids[5] || AID3 == nids[4])
                   && (AID4 == nids[0] || AID4 == nids[1] || AID4 == nids[5] || AID4 == nids[4])) {
            // front
            numFacet = 4;

            *AImin = 1;
            *AImax = InterfaceNodesPosSmoothVF_NBSUB - 1;
            *AJmin = 0;
            *AJmax = 0;
            *AKmin = 1;
            *AKmax = InterfaceNodesPosSmoothVF_NBSUB - 1;

        } else if ((AID1 == nids[3] || AID1 == nids[2] || AID1 == nids[6] || AID1 == nids[7])
                   && (AID2 == nids[3] || AID2 == nids[2] || AID2 == nids[6] || AID2 == nids[7])
                   && (AID3 == nids[3] || AID3 == nids[2] || AID3 == nids[6] || AID3 == nids[7])
                   && (AID4 == nids[3] || AID4 == nids[2] || AID4 == nids[6] || AID4 == nids[7])) {
            // back
            numFacet = 5;

            *AImin = 1;
            *AImax = InterfaceNodesPosSmoothVF_NBSUB - 1;
            *AJmin = InterfaceNodesPosSmoothVF_NBSUB;
            *AJmax = InterfaceNodesPosSmoothVF_NBSUB;
            *AKmin = 1;
            *AKmax = InterfaceNodesPosSmoothVF_NBSUB - 1;
        }

        // find the relation between the face and its corresponding region facet
        const std::vector<std::vector<kmds::TCellID> > ijkNodes = AR.getIJKNodesFacesIDs();

        bool direct = true;
        int offset = -1;

        for (int i=0; i<ijkNodes[numFacet].size(); i++) {
            if(ijkNodes[numFacet][i] == AID1) {
                offset = i;

                if(ijkNodes[numFacet][(i + 1) % ijkNodes[numFacet].size()] != AID2) {
                    direct = false;
                    break;
                }
            }
        }

        *AFaceRel = offset;
        if(!direct) {
            *AFaceRel += 4;
        }

    }

    /*----------------------------------------------------------------------------*/
    struct InterfaceNodesPosSmoothVF_identifyCells_xD {
        const kmds::GrowingView<kmds::TCellID> *selection;
        kmds::Mesh *mesh;
        const kmds::Connectivity *c_N2C;
        const FracPres *fp;

        Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> cellsSet;

        kmds::GrowingView<kmds::TCellID> *selectionCells;

        InterfaceNodesPosSmoothVF_identifyCells_xD(

                const kmds::GrowingView<kmds::TCellID> *selection_,
                kmds::Mesh *mesh_,
                const kmds::Connectivity *c_N2C_,
                const FracPres *fp_,

                Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> cellsSet_,

                kmds::GrowingView<kmds::TCellID> *selectionCells_

        )
                : selection(selection_), mesh(mesh_), c_N2C(c_N2C_), fp(fp_), cellsSet(cellsSet_),
                  selectionCells(selectionCells_) {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            const int nid = selection->get(i);

            Kokkos::View<kmds::TCellID *> cids;
            c_N2C->get(nid, cids);

            bool toKeep = false;
            for (int ic = 0; ic < cids.size(); ic++) {
                if (fp->isMixedCell(cids[ic])) {
                    toKeep = true;
                }
            }

            // keep the adjacent cells if there is a mixed cell
            if (toKeep) {

                for (int ic = 0; ic < cids.size(); ic++) {
                    Kokkos::UnorderedMapInsertResult res = cellsSet.insert(cids[ic]);
                    bool success = res.success();

                    if (success) {
                        selectionCells->push_back(cids[ic]);
                    }
                }

            }
        }

    };

/*----------------------------------------------------------------------------*/
    struct InterfaceNodesPosSmoothVF_identifyNodesFromN2C_xD {

        const kmds::GrowingView<kmds::TCellID> *selection;
        const kmds::Connectivity *c_N2C;
        const kmds::Variable<bool> *varCellsKept;

        kmds::GrowingView<kmds::TCellID> *selectionNodes;

        InterfaceNodesPosSmoothVF_identifyNodesFromN2C_xD(

                const kmds::GrowingView<kmds::TCellID> *selection_,
                const kmds::Connectivity *c_N2C_,
                const kmds::Variable<bool> *varCellsKept_,

                kmds::GrowingView<kmds::TCellID> *selectionNodes_

        )
                : selection(selection_), c_N2C(c_N2C_), varCellsKept(varCellsKept_), selectionNodes(selectionNodes_) {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            const kmds::TCellID nid = selection->get(i);

            Kokkos::View<kmds::TCellID *> cids;
            c_N2C->get(nid, cids);

            bool toKeep = false;
            for (int ic = 0; ic < cids.size(); ic++) {
                if ((*varCellsKept)[cids[ic]]) {
                    toKeep = true;
                }
            }

            if (toKeep) {
                selectionNodes->push_back(nid);
            }
        }

    };

    /*----------------------------------------------------------------------------*/
    struct InterfaceNodesPosSmoothVF_identifyEdgesFromE2C_xD {

        const kmds::GrowingView<kmds::TCellID> *selection;
        const kmds::Connectivity *c_E2C;
        const kmds::Variable<bool> *varCellsKept;

        kmds::GrowingView<kmds::TCellID> *selectionEdges;

        InterfaceNodesPosSmoothVF_identifyEdgesFromE2C_xD(

                const kmds::GrowingView<kmds::TCellID> *selection_,
                const kmds::Connectivity *c_E2C_,
                const kmds::Variable<bool> *varCellsKept_,

                kmds::GrowingView<kmds::TCellID> *selectionEdges_

        )
                : selection(selection_), c_E2C(c_E2C_), varCellsKept(varCellsKept_), selectionEdges(selectionEdges_) {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            const kmds::TCellID eid = selection->get(i);

            Kokkos::View<kmds::TCellID *> cids;
            c_E2C->get(eid, cids);

            for(int ic=0; ic<cids.size(); ic++) {
                if ((*varCellsKept)[cids[ic]]) {
                    selectionEdges->push_back(eid);
                    break;
                }
            }

        }

    };

    /*----------------------------------------------------------------------------*/
    struct InterfaceNodesPosSmoothVF_identifyFacesFromF2C_xD {

        const kmds::GrowingView<kmds::TCellID> *selection;
        const kmds::Connectivity *c_F2C;
        const kmds::Variable<bool> *varCellsKept;

        kmds::GrowingView<kmds::TCellID> *selectionFaces;

        InterfaceNodesPosSmoothVF_identifyFacesFromF2C_xD(

                const kmds::GrowingView<kmds::TCellID> *selection_,
                const kmds::Connectivity *c_F2C_,
                const kmds::Variable<bool> *varCellsKept_,

                kmds::GrowingView<kmds::TCellID> *selectionFaces_

        )
                : selection(selection_), c_F2C(c_F2C_), varCellsKept(varCellsKept_), selectionFaces(selectionFaces_) {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            const kmds::TCellID fid = selection->get(i);

            Kokkos::View<kmds::TCellID *> cids;
            c_F2C->get(fid, cids);

            for(int ic=0; ic<cids.size(); ic++) {
                if ((*varCellsKept)[cids[ic]]) {
                    selectionFaces->push_back(fid);
                    break;
                }
            }

        }

    };

/*----------------------------------------------------------------------------*/
    struct InterfaceNodesPosSmoothVF_createNodesFromNodes_2D {
        const kmds::GrowingView<kmds::TCellID> *selection;
        kmds::Mesh *mesh_source;
        kmds::Mesh *mesh_target;
        const kmds::Connectivity *c_N2C;
        const kmds::Variable<kmds::TCellID> *varOldCells2cellsSelection;

        Kokkos::View<kmds::TCellID *[InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB +
                                                                          1]> oldCell2subNodes;

        InterfaceNodesPosSmoothVF_createNodesFromNodes_2D(

                const kmds::GrowingView<kmds::TCellID> *selection_,
                kmds::Mesh *mesh_source_,
                kmds::Mesh *mesh_target_,
                const kmds::Connectivity *c_N2C_,
                const kmds::Variable<kmds::TCellID> *varOldCells2cellsSelection_,
                Kokkos::View<kmds::TCellID *[InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB +
                                                                                  1]> oldCell2subNodes_
        )
                : selection(selection_), mesh_source(mesh_source_), mesh_target(mesh_target_), c_N2C(c_N2C_),
                  varOldCells2cellsSelection(varOldCells2cellsSelection_), oldCell2subNodes(oldCell2subNodes_) {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            const int nid = selection->get(i);

            double xyz[3];
            mesh_source->getNodeLocation(nid, xyz[0], xyz[1], xyz[2]);
            const kmds::TCellID subID = mesh_target->addNode();
            mesh_target->setNodeLocation(subID, xyz[0], xyz[1], xyz[2]);

            Kokkos::View<kmds::TCellID *> cids;
            c_N2C->get(nid, cids);

            for (int ic = 0; ic < cids.size(); ic++) {

                const kmds::TCellID oldcid = cids[ic];
                if ((*varOldCells2cellsSelection)[oldcid] != kmds::NullID) {

                    const kmds::Face f = mesh_source->getFace(oldcid);

                    const kmds::TCellID subCid = (*varOldCells2cellsSelection)[oldcid];
                    int ix = -1;
                    int iy = -1;
                    InterfaceNodesPosSmoothVF_getNodeSubIndices_2D(f, nid, &ix, &iy);
                    oldCell2subNodes(subCid, ix, iy) = subID;
                }
            }

        }

    };

    /*----------------------------------------------------------------------------*/
    struct InterfaceNodesPosSmoothVF_createNodesFromNodes_3D {
        const kmds::GrowingView<kmds::TCellID> *selection;
        kmds::Mesh *mesh_source;
        kmds::Mesh *mesh_target;
        const kmds::Connectivity *c_N2C;
        const kmds::Variable<kmds::TCellID> *varOldCells2cellsSelection;

        Kokkos::View<kmds::TCellID *[InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB +
                                                                          1]> oldCell2subNodes;

        InterfaceNodesPosSmoothVF_createNodesFromNodes_3D(

                const kmds::GrowingView<kmds::TCellID> *selection_,
                kmds::Mesh *mesh_source_,
                kmds::Mesh *mesh_target_,
                const kmds::Connectivity *c_N2C_,
                const kmds::Variable<kmds::TCellID> *varOldCells2cellsSelection_,
                Kokkos::View<kmds::TCellID *[InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB +
                                                                                  1]> oldCell2subNodes_
        )
                : selection(selection_), mesh_source(mesh_source_), mesh_target(mesh_target_), c_N2C(c_N2C_),
                  varOldCells2cellsSelection(varOldCells2cellsSelection_), oldCell2subNodes(oldCell2subNodes_) {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            const int nid = selection->get(i);

            double xyz[3];
            mesh_source->getNodeLocation(nid, xyz[0], xyz[1], xyz[2]);
            const kmds::TCellID subID = mesh_target->addNode();
            mesh_target->setNodeLocation(subID, xyz[0], xyz[1], xyz[2]);

            Kokkos::View<kmds::TCellID *> cids;
            c_N2C->get(nid, cids);

            for (int ic = 0; ic < cids.size(); ic++) {

                const kmds::TCellID oldcid = cids[ic];
                if ((*varOldCells2cellsSelection)[oldcid] != kmds::NullID) {

                    const kmds::Region r = mesh_source->getRegion(oldcid);

                    const kmds::TCellID subCid = (*varOldCells2cellsSelection)[oldcid];
                    int ix = -1;
                    int iy = -1;
                    int iz = -1;
                    InterfaceNodesPosSmoothVF_getNodeSubIndices_3D(r, nid, &ix, &iy, &iz);
                    oldCell2subNodes(subCid, ix, iy, iz) = subID;

                }
            }

        }

    };

    /*----------------------------------------------------------------------------*/
    struct InterfaceNodesPosSmoothVF_createNodesFromEdges_2D {
        const kmds::GrowingView<kmds::TCellID> *selection;
        kmds::Mesh *mesh_source;
        kmds::Mesh *mesh_target;
        const kmds::Connectivity *c_E2C;
        const kmds::Variable<kmds::TCellID> *varOldCells2cellsSelection;

        Kokkos::View<kmds::TCellID *[InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB +
                                                                          1]> oldCell2subNodes;

        InterfaceNodesPosSmoothVF_createNodesFromEdges_2D(

                const kmds::GrowingView<kmds::TCellID> *selection_,
                kmds::Mesh *mesh_source_,
                kmds::Mesh *mesh_target_,
                const kmds::Connectivity *c_E2C_,
                const kmds::Variable<kmds::TCellID> *varOldCells2cellsSelection_,
                Kokkos::View<kmds::TCellID *[InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB +
                                                                                  1]> oldCell2subNodes_
        )
                : selection(selection_), mesh_source(mesh_source_), mesh_target(mesh_target_), c_E2C(c_E2C_),
                  varOldCells2cellsSelection(varOldCells2cellsSelection_), oldCell2subNodes(oldCell2subNodes_) {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            const int eid = selection->get(i);

            const int nbSub = InterfaceNodesPosSmoothVF_NBSUB - 1;
            if (nbSub != 0) {

                kmds::Edge e = mesh_source->getEdge(eid);
                Kokkos::View<kmds::TCellID *> nids;
                e.nodeIds(nids);

                double xyz_0[3];
                mesh_source->getNodeLocation(nids[0], xyz_0[0], xyz_0[1], xyz_0[2]);
                double xyz_1[3];
                mesh_source->getNodeLocation(nids[1], xyz_1[0], xyz_1[1], xyz_1[2]);

                const kmds::TCellID subID_0 = mesh_target->addNodes(nbSub);

                for (int isub = 1; isub <= nbSub; isub++) {
                    mesh_target->setNodeLocation(subID_0 + isub - 1,
                                                 xyz_0[0] + ((double) isub / (double) InterfaceNodesPosSmoothVF_NBSUB) *
                                                            (xyz_1[0] - xyz_0[0]),
                                                 xyz_0[1] + ((double) isub / (double) InterfaceNodesPosSmoothVF_NBSUB) *
                                                            (xyz_1[1] - xyz_0[1]),
                                                 xyz_0[2] + ((double) isub / (double) InterfaceNodesPosSmoothVF_NBSUB) *
                                                            (xyz_1[2] - xyz_0[2]));
                }

                Kokkos::View<kmds::TCellID *> cids;
                c_E2C->get(eid, cids);

                for (int ic = 0; ic < cids.size(); ic++) {

                    const kmds::TCellID oldcid = cids[ic];

                    if ((*varOldCells2cellsSelection)[oldcid] != kmds::NullID) {

                        const kmds::Face f = mesh_source->getFace(oldcid);

                        const kmds::TCellID subCid = (*varOldCells2cellsSelection)[oldcid];
                        int imin = -1;
                        int imax = -1;
                        int jmin = -1;
                        int jmax = -1;
                        bool direct = false;
                        InterfaceNodesPosSmoothVF_getEdgeSubIndices_2D(f,
                                                                       nids[0],
                                                                       nids[1],
                                                                       &imin,
                                                                       &imax,
                                                                       &jmin,
                                                                       &jmax,
                                                                       &direct);

                        int subID_first = subID_0;
                        int subID_incr = 1;
                        if (!direct) {
                            subID_first = subID_0 + InterfaceNodesPosSmoothVF_NBSUB - 2;
                            subID_incr = -1;
                        }
                        for (int ix = imin; ix <= imax; ix++) {
                            for (int iy = jmin; iy <= jmax; iy++) {
                                oldCell2subNodes(subCid, ix, iy) = subID_first;
                                subID_first += subID_incr;
                            }
                        }
                    }
                }
            }
        }

    };

    /*----------------------------------------------------------------------------*/
    struct InterfaceNodesPosSmoothVF_createNodesFromEdges_3D {
        const kmds::GrowingView<kmds::TCellID> *selection;
        kmds::Mesh *mesh_source;
        kmds::Mesh *mesh_target;
        const kmds::Connectivity *c_E2C;
        const kmds::Variable<kmds::TCellID> *varOldCells2cellsSelection;

        Kokkos::View<kmds::TCellID *[InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB +
                                                                          1]> oldCell2subNodes;

        InterfaceNodesPosSmoothVF_createNodesFromEdges_3D(

                const kmds::GrowingView<kmds::TCellID> *selection_,
                kmds::Mesh *mesh_source_,
                kmds::Mesh *mesh_target_,
                const kmds::Connectivity *c_E2C_,
                const kmds::Variable<kmds::TCellID> *varOldCells2cellsSelection_,
                Kokkos::View<kmds::TCellID *[InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB +
                                                                                  1]> oldCell2subNodes_
        )
                : selection(selection_), mesh_source(mesh_source_), mesh_target(mesh_target_), c_E2C(c_E2C_),
                  varOldCells2cellsSelection(varOldCells2cellsSelection_), oldCell2subNodes(oldCell2subNodes_) {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            const int eid = selection->get(i);

            const int nbSub = InterfaceNodesPosSmoothVF_NBSUB - 1;
            if (nbSub != 0) {

                kmds::Edge e = mesh_source->getEdge(eid);
                Kokkos::View<kmds::TCellID *> nids;
                e.nodeIds(nids);

                double xyz_0[3];
                mesh_source->getNodeLocation(nids[0], xyz_0[0], xyz_0[1], xyz_0[2]);
                double xyz_1[3];
                mesh_source->getNodeLocation(nids[1], xyz_1[0], xyz_1[1], xyz_1[2]);

                const kmds::TCellID subID_0 = mesh_target->addNodes(nbSub);

                for (int isub = 1; isub <= nbSub; isub++) {
                    mesh_target->setNodeLocation(subID_0 + isub - 1,
                                                 xyz_0[0] + ((double) isub / (double) InterfaceNodesPosSmoothVF_NBSUB) *
                                                            (xyz_1[0] - xyz_0[0]),
                                                 xyz_0[1] + ((double) isub / (double) InterfaceNodesPosSmoothVF_NBSUB) *
                                                            (xyz_1[1] - xyz_0[1]),
                                                 xyz_0[2] + ((double) isub / (double) InterfaceNodesPosSmoothVF_NBSUB) *
                                                            (xyz_1[2] - xyz_0[2]));
                }

                Kokkos::View<kmds::TCellID *> cids;
                c_E2C->get(eid, cids);

                for (int ic = 0; ic < cids.size(); ic++) {

                    const kmds::TCellID oldcid = cids[ic];

                    if ((*varOldCells2cellsSelection)[oldcid] != kmds::NullID) {

                        const kmds::Region r = mesh_source->getRegion(oldcid);

                        const kmds::TCellID subCid = (*varOldCells2cellsSelection)[oldcid];
                        int imin = -1;
                        int imax = -1;
                        int jmin = -1;
                        int jmax = -1;
                        int kmin = -1;
                        int kmax = -1;
                        bool direct = false;
                        InterfaceNodesPosSmoothVF_getEdgeSubIndices_3D(r,
                                                                       nids[0],
                                                                       nids[1],
                                                                       &imin,
                                                                       &imax,
                                                                       &jmin,
                                                                       &jmax,
                                                                       &kmin,
                                                                       &kmax,
                                                                       &direct);

                        int subID_first = subID_0;
                        int subID_incr = 1;
                        if (!direct) {
                            subID_first = subID_0 + InterfaceNodesPosSmoothVF_NBSUB - 2;
                            subID_incr = -1;
                        }
                        for (int ix = imin; ix <= imax; ix++) {
                            for (int iy = jmin; iy <= jmax; iy++) {
                                for (int iz = kmin; iz <= kmax; iz++) {
                                    oldCell2subNodes(subCid, ix, iy, iz) = subID_first;
                                    subID_first += subID_incr;
                                }
                            }
                        }
                    }
                }
            }
        }

    };

    /*----------------------------------------------------------------------------*/
    struct InterfaceNodesPosSmoothVF_createNodesFromFaces_3D {
        const kmds::GrowingView<kmds::TCellID> *selection;
        kmds::Mesh *mesh_source;
        kmds::Mesh *mesh_target;
        const kmds::Connectivity *c_F2C;
        const kmds::Variable<kmds::TCellID> *varOldCells2cellsSelection;

        Kokkos::View<kmds::TCellID *[InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB +
                                                                                                               1]> oldCell2subNodes;

        InterfaceNodesPosSmoothVF_createNodesFromFaces_3D(

                const kmds::GrowingView<kmds::TCellID> *selection_,
                kmds::Mesh *mesh_source_,
                kmds::Mesh *mesh_target_,
                const kmds::Connectivity *c_F2C_,
                const kmds::Variable<kmds::TCellID> *varOldCells2cellsSelection_,
                Kokkos::View<kmds::TCellID *[InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB +
                                                                                                                       1]> oldCell2subNodes_
        )
                : selection(selection_), mesh_source(mesh_source_), mesh_target(mesh_target_), c_F2C(c_F2C_),
                  varOldCells2cellsSelection(varOldCells2cellsSelection_), oldCell2subNodes(oldCell2subNodes_) {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            const int fid = selection->get(i);

            const int nbSub = (InterfaceNodesPosSmoothVF_NBSUB - 1);
            if (nbSub != 0) {

                kmds::Face f = mesh_source->getFace(fid);
                Kokkos::View<kmds::TCellID *> nids;
                f.nodeIds(nids);

                double xyz[4][3];
                for (int in = 0; in < nids.size(); in++) {
                    mesh_source->getNodeLocation(nids[in], xyz[in][0], xyz[in][1], xyz[in][2]);
                }
                gmds::math::Point pt0(xyz[0][0], xyz[0][1], xyz[0][2]);
                gmds::math::Point pt1(xyz[1][0], xyz[1][1], xyz[1][2]);
                gmds::math::Point pt2(xyz[2][0], xyz[2][1], xyz[2][2]);
                gmds::math::Point pt3(xyz[3][0], xyz[3][1], xyz[3][2]);

                const kmds::TCellID subID_0 = mesh_target->addNodes(nbSub*nbSub);

                int localindex = 0;
                for (int isub = 1; isub <= nbSub; isub++) {

                    gmds::math::Point pt_tmp0 =
                            pt0 + ((double) isub / (double) InterfaceNodesPosSmoothVF_NBSUB) * (pt1 - pt0);
                    gmds::math::Point pt_tmp1 =
                            pt3 + ((double) isub / (double) InterfaceNodesPosSmoothVF_NBSUB) * (pt2 - pt3);

                    for (int jsub = 1; jsub <= nbSub; jsub++) {

                        gmds::math::Point pt = pt_tmp0 + ((double) jsub / (double) InterfaceNodesPosSmoothVF_NBSUB) *
                                                         (pt_tmp1 - pt_tmp0);

                        mesh_target->setNodeLocation(subID_0 + localindex,
                                                     pt.X(),
                                                     pt.Y(),
                                                     pt.Z()
                        );

                        localindex++;
                    }
                }

                Kokkos::View<kmds::TCellID *> cids;
                c_F2C->get(fid, cids);

                for (int ic = 0; ic < cids.size(); ic++) {

                    const kmds::TCellID oldcid = cids[ic];

                    if ((*varOldCells2cellsSelection)[oldcid] != kmds::NullID) {

                        const kmds::Region r = mesh_source->getRegion(oldcid);

                        const kmds::TCellID subCid = (*varOldCells2cellsSelection)[oldcid];
                        int imin = -1;
                        int imax = -1;
                        int jmin = -1;
                        int jmax = -1;
                        int kmin = -1;
                        int kmax = -1;
                        int faceRel = -1;
                        InterfaceNodesPosSmoothVF_getFaceSubIndices_3D(r,
                                                                       nids[0],
                                                                       nids[1],
                                                                       nids[2],
                                                                       nids[3],
                                                                       &imin,
                                                                       &imax,
                                                                       &jmin,
                                                                       &jmax,
                                                                       &kmin,
                                                                       &kmax,
                                                                       &faceRel);

                        int begin;
                        int istride;
                        int jstride;
                        int kstride = 0;

                        switch(faceRel) {
                            case 0 :
                                begin = 0;
                                istride = nbSub;
                                jstride = 1;
                                break;
                            case 1 :
                                begin = nbSub - 1;
                                istride = -1;
                                jstride = nbSub;
                                break;
                            case 2 :
                                begin = nbSub * nbSub - 1;
                                istride = -nbSub;
                                jstride = -1;
                                break;
                            case 3 :
                                begin = (nbSub - 1) * nbSub;
                                istride = 1;
                                jstride = -nbSub;
                                break;
                            case 4 :
                                begin = 0;
                                istride = 1;
                                jstride = nbSub;
                                break;
                            case 5 :
                                begin = (nbSub - 1) * nbSub;
                                istride = -nbSub;
                                jstride = 1;
                                break;
                            case 6 :
                                begin = nbSub * nbSub - 1;
                                istride = -1;
                                jstride = -nbSub;
                                break;
                            case 7 :
                                begin = (nbSub - 1);
                                istride = nbSub;
                                jstride = -1;
                                break;
                            default:
                                break;
                        }

                        // left or right side
                        if(imin == 0 || imax == InterfaceNodesPosSmoothVF_NBSUB) {
                            kstride = jstride;
                            jstride = istride;
                            istride = 0;
                        }

                        // front or back
                        if(jmin == 0 || jmax == InterfaceNodesPosSmoothVF_NBSUB) {
                            kstride = jstride;
                            jstride = 0;
                        }

                        for (int ix = imin; ix <= imax; ix++) {
                            for (int iy = jmin; iy <= jmax; iy++) {
                                for (int iz = kmin; iz <= kmax; iz++) {
                                    oldCell2subNodes(subCid, ix, iy, iz) = subID_0 + begin + (ix-1)*istride + (iy-1)*jstride + (iz-1)*kstride;
                                }
                            }
                        }
                    }
                }
            }
        }

    };

    /*----------------------------------------------------------------------------*/
    struct InterfaceNodesPosSmoothVF_createNodesFromCells_2D {
        const kmds::GrowingView<kmds::TCellID> *selection;
        kmds::Mesh *mesh_source;
        kmds::Mesh *mesh_target;
        const kmds::Variable<kmds::TCellID> *varOldCells2cellsSelection;

        Kokkos::View<kmds::TCellID *[InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB +
                                                                          1]> oldCell2subNodes;

        InterfaceNodesPosSmoothVF_createNodesFromCells_2D(

                const kmds::GrowingView<kmds::TCellID> *selection_,
                kmds::Mesh *mesh_source_,
                kmds::Mesh *mesh_target_,
                const kmds::Variable<kmds::TCellID> *varOldCells2cellsSelection_,
                Kokkos::View<kmds::TCellID *[InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB +
                                                                                  1]> oldCell2subNodes_
        )
                : selection(selection_), mesh_source(mesh_source_), mesh_target(mesh_target_),
                  varOldCells2cellsSelection(varOldCells2cellsSelection_), oldCell2subNodes(oldCell2subNodes_) {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            const int cid = selection->get(i);

            const int nbSub = InterfaceNodesPosSmoothVF_NBSUB - 1;
            if (nbSub != 0) {

                const kmds::Face f = mesh_source->getFace(cid);

                Kokkos::View<kmds::TCellID *> nids;
                f.nodeIds(nids);

                double xyz[4][3];
                for (int in = 0; in < nids.size(); in++) {
                    mesh_source->getNodeLocation(nids[in], xyz[in][0], xyz[in][1], xyz[in][2]);
                }
                gmds::math::Point pt0(xyz[0][0], xyz[0][1], xyz[0][2]);
                gmds::math::Point pt1(xyz[1][0], xyz[1][1], xyz[1][2]);
                gmds::math::Point pt2(xyz[2][0], xyz[2][1], xyz[2][2]);
                gmds::math::Point pt3(xyz[3][0], xyz[3][1], xyz[3][2]);

                const kmds::TCellID subID_0 = mesh_target->addNodes(
                        (InterfaceNodesPosSmoothVF_NBSUB - 1) * (InterfaceNodesPosSmoothVF_NBSUB - 1));

                const kmds::TCellID subCid = (*varOldCells2cellsSelection)[cid];

                int localindex = 0;
                for (int isub = 1; isub <= nbSub; isub++) {

                    gmds::math::Point pt_tmp0 =
                            pt0 + ((double) isub / (double) InterfaceNodesPosSmoothVF_NBSUB) * (pt1 - pt0);
                    gmds::math::Point pt_tmp1 =
                            pt3 + ((double) isub / (double) InterfaceNodesPosSmoothVF_NBSUB) * (pt2 - pt3);

                    for (int jsub = 1; jsub <= nbSub; jsub++) {

                        gmds::math::Point pt = pt_tmp0 + ((double) jsub / (double) InterfaceNodesPosSmoothVF_NBSUB) *
                                                         (pt_tmp1 - pt_tmp0);

                        mesh_target->setNodeLocation(subID_0 + localindex,
                                                     pt.X(),
                                                     pt.Y(),
                                                     pt.Z()
                        );

                        oldCell2subNodes(subCid, isub, jsub) = subID_0 + localindex;

                        localindex++;
                    }
                }

            }

        }

    };

    /*----------------------------------------------------------------------------*/
    struct InterfaceNodesPosSmoothVF_createNodesFromCells_3D {
        const kmds::GrowingView<kmds::TCellID> *selection;
        kmds::Mesh *mesh_source;
        kmds::Mesh *mesh_target;
        const kmds::Variable<kmds::TCellID> *varOldCells2cellsSelection;

        Kokkos::View<kmds::TCellID *[InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB +
                                                                          1]> oldCell2subNodes;

        InterfaceNodesPosSmoothVF_createNodesFromCells_3D(

                const kmds::GrowingView<kmds::TCellID> *selection_,
                kmds::Mesh *mesh_source_,
                kmds::Mesh *mesh_target_,
                const kmds::Variable<kmds::TCellID> *varOldCells2cellsSelection_,
                Kokkos::View<kmds::TCellID *[InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB +
                                                                                  1]> oldCell2subNodes_
        )
                : selection(selection_), mesh_source(mesh_source_), mesh_target(mesh_target_),
                  varOldCells2cellsSelection(varOldCells2cellsSelection_), oldCell2subNodes(oldCell2subNodes_) {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            const int cid = selection->get(i);

            const int nbSub = InterfaceNodesPosSmoothVF_NBSUB - 1;
            if (nbSub != 0) {

                const kmds::Region r = mesh_source->getRegion(cid);

                Kokkos::View<kmds::TCellID *> nids;
                r.nodeIds(nids);

                double xyz[8][3];
                for (int in = 0; in < nids.size(); in++) {
                    mesh_source->getNodeLocation(nids[in], xyz[in][0], xyz[in][1], xyz[in][2]);
                }
                gmds::math::Point pt0(xyz[0][0], xyz[0][1], xyz[0][2]);
                gmds::math::Point pt1(xyz[1][0], xyz[1][1], xyz[1][2]);
                gmds::math::Point pt2(xyz[2][0], xyz[2][1], xyz[2][2]);
                gmds::math::Point pt3(xyz[3][0], xyz[3][1], xyz[3][2]);
                gmds::math::Point pt4(xyz[4][0], xyz[4][1], xyz[4][2]);
                gmds::math::Point pt5(xyz[5][0], xyz[5][1], xyz[5][2]);
                gmds::math::Point pt6(xyz[6][0], xyz[6][1], xyz[6][2]);
                gmds::math::Point pt7(xyz[7][0], xyz[7][1], xyz[7][2]);

                const kmds::TCellID subID_0 = mesh_target->addNodes(
                        (InterfaceNodesPosSmoothVF_NBSUB - 1) * (InterfaceNodesPosSmoothVF_NBSUB - 1) * (InterfaceNodesPosSmoothVF_NBSUB - 1));

                const kmds::TCellID subCid = (*varOldCells2cellsSelection)[cid];

                int localindex = 0;
                for (int isub = 1; isub <= nbSub; isub++) {

                    gmds::math::Point pt_tmp0 =
                            pt0 + ((double) isub / (double) InterfaceNodesPosSmoothVF_NBSUB) * (pt1 - pt0);
                    gmds::math::Point pt_tmp1 =
                            pt3 + ((double) isub / (double) InterfaceNodesPosSmoothVF_NBSUB) * (pt2 - pt3);

                    gmds::math::Point pt_tmp2 =
                            pt4 + ((double) isub / (double) InterfaceNodesPosSmoothVF_NBSUB) * (pt5 - pt4);
                    gmds::math::Point pt_tmp3 =
                            pt7 + ((double) isub / (double) InterfaceNodesPosSmoothVF_NBSUB) * (pt6 - pt7);

                    for (int jsub = 1; jsub <= nbSub; jsub++) {

                        gmds::math::Point pt_tmp4 =
                                pt_tmp0 + ((double) jsub / (double) InterfaceNodesPosSmoothVF_NBSUB) *
                                          (pt_tmp1 - pt_tmp0);

                        gmds::math::Point pt_tmp5 =
                                pt_tmp2 + ((double) jsub / (double) InterfaceNodesPosSmoothVF_NBSUB) *
                                          (pt_tmp3 - pt_tmp2);

                        for (int ksub = 1; ksub <= nbSub; ksub++) {

                            gmds::math::Point pt =
                                    pt_tmp4 + ((double) ksub / (double) InterfaceNodesPosSmoothVF_NBSUB) *
                                              (pt_tmp5 - pt_tmp4);



                            mesh_target->setNodeLocation(subID_0 + localindex,
                                                         pt.X(),
                                                         pt.Y(),
                                                         pt.Z()
                            );

                            oldCell2subNodes(subCid, isub, jsub, ksub) = subID_0 + localindex;

                            localindex++;
                        }
                    }
                }

            }

        }

    };

    /*----------------------------------------------------------------------------*/
    struct InterfaceNodesPosSmoothVF_createSubCells_2D {
        const kmds::GrowingView<kmds::TCellID> *selection;
        kmds::Mesh *mesh_source;
        kmds::Mesh *mesh_target;
        const kmds::Variable<kmds::TCellID> *varOldCells2cellsSelection;

        Kokkos::View<kmds::TCellID *[InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB +
                                                                          1]> oldCell2subNodes;
        kmds::Variable<kmds::TCellID> *varSubCells2oldCells;
        kmds::Variable<kmds::TCellID> *varOldCells2firstSubCells;

        InterfaceNodesPosSmoothVF_createSubCells_2D(

                const kmds::GrowingView<kmds::TCellID> *selection_,
                kmds::Mesh *mesh_source_,
                kmds::Mesh *mesh_target_,
                const kmds::Variable<kmds::TCellID> *varOldCells2cellsSelection_,
                Kokkos::View<kmds::TCellID *[InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB +
                                                                                  1]> oldCell2subNodes_,
                kmds::Variable<kmds::TCellID> *varSubCells2oldCells_,
                kmds::Variable<kmds::TCellID> *varOldCells2firstSubCells_
        )
                : selection(selection_), mesh_source(mesh_source_), mesh_target(mesh_target_),
                  varOldCells2cellsSelection(varOldCells2cellsSelection_), oldCell2subNodes(oldCell2subNodes_),
                  varSubCells2oldCells(varSubCells2oldCells_), varOldCells2firstSubCells(varOldCells2firstSubCells_) {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            const int cid = selection->get(i);

            const kmds::TCellID subCid = (*varOldCells2cellsSelection)[cid];

            kmds::TCellID subID_0 = mesh_target->addQuads(
                    InterfaceNodesPosSmoothVF_NBSUB2);

            (*varOldCells2firstSubCells)[cid] = subID_0;

            kmds::TCellID localIndex = 0;
            for (int isub = 0; isub < InterfaceNodesPosSmoothVF_NBSUB; isub++) {
                for (int jsub = 0; jsub < InterfaceNodesPosSmoothVF_NBSUB; jsub++) {
                    kmds::TCellID nids[4];
                    nids[0] = oldCell2subNodes(subCid, isub, jsub);
                    nids[1] = oldCell2subNodes(subCid, isub + 1, jsub);
                    nids[2] = oldCell2subNodes(subCid, isub + 1, jsub + 1);
                    nids[3] = oldCell2subNodes(subCid, isub, jsub + 1);
                    kmds::Face f = mesh_target->getFace(subID_0 + localIndex);
                    f.setNodes(nids, 4);

                    (*varSubCells2oldCells)[subID_0 + localIndex] = cid;
                    localIndex++;
                }
            }

        }

    };

    /*----------------------------------------------------------------------------*/
    struct InterfaceNodesPosSmoothVF_createSubCells_3D {
        const kmds::GrowingView<kmds::TCellID> *selection;
        kmds::Mesh *mesh_source;
        kmds::Mesh *mesh_target;
        const kmds::Variable<kmds::TCellID> *varOldCells2cellsSelection;

        Kokkos::View<kmds::TCellID *[InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB +
                                                                          1]> oldCell2subNodes;
        kmds::Variable<kmds::TCellID> *varSubCells2oldCells;
        kmds::Variable<kmds::TCellID> *varOldCells2firstSubCells;

        InterfaceNodesPosSmoothVF_createSubCells_3D(

                const kmds::GrowingView<kmds::TCellID> *selection_,
                kmds::Mesh *mesh_source_,
                kmds::Mesh *mesh_target_,
                const kmds::Variable<kmds::TCellID> *varOldCells2cellsSelection_,
                Kokkos::View<kmds::TCellID *[InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB +
                                                                                  1]> oldCell2subNodes_,
                kmds::Variable<kmds::TCellID> *varSubCells2oldCells_,
                kmds::Variable<kmds::TCellID> *varOldCells2firstSubCells_
        )
                : selection(selection_), mesh_source(mesh_source_), mesh_target(mesh_target_),
                  varOldCells2cellsSelection(varOldCells2cellsSelection_), oldCell2subNodes(oldCell2subNodes_),
                  varSubCells2oldCells(varSubCells2oldCells_), varOldCells2firstSubCells(varOldCells2firstSubCells_) {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            const int cid = selection->get(i);

            const kmds::TCellID subCid = (*varOldCells2cellsSelection)[cid];

            kmds::TCellID subID_0 = mesh_target->addHexahedra(
                    InterfaceNodesPosSmoothVF_NBSUB3);

            (*varOldCells2firstSubCells)[cid] = subID_0;

            kmds::TCellID localIndex = 0;
            for (int isub = 0; isub < InterfaceNodesPosSmoothVF_NBSUB; isub++) {
                for (int jsub = 0; jsub < InterfaceNodesPosSmoothVF_NBSUB; jsub++) {
                    for (int ksub = 0; ksub < InterfaceNodesPosSmoothVF_NBSUB; ksub++) {
                        kmds::TCellID nids[8];
                        nids[0] = oldCell2subNodes(subCid, isub, jsub, ksub);
                        nids[1] = oldCell2subNodes(subCid, isub + 1, jsub, ksub);
                        nids[2] = oldCell2subNodes(subCid, isub + 1, jsub + 1, ksub);
                        nids[3] = oldCell2subNodes(subCid, isub, jsub + 1, ksub);
                        nids[4] = oldCell2subNodes(subCid, isub, jsub, ksub + 1);
                        nids[5] = oldCell2subNodes(subCid, isub + 1, jsub, ksub + 1);
                        nids[6] = oldCell2subNodes(subCid, isub + 1, jsub + 1, ksub + 1);
                        nids[7] = oldCell2subNodes(subCid, isub, jsub + 1, ksub + 1);
                        kmds::Region r = mesh_target->getRegion(subID_0 + localIndex);
                        r.setNodes(nids, 8);

                        (*varSubCells2oldCells)[subID_0 + localIndex] = cid;
                        localIndex++;
                    }
                }
            }

        }

    };

    /*----------------------------------------------------------------------------*/
    struct InterfaceNodesPosSmoothVF_keepCells_xD {
        const kmds::GrowingView<kmds::TCellID> *selection;
        const kmds::Connectivity *c_C2C;
        const kmds::Variable<bool> *varMixedCells;

        kmds::Variable<kmds::TCellID> *varCells2vertices;
        kmds::GrowingView<kmds::TCellID> *selectionCells;

        InterfaceNodesPosSmoothVF_keepCells_xD(

                const kmds::GrowingView<kmds::TCellID> *selection_,
                const kmds::Connectivity *c_C2C_,
                const kmds::Variable<bool> *varMixedCells_,

                kmds::Variable<kmds::TCellID> *varCells2vertices_,
                kmds::GrowingView<kmds::TCellID> *selectionCells_

        )
                : selection(selection_), c_C2C(c_C2C_), varMixedCells(varMixedCells_),
                  varCells2vertices(varCells2vertices_), selectionCells(selectionCells_) {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            const int cid = selection->get(i);

            Kokkos::View<kmds::TCellID *> cids;
            c_C2C->get(cid, cids);

            bool toKeep = false;
            for (int ic = 0; ic < cids.size(); ic++) {
                if ((*varMixedCells)[cids[ic]]) {
                    toKeep = true;
                    break;
                }
            }

            if (toKeep) {

                kmds::TCellID index = selectionCells->push_back(cid);
                (*varCells2vertices)[cid] = index;
            }
        }

    };

    /*----------------------------------------------------------------------------*/
    struct InterfaceNodesPosSmoothVF_buildCellNeighbors_xD {
        const kmds::GrowingView<kmds::TCellID> *selection;

        const kmds::Variable<kmds::TCellID> *varCell2vertices;
        const kmds::Connectivity *c_C2C;

        kmds::Graph *graph;


        InterfaceNodesPosSmoothVF_buildCellNeighbors_xD(
                const kmds::GrowingView<kmds::TCellID> *Selection_,
                const kmds::Variable<kmds::TCellID> *varCell2vertices_,
                const kmds::Connectivity *c_C2C_,
                kmds::Graph *graph_
        )
                : selection(Selection_), varCell2vertices(varCell2vertices_), c_C2C(c_C2C_), graph(graph_) {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int cid = selection->get(i);

            Kokkos::View<kmds::TCellID *> cids;
            c_C2C->get(cid, cids);

            std::set<kmds::TCellID> dist_N2C2N_nodes;

            for (int ic = 0; ic < cids.size(); ic++) {

                // we keep as neighbors the nodes in the selection only
                kmds::TCellID index = (*varCell2vertices)[cids[ic]];
                if (index != kmds::NullID) {
                    dist_N2C2N_nodes.insert(index);
                }
            }

            graph->setNeighbors(i, dist_N2C2N_nodes);
        }
    };

    /*----------------------------------------------------------------------------*/
    struct InterfaceNodesPosSmoothVF_averagerVF_xD {

        const Kokkos::View<bool *> isfixed;
        const int nbMat;
        const double threshold;

        std::vector<gmds::math::VectorDyn> *vf_next;
        std::vector<gmds::math::VectorDyn> *vf_current;
//        std::vector<gmds::math::Vector> *vf_next;
//        std::vector<gmds::math::Vector> *vf_current;


        const kmds::Graph *graph;


        InterfaceNodesPosSmoothVF_averagerVF_xD(
                const Kokkos::View<bool *> isfixed_,

                const int nbMat_,
                const double threshold_,

                std::vector<gmds::math::VectorDyn> *vf_next_,
                std::vector<gmds::math::VectorDyn> *vf_current_,
//                std::vector<gmds::math::Vector> *vf_next_,
//                std::vector<gmds::math::Vector> *vf_current_,

                const kmds::Graph *graph_
        )
                : isfixed(isfixed_),
                  nbMat(nbMat_),
                  threshold(threshold_),
                  vf_next(vf_next_),
                  vf_current(vf_current_),
                  graph(graph_)

                {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int p) const {
            if (!isfixed[p]) {
                gmds::math::VectorDyn v((*vf_current)[p]);
//                gmds::math::Vector v((*vf_current)[p]);

                // do not smooth if one component is higher than threshold
                bool isAboveThres = false;
                for (int imat = 0; imat < nbMat; imat++) {
                    if (v[imat] >= threshold) {
                        isAboveThres = true;
                        break;
                    }
                }
                if (isAboveThres) {
                    (*vf_next)[p] = v;
                    return;
                }

                Kokkos::View<kmds::TCellID *> nids;
                const int nbNeighbors = graph->getNeighbors(p, nids);

                for (int i = 0; i < nbNeighbors; i++) {
                    v = v + (*vf_current)[nids[i]];
                }

                v = 1. / (nbNeighbors + 1) * v;

// no need for normalization, the norm1(average) = 1
                double norm1 = v.normL1();
                (*vf_next)[p] = 1. / norm1 * v;
//                        vf_next[p] = v;
            }
        }

    };

/*----------------------------------------------------------------------------*/
    void
    InterfaceNodesPosSmoothVF_buildSubMesh_2D(kmds::Mesh *AMesh_source,
                                              kmds::Mesh *AMesh_target,
                                              const kmds::Connectivity *c_N2F,
                                              const kmds::Connectivity *c_E2F,
                                              const FracPres *Afp_source,
                                              MaterialAssignment *Ama_target,
                                              kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                                              kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                                              kmds::Variable<bool> *AVarMixedCells_source,
                                              kmds::Variable<bool> *AVarMixedCells_target,
                                              kmds::Variable<double> *AVarSurfvol_source,
                                              kmds::Variable<double> *AVarSurfvol_target)
    {
        kmds::GrowingView<kmds::TCellID> nodesIDs("NODES", AMesh_source->getNbNodes());
        AMesh_source->getNodeIDs(&nodesIDs);

        kmds::GrowingView<kmds::TCellID> selectionCells("CELLS", AMesh_source->getNbFaces());

        Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> cellsSet(AMesh_source->getNbFaces());

        // select nodes and cells to keep
        Kokkos::parallel_for(nodesIDs.getNbElems(),
                             InterfaceNodesPosSmoothVF_identifyCells_xD(
                                     &nodesIDs,
                                     AMesh_source,
                                     c_N2F,
                                     Afp_source,
                                     cellsSet,
                                     &selectionCells)
        );


        std::cout << "selectionCells.getNbElems() " << selectionCells.getNbElems() << std::endl;

        // mark kept cells
        kmds::Variable<bool> *varCellsKept =
                AMesh_source->createVariable<bool>(false, kmds::KMDS_FACE, "cellsKept");

        Kokkos::parallel_for(selectionCells.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID cid = selectionCells.get(i);
                                 (*varCellsKept)[cid] = true;
                             });

        // select nodes to keep
        kmds::GrowingView<kmds::TCellID> selectionNodes("NODES", AMesh_source->getNbNodes());

        Kokkos::parallel_for(AMesh_source->getNbNodes(),
                             InterfaceNodesPosSmoothVF_identifyNodesFromN2C_xD(
                                     &nodesIDs,
                                     c_N2F,
                                     varCellsKept,
                                     &selectionNodes
                             ));

        std::cout << "selectionNodes.getNbElems() " << selectionNodes.getNbElems() << std::endl;

        // select edges to keep
        kmds::GrowingView<kmds::TCellID> edgesIDs("EDGES", AMesh_source->getNbEdges());
        AMesh_source->getEdgeIDs(&edgesIDs);

        kmds::GrowingView<kmds::TCellID> selectionEdges("EDGES", AMesh_source->getNbEdges());

        Kokkos::parallel_for(AMesh_source->getNbEdges(),
                             InterfaceNodesPosSmoothVF_identifyEdgesFromE2C_xD(
                                     &edgesIDs,
                                     c_E2F,
                                     varCellsKept,
                                     &selectionEdges
                             ));

        std::cout << "selectionEdges.getNbElems() " << selectionEdges.getNbElems() << std::endl;

        // build subnodes

        // update target mesh containers capacity
        const int nbOldCellsSelection = selectionCells.getNbElems();
        AMesh_target->updateNodeCapacity(
                nbOldCellsSelection * (InterfaceNodesPosSmoothVF_NBSUB + 1) * (InterfaceNodesPosSmoothVF_NBSUB + 1));
        AMesh_target->updateFaceCapacity(
                nbOldCellsSelection * InterfaceNodesPosSmoothVF_NBSUB2);

        // prepare for oldcell to subNodes connectivity
        kmds::Variable<kmds::TCellID> *varOldCells2cellsSelection =
                AMesh_source->createVariable<kmds::TCellID>(kmds::NullID, kmds::KMDS_FACE, "oldCells2cellsSelection");
        Kokkos::View<kmds::TCellID *[InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB +
                                                                          1]> oldCell2subNodes("oldCell2subNodes",
                                                                                               nbOldCellsSelection);

        Kokkos::parallel_for(nbOldCellsSelection,
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID cid = selectionCells.get(i);
                                 (*varOldCells2cellsSelection)[cid] = i;
                             });

        // create nodes from nodes
        Kokkos::parallel_for(selectionNodes.getNbElems(),
                             InterfaceNodesPosSmoothVF_createNodesFromNodes_2D(
                                     &selectionNodes,
                                     AMesh_source,
                                     AMesh_target,
                                     c_N2F,
                                     varOldCells2cellsSelection,
                                     oldCell2subNodes
                             ));

        // create nodes from edges
        Kokkos::parallel_for(selectionEdges.getNbElems(),
                             InterfaceNodesPosSmoothVF_createNodesFromEdges_2D(
                                     &selectionEdges,
                                     AMesh_source,
                                     AMesh_target,
                                     c_E2F,
                                     varOldCells2cellsSelection,
                                     oldCell2subNodes
                             ));

        // create nodes from cells
        Kokkos::parallel_for(selectionCells.getNbElems(),
                             InterfaceNodesPosSmoothVF_createNodesFromCells_2D(
                                     &selectionCells,
                                     AMesh_source,
                                     AMesh_target,
                                     varOldCells2cellsSelection,
                                     oldCell2subNodes
                             ));

        // create the sub cells
        Kokkos::parallel_for(selectionCells.getNbElems(),
                             InterfaceNodesPosSmoothVF_createSubCells_2D(
                                     &selectionCells,
                                     AMesh_source,
                                     AMesh_target,
                                     varOldCells2cellsSelection,
                                     oldCell2subNodes,
                                     AVarSubCells2oldCells,
                                     AVarOldCells2firstSubCells
                             ));


        // fill data needed for computation
        kmds::GrowingView<kmds::TCellID> cellsIDs_source("CELLS", AMesh_source->getNbFaces());
        AMesh_source->getFaceIDs(&cellsIDs_source);

        // mark source mixed cells
        Kokkos::parallel_for(cellsIDs_source.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID cid = cellsIDs_source.get(i);
                                 if (Afp_source->isMixedCell(cid)) {
                                     (*AVarMixedCells_source)[cid] = true;
                                 }
                             });

        // compute source cells surfvol
        Kokkos::parallel_for(cellsIDs_source.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID cid = cellsIDs_source.get(i);
                                 kmds::Face f = AMesh_source->getFace(cid);
                                 (*AVarSurfvol_source)[cid] = f.surfvol();
                             });

        kmds::GrowingView<kmds::TCellID> cellsIDs_target("CELLS", AMesh_target->getNbFaces());
        AMesh_target->getFaceIDs(&cellsIDs_target);

        // mark target mixed cells
        Kokkos::parallel_for(cellsIDs_target.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID cid = cellsIDs_target.get(i);
                                 if ((*AVarMixedCells_source)[(*AVarSubCells2oldCells)[cid]]) {
                                     (*AVarMixedCells_target)[cid] = true;
                                 }
                             });

        // compute target cells surfvol
        Kokkos::parallel_for(cellsIDs_target.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID cid = cellsIDs_target.get(i);
                                 kmds::Face f = AMesh_target->getFace(cid);
                                 (*AVarSurfvol_target)[cid] = f.surfvol();
                             });


        // assign subcells created from pure source cells

        // initialize material assignment in target mesh
        Ama_target->updateCapacity(AMesh_target->getNbFaces());
        Ama_target->createMaterials(Afp_source->getMaterialList());

        // create a "fake" material that regroups the unassigned pixels
        // and set the pixels to unassigned
        int matUndefID = Ama_target->createMaterial("unassigned");

        Kokkos::parallel_for(cellsIDs_target.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID cid = cellsIDs_target.get(i);
                                 Ama_target->setMaterial(matUndefID, cid);
                             });



        Kokkos::parallel_for(cellsIDs_source.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID cid = cellsIDs_source.get(i);
                                 if (!Afp_source->isMixedCell(cid)) {
                                     kmds::TCellID subID_0 = (*AVarOldCells2firstSubCells)[cid];
                                     if(subID_0 != kmds::NullID) {
                                         for(int isub=0; isub<InterfaceNodesPosSmoothVF_NBSUB2; isub++) {
                                             Ama_target->setMaterial(Afp_source->getMaxMatFracPresIndex(cid), subID_0 + isub);
                                         }
                                     }
                                 }
                             });


        // free temporary resources
        AMesh_source->deleteVariable(kmds::KMDS_FACE, "cellsKept");
        AMesh_source->deleteVariable(kmds::KMDS_FACE, "oldCells2cellsSelection");
    }

/*----------------------------------------------------------------------------*/
    void
    InterfaceNodesPosSmoothVF_buildSubMesh_3D(kmds::Mesh *AMesh_source,
                                              kmds::Mesh *AMesh_target,
                                              const kmds::Connectivity *c_N2C,
                                              const kmds::Connectivity *c_E2C,
                                              const kmds::Connectivity *c_F2C,
                                              const FracPres *Afp_source,
                                              MaterialAssignment *Ama_target,
                                              kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                                              kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                                              kmds::Variable<bool> *AVarMixedCells_source,
                                              kmds::Variable<bool> *AVarMixedCells_target,
                                              kmds::Variable<double> *AVarSurfvol_source,
                                              kmds::Variable<double> *AVarSurfvol_target) {
        kmds::GrowingView<kmds::TCellID> nodesIDs("NODES", AMesh_source->getNbNodes());
        AMesh_source->getNodeIDs(&nodesIDs);

        kmds::GrowingView<kmds::TCellID> selectionCells("CELLS", AMesh_source->getNbFaces());

        Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> cellsSet(AMesh_source->getNbFaces());

        // select nodes and cells to keep
        Kokkos::parallel_for(nodesIDs.getNbElems(),
                             InterfaceNodesPosSmoothVF_identifyCells_xD(
                                     &nodesIDs,
                                     AMesh_source,
                                     c_N2C,
                                     Afp_source,
                                     cellsSet,
                                     &selectionCells)
        );


        std::cout << "selectionCells.getNbElems() " << selectionCells.getNbElems() << std::endl;

        // mark kept cells
        kmds::Variable<bool> *varCellsKept =
                AMesh_source->createVariable<bool>(false, kmds::KMDS_REGION, "cellsKept");

        Kokkos::parallel_for(selectionCells.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID cid = selectionCells.get(i);
                                 (*varCellsKept)[cid] = true;
                             });

        // select nodes to keep
        kmds::GrowingView<kmds::TCellID> selectionNodes("NODES", AMesh_source->getNbNodes());

        Kokkos::parallel_for(AMesh_source->getNbNodes(),
                             InterfaceNodesPosSmoothVF_identifyNodesFromN2C_xD(
                                     &nodesIDs,
                                     c_N2C,
                                     varCellsKept,
                                     &selectionNodes
                             ));

        std::cout << "selectionNodes.getNbElems() " << selectionNodes.getNbElems() << std::endl;

        // select edges to keep
        kmds::GrowingView<kmds::TCellID> edgesIDs("EDGES", AMesh_source->getNbEdges());
        AMesh_source->getEdgeIDs(&edgesIDs);

        kmds::GrowingView<kmds::TCellID> selectionEdges("EDGES", AMesh_source->getNbEdges());

        Kokkos::parallel_for(AMesh_source->getNbEdges(),
                             InterfaceNodesPosSmoothVF_identifyEdgesFromE2C_xD(
                                     &edgesIDs,
                                     c_E2C,
                                     varCellsKept,
                                     &selectionEdges
                             ));

        std::cout << "selectionEdges.getNbElems() " << selectionEdges.getNbElems() << std::endl;

        // select faces to keep
        kmds::GrowingView<kmds::TCellID> facesIDs("FACES", AMesh_source->getNbFaces());
        AMesh_source->getFaceIDs(&facesIDs);

        kmds::GrowingView<kmds::TCellID> selectionFaces("FACES", AMesh_source->getNbFaces());

        Kokkos::parallel_for(AMesh_source->getNbFaces(),
                             InterfaceNodesPosSmoothVF_identifyFacesFromF2C_xD(
                                     &facesIDs,
                                     c_F2C,
                                     varCellsKept,
                                     &selectionFaces
                             ));

        std::cout << "selectionFaces.getNbElems() " << selectionFaces.getNbElems() << std::endl;

        // build subnodes

        // update target mesh containers capacity
        const int nbOldCellsSelection = selectionCells.getNbElems();
        AMesh_target->updateNodeCapacity(
                nbOldCellsSelection * (InterfaceNodesPosSmoothVF_NBSUB + 1) * (InterfaceNodesPosSmoothVF_NBSUB + 1) * (InterfaceNodesPosSmoothVF_NBSUB + 1));
        AMesh_target->updateRegionCapacity(
                nbOldCellsSelection * InterfaceNodesPosSmoothVF_NBSUB3);

        // prepare for oldcell to subNodes connectivity
        kmds::Variable<kmds::TCellID> *varOldCells2cellsSelection =
                AMesh_source->createVariable<kmds::TCellID>(kmds::NullID, kmds::KMDS_REGION, "oldCells2cellsSelection");
        Kokkos::View<kmds::TCellID *[InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB + 1][InterfaceNodesPosSmoothVF_NBSUB + 1]> oldCell2subNodes("oldCell2subNodes",
                                                                                                                                                                      nbOldCellsSelection);

        Kokkos::parallel_for(nbOldCellsSelection,
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID cid = selectionCells.get(i);
                                 (*varOldCells2cellsSelection)[cid] = i;
                             });

        // create nodes from nodes
        Kokkos::parallel_for(selectionNodes.getNbElems(),
                             InterfaceNodesPosSmoothVF_createNodesFromNodes_3D(
                                     &selectionNodes,
                                     AMesh_source,
                                     AMesh_target,
                                     c_N2C,
                                     varOldCells2cellsSelection,
                                     oldCell2subNodes
                             ));

        // create nodes from edges
        Kokkos::parallel_for(selectionEdges.getNbElems(),
                             InterfaceNodesPosSmoothVF_createNodesFromEdges_3D(
                                     &selectionEdges,
                                     AMesh_source,
                                     AMesh_target,
                                     c_E2C,
                                     varOldCells2cellsSelection,
                                     oldCell2subNodes
                             ));

        // create nodes from faces
        Kokkos::parallel_for(selectionFaces.getNbElems(),
                             InterfaceNodesPosSmoothVF_createNodesFromFaces_3D(
                                     &selectionFaces,
                                     AMesh_source,
                                     AMesh_target,
                                     c_F2C,
                                     varOldCells2cellsSelection,
                                     oldCell2subNodes
                             ));

        // create nodes from cells
        Kokkos::parallel_for(selectionCells.getNbElems(),
                             InterfaceNodesPosSmoothVF_createNodesFromCells_3D(
                                     &selectionCells,
                                     AMesh_source,
                                     AMesh_target,
                                     varOldCells2cellsSelection,
                                     oldCell2subNodes
                             ));

//        kmds::TCellID cid_source = 58545;
//        kmds::TCellID cid_select = (*varOldCells2cellsSelection)[cid_source];
//        std::cout<<"cell_58545 cid_select "<<cid_select<<std::endl;
//
//        for (int isub = 0; isub < InterfaceNodesPosSmoothVF_NBSUB+1; isub++) {
//            for (int jsub = 0; jsub < InterfaceNodesPosSmoothVF_NBSUB + 1; jsub++) {
//                for (int ksub = 0; ksub < InterfaceNodesPosSmoothVF_NBSUB + 1; ksub++) {
//
//                    std::cout<<isub<<" "<<jsub<<" "<<ksub<<" "<<oldCell2subNodes(cid_select, isub, jsub, ksub)<<std::endl;
//                }
//            }
//        }
//
//
//exit(-1);

        // create the sub cells
        Kokkos::parallel_for(selectionCells.getNbElems(),
                             InterfaceNodesPosSmoothVF_createSubCells_3D(
                                     &selectionCells,
                                     AMesh_source,
                                     AMesh_target,
                                     varOldCells2cellsSelection,
                                     oldCell2subNodes,
                                     AVarSubCells2oldCells,
                                     AVarOldCells2firstSubCells
                             ));


        // fill data needed for computation
        kmds::GrowingView<kmds::TCellID> cellsIDs_source("CELLS", AMesh_source->getNbRegions());
        AMesh_source->getRegionIDs(&cellsIDs_source);

        // mark source mixed cells
        Kokkos::parallel_for(cellsIDs_source.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID cid = cellsIDs_source.get(i);
                                 if (Afp_source->isMixedCell(cid)) {
                                     (*AVarMixedCells_source)[cid] = true;
                                 }
                             });

        // compute source cells surfvol
        Kokkos::parallel_for(cellsIDs_source.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID cid = cellsIDs_source.get(i);
                                 kmds::Region r = AMesh_source->getRegion(cid);
                                 (*AVarSurfvol_source)[cid] = r.surfvol();
                             });

        kmds::GrowingView<kmds::TCellID> cellsIDs_target("CELLS", AMesh_target->getNbRegions());
        AMesh_target->getRegionIDs(&cellsIDs_target);

        // mark target mixed cells
        Kokkos::parallel_for(cellsIDs_target.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID cid = cellsIDs_target.get(i);
                                 if ((*AVarMixedCells_source)[(*AVarSubCells2oldCells)[cid]]) {
                                     (*AVarMixedCells_target)[cid] = true;
                                 }
                             });

        // compute target cells surfvol
        Kokkos::parallel_for(cellsIDs_target.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID cid = cellsIDs_target.get(i);
                                 kmds::Region r = AMesh_target->getRegion(cid);
                                 (*AVarSurfvol_target)[cid] = r.surfvol();
                             });


        // assign subcells created from pure source cells

        // initialize material assignment in target mesh
        Ama_target->updateCapacity(AMesh_target->getNbRegions());
        Ama_target->createMaterials(Afp_source->getMaterialList());

        Kokkos::parallel_for(cellsIDs_source.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID cid = cellsIDs_source.get(i);
                                 if (!Afp_source->isMixedCell(cid)) {
                                     kmds::TCellID subID_0 = (*AVarOldCells2firstSubCells)[cid];
                                     if(subID_0 != kmds::NullID) {
                                         for(int isub=0; isub<InterfaceNodesPosSmoothVF_NBSUB3; isub++) {
                                             Ama_target->setMaterial(Afp_source->getMaxMatFracPresIndex(cid), subID_0 + isub);
                                         }
                                     }
                                 }
                             });


        AMesh_source->deleteVariable(kmds::KMDS_REGION, "cellsKept");
        AMesh_source->deleteVariable(kmds::KMDS_REGION, "oldCells2cellsSelection");
    }

    /*----------------------------------------------------------------------------*/
    void
    InterfaceNodesPosSmoothVF_selectCells_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                                             const kmds::Connectivity *c_C2C_byN,
                                             const kmds::Variable<bool> *AVarMixedCells_target,
                                             kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                             kmds::GrowingView<kmds::TCellID> *ASelectionCellsIDs_target)
    {
        Kokkos::parallel_for(ACellsIDs_target->getNbElems(),
                             InterfaceNodesPosSmoothVF_keepCells_xD(
                                     ACellsIDs_target,
                                     c_C2C_byN,
                                     AVarMixedCells_target,
                                     AVarCells2vertices,
                                     ASelectionCellsIDs_target
                             ));

        std::cout << "SubMesh ASelectionCellsIDs_target->getNbElems() " << ASelectionCellsIDs_target->getNbElems() << std::endl;
    }

    /*----------------------------------------------------------------------------*/
    void
    InterfaceNodesPosSmoothVF_buildGraph_xD(const kmds::GrowingView<kmds::TCellID> *ASelectionCellsIDs_target,
                                            const kmds::Connectivity *c_C2C_byN,
                                            const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                            kmds::Graph *AGraph)
    {
        Kokkos::parallel_for(ASelectionCellsIDs_target->getNbElems(),
                             InterfaceNodesPosSmoothVF_buildCellNeighbors_xD(
                                     ASelectionCellsIDs_target,
                                     AVarCells2vertices,
                                     c_C2C_byN,
                                     AGraph));
    }

    /*----------------------------------------------------------------------------*/
    void
    InterfaceNodesPosSmoothVF_assignMIP_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
                                           const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                                           const FracPres *Afp_source,
                                           MaterialAssignment *Ama_target,
                                           const kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                                           const kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                                           const kmds::Variable<bool> *AVarMixedCells_source,
                                           const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                           const kmds::Graph *AGraph,
                                           const int ANbSubPixels)
    {
        // =======================================
        const int nbVert = AGraph->getNbVec();
        const int nbMat = Afp_source->getNbMaterials();


        const int nbCells_target = ACellsIDs_target->getNbElems();

        std::map<kmds::TCellID, int> vert2puremat;
        for (int i = 0; i < nbCells_target; i++) {
            kmds::TCellID cid = ACellsIDs_target->get(i);

            kmds::TCellID vert = (*AVarCells2vertices)[cid];
            if (vert != kmds::NullID) {
                kmds::TCellID cid_old = (*AVarSubCells2oldCells)[cid];
                if(!Afp_source->isMixedCell(cid_old)) {
                    for(int imat=0; imat<nbMat; imat++) {
                        if(Afp_source->getFracPres(imat, cid_old) == 1.) {
                            vert2puremat.emplace(vert, imat);
                            break;
                        }
                    }
                }
            }
        }

        int nbMixedCells_source = 0;
        Kokkos::parallel_reduce(ACellsIDs_source->getNbElems(),
                                KOKKOS_LAMBDA(const int i, int& sum) {
                                    kmds::TCellID cid = ACellsIDs_source->get(i);
                                    if ((*AVarMixedCells_source)[cid]) {
                                        sum++;
                                    }
                                },
                                nbMixedCells_source);


        // =======================================
        // GLPK
        // =======================================

        glp_prob *lp;
        int *ia, *ja;
        double *ar;

        lp = glp_create_prob();
        glp_set_obj_dir(lp, GLP_MIN);

        // compute problem size
        int nbRows = 0;
        int nbCols = 0;
        int nnz = 0;

        // pixel color variables
        nbCols += nbVert * nbMat;

        // color unicity per pixel
        nbRows += nbVert;
        nnz += nbVert * nbMat;

        // vf per coarse mixed cell
        nbRows += nbMixedCells_source * nbMat;
        nnz += nbMixedCells_source * nbMat *ANbSubPixels;

        // neighboring similarity
        nbRows += nbVert * nbMat;
        nbCols += 2 * nbVert * nbMat;
        nnz += 3 * nbVert * nbMat + nbMat * AGraph->getNbEdges(0);

        // pure pixels
        nbRows += vert2puremat.size();
        nnz += vert2puremat.size();


        // enforce constraints
        glp_add_rows(lp, nbRows);
        glp_add_cols(lp, nbCols);

        ia = (int *) calloc(1 + nnz, sizeof(int));
        ja = (int *) calloc(1 + nnz, sizeof(int));
        ar = (double *) calloc(1 + nnz, sizeof(double));


        int irow = 1; // start at 1, not 0 !
        int icol = 1; // start at 1, not 0 !
        int innz = 1; // start at 1, not 0 !

        // pixel color variables
        // represented as booleans
        for(int p=0; p<nbVert; p++) {
            for(int imat=0; imat<nbMat; imat++) {
                std::string col_name = "pixelcolor_" + std::to_string(p) + "_" + std::to_string(imat);
                glp_set_col_name(lp, icol, col_name.c_str());
                glp_set_col_kind(lp, icol, GLP_BV);
                // glp_set_col_bnds(lp, icol, GLP_DB, 0, 1); // not necessary because set by GLP_BV

                glp_set_obj_coef(lp, icol, 0.);

                icol++;
            }
        }

        // color unicity per pixel
        for(int p=0; p<nbVert; p++) {
            std::string row_name = "colorunicity_" + std::to_string(p);
            glp_set_row_name(lp, irow, row_name.c_str());
            glp_set_row_bnds(lp, irow, GLP_FX, 1, 0);

            for(int imat=0; imat<nbMat; imat++) {
                ia[innz] = irow;
                ja[innz] = p * nbMat + imat;
                ar[innz] = 1.;

                innz++;
            }

            irow++;
        }

        // vf per coarse mixed cell per material
        const int nbCells_source = ACellsIDs_source->getNbElems();
        for (int i = 0; i < nbCells_source; i++) {
            kmds::TCellID cid = ACellsIDs_source->get(i);
            if ((*AVarMixedCells_source)[cid]) {
                kmds::TCellID subCid_first = (*AVarOldCells2firstSubCells)[cid];

                int sum_vf2int = 0;
                for(int imat=0; imat<nbMat; imat++) {
                    std::string row_name = "vfpreserve_" + std::to_string(cid) + "_" + std::to_string(imat);
                    glp_set_row_name(lp, irow, row_name.c_str());

                    int vf2int = ANbSubPixels * Afp_source->getFracPres(imat, cid);
                    glp_set_row_bnds(lp, irow, GLP_FX, vf2int, 0);

                    for (int isub = 0;
                         isub < ANbSubPixels; isub++) {
                        kmds::TCellID vert = (*AVarCells2vertices)[subCid_first + isub];

                        ia[innz] = irow;
                        ja[innz] = vert * nbMat + imat;
                        ar[innz] = 1.;

                        innz++;
                    }

                    sum_vf2int += vf2int;
                    irow++;
                }

                // check correctness for vf * nbPixels treated as an integer
                if(sum_vf2int != ANbSubPixels) {
                    std::cout<<"sum_vf2int "<<sum_vf2int<<std::endl;
                    throw kmds::KException("InterfaceNodesPosSmoothVF_buildSubGraph_glpk_2D : sum_vf2int not equal to number of pixels.");
                }
            }
        }

        // neighboring similarity
        for(int p=0; p<nbVert; p++) {
            for(int imat=0; imat<nbMat; imat++) {
                std::string col_pos_name = "neighbourcolor_pos_" + std::to_string(p) + "_" + std::to_string(imat);
                glp_set_col_name(lp, icol, col_pos_name.c_str());
                glp_set_col_bnds(lp, icol, GLP_LO, 0, 0);

                glp_set_obj_coef(lp, icol, 1.);

                icol++;

                std::string col_neg_name = "neighbourcolor_neg_" + std::to_string(p) + "_" + std::to_string(imat);
                glp_set_col_name(lp, icol, col_neg_name.c_str());
                glp_set_col_bnds(lp, icol, GLP_LO, 0, 0);

                glp_set_obj_coef(lp, icol, 1.);

                icol++;

                std::string row_name = "neighbourcolor_abs_" + std::to_string(p) + "_" + std::to_string(imat);
                glp_set_row_name(lp, irow, row_name.c_str());
                glp_set_row_bnds(lp, irow, GLP_FX, 0, 0);

                // Cpk
                ia[innz] = irow;
                ja[innz] = p * nbMat + imat;
                ar[innz] = 1.;
                innz++;

                // De
                ia[innz] = irow;
                ja[innz] = nbVert * nbMat + 2 * (p * nbMat + imat) + 0;
                ar[innz] = 1.;
                innz++;

                // de
                ia[innz] = irow;
                ja[innz] = nbVert * nbMat + 2 * (p * nbMat + imat) + 1;
                ar[innz] = - 1.;
                innz++;

                Kokkos::View<kmds::TCellID *> neighbors;
                int nb = AGraph->getNeighbors(p, neighbors);

                for (int in=0; in<neighbors.size(); in++) {

                    ia[innz] = irow;
                    ja[innz] = neighbors[in] * nbMat + imat;
                    ar[innz] = - (1. / neighbors.size());

                    innz++;
                }

                irow++;
            }
        }

        // known color for subpixels of pure cells
        for(int p=0; p<nbVert; p++) {
            if(vert2puremat.find(p) != vert2puremat.end()) {
                int imat = vert2puremat[p];
                std::string row_name = "purepixels_" + std::to_string(p);
                glp_set_row_name(lp, irow, row_name.c_str());
                glp_set_row_bnds(lp, irow, GLP_FX, 1, 0);

                ia[innz] = irow;
                ja[innz] = p * nbMat + imat;
                ar[innz] = 1.;

                innz++;


                irow++;
            }
        }


        std::cout<<"row "<<nbRows<<" irow "<<irow<<std::endl;
        std::cout<<"col "<<nbCols<<" icol "<<icol<<std::endl;
        std::cout<<"nnz "<<nnz<<" innz "<<innz<<std::endl;

        // add one because that is the indices starting index
        for(int j=1; j<=nnz; j++) {
            ja[j]++;
        }

        glp_load_matrix(lp, nnz, ia, ja, ar);

        // glp_term_out( GLP_OFF );
        glp_term_out(GLP_ON);

        glp_iocp glpParams;
        glp_init_iocp(&glpParams);
        glpParams.presolve = GLP_ON;

        glpParams.tm_lim = 600000;
//        glpParams.ps_heur = GLP_ON;

        glp_write_lp(lp, NULL, "cplex.txt");

        glp_write_mps(lp, GLP_MPS_FILE, NULL, "mps.txt");

        int glpErr = 0;
        glpErr = glp_intopt(lp, &glpParams);
        switch (glpErr) {
            case 0:
                std::cout << "GLP OK" << std::endl;
                break;
            default:
                std::cout << "pb solving in GLP." << std::endl;
                break;
        }

        glpErr = glp_mip_status(lp);
        switch (glpErr) {
            case GLP_UNDEF:
                std::cout << " MIP solution is undefined" << std::endl;
                throw kmds::KException("SubMapping::boundaryDiscretization MIP solution is undefined");
                break;
            case GLP_OPT:
                std::cout << " MIP solution is integer optimal" << std::endl;
                break;
            case GLP_FEAS:
                std::cout << " MIP solution is integer feasible" << std::endl;
                break;
            case GLP_NOFEAS:
                std::cout << "problem has no integer feasible solution" << std::endl;
                throw kmds::KException("SubMapping::boundaryDiscretization problem has no integer feasible solution");
                break;
            default:
                throw kmds::KException("SubMapping::boundaryDiscretization glp_intopt unknown return code.");
        }

        glp_print_mip(lp, "pixelAssignment.txt");

        for (int i = 0; i < nbCells_target; i++) {
            kmds::TCellID cid = ACellsIDs_target->get(i);

            kmds::TCellID vert = (*AVarCells2vertices)[cid];
            if (vert != kmds::NullID) {
                int mat = nbMat; // default;
                bool found = false;
                for(int imat=0; imat<nbMat; imat++) {
                    int color = glp_mip_col_val(lp, vert * nbMat + imat + 1);
                    if(color == 1) {
                        if(found) {
                            std::cout<<"QQQQQQQQQQQQQQQQ"<<std::endl;
                        }
                        mat = imat;
                        found = true;
                    }
                }
                Ama_target->setMaterial(mat, cid);
            } else {
//                Ama_target->setMaterial(nbMat, cid);
            }
        }

        free(ia);
        free(ja);
        free(ar);

        glp_delete_prob(lp);

    }

    /*----------------------------------------------------------------------------*/
    void
    InterfaceNodesPosSmoothVF_summary_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
                                         const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                                         kmds::Mesh* AMesh_source,
                                         kmds::Mesh* AMesh_target,
                                         const FracPres *Afp_source,
                                         MaterialAssignment *Ama_target,
                                         const kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                                         const kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                                         const kmds::Variable<bool> *AVarMixedCells_source,
                                         const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                         const kmds::Graph *AGraph,
                                         const kmds::Variable<double> *AVarSurfVol_source,
                                         const kmds::Variable<double> *AVarSurfVol_target,
                                         const int ANbSubPixels)
    {
        // =======================================
        // initialize
        const int dim = (AMesh_source->getNbRegions() > 0) ? 3 : 2;
        const int nbMat = Afp_source->getNbMaterials();

        // initialize discrepancy
        double discrepancy_d = 0.;
        double discrepancy_i = 0.;

        kmds::Variable<double> *vartmp_discrepancy_d;
        kmds::Variable<int> *vartmp_discrepancy_i;

        if(3 == dim) {
            vartmp_discrepancy_d = AMesh_source->createVariable<double>(-1, kmds::KMDS_REGION, "vartmp_discrepancy_d");
            vartmp_discrepancy_i = AMesh_source->createVariable<int>(-1, kmds::KMDS_REGION, "vartmp_discrepancy_i");
        } else {
            vartmp_discrepancy_d = AMesh_source->createVariable<double>(-1, kmds::KMDS_FACE, "vartmp_discrepancy_d");
            vartmp_discrepancy_i = AMesh_source->createVariable<int>(-1, kmds::KMDS_FACE, "vartmp_discrepancy_i");
        }


        // =======================================
        // compute discrepancy
        for(int ic=0; ic<ACellsIDs_source->getNbElems(); ic++) {
            kmds::TCellID cid_source = ACellsIDs_source->get(ic);

            if ((*AVarMixedCells_source)[cid_source]) {

                kmds::TCellID subID_0 = (*AVarOldCells2firstSubCells)[cid_source];

                // WARNING not valid in case of a mesh with several types of cells
                const double vol_d = (*AVarSurfVol_source)[cid_source];
                const int vol_i = ANbSubPixels;

                std::vector<double> volMat_d(nbMat, 0.);
                std::vector<int> volMat_i(nbMat, 0);

                for (int p=0; p<ANbSubPixels; p++) {

                    const kmds::TCellID cid = subID_0 + p;

                    int mat = Ama_target->getMaterial(cid);

                    volMat_d[mat] += (*AVarSurfVol_target)[cid];
                    volMat_i[mat] += 1.;
                }


                double discrepancy_d_loc = 0.;
                double discrepancy_i_loc = 0;
                for (int imat = 0; imat < nbMat; imat++) {
                    discrepancy_d_loc += std::fabs(Afp_source->getFracPres(imat, cid_source) * vol_d - volMat_d[imat]);
                    discrepancy_i_loc += std::fabs(Afp_source->getFracPres(imat, cid_source) * vol_i - volMat_i[imat]);
                }

                (*vartmp_discrepancy_d)[cid_source] = discrepancy_d_loc;
                (*vartmp_discrepancy_i)[cid_source] = discrepancy_i_loc;

                discrepancy_d += discrepancy_d_loc;
                discrepancy_i += discrepancy_i_loc;
            }
        }

        // =======================================
        // initialize edgecut
        double edgecut_d = 0;
        int edgecut_i = 0;

        Kokkos::View<int *> vertAssign ("vertices_assignment", AGraph->getNbVec());
        Kokkos::View<gmds::math::Point *> vertCoords ("vertices_xyz", AGraph->getNbVec());

        for(int ic=0; ic<ACellsIDs_target->getNbElems(); ic++) {
            kmds::TCellID cid_target = ACellsIDs_target->get(ic);

            kmds::TCellID vert = (*AVarCells2vertices)[cid_target];

            if (kmds::NullID != vert) {
                vertAssign[vert] = Ama_target->getMaterial(cid_target);

                if(3 == dim) {
                    vertCoords[vert] = AMesh_target->getRegion(cid_target).midpoint();
                } else {
                    vertCoords[vert] = AMesh_target->getFace(cid_target).midpoint();
                }
            }
        }

        // =======================================
        // compute edgecut
        for(int ic=0; ic<ACellsIDs_target->getNbElems(); ic++) {
            kmds::TCellID cid_target = ACellsIDs_target->get(ic);

            kmds::TCellID vert = (*AVarCells2vertices)[cid_target];

            if(kmds::NullID != vert) {

                int mat = vertAssign[vert];
                gmds::math::Point pt = vertCoords[vert];

                Kokkos::View<kmds::TCellID *> neighbors;
                int nb_0 = AGraph->getNeighbors(vert, neighbors);

                for(int in=0; in<nb_0; in++) {
                    if(mat != vertAssign[neighbors[in]]) {
                        edgecut_i++;

                        edgecut_d += 1. / pt.distance2(vertCoords[neighbors[in]]);
                    }
                }
            }

        }


        std::cout<<"=================="<<std::endl;
        std::cout<<"Pixelated summary "<<std::endl;
        std::cout<<"discrepancy_i "<<discrepancy_i<<std::endl;
        std::cout<<"discrepancy_d "<<discrepancy_d<<std::endl;
        std::cout<<"edgecut_i "<<edgecut_i<<std::endl;
        std::cout<<"edgecut_d "<<edgecut_d<<std::endl;


        // writer
        kmds::VTKWriter<kmds::Mesh> w(*AMesh_source);
        if(3 == dim) {
            w.write("poyop", kmds::R);
        } else {
            w.write("poyop", kmds::F);
        }

        // clean-up
        if(3 == dim) {
            AMesh_source->deleteVariable(kmds::KMDS_REGION, "vartmp_discrepancy_d");
            AMesh_source->deleteVariable(kmds::KMDS_REGION, "vartmp_discrepancy_i");
        } else {
            AMesh_source->deleteVariable(kmds::KMDS_FACE, "vartmp_discrepancy_d");
            AMesh_source->deleteVariable(kmds::KMDS_FACE, "vartmp_discrepancy_i");
        }

    }

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
