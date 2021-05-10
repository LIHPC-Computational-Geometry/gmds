/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    MoveToNewPos.h
 *  \author  legoff
 *  \date    05/18/2018
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_MOVETONEWPOS_H_
#define ELG3D_MOVETONEWPOS_H_
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <KM/DS/Mesh.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/DATACMPNT/FracPres.h"
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/

    const int MoveToNewPos_NBMAXCOLORS = 20;


    /*------------------------------------------------------------------------*/
    /** \brief  Move the nodes from the selection to their provided new position.
     *          Additionally projects unto the geometric association of the node, if any,
     *          and does not move if the new position provided is incorrect (ie NaN).
     *
     *  \param[in] ASelectionNodes the selection of the nodes that will be moved
     *  \param[in,out]  AMesh the mesh
     *  \param[in] AVarNodeDestination the new prospective positions
     *  \param[in] AVarNodeGeomAssociation the geometric association of the nodes, for those which have one
     *  \param[out] AVarNodeNbMoves the number of effective displacements
     */
    void moveToNewPos_basicMove(const kmds::GrowingView<kmds::TCellID>* ASelectionInterfaceNodes,
                                kmds::Mesh* AMesh,
                                const kmds::Variable<gmds::math::Point>* AVarNewPos,
                                const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation);


    /*------------------------------------------------------------------------*/
    /** \brief  Move the nodes from the selection to their provided new position, by increments.
     *          Additionally projects unto the geometric association of the node along the way, if any.
     *
     *  \param[in] ANbIter the displacement will be done by this number of increments
     *  \param[in] AQualThreshold a node will stay in place if the minimum scaled jacobian of
     *             the adjacent cells falls below this threshold
     *  \param[in] ASelectionNodes the selection of the nodes that will be moved
     *  \param[in,out]  AMesh the mesh
     *  \param[in] AVarNodeDestination the new prospective positions
     *  \param[in] AVarNodeGeomAssociation the geometric association of the nodes, for those which have one
     *  \param[out] AVarNodeNbMoves the number of effective displacements for each node (for debug purposes)
     *
     */
    void moveToNewPos_noBadMove_2D(const int ANbIter,
                                   const kmds::TCoord qualThreshold,
                                   const kmds::GrowingView<kmds::TCellID>* ASelectionNodes,
                                   kmds::Mesh* AMesh,
                                   const kmds::Connectivity* c_N2F,
                                   const kmds::Variable<gmds::math::Point>* AVarNodeDestination,
                                   const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                                   kmds::Variable<int>* AVarNodeNbMoves);

    void moveToNewPos_noBadMove_3D(const int ANbIter,
                                   const kmds::TCoord qualThreshold,
                                   const kmds::GrowingView<kmds::TCellID>* ASelectionNodes,
                                   kmds::Mesh* AMesh,
                                   const kmds::Connectivity* c_N2R,
                                   const kmds::Variable<gmds::math::Point>* AVarNodeDestination,
                                   const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                                   kmds::Variable<int>* AVarNodeNbMoves);


/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_MOVETONEWPOS_H_ */
