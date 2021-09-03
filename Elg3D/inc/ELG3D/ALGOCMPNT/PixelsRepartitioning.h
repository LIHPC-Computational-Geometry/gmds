/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    PixelsRepartitioning.h
 *  \author  legoff
 *  \date    10/29/2018
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_PIXELSREPARTITIONING_H_
#define ELG3D_PIXELSREPARTITIONING_H_
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
#include <KM/Utils/Graph.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/MaterialGradientComputation.h"
#include "ELG3D/DATACMPNT/FracPres.h"
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/

    const int PixelsRepartitioning_NBSUB = 3;

    const int PixelsRepartitioning_NBSUB2 = PixelsRepartitioning_NBSUB * PixelsRepartitioning_NBSUB;

    const int PixelsRepartitioning_NBSUB3 = PixelsRepartitioning_NBSUB * PixelsRepartitioning_NBSUB * PixelsRepartitioning_NBSUB;

    const double PixelsRepartitioning_IMBALANCE = 0.011;
//    const double PixelsRepartitioning_IMBALANCE = 0.0011;

    const double PixelsRepartitioning_MAXNBCONSECUTIVENEGGAIN = 50;

    const double PixelsRepartitioning_ZEROGAIN = 1e-16;

/*----------------------------------------------------------------------------*/
    void
    PixelsRepartitioning_KL_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
                               const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                               const FracPres *Afp_source,
                               MaterialAssignment *Ama_target,
                               const kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                               const kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                               const kmds::Variable<bool> *AVarMixedCells_source,
                               const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                               const kmds::Variable<double> *AVarSurfvol_source,
                               const kmds::Variable<double> *AVarSurfvol_target,
                               const kmds::Graph *AGraph,
                               const kmds::Variable<gmds::math::Point> *AVarMidPoints_target,
                               const int ANbSubPixels,
                               kmds::Mesh* AMesh_target,
                               const FracPres *Afp_target);

    void
    PixelsRepartitioning_KL_grad_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
                                    const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                                    const FracPres *Afp_source,
                                    MaterialAssignment *Ama_target,
                                    const kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                                    const kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                                    const kmds::Variable<bool> *AVarMixedCells_source,
                                    const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                    const kmds::Variable<double> *AVarSurfvol_source,
                                    const kmds::Variable<double> *AVarSurfvol_target,
                                    const kmds::Graph *AGraph,
                                    const kmds::Variable<gmds::math::Point> *AVarMidPoints_target,
                                    const int ANbSubPixels,
                                    kmds::Mesh* AMesh_target,
                                    const FracPres *Afp_target,
                                    const kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads_source);

    void
    PixelsRepartitioning_FM_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
                               const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                               const FracPres *Afp_source,
                               MaterialAssignment *Ama_target,
                               const kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                               const kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                               const kmds::Variable<bool> *AVarMixedCells_source,
                               const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                               const kmds::Variable<double> *AVarSurfvol_source,
                               const kmds::Variable<double> *AVarSurfvol_target,
                               const kmds::Graph *AGraph,
                               const kmds::Variable<gmds::math::Point> *AVarMidPoints_target,
                               const int ANbSubPixels,
                               kmds::Mesh* AMesh_target,
                               const FracPres *Afp_target);

    void
    PixelsRepartitioning_FM_grad_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
                                    const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                                    const FracPres *Afp_source,
                                    MaterialAssignment *Ama_target,
                                    const kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                                    const kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                                    const kmds::Variable<bool> *AVarMixedCells_source,
                                    const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                    const kmds::Variable<double> *AVarSurfvol_source,
                                    const kmds::Variable<double> *AVarSurfvol_target,
                                    const kmds::Graph *AGraph,
                                    const kmds::Variable<gmds::math::Point> *AVarMidPoints_target,
                                    const int ANbSubPixels,
                                    kmds::Mesh* AMesh_target,
                                    const FracPres *Afp_target,
                                    const kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads_source);

    void
    PixelsRepartitioning_summary_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
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
                                    const int ANbSubPixels);

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_PIXELSREPARTITIONING_H_ */
