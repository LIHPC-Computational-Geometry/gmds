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
