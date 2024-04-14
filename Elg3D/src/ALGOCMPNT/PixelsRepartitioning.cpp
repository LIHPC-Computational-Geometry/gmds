/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    PixelsRepartitioning.cpp
 *  \author  legoff
 *  \date    05/20/2020
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/PixelsRepartitioning.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <map>
#include <set>
#include <random>

#include <algorithm>
#include <numeric>

#include <unistd.h>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <Kokkos_Core.hpp>
#include <Kokkos_UnorderedMap.hpp>
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
#include "ELG3D/DATACMPNT/FracPres.h"
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
#include "ELG3D/DATACMPNT/Parameters.h"
#include <ELG3D/ALGOCMPNT/Tools.h>
/*----------------------------------------------------------------------------*/
namespace elg3d {

/*----------------------------------------------------------------------------*/
    void
    PixelsRepartitioning_KL_computeGain(const kmds::TCellID AVert,
                                        const kmds::Graph *AGraph,
                                        const Kokkos::View<gmds::math::Point *> AMidpoints,
                                        const std::vector<int>* AMat,
                                        const int ANbMat,
                                        Kokkos::View<double **> AllGains,
                                        const double AEdgelength_max_tot
                                         )
    {
        for(int imat=0; imat<ANbMat; imat++) {
            AllGains(AVert,imat) = 0;
        }

        Kokkos::View<kmds::TCellID *> nids;
        const int nbNeighbors = AGraph->getNeighbors(AVert, nids);

        double innerGain = 0.;
        double outerGain = 0.;

        const gmds::math::Point pt_v = AMidpoints[AVert];

        double minEdgeLength = HUGE_VALF;

        for(int i_n=0; i_n<nbNeighbors; i_n++) {

            const gmds::math::Point pt_n = AMidpoints[nids[i_n]];
            const double dist = pt_v.distance(pt_n);

            if(minEdgeLength > dist) {
                minEdgeLength = dist;
            }
        }

//            double gain = 1 - pt_v.distance(pt_n) / AEdgelength_max_tot;
//            double gain = 1;

        for(int i_n=0; i_n<nbNeighbors; i_n++) {
            const gmds::math::Point pt_n = AMidpoints[nids[i_n]];
            const double dist = pt_v.distance(pt_n);
            double gain = minEdgeLength / dist;

            // horizontal weight
            const double weight = (1./1.1) * (.1 + std::abs(gmds::math::Vector3d(pt_n - pt_v).getNormalize().dot(gmds::math::Vector3d ({1., 0., 0.}))));
//            gain *= weight;
//            gmds::math::Vector dir_hor(0., 1., 0.);
//            gmds::math::Vector dir(0., 1., 0.);

            AllGains(AVert,(*AMat)[nids[i_n]]) = AllGains(AVert,(*AMat)[nids[i_n]]) + gain;
        }

        for(int imat=0; imat<ANbMat; imat++) {
            AllGains(AVert,imat) = AllGains(AVert,imat) / nbNeighbors;
        }
    }

    /*----------------------------------------------------------------------------*/
    void
    PixelsRepartitioning_KL_computeGain_grad(const kmds::TCellID AVert,
                                             const kmds::Graph *AGraph,
                                             const Kokkos::View<gmds::math::Point *> AMidpoints,
                                             const std::vector<int>* AMat,
                                             const int ANbMat,
                                             Kokkos::View<double **> AllGains,
                                             const double AEdgelength_max_tot,
                                             const kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads_source,
                                             const kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                                             Kokkos::View<kmds::TCellID *> AVert2Cell_target
    )
    {
        const kmds::TCellID cid_target = AVert2Cell_target(AVert);
        const kmds::TCellID cid_source = (*AVarSubCells2oldCells)[cid_target];

        for(int imat=0; imat<ANbMat; imat++) {
            AllGains(AVert,imat) = 0;
        }

        Kokkos::View<kmds::TCellID *> nids;
        const int nbNeighbors = AGraph->getNeighbors(AVert, nids);

        double innerGain = 0.;
        double outerGain = 0.;

        const gmds::math::Point pt_v = AMidpoints[AVert];

        double minEdgeLength = HUGE_VALF;

        for(int i_n=0; i_n<nbNeighbors; i_n++) {

            const gmds::math::Point pt_n = AMidpoints[nids[i_n]];
            const double dist = pt_v.distance(pt_n);

            if(minEdgeLength > dist) {
                minEdgeLength = dist;
            }
        }

//            double gain = 1 - pt_v.distance(pt_n) / AEdgelength_max_tot;
//            double gain = 1;

        for(int i_n=0; i_n<nbNeighbors; i_n++) {
            const gmds::math::Point pt_n = AMidpoints[nids[i_n]];
            const double dist = pt_v.distance(pt_n);
            double gain = minEdgeLength / dist;

            // gradient weight
            const int mat = (*AMat)[nids[i_n]];
            const gmds::math::Vector dir = (*AVarGrads_source)[cid_source].m_grad[mat];

            double weight = 1.;
            if(!dir.isZero()) {
                weight = (1./1.1) * (1.1 - std::abs(gmds::math::Vector3d(pt_n - pt_v).getNormalize().dot(dir)));
            }
            gain *= weight;

            AllGains(AVert,mat) = AllGains(AVert,mat) + gain;
        }

        for(int imat=0; imat<ANbMat; imat++) {
            AllGains(AVert,imat) = AllGains(AVert,imat) / nbNeighbors;
        }
    }

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
                               const FracPres *Afp_target)
    {
        // =======================================
        const int dimension = (AMesh_target->getNbRegions() > 0) ? 3 : 2;
        const kmds::ECellType celltype = (AMesh_target->getNbRegions() > 0) ? kmds::KMDS_REGION : kmds::KMDS_FACE;

        const int nbVert = AGraph->getNbVec();
        const int nbMat = Afp_source->getNbMaterials();

        int nbFixed = 0;
        Kokkos::View<bool *> isfixed("isFixed", nbVert);
        std::vector<int> mat(nbVert, -1);


        const int nbCells_source = ACellsIDs_source->getNbElems();
        Kokkos::View<double **> remainingMatVol("remainingMatVol", nbCells_source, nbMat);

        for (int i = 0; i < nbCells_source; i++) {
            for (int imat = 0; imat < nbMat; imat++) {

                const kmds::TCellID cid = ACellsIDs_source->get(i);
                remainingMatVol(i,imat) = Afp_source->getFracPres(imat, cid) * (*AVarSurfvol_source)[cid];

                // TODO : investigate
                remainingMatVol(i,imat) *= 0.99999999;
            }
        }

        const int nbCells_target = ACellsIDs_target->getNbElems();

        Kokkos::View<kmds::TCellID *> vert2Cell_target("vert2Cell_target", nbVert);

        for (int i = 0; i < nbCells_target; i++) {
            kmds::TCellID cid = ACellsIDs_target->get(i);

            kmds::TCellID vert = (*AVarCells2vertices)[cid];
            if (vert != kmds::NullID) {
                vert2Cell_target(vert) = cid;
            }
        }

        // associate cell midpoints to graph vertices
        Kokkos::View<gmds::math::Point *> midpoints("midPoints", nbVert);

        for (int i = 0; i < nbCells_target; i++) {
            kmds::TCellID cid = ACellsIDs_target->get(i);

            kmds::TCellID vert = (*AVarCells2vertices)[cid];
            if (vert != kmds::NullID) {
                midpoints[vert] = (*AVarMidPoints_target)[cid];
            }
        }

        Kokkos::View<bool *> visitedVertices("visitedVertices", nbVert);
        Kokkos::View<double **> allGains("allGains", nbVert, nbMat);
//        Kokkos::View<double *> innerGains("innerGains", nbVert);
//        Kokkos::View<double *> outerGains("outerGains", nbVert);


        // fix vertices issued from pure source cells
        // initialize mat
        for (int p = 0; p < nbVert; p++) {
            const kmds::TCellID cid_target = vert2Cell_target(p);
            const kmds::TCellID cid_source = (*AVarSubCells2oldCells)[cid_target];

            if(!Afp_source->isMixedCell(cid_source)) {
                isfixed[p] = true;
            }

            mat[p] = Ama_target->getMaterial(cid_target);
        }

        // initial gains computation
        kmds::TCoord edgelength_max_tot = -HUGE_VALF;
        {
            Kokkos::parallel_reduce("Maxreduce",
                                    nbVert,
                                    KOKKOS_LAMBDA(const kmds::TSize p, kmds::TCoord& maxdist) {
                                        Kokkos::View<kmds::TCellID *> nids;
                                        const int nbNeighbors = AGraph->getNeighbors(p, nids);

                                        const gmds::math::Point pt_v = midpoints[p];

                                        for(int i_n=0; i_n<nbNeighbors; i_n++) {

                                            const gmds::math::Point pt_n = midpoints[nids[i_n]];
                                            kmds::TCoord dist = pt_v.distance(pt_n);

                                            maxdist = std::max(maxdist, dist);
                                        }
                                    },
                                    Kokkos::Max<kmds::TCoord> (edgelength_max_tot));

            Kokkos::parallel_for(nbVert,
                                 KOKKOS_LAMBDA(const kmds::TSize p) {

                                     PixelsRepartitioning_KL_computeGain(p,
                                                                         AGraph,
                                                                         midpoints,
                                                                         &mat,
                                                                         nbMat,
                                                                         allGains,
                                                                         edgelength_max_tot);

                                 });
        }

        std::cout<<"edgelength_max_tot "<<edgelength_max_tot<<std::endl;

        ////// DEBUG
        if (Parameters::debug_ouputfile_pixels) {
            std::vector<kmds::Variable<double>* > varGainMat_tmp(nbMat, nullptr);
            for(int imat=0; imat<nbMat; imat++) {
                varGainMat_tmp[imat] = AMesh_target->createVariable<double>(-1, Parameters::celltype,
                                                                            std::string("gainMat_") +
                                                                            std::to_string(imat));

            }

            for(int p=0; p<nbVert; p++) {
                const kmds::TCellID cid_target = vert2Cell_target(p);

                for(int imat=0; imat<nbMat; imat++) {
                    (*(varGainMat_tmp[imat]))[cid_target] = allGains(p, imat);
                }
            }
            if (Parameters::dimension == 2) {
                elg3d::Tools_write_2D(AMesh_target, Afp_target, Ama_target, "gains_KL_iter_0");
            } else {
                elg3d::Tools_write_3D(AMesh_target, Afp_target, Ama_target, "gains_KL_iter_0");
            }

            for(int imat=0; imat<nbMat; imat++) {
                AMesh_target->deleteVariable(Parameters::celltype, std::string("gainMat_") + std::to_string(imat));
            }
        }

        double gain_max = 0.;
        bool hasChanged = true;

        int numOuterIter = 0;

        Kokkos::Timer timer_loop;
        Kokkos::Timer timer_smooth;

        // iterative loop

        while (hasChanged) {

            timer_loop.reset();

            hasChanged = false;

            // fill in the swaps and their cumulative gains
            std::vector<double> swap_sequence_gains;
            std::vector<std::pair<kmds::TCellID, kmds::TCellID> > swap_sequence;

            std::vector<int> mat_tmp = mat;

            Kokkos::View<bool *> isVisited("isVisited", nbVert);
            Kokkos::parallel_for(nbVert,
                                 KOKKOS_LAMBDA(const kmds::TSize p) {
                                     isVisited[p] = isfixed[p];
                                 });


            int iter = 0;
            int nbConsecutiveNegIter = 0;
            double current_cumul_gain = 0;
            double current_cumul_gain_max = -HUGE_VALF;
            double current_cumul_gain_retry = 0;
            while (true) {

                // find the best gain
                double bestGain = -HUGE_VALF;
                kmds::TCellID v0 = kmds::NullID;
                kmds::TCellID v1 = kmds::NullID;

                for (int p = 0; p < nbVert; p++) {

                    if (!isVisited[p]) {

                        const gmds::math::Point pt_v = midpoints[p];

                        const int matp = mat_tmp[p];

                        const kmds::TCellID cid_target = vert2Cell_target(p);
                        const kmds::TCellID cid_source = (*AVarSubCells2oldCells)[cid_target];

                        const kmds::TCellID cid_target_bis = (*AVarOldCells2firstSubCells)[cid_source];
                        for(int i_n=0; i_n<ANbSubPixels; i_n++) {
                            const int p1 = (*AVarCells2vertices)[cid_target_bis + i_n];

                            // no need to check whether p1 exists (ie source cell is mixed); p would be marked as
                            // visited otherwise
                            if(p1 > p) {
                                if (!isVisited[p1]) {
                                    if (mat_tmp[p1] != matp) {

                                        // compute gain of swap between p and p1
                                        double gain = allGains(p,mat_tmp[p1]) + allGains(p1,matp)
                                                      - allGains(p,matp) - allGains(p1,mat_tmp[p1]);

                                        bool areNeighbors = AGraph->isNeighbor(p, p1);
                                        if (areNeighbors) {
                                            const gmds::math::Point pt_n = midpoints[p1];
//                                            double gain_edge = 1 - pt_v.distance(pt_n) / edgelength_max_tot;
                                            double gain_edge = 1;

                                            gain -= 2 * gain_edge;
                                        }

                                        if (bestGain < gain) {
                                            bestGain = gain;
                                            v0 = p;
                                            v1 = p1;
                                        }
                                    }
                                }
                            }
                        }

//                        for (int p1 = p + 1; p1 < nbVert; p1++) {
//
//                            if (!isVisited[p1]) {
//
////                                std::cout<<"matp "<<matp<<" "<<mat_tmp[p1]<<std::endl;
//
//                                if (mat_tmp[p1] != matp) {
//
////                                    std::cout<<"p "<<p<<" p1 "<<p1<<std::endl;
//
//                                    // swap only authorized inside one source cell
//                                    const kmds::TCellID c1_target = vert2Cell_target(p1);
//                                    const kmds::TCellID c1_source = (*AVarSubCells2oldCells)[c1_target];
//                                    if (c1_source == cid_source) {
//
//                                        // compute gain of swap between p and p1
//                                        double gain = allGains(p,mat_tmp[p1]) + allGains(p1,matp)
//                                                      - allGains(p,matp) - allGains(p1,mat_tmp[p1]);
//
//                                        bool areNeighbors = AGraph->isNeighbor(p, p1);
//                                        if (areNeighbors) {
//                                            const gmds::math::Point pt_n = midpoints[p1];
////                                            double gain_edge = 1 - pt_v.distance(pt_n) / edgelength_max_tot;
//                                            double gain_edge = 1;
//
//                                            gain -= 2 * gain_edge;
//                                        }
//
//                                        if (bestGain < gain) {
//                                            bestGain = gain;
//                                            v0 = p;
//                                            v1 = p1;
//                                        }
//                                    }
//                                }
//                            }
//
//                        }  // for (int p1 = p + 1; p1 < nbVert; p1++)
                    }
                }  // for (int p = 0; p < nbVert; p++) {

                // we exit the loop when no swap can be found
                if(bestGain == -HUGE_VALF) {
                    break;
                }

                current_cumul_gain += bestGain;
                current_cumul_gain_retry += bestGain;
                if(current_cumul_gain_max < current_cumul_gain) {
                    current_cumul_gain_max = current_cumul_gain;
                }
                if(bestGain < 0.) {
                    nbConsecutiveNegIter++;
                    current_cumul_gain_retry = bestGain;
                } else {
                    if((bestGain > 0.01 * current_cumul_gain_max) || (current_cumul_gain_retry > 0.1 * current_cumul_gain_max)) {
                        nbConsecutiveNegIter = 0;
                    }
                }
                if(nbConsecutiveNegIter > PixelsRepartitioning_MAXNBCONSECUTIVENEGGAIN) {
                    break;
                }

                std::cout<<"bestgain "<<bestGain<<" "<<v0<<" "<<v1<<" nbConsecutiveNegIter "<<nbConsecutiveNegIter<<" "<<current_cumul_gain_retry<<" "<<current_cumul_gain_max<<std::endl;
                std::swap(mat_tmp[v0], mat_tmp[v1]);


                // update the gains for the not yet visited neighbors of v0 and v1
                // no need to update v0 and v1 as they are fixed at te moment
                {
                    // v0
                    Kokkos::View<kmds::TCellID *> nids;
                    const int nbNeighbors = AGraph->getNeighbors(v0, nids);

                    for(int i_n=0; i_n<nbNeighbors; i_n++) {

                        if(!isVisited[nids[i_n]]) {
                            PixelsRepartitioning_KL_computeGain(nids[i_n],
                                                                AGraph,
                                                                midpoints,
                                                                &mat_tmp,
                                                                nbMat,
                                                                allGains,
                                                                edgelength_max_tot);
                        }
                    }
                }
                {
                    // v1
                    Kokkos::View<kmds::TCellID *> nids;
                    const int nbNeighbors = AGraph->getNeighbors(v1, nids);

                    for(int i_n=0; i_n<nbNeighbors; i_n++) {

                        if(!isVisited[nids[i_n]]) {
                            PixelsRepartitioning_KL_computeGain(nids[i_n],
                                                                AGraph,
                                                                midpoints,
                                                                &mat_tmp,
                                                                nbMat,
                                                                allGains,
                                                                edgelength_max_tot);
                        }
                    }
                }

                // store the best swap in the sequence
                swap_sequence_gains.push_back(bestGain);
                swap_sequence.push_back(std::pair<kmds::TCellID, kmds::TCellID> (v0, v1));

                isVisited[v0] = true;
                isVisited[v1] = true;

                iter++;
            }  // while (true)


            // retrieve the maximum cumulative gain
            std::vector<double> swap_sequence_gains_cumul(swap_sequence_gains.size());
            std::partial_sum(swap_sequence_gains.begin(), swap_sequence_gains.end(), swap_sequence_gains_cumul.begin());
            auto gain_cumul_max_it = std::max_element(swap_sequence_gains_cumul.begin(), swap_sequence_gains_cumul.end());

            // check whether the maximum gain is positive
            if (*gain_cumul_max_it > PixelsRepartitioning_ZEROGAIN) {
                hasChanged = true;

                int nbSwaps = std::distance(swap_sequence_gains_cumul.begin(), gain_cumul_max_it) + 1;

                std::cout << "swap sequence cumul gain " << *gain_cumul_max_it << " nbswaps " << nbSwaps << std::endl;

                for (int iswap = 0; iswap < nbSwaps; iswap++) {
                    std::swap(mat[swap_sequence[iswap].first], mat[swap_sequence[iswap].second]);
                }

                // update the inner and outer gains
                Kokkos::parallel_for(nbVert,
                                     KOKKOS_LAMBDA(const kmds::TSize p) {

                                         PixelsRepartitioning_KL_computeGain(p,
                                                                             AGraph,
                                                                             midpoints,
                                                                             &mat,
                                                                             nbMat,
                                                                             allGains,
                                                                             edgelength_max_tot);
                                     });


                // output for debug
                // update material assignment
                for (int i = 0; i < nbCells_target; i++) {
                    kmds::TCellID cid = ACellsIDs_target->get(i);

                    kmds::TCellID vert = (*AVarCells2vertices)[cid];
                    if (vert != kmds::NullID) {
                        Ama_target->setMaterial(mat[vert], cid);
                    } else {
//                Ama_target->setMaterial(nbMat, cid);
                    }
                }

                ////// DEBUG
                if (Parameters::debug_ouputfile_pixels) {
                    std::vector<kmds::Variable<double> *> varGainMat_tmp(nbMat, nullptr);
                    for (int imat = 0; imat < nbMat; imat++) {
                        varGainMat_tmp[imat] = AMesh_target->createVariable<double>(-1, Parameters::celltype,
                                                                                    std::string("gainMat_") +
                                                                                    std::to_string(imat));
                    }

                    for (int p = 0; p < nbVert; p++) {
                        const kmds::TCellID cid_target = vert2Cell_target(p);

                        for (int imat = 0; imat < nbMat; imat++) {
                            (*(varGainMat_tmp[imat]))[cid_target] = allGains(p, imat);
                        }
                    }
                    std::string outputfilename = std::string("gains_KL_iter_") + std::to_string(numOuterIter + 1);
                    if (Parameters::dimension == 2) {
                        elg3d::Tools_write_2D(AMesh_target, Afp_target, Ama_target, outputfilename);
                    } else {
                        elg3d::Tools_write_3D(AMesh_target, Afp_target, Ama_target, outputfilename);
                    }

                    for (int imat = 0; imat < nbMat; imat++) {
                        AMesh_target->deleteVariable(Parameters::celltype,
                                                     std::string("gainMat_") + std::to_string(imat));
                    }
                }

            }  // if (*gain_cumul_max_it > 0)

            numOuterIter++;
//            break;
        }  // while (hasChanged)


        // update material assignment
        for (int i = 0; i < nbCells_target; i++) {
            kmds::TCellID cid = ACellsIDs_target->get(i);

            kmds::TCellID vert = (*AVarCells2vertices)[cid];
            if (vert != kmds::NullID) {
                Ama_target->setMaterial(mat[vert], cid);
            } else {
//                Ama_target->setMaterial(nbMat, cid);
            }
        }

    }

    /*----------------------------------------------------------------------------*/
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
                                    const kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads_source)
    {
        // =======================================
        const int dimension = (AMesh_target->getNbRegions() > 0) ? 3 : 2;
        const kmds::ECellType celltype = (AMesh_target->getNbRegions() > 0) ? kmds::KMDS_REGION : kmds::KMDS_FACE;

        const int nbVert = AGraph->getNbVec();
        const int nbMat = Afp_source->getNbMaterials();

        int nbFixed = 0;
        Kokkos::View<bool *> isfixed("isFixed", nbVert);
        std::vector<int> mat(nbVert, -1);


        const int nbCells_source = ACellsIDs_source->getNbElems();
        Kokkos::View<double **> remainingMatVol("remainingMatVol", nbCells_source, nbMat);

        for (int i = 0; i < nbCells_source; i++) {
            for (int imat = 0; imat < nbMat; imat++) {

                const kmds::TCellID cid = ACellsIDs_source->get(i);
                remainingMatVol(i,imat) = Afp_source->getFracPres(imat, cid) * (*AVarSurfvol_source)[cid];

                // TODO : investigate
                remainingMatVol(i,imat) *= 0.99999999;
            }
        }

        const int nbCells_target = ACellsIDs_target->getNbElems();

        Kokkos::View<kmds::TCellID *> vert2Cell_target("vert2Cell_target", nbVert);

        for (int i = 0; i < nbCells_target; i++) {
            kmds::TCellID cid = ACellsIDs_target->get(i);

            kmds::TCellID vert = (*AVarCells2vertices)[cid];
            if (vert != kmds::NullID) {
                vert2Cell_target(vert) = cid;
            }
        }

        // associate cell midpoints to graph vertices
        Kokkos::View<gmds::math::Point *> midpoints("midPoints", nbVert);

        for (int i = 0; i < nbCells_target; i++) {
            kmds::TCellID cid = ACellsIDs_target->get(i);

            kmds::TCellID vert = (*AVarCells2vertices)[cid];
            if (vert != kmds::NullID) {
                midpoints[vert] = (*AVarMidPoints_target)[cid];
            }
        }

        Kokkos::View<bool *> visitedVertices("visitedVertices", nbVert);
        Kokkos::View<double **> allGains("allGains", nbVert, nbMat);
//        Kokkos::View<double *> innerGains("innerGains", nbVert);
//        Kokkos::View<double *> outerGains("outerGains", nbVert);


        // fix vertices issued from pure source cells
        // initialize mat
        for (int p = 0; p < nbVert; p++) {
            const kmds::TCellID cid_target = vert2Cell_target(p);
            const kmds::TCellID cid_source = (*AVarSubCells2oldCells)[cid_target];

            if(!Afp_source->isMixedCell(cid_source)) {
                isfixed[p] = true;
            }

            mat[p] = Ama_target->getMaterial(cid_target);
        }

        // initial gains computation
        kmds::TCoord edgelength_max_tot = -HUGE_VALF;
        {
            Kokkos::parallel_reduce("Maxreduce",
                                    nbVert,
                                    KOKKOS_LAMBDA(const kmds::TSize p, kmds::TCoord& maxdist) {
                                        Kokkos::View<kmds::TCellID *> nids;
                                        const int nbNeighbors = AGraph->getNeighbors(p, nids);

                                        const gmds::math::Point pt_v = midpoints[p];

                                        for(int i_n=0; i_n<nbNeighbors; i_n++) {

                                            const gmds::math::Point pt_n = midpoints[nids[i_n]];
                                            kmds::TCoord dist = pt_v.distance(pt_n);

                                            maxdist = std::max(maxdist, dist);
                                        }
                                    },
                                    Kokkos::Max<kmds::TCoord> (edgelength_max_tot));

            Kokkos::parallel_for(nbVert,
                                 KOKKOS_LAMBDA(const kmds::TSize p) {

                                     PixelsRepartitioning_KL_computeGain_grad(p,
                                                                              AGraph,
                                                                              midpoints,
                                                                              &mat,
                                                                              nbMat,
                                                                              allGains,
                                                                              edgelength_max_tot,
                                                                              AVarGrads_source,
                                                                              AVarSubCells2oldCells,
                                                                              vert2Cell_target);

                                 });
        }

        std::cout<<"edgelength_max_tot "<<edgelength_max_tot<<std::endl;

        ////// DEBUG
        if (Parameters::debug_ouputfile_pixels) {
            std::vector<kmds::Variable<double>* > varGainMat_tmp(nbMat, nullptr);
            for(int imat=0; imat<nbMat; imat++) {
                varGainMat_tmp[imat] = AMesh_target->createVariable<double>(-1, Parameters::celltype, std::string("gainMat_") + std::to_string(imat));
            }

            for(int p=0; p<nbVert; p++) {
                const kmds::TCellID cid_target = vert2Cell_target(p);

                for(int imat=0; imat<nbMat; imat++) {
                    (*(varGainMat_tmp[imat]))[cid_target] = allGains(p, imat);
                }
            }
            if(Parameters::dimension == 2) {
                elg3d::Tools_write_2D(AMesh_target, Afp_target, Ama_target, "gains_KL_iter_0");
            } else {
                elg3d::Tools_write_3D(AMesh_target, Afp_target, Ama_target, "gains_KL_iter_0");
            }

            for(int imat=0; imat<nbMat; imat++) {
                AMesh_target->deleteVariable(Parameters::celltype, std::string("gainMat_") + std::to_string(imat));
            }
        }

        double gain_max = 0.;
        bool hasChanged = true;

        int numOuterIter = 0;

        Kokkos::Timer timer_loop;
        Kokkos::Timer timer_smooth;

        // iterative loop

        while (hasChanged) {

            timer_loop.reset();

            hasChanged = false;

            // fill in the swaps and their cumulative gains
            std::vector<double> swap_sequence_gains;
            std::vector<std::pair<kmds::TCellID, kmds::TCellID> > swap_sequence;

            std::vector<int> mat_tmp = mat;

            Kokkos::View<bool *> isVisited("isVisited", nbVert);
            Kokkos::parallel_for(nbVert,
                                 KOKKOS_LAMBDA(const kmds::TSize p) {
                                     isVisited[p] = isfixed[p];
                                 });


            int iter = 0;
            while (true) {

                // find the best gain
                double bestGain = -HUGE_VALF;
                kmds::TCellID v0 = kmds::NullID;
                kmds::TCellID v1 = kmds::NullID;

                for (int p = 0; p < nbVert; p++) {

                    if (!isVisited[p]) {

                        const gmds::math::Point pt_v = midpoints[p];

                        const int matp = mat_tmp[p];

                        const kmds::TCellID cid_target = vert2Cell_target(p);
                        const kmds::TCellID cid_source = (*AVarSubCells2oldCells)[cid_target];

                        const kmds::TCellID cid_target_bis = (*AVarOldCells2firstSubCells)[cid_source];
                        for(int i_n=0; i_n<ANbSubPixels; i_n++) {
                            const int p1 = (*AVarCells2vertices)[cid_target_bis + i_n];

                            // no need to check whether p1 exists (ie source cell is mixed); p would be marked as
                            // visited otherwise
                            if (p1 > p) {
                                if (!isVisited[p1]) {
                                    if (mat_tmp[p1] != matp) {

                                        // compute gain of swap between p and p1
                                        double gain = allGains(p, mat_tmp[p1]) + allGains(p1, matp)
                                                      - allGains(p, matp) - allGains(p1, mat_tmp[p1]);

                                        bool areNeighbors = AGraph->isNeighbor(p, p1);
                                        if (areNeighbors) {
                                            const gmds::math::Point pt_n = midpoints[p1];
//                                            double gain_edge = 1 - pt_v.distance(pt_n) / edgelength_max_tot;
                                            double gain_edge = 1;

                                            gain -= 2 * gain_edge;
                                        }

                                        if (bestGain < gain) {
                                            bestGain = gain;
                                            v0 = p;
                                            v1 = p1;
                                        }
                                    }
                                }
                            }
                        }

//                        for (int p1 = p + 1; p1 < nbVert; p1++) {
//
//                            if (!isVisited[p1]) {
//
////                                std::cout<<"matp "<<matp<<" "<<mat_tmp[p1]<<std::endl;
//
//                                if (mat_tmp[p1] != matp) {
//
////                                    std::cout<<"p "<<p<<" p1 "<<p1<<std::endl;
//
//                                    // swap only authorized inside one source cell
//                                    const kmds::TCellID c1_target = vert2Cell_target(p1);
//                                    const kmds::TCellID c1_source = (*AVarSubCells2oldCells)[c1_target];
//                                    if (c1_source == cid_source) {
//
//                                        // compute gain of swap between p and p1
//                                        double gain = allGains(p,mat_tmp[p1]) + allGains(p1,matp)
//                                                      - allGains(p,matp) - allGains(p1,mat_tmp[p1]);
//
//                                        bool areNeighbors = AGraph->isNeighbor(p, p1);
//                                        if (areNeighbors) {
//                                            const gmds::math::Point pt_n = midpoints[p1];
////                                            double gain_edge = 1 - pt_v.distance(pt_n) / edgelength_max_tot;
//                                            double gain_edge = 1;
//
//                                            gain -= 2 * gain_edge;
//                                        }
//
//                                        if (bestGain < gain) {
//                                            bestGain = gain;
//                                            v0 = p;
//                                            v1 = p1;
//                                        }
//                                    }
//                                }
//                            }
//
//                        }  // for (int p1 = p + 1; p1 < nbVert; p1++)
                    }
                }  // for (int p = 0; p < nbVert; p++) {

                // we exit the loop when no swap can be found
                if(bestGain == -HUGE_VALF) {
                    break;
                }
                std::cout<<"bestgain "<<bestGain<<" "<<v0<<" "<<v1<<std::endl;
                std::swap(mat_tmp[v0], mat_tmp[v1]);


                // update the gains for the not yet visited neighbors of v0 and v1
                // no need to update v0 and v1 as they are fixed at te moment
                {
                    // v0
                    Kokkos::View<kmds::TCellID *> nids;
                    const int nbNeighbors = AGraph->getNeighbors(v0, nids);

                    for(int i_n=0; i_n<nbNeighbors; i_n++) {

                        if(!isVisited[nids[i_n]]) {
                            PixelsRepartitioning_KL_computeGain_grad(nids[i_n],
                                                                     AGraph,
                                                                     midpoints,
                                                                     &mat_tmp,
                                                                     nbMat,
                                                                     allGains,
                                                                     edgelength_max_tot,
                                                                     AVarGrads_source,
                                                                     AVarSubCells2oldCells,
                                                                     vert2Cell_target);
                        }
                    }
                }
                {
                    // v1
                    Kokkos::View<kmds::TCellID *> nids;
                    const int nbNeighbors = AGraph->getNeighbors(v1, nids);

                    for(int i_n=0; i_n<nbNeighbors; i_n++) {

                        if(!isVisited[nids[i_n]]) {
                            PixelsRepartitioning_KL_computeGain_grad(nids[i_n],
                                                                     AGraph,
                                                                     midpoints,
                                                                     &mat_tmp,
                                                                     nbMat,
                                                                     allGains,
                                                                     edgelength_max_tot,
                                                                     AVarGrads_source,
                                                                     AVarSubCells2oldCells,
                                                                     vert2Cell_target);
                        }
                    }
                }

                // store the best swap in the sequence
                swap_sequence_gains.push_back(bestGain);
                swap_sequence.push_back(std::pair<kmds::TCellID, kmds::TCellID> (v0, v1));

                isVisited[v0] = true;
                isVisited[v1] = true;

                iter++;
            }  // while (true)


            // retrieve the maximum cumulative gain
            std::vector<double> swap_sequence_gains_cumul(swap_sequence_gains.size());
            std::partial_sum(swap_sequence_gains.begin(), swap_sequence_gains.end(), swap_sequence_gains_cumul.begin());
            auto gain_cumul_max_it = std::max_element(swap_sequence_gains_cumul.begin(), swap_sequence_gains_cumul.end());

            // check whether the maximum gain is positive
            if (*gain_cumul_max_it > PixelsRepartitioning_ZEROGAIN) {
                hasChanged = true;

                int nbSwaps = std::distance(swap_sequence_gains_cumul.begin(), gain_cumul_max_it) + 1;

                std::cout<<"swap sequence cumul gain "<<*gain_cumul_max_it<<" nbswaps "<<nbSwaps<<std::endl;

                for(int iswap=0; iswap<nbSwaps; iswap++) {
                    std::swap(mat[swap_sequence[iswap].first], mat[swap_sequence[iswap].second]);
                }

                // update the inner and outer gains
                Kokkos::parallel_for(nbVert,
                                     KOKKOS_LAMBDA(const kmds::TSize p) {

                                         PixelsRepartitioning_KL_computeGain_grad(p,
                                                                                  AGraph,
                                                                                  midpoints,
                                                                                  &mat,
                                                                                  nbMat,
                                                                                  allGains,
                                                                                  edgelength_max_tot,
                                                                                  AVarGrads_source,
                                                                                  AVarSubCells2oldCells,
                                                                                  vert2Cell_target);
                                     });

                // output for debug
                // update material assignment
                for (int i = 0; i < nbCells_target; i++) {
                    kmds::TCellID cid = ACellsIDs_target->get(i);

                    kmds::TCellID vert = (*AVarCells2vertices)[cid];
                    if (vert != kmds::NullID) {
                        Ama_target->setMaterial(mat[vert], cid);
                    } else {
//                Ama_target->setMaterial(nbMat, cid);
                    }
                }

                ////// DEBUG
                if (Parameters::debug_ouputfile_pixels) {
                    std::vector<kmds::Variable<double>* > varGainMat_tmp(nbMat, nullptr);
                    for(int imat=0; imat<nbMat; imat++) {
                        varGainMat_tmp[imat] = AMesh_target->createVariable<double>(-1, Parameters::celltype, std::string("gainMat_") + std::to_string(imat));
                    }

                    for(int p=0; p<nbVert; p++) {
                        const kmds::TCellID cid_target = vert2Cell_target(p);

                        for(int imat=0; imat<nbMat; imat++) {
                            (*(varGainMat_tmp[imat]))[cid_target] = allGains(p, imat);
                        }
                    }
                    std::string outputfilename = std::string("gains_KL_iter_") + std::to_string(numOuterIter+1);
                    if(Parameters::dimension == 2) {
                        elg3d::Tools_write_2D(AMesh_target, Afp_target, Ama_target, outputfilename);
                    } else {
                        elg3d::Tools_write_3D(AMesh_target, Afp_target, Ama_target, outputfilename);
                    }

                    for(int imat=0; imat<nbMat; imat++) {
                        AMesh_target->deleteVariable(Parameters::celltype, std::string("gainMat_") + std::to_string(imat));
                    }
                }

            }  // if (*gain_cumul_max_it > 0)

            numOuterIter++;
//            break;
        }  // while (hasChanged)


        // update material assignment
        for (int i = 0; i < nbCells_target; i++) {
            kmds::TCellID cid = ACellsIDs_target->get(i);

            kmds::TCellID vert = (*AVarCells2vertices)[cid];
            if (vert != kmds::NullID) {
                Ama_target->setMaterial(mat[vert], cid);
            } else {
//                Ama_target->setMaterial(nbMat, cid);
            }
        }

    }

    /*----------------------------------------------------------------------------*/
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
                               const FracPres *Afp_target)
    {
        // =======================================
        const int dimension = (AMesh_target->getNbRegions() > 0) ? 3 : 2;
        const kmds::ECellType celltype = (AMesh_target->getNbRegions() > 0) ? kmds::KMDS_REGION : kmds::KMDS_FACE;

        const int nbVert = AGraph->getNbVec();
        const int nbMat = Afp_source->getNbMaterials();

        int nbFixed = 0;
        Kokkos::View<bool *> isfixed("isFixed", nbVert);
        std::vector<int> mat(nbVert, -1);


        const int nbCells_source = ACellsIDs_source->getNbElems();
        Kokkos::View<double **> remainingMatVol("remainingMatVol", nbCells_source, nbMat);

        const int nbCells_target = ACellsIDs_target->getNbElems();

        Kokkos::View<kmds::TCellID *> vert2Cell_target("vert2Cell_target", nbVert);

        for (int i = 0; i < nbCells_target; i++) {
            kmds::TCellID cid = ACellsIDs_target->get(i);

            kmds::TCellID vert = (*AVarCells2vertices)[cid];
            if (vert != kmds::NullID) {
                vert2Cell_target(vert) = cid;
            }
        }

        for (int i = 0; i < nbCells_source; i++) {
            for (int imat = 0; imat < nbMat; imat++) {

                const kmds::TCellID cid = ACellsIDs_source->get(i);
                remainingMatVol(cid,imat) = Afp_source->getFracPres(imat, cid) * (*AVarSurfvol_source)[cid];

//                // TODO : investigate
//                remainingMatVol(cid,imat) *= 0.99999999;
//                remainingMatVol(cid,imat) = 0;
            }
        }
        for (int i = 0; i < nbCells_target; i++) {
            kmds::TCellID cid_target = ACellsIDs_target->get(i);
            const kmds::TCellID cid_source = (*AVarSubCells2oldCells)[cid_target];
            const int mat_loc = Ama_target->getMaterial(cid_target);

            remainingMatVol(cid_source,mat_loc) = remainingMatVol(cid_source,mat_loc) - (*AVarSurfvol_target)[cid_target];
        }


        // associate cell midpoints to graph vertices
        Kokkos::View<gmds::math::Point *> midpoints("midPoints", nbVert);

        for (int i = 0; i < nbCells_target; i++) {
            kmds::TCellID cid = ACellsIDs_target->get(i);

            kmds::TCellID vert = (*AVarCells2vertices)[cid];
            if (vert != kmds::NullID) {
                midpoints[vert] = (*AVarMidPoints_target)[cid];
            }
        }

        Kokkos::View<bool *> visitedVertices("visitedVertices", nbVert);
        Kokkos::View<double **> allGains("allGains", nbVert, nbMat);
//        Kokkos::View<double *> innerGains("innerGains", nbVert);
//        Kokkos::View<double *> outerGains("outerGains", nbVert);


        // fix vertices issued from pure source cells
        // initialize mat
        for (int p = 0; p < nbVert; p++) {
            const kmds::TCellID cid_target = vert2Cell_target(p);
            const kmds::TCellID cid_source = (*AVarSubCells2oldCells)[cid_target];

            if(!Afp_source->isMixedCell(cid_source)) {
                isfixed[p] = true;
            }

            mat[p] = Ama_target->getMaterial(cid_target);
        }

        // initial gains computation
        kmds::TCoord edgelength_max_tot = -HUGE_VALF;
        {
            Kokkos::parallel_reduce("Maxreduce",
                                    nbVert,
                                    KOKKOS_LAMBDA(const kmds::TSize p, kmds::TCoord& maxdist) {
                                        Kokkos::View<kmds::TCellID *> nids;
                                        const int nbNeighbors = AGraph->getNeighbors(p, nids);

                                        const gmds::math::Point pt_v = midpoints[p];

                                        for(int i_n=0; i_n<nbNeighbors; i_n++) {

                                            const gmds::math::Point pt_n = midpoints[nids[i_n]];
                                            kmds::TCoord dist = pt_v.distance(pt_n);

                                            maxdist = std::max(maxdist, dist);
                                        }
                                    },
                                    Kokkos::Max<kmds::TCoord> (edgelength_max_tot));

            Kokkos::parallel_for(nbVert,
                                 KOKKOS_LAMBDA(const kmds::TSize p) {

                                     PixelsRepartitioning_KL_computeGain(p,
                                                                         AGraph,
                                                                         midpoints,
                                                                         &mat,
                                                                         nbMat,
                                                                         allGains,
                                                                         edgelength_max_tot);
                                 });
        }

        std::cout<<"edgelength_max_tot "<<edgelength_max_tot<<std::endl;

        ////// DEBUG
        if (Parameters::debug_ouputfile_pixels) {
            std::vector<kmds::Variable<double>* > varGainMat_tmp(nbMat, nullptr);
            for(int imat=0; imat<nbMat; imat++) {
                varGainMat_tmp[imat] = AMesh_target->createVariable<double>(-1, Parameters::celltype, std::string("gainMat_") + std::to_string(imat));
            }

            for(int p=0; p<nbVert; p++) {
                const kmds::TCellID cid_target = vert2Cell_target(p);

                for(int imat=0; imat<nbMat; imat++) {
                    (*(varGainMat_tmp[imat]))[cid_target] = allGains(p, imat);
                }
            }
            if(Parameters::dimension == 2) {
                elg3d::Tools_write_2D(AMesh_target, Afp_target, Ama_target, "gains_FM_iter_0");
            } else {
                elg3d::Tools_write_3D(AMesh_target, Afp_target, Ama_target, "gains_FM_iter_0");
            }

            for(int imat=0; imat<nbMat; imat++) {
                AMesh_target->deleteVariable(Parameters::celltype, std::string("gainMat_") + std::to_string(imat));
            }
        }

        double gain_max = 0.;
        bool hasChanged = true;

        int numOuterIter = 0;

        Kokkos::Timer timer_loop;
        Kokkos::Timer timer_smooth;

        // iterative loop

        while (hasChanged) {

            timer_loop.reset();

            hasChanged = false;

            // fill in the swaps and their cumulative gains
            std::vector<double> swap_sequence_gains;
            std::vector<std::pair<kmds::TCellID, kmds::TCellID> > swap_sequence;

            std::vector<int> mat_tmp = mat;

            Kokkos::View<double **> remainingMatVol_tmp("remainingMatVol_tmp", nbCells_source, nbMat);
            for (int i = 0; i < nbCells_source; i++) {

                const kmds::TCellID cid = ACellsIDs_source->get(i);

                for(int imat=0; imat<nbMat; imat++) {
                    remainingMatVol_tmp(cid,imat) = remainingMatVol(cid,imat);
                }
            }

            // necessary for fixed pixels (those spawned from pure source cells)
            Kokkos::View<bool *> isVisited("isVisited", nbVert);
            Kokkos::parallel_for(nbVert,
                                 KOKKOS_LAMBDA(const kmds::TSize p) {
                                     isVisited[p] = isfixed[p];
                                 });


            int iter = 0;
            while (true) {

                // find the best gain
                double bestGain = -HUGE_VALF;
                kmds::TCellID v0 = kmds::NullID;
                int mat_change = -1;

                for (int p = 0; p < nbVert; p++) {

                    if (!isVisited[p]) {

                        const gmds::math::Point pt_v = midpoints[p];

                        const int matp = mat_tmp[p];

                        // retrieve the source cell
                        const kmds::TCellID cid_target = vert2Cell_target(p);
                        const kmds::TCellID cid_source = (*AVarSubCells2oldCells)[cid_target];

                        // compute the gain of a change of material for p
                        double gain_loc = -HUGE_VALF;
                        int mat_loc = -1;
                        for(int imat=0; imat<nbMat; imat++) {
                            if(imat != matp) {
                                // only allow change to materials present in the source cell
                                if(Afp_source->getFracPres(imat, cid_source) > 0.) {

                                    // do not stray too far from the frac pres of the source cell
//                                    if((remainingMatVol_tmp(cid_source, imat) > -1)
//                                    && (remainingMatVol_tmp(cid_source, matp) <  1)) {
                                    if((remainingMatVol_tmp(cid_source, imat) - (*AVarSurfvol_target)[cid_target]>  -PixelsRepartitioning_IMBALANCE * (*AVarSurfvol_source)[cid_source])
                                       && (remainingMatVol_tmp(cid_source, matp) + (*AVarSurfvol_target)[cid_target] <  PixelsRepartitioning_IMBALANCE * (*AVarSurfvol_source)[cid_source])) {

                                        double gain = allGains(p, imat) - allGains(p, matp);
                                        if (gain_loc < gain) {
                                            gain_loc = gain;
                                            mat_loc = imat;
                                        }
                                    }
                                }
                            }
                        }


//                        double gain = allGains(p, mat_tmp[p1]) + allGains(p1, matp)
//                                      - allGains(p, matp) - allGains(p1, mat_tmp[p1]);
//
//                        bool areNeighbors = AGraph->isNeighbor(p, p1);
//                        if (areNeighbors) {
//                            const gmds::math::Point pt_n = midpoints[p1];
////                                            double gain_edge = 1 - pt_v.distance(pt_n) / edgelength_max_tot;
//                            double gain_edge = 1;
//
//                            gain -= 2 * gain_edge;
//                        }

                        if (bestGain < gain_loc) {
                            bestGain = gain_loc;
                            v0 = p;
                            mat_change = mat_loc;
                        }
                    }
            }  // for (int p = 0; p < nbVert; p++) {

                // we exit the loop when no swap can be found
                if(bestGain == -HUGE_VALF) {
                    break;
                }

                // update the material frac pres for the source cell of v0
                {
                    const kmds::TCellID cid_target = vert2Cell_target(v0);
                    const kmds::TCellID cid_source = (*AVarSubCells2oldCells)[cid_target];

//                    remainingMatVol_tmp(cid_source, mat_tmp[v0]) = remainingMatVol_tmp(cid_source, mat_tmp[v0]) + 1;
//                    remainingMatVol_tmp(cid_source, mat_change) = remainingMatVol_tmp(cid_source, mat_change) - 1;
                    remainingMatVol_tmp(cid_source, mat_tmp[v0]) = remainingMatVol_tmp(cid_source, mat_tmp[v0]) + (*AVarSurfvol_target)[cid_target];
                    remainingMatVol_tmp(cid_source, mat_change) = remainingMatVol_tmp(cid_source, mat_change) - (*AVarSurfvol_target)[cid_target];

//                    std::cout << "bestgain " << bestGain << " " << v0 << " " << mat_change << " "
//                              << remainingMatVol_tmp(cid_source, mat_change) << std::endl;
                }
                mat_tmp[v0] = mat_change;

                // update the gains for the not yet visited neighbors of v0
                // no need to update v0  as it is fixed at the moment
                {
                    // v0
                    Kokkos::View<kmds::TCellID *> nids;
                    const int nbNeighbors = AGraph->getNeighbors(v0, nids);

                    for(int i_n=0; i_n<nbNeighbors; i_n++) {

                        if(!isVisited[nids[i_n]]) {
                            PixelsRepartitioning_KL_computeGain(nids[i_n],
                                                                AGraph,
                                                                midpoints,
                                                                &mat_tmp,
                                                                nbMat,
                                                                allGains,
                                                                edgelength_max_tot);
                        }
                    }
                }

                // store the best swap in the sequence
                swap_sequence_gains.push_back(bestGain);
                swap_sequence.push_back(std::pair<kmds::TCellID, kmds::TCellID> (v0, mat_change));

                isVisited[v0] = true;

                iter++;
            }  // while (true)


            // retrieve the maximum cumulative gain
            std::vector<double> swap_sequence_gains_cumul(swap_sequence_gains.size());
            std::partial_sum(swap_sequence_gains.begin(), swap_sequence_gains.end(), swap_sequence_gains_cumul.begin());
            auto gain_cumul_max_it = std::max_element(swap_sequence_gains_cumul.begin(), swap_sequence_gains_cumul.end());

            // check whether the maximum gain is positive
            if ((gain_cumul_max_it != swap_sequence_gains_cumul.end()) && (*gain_cumul_max_it > 0)) {
                hasChanged = true;

                int nbSwaps = std::distance(swap_sequence_gains_cumul.begin(), gain_cumul_max_it) + 1;

                std::cout<<"swap sequence cumul gain "<<*gain_cumul_max_it<<" nbswaps "<<nbSwaps<<std::endl;

                for(int iswap=0; iswap<nbSwaps; iswap++) {

                    // update the material frac pres for the source cell of v0
                    {
                        kmds::TCellID v0 = swap_sequence[iswap].first;
                        int mat_change = swap_sequence[iswap].second;
                        const kmds::TCellID cid_target = vert2Cell_target(v0);
                        const kmds::TCellID cid_source = (*AVarSubCells2oldCells)[cid_target];

//                        remainingMatVol(cid_source, mat[v0]) = remainingMatVol(cid_source, mat[v0]) + 1;
//                        remainingMatVol(cid_source, mat_change)  = remainingMatVol(cid_source, mat_change) - 1;
                        remainingMatVol(cid_source, mat[v0]) = remainingMatVol(cid_source, mat[v0]) + (*AVarSurfvol_target)[cid_target];
                        remainingMatVol(cid_source, mat_change)  = remainingMatVol(cid_source, mat_change) - (*AVarSurfvol_target)[cid_target];
                    }

                    mat[swap_sequence[iswap].first] = swap_sequence[iswap].second;
                }

                // update the inner and outer gains
                Kokkos::parallel_for(nbVert,
                                     KOKKOS_LAMBDA(const kmds::TSize p) {

                                         PixelsRepartitioning_KL_computeGain(p,
                                                                             AGraph,
                                                                             midpoints,
                                                                             &mat,
                                                                             nbMat,
                                                                             allGains,
                                                                             edgelength_max_tot);
                                     });

                // output for debug
                // update material assignment
                for (int i = 0; i < nbCells_target; i++) {
                    kmds::TCellID cid = ACellsIDs_target->get(i);

                    kmds::TCellID vert = (*AVarCells2vertices)[cid];
                    if (vert != kmds::NullID) {
                        Ama_target->setMaterial(mat[vert], cid);
                    } else {
//                Ama_target->setMaterial(nbMat, cid);
                    }
                }

                ////// DEBUG
                if (Parameters::debug_ouputfile_pixels) {
                    std::vector<kmds::Variable<double>* > varGainMat_tmp(nbMat, nullptr);
                    for(int imat=0; imat<nbMat; imat++) {
                        varGainMat_tmp[imat] = AMesh_target->createVariable<double>(-1, Parameters::celltype, std::string("gainMat_") + std::to_string(imat));
                    }

                    for(int p=0; p<nbVert; p++) {
                        const kmds::TCellID cid_target = vert2Cell_target(p);

                        for(int imat=0; imat<nbMat; imat++) {
                            (*(varGainMat_tmp[imat]))[cid_target] = allGains(p, imat);
                        }
                    }
                    std::string outputfilename = std::string("gains_FM_iter_") + std::to_string(numOuterIter+1);
                    if(Parameters::dimension == 2) {
                        elg3d::Tools_write_2D(AMesh_target, Afp_target, Ama_target, outputfilename);
                    } else {
                        elg3d::Tools_write_3D(AMesh_target, Afp_target, Ama_target, outputfilename);
                    }

                    for(int imat=0; imat<nbMat; imat++) {
                        AMesh_target->deleteVariable(Parameters::celltype, std::string("gainMat_") + std::to_string(imat));
                    }
                }

            }  // if ((gain_cumul_max_it != swap_sequence_gains_cumul.end()) && (*gain_cumul_max_it > 0))

            numOuterIter++;
//            break;
        }  // while (hasChanged)


        // update material assignment
        for (int i = 0; i < nbCells_target; i++) {
            kmds::TCellID cid = ACellsIDs_target->get(i);

            kmds::TCellID vert = (*AVarCells2vertices)[cid];
            if (vert != kmds::NullID) {
                Ama_target->setMaterial(mat[vert], cid);
            } else {
//                Ama_target->setMaterial(nbMat, cid);
            }
        }

    }

    /*----------------------------------------------------------------------------*/
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
                                    const kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads_source)
    {
        // =======================================
        const int dimension = (AMesh_target->getNbRegions() > 0) ? 3 : 2;
        const kmds::ECellType celltype = (AMesh_target->getNbRegions() > 0) ? kmds::KMDS_REGION : kmds::KMDS_FACE;

        const int nbVert = AGraph->getNbVec();
        const int nbMat = Afp_source->getNbMaterials();

        int nbFixed = 0;
        Kokkos::View<bool *> isfixed("isFixed", nbVert);
        std::vector<int> mat(nbVert, -1);


        const int nbCells_source = ACellsIDs_source->getNbElems();
        Kokkos::View<double **> remainingMatVol("remainingMatVol", nbCells_source, nbMat);

        const int nbCells_target = ACellsIDs_target->getNbElems();

        Kokkos::View<kmds::TCellID *> vert2Cell_target("vert2Cell_target", nbVert);

        for (int i = 0; i < nbCells_target; i++) {
            kmds::TCellID cid = ACellsIDs_target->get(i);

            kmds::TCellID vert = (*AVarCells2vertices)[cid];
            if (vert != kmds::NullID) {
                vert2Cell_target(vert) = cid;
            }
        }

        for (int i = 0; i < nbCells_source; i++) {
            for (int imat = 0; imat < nbMat; imat++) {

                const kmds::TCellID cid = ACellsIDs_source->get(i);
                remainingMatVol(cid,imat) = Afp_source->getFracPres(imat, cid) * (*AVarSurfvol_source)[cid];

//                // TODO : investigate
//                remainingMatVol(cid,imat) *= 0.99999999;
//                remainingMatVol(cid,imat) = 0;
            }
        }
        for (int i = 0; i < nbCells_target; i++) {
            kmds::TCellID cid_target = ACellsIDs_target->get(i);
            const kmds::TCellID cid_source = (*AVarSubCells2oldCells)[cid_target];
            const int mat_loc = Ama_target->getMaterial(cid_target);

            remainingMatVol(cid_source,mat_loc) = remainingMatVol(cid_source,mat_loc) - (*AVarSurfvol_target)[cid_target];
        }


        // associate cell midpoints to graph vertices
        Kokkos::View<gmds::math::Point *> midpoints("midPoints", nbVert);

        for (int i = 0; i < nbCells_target; i++) {
            kmds::TCellID cid = ACellsIDs_target->get(i);

            kmds::TCellID vert = (*AVarCells2vertices)[cid];
            if (vert != kmds::NullID) {
                midpoints[vert] = (*AVarMidPoints_target)[cid];
            }
        }

        Kokkos::View<bool *> visitedVertices("visitedVertices", nbVert);
        Kokkos::View<double **> allGains("allGains", nbVert, nbMat);
//        Kokkos::View<double *> innerGains("innerGains", nbVert);
//        Kokkos::View<double *> outerGains("outerGains", nbVert);


        // fix vertices issued from pure source cells
        // initialize mat
        for (int p = 0; p < nbVert; p++) {
            const kmds::TCellID cid_target = vert2Cell_target(p);
            const kmds::TCellID cid_source = (*AVarSubCells2oldCells)[cid_target];

            if(!Afp_source->isMixedCell(cid_source)) {
                isfixed[p] = true;
            }

            mat[p] = Ama_target->getMaterial(cid_target);
        }

        // initial gains computation
        kmds::TCoord edgelength_max_tot = -HUGE_VALF;
        {
            Kokkos::parallel_reduce("Maxreduce",
                                    nbVert,
                                    KOKKOS_LAMBDA(const kmds::TSize p, kmds::TCoord& maxdist) {
                                        Kokkos::View<kmds::TCellID *> nids;
                                        const int nbNeighbors = AGraph->getNeighbors(p, nids);

                                        const gmds::math::Point pt_v = midpoints[p];

                                        for(int i_n=0; i_n<nbNeighbors; i_n++) {

                                            const gmds::math::Point pt_n = midpoints[nids[i_n]];
                                            kmds::TCoord dist = pt_v.distance(pt_n);

                                            maxdist = std::max(maxdist, dist);
                                        }
                                    },
                                    Kokkos::Max<kmds::TCoord> (edgelength_max_tot));

            Kokkos::parallel_for(nbVert,
                                 KOKKOS_LAMBDA(const kmds::TSize p) {

                                     PixelsRepartitioning_KL_computeGain_grad(p,
                                                                              AGraph,
                                                                              midpoints,
                                                                              &mat,
                                                                              nbMat,
                                                                              allGains,
                                                                              edgelength_max_tot,
                                                                              AVarGrads_source,
                                                                              AVarSubCells2oldCells,
                                                                              vert2Cell_target);
                                 });
        }

        std::cout<<"edgelength_max_tot "<<edgelength_max_tot<<std::endl;

        ////// DEBUG
        if (Parameters::debug_ouputfile_pixels) {
            std::vector<kmds::Variable<double>* > varGainMat_tmp(nbMat, nullptr);
            for(int imat=0; imat<nbMat; imat++) {
                varGainMat_tmp[imat] = AMesh_target->createVariable<double>(-1, Parameters::celltype, std::string("gainMat_") + std::to_string(imat));
            }

            for(int p=0; p<nbVert; p++) {
                const kmds::TCellID cid_target = vert2Cell_target(p);

                for(int imat=0; imat<nbMat; imat++) {
                    (*(varGainMat_tmp[imat]))[cid_target] = allGains(p, imat);
                }
            }
            if(Parameters::dimension == 2) {
                elg3d::Tools_write_2D(AMesh_target, Afp_target, Ama_target, "gains_FM_iter_0");
            } else {
                elg3d::Tools_write_3D(AMesh_target, Afp_target, Ama_target, "gains_FM_iter_0");
            }

            for(int imat=0; imat<nbMat; imat++) {
                AMesh_target->deleteVariable(Parameters::celltype, std::string("gainMat_") + std::to_string(imat));
            }
        }

        double gain_max = 0.;
        bool hasChanged = true;

        int numOuterIter = 0;

        Kokkos::Timer timer_loop;
        Kokkos::Timer timer_smooth;

        // iterative loop

        while (hasChanged) {

            timer_loop.reset();

            hasChanged = false;

            // fill in the swaps and their cumulative gains
            std::vector<double> swap_sequence_gains;
            std::vector<std::pair<kmds::TCellID, kmds::TCellID> > swap_sequence;

            std::vector<int> mat_tmp = mat;

            Kokkos::View<double **> remainingMatVol_tmp("remainingMatVol_tmp", nbCells_source, nbMat);
            for (int i = 0; i < nbCells_source; i++) {

                const kmds::TCellID cid = ACellsIDs_source->get(i);

                for(int imat=0; imat<nbMat; imat++) {
                    remainingMatVol_tmp(cid,imat) = remainingMatVol(cid,imat);
                }
            }

            // necessary for fixed pixels (those spawned from pure source cells)
            Kokkos::View<bool *> isVisited("isVisited", nbVert);
            Kokkos::parallel_for(nbVert,
                                 KOKKOS_LAMBDA(const kmds::TSize p) {
                                     isVisited[p] = isfixed[p];
                                 });


            int iter = 0;
            while (true) {

                // find the best gain
                double bestGain = -HUGE_VALF;
                kmds::TCellID v0 = kmds::NullID;
                int mat_change = -1;

                for (int p = 0; p < nbVert; p++) {

                    if (!isVisited[p]) {

                        const gmds::math::Point pt_v = midpoints[p];

                        const int matp = mat_tmp[p];

                        // retrieve the source cell
                        const kmds::TCellID cid_target = vert2Cell_target(p);
                        const kmds::TCellID cid_source = (*AVarSubCells2oldCells)[cid_target];

                        // compute the gain of a change of material for p
                        double gain_loc = -HUGE_VALF;
                        int mat_loc = -1;
                        for(int imat=0; imat<nbMat; imat++) {
                            if(imat != matp) {
                                // only allow change to materials present in the source cell
                                if(Afp_source->getFracPres(imat, cid_source) > 0.) {

                                    // do not stray too far from the frac pres of the source cell
//                                    if((remainingMatVol_tmp(cid_source, imat) > -1)
//                                    && (remainingMatVol_tmp(cid_source, matp) <  1)) {
                                    if((remainingMatVol_tmp(cid_source, imat) - (*AVarSurfvol_target)[cid_target]>  -PixelsRepartitioning_IMBALANCE * (*AVarSurfvol_source)[cid_source])
                                       && (remainingMatVol_tmp(cid_source, matp) + (*AVarSurfvol_target)[cid_target] <  PixelsRepartitioning_IMBALANCE * (*AVarSurfvol_source)[cid_source])) {

                                        double gain = allGains(p, imat) - allGains(p, matp);
                                        if (gain_loc < gain) {
                                            gain_loc = gain;
                                            mat_loc = imat;
                                        }
                                    }
                                }
                            }
                        }


//                        double gain = allGains(p, mat_tmp[p1]) + allGains(p1, matp)
//                                      - allGains(p, matp) - allGains(p1, mat_tmp[p1]);
//
//                        bool areNeighbors = AGraph->isNeighbor(p, p1);
//                        if (areNeighbors) {
//                            const gmds::math::Point pt_n = midpoints[p1];
////                                            double gain_edge = 1 - pt_v.distance(pt_n) / edgelength_max_tot;
//                            double gain_edge = 1;
//
//                            gain -= 2 * gain_edge;
//                        }

                        if (bestGain < gain_loc) {
                            bestGain = gain_loc;
                            v0 = p;
                            mat_change = mat_loc;
                        }
                    }
                }  // for (int p = 0; p < nbVert; p++) {

                // we exit the loop when no swap can be found
                if(bestGain == -HUGE_VALF) {
                    break;
                }

                // update the material frac pres for the source cell of v0
                {
                    const kmds::TCellID cid_target = vert2Cell_target(v0);
                    const kmds::TCellID cid_source = (*AVarSubCells2oldCells)[cid_target];

//                    remainingMatVol_tmp(cid_source, mat_tmp[v0]) = remainingMatVol_tmp(cid_source, mat_tmp[v0]) + 1;
//                    remainingMatVol_tmp(cid_source, mat_change) = remainingMatVol_tmp(cid_source, mat_change) - 1;
                    remainingMatVol_tmp(cid_source, mat_tmp[v0]) = remainingMatVol_tmp(cid_source, mat_tmp[v0]) + (*AVarSurfvol_target)[cid_target];
                    remainingMatVol_tmp(cid_source, mat_change) = remainingMatVol_tmp(cid_source, mat_change) - (*AVarSurfvol_target)[cid_target];

//                    std::cout << "bestgain " << bestGain << " " << v0 << " " << mat_change << " "
//                              << remainingMatVol_tmp(cid_source, mat_change) << std::endl;
                }
                mat_tmp[v0] = mat_change;

                // update the gains for the not yet visited neighbors of v0
                // no need to update v0  as it is fixed at the moment
                {
                    // v0
                    Kokkos::View<kmds::TCellID *> nids;
                    const int nbNeighbors = AGraph->getNeighbors(v0, nids);

                    for(int i_n=0; i_n<nbNeighbors; i_n++) {

                        if(!isVisited[nids[i_n]]) {
                            PixelsRepartitioning_KL_computeGain_grad(nids[i_n],
                                                                     AGraph,
                                                                     midpoints,
                                                                     &mat_tmp,
                                                                     nbMat,
                                                                     allGains,
                                                                     edgelength_max_tot,
                                                                     AVarGrads_source,
                                                                     AVarSubCells2oldCells,
                                                                     vert2Cell_target);
                        }
                    }
                }

                // store the best swap in the sequence
                swap_sequence_gains.push_back(bestGain);
                swap_sequence.push_back(std::pair<kmds::TCellID, kmds::TCellID> (v0, mat_change));

                isVisited[v0] = true;

                iter++;
            }  // while (true)


            // retrieve the maximum cumulative gain
            std::vector<double> swap_sequence_gains_cumul(swap_sequence_gains.size());
            std::partial_sum(swap_sequence_gains.begin(), swap_sequence_gains.end(), swap_sequence_gains_cumul.begin());
            auto gain_cumul_max_it = std::max_element(swap_sequence_gains_cumul.begin(), swap_sequence_gains_cumul.end());

            // check whether the maximum gain is positive
            if ((gain_cumul_max_it != swap_sequence_gains_cumul.end()) && (*gain_cumul_max_it > 0)) {
                hasChanged = true;

                int nbSwaps = std::distance(swap_sequence_gains_cumul.begin(), gain_cumul_max_it) + 1;

                std::cout<<"swap sequence cumul gain "<<*gain_cumul_max_it<<" nbswaps "<<nbSwaps<<std::endl;

                for(int iswap=0; iswap<nbSwaps; iswap++) {

                    // update the material frac pres for the source cell of v0
                    {
                        kmds::TCellID v0 = swap_sequence[iswap].first;
                        int mat_change = swap_sequence[iswap].second;
                        const kmds::TCellID cid_target = vert2Cell_target(v0);
                        const kmds::TCellID cid_source = (*AVarSubCells2oldCells)[cid_target];

//                        remainingMatVol(cid_source, mat[v0]) = remainingMatVol(cid_source, mat[v0]) + 1;
//                        remainingMatVol(cid_source, mat_change)  = remainingMatVol(cid_source, mat_change) - 1;
                        remainingMatVol(cid_source, mat[v0]) = remainingMatVol(cid_source, mat[v0]) + (*AVarSurfvol_target)[cid_target];
                        remainingMatVol(cid_source, mat_change)  = remainingMatVol(cid_source, mat_change) - (*AVarSurfvol_target)[cid_target];
                    }

                    mat[swap_sequence[iswap].first] = swap_sequence[iswap].second;
                }

                // update the inner and outer gains
                Kokkos::parallel_for(nbVert,
                                     KOKKOS_LAMBDA(const kmds::TSize p) {

                                         PixelsRepartitioning_KL_computeGain_grad(p,
                                                                                  AGraph,
                                                                                  midpoints,
                                                                                  &mat,
                                                                                  nbMat,
                                                                                  allGains,
                                                                                  edgelength_max_tot,
                                                                                  AVarGrads_source,
                                                                                  AVarSubCells2oldCells,
                                                                                  vert2Cell_target);
                                     });

                // output for debug
                // update material assignment
                for (int i = 0; i < nbCells_target; i++) {
                    kmds::TCellID cid = ACellsIDs_target->get(i);

                    kmds::TCellID vert = (*AVarCells2vertices)[cid];
                    if (vert != kmds::NullID) {
                        Ama_target->setMaterial(mat[vert], cid);
                    } else {
//                Ama_target->setMaterial(nbMat, cid);
                    }
                }

                ////// DEBUG
                if (Parameters::debug_ouputfile_pixels) {
                    std::vector<kmds::Variable<double>* > varGainMat_tmp(nbMat, nullptr);
                    for(int imat=0; imat<nbMat; imat++) {
                        varGainMat_tmp[imat] = AMesh_target->createVariable<double>(-1, Parameters::celltype, std::string("gainMat_") + std::to_string(imat));
                    }

                    for(int p=0; p<nbVert; p++) {
                        const kmds::TCellID cid_target = vert2Cell_target(p);

                        for(int imat=0; imat<nbMat; imat++) {
                            (*(varGainMat_tmp[imat]))[cid_target] = allGains(p, imat);
                        }
                    }
                    std::string outputfilename = std::string("gains_FM_iter_") + std::to_string(numOuterIter+1);
                    if(Parameters::dimension == 2) {
                        elg3d::Tools_write_2D(AMesh_target, Afp_target, Ama_target, outputfilename);
                    } else {
                        elg3d::Tools_write_3D(AMesh_target, Afp_target, Ama_target, outputfilename);
                    }

                    for(int imat=0; imat<nbMat; imat++) {
                        AMesh_target->deleteVariable(Parameters::celltype, std::string("gainMat_") + std::to_string(imat));
                    }
                }

            }  // if ((gain_cumul_max_it != swap_sequence_gains_cumul.end()) && (*gain_cumul_max_it > 0))

            numOuterIter++;
//            break;
        }  // while (hasChanged)


        // update material assignment
        for (int i = 0; i < nbCells_target; i++) {
            kmds::TCellID cid = ACellsIDs_target->get(i);

            kmds::TCellID vert = (*AVarCells2vertices)[cid];
            if (vert != kmds::NullID) {
                Ama_target->setMaterial(mat[vert], cid);
            } else {
//                Ama_target->setMaterial(nbMat, cid);
            }
        }

    }

    /*----------------------------------------------------------------------------*/
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
