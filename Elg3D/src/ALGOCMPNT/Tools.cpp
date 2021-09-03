/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    Tools.cpp
 *  \author  legoff
 *  \date    05/03/2018
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/Tools.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <map>
#include <string>
#include <vector>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <Kokkos_Core.hpp>

#include <gts.h>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/r2d.h"
#include "ELG3D/ALGOCMPNT/r3d.h"

#include <gmds/ig/Mesh.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/IGMeshIOService.h>

#include <KM/IO/VTKWriter.h>

//#include <GMDS/Algo/GETMe.h>

#include <gmds/math/Point.h>
#include <gmds/math/Triangle.h>
#include <gmds/math/Quadrilateral.h>
#include <gmds/math/Tetrahedron.h>
#include <gmds/math/Hexahedron.h>
#include <gmds/math/Numerics.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
    void
    Tools_computeMidPointVar_2D(const kmds::GrowingView<kmds::TCellID>* ASelectionCells,
                                kmds::Mesh* AMesh,
                                kmds::Variable<gmds::math::Point>* AVarMidpoints)
    {
        Kokkos::parallel_for(ASelectionCells->getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 (*AVarMidpoints)[ASelectionCells->get(i)] = AMesh->getFace(ASelectionCells->get(i)).midpoint();
                             });
    }
/*----------------------------------------------------------------------------*/
    void
    Tools_computeMidPointVar_3D(const kmds::GrowingView<kmds::TCellID>* ASelectionCells,
                                kmds::Mesh* AMesh,
                                kmds::Variable<gmds::math::Point>* AVarMidpoints)
    {
        Kokkos::parallel_for(ASelectionCells->getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 (*AVarMidpoints)[ASelectionCells->get(i)] = AMesh->getRegion(ASelectionCells->get(i)).midpoint();
                             });
    }
/*----------------------------------------------------------------------------*/
    void
    Tools_storeNodePos_xD(kmds::Mesh* AMesh,
                          kmds::Variable<gmds::math::Point>* AVarNodePos)
    {
        kmds::GrowingView <kmds::TCellID> nodesIDs("NODESIDS", AMesh->getNbNodes());
        AMesh->getNodeIDs(&nodesIDs);

        Kokkos::parallel_for(nodesIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID nid = nodesIDs.get(i);
                                 (*AVarNodePos)[nid] = AMesh->getNodeLocation(nid);
                             });
    }
/*----------------------------------------------------------------------------*/
    void
    Tools_loadNodePos_xD(kmds::Mesh* AMesh,
                         const kmds::Variable<gmds::math::Point>* AVarNodePos)
    {
        kmds::GrowingView <kmds::TCellID> nodesIDs("NODESIDS", AMesh->getNbNodes());
        AMesh->getNodeIDs(&nodesIDs);

        Kokkos::parallel_for(nodesIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i) {
                                 kmds::TCellID nid = nodesIDs.get(i);
                                 AMesh->setNodeLocation(nid, (*AVarNodePos)[nid]);
                             });
    }

    /*----------------------------------------------------------------------------*/
    kmds::TCoord
    Tools_computeScaledJacobian_2D(const kmds::Mesh* AMesh)
    {
        kmds::GrowingView <kmds::TCellID> cellsIDs("CELLSIDS", AMesh->getNbFaces());
        AMesh->getFaceIDs(&cellsIDs);

        kmds::TCoord minQual = HUGE_VALF;
        Kokkos::parallel_reduce("MinReduce",
                                cellsIDs.getNbElems(),
                                KOKKOS_LAMBDA(const int i, kmds::TCoord& minQ) {
                                    kmds::TCellID cid = cellsIDs.get(i);

                                    kmds::Face c = AMesh->getFace(cid);
                                    kmds::TCoord qual = c.scaledJacobian();

                                    minQ = std::min(minQ, qual);
                                },
                                Kokkos::Min<kmds::TCoord> (minQual));

        return minQual;
    }

    /*----------------------------------------------------------------------------*/
    kmds::TCoord
    Tools_computeScaledJacobian_3D(const kmds::Mesh* AMesh)
    {
        kmds::GrowingView <kmds::TCellID> cellsIDs("CELLSIDS", AMesh->getNbRegions());
        AMesh->getRegionIDs(&cellsIDs);

        kmds::TCoord minQual = HUGE_VALF;
        Kokkos::parallel_reduce("MinReduce",
                                cellsIDs.getNbElems(),
                                KOKKOS_LAMBDA(const int i, kmds::TCoord& minQ) {
                                    kmds::TCellID cid = cellsIDs.get(i);

                                    kmds::Region c = AMesh->getRegion(cid);
                                    kmds::TCoord qual = c.scaledJacobian();

                                    minQ = std::min(minQ, qual);
                                },
                                Kokkos::Min<kmds::TCoord> (minQual));

        return minQual;
    }

/*----------------------------------------------------------------------------*/
    kmds::TCoord
    Tools_computeScaledJacobian_2D(kmds::Mesh* AMesh,
                                   kmds::Variable<kmds::TCoord>* AVarQuality)
    {
        kmds::GrowingView <kmds::TCellID> cellsIDs("CELLSIDS", AMesh->getNbFaces());
        AMesh->getFaceIDs(&cellsIDs);

        kmds::TCoord minQual = HUGE_VALF;
        Kokkos::parallel_reduce("MinReduce",
                                cellsIDs.getNbElems(),
                                KOKKOS_LAMBDA(const int i, kmds::TCoord& minQ) {
                                    kmds::TCellID cid = cellsIDs.get(i);

                                    kmds::Face c = AMesh->getFace(cid);
                                    kmds::TCoord qual = c.scaledJacobian();

                                    (*AVarQuality)[cid] = qual;

                                    minQ = std::min(minQ, qual);
                                },
                                Kokkos::Min<kmds::TCoord> (minQual));

        return minQual;
    }

/*----------------------------------------------------------------------------*/
    kmds::TCoord
    Tools_computeScaledJacobian_3D(kmds::Mesh* AMesh,
                                   kmds::Variable<kmds::TCoord>* AVarQuality)
    {
        kmds::GrowingView <kmds::TCellID> cellsIDs("CELLSIDS", AMesh->getNbRegions());
        AMesh->getRegionIDs(&cellsIDs);

        kmds::TCoord minQual = HUGE_VALF;
        Kokkos::parallel_reduce("MinReduce",
                                cellsIDs.getNbElems(),
                                KOKKOS_LAMBDA(const int i, kmds::TCoord& minQ) {
                                    kmds::TCellID cid = cellsIDs.get(i);

                                    kmds::Region c = AMesh->getRegion(cid);
                                    kmds::TCoord qual = c.scaledJacobian();

                                    (*AVarQuality)[cid] = qual;

                                    minQ = std::min(minQ, qual);

                                },
                                Kokkos::Min<kmds::TCoord> (minQual));

        return minQual;
    }
/*----------------------------------------------------------------------------*/
    kmds::TCoord
    Tools_computeDistance_xD(const kmds::GrowingView<kmds::TCellID>* ASelectionNodes,
                             kmds::Mesh* AMesh,
                             const kmds::Variable<gmds::math::Point>* AVarTargetPosition,
                             kmds::Variable<kmds::TCoord>* AVarDist)
    {
        AVarDist->setValuesTo(0.);

        kmds::TCoord maxDist = -HUGE_VALF;
        Kokkos::parallel_reduce("MaxReduce",
                                ASelectionNodes->getNbElems(),
                                KOKKOS_LAMBDA(const int i, kmds::TCoord& maxD) {
                                    kmds::TCellID nid = ASelectionNodes->get(i);

                                    gmds::math::Point pt = AMesh->getNodeLocation(nid);
                                    kmds::TCoord dist = pt.distance((*AVarTargetPosition)[nid]);

                                    (*AVarDist)[nid] = dist;

//                                    maxD = std::max(maxD, dist);
maxD += dist;

                                },
//                                Kokkos::Max<kmds::TCoord> (maxDist));
maxDist);

        return maxDist;
    }
/*----------------------------------------------------------------------------*/
    void
    Tools_callGETME_2D(kmds::Mesh* AMesh,
                       const kmds::GrowingView<kmds::TCellID>* ASelectionFixedNodes,
                       const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation)
    {
        throw gmds::GMDSException("Not yet implemented!");

        /*

        kmds::TSize nbNodes = AMesh->getNbNodes();
        kmds::GrowingView<kmds::TCellID> nodeIDs("NODES", nbNodes);
        AMesh->getNodeIDsSeq(&nodeIDs);

        kmds::Variable<int >* varkmds2gmdsID =
                AMesh->createVariable<int >(gmds::NullID, kmds::KMDS_NODE, "tmp_varkmds2gmds");

        kmds::TSize nbCells = AMesh->getNbFaces();
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", nbCells);
        AMesh->getFaceIDsSeq(&cellIDs);

        gmds::Mesh gmdsMesh(gmds::DIM2|gmds::N|gmds::F|gmds::F2N);
        gmdsMesh.initializeGeometryClassification();
        gmds::Variable<gmds::cad::GeomEntity* >* nodesClassification = gmdsMesh.getGeometricClassification(0);

        for(int i=0; i<nbNodes; i++) {

            kmds::TCellID nid = nodeIDs.get(i);
            gmds::math::Point pt = AMesh->getNodeLocation(nid);

            gmds::Node n = gmdsMesh.newNode(pt);
            (*varkmds2gmdsID)[nid] = n.id();

            if ((*AVarNodeGeomAssociation)[nid] != reinterpret_cast<std::uintptr_t>(nullptr)) {
                (*nodesClassification)[n.id()] = reinterpret_cast<gmds::cad::GeomEntity *>((*AVarNodeGeomAssociation)[nid]);
            }
        }

        int markFixedNodes = gmdsMesh.newMark<gmds::Node> ();
        for(int i=0; i<ASelectionFixedNodes->getNbElems(); i++) {
            kmds::TCellID nid = ASelectionFixedNodes->get(i);

            gmds::Node n = gmdsMesh.get<gmds::Node> ((*varkmds2gmdsID)[nid]);

            gmdsMesh.mark(n, markFixedNodes);
        }

        for(int i=0; i<nbCells; i++) {

            kmds::TCellID cid = cellIDs.get(i);
            kmds::Face c = AMesh->getFace(cid);

            Kokkos::View<kmds::TCellID*> nids;
            c.nodeIds(nids);

            gmdsMesh.newQuad((*varkmds2gmdsID)[nids(0)],
                             (*varkmds2gmdsID)[nids(1)],
                             (*varkmds2gmdsID)[nids(2)],
                             (*varkmds2gmdsID)[nids(3)]
                    );
        }

        gmds::GETMe getme(gmdsMesh);
        getme.execSimult2D(
                10,
                0.8,
                0.5,
                0.,
                1.,
                1.,
                true,
                true,
                false,
                true,
//false,
                markFixedNodes
        );


        for(int i=0; i<nbNodes; i++) {

            kmds::TCellID nid = nodeIDs.get(i);
            gmds::Node n = gmdsMesh.get<gmds::Node> ((*varkmds2gmdsID)[nid]);

            AMesh->setNodeLocation(nid, n.getPoint());
        }


        // clean-up
        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varkmds2gmds");

         */
    }
/*----------------------------------------------------------------------------*/
    void
    Tools_callGETME_3D(kmds::Mesh* AMesh,
                       const kmds::GrowingView<kmds::TCellID>* ASelectionFixedNodes,
                       const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation)
    {
        throw gmds::GMDSException("Not yet implemented!");

        /*

        kmds::TSize nbNodes = AMesh->getNbNodes();
        kmds::GrowingView<kmds::TCellID> nodeIDs("NODES", nbNodes);
        AMesh->getNodeIDsSeq(&nodeIDs);

        kmds::Variable<int >* varkmds2gmdsID =
                AMesh->createVariable<int >(gmds::NullID, kmds::KMDS_NODE, "tmp_varkmds2gmds");

        kmds::TSize nbCells = AMesh->getNbRegions();
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", nbCells);
        AMesh->getRegionIDsSeq(&cellIDs);

        gmds::Mesh gmdsMesh(gmds::DIM3|gmds::N|gmds::R|gmds::R2N);
        gmdsMesh.initializeGeometryClassification();
        gmds::Variable<gmds::cad::GeomEntity* >* nodesClassification = gmdsMesh.getGeometricClassification(0);

        for(int i=0; i<nbNodes; i++) {

            kmds::TCellID nid = nodeIDs.get(i);
            gmds::math::Point pt = AMesh->getNodeLocation(nid);

            gmds::Node n = gmdsMesh.newNode(pt);
            (*varkmds2gmdsID)[nid] = n.id();

            if ((*AVarNodeGeomAssociation)[nid] != reinterpret_cast<std::uintptr_t>(nullptr)) {
                (*nodesClassification)[n.id()] = reinterpret_cast<gmds::cad::GeomEntity *>((*AVarNodeGeomAssociation)[nid]);
            }
        }

        int markFixedNodes = gmdsMesh.newMark<gmds::Node> ();
        for(int i=0; i<ASelectionFixedNodes->getNbElems(); i++) {
            kmds::TCellID nid = ASelectionFixedNodes->get(i);

            gmds::Node n = gmdsMesh.get<gmds::Node> ((*varkmds2gmdsID)[nid]);

            gmdsMesh.mark(n, markFixedNodes);
        }

        for(int i=0; i<nbCells; i++) {

            kmds::TCellID cid = cellIDs.get(i);
            kmds::Region c = AMesh->getRegion(cid);

            Kokkos::View<kmds::TCellID*> nids;
            c.nodeIds(nids);

            gmdsMesh.newHex((*varkmds2gmdsID)[nids(0)],
                            (*varkmds2gmdsID)[nids(1)],
                            (*varkmds2gmdsID)[nids(2)],
                            (*varkmds2gmdsID)[nids(3)],
                            (*varkmds2gmdsID)[nids(4)],
                            (*varkmds2gmdsID)[nids(5)],
                            (*varkmds2gmdsID)[nids(6)],
                            (*varkmds2gmdsID)[nids(7)]
            );
        }

        gmds::GETMe getme(gmdsMesh);
        getme.execSimult(
                100,
                0.8,
                0.5,
                0.,
                1.,
                1.,
                true,
                true,
                false,
                true,
                markFixedNodes
        );


        for(int i=0; i<nbNodes; i++) {

            kmds::TCellID nid = nodeIDs.get(i);
            gmds::Node n = gmdsMesh.get<gmds::Node> ((*varkmds2gmdsID)[nid]);

            AMesh->setNodeLocation(nid, n.getPoint());
        }


        // clean-up
        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varkmds2gmds");

         */
    }
/*----------------------------------------------------------------------------*/
    void
    Tools_write_2D(kmds::Mesh* AMesh,
                   const elg3d::FracPres* Afp,
                   const elg3d::MaterialAssignment *Ama,
                   const kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads,
                   const std::string AFileName)
    {
        kmds::TSize nbCells = AMesh->getNbFaces();
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", nbCells);
        AMesh->getFaceIDsSeq(&cellIDs);

        // build variable based on fracpres
        std::map<int, std::string> materialsfp = Afp->getMaterialList();
        for(auto mat: materialsfp) {

            int matindex = mat.first;

            std::string varname = std::string("fp_mat_") + mat.second;
            kmds::Variable<double>* var =
                    AMesh->createVariable<double>(-HUGE_VALF, kmds::KMDS_FACE, varname);

            for(int icell=0; icell<nbCells; icell++) {
                kmds::TCellID cid = cellIDs.get(icell);
                (*var)[cid] = Afp->getFracPres(matindex, cid);
            }
        }

        // build variable based on materialassignment
        kmds::Variable<int>* varma =
                AMesh->createVariable<int>(-1, kmds::KMDS_FACE, "materialassignment");

        for(int icell=0; icell<nbCells; icell++) {
            kmds::TCellID cid = cellIDs.get(icell);
            (*varma)[cid] = Ama->getMaterial(cid);
        }

        std::map<int, std::string> materialsma = Ama->getMaterialList();
        for(auto mat: materialsfp) {

            int matindex = mat.first;

            std::string varname = std::string("ma_mat_") + mat.second;
            kmds::Variable<int>* var =
                    AMesh->createVariable<int>(0, kmds::KMDS_FACE, varname);

            for(int icell=0; icell<nbCells; icell++) {
                kmds::TCellID cid = cellIDs.get(icell);
                if(Ama->getMaterial(cid) == matindex) {
                    (*var)[cid] = 1;
                }
            }
        }

        // build variables of vector based on gradients
        for(auto mat: materialsfp) {

            int matindex = mat.first;

            std::string varname = std::string("gradvec_mat_") + mat.second;
            kmds::Variable<gmds::math::Vector>* var =
                    AMesh->createVariable<gmds::math::Vector>(gmds::math::Vector (0., 0., 0.), kmds::KMDS_FACE, varname);

            for(int icell=0; icell<nbCells; icell++) {
                kmds::TCellID cid = cellIDs.get(icell);

                int indexfound = -1;

                for(int imat=0; imat<MaterialGradientComputation_MAXNBMATPERCELLBALL; imat++) {
                    if ((*AVarGrads)[cid].m_matindex[imat] == matindex) {
                        indexfound = imat;
                    }
                }
                if(indexfound != -1) {
                    (*var)[cid] = (*AVarGrads)[cid].m_grad[indexfound];
                }
            }
        }

        kmds::VTKWriter<kmds::Mesh> w(*AMesh);
        w.write(AFileName, kmds::F);

        // delete the created variables
        for(auto mat: materialsfp) {

            std::string varnamefp = std::string("fp_mat_") + mat.second;
            AMesh->deleteVariable(kmds::KMDS_FACE, varnamefp);
            std::string varnamegrad = std::string("gradvec_mat_") + mat.second;
            AMesh->deleteVariable(kmds::KMDS_FACE, varnamegrad);
        }

        AMesh->deleteVariable(kmds::KMDS_FACE, varma->getName());

        for(auto mat: materialsfp) {

            std::string varname = std::string("ma_mat_") + mat.second;
            AMesh->deleteVariable(kmds::KMDS_FACE, varname);
        }
    }
    /*----------------------------------------------------------------------------*/
    void
    Tools_write_2D(kmds::Mesh* AMesh,
                   const elg3d::FracPres* Afp,
                   const kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads,
                   const std::string AFileName)
    {
        kmds::TSize nbCells = AMesh->getNbFaces();
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", nbCells);
        AMesh->getFaceIDsSeq(&cellIDs);

        // build variable based on fracpres
        std::map<int, std::string> materialsfp = Afp->getMaterialList();
        for(auto mat: materialsfp) {

            int matindex = mat.first;

            std::string varname = std::string("fp_mat_") + mat.second;
            kmds::Variable<double>* var =
                    AMesh->createVariable<double>(-HUGE_VALF, kmds::KMDS_FACE, varname);

            for(int icell=0; icell<nbCells; icell++) {
                kmds::TCellID cid = cellIDs.get(icell);
                (*var)[cid] = Afp->getFracPres(matindex, cid);
            }
        }

        // build variables of vector based on gradients
        for(auto mat: materialsfp) {

            int matindex = mat.first;

            std::string varname = std::string("gradvec_mat_") + mat.second;
            kmds::Variable<gmds::math::Vector>* var =
                    AMesh->createVariable<gmds::math::Vector>(gmds::math::Vector (0., 0., 0.), kmds::KMDS_FACE, varname);

            for(int icell=0; icell<nbCells; icell++) {
                kmds::TCellID cid = cellIDs.get(icell);

                int indexfound = -1;

                for(int imat=0; imat<MaterialGradientComputation_MAXNBMATPERCELLBALL; imat++) {
                    if ((*AVarGrads)[cid].m_matindex[imat] == matindex) {
                        indexfound = imat;
                    }
                }
                if(indexfound != -1) {
                    (*var)[cid] = (*AVarGrads)[cid].m_grad[indexfound];
                }
            }
        }

        kmds::VTKWriter<kmds::Mesh> w(*AMesh);
        w.write(AFileName, kmds::F);

        // delete the created variables
        for(auto mat: materialsfp) {

            std::string varnamefp = std::string("fp_mat_") + mat.second;
            AMesh->deleteVariable(kmds::KMDS_FACE, varnamefp);
            std::string varnamegrad = std::string("gradvec_mat_") + mat.second;
            AMesh->deleteVariable(kmds::KMDS_FACE, varnamegrad);
        }

//        for(auto mat: materialsfp) {
//
//            std::string varname = std::string("ma_mat_") + mat.second;
//            AMesh->deleteVariable(kmds::KMDS_FACE, varname);
//        }
    }

/*----------------------------------------------------------------------------*/
    void
    Tools_write_3D(kmds::Mesh* AMesh,
                                         const elg3d::FracPres* Afp,
                                         const elg3d::MaterialAssignment *Ama,
                                         const kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads,
                                         const std::string AFileName)
    {
        kmds::TSize nbCells = AMesh->getNbRegions();
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", nbCells);
        AMesh->getRegionIDsSeq(&cellIDs);

        // build variable based on fracpres
        std::map<int, std::string> materialsfp = Afp->getMaterialList();
        for(auto mat: materialsfp) {

            int matindex = mat.first;

            std::string varname = std::string("fp_mat_") + mat.second;
            kmds::Variable<double>* var =
                    AMesh->createVariable<double>(-HUGE_VALF, kmds::KMDS_REGION, varname);

            for(int icell=0; icell<nbCells; icell++) {
                kmds::TCellID cid = cellIDs.get(icell);
                (*var)[cid] = Afp->getFracPres(matindex, cid);
            }
        }

        // build variable based on materialassignment
        kmds::Variable<int>* varma =
                AMesh->createVariable<int>(-1, kmds::KMDS_REGION, "materialassignment");

        for(int icell=0; icell<nbCells; icell++) {
            kmds::TCellID cid = cellIDs.get(icell);
            (*varma)[cid] = Ama->getMaterial(cid);
        }

        std::map<int, std::string> materialsma = Ama->getMaterialList();
        for(auto mat: materialsfp) {

            int matindex = mat.first;

            std::string varname = std::string("ma_mat_") + mat.second;
            kmds::Variable<int>* var =
                    AMesh->createVariable<int>(0, kmds::KMDS_REGION, varname);

            for(int icell=0; icell<nbCells; icell++) {
                kmds::TCellID cid = cellIDs.get(icell);
                if(Ama->getMaterial(cid) == matindex) {
                    (*var)[cid] = 1;
                }
            }
        }

        // build variables of vector based on gradients
        for(auto mat: materialsfp) {

            int matindex = mat.first;

            std::string varname = std::string("gradvec_mat_") + mat.second;
            kmds::Variable<gmds::math::Vector>* var =
                    AMesh->createVariable<gmds::math::Vector>(gmds::math::Vector (0., 0., 0.), kmds::KMDS_REGION, varname);

            for(int icell=0; icell<nbCells; icell++) {
                kmds::TCellID cid = cellIDs.get(icell);

                int indexfound = -1;

                for(int imat=0; imat<MaterialGradientComputation_MAXNBMATPERCELLBALL; imat++) {
                    if ((*AVarGrads)[cid].m_matindex[imat] == matindex) {
                        indexfound = imat;
                    }
                }
                if(indexfound != -1) {

                    (*var)[cid] = (*AVarGrads)[cid].m_grad[indexfound];
                }
            }
        }


        kmds::VTKWriter<kmds::Mesh> w(*AMesh);
        w.write(AFileName, kmds::R);

        // delete the created variables
        for(auto mat: materialsfp) {

            std::string varnamefp = std::string("fp_mat_") + mat.second;
            AMesh->deleteVariable(kmds::KMDS_REGION, varnamefp);
            std::string varnamegrad = std::string("gradvec_mat_") + mat.second;
            AMesh->deleteVariable(kmds::KMDS_REGION, varnamegrad);
        }

        AMesh->deleteVariable(kmds::KMDS_REGION, varma->getName());

        for(auto mat: materialsfp) {

            std::string varname = std::string("ma_mat_") + mat.second;
            AMesh->deleteVariable(kmds::KMDS_REGION, varname);
        }
    }

    /*----------------------------------------------------------------------------*/
    void
    Tools_write_3D(kmds::Mesh* AMesh,
                   const elg3d::FracPres* Afp,
                   const kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads,
                   const std::string AFileName)
    {
        kmds::TSize nbCells = AMesh->getNbRegions();
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", nbCells);
        AMesh->getRegionIDsSeq(&cellIDs);

        // build variable based on fracpres
        std::map<int, std::string> materialsfp = Afp->getMaterialList();
        for(auto mat: materialsfp) {

            int matindex = mat.first;

            std::string varname = std::string("fp_mat_") + mat.second;
            kmds::Variable<double>* var =
                    AMesh->createVariable<double>(-HUGE_VALF, kmds::KMDS_REGION, varname);

            for(int icell=0; icell<nbCells; icell++) {
                kmds::TCellID cid = cellIDs.get(icell);
                (*var)[cid] = Afp->getFracPres(matindex, cid);
            }
        }

        // build variables of vector based on gradients
        for(auto mat: materialsfp) {

            int matindex = mat.first;

            std::string varname = std::string("gradvec_mat_") + mat.second;
            kmds::Variable<gmds::math::Vector>* var =
                    AMesh->createVariable<gmds::math::Vector>(gmds::math::Vector (0., 0., 0.), kmds::KMDS_REGION, varname);

            for(int icell=0; icell<nbCells; icell++) {
                kmds::TCellID cid = cellIDs.get(icell);

                int indexfound = -1;

                for(int imat=0; imat<MaterialGradientComputation_MAXNBMATPERCELLBALL; imat++) {
                    if ((*AVarGrads)[cid].m_matindex[imat] == matindex) {
                        indexfound = imat;
                    }
                }
                if(indexfound != -1) {

                    (*var)[cid] = (*AVarGrads)[cid].m_grad[indexfound];
                }
            }
        }


        kmds::VTKWriter<kmds::Mesh> w(*AMesh);
        w.write(AFileName, kmds::R);

        // delete the created variables
        for(auto mat: materialsfp) {

            std::string varnamefp = std::string("fp_mat_") + mat.second;
            AMesh->deleteVariable(kmds::KMDS_REGION, varnamefp);
            std::string varnamegrad = std::string("gradvec_mat_") + mat.second;
            AMesh->deleteVariable(kmds::KMDS_REGION, varnamegrad);
        }
    }

/*----------------------------------------------------------------------------*/
    void
    Tools_write_2D(kmds::Mesh* AMesh,
                   const elg3d::FracPres* Afp,
                   const elg3d::MaterialAssignment *Ama,
                   const std::string AFileName)
    {
//        return;
        kmds::TSize nbCells = AMesh->getNbFaces();
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", nbCells);
        AMesh->getFaceIDsSeq(&cellIDs);

        // build variable based on fracpres
        std::map<int, std::string> materialsfp = Afp->getMaterialList();
        for(auto mat: materialsfp) {

            int matindex = mat.first;

            std::string varname = std::string("fp_mat_") + mat.second;
            kmds::Variable<double>* var =
                    AMesh->createVariable<double>(-HUGE_VALF, kmds::KMDS_FACE, varname);

            for(int icell=0; icell<nbCells; icell++) {
                kmds::TCellID cid = cellIDs.get(icell);
                (*var)[cid] = Afp->getFracPres(matindex, cid);
            }
        }

        // build variable based on materialassignment
        kmds::Variable<int>* varma =
                AMesh->createVariable<int>(-1, kmds::KMDS_FACE, "materialassignment");

        for(int icell=0; icell<nbCells; icell++) {
            kmds::TCellID cid = cellIDs.get(icell);
            (*varma)[cid] = Ama->getMaterial(cid);
        }

        std::map<int, std::string> materialsma = Ama->getMaterialList();
        for(auto mat: materialsfp) {

            int matindex = mat.first;

            std::string varname = std::string("ma_mat_") + mat.second;
            kmds::Variable<int>* var =
                    AMesh->createVariable<int>(0, kmds::KMDS_FACE, varname);

            for(int icell=0; icell<nbCells; icell++) {
                kmds::TCellID cid = cellIDs.get(icell);
                if(Ama->getMaterial(cid) == matindex) {
                    (*var)[cid] = 1;
                }
            }
        }


        kmds::VTKWriter<kmds::Mesh> w(*AMesh);
        w.write(AFileName, kmds::F);


        // delete the created variables
        for(auto mat: materialsfp) {

            std::string varnamefp = std::string("fp_mat_") + mat.second;
            AMesh->deleteVariable(kmds::KMDS_FACE, varnamefp);
        }

        AMesh->deleteVariable(kmds::KMDS_FACE, varma->getName());

        for(auto mat: materialsfp) {

            std::string varname = std::string("ma_mat_") + mat.second;
            AMesh->deleteVariable(kmds::KMDS_FACE, varname);
        }
    }
    /*----------------------------------------------------------------------------*/
    void
    Tools_write_lite_2D(kmds::Mesh* AMesh,
                        const elg3d::MaterialAssignment *Ama,
                        const std::string AFileName)
    {
//        return;

        kmds::TSize nbCells = AMesh->getNbFaces();
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", nbCells);
        AMesh->getFaceIDsSeq(&cellIDs);

        // build variable based on materialassignment
        kmds::Variable<int>* varma =
                AMesh->createVariable<int>(-1, kmds::KMDS_FACE, "materialassignment");

        for(int icell=0; icell<nbCells; icell++) {
            kmds::TCellID cid = cellIDs.get(icell);
            (*varma)[cid] = Ama->getMaterial(cid);
        }


        kmds::VTKWriter<kmds::Mesh> w(*AMesh);
        w.write(AFileName, kmds::F);


        // delete the created variable
        AMesh->deleteVariable(kmds::KMDS_FACE, varma->getName());
    }
/*----------------------------------------------------------------------------*/
    void
    Tools_write_3D(kmds::Mesh* AMesh,
                   const elg3d::FracPres* Afp,
                   const elg3d::MaterialAssignment *Ama,
                   const std::string AFileName)
    {
//        return;

        kmds::TSize nbCells = AMesh->getNbRegions();
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", nbCells);
        AMesh->getRegionIDsSeq(&cellIDs);

        // build variable based on fracpres
        std::map<int, std::string> materialsfp = Afp->getMaterialList();
        for(auto mat: materialsfp) {

            int matindex = mat.first;

            std::string varname = std::string("fp_mat_") + mat.second;
            kmds::Variable<double>* var =
                    AMesh->createVariable<double>(-HUGE_VALF, kmds::KMDS_REGION, varname);

            for(int icell=0; icell<nbCells; icell++) {
                kmds::TCellID cid = cellIDs.get(icell);
                (*var)[cid] = Afp->getFracPres(matindex, cid);
            }
        }

        // build variable based on materialassignment
        kmds::Variable<int>* varma =
                AMesh->createVariable<int>(-1, kmds::KMDS_REGION, "materialassignment");

        for(int icell=0; icell<nbCells; icell++) {
            kmds::TCellID cid = cellIDs.get(icell);
            (*varma)[cid] = Ama->getMaterial(cid);
        }

        std::map<int, std::string> materialsma = Ama->getMaterialList();
        for(auto mat: materialsfp) {

            int matindex = mat.first;

            std::string varname = std::string("ma_mat_") + mat.second;
            kmds::Variable<int>* var =
                    AMesh->createVariable<int>(0, kmds::KMDS_REGION, varname);

            for(int icell=0; icell<nbCells; icell++) {
                kmds::TCellID cid = cellIDs.get(icell);
                if(Ama->getMaterial(cid) == matindex) {
                    (*var)[cid] = 1;
                }
            }
        }


        kmds::VTKWriter<kmds::Mesh> w(*AMesh);
        w.write(AFileName, kmds::R);

        // delete the created variables
        for(auto mat: materialsfp) {

            std::string varnamefp = std::string("fp_mat_") + mat.second;
            AMesh->deleteVariable(kmds::KMDS_REGION, varnamefp);
        }

        AMesh->deleteVariable(kmds::KMDS_REGION, varma->getName());

        for(auto mat: materialsfp) {

            std::string varname = std::string("ma_mat_") + mat.second;
            AMesh->deleteVariable(kmds::KMDS_REGION, varname);
        }
    }
    /*----------------------------------------------------------------------------*/
    void
    Tools_write_lite_3D(kmds::Mesh* AMesh,
                        const elg3d::MaterialAssignment *Ama,
                        const std::string AFileName)
    {
//        return;

        kmds::TSize nbCells = AMesh->getNbRegions();
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", nbCells);
        AMesh->getRegionIDsSeq(&cellIDs);

        // build variable based on materialassignment
        kmds::Variable<int>* varma =
                AMesh->createVariable<int>(-1, kmds::KMDS_REGION, "materialassignment");

        for(int icell=0; icell<nbCells; icell++) {
            kmds::TCellID cid = cellIDs.get(icell);
            (*varma)[cid] = Ama->getMaterial(cid);
        }


        kmds::VTKWriter<kmds::Mesh> w(*AMesh);
        w.write(AFileName, kmds::R);


        // delete the created variable
        AMesh->deleteVariable(kmds::KMDS_REGION, varma->getName());
    }
    /*----------------------------------------------------------------------------*/
    void
    Tools_write_interfaces_3D(const kmds::GrowingView<kmds::TCellID>* AInterfaceFaces,
                              const kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
                              kmds::Mesh* AMesh,
                              const kmds::Connectivity* c_F2C,
                              const elg3d::MaterialAssignment *Ama,
                              const std::string AFileName)
    {
//        return;

        const kmds::TCellID nbNodes = AInterfaceNodes->getNbElems();
        const kmds::TCellID nbFaces = AInterfaceFaces->getNbElems();

        kmds::Mesh mesh;
        mesh.updateNodeCapacity(nbNodes);
        mesh.updateFaceCapacity(nbFaces);

        std::map<kmds::TCellID, kmds::TCellID> current2tmpNodeIDs;

        kmds::TCellID nid0 = mesh.addNodes(nbNodes);
        for(int i=0; i<nbNodes; i++) {
            kmds::TCellID nid = AInterfaceNodes->get(i);

            kmds::TCoord xyz[3];
            AMesh->getNodeLocation(nid, xyz[0], xyz[1], xyz[2]);

            mesh.setNodeLocation(nid0 + i, xyz[0], xyz[1], xyz[2]);

            current2tmpNodeIDs.emplace(nid, nid0 + i);
        }

        std::map<kmds::TCellID, kmds::TCellID> current2tmpFaceIDs;

        kmds::TCellID fid0 = mesh.addQuads(nbFaces);
        for(int i=0; i<nbFaces; i++) {
            kmds::TCellID fid = AInterfaceFaces->get(i);

            kmds::Face f = AMesh->getFace(fid);
            Kokkos::View<kmds::TCellID *> nids;
            f.nodeIds(nids);

            kmds::TCellID nids_tmp[4];
            for(int n=0; n<nids.extent(0); n++) {
                nids_tmp[n] = current2tmpNodeIDs[nids(n)];
            }

            kmds::Face f_tmp = mesh.getFace(fid0 + i);

            f_tmp.setNodes(nids_tmp, nids.extent(0));

            current2tmpNodeIDs.emplace(fid, fid0 + i);
        }

        kmds::VTKWriter<kmds::Mesh> w(mesh);
        w.write(AFileName, kmds::F);

    }
/*----------------------------------------------------------------------------*/
    void
    Tools_read_fracpres_vtk_3D(kmds::Mesh *AMesh,
                               elg3d::FracPres *Afp,
                               const std::string AFileName,
                               bool ADeduceVoid)
    {
        gmds::Mesh mesh_tmp(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::R|gmds::R2N));
        gmds::IGMeshIOService ioService(&mesh_tmp);
        gmds::VTKReader vtkReader(&ioService);
        vtkReader.setCellOptions(gmds::N|gmds::R);
        vtkReader.setDataOptions(gmds::R);

        vtkReader.read(AFileName);

        kmds::TSize nbNodes = mesh_tmp.getNbNodes();
        kmds::TSize nbCells = mesh_tmp.getNbRegions();

        AMesh->updateNodeCapacity(nbNodes);
        AMesh->updateRegionCapacity(nbCells);

        for (gmds::Mesh::nodes_iterator it_nodes = mesh_tmp.nodes_begin(); it_nodes != mesh_tmp.nodes_end(); ++it_nodes) {
            gmds::Node n = mesh_tmp.get<gmds::Node>(*it_nodes);

            kmds::TCellID firstid = AMesh->newNode(n.X(), n.Y(), n.Z());
        }

        for (gmds::Mesh::regions_iterator it_regions = mesh_tmp.regions_begin(); it_regions != mesh_tmp.regions_end(); ++it_regions) {
            gmds::Region r = mesh_tmp.get<gmds::Region>(*it_regions);

            std::vector<gmds::TCellID> nids = r.getIDs<gmds::Node>();

            kmds::TCellID firstid = AMesh->newHexahedron(nids[0],
                                                         nids[1],
                                                         nids[2],
                                                         nids[3],
                                                         nids[4],
                                                         nids[5],
                                                         nids[6],
                                                         nids[7]);
        }

        kmds::Variable<double>* AVoidSurfVol =
                AMesh->createVariable<double>(1, kmds::KMDS_REGION, "voidSurfVol");

        std::vector<gmds::VariableItf*> region_variables = mesh_tmp.getAllVariables(gmds::GMDS_REGION);

        std::cout<<"region_variables.size() "<<region_variables.size()<<std::endl;

        for (unsigned int i = 0; i < region_variables.size(); i++) {

            gmds::VariableItf *current_var = region_variables[i];
            switch (ioService.getType(current_var)) {
                case (gmds::IMeshIOService::var_double) : {
                    gmds::Variable<double> *v_double = dynamic_cast<gmds::Variable<double> *> (current_var);

                    int mat = Afp->createMaterial(v_double->getName().c_str());

                    for (gmds::Mesh::regions_iterator it = mesh_tmp.regions_begin(); it != mesh_tmp.regions_end(); ++it) {

                        Afp->setFracPres(mat, *it, (*v_double)[*it]);

                        (*AVoidSurfVol)[*it] = (*AVoidSurfVol)[*it] - (*v_double)[*it];
                    }
                }
                    break;
                default:
                    break;
            }
        }

        if(ADeduceVoid) {
            int mat = Afp->createMaterial("void");

            for (gmds::Mesh::regions_iterator it = mesh_tmp.regions_begin(); it != mesh_tmp.regions_end(); ++it) {
                Afp->setFracPres(mat, *it, (*AVoidSurfVol)[*it]);
            }
        }

        AMesh->deleteVariable(kmds::KMDS_REGION, "voidSurfVol");

    }

    /*----------------------------------------------------------------------------*/
    void
    Tools_read_fracpres_vtk_2D(kmds::Mesh *AMesh,
                               elg3d::FracPres *Afp,
                               const std::string AFileName,
                               bool ADeduceVoid)
    {
        gmds::Mesh mesh_tmp(gmds::MeshModel(gmds::DIM2|gmds::N|gmds::F|gmds::F2N));
        gmds::IGMeshIOService ioService(&mesh_tmp);
        gmds::VTKReader vtkReader(&ioService);


        kmds::TSize nbNodes = mesh_tmp.getNbNodes();
        kmds::TSize nbCells = mesh_tmp.getNbFaces();

        AMesh->updateNodeCapacity(nbNodes);
        AMesh->updateFaceCapacity(nbCells);

        for (gmds::Mesh::nodes_iterator it_nodes = mesh_tmp.nodes_begin(); it_nodes != mesh_tmp.nodes_end(); ++it_nodes) {
            gmds::Node n = mesh_tmp.get<gmds::Node>(*it_nodes);

            kmds::TCellID firstid = AMesh->newNode(n.X(), n.Y(), 0.);
        }

        for (gmds::Mesh::faces_iterator it_faces = mesh_tmp.faces_begin(); it_faces != mesh_tmp.faces_end(); ++it_faces) {
            gmds::Face f = mesh_tmp.get<gmds::Face>(*it_faces);
            std::vector<gmds::TCellID> nids = f.getIDs<gmds::Node>();

            kmds::TCellID firstid = AMesh->newQuad(nids[0],
                                                   nids[1],
                                                   nids[2],
                                                   nids[3]);
        }

        kmds::Variable<double>* AVoidSurfVol =
                AMesh->createVariable<double>(1, kmds::KMDS_FACE, "voidSurfVol");

        std::vector<gmds::VariableItf*> face_variables = mesh_tmp.getAllVariables(gmds::GMDS_FACE);
        for (unsigned int i = 0; i < face_variables.size(); i++) {

            gmds::VariableItf *current_var = face_variables[i];
            switch (ioService.getType(current_var)) {
                case (gmds::IMeshIOService::var_double) : {
                    gmds::Variable<double> *v_double = dynamic_cast<gmds::Variable<double> *> (current_var);

                    int mat = Afp->createMaterial(v_double->getName().c_str());

                    for (gmds::Mesh::faces_iterator it = mesh_tmp.faces_begin(); it != mesh_tmp.faces_end(); ++it) {

                        Afp->setFracPres(mat, *it, (*v_double)[*it]);

                        (*AVoidSurfVol)[*it] = (*AVoidSurfVol)[*it] - (*v_double)[*it];
                    }
                }
                    break;
                default:
                    break;
            }
        }

        if(ADeduceVoid) {
            int mat = Afp->createMaterial("void");

            for (gmds::Mesh::faces_iterator it = mesh_tmp.faces_begin(); it != mesh_tmp.faces_end(); ++it) {
                Afp->setFracPres(mat, *it, (*AVoidSurfVol)[*it]);
            }
        }

        AMesh->deleteVariable(kmds::KMDS_FACE, "voidSurfVol");

    }

    /*----------------------------------------------------------------------------*/
    void
    Tools_read_ma_vtk_3D(kmds::Mesh *AMesh,
                         elg3d::MaterialAssignment *Ama,
                         const std::string AFileName)
    {
        gmds::Mesh mesh_tmp(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::R|gmds::R2N));
        gmds::IGMeshIOService ioService(&mesh_tmp);
        gmds::VTKReader vtkReader(&ioService);
        vtkReader.setCellOptions(gmds::N|gmds::R);
        vtkReader.setDataOptions(gmds::R);

        vtkReader.read(AFileName);

        kmds::TSize nbNodes = mesh_tmp.getNbNodes();
        kmds::TSize nbCells = mesh_tmp.getNbRegions();

        AMesh->updateNodeCapacity(nbNodes);
        AMesh->updateRegionCapacity(nbCells);

        for (gmds::Mesh::nodes_iterator it_nodes = mesh_tmp.nodes_begin(); it_nodes != mesh_tmp.nodes_end(); ++it_nodes) {
            gmds::Node n = mesh_tmp.get<gmds::Node>(*it_nodes);

            kmds::TCellID firstid = AMesh->newNode(n.X(), n.Y(), n.Z());
        }

        for (gmds::Mesh::regions_iterator it_regions = mesh_tmp.regions_begin(); it_regions != mesh_tmp.regions_end(); ++it_regions) {
            gmds::Region r = mesh_tmp.get<gmds::Region>(*it_regions);

            std::vector<gmds::TCellID> nids = r.getIDs<gmds::Node>();

            kmds::TCellID firstid = AMesh->newHexahedron(nids[0],
                                                         nids[1],
                                                         nids[2],
                                                         nids[3],
                                                         nids[4],
                                                         nids[5],
                                                         nids[6],
                                                         nids[7]);
        }

        Ama->updateCapacity(nbCells);

        std::vector<gmds::VariableItf*> region_variables = mesh_tmp.getAllVariables(gmds::GMDS_REGION);
        std::cout<<"region_variables.size() "<<region_variables.size()<<std::endl;

        for (unsigned int i = 0; i < region_variables.size(); i++) {

            gmds::VariableItf *current_var = region_variables[i];
            switch (ioService.getType(current_var)) {
                case (gmds::IMeshIOService::var_int) : {
                    gmds::Variable<int> *v_int = dynamic_cast<gmds::Variable<int> *> (current_var);

//                    if((v_int->getName().length() > 3) && (v_int->getName().substr(0, 3) == std::string("ma_"))) {
//                        int mat = Ama->createMaterial(v_int->getName().substr(3, std::string::npos));
//
//                        for (gmds::Mesh::regions_iterator it = mesh_tmp.regions_begin();
//                             it != mesh_tmp.regions_end(); ++it) {
//
//                            if((*v_int)[*it] == 1) {
//
//                                Ama->setMaterial(mat, *it);
//                            }
//                        }
//                    }

                    if (v_int->getName() == std::string("materialassignment")) {

                        for (gmds::Mesh::regions_iterator it = mesh_tmp.regions_begin();
                             it != mesh_tmp.regions_end(); ++it) {

                            Ama->setMaterial((*v_int)[*it], *it);
                        }
                    }
                }
                    break;
                default:
                    break;
            }
        }
    }

    /*----------------------------------------------------------------------------*/
    void Tools_compute_fracpres_2D(kmds::Mesh* AMesh,
                                   gmds::Mesh* AgmdsMesh,
                                   bool ADeduceVoid,
                                   elg3d::FracPres* Afp)
    {
        // create the materials
        for (gmds::Mesh::group_iterator<gmds::Face> its = AgmdsMesh->groups_begin<gmds::Face>(); its != AgmdsMesh->groups_end<gmds::Face>(); its++)
        {
            std::string matName = (*its)->name();
            Afp->createMaterial(matName);

            std::cout<<"create mat "<<matName<<std::endl;
        }

        // is the void deduced ?
        int voidID = 0;
        if(ADeduceVoid) {
            std::string matName = std::string("void");
            Afp->createMaterial(matName);
            voidID = Afp->getMaterialID("void");
        }

        std::cout<<"number of materials "<<Afp->getNbMaterials()<<std::endl;

        // initialize the fracpres of the kmds mesh
        int nbCells = AMesh->getNbFaces();

        // this is the current fracpres, useful in the case where we have overlapping material meshes
        std::map<int, double> sumFracPres;

        for(int i=0; i<nbCells; i++) {
            if(ADeduceVoid) {
                Afp->setFracPres(voidID, i, 1.);
            }

            sumFracPres[i] = 0.;
        }


        // populate the volume fractions
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", nbCells);
        AMesh->getFaceIDs_dummy(&cellIDs);

        for (gmds::Mesh::group_iterator<gmds::Face> its = AgmdsMesh->groups_begin<gmds::Face>(); its != AgmdsMesh->groups_end<gmds::Face>(); its++)
        {
            std::string matName = (*its)->name();
            int matIndex = Afp->getMaterialID(matName);

            const std::vector<gmds::TCellID> cells = (*its)->cells();

            for(int i_gmds=0; i_gmds<cells.size(); i_gmds++) {
                gmds::Face c = AgmdsMesh->get<gmds::Face>(cells[i_gmds]);
                std::vector<gmds::Node> n = c.get<gmds::Node>();
                gmds::math::Triangle tri(n[0].getPoint(), n[1].getPoint(), n[2].getPoint());

                double xyz_min[3];
                double xyz_max[3];
                tri.computeBoundingBox(xyz_min, xyz_max);

                r2d_rvec2 vertices[3];
                vertices[0].x = n[0].X();
                vertices[0].y = n[0].Y();
                vertices[1].x = n[1].X();
                vertices[1].y = n[1].Y();
                vertices[2].x = n[2].X();
                vertices[2].y = n[2].Y();
                r2d_int numverts = 3;
                r2d_plane planes[3];
                r2d_poly_faces_from_verts(planes,  vertices, numverts);


                for(int i_kmds=0; i_kmds<nbCells; i_kmds++) {

                    kmds::Face f = AMesh->getFace(cellIDs.get(i_kmds));

                    Kokkos::View<kmds::TCellID*> nids;
                    f.nodeIds(nids);

                    kmds::TCoord xyz[4][3];
                    AMesh->getNodeLocation(nids[0], xyz[0][0], xyz[0][1], xyz[0][2]);
                    AMesh->getNodeLocation(nids[1], xyz[1][0], xyz[1][1], xyz[1][2]);
                    AMesh->getNodeLocation(nids[2], xyz[2][0], xyz[2][1], xyz[2][2]);
                    AMesh->getNodeLocation(nids[3], xyz[3][0], xyz[3][1], xyz[3][2]);

                    gmds::math::Point pt0(xyz[0][0], xyz[0][1], xyz[0][2]);
                    gmds::math::Point pt1(xyz[1][0], xyz[1][1], xyz[1][2]);
                    gmds::math::Point pt2(xyz[2][0], xyz[2][1], xyz[2][2]);
                    gmds::math::Point pt3(xyz[3][0], xyz[3][1], xyz[3][2]);
                    gmds::math::Quadrilateral q(pt0, pt1, pt2, pt3);

                    const double volsurf = q.area();

                    r2d_poly poly;
                    r2d_rvec2 verts[4];

                    verts[0].x = xyz[0][0];
                    verts[0].y = xyz[0][1];
                    verts[1].x = xyz[1][0];
                    verts[1].y = xyz[1][1];
                    verts[2].x = xyz[2][0];
                    verts[2].y = xyz[2][1];
                    verts[3].x = xyz[3][0];
                    verts[3].y = xyz[3][1];

                    r2d_init_poly(&poly, verts, 4);

                    r2d_clip(&poly, planes, 3);
                    r2d_int POLY_ORDER = 2;
                    r2d_real om[R2D_NUM_MOMENTS(POLY_ORDER)];
                    r2d_reduce(&poly, om, POLY_ORDER);

                    double fracpres = Afp->getFracPres(matIndex, i_kmds);

                    // add but do not forget to divide by cell area
                    double fracpresadd = om[0]/volsurf;

                    if(fracpresadd < 1.e-12) {
                        fracpresadd = 0.;
                    }

                    if(sumFracPres[i_kmds] + fracpresadd > 1.) {
                        fracpresadd = 1.- sumFracPres[i_kmds];
                        sumFracPres[i_kmds] = 1.;
                    } else {
                        sumFracPres[i_kmds] += fracpresadd;
                    }

                    fracpres += fracpresadd;

                    Afp->setFracPres(matIndex, i_kmds, fracpres);

                    if(ADeduceVoid) {
                        double fracpresVoid = Afp->getFracPres(voidID, i_kmds);
                        fracpresVoid -= fracpresadd;
                        if(fracpresVoid < 0.) {
                            fracpresVoid = 0.;
                        }

                        Afp->setFracPres(voidID,i_kmds,fracpresVoid);
                    }
                }
            }

        }

        Afp->removeTraceAmounts();
    }

    /*----------------------------------------------------------------------------*/
    void Tools_compute_fracpres_3D(kmds::Mesh* AMesh,
                                   gmds::Mesh* AgmdsMesh,
                                   bool ADeduceVoid,
                                   elg3d::FracPres* Afp)
    {
        // create the materials

        for (gmds::Mesh::group_iterator<gmds::Region> itv = AgmdsMesh->groups_begin<gmds::Region>(); itv != AgmdsMesh->groups_end<gmds::Region>(); itv++)
        {
            std::string matName = (*itv)->name();
            Afp->createMaterial(matName);

            std::cout<<"create mat "<<matName<<std::endl;
        }

        // is the void deduced ?
        int voidID = 0;
        if(ADeduceVoid) {
            std::string matName = std::string("void");
            Afp->createMaterial(matName);
            voidID = Afp->getMaterialID("void");
        }

        std::cout<<"number of materials "<<Afp->getNbMaterials()<<std::endl;

        // initialize the fracpres of the kmds mesh
        int nbCells = AMesh->getNbRegions();

        // this is the current fracpres, useful in the case where we have overlapping material meshes
        std::map<int, double> sumFracPres;

        for(int i=0; i<nbCells; i++) {
            if(ADeduceVoid) {
                Afp->setFracPres(voidID, i, 1.);
            }

            sumFracPres[i] = 0.;
        }


        // prepare the target mesh in the AABBtree
        // Axis-Aligned Bounding Box tree for the cells
        // The bounded value is the CellID
        GNode *aabbCellsTree;
        {
            GSList *list = NULL;


            kmds::GrowingView<kmds::TCellID> cellIDs_target("CELLS", AMesh->getNbRegions());
            AMesh->getRegionIDs(&cellIDs_target);

            for (int i = 0; i < cellIDs_target.getNbElems(); i++) {
                kmds::TCellID cid = cellIDs_target.get(i);

                kmds::Region c = AMesh->getRegion(cid);

                double minXYZ[3];
                double maxXYZ[3];

                c.computeBoundingBox(minXYZ, maxXYZ);

                gpointer pointer = GINT_TO_POINTER(cid);
                GtsBBox *bbox = gts_bbox_new(
                        gts_bbox_class(),
                        pointer,
                        minXYZ[0], minXYZ[1], minXYZ[2],
                        maxXYZ[0], maxXYZ[1], maxXYZ[2]);

                list = g_slist_prepend(list, bbox);
            }
            aabbCellsTree = gts_bb_tree_new(list);

            g_slist_free(list);
        }

        gpointer pointer = NULL;
        GtsBBox* bbox = gts_bbox_new(
                gts_bbox_class (),
                pointer,
                0.,0.,0.,
                0.,0.,0.);



        // populate the volume fractions
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", nbCells);
        AMesh->getRegionIDs_dummy(&cellIDs);

        std::cout<<"nbVolumes "<<AgmdsMesh->getNbGroups<gmds::Region>()<<std::endl;

        for (gmds::Mesh::group_iterator<gmds::Region> itv = AgmdsMesh->groups_begin<gmds::Region>(); itv != AgmdsMesh->groups_end<gmds::Region>(); itv++)
        {

            std::cout<<"vol frac for mat "<<(*itv)->name()<<std::endl;

            std::string matName = (*itv)->name();
            int matIndex = Afp->getMaterialID(matName);

            const std::vector<gmds::TCellID > cells = (*itv)->cells();

            for(int i_gmds=0; i_gmds<cells.size(); i_gmds++) {

                std::cout<<"i_gmds "<<i_gmds<<" of "<<cells.size()<<std::endl;

                gmds::Region c = AgmdsMesh->get<gmds::Region>(cells[i_gmds]);
                std::vector<gmds::Node> n = c.get<gmds::Node>();
                gmds::math::Tetrahedron tet(n[0].getPoint(), n[1].getPoint(), n[2].getPoint(), n[3].getPoint());

                double xyz_min_tet[3];
                double xyz_max_tet[3];
                tet.computeBoundingBox(xyz_min_tet, xyz_max_tet);

                r3d_rvec3 vertices[4];
                for(int i_n=0; i_n<4; i_n++) {
                    vertices[i_n].x = n[i_n].X();
                    vertices[i_n].y = n[i_n].Y();
                    vertices[i_n].z = n[i_n].Z();
                }

//                r3d_poly poly_tet;
//                r3d_init_tet(&poly_tet, vertices);
//
//                std::cout<<"isgood tet "<<r3d_is_good(&poly_tet)<<std::endl;
//                r3d_print(&poly_tet);

                r3d_plane planes[4];
                r3d_tet_faces_from_verts(planes,  vertices);


                // search for the cells of kmds mesh in the tree
                gts_bbox_set(
                    bbox,
                    pointer,
                    xyz_min_tet[0],xyz_min_tet[1],xyz_min_tet[2],
                    xyz_max_tet[0],xyz_max_tet[1],xyz_max_tet[2]);

                GSList* boxList = gts_bb_tree_overlap(aabbCellsTree,bbox);
                if(boxList == NULL) {
//                    throw kmds::KException("Tools_compute_fracpres_3D : a cell from source mesh does not intersect target mesh.");
                }

//                for(int i_kmds=0; i_kmds<nbCells; i_kmds++) {
//
//                    kmds::Region r = AMesh->getRegion(cellIDs.get(i_kmds));

                while (boxList != NULL) {

                    GtsBBox *box = (GtsBBox *) (boxList->data);
                    kmds::TCellID cid_target = GPOINTER_TO_INT(box->bounded);

                    kmds::Region c_target = AMesh->getRegion(cid_target);

                    Kokkos::View<kmds::TCellID*> nids;
                    c_target.nodeIds(nids);

                    kmds::TCoord xyz[8][3];
                    for(int i_n=0; i_n<8; i_n++) {
                        AMesh->getNodeLocation(nids[i_n], xyz[i_n][0], xyz[i_n][1], xyz[i_n][2]);
                    }

                    gmds::math::Point pts[8];
                    for(int i_n=0; i_n<8; i_n++) {
                        pts[i_n].setXYZ(xyz[i_n][0], xyz[i_n][1], xyz[i_n][2]);
                    }

                    gmds::math::Hexahedron h(pts[0], pts[1], pts[2], pts[3], pts[4], pts[5], pts[6], pts[7]);

                    double xyz_min_h[3];
                    double xyz_max_h[3];
                    h.computeBoundingBox(xyz_min_h, xyz_max_h);

                    if(!gmds::math::intersectBoundingBox(xyz_min_tet, xyz_max_tet, xyz_min_h, xyz_max_h)) {
                        continue;
                    }

                    const double volsurf = h.getVolume();

                    r3d_poly poly;
                    poly.nverts = 8;
//                    r3d_rvec3 verts[8];

                    for(int i_n=0; i_n<8; i_n++) {
                        poly.verts[i_n].pos.xyz[0] = xyz[i_n][0];
                        poly.verts[i_n].pos.xyz[1] = xyz[i_n][1];
                        poly.verts[i_n].pos.xyz[2] = xyz[i_n][2];
                    }

                    poly.verts[0].pnbrs[0] = 1;
                    poly.verts[0].pnbrs[1] = 4;
                    poly.verts[0].pnbrs[2] = 3;
                    poly.verts[1].pnbrs[0] = 2;
                    poly.verts[1].pnbrs[1] = 5;
                    poly.verts[1].pnbrs[2] = 0;
                    poly.verts[2].pnbrs[0] = 3;
                    poly.verts[2].pnbrs[1] = 6;
                    poly.verts[2].pnbrs[2] = 1;
                    poly.verts[3].pnbrs[0] = 0;
                    poly.verts[3].pnbrs[1] = 7;
                    poly.verts[3].pnbrs[2] = 2;

                    poly.verts[4].pnbrs[0] = 5;
                    poly.verts[4].pnbrs[1] = 7;
                    poly.verts[4].pnbrs[2] = 0;
                    poly.verts[5].pnbrs[0] = 6;
                    poly.verts[5].pnbrs[1] = 4;
                    poly.verts[5].pnbrs[2] = 1;
                    poly.verts[6].pnbrs[0] = 7;
                    poly.verts[6].pnbrs[1] = 5;
                    poly.verts[6].pnbrs[2] = 2;
                    poly.verts[7].pnbrs[0] = 4;
                    poly.verts[7].pnbrs[1] = 6;
                    poly.verts[7].pnbrs[2] = 3;

//                    std::cout<<"isgood "<<r3d_is_good(&poly)<<std::endl;
//                    r3d_print(&poly);
//                    exit(0);


                    r3d_clip(&poly, planes, 4);
                    r3d_int POLY_ORDER = 2;
                    r3d_real om[R3D_NUM_MOMENTS(POLY_ORDER)];
                    r3d_reduce(&poly, om, POLY_ORDER);

                    double fracpres = Afp->getFracPres(matIndex, cid_target);

                    // add but do not forget to divide by cell area
                    double fracpresadd = om[0]/volsurf;

//                    std::cout<<"INTERSECT "<<fracpresadd<<std::endl;

                    if(fracpresadd < 1.e-10) {
                        fracpresadd = 0.;
                    }

                    if(sumFracPres[cid_target] + fracpresadd > 1.) {
                        fracpresadd = 1.- sumFracPres[cid_target];
                        sumFracPres[cid_target] = 1.;
                    } else {
                        sumFracPres[cid_target] += fracpresadd;
                    }

                    fracpres += fracpresadd;

                    Afp->setFracPres(matIndex, cid_target, fracpres);

                    if(ADeduceVoid) {
                        double fracpresVoid = Afp->getFracPres(voidID, cid_target);
                        fracpresVoid -= fracpresadd;
                        if(fracpresVoid < 0.) {
                            fracpresVoid = 0.;
                        }

                        Afp->setFracPres(voidID,cid_target,fracpresVoid);
                    }

                    boxList = g_slist_next(boxList);
                }

                g_slist_free(boxList);
            }

        }

        // clean-up
        gts_bb_tree_destroy(aabbCellsTree, true);

    }

    /*----------------------------------------------------------------------------*/
    void Tools_copy_mesh_2D(kmds::Mesh* AMesh,
                            kmds::Mesh* AMesh_copy)
    {
        // TODO check that AMesh is compact
        // TODO check that AMesh is full quad
        assert(AMesh->getNbFaces() == AMesh->getNbQuads());


        const int nbNodes = AMesh->getNbNodes();

        AMesh_copy->updateNodeCapacity(nbNodes);

        AMesh_copy->addNodes(nbNodes);

        Kokkos::parallel_for(nbNodes,
                             KOKKOS_LAMBDA(const int i) {
                                 gmds::math::Point pt = AMesh->getNodeLocation(i);
                                 AMesh_copy->setNodeLocation(i, pt);
                             });



        const int nbCells = AMesh->getNbFaces();

        AMesh_copy->updateFaceCapacity(nbCells);

        AMesh_copy->addQuads(nbCells);

        Kokkos::parallel_for(nbCells,
                             KOKKOS_LAMBDA(const int i) {

                                 kmds::Face c = AMesh->getFace(i);

                                 Kokkos::View<kmds::TCellID*> nids;
                                 c.nodeIds(nids);

                                 kmds::Face c_copy = AMesh_copy->getFace(i);
                                 c_copy.setNodes(nids);
                             });

    }

    /*----------------------------------------------------------------------------*/
    void Tools_copy_mesh_3D(kmds::Mesh* AMesh,
                            kmds::Mesh* AMesh_copy)
    {
        // TODO check that AMesh is compact
        // TODO check that AMesh is full hex
        assert(AMesh->getNbRegions() == AMesh->getNbHexahedra());


        const int nbNodes = AMesh->getNbNodes();

        AMesh_copy->updateNodeCapacity(nbNodes);

        AMesh_copy->addNodes(nbNodes);

        Kokkos::parallel_for(nbNodes,
                             KOKKOS_LAMBDA(const int i) {
                                 gmds::math::Point pt = AMesh->getNodeLocation(i);
                                 AMesh_copy->setNodeLocation(i, pt);
                             });



        const int nbCells = AMesh->getNbRegions();

        AMesh_copy->updateRegionCapacity(nbCells);

        AMesh_copy->addHexahedra(nbCells);

        Kokkos::parallel_for(nbCells,
                             KOKKOS_LAMBDA(const int i) {

                                 kmds::Region c = AMesh->getRegion(i);

                                 Kokkos::View<kmds::TCellID*> nids;
                                 c.nodeIds(nids);

                                 kmds::Region c_copy = AMesh_copy->getRegion(i);
                                 c_copy.setNodes(nids);
                             });

    }

    /*----------------------------------------------------------------------------*/
    void Tools_compute_fracpres_2D(kmds::Mesh* AMesh_source,
                                   kmds::Mesh* AMesh_target,
                                   const elg3d::MaterialAssignment* Ama_source,
                                   elg3d::FracPres* Afp_target)
    {

        const kmds::TCoord minScaledJacobian_source = Tools_computeScaledJacobian_2D(AMesh_source);
        if(minScaledJacobian_source < 0.) {
            throw kmds::KException("Tools_compute_fracpres_2D : source mesh has an inverted cell");
        }
        const kmds::TCoord minScaledJacobian_target = Tools_computeScaledJacobian_2D(AMesh_target);
        if(minScaledJacobian_target < 0.) {
            throw kmds::KException("Tools_compute_fracpres_2D : target mesh has an inverted cell");
        }

        // prepare the materials
        std::map<int, std::string> matList = Ama_source->getMaterialList();
        Afp_target->setMaterialList(matList);


        // prepare the target mesh in the AABBtree
        // Axis-Aligned Bounding Box tree for the cells
        // The bounded value is the CellID
        GSList* list = NULL;


        kmds::GrowingView<kmds::TCellID> cellIDs_target("FACES", AMesh_target->getNbFaces());
        AMesh_target->getFaceIDs(&cellIDs_target);

        for(int i=0; i<cellIDs_target.getNbElems(); i++) {
            kmds::TCellID cid = cellIDs_target.get(i);

            kmds::Face c = AMesh_target->getFace(cid);

            double minXYZ[3];
            double maxXYZ[3];

            c.computeBoundingBox(minXYZ,maxXYZ);

            gpointer pointer = GINT_TO_POINTER(cid);
            GtsBBox* bbox = gts_bbox_new(
                    gts_bbox_class (),
                    pointer,
                    minXYZ[0],minXYZ[1],minXYZ[2],
                    maxXYZ[0],maxXYZ[1],maxXYZ[2]);

            list = g_slist_prepend(list,bbox);
        }
        GNode* aabbCellsTree = gts_bb_tree_new(list);

        g_slist_free(list);


        // compute discrepancy
        gpointer pointer = NULL;
        GtsBBox* bbox = gts_bbox_new(
                gts_bbox_class (),
                pointer,
                0.,0.,0.,
                0.,0.,0.);

        kmds::Variable<double>* varSumFracpres = AMesh_target->createVariable<double>(0., kmds::KMDS_FACE, "tmp_varSumFracpres");



        kmds::GrowingView<kmds::TCellID> cellIDs_source("FACES", AMesh_source->getNbFaces());
        AMesh_source->getFaceIDs(&cellIDs_source);

        for(int i=0; i<cellIDs_source.getNbElems(); i++) {
            kmds::TCellID cid_source = cellIDs_source.get(i);

            kmds::Face c_source = AMesh_source->getFace(cid_source);

            double minXYZ[3];
            double maxXYZ[3];

            c_source.computeBoundingBox(minXYZ, maxXYZ);

            gts_bbox_set(
                    bbox,
                    pointer,
                    minXYZ[0],minXYZ[1],minXYZ[2],
                    maxXYZ[0],maxXYZ[1],maxXYZ[2]);

            GSList* boxList = gts_bb_tree_overlap(aabbCellsTree,bbox);
            if(boxList == NULL) {
                throw kmds::KException("Tools_compute_fracpres_2D : a cell from source mesh does not intersect target mesh.");
            }

            int matIndex = Ama_source->getMaterial(cid_source);


            switch(c_source.computeType()) {
                case kmds::KMDS_QUAD :
                {
                    Kokkos::View<kmds::TCellID *> nids_source;
                    c_source.nodeIds(nids_source);

                    kmds::TCoord xyz_source[4][3];
                    AMesh_source->getNodeLocation(nids_source[0], xyz_source[0][0], xyz_source[0][1], xyz_source[0][2]);
                    AMesh_source->getNodeLocation(nids_source[1], xyz_source[1][0], xyz_source[1][1], xyz_source[1][2]);
                    AMesh_source->getNodeLocation(nids_source[2], xyz_source[2][0], xyz_source[2][1], xyz_source[2][2]);
                    AMesh_source->getNodeLocation(nids_source[3], xyz_source[3][0], xyz_source[3][1], xyz_source[3][2]);

                    r2d_rvec2 vertices[4];
                    vertices[0].x = xyz_source[0][0];
                    vertices[0].y = xyz_source[0][1];
                    vertices[1].x = xyz_source[1][0];
                    vertices[1].y = xyz_source[1][1];
                    vertices[2].x = xyz_source[2][0];
                    vertices[2].y = xyz_source[2][1];
                    vertices[3].x = xyz_source[3][0];
                    vertices[3].y = xyz_source[3][1];
                    r2d_int numverts = 4;
                    r2d_plane planes[4];
                    r2d_poly_faces_from_verts(planes,  vertices, numverts);


                    while (boxList != NULL) {

                        GtsBBox *box = (GtsBBox *) (boxList->data);
                        kmds::TCellID cid_target = GPOINTER_TO_INT(box->bounded);

                        kmds::Face c_target = AMesh_target->getFace(cid_target);


                        Kokkos::View<kmds::TCellID *> nids_target;
                        c_target.nodeIds(nids_target);

                        kmds::TCoord xyz_target[4][3];
                        AMesh_target->getNodeLocation(nids_target[0], xyz_target[0][0], xyz_target[0][1],
                                                      xyz_target[0][2]);
                        AMesh_target->getNodeLocation(nids_target[1], xyz_target[1][0], xyz_target[1][1],
                                                      xyz_target[1][2]);
                        AMesh_target->getNodeLocation(nids_target[2], xyz_target[2][0], xyz_target[2][1],
                                                      xyz_target[2][2]);
                        AMesh_target->getNodeLocation(nids_target[3], xyz_target[3][0], xyz_target[3][1],
                                                      xyz_target[3][2]);

                        gmds::math::Point pt0(xyz_target[0][0], xyz_target[0][1], xyz_target[0][2]);
                        gmds::math::Point pt1(xyz_target[1][0], xyz_target[1][1], xyz_target[1][2]);
                        gmds::math::Point pt2(xyz_target[2][0], xyz_target[2][1], xyz_target[2][2]);
                        gmds::math::Point pt3(xyz_target[3][0], xyz_target[3][1], xyz_target[3][2]);
                        gmds::math::Quadrilateral q(pt0, pt1, pt2, pt3);

                        const double volsurf = q.area();

                        r2d_poly poly;
                        r2d_rvec2 verts[4];

                        verts[0].x = xyz_target[0][0];
                        verts[0].y = xyz_target[0][1];
                        verts[1].x = xyz_target[1][0];
                        verts[1].y = xyz_target[1][1];
                        verts[2].x = xyz_target[2][0];
                        verts[2].y = xyz_target[2][1];
                        verts[3].x = xyz_target[3][0];
                        verts[3].y = xyz_target[3][1];

                        r2d_init_poly(&poly, verts, 4);

                        r2d_clip(&poly, planes, 3);
                        r2d_int POLY_ORDER = 2;
                        r2d_real om[R2D_NUM_MOMENTS(POLY_ORDER)];
                        r2d_reduce(&poly, om, POLY_ORDER);

                        double fracpres = Afp_target->getFracPres(matIndex, cid_target);

                        // add but do not forget to divide by cell area
                        double fracpresadd = om[0] / volsurf;

                        if (fracpresadd < 1.e-10) {
                            fracpresadd = 0.;
                        }

                        if((*varSumFracpres)[cid_target] + fracpresadd > 1.) {
                            fracpresadd = 1.- (*varSumFracpres)[cid_target];
                            (*varSumFracpres)[cid_target] = 1.;
                        } else {
                            (*varSumFracpres)[cid_target] += fracpresadd;
                        }

                        fracpres += fracpresadd;

                        Afp_target->setFracPres(matIndex, cid_target, fracpres);


                        boxList = g_slist_next(boxList);

                    }  // while (boxList != NULL)

//                    g_slist_free(boxList);

                }
                    break;
                default:
                    throw kmds::KException("Tools_compute_fracpres_2D : cell type not handled.");
                    break;
            }

        }


        // clean-up
        gts_bb_tree_destroy(aabbCellsTree, true);

        AMesh_target->deleteVariable(kmds::KMDS_FACE, "tmp_varSumFracpres");
    }

    /*----------------------------------------------------------------------------*/
    void Tools_compute_fracpres_3D(kmds::Mesh* AMesh_source,
                                   kmds::Mesh* AMesh_target,
                                   const elg3d::MaterialAssignment* Ama_source,
                                   elg3d::FracPres* Afp_target)
    {

        const kmds::TCoord minScaledJacobian_source = Tools_computeScaledJacobian_3D(AMesh_source);
        if(minScaledJacobian_source < 0.) {
            throw kmds::KException("Tools_compute_fracpres_3D : source mesh has an inverted cell");
        }
        const kmds::TCoord minScaledJacobian_target = Tools_computeScaledJacobian_3D(AMesh_target);
        if(minScaledJacobian_target < 0.) {
            throw kmds::KException("Tools_compute_fracpres_3D : target mesh has an inverted cell");
        }

        // prepare the materials
        std::map<int, std::string> matList = Ama_source->getMaterialList();
        Afp_target->setMaterialList(matList);


        // prepare the target mesh in the AABBtree
        // Axis-Aligned Bounding Box tree for the cells
        // The bounded value is the CellID
        GSList* list = NULL;


        kmds::GrowingView<kmds::TCellID> cellIDs_target("CELLS", AMesh_target->getNbRegions());
        AMesh_target->getRegionIDs(&cellIDs_target);

        for(int i=0; i<cellIDs_target.getNbElems(); i++) {
            kmds::TCellID cid = cellIDs_target.get(i);

            kmds::Region c = AMesh_target->getRegion(cid);

            double minXYZ[3];
            double maxXYZ[3];

            c.computeBoundingBox(minXYZ,maxXYZ);

            gpointer pointer = GINT_TO_POINTER(cid);
            GtsBBox* bbox = gts_bbox_new(
                    gts_bbox_class (),
                    pointer,
                    minXYZ[0],minXYZ[1],minXYZ[2],
                    maxXYZ[0],maxXYZ[1],maxXYZ[2]);

            list = g_slist_prepend(list,bbox);
        }
        GNode* aabbCellsTree = gts_bb_tree_new(list);

        g_slist_free(list);


        // compute discrepancy
        gpointer pointer = NULL;
        GtsBBox* bbox = gts_bbox_new(
                gts_bbox_class (),
                pointer,
                0.,0.,0.,
                0.,0.,0.);

        kmds::Variable<double>* varSumFracpres = AMesh_target->createVariable<double>(0., kmds::KMDS_REGION, "tmp_varSumFracpres");



        kmds::GrowingView<kmds::TCellID> cellIDs_source("CELLS", AMesh_source->getNbRegions());
        AMesh_source->getRegionIDs(&cellIDs_source);
        const kmds::TCellID nbCells_source = cellIDs_source.getNbElems();

        for(int i=0; i<nbCells_source; i++) {

            kmds::TCellID cid_source = cellIDs_source.get(i);

            if(cid_source % 10000 == 0) {
                std::cout << "cid_source " << cid_source << " of " << nbCells_source << std::endl;
            }

            kmds::Region c_source = AMesh_source->getRegion(cid_source);

            double minXYZ[3];
            double maxXYZ[3];

            c_source.computeBoundingBox(minXYZ, maxXYZ);

            gts_bbox_set(
                    bbox,
                    pointer,
                    minXYZ[0],minXYZ[1],minXYZ[2],
                    maxXYZ[0],maxXYZ[1],maxXYZ[2]);

            GSList* boxList = gts_bb_tree_overlap(aabbCellsTree,bbox);
            if(boxList == NULL) {
                throw kmds::KException("Tools_compute_fracpres_3D : a cell from source mesh does not intersect target mesh.");
            }

            int matIndex = Ama_source->getMaterial(cid_source);


            switch(c_source.computeType()) {
                case kmds::KMDS_HEX :
                {
                    Kokkos::View<kmds::TCellID *> nids_source;
                    c_source.nodeIds(nids_source);

                    kmds::TCoord xyz_source[8][3];
                    AMesh_source->getNodeLocation(nids_source[0], xyz_source[0][0], xyz_source[0][1], xyz_source[0][2]);
                    AMesh_source->getNodeLocation(nids_source[1], xyz_source[1][0], xyz_source[1][1], xyz_source[1][2]);
                    AMesh_source->getNodeLocation(nids_source[2], xyz_source[2][0], xyz_source[2][1], xyz_source[2][2]);
                    AMesh_source->getNodeLocation(nids_source[3], xyz_source[3][0], xyz_source[3][1], xyz_source[3][2]);
                    AMesh_source->getNodeLocation(nids_source[4], xyz_source[4][0], xyz_source[4][1], xyz_source[4][2]);
                    AMesh_source->getNodeLocation(nids_source[5], xyz_source[5][0], xyz_source[5][1], xyz_source[5][2]);
                    AMesh_source->getNodeLocation(nids_source[6], xyz_source[6][0], xyz_source[6][1], xyz_source[6][2]);
                    AMesh_source->getNodeLocation(nids_source[7], xyz_source[7][0], xyz_source[7][1], xyz_source[7][2]);

//                    r3d_rvec3 verts[8];
//                    verts[0].x = xyz_source[0][0];
//                    verts[0].y = xyz_source[0][1];
//                    verts[0].z = xyz_source[0][2];
//                    verts[1].x = xyz_source[1][0];
//                    verts[1].y = xyz_source[1][1];
//                    verts[1].z = xyz_source[1][2];
//                    verts[2].x = xyz_source[2][0];
//                    verts[2].y = xyz_source[2][1];
//                    verts[2].z = xyz_source[2][2];
//                    verts[3].x = xyz_source[3][0];
//                    verts[3].y = xyz_source[3][1];
//                    verts[3].z = xyz_source[3][2];
//                    verts[4].x = xyz_source[4][0];
//                    verts[4].y = xyz_source[4][1];
//                    verts[4].z = xyz_source[4][2];
//                    verts[5].x = xyz_source[5][0];
//                    verts[5].y = xyz_source[5][1];
//                    verts[5].z = xyz_source[5][2];
//                    verts[6].x = xyz_source[6][0];
//                    verts[6].y = xyz_source[6][1];
//                    verts[6].z = xyz_source[6][2];
//                    verts[7].x = xyz_source[7][0];
//                    verts[7].y = xyz_source[7][1];
//                    verts[7].z = xyz_source[7][2];

//                    // create polyhedron
//                    r3d_poly src_r3dpoly;
//                    const int num_verts = 14;
//                    const int num_faces = 24;
//
//                    r3d_int face_num_verts[num_faces];
//                    for(int i_f = 0; i_f<num_faces; i_f++) {
//                        face_num_verts[i_f] = 3;
//                    }
//
//                    r3d_int **face_vert_ids = new r3d_int *[num_faces];
//                    for (int i_f = 0; i_f < num_faces; i_f++) {
//                        face_vert_ids[i_f] = new r3d_int[face_num_verts[i_f]];
//                    }
//                    face_vert_ids[0][0] = 0;  // bottom
//                    face_vert_ids[0][1] = 1;
//                    face_vert_ids[0][2] = 8;
//                    face_vert_ids[1][0] = 1;
//                    face_vert_ids[1][1] = 2;
//                    face_vert_ids[1][2] = 8;
//                    face_vert_ids[2][0] = 2;
//                    face_vert_ids[2][1] = 3;
//                    face_vert_ids[2][2] = 8;
//                    face_vert_ids[3][0] = 3;
//                    face_vert_ids[3][1] = 0;
//                    face_vert_ids[3][2] = 8;
//                    face_vert_ids[0][0] = 0;  // top
//                    face_vert_ids[0][1] = 1;
//                    face_vert_ids[0][2] = 9;
//                    face_vert_ids[1][0] = 1;
//                    face_vert_ids[1][1] = 2;
//                    face_vert_ids[1][2] = 9;
//                    face_vert_ids[2][0] = 2;
//                    face_vert_ids[2][1] = 3;
//                    face_vert_ids[2][2] = 9;
//                    face_vert_ids[3][0] = 3;
//                    face_vert_ids[3][1] = 0;
//                    face_vert_ids[3][2] = 9;
//
//                    r3d_init_poly(&src_r3dpoly, verts, num_verts, face_vert_ids, face_num_verts,
//                num_faces);

//                    gmds::math::Point pt_center = c_source.midpoint();
//                    r3d_rvec3 vec_center;
//                    vec_center.x = pt_center.X();
//                    vec_center.y = pt_center.Y();
//                    vec_center.z = pt_center.Z();
//

                    // create polyhedron
                    r3d_poly src_cell;
                    src_cell.nverts = 8;

                    src_cell.verts[0].pos.xyz[0] = xyz_source[0][0];
                    src_cell.verts[0].pos.xyz[1] = xyz_source[0][1];
                    src_cell.verts[0].pos.xyz[2] = xyz_source[0][2];
                    src_cell.verts[1].pos.xyz[0] = xyz_source[1][0];
                    src_cell.verts[1].pos.xyz[1] = xyz_source[1][1];
                    src_cell.verts[1].pos.xyz[2] = xyz_source[1][2];
                    src_cell.verts[2].pos.xyz[0] = xyz_source[2][0];
                    src_cell.verts[2].pos.xyz[1] = xyz_source[2][1];
                    src_cell.verts[2].pos.xyz[2] = xyz_source[2][2];
                    src_cell.verts[3].pos.xyz[0] = xyz_source[3][0];
                    src_cell.verts[3].pos.xyz[1] = xyz_source[3][1];
                    src_cell.verts[3].pos.xyz[2] = xyz_source[3][2];
                    src_cell.verts[4].pos.xyz[0] = xyz_source[4][0];
                    src_cell.verts[4].pos.xyz[1] = xyz_source[4][1];
                    src_cell.verts[4].pos.xyz[2] = xyz_source[4][2];
                    src_cell.verts[5].pos.xyz[0] = xyz_source[5][0];
                    src_cell.verts[5].pos.xyz[1] = xyz_source[5][1];
                    src_cell.verts[5].pos.xyz[2] = xyz_source[5][2];
                    src_cell.verts[6].pos.xyz[0] = xyz_source[6][0];
                    src_cell.verts[6].pos.xyz[1] = xyz_source[6][1];
                    src_cell.verts[6].pos.xyz[2] = xyz_source[6][2];
                    src_cell.verts[7].pos.xyz[0] = xyz_source[7][0];
                    src_cell.verts[7].pos.xyz[1] = xyz_source[7][1];
                    src_cell.verts[7].pos.xyz[2] = xyz_source[7][2];

                    src_cell.verts[0].pnbrs[0] = 1;
                    src_cell.verts[0].pnbrs[1] = 4;
                    src_cell.verts[0].pnbrs[2] = 3;
                    src_cell.verts[1].pnbrs[0] = 2;
                    src_cell.verts[1].pnbrs[1] = 5;
                    src_cell.verts[1].pnbrs[2] = 0;
                    src_cell.verts[2].pnbrs[0] = 3;
                    src_cell.verts[2].pnbrs[1] = 6;
                    src_cell.verts[2].pnbrs[2] = 1;
                    src_cell.verts[3].pnbrs[0] = 0;
                    src_cell.verts[3].pnbrs[1] = 7;
                    src_cell.verts[3].pnbrs[2] = 2;
                    src_cell.verts[4].pnbrs[0] = 5;
                    src_cell.verts[4].pnbrs[1] = 7;
                    src_cell.verts[4].pnbrs[2] = 0;
                    src_cell.verts[5].pnbrs[0] = 6;
                    src_cell.verts[5].pnbrs[1] = 4;
                    src_cell.verts[5].pnbrs[2] = 1;
                    src_cell.verts[6].pnbrs[0] = 7;
                    src_cell.verts[6].pnbrs[1] = 5;
                    src_cell.verts[6].pnbrs[2] = 2;
                    src_cell.verts[7].pnbrs[0] = 4;
                    src_cell.verts[7].pnbrs[1] = 6;
                    src_cell.verts[7].pnbrs[2] = 3;

                    r3d_int isgood = r3d_is_good(&src_cell);
                    if(isgood == 0) {
                        std::cout<<"NOT GOOD R3D"<<std::endl;
                    }

                    while (boxList != NULL) {

                        GtsBBox *box = (GtsBBox *) (boxList->data);
                        kmds::TCellID cid_target = GPOINTER_TO_INT(box->bounded);

                        kmds::Region c_target = AMesh_target->getRegion(cid_target);


                        Kokkos::View<kmds::TCellID *> nids_target;
                        c_target.nodeIds(nids_target);

                        kmds::TCoord xyz_target[15][3];
                        AMesh_target->getNodeLocation(nids_target[0], xyz_target[0][0], xyz_target[0][1],
                                                      xyz_target[0][2]);
                        AMesh_target->getNodeLocation(nids_target[1], xyz_target[1][0], xyz_target[1][1],
                                                      xyz_target[1][2]);
                        AMesh_target->getNodeLocation(nids_target[2], xyz_target[2][0], xyz_target[2][1],
                                                      xyz_target[2][2]);
                        AMesh_target->getNodeLocation(nids_target[3], xyz_target[3][0], xyz_target[3][1],
                                                      xyz_target[3][2]);
                        AMesh_target->getNodeLocation(nids_target[4], xyz_target[4][0], xyz_target[4][1],
                                                      xyz_target[4][2]);
                        AMesh_target->getNodeLocation(nids_target[5], xyz_target[5][0], xyz_target[5][1],
                                                      xyz_target[5][2]);
                        AMesh_target->getNodeLocation(nids_target[6], xyz_target[6][0], xyz_target[6][1],
                                                      xyz_target[6][2]);
                        AMesh_target->getNodeLocation(nids_target[7], xyz_target[7][0], xyz_target[7][1],
                                                      xyz_target[7][2]);

                        gmds::math::Point pt0(xyz_target[0][0], xyz_target[0][1], xyz_target[0][2]);
                        gmds::math::Point pt1(xyz_target[1][0], xyz_target[1][1], xyz_target[1][2]);
                        gmds::math::Point pt2(xyz_target[2][0], xyz_target[2][1], xyz_target[2][2]);
                        gmds::math::Point pt3(xyz_target[3][0], xyz_target[3][1], xyz_target[3][2]);
                        gmds::math::Point pt4(xyz_target[4][0], xyz_target[4][1], xyz_target[4][2]);
                        gmds::math::Point pt5(xyz_target[5][0], xyz_target[5][1], xyz_target[5][2]);
                        gmds::math::Point pt6(xyz_target[6][0], xyz_target[6][1], xyz_target[6][2]);
                        gmds::math::Point pt7(xyz_target[7][0], xyz_target[7][1], xyz_target[7][2]);

                        gmds::math::Hexahedron h(pt0, pt1, pt2, pt3, pt4, pt5, pt6, pt7);
                        const double volsurf = h.getVolume();

                        r3d_rvec3 verts2[15];
                        verts2[0].x = xyz_target[0][0];
                        verts2[0].y = xyz_target[0][1];
                        verts2[0].z = xyz_target[0][2];
                        verts2[1].x = xyz_target[1][0];
                        verts2[1].y = xyz_target[1][1];
                        verts2[1].z = xyz_target[1][2];
                        verts2[2].x = xyz_target[2][0];
                        verts2[2].y = xyz_target[2][1];
                        verts2[2].z = xyz_target[2][2];
                        verts2[3].x = xyz_target[3][0];
                        verts2[3].y = xyz_target[3][1];
                        verts2[3].z = xyz_target[3][2];
                        verts2[4].x = xyz_target[4][0];
                        verts2[4].y = xyz_target[4][1];
                        verts2[4].z = xyz_target[4][2];
                        verts2[5].x = xyz_target[5][0];
                        verts2[5].y = xyz_target[5][1];
                        verts2[5].z = xyz_target[5][2];
                        verts2[6].x = xyz_target[6][0];
                        verts2[6].y = xyz_target[6][1];
                        verts2[6].z = xyz_target[6][2];
                        verts2[7].x = xyz_target[7][0];
                        verts2[7].y = xyz_target[7][1];
                        verts2[7].z = xyz_target[7][2];

                        gmds::math::Point pt_bottom = 0.25 * (pt0 + pt1 + pt2 + pt3);
                        gmds::math::Point pt_top    = 0.25 * (pt4 + pt5 + pt6 + pt7);
                        gmds::math::Point pt_left   = 0.25 * (pt0 + pt3 + pt7 + pt4);
                        gmds::math::Point pt_right  = 0.25 * (pt1 + pt2 + pt6 + pt5);
                        gmds::math::Point pt_front  = 0.25 * (pt0 + pt1 + pt5 + pt4);
                        gmds::math::Point pt_back   = 0.25 * (pt3 + pt2 + pt6 + pt7);
                        verts2[8].x  = pt_bottom.X();
                        verts2[8].y  = pt_bottom.Y();
                        verts2[8].z  = pt_bottom.Z();
                        verts2[9].x  = pt_top.X();
                        verts2[9].y  = pt_top.Y();
                        verts2[9].z  = pt_top.Z();
                        verts2[10].x = pt_left.X();
                        verts2[10].y = pt_left.Y();
                        verts2[10].z = pt_left.Z();
                        verts2[11].x = pt_right.X();
                        verts2[11].y = pt_right.Y();
                        verts2[11].z = pt_right.Z();
                        verts2[12].x = pt_front.X();
                        verts2[12].y = pt_front.Y();
                        verts2[12].z = pt_front.Z();
                        verts2[13].x = pt_back.X();
                        verts2[13].y = pt_back.Y();
                        verts2[13].z = pt_back.Z();

                        gmds::math::Point pt_center = c_target.midpoint();
                        verts2[14].x = pt_center.X();
                        verts2[14].y = pt_center.Y();
                        verts2[14].z = pt_center.Z();

                        r3d_rvec3 verts_tet[4];
                        r3d_plane faces[4];
                        r3d_poly src_cell_copy;
                        const int POLY_ORDER = 1;
                        r3d_real intersectVol = 0.;
                        r3d_real om[R3D_NUM_MOMENTS(POLY_ORDER)];

                        // bottom
                        verts_tet[0] = verts2[0];
                        verts_tet[1] = verts2[1];
                        verts_tet[2] = verts2[8];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[1];
                        verts_tet[1] = verts2[2];
                        verts_tet[2] = verts2[8];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[2];
                        verts_tet[1] = verts2[3];
                        verts_tet[2] = verts2[8];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[3];
                        verts_tet[1] = verts2[0];
                        verts_tet[2] = verts2[8];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];

                        // top
                        verts_tet[0] = verts2[4];
                        verts_tet[1] = verts2[7];
                        verts_tet[2] = verts2[9];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[7];
                        verts_tet[1] = verts2[6];
                        verts_tet[2] = verts2[9];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[6];
                        verts_tet[1] = verts2[5];
                        verts_tet[2] = verts2[9];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[5];
                        verts_tet[1] = verts2[4];
                        verts_tet[2] = verts2[9];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];

                        // left
                        verts_tet[0] = verts2[0];
                        verts_tet[1] = verts2[3];
                        verts_tet[2] = verts2[10];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[3];
                        verts_tet[1] = verts2[7];
                        verts_tet[2] = verts2[10];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[7];
                        verts_tet[1] = verts2[4];
                        verts_tet[2] = verts2[10];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[4];
                        verts_tet[1] = verts2[0];
                        verts_tet[2] = verts2[10];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];

                        // right
                        verts_tet[0] = verts2[2];
                        verts_tet[1] = verts2[1];
                        verts_tet[2] = verts2[11];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[1];
                        verts_tet[1] = verts2[5];
                        verts_tet[2] = verts2[11];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[5];
                        verts_tet[1] = verts2[6];
                        verts_tet[2] = verts2[11];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[6];
                        verts_tet[1] = verts2[2];
                        verts_tet[2] = verts2[11];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];

                        // front
                        verts_tet[0] = verts2[1];
                        verts_tet[1] = verts2[0];
                        verts_tet[2] = verts2[12];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[0];
                        verts_tet[1] = verts2[4];
                        verts_tet[2] = verts2[12];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[4];
                        verts_tet[1] = verts2[5];
                        verts_tet[2] = verts2[12];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[5];
                        verts_tet[1] = verts2[1];
                        verts_tet[2] = verts2[12];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];

                        // back
                        verts_tet[0] = verts2[3];
                        verts_tet[1] = verts2[2];
                        verts_tet[2] = verts2[13];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[2];
                        verts_tet[1] = verts2[6];
                        verts_tet[2] = verts2[13];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[6];
                        verts_tet[1] = verts2[7];
                        verts_tet[2] = verts2[13];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[7];
                        verts_tet[1] = verts2[3];
                        verts_tet[2] = verts2[13];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];

                        double fracpres = Afp_target->getFracPres(matIndex, cid_target);

                        // add but do not forget to divide by cell area
                        double fracpresadd = intersectVol / volsurf;

                        if (fracpresadd < 1.e-10) {
                            fracpresadd = 0.;
                        }

                        if((*varSumFracpres)[cid_target] + fracpresadd > 1.) {
                            fracpresadd = 1.- (*varSumFracpres)[cid_target];
                            (*varSumFracpres)[cid_target] = 1.;
                        } else {
                            (*varSumFracpres)[cid_target] += fracpresadd;
                        }

                        fracpres += fracpresadd;

                        Afp_target->setFracPres(matIndex, cid_target, fracpres);


                        boxList = g_slist_next(boxList);

                    }  // while (boxList != NULL)

//                    g_slist_free(boxList);

                }
                    break;
                default:
                    throw kmds::KException("Tools_compute_fracpres_3D : cell type not handled.");
                    break;
            }

        }


        // clean-up
        gts_bb_tree_destroy(aabbCellsTree, true);

        AMesh_target->deleteVariable(kmds::KMDS_REGION, "tmp_varSumFracpres");
    }

    /*----------------------------------------------------------------------------*/
    double
    Tools_compute_discrepancy_2D(kmds::Mesh* AMesh,
                                 const elg3d::FracPres* Afp_ref,
                                 const elg3d::FracPres* Afp_new,
                                 kmds::Variable<double>* AVarDiscrepancy)
    {

        const int nbMat = Afp_ref->getNbMaterials();

        double discrepancy_glob = 0.;

        kmds::GrowingView<kmds::TCellID> cellIDs("FACES", AMesh->getNbFaces());
        AMesh->getFaceIDs(&cellIDs);

        Kokkos::parallel_reduce(cellIDs.getNbElems(),
                             KOKKOS_LAMBDA(const int i, double& sum) {

                                 kmds::TCellID cid = cellIDs.get(i);
                                 kmds::Face c = AMesh->getFace(i);

                                 double vol = c.surfvol();

                                 double discrepancy_loc = 0.;
                                 for(int imat=0; imat<nbMat; imat++) {

                                     discrepancy_loc += std::fabs(Afp_ref->getFracPres(imat, cid) - Afp_new->getFracPres(imat, cid));
                                 }

                                 discrepancy_loc *= vol;

                                 (*AVarDiscrepancy)[cid] = discrepancy_loc;
                                 sum += discrepancy_loc;
                             },
                             discrepancy_glob);

    }

    /*----------------------------------------------------------------------------*/
    double
    Tools_compute_discrepancy_3D(kmds::Mesh* AMesh,
                                 const elg3d::FracPres* Afp_ref,
                                 const elg3d::FracPres* Afp_new,
                                 kmds::Variable<double>* AVarDiscrepancy)
    {

        const int nbMat = Afp_ref->getNbMaterials();

        double discrepancy_glob = 0.;

        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", AMesh->getNbRegions());
        AMesh->getRegionIDs(&cellIDs);

        Kokkos::parallel_reduce(cellIDs.getNbElems(),
                                KOKKOS_LAMBDA(const int i, double& sum) {

                                    kmds::TCellID cid = cellIDs.get(i);
                                    kmds::Region c = AMesh->getRegion(i);

                                    double vol = c.surfvol();

                                    double discrepancy_loc = 0.;
                                    for(int imat=0; imat<nbMat; imat++) {

                                        discrepancy_loc += std::fabs(Afp_ref->getFracPres(imat, cid) - Afp_new->getFracPres(imat, cid));
                                    }

                                    discrepancy_loc *= vol;

                                    (*AVarDiscrepancy)[cid] = discrepancy_loc;
                                    sum += discrepancy_loc;
                                },
                                discrepancy_glob);

    }

    /*----------------------------------------------------------------------------*/
    void Tools_compute_fracpres_source_and_submesh_2D(const kmds::Mesh* AMesh_source,
                                                      const elg3d::FracPres* Afp_source,
                                                      const kmds::Variable<kmds::TCellID>* AVarOldCells2firstSubCells,
                                                      const kmds::Mesh* AMesh_submesh,
                                                      const elg3d::MaterialAssignment* Ama_submesh,
                                                      kmds::Mesh* AMesh_target,
                                                      elg3d::FracPres* Afp_target)
    {

        const kmds::TCoord minScaledJacobian_source = Tools_computeScaledJacobian_2D(AMesh_source);
        if(minScaledJacobian_source < 0.) {
            throw kmds::KException("Tools_compute_fracpres_2D : source mesh has an inverted cell");
        }
        const kmds::TCoord minScaledJacobian_target = Tools_computeScaledJacobian_2D(AMesh_target);
        if(minScaledJacobian_target < 0.) {
            throw kmds::KException("Tools_compute_fracpres_2D : target mesh has an inverted cell");
        }

        // prepare the materials
        std::map<int, std::string> matList = Afp_source->getMaterialList();
        Afp_target->setMaterialList(matList);


        // prepare the target mesh in the AABBtree
        // Axis-Aligned Bounding Box tree for the cells
        // The bounded value is the CellID
        GSList* list = NULL;


        kmds::GrowingView<kmds::TCellID> cellIDs_target("FACES", AMesh_target->getNbFaces());
        AMesh_target->getFaceIDs(&cellIDs_target);

        for(int i=0; i<cellIDs_target.getNbElems(); i++) {
            kmds::TCellID cid = cellIDs_target.get(i);

            kmds::Face c = AMesh_target->getFace(cid);

            double minXYZ[3];
            double maxXYZ[3];

            c.computeBoundingBox(minXYZ,maxXYZ);

            gpointer pointer = GINT_TO_POINTER(cid);
            GtsBBox* bbox = gts_bbox_new(
                    gts_bbox_class (),
                    pointer,
                    minXYZ[0],minXYZ[1],minXYZ[2],
                    maxXYZ[0],maxXYZ[1],maxXYZ[2]);

            list = g_slist_prepend(list,bbox);
        }
        GNode* aabbCellsTree = gts_bb_tree_new(list);

        g_slist_free(list);


        // compute discrepancy
        gpointer pointer = NULL;
        GtsBBox* bbox = gts_bbox_new(
                gts_bbox_class (),
                pointer,
                0.,0.,0.,
                0.,0.,0.);

        kmds::Variable<double>* varSumFracpres = AMesh_target->createVariable<double>(0., kmds::KMDS_FACE, "tmp_varSumFracpres");


        kmds::GrowingView<kmds::TCellID> cellIDs_source("FACES", AMesh_source->getNbFaces());
        AMesh_source->getFaceIDs(&cellIDs_source);

        for(int i=0; i<cellIDs_source.getNbElems(); i++) {
            kmds::TCellID cid_source = cellIDs_source.get(i);

            // only project source cells not in the submesh
            if((*AVarOldCells2firstSubCells)[cid_source] == kmds::NullID) {

                kmds::Face c_source = AMesh_source->getFace(cid_source);

                double minXYZ[3];
                double maxXYZ[3];

                c_source.computeBoundingBox(minXYZ, maxXYZ);

                gts_bbox_set(
                        bbox,
                        pointer,
                        minXYZ[0], minXYZ[1], minXYZ[2],
                        maxXYZ[0], maxXYZ[1], maxXYZ[2]);

                GSList *boxList = gts_bb_tree_overlap(aabbCellsTree, bbox);
                if (boxList == NULL) {
                    throw kmds::KException(
                            "Tools_compute_fracpres_2D : a cell from source mesh does not intersect target mesh.");
                }

                int matIndex = Afp_source->getMaxMatFracPresIndex(cid_source);

                switch (c_source.computeType()) {
                    case kmds::KMDS_QUAD : {
                        Kokkos::View<kmds::TCellID *> nids_source;
                        c_source.nodeIds(nids_source);

                        kmds::TCoord xyz_source[4][3];
                        AMesh_source->getNodeLocation(nids_source[0], xyz_source[0][0], xyz_source[0][1],
                                                      xyz_source[0][2]);
                        AMesh_source->getNodeLocation(nids_source[1], xyz_source[1][0], xyz_source[1][1],
                                                      xyz_source[1][2]);
                        AMesh_source->getNodeLocation(nids_source[2], xyz_source[2][0], xyz_source[2][1],
                                                      xyz_source[2][2]);
                        AMesh_source->getNodeLocation(nids_source[3], xyz_source[3][0], xyz_source[3][1],
                                                      xyz_source[3][2]);

                        r2d_rvec2 vertices[4];
                        vertices[0].x = xyz_source[0][0];
                        vertices[0].y = xyz_source[0][1];
                        vertices[1].x = xyz_source[1][0];
                        vertices[1].y = xyz_source[1][1];
                        vertices[2].x = xyz_source[2][0];
                        vertices[2].y = xyz_source[2][1];
                        vertices[3].x = xyz_source[3][0];
                        vertices[3].y = xyz_source[3][1];
                        r2d_int numverts = 4;
                        r2d_plane planes[4];
                        r2d_poly_faces_from_verts(planes, vertices, numverts);


                        while (boxList != NULL) {

                            GtsBBox *box = (GtsBBox *) (boxList->data);
                            kmds::TCellID cid_target = GPOINTER_TO_INT(box->bounded);

                            kmds::Face c_target = AMesh_target->getFace(cid_target);


                            Kokkos::View<kmds::TCellID *> nids_target;
                            c_target.nodeIds(nids_target);

                            kmds::TCoord xyz_target[4][3];
                            AMesh_target->getNodeLocation(nids_target[0], xyz_target[0][0], xyz_target[0][1],
                                                          xyz_target[0][2]);
                            AMesh_target->getNodeLocation(nids_target[1], xyz_target[1][0], xyz_target[1][1],
                                                          xyz_target[1][2]);
                            AMesh_target->getNodeLocation(nids_target[2], xyz_target[2][0], xyz_target[2][1],
                                                          xyz_target[2][2]);
                            AMesh_target->getNodeLocation(nids_target[3], xyz_target[3][0], xyz_target[3][1],
                                                          xyz_target[3][2]);

                            gmds::math::Point pt0(xyz_target[0][0], xyz_target[0][1], xyz_target[0][2]);
                            gmds::math::Point pt1(xyz_target[1][0], xyz_target[1][1], xyz_target[1][2]);
                            gmds::math::Point pt2(xyz_target[2][0], xyz_target[2][1], xyz_target[2][2]);
                            gmds::math::Point pt3(xyz_target[3][0], xyz_target[3][1], xyz_target[3][2]);
                            gmds::math::Quadrilateral q(pt0, pt1, pt2, pt3);

                            const double volsurf = q.area();

                            r2d_poly poly;
                            r2d_rvec2 verts[4];

                            verts[0].x = xyz_target[0][0];
                            verts[0].y = xyz_target[0][1];
                            verts[1].x = xyz_target[1][0];
                            verts[1].y = xyz_target[1][1];
                            verts[2].x = xyz_target[2][0];
                            verts[2].y = xyz_target[2][1];
                            verts[3].x = xyz_target[3][0];
                            verts[3].y = xyz_target[3][1];

                            r2d_init_poly(&poly, verts, 4);

                            r2d_clip(&poly, planes, 3);
                            r2d_int POLY_ORDER = 2;
                            r2d_real om[R2D_NUM_MOMENTS(POLY_ORDER)];
                            r2d_reduce(&poly, om, POLY_ORDER);

                            double fracpres = Afp_target->getFracPres(matIndex, cid_target);

                            // add but do not forget to divide by cell area
                            double fracpresadd = om[0] / volsurf;

                            if (fracpresadd < 1.e-10) {
                                fracpresadd = 0.;
                            }

                            if ((*varSumFracpres)[cid_target] + fracpresadd > 1.) {
                                fracpresadd = 1. - (*varSumFracpres)[cid_target];
                                (*varSumFracpres)[cid_target] = 1.;
                            } else {
                                (*varSumFracpres)[cid_target] += fracpresadd;
                            }

                            fracpres += fracpresadd;

                            Afp_target->setFracPres(matIndex, cid_target, fracpres);


                            boxList = g_slist_next(boxList);

                        }  // while (boxList != NULL)

//                    g_slist_free(boxList);

                    }
                        break;
                    default:
                        throw kmds::KException("Tools_compute_fracpres_2D : cell type not handled.");
                        break;
                }

            }  // if((*AVarOldCells2firstSubCells)[cid_source] == kmds::NullID)
        }

        kmds::GrowingView<kmds::TCellID> cellIDs_submesh("FACES", AMesh_submesh->getNbFaces());
        AMesh_submesh->getFaceIDs(&cellIDs_submesh);

        for(int i=0; i<cellIDs_submesh.getNbElems(); i++) {
            kmds::TCellID cid_source = cellIDs_submesh.get(i);


            kmds::Face c_source = AMesh_submesh->getFace(cid_source);

            double minXYZ[3];
            double maxXYZ[3];

            c_source.computeBoundingBox(minXYZ, maxXYZ);

            gts_bbox_set(
                    bbox,
                    pointer,
                    minXYZ[0], minXYZ[1], minXYZ[2],
                    maxXYZ[0], maxXYZ[1], maxXYZ[2]);

            GSList *boxList = gts_bb_tree_overlap(aabbCellsTree, bbox);
            if (boxList == NULL) {
                throw kmds::KException(
                        "Tools_compute_fracpres_2D : a cell from source mesh does not intersect target mesh.");
            }

            int matIndex = Ama_submesh->getMaterial(cid_source);

            switch (c_source.computeType()) {
                case kmds::KMDS_QUAD : {
                    Kokkos::View<kmds::TCellID *> nids_source;
                    c_source.nodeIds(nids_source);

                    kmds::TCoord xyz_source[4][3];
                    AMesh_submesh->getNodeLocation(nids_source[0], xyz_source[0][0], xyz_source[0][1],
                                                  xyz_source[0][2]);
                    AMesh_submesh->getNodeLocation(nids_source[1], xyz_source[1][0], xyz_source[1][1],
                                                  xyz_source[1][2]);
                    AMesh_submesh->getNodeLocation(nids_source[2], xyz_source[2][0], xyz_source[2][1],
                                                  xyz_source[2][2]);
                    AMesh_submesh->getNodeLocation(nids_source[3], xyz_source[3][0], xyz_source[3][1],
                                                  xyz_source[3][2]);

                    r2d_rvec2 vertices[4];
                    vertices[0].x = xyz_source[0][0];
                    vertices[0].y = xyz_source[0][1];
                    vertices[1].x = xyz_source[1][0];
                    vertices[1].y = xyz_source[1][1];
                    vertices[2].x = xyz_source[2][0];
                    vertices[2].y = xyz_source[2][1];
                    vertices[3].x = xyz_source[3][0];
                    vertices[3].y = xyz_source[3][1];
                    r2d_int numverts = 4;
                    r2d_plane planes[4];
                    r2d_poly_faces_from_verts(planes, vertices, numverts);


                    while (boxList != NULL) {

                        GtsBBox *box = (GtsBBox *) (boxList->data);
                        kmds::TCellID cid_target = GPOINTER_TO_INT(box->bounded);

                        kmds::Face c_target = AMesh_target->getFace(cid_target);


                        Kokkos::View<kmds::TCellID *> nids_target;
                        c_target.nodeIds(nids_target);

                        kmds::TCoord xyz_target[4][3];
                        AMesh_target->getNodeLocation(nids_target[0], xyz_target[0][0], xyz_target[0][1],
                                                      xyz_target[0][2]);
                        AMesh_target->getNodeLocation(nids_target[1], xyz_target[1][0], xyz_target[1][1],
                                                      xyz_target[1][2]);
                        AMesh_target->getNodeLocation(nids_target[2], xyz_target[2][0], xyz_target[2][1],
                                                      xyz_target[2][2]);
                        AMesh_target->getNodeLocation(nids_target[3], xyz_target[3][0], xyz_target[3][1],
                                                      xyz_target[3][2]);

                        gmds::math::Point pt0(xyz_target[0][0], xyz_target[0][1], xyz_target[0][2]);
                        gmds::math::Point pt1(xyz_target[1][0], xyz_target[1][1], xyz_target[1][2]);
                        gmds::math::Point pt2(xyz_target[2][0], xyz_target[2][1], xyz_target[2][2]);
                        gmds::math::Point pt3(xyz_target[3][0], xyz_target[3][1], xyz_target[3][2]);
                        gmds::math::Quadrilateral q(pt0, pt1, pt2, pt3);

                        const double volsurf = q.area();

                        r2d_poly poly;
                        r2d_rvec2 verts[4];

                        verts[0].x = xyz_target[0][0];
                        verts[0].y = xyz_target[0][1];
                        verts[1].x = xyz_target[1][0];
                        verts[1].y = xyz_target[1][1];
                        verts[2].x = xyz_target[2][0];
                        verts[2].y = xyz_target[2][1];
                        verts[3].x = xyz_target[3][0];
                        verts[3].y = xyz_target[3][1];

                        r2d_init_poly(&poly, verts, 4);

                        r2d_clip(&poly, planes, 3);
                        r2d_int POLY_ORDER = 2;
                        r2d_real om[R2D_NUM_MOMENTS(POLY_ORDER)];
                        r2d_reduce(&poly, om, POLY_ORDER);

                        double fracpres = Afp_target->getFracPres(matIndex, cid_target);

                        // add but do not forget to divide by cell area
                        double fracpresadd = om[0] / volsurf;

                        if (fracpresadd < 1.e-10) {
                            fracpresadd = 0.;
                        }

                        if ((*varSumFracpres)[cid_target] + fracpresadd > 1.) {
                            fracpresadd = 1. - (*varSumFracpres)[cid_target];
                            (*varSumFracpres)[cid_target] = 1.;
                        } else {
                            (*varSumFracpres)[cid_target] += fracpresadd;
                        }

                        fracpres += fracpresadd;

                        Afp_target->setFracPres(matIndex, cid_target, fracpres);


                        boxList = g_slist_next(boxList);

                    }  // while (boxList != NULL)

//                    g_slist_free(boxList);

                }
                    break;
                default:
                    throw kmds::KException("Tools_compute_fracpres_2D : cell type not handled.");
                    break;
            }

        }


        // clean-up
        gts_bb_tree_destroy(aabbCellsTree, true);

        AMesh_target->deleteVariable(kmds::KMDS_FACE, "tmp_varSumFracpres");
    }

    /*----------------------------------------------------------------------------*/
    void Tools_compute_fracpres_source_and_submesh_3D(const kmds::Mesh* AMesh_source,
                                                      const elg3d::FracPres* Afp_source,
                                                      const kmds::Variable<kmds::TCellID>* AVarOldCells2firstSubCells,
                                                      const kmds::Mesh* AMesh_submesh,
                                                      const elg3d::MaterialAssignment* Ama_submesh,
                                                      kmds::Mesh* AMesh_target,
                                                      elg3d::FracPres* Afp_target)
    {

        const kmds::TCoord minScaledJacobian_source = Tools_computeScaledJacobian_3D(AMesh_source);
        if(minScaledJacobian_source < 0.) {
            throw kmds::KException("Tools_compute_fracpres_3D : source mesh has an inverted cell");
        }
        const kmds::TCoord minScaledJacobian_target = Tools_computeScaledJacobian_3D(AMesh_target);
        if(minScaledJacobian_target < 0.) {
            throw kmds::KException("Tools_compute_fracpres_3D : target mesh has an inverted cell");
        }

        // prepare the materials
        std::map<int, std::string> matList = Afp_source->getMaterialList();
        Afp_target->setMaterialList(matList);


        // prepare the target mesh in the AABBtree
        // Axis-Aligned Bounding Box tree for the cells
        // The bounded value is the CellID
        GSList* list = NULL;


        kmds::GrowingView<kmds::TCellID> cellIDs_target("CELLS", AMesh_target->getNbRegions());
        AMesh_target->getRegionIDs(&cellIDs_target);

        for(int i=0; i<cellIDs_target.getNbElems(); i++) {
            kmds::TCellID cid = cellIDs_target.get(i);

            kmds::Region c = AMesh_target->getRegion(cid);

            double minXYZ[3];
            double maxXYZ[3];

            c.computeBoundingBox(minXYZ,maxXYZ);

            gpointer pointer = GINT_TO_POINTER(cid);
            GtsBBox* bbox = gts_bbox_new(
                    gts_bbox_class (),
                    pointer,
                    minXYZ[0],minXYZ[1],minXYZ[2],
                    maxXYZ[0],maxXYZ[1],maxXYZ[2]);

            list = g_slist_prepend(list,bbox);
        }
        GNode* aabbCellsTree = gts_bb_tree_new(list);

        g_slist_free(list);


        // compute discrepancy
        gpointer pointer = NULL;
        GtsBBox* bbox = gts_bbox_new(
                gts_bbox_class (),
                pointer,
                0.,0.,0.,
                0.,0.,0.);

        kmds::Variable<double>* varSumFracpres = AMesh_target->createVariable<double>(0., kmds::KMDS_REGION, "tmp_varSumFracpres");



        kmds::GrowingView<kmds::TCellID> cellIDs_source("CELLS", AMesh_source->getNbRegions());
        AMesh_source->getRegionIDs(&cellIDs_source);
        const kmds::TCellID nbCells_source = cellIDs_source.getNbElems();

        for(int i=0; i<nbCells_source; i++) {

            kmds::TCellID cid_source = cellIDs_source.get(i);

            std::cout<<"cid_source "<<cid_source<<" of "<<nbCells_source<<std::endl;

            if((*AVarOldCells2firstSubCells)[cid_source] == kmds::NullID) {

                kmds::Region c_source = AMesh_source->getRegion(cid_source);

                double minXYZ[3];
                double maxXYZ[3];

                c_source.computeBoundingBox(minXYZ, maxXYZ);

                gts_bbox_set(
                        bbox,
                        pointer,
                        minXYZ[0], minXYZ[1], minXYZ[2],
                        maxXYZ[0], maxXYZ[1], maxXYZ[2]);

                GSList *boxList = gts_bb_tree_overlap(aabbCellsTree, bbox);
                if (boxList == NULL) {
                    throw kmds::KException(
                            "Tools_compute_fracpres_3D : a cell from source mesh does not intersect target mesh.");
                }

                int matIndex = Afp_source->getMaxMatFracPresIndex(cid_source);


                switch (c_source.computeType()) {
                    case kmds::KMDS_HEX : {
                        Kokkos::View<kmds::TCellID *> nids_source;
                        c_source.nodeIds(nids_source);

                        kmds::TCoord xyz_source[8][3];
                        AMesh_source->getNodeLocation(nids_source[0], xyz_source[0][0], xyz_source[0][1],
                                                      xyz_source[0][2]);
                        AMesh_source->getNodeLocation(nids_source[1], xyz_source[1][0], xyz_source[1][1],
                                                      xyz_source[1][2]);
                        AMesh_source->getNodeLocation(nids_source[2], xyz_source[2][0], xyz_source[2][1],
                                                      xyz_source[2][2]);
                        AMesh_source->getNodeLocation(nids_source[3], xyz_source[3][0], xyz_source[3][1],
                                                      xyz_source[3][2]);
                        AMesh_source->getNodeLocation(nids_source[4], xyz_source[4][0], xyz_source[4][1],
                                                      xyz_source[4][2]);
                        AMesh_source->getNodeLocation(nids_source[5], xyz_source[5][0], xyz_source[5][1],
                                                      xyz_source[5][2]);
                        AMesh_source->getNodeLocation(nids_source[6], xyz_source[6][0], xyz_source[6][1],
                                                      xyz_source[6][2]);
                        AMesh_source->getNodeLocation(nids_source[7], xyz_source[7][0], xyz_source[7][1],
                                                      xyz_source[7][2]);

//                    r3d_rvec3 verts[8];
//                    verts[0].x = xyz_source[0][0];
//                    verts[0].y = xyz_source[0][1];
//                    verts[0].z = xyz_source[0][2];
//                    verts[1].x = xyz_source[1][0];
//                    verts[1].y = xyz_source[1][1];
//                    verts[1].z = xyz_source[1][2];
//                    verts[2].x = xyz_source[2][0];
//                    verts[2].y = xyz_source[2][1];
//                    verts[2].z = xyz_source[2][2];
//                    verts[3].x = xyz_source[3][0];
//                    verts[3].y = xyz_source[3][1];
//                    verts[3].z = xyz_source[3][2];
//                    verts[4].x = xyz_source[4][0];
//                    verts[4].y = xyz_source[4][1];
//                    verts[4].z = xyz_source[4][2];
//                    verts[5].x = xyz_source[5][0];
//                    verts[5].y = xyz_source[5][1];
//                    verts[5].z = xyz_source[5][2];
//                    verts[6].x = xyz_source[6][0];
//                    verts[6].y = xyz_source[6][1];
//                    verts[6].z = xyz_source[6][2];
//                    verts[7].x = xyz_source[7][0];
//                    verts[7].y = xyz_source[7][1];
//                    verts[7].z = xyz_source[7][2];

//                    // create polyhedron
//                    r3d_poly src_r3dpoly;
//                    const int num_verts = 14;
//                    const int num_faces = 24;
//
//                    r3d_int face_num_verts[num_faces];
//                    for(int i_f = 0; i_f<num_faces; i_f++) {
//                        face_num_verts[i_f] = 3;
//                    }
//
//                    r3d_int **face_vert_ids = new r3d_int *[num_faces];
//                    for (int i_f = 0; i_f < num_faces; i_f++) {
//                        face_vert_ids[i_f] = new r3d_int[face_num_verts[i_f]];
//                    }
//                    face_vert_ids[0][0] = 0;  // bottom
//                    face_vert_ids[0][1] = 1;
//                    face_vert_ids[0][2] = 8;
//                    face_vert_ids[1][0] = 1;
//                    face_vert_ids[1][1] = 2;
//                    face_vert_ids[1][2] = 8;
//                    face_vert_ids[2][0] = 2;
//                    face_vert_ids[2][1] = 3;
//                    face_vert_ids[2][2] = 8;
//                    face_vert_ids[3][0] = 3;
//                    face_vert_ids[3][1] = 0;
//                    face_vert_ids[3][2] = 8;
//                    face_vert_ids[0][0] = 0;  // top
//                    face_vert_ids[0][1] = 1;
//                    face_vert_ids[0][2] = 9;
//                    face_vert_ids[1][0] = 1;
//                    face_vert_ids[1][1] = 2;
//                    face_vert_ids[1][2] = 9;
//                    face_vert_ids[2][0] = 2;
//                    face_vert_ids[2][1] = 3;
//                    face_vert_ids[2][2] = 9;
//                    face_vert_ids[3][0] = 3;
//                    face_vert_ids[3][1] = 0;
//                    face_vert_ids[3][2] = 9;
//
//                    r3d_init_poly(&src_r3dpoly, verts, num_verts, face_vert_ids, face_num_verts,
//                num_faces);

//                    gmds::math::Point pt_center = c_source.midpoint();
//                    r3d_rvec3 vec_center;
//                    vec_center.x = pt_center.X();
//                    vec_center.y = pt_center.Y();
//                    vec_center.z = pt_center.Z();
//

                        // create polyhedron
                        r3d_poly src_cell;
                        src_cell.nverts = 8;

                        src_cell.verts[0].pos.xyz[0] = xyz_source[0][0];
                        src_cell.verts[0].pos.xyz[1] = xyz_source[0][1];
                        src_cell.verts[0].pos.xyz[2] = xyz_source[0][2];
                        src_cell.verts[1].pos.xyz[0] = xyz_source[1][0];
                        src_cell.verts[1].pos.xyz[1] = xyz_source[1][1];
                        src_cell.verts[1].pos.xyz[2] = xyz_source[1][2];
                        src_cell.verts[2].pos.xyz[0] = xyz_source[2][0];
                        src_cell.verts[2].pos.xyz[1] = xyz_source[2][1];
                        src_cell.verts[2].pos.xyz[2] = xyz_source[2][2];
                        src_cell.verts[3].pos.xyz[0] = xyz_source[3][0];
                        src_cell.verts[3].pos.xyz[1] = xyz_source[3][1];
                        src_cell.verts[3].pos.xyz[2] = xyz_source[3][2];
                        src_cell.verts[4].pos.xyz[0] = xyz_source[4][0];
                        src_cell.verts[4].pos.xyz[1] = xyz_source[4][1];
                        src_cell.verts[4].pos.xyz[2] = xyz_source[4][2];
                        src_cell.verts[5].pos.xyz[0] = xyz_source[5][0];
                        src_cell.verts[5].pos.xyz[1] = xyz_source[5][1];
                        src_cell.verts[5].pos.xyz[2] = xyz_source[5][2];
                        src_cell.verts[6].pos.xyz[0] = xyz_source[6][0];
                        src_cell.verts[6].pos.xyz[1] = xyz_source[6][1];
                        src_cell.verts[6].pos.xyz[2] = xyz_source[6][2];
                        src_cell.verts[7].pos.xyz[0] = xyz_source[7][0];
                        src_cell.verts[7].pos.xyz[1] = xyz_source[7][1];
                        src_cell.verts[7].pos.xyz[2] = xyz_source[7][2];

                        src_cell.verts[0].pnbrs[0] = 1;
                        src_cell.verts[0].pnbrs[1] = 4;
                        src_cell.verts[0].pnbrs[2] = 3;
                        src_cell.verts[1].pnbrs[0] = 2;
                        src_cell.verts[1].pnbrs[1] = 5;
                        src_cell.verts[1].pnbrs[2] = 0;
                        src_cell.verts[2].pnbrs[0] = 3;
                        src_cell.verts[2].pnbrs[1] = 6;
                        src_cell.verts[2].pnbrs[2] = 1;
                        src_cell.verts[3].pnbrs[0] = 0;
                        src_cell.verts[3].pnbrs[1] = 7;
                        src_cell.verts[3].pnbrs[2] = 2;
                        src_cell.verts[4].pnbrs[0] = 5;
                        src_cell.verts[4].pnbrs[1] = 7;
                        src_cell.verts[4].pnbrs[2] = 0;
                        src_cell.verts[5].pnbrs[0] = 6;
                        src_cell.verts[5].pnbrs[1] = 4;
                        src_cell.verts[5].pnbrs[2] = 1;
                        src_cell.verts[6].pnbrs[0] = 7;
                        src_cell.verts[6].pnbrs[1] = 5;
                        src_cell.verts[6].pnbrs[2] = 2;
                        src_cell.verts[7].pnbrs[0] = 4;
                        src_cell.verts[7].pnbrs[1] = 6;
                        src_cell.verts[7].pnbrs[2] = 3;

                        r3d_int isgood = r3d_is_good(&src_cell);
                        if (isgood == 0) {
                            std::cout << "NOT GOOD R3D" << std::endl;
                        }

                        while (boxList != NULL) {

                            GtsBBox *box = (GtsBBox *) (boxList->data);
                            kmds::TCellID cid_target = GPOINTER_TO_INT(box->bounded);

                            kmds::Region c_target = AMesh_target->getRegion(cid_target);


                            Kokkos::View<kmds::TCellID *> nids_target;
                            c_target.nodeIds(nids_target);

                            kmds::TCoord xyz_target[15][3];
                            AMesh_target->getNodeLocation(nids_target[0], xyz_target[0][0], xyz_target[0][1],
                                                          xyz_target[0][2]);
                            AMesh_target->getNodeLocation(nids_target[1], xyz_target[1][0], xyz_target[1][1],
                                                          xyz_target[1][2]);
                            AMesh_target->getNodeLocation(nids_target[2], xyz_target[2][0], xyz_target[2][1],
                                                          xyz_target[2][2]);
                            AMesh_target->getNodeLocation(nids_target[3], xyz_target[3][0], xyz_target[3][1],
                                                          xyz_target[3][2]);
                            AMesh_target->getNodeLocation(nids_target[4], xyz_target[4][0], xyz_target[4][1],
                                                          xyz_target[4][2]);
                            AMesh_target->getNodeLocation(nids_target[5], xyz_target[5][0], xyz_target[5][1],
                                                          xyz_target[5][2]);
                            AMesh_target->getNodeLocation(nids_target[6], xyz_target[6][0], xyz_target[6][1],
                                                          xyz_target[6][2]);
                            AMesh_target->getNodeLocation(nids_target[7], xyz_target[7][0], xyz_target[7][1],
                                                          xyz_target[7][2]);

                            gmds::math::Point pt0(xyz_target[0][0], xyz_target[0][1], xyz_target[0][2]);
                            gmds::math::Point pt1(xyz_target[1][0], xyz_target[1][1], xyz_target[1][2]);
                            gmds::math::Point pt2(xyz_target[2][0], xyz_target[2][1], xyz_target[2][2]);
                            gmds::math::Point pt3(xyz_target[3][0], xyz_target[3][1], xyz_target[3][2]);
                            gmds::math::Point pt4(xyz_target[4][0], xyz_target[4][1], xyz_target[4][2]);
                            gmds::math::Point pt5(xyz_target[5][0], xyz_target[5][1], xyz_target[5][2]);
                            gmds::math::Point pt6(xyz_target[6][0], xyz_target[6][1], xyz_target[6][2]);
                            gmds::math::Point pt7(xyz_target[7][0], xyz_target[7][1], xyz_target[7][2]);

                            gmds::math::Hexahedron h(pt0, pt1, pt2, pt3, pt4, pt5, pt6, pt7);
                            const double volsurf = h.getVolume();

                            r3d_rvec3 verts2[15];
                            verts2[0].x = xyz_target[0][0];
                            verts2[0].y = xyz_target[0][1];
                            verts2[0].z = xyz_target[0][2];
                            verts2[1].x = xyz_target[1][0];
                            verts2[1].y = xyz_target[1][1];
                            verts2[1].z = xyz_target[1][2];
                            verts2[2].x = xyz_target[2][0];
                            verts2[2].y = xyz_target[2][1];
                            verts2[2].z = xyz_target[2][2];
                            verts2[3].x = xyz_target[3][0];
                            verts2[3].y = xyz_target[3][1];
                            verts2[3].z = xyz_target[3][2];
                            verts2[4].x = xyz_target[4][0];
                            verts2[4].y = xyz_target[4][1];
                            verts2[4].z = xyz_target[4][2];
                            verts2[5].x = xyz_target[5][0];
                            verts2[5].y = xyz_target[5][1];
                            verts2[5].z = xyz_target[5][2];
                            verts2[6].x = xyz_target[6][0];
                            verts2[6].y = xyz_target[6][1];
                            verts2[6].z = xyz_target[6][2];
                            verts2[7].x = xyz_target[7][0];
                            verts2[7].y = xyz_target[7][1];
                            verts2[7].z = xyz_target[7][2];

                            gmds::math::Point pt_bottom = 0.25 * (pt0 + pt1 + pt2 + pt3);
                            gmds::math::Point pt_top = 0.25 * (pt4 + pt5 + pt6 + pt7);
                            gmds::math::Point pt_left = 0.25 * (pt0 + pt3 + pt7 + pt4);
                            gmds::math::Point pt_right = 0.25 * (pt1 + pt2 + pt6 + pt5);
                            gmds::math::Point pt_front = 0.25 * (pt0 + pt1 + pt5 + pt4);
                            gmds::math::Point pt_back = 0.25 * (pt3 + pt2 + pt6 + pt7);
                            verts2[8].x = pt_bottom.X();
                            verts2[8].y = pt_bottom.Y();
                            verts2[8].z = pt_bottom.Z();
                            verts2[9].x = pt_top.X();
                            verts2[9].y = pt_top.Y();
                            verts2[9].z = pt_top.Z();
                            verts2[10].x = pt_left.X();
                            verts2[10].y = pt_left.Y();
                            verts2[10].z = pt_left.Z();
                            verts2[11].x = pt_right.X();
                            verts2[11].y = pt_right.Y();
                            verts2[11].z = pt_right.Z();
                            verts2[12].x = pt_front.X();
                            verts2[12].y = pt_front.Y();
                            verts2[12].z = pt_front.Z();
                            verts2[13].x = pt_back.X();
                            verts2[13].y = pt_back.Y();
                            verts2[13].z = pt_back.Z();

                            gmds::math::Point pt_center = c_target.midpoint();
                            verts2[14].x = pt_center.X();
                            verts2[14].y = pt_center.Y();
                            verts2[14].z = pt_center.Z();

                            r3d_rvec3 verts_tet[4];
                            r3d_plane faces[4];
                            r3d_poly src_cell_copy;
                            const int POLY_ORDER = 1;
                            r3d_real intersectVol = 0.;
                            r3d_real om[R3D_NUM_MOMENTS(POLY_ORDER)];

                            // bottom
                            verts_tet[0] = verts2[0];
                            verts_tet[1] = verts2[1];
                            verts_tet[2] = verts2[8];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];
                            verts_tet[0] = verts2[1];
                            verts_tet[1] = verts2[2];
                            verts_tet[2] = verts2[8];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];
                            verts_tet[0] = verts2[2];
                            verts_tet[1] = verts2[3];
                            verts_tet[2] = verts2[8];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];
                            verts_tet[0] = verts2[3];
                            verts_tet[1] = verts2[0];
                            verts_tet[2] = verts2[8];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];

                            // top
                            verts_tet[0] = verts2[4];
                            verts_tet[1] = verts2[7];
                            verts_tet[2] = verts2[9];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];
                            verts_tet[0] = verts2[7];
                            verts_tet[1] = verts2[6];
                            verts_tet[2] = verts2[9];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];
                            verts_tet[0] = verts2[6];
                            verts_tet[1] = verts2[5];
                            verts_tet[2] = verts2[9];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];
                            verts_tet[0] = verts2[5];
                            verts_tet[1] = verts2[4];
                            verts_tet[2] = verts2[9];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];

                            // left
                            verts_tet[0] = verts2[0];
                            verts_tet[1] = verts2[3];
                            verts_tet[2] = verts2[10];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];
                            verts_tet[0] = verts2[3];
                            verts_tet[1] = verts2[7];
                            verts_tet[2] = verts2[10];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];
                            verts_tet[0] = verts2[7];
                            verts_tet[1] = verts2[4];
                            verts_tet[2] = verts2[10];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];
                            verts_tet[0] = verts2[4];
                            verts_tet[1] = verts2[0];
                            verts_tet[2] = verts2[10];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];

                            // right
                            verts_tet[0] = verts2[2];
                            verts_tet[1] = verts2[1];
                            verts_tet[2] = verts2[11];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];
                            verts_tet[0] = verts2[1];
                            verts_tet[1] = verts2[5];
                            verts_tet[2] = verts2[11];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];
                            verts_tet[0] = verts2[5];
                            verts_tet[1] = verts2[6];
                            verts_tet[2] = verts2[11];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];
                            verts_tet[0] = verts2[6];
                            verts_tet[1] = verts2[2];
                            verts_tet[2] = verts2[11];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];

                            // front
                            verts_tet[0] = verts2[1];
                            verts_tet[1] = verts2[0];
                            verts_tet[2] = verts2[12];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];
                            verts_tet[0] = verts2[0];
                            verts_tet[1] = verts2[4];
                            verts_tet[2] = verts2[12];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];
                            verts_tet[0] = verts2[4];
                            verts_tet[1] = verts2[5];
                            verts_tet[2] = verts2[12];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];
                            verts_tet[0] = verts2[5];
                            verts_tet[1] = verts2[1];
                            verts_tet[2] = verts2[12];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];

                            // back
                            verts_tet[0] = verts2[3];
                            verts_tet[1] = verts2[2];
                            verts_tet[2] = verts2[13];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];
                            verts_tet[0] = verts2[2];
                            verts_tet[1] = verts2[6];
                            verts_tet[2] = verts2[13];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];
                            verts_tet[0] = verts2[6];
                            verts_tet[1] = verts2[7];
                            verts_tet[2] = verts2[13];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];
                            verts_tet[0] = verts2[7];
                            verts_tet[1] = verts2[3];
                            verts_tet[2] = verts2[13];
                            verts_tet[3] = verts2[14];
                            r3d_tet_faces_from_verts(faces, verts_tet);
                            src_cell_copy = src_cell;
                            r3d_clip(&src_cell_copy, &faces[0], 4);
                            r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                            intersectVol += om[0];

                            double fracpres = Afp_target->getFracPres(matIndex, cid_target);

                            // add but do not forget to divide by cell area
                            double fracpresadd = intersectVol / volsurf;

                            if (fracpresadd < 1.e-10) {
                                fracpresadd = 0.;
                            }

                            if ((*varSumFracpres)[cid_target] + fracpresadd > 1.) {
                                fracpresadd = 1. - (*varSumFracpres)[cid_target];
                                (*varSumFracpres)[cid_target] = 1.;
                            } else {
                                (*varSumFracpres)[cid_target] += fracpresadd;
                            }

                            fracpres += fracpresadd;

                            Afp_target->setFracPres(matIndex, cid_target, fracpres);


                            boxList = g_slist_next(boxList);

                        }  // while (boxList != NULL)

//                    g_slist_free(boxList);

                    }
                        break;
                    default:
                        throw kmds::KException("Tools_compute_fracpres_3D : cell type not handled.");
                        break;
                }

            }  // if((*AVarOldCells2firstSubCells)[cid_source] == kmds::NullID)

        }


        kmds::GrowingView<kmds::TCellID> cellIDs_submesh("CELLS", AMesh_submesh->getNbRegions());
        AMesh_submesh->getRegionIDs(&cellIDs_submesh);
        const kmds::TCellID nbCells_submesh = cellIDs_submesh.getNbElems();

        for(int i=0; i<nbCells_submesh; i++) {

            kmds::TCellID cid_source = cellIDs_submesh.get(i);

            std::cout<<"cid_source "<<cid_source<<" of "<<nbCells_submesh<<std::endl;

            kmds::Region c_source = AMesh_submesh->getRegion(cid_source);

            double minXYZ[3];
            double maxXYZ[3];

            c_source.computeBoundingBox(minXYZ, maxXYZ);

            gts_bbox_set(
                    bbox,
                    pointer,
                    minXYZ[0], minXYZ[1], minXYZ[2],
                    maxXYZ[0], maxXYZ[1], maxXYZ[2]);

            GSList *boxList = gts_bb_tree_overlap(aabbCellsTree, bbox);
            if (boxList == NULL) {
                throw kmds::KException(
                        "Tools_compute_fracpres_3D : a cell from source mesh does not intersect target mesh.");
            }

            int matIndex = Ama_submesh->getMaterial(cid_source);


            switch (c_source.computeType()) {
                case kmds::KMDS_HEX : {
                    Kokkos::View<kmds::TCellID *> nids_source;
                    c_source.nodeIds(nids_source);

                    kmds::TCoord xyz_source[8][3];
                    AMesh_submesh->getNodeLocation(nids_source[0], xyz_source[0][0], xyz_source[0][1],
                                                   xyz_source[0][2]);
                    AMesh_submesh->getNodeLocation(nids_source[1], xyz_source[1][0], xyz_source[1][1],
                                                   xyz_source[1][2]);
                    AMesh_submesh->getNodeLocation(nids_source[2], xyz_source[2][0], xyz_source[2][1],
                                                   xyz_source[2][2]);
                    AMesh_submesh->getNodeLocation(nids_source[3], xyz_source[3][0], xyz_source[3][1],
                                                   xyz_source[3][2]);
                    AMesh_submesh->getNodeLocation(nids_source[4], xyz_source[4][0], xyz_source[4][1],
                                                   xyz_source[4][2]);
                    AMesh_submesh->getNodeLocation(nids_source[5], xyz_source[5][0], xyz_source[5][1],
                                                   xyz_source[5][2]);
                    AMesh_submesh->getNodeLocation(nids_source[6], xyz_source[6][0], xyz_source[6][1],
                                                   xyz_source[6][2]);
                    AMesh_submesh->getNodeLocation(nids_source[7], xyz_source[7][0], xyz_source[7][1],
                                                   xyz_source[7][2]);

//                    r3d_rvec3 verts[8];
//                    verts[0].x = xyz_source[0][0];
//                    verts[0].y = xyz_source[0][1];
//                    verts[0].z = xyz_source[0][2];
//                    verts[1].x = xyz_source[1][0];
//                    verts[1].y = xyz_source[1][1];
//                    verts[1].z = xyz_source[1][2];
//                    verts[2].x = xyz_source[2][0];
//                    verts[2].y = xyz_source[2][1];
//                    verts[2].z = xyz_source[2][2];
//                    verts[3].x = xyz_source[3][0];
//                    verts[3].y = xyz_source[3][1];
//                    verts[3].z = xyz_source[3][2];
//                    verts[4].x = xyz_source[4][0];
//                    verts[4].y = xyz_source[4][1];
//                    verts[4].z = xyz_source[4][2];
//                    verts[5].x = xyz_source[5][0];
//                    verts[5].y = xyz_source[5][1];
//                    verts[5].z = xyz_source[5][2];
//                    verts[6].x = xyz_source[6][0];
//                    verts[6].y = xyz_source[6][1];
//                    verts[6].z = xyz_source[6][2];
//                    verts[7].x = xyz_source[7][0];
//                    verts[7].y = xyz_source[7][1];
//                    verts[7].z = xyz_source[7][2];

//                    // create polyhedron
//                    r3d_poly src_r3dpoly;
//                    const int num_verts = 14;
//                    const int num_faces = 24;
//
//                    r3d_int face_num_verts[num_faces];
//                    for(int i_f = 0; i_f<num_faces; i_f++) {
//                        face_num_verts[i_f] = 3;
//                    }
//
//                    r3d_int **face_vert_ids = new r3d_int *[num_faces];
//                    for (int i_f = 0; i_f < num_faces; i_f++) {
//                        face_vert_ids[i_f] = new r3d_int[face_num_verts[i_f]];
//                    }
//                    face_vert_ids[0][0] = 0;  // bottom
//                    face_vert_ids[0][1] = 1;
//                    face_vert_ids[0][2] = 8;
//                    face_vert_ids[1][0] = 1;
//                    face_vert_ids[1][1] = 2;
//                    face_vert_ids[1][2] = 8;
//                    face_vert_ids[2][0] = 2;
//                    face_vert_ids[2][1] = 3;
//                    face_vert_ids[2][2] = 8;
//                    face_vert_ids[3][0] = 3;
//                    face_vert_ids[3][1] = 0;
//                    face_vert_ids[3][2] = 8;
//                    face_vert_ids[0][0] = 0;  // top
//                    face_vert_ids[0][1] = 1;
//                    face_vert_ids[0][2] = 9;
//                    face_vert_ids[1][0] = 1;
//                    face_vert_ids[1][1] = 2;
//                    face_vert_ids[1][2] = 9;
//                    face_vert_ids[2][0] = 2;
//                    face_vert_ids[2][1] = 3;
//                    face_vert_ids[2][2] = 9;
//                    face_vert_ids[3][0] = 3;
//                    face_vert_ids[3][1] = 0;
//                    face_vert_ids[3][2] = 9;
//
//                    r3d_init_poly(&src_r3dpoly, verts, num_verts, face_vert_ids, face_num_verts,
//                num_faces);

//                    gmds::math::Point pt_center = c_source.midpoint();
//                    r3d_rvec3 vec_center;
//                    vec_center.x = pt_center.X();
//                    vec_center.y = pt_center.Y();
//                    vec_center.z = pt_center.Z();
//

                    // create polyhedron
                    r3d_poly src_cell;
                    src_cell.nverts = 8;

                    src_cell.verts[0].pos.xyz[0] = xyz_source[0][0];
                    src_cell.verts[0].pos.xyz[1] = xyz_source[0][1];
                    src_cell.verts[0].pos.xyz[2] = xyz_source[0][2];
                    src_cell.verts[1].pos.xyz[0] = xyz_source[1][0];
                    src_cell.verts[1].pos.xyz[1] = xyz_source[1][1];
                    src_cell.verts[1].pos.xyz[2] = xyz_source[1][2];
                    src_cell.verts[2].pos.xyz[0] = xyz_source[2][0];
                    src_cell.verts[2].pos.xyz[1] = xyz_source[2][1];
                    src_cell.verts[2].pos.xyz[2] = xyz_source[2][2];
                    src_cell.verts[3].pos.xyz[0] = xyz_source[3][0];
                    src_cell.verts[3].pos.xyz[1] = xyz_source[3][1];
                    src_cell.verts[3].pos.xyz[2] = xyz_source[3][2];
                    src_cell.verts[4].pos.xyz[0] = xyz_source[4][0];
                    src_cell.verts[4].pos.xyz[1] = xyz_source[4][1];
                    src_cell.verts[4].pos.xyz[2] = xyz_source[4][2];
                    src_cell.verts[5].pos.xyz[0] = xyz_source[5][0];
                    src_cell.verts[5].pos.xyz[1] = xyz_source[5][1];
                    src_cell.verts[5].pos.xyz[2] = xyz_source[5][2];
                    src_cell.verts[6].pos.xyz[0] = xyz_source[6][0];
                    src_cell.verts[6].pos.xyz[1] = xyz_source[6][1];
                    src_cell.verts[6].pos.xyz[2] = xyz_source[6][2];
                    src_cell.verts[7].pos.xyz[0] = xyz_source[7][0];
                    src_cell.verts[7].pos.xyz[1] = xyz_source[7][1];
                    src_cell.verts[7].pos.xyz[2] = xyz_source[7][2];

                    src_cell.verts[0].pnbrs[0] = 1;
                    src_cell.verts[0].pnbrs[1] = 4;
                    src_cell.verts[0].pnbrs[2] = 3;
                    src_cell.verts[1].pnbrs[0] = 2;
                    src_cell.verts[1].pnbrs[1] = 5;
                    src_cell.verts[1].pnbrs[2] = 0;
                    src_cell.verts[2].pnbrs[0] = 3;
                    src_cell.verts[2].pnbrs[1] = 6;
                    src_cell.verts[2].pnbrs[2] = 1;
                    src_cell.verts[3].pnbrs[0] = 0;
                    src_cell.verts[3].pnbrs[1] = 7;
                    src_cell.verts[3].pnbrs[2] = 2;
                    src_cell.verts[4].pnbrs[0] = 5;
                    src_cell.verts[4].pnbrs[1] = 7;
                    src_cell.verts[4].pnbrs[2] = 0;
                    src_cell.verts[5].pnbrs[0] = 6;
                    src_cell.verts[5].pnbrs[1] = 4;
                    src_cell.verts[5].pnbrs[2] = 1;
                    src_cell.verts[6].pnbrs[0] = 7;
                    src_cell.verts[6].pnbrs[1] = 5;
                    src_cell.verts[6].pnbrs[2] = 2;
                    src_cell.verts[7].pnbrs[0] = 4;
                    src_cell.verts[7].pnbrs[1] = 6;
                    src_cell.verts[7].pnbrs[2] = 3;

                    r3d_int isgood = r3d_is_good(&src_cell);
                    if (isgood == 0) {
                        std::cout << "NOT GOOD R3D" << std::endl;
                    }

                    while (boxList != NULL) {

                        GtsBBox *box = (GtsBBox *) (boxList->data);
                        kmds::TCellID cid_target = GPOINTER_TO_INT(box->bounded);

                        kmds::Region c_target = AMesh_target->getRegion(cid_target);


                        Kokkos::View<kmds::TCellID *> nids_target;
                        c_target.nodeIds(nids_target);

                        kmds::TCoord xyz_target[15][3];
                        AMesh_target->getNodeLocation(nids_target[0], xyz_target[0][0], xyz_target[0][1],
                                                      xyz_target[0][2]);
                        AMesh_target->getNodeLocation(nids_target[1], xyz_target[1][0], xyz_target[1][1],
                                                      xyz_target[1][2]);
                        AMesh_target->getNodeLocation(nids_target[2], xyz_target[2][0], xyz_target[2][1],
                                                      xyz_target[2][2]);
                        AMesh_target->getNodeLocation(nids_target[3], xyz_target[3][0], xyz_target[3][1],
                                                      xyz_target[3][2]);
                        AMesh_target->getNodeLocation(nids_target[4], xyz_target[4][0], xyz_target[4][1],
                                                      xyz_target[4][2]);
                        AMesh_target->getNodeLocation(nids_target[5], xyz_target[5][0], xyz_target[5][1],
                                                      xyz_target[5][2]);
                        AMesh_target->getNodeLocation(nids_target[6], xyz_target[6][0], xyz_target[6][1],
                                                      xyz_target[6][2]);
                        AMesh_target->getNodeLocation(nids_target[7], xyz_target[7][0], xyz_target[7][1],
                                                      xyz_target[7][2]);

                        gmds::math::Point pt0(xyz_target[0][0], xyz_target[0][1], xyz_target[0][2]);
                        gmds::math::Point pt1(xyz_target[1][0], xyz_target[1][1], xyz_target[1][2]);
                        gmds::math::Point pt2(xyz_target[2][0], xyz_target[2][1], xyz_target[2][2]);
                        gmds::math::Point pt3(xyz_target[3][0], xyz_target[3][1], xyz_target[3][2]);
                        gmds::math::Point pt4(xyz_target[4][0], xyz_target[4][1], xyz_target[4][2]);
                        gmds::math::Point pt5(xyz_target[5][0], xyz_target[5][1], xyz_target[5][2]);
                        gmds::math::Point pt6(xyz_target[6][0], xyz_target[6][1], xyz_target[6][2]);
                        gmds::math::Point pt7(xyz_target[7][0], xyz_target[7][1], xyz_target[7][2]);

                        gmds::math::Hexahedron h(pt0, pt1, pt2, pt3, pt4, pt5, pt6, pt7);
                        const double volsurf = h.getVolume();

                        r3d_rvec3 verts2[15];
                        verts2[0].x = xyz_target[0][0];
                        verts2[0].y = xyz_target[0][1];
                        verts2[0].z = xyz_target[0][2];
                        verts2[1].x = xyz_target[1][0];
                        verts2[1].y = xyz_target[1][1];
                        verts2[1].z = xyz_target[1][2];
                        verts2[2].x = xyz_target[2][0];
                        verts2[2].y = xyz_target[2][1];
                        verts2[2].z = xyz_target[2][2];
                        verts2[3].x = xyz_target[3][0];
                        verts2[3].y = xyz_target[3][1];
                        verts2[3].z = xyz_target[3][2];
                        verts2[4].x = xyz_target[4][0];
                        verts2[4].y = xyz_target[4][1];
                        verts2[4].z = xyz_target[4][2];
                        verts2[5].x = xyz_target[5][0];
                        verts2[5].y = xyz_target[5][1];
                        verts2[5].z = xyz_target[5][2];
                        verts2[6].x = xyz_target[6][0];
                        verts2[6].y = xyz_target[6][1];
                        verts2[6].z = xyz_target[6][2];
                        verts2[7].x = xyz_target[7][0];
                        verts2[7].y = xyz_target[7][1];
                        verts2[7].z = xyz_target[7][2];

                        gmds::math::Point pt_bottom = 0.25 * (pt0 + pt1 + pt2 + pt3);
                        gmds::math::Point pt_top = 0.25 * (pt4 + pt5 + pt6 + pt7);
                        gmds::math::Point pt_left = 0.25 * (pt0 + pt3 + pt7 + pt4);
                        gmds::math::Point pt_right = 0.25 * (pt1 + pt2 + pt6 + pt5);
                        gmds::math::Point pt_front = 0.25 * (pt0 + pt1 + pt5 + pt4);
                        gmds::math::Point pt_back = 0.25 * (pt3 + pt2 + pt6 + pt7);
                        verts2[8].x = pt_bottom.X();
                        verts2[8].y = pt_bottom.Y();
                        verts2[8].z = pt_bottom.Z();
                        verts2[9].x = pt_top.X();
                        verts2[9].y = pt_top.Y();
                        verts2[9].z = pt_top.Z();
                        verts2[10].x = pt_left.X();
                        verts2[10].y = pt_left.Y();
                        verts2[10].z = pt_left.Z();
                        verts2[11].x = pt_right.X();
                        verts2[11].y = pt_right.Y();
                        verts2[11].z = pt_right.Z();
                        verts2[12].x = pt_front.X();
                        verts2[12].y = pt_front.Y();
                        verts2[12].z = pt_front.Z();
                        verts2[13].x = pt_back.X();
                        verts2[13].y = pt_back.Y();
                        verts2[13].z = pt_back.Z();

                        gmds::math::Point pt_center = c_target.midpoint();
                        verts2[14].x = pt_center.X();
                        verts2[14].y = pt_center.Y();
                        verts2[14].z = pt_center.Z();

                        r3d_rvec3 verts_tet[4];
                        r3d_plane faces[4];
                        r3d_poly src_cell_copy;
                        const int POLY_ORDER = 1;
                        r3d_real intersectVol = 0.;
                        r3d_real om[R3D_NUM_MOMENTS(POLY_ORDER)];

                        // bottom
                        verts_tet[0] = verts2[0];
                        verts_tet[1] = verts2[1];
                        verts_tet[2] = verts2[8];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[1];
                        verts_tet[1] = verts2[2];
                        verts_tet[2] = verts2[8];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[2];
                        verts_tet[1] = verts2[3];
                        verts_tet[2] = verts2[8];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[3];
                        verts_tet[1] = verts2[0];
                        verts_tet[2] = verts2[8];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];

                        // top
                        verts_tet[0] = verts2[4];
                        verts_tet[1] = verts2[7];
                        verts_tet[2] = verts2[9];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[7];
                        verts_tet[1] = verts2[6];
                        verts_tet[2] = verts2[9];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[6];
                        verts_tet[1] = verts2[5];
                        verts_tet[2] = verts2[9];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[5];
                        verts_tet[1] = verts2[4];
                        verts_tet[2] = verts2[9];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];

                        // left
                        verts_tet[0] = verts2[0];
                        verts_tet[1] = verts2[3];
                        verts_tet[2] = verts2[10];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[3];
                        verts_tet[1] = verts2[7];
                        verts_tet[2] = verts2[10];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[7];
                        verts_tet[1] = verts2[4];
                        verts_tet[2] = verts2[10];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[4];
                        verts_tet[1] = verts2[0];
                        verts_tet[2] = verts2[10];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];

                        // right
                        verts_tet[0] = verts2[2];
                        verts_tet[1] = verts2[1];
                        verts_tet[2] = verts2[11];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[1];
                        verts_tet[1] = verts2[5];
                        verts_tet[2] = verts2[11];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[5];
                        verts_tet[1] = verts2[6];
                        verts_tet[2] = verts2[11];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[6];
                        verts_tet[1] = verts2[2];
                        verts_tet[2] = verts2[11];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];

                        // front
                        verts_tet[0] = verts2[1];
                        verts_tet[1] = verts2[0];
                        verts_tet[2] = verts2[12];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[0];
                        verts_tet[1] = verts2[4];
                        verts_tet[2] = verts2[12];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[4];
                        verts_tet[1] = verts2[5];
                        verts_tet[2] = verts2[12];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[5];
                        verts_tet[1] = verts2[1];
                        verts_tet[2] = verts2[12];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];

                        // back
                        verts_tet[0] = verts2[3];
                        verts_tet[1] = verts2[2];
                        verts_tet[2] = verts2[13];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[2];
                        verts_tet[1] = verts2[6];
                        verts_tet[2] = verts2[13];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[6];
                        verts_tet[1] = verts2[7];
                        verts_tet[2] = verts2[13];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];
                        verts_tet[0] = verts2[7];
                        verts_tet[1] = verts2[3];
                        verts_tet[2] = verts2[13];
                        verts_tet[3] = verts2[14];
                        r3d_tet_faces_from_verts(faces, verts_tet);
                        src_cell_copy = src_cell;
                        r3d_clip(&src_cell_copy, &faces[0], 4);
                        r3d_reduce(&src_cell_copy, om, POLY_ORDER);
                        intersectVol += om[0];

                        double fracpres = Afp_target->getFracPres(matIndex, cid_target);

                        // add but do not forget to divide by cell area
                        double fracpresadd = intersectVol / volsurf;

                        if (fracpresadd < 1.e-10) {
                            fracpresadd = 0.;
                        }

                        if ((*varSumFracpres)[cid_target] + fracpresadd > 1.) {
                            fracpresadd = 1. - (*varSumFracpres)[cid_target];
                            (*varSumFracpres)[cid_target] = 1.;
                        } else {
                            (*varSumFracpres)[cid_target] += fracpresadd;
                        }

                        fracpres += fracpresadd;

                        Afp_target->setFracPres(matIndex, cid_target, fracpres);


                        boxList = g_slist_next(boxList);

                    }  // while (boxList != NULL)

//                    g_slist_free(boxList);

                }
                    break;
                default:
                    throw kmds::KException("Tools_compute_fracpres_3D : cell type not handled.");
                    break;
            }

        }

        // clean-up
        gts_bb_tree_destroy(aabbCellsTree, true);

        AMesh_target->deleteVariable(kmds::KMDS_REGION, "tmp_varSumFracpres");
    }

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
