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
/** \file    InterfaceNodesPos.cpp
 *  \author  legoff
 *  \date    03/26/2018
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/InterfaceNodesPos.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <set>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
// Usage of the Eigen Template library
#include <Eigen/Dense>
#include <Kokkos_UnorderedMap.hpp>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/math/Line.h>
#include <gmds/math/Plane.h>
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/MaterialInterfaces.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {

/*----------------------------------------------------------------------------*/
    struct InterfaceNodesPos_computePos_onenode_2D {

        const kmds::GrowingView<kmds::TCellID> *nodeIDsAccessor;

        kmds::Mesh *mesh;
        const kmds::Connectivity *c_N2F;
        const kmds::Variable<gmds::math::Point> *varPCrossPt;
        const kmds::Variable<gmds::math::Vector> *varPCrossVec;

        kmds::Variable<gmds::math::Point> *varNewPos;


        InterfaceNodesPos_computePos_onenode_2D(const kmds::GrowingView<kmds::TCellID> *ANodeIDsAccessor_,
                                                kmds::Mesh *AMesh_,
                                                const kmds::Connectivity *c_N2F_,
                                                const kmds::Variable<gmds::math::Point> *AVarPCrossPt_,
                                                const kmds::Variable<gmds::math::Vector> *AVarPCrossVec_,
                                                kmds::Variable<gmds::math::Point> *AVarNewPos
        )
                : nodeIDsAccessor(ANodeIDsAccessor_)
                , mesh(AMesh_)
                , c_N2F(c_N2F_)
                , varPCrossPt(AVarPCrossPt_)
                , varPCrossVec(AVarPCrossVec_)
                , varNewPos(AVarNewPos) {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int nid = nodeIDsAccessor->get(i);

            kmds::TCoord coords[3];
            mesh->getNodeLocation(nid, coords[0], coords[1], coords[2]);
            gmds::math::Point oldPos(coords[0], coords[1], coords[2]);

            Kokkos::View<kmds::TCellID*> faces;
            c_N2F->get(nid, faces);


            gmds::math::Point newPos(0.,0.,0.);
            int nbplanes = 0;
            for(int i_f=0; i_f<faces.size(); i_f++) {
                gmds::math::Point pcrosspt = (*varPCrossPt)[faces(i_f)];
                gmds::math::Vector pcrossvec = (*varPCrossVec)[faces(i_f)];

                if(pcrosspt != gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF)) {

                    gmds::math::Vector pcrossvec_orth(-pcrossvec.Y(), pcrossvec.X(), 0.);

                    gmds::math::Line pl(pcrosspt, pcrossvec_orth);
                    gmds::math::Point proj = pl.project(oldPos);
                    newPos = newPos + proj;
                    nbplanes++;

                }
            }

            newPos = (1./nbplanes) * newPos;

            (*varNewPos)[nid] = newPos;

//            exit(-1);
        }
    };
/*----------------------------------------------------------------------------*/
    struct InterfaceNodesPos_computePos_onenode_3D {

        const kmds::GrowingView<kmds::TCellID> *nodeIDsAccessor;

        kmds::Mesh *mesh;
        const kmds::Connectivity *c_N2F;
        const kmds::Variable<gmds::math::Point> *varPCrossPt;
        const kmds::Variable<gmds::math::Vector> *varPCrossVec;

        kmds::Variable<gmds::math::Point> *varNewPos;


        InterfaceNodesPos_computePos_onenode_3D(const kmds::GrowingView<kmds::TCellID> *ANodeIDsAccessor_,
                                                kmds::Mesh *AMesh_,
                                                const kmds::Connectivity *c_N2F_,
                                                const kmds::Variable<gmds::math::Point> *AVarPCrossPt_,
                                                const kmds::Variable<gmds::math::Vector> *AVarPCrossVec_,
                                                kmds::Variable<gmds::math::Point> *AVarNewPos
        )
                : nodeIDsAccessor(ANodeIDsAccessor_)
                , mesh(AMesh_)
                , c_N2F(c_N2F_)
                , varPCrossPt(AVarPCrossPt_)
                , varPCrossVec(AVarPCrossVec_)
                , varNewPos(AVarNewPos) {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int nid = nodeIDsAccessor->get(i);

            kmds::TCoord coords[3];
            mesh->getNodeLocation(nid, coords[0], coords[1], coords[2]);
            gmds::math::Point oldPos(coords[0], coords[1], coords[2]);

            Kokkos::View<kmds::TCellID*> faces;
            c_N2F->get(nid, faces);

            gmds::math::Point newPos(0.,0.,0.);
            int nbplanes = 0;
            for(int i_f=0; i_f<faces.size(); i_f++) {
                gmds::math::Point pcrosspt = (*varPCrossPt)[faces(i_f)];
                gmds::math::Vector pcrossvec = (*varPCrossVec)[faces(i_f)];

                if(pcrosspt != gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF)) {

                    gmds::math::Plane pl(pcrosspt, pcrossvec);
                    gmds::math::Point proj = pl.project(oldPos);
                    newPos = newPos + proj;
                    nbplanes++;

//                    if(nid == 342) {
//                        std::cout << "nid " << "oldPos " << oldPos << std::endl;
//                        std::cout << "faces(i_f) " << faces(i_f) << std::endl;
//                        std::cout << pcrosspt << " " << pcrossvec << std::endl;
//                        std::cout << "proj " << proj << std::endl;
//
//                        kmds::Face f = mesh->getFace(faces(i_f));
//                        std::cout<<"midpoint "<<f.midpoint()<<std::endl;
////                        Kokkos::View<kmds::TCellID *> nodesids;
////                        f.nodeIds(nodesids);
////
////                        for()
//                    }

                }
            }

            assert(nbplanes > 0);
//            std::cout<<"newPos "<<nid<<" "<<newPos<<" "<<nbplanes<<std::endl;

            newPos = (1./nbplanes) * newPos;

            (*varNewPos)[nid] = newPos;

//            exit(-1);
        }
    };
/*----------------------------------------------------------------------------*/
    struct InterfaceNodesPos_computePcross_oneface_xD {

        const kmds::GrowingView<kmds::TCellID> *faceIDsAccessor;

        kmds::Mesh *mesh;
        const kmds::Connectivity *c_F2C;
        const elg3d::FracPres *fp;
        const elg3d::MaterialAssignment *ma;
        const kmds::Variable<gmds::math::Point> *varMidpoints;
        const kmds::Variable<MaterialGradientComputation_midcellGradients>* varGrad;
        kmds::Variable<gmds::math::Point> *varPCrossPt;
        kmds::Variable<gmds::math::Vector> *varPCrossVec;


        InterfaceNodesPos_computePcross_oneface_xD(const kmds::GrowingView<kmds::TCellID> *AFaceIDsAccessor_,
                                                   kmds::Mesh *AMesh_,
                                                   const kmds::Connectivity *c_F2C_,
                                                   const elg3d::FracPres *fp_,
                                                   const elg3d::MaterialAssignment *ma_,
                                                   const kmds::Variable<gmds::math::Point> *AVarMidpoints_,
                                                   const kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrad_,
                                                   kmds::Variable<gmds::math::Point> *AVarPCrossPt_,
                                                   kmds::Variable<gmds::math::Vector> *AVarPCrossVec_
        )
        : faceIDsAccessor(AFaceIDsAccessor_)
        , mesh(AMesh_)
        , c_F2C(c_F2C_)
        , fp(fp_)
        , ma(ma_)
        , varMidpoints(AVarMidpoints_)
        , varGrad(AVarGrad_)
        , varPCrossPt(AVarPCrossPt_)
        , varPCrossVec(AVarPCrossVec_){
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int fid = faceIDsAccessor->get(i);

            Kokkos::View<kmds::TCellID*> cells;
            c_F2C->get(fid, cells);

            kmds::TCellID cid0 = cells(0);
            kmds::TCellID cid1 = cells(1);

            int mat0 = ma->getMaterial(cid0);
            int mat1 = ma->getMaterial(cid1);

            double fp0m0 = fp->getFracPres(mat0, cid0);
            double fp0m1 = fp->getFracPres(mat1, cid0);
            double fp1m0 = fp->getFracPres(mat0, cid1);
            double fp1m1 = fp->getFracPres(mat1, cid1);

            double d = fp1m0 - fp0m0 - fp1m1 + fp0m1;
            double tt = 0.5;

            // TODO hard-coded value
            if (fabs(d) >= 1e-5)
                tt = (fp0m1 - fp0m0) / d;

            // TODO check domain of validity
            if(tt < 0.) tt = 0.;
            if(tt > 1.) tt = 1.;

            gmds::math::Point p0 = (*varMidpoints)[cid0];
            gmds::math::Point p1 = (*varMidpoints)[cid1];
            gmds::math::Point pt(p0 + tt * (p1 - p0));
            (*varPCrossPt)[fid] = pt;

            // TODO gradient vectors
            int indexc0m0 = -1;
            int indexc1m0 = -1;
            int indexc0m1 = -1;
            int indexc1m1 = -1;
            for(int imat=0; imat<MaterialGradientComputation_MAXNBMATPERCELLBALL; imat++) {
                if((*varGrad)[cid0].m_matindex[imat] == mat0) {
                    indexc0m0 = imat;
                }
                if((*varGrad)[cid1].m_matindex[imat] == mat0) {
                    indexc1m0 = imat;
                }
                if((*varGrad)[cid0].m_matindex[imat] == mat1) {
                    indexc0m1 = imat;
                }
                if((*varGrad)[cid1].m_matindex[imat] == mat1) {
                    indexc1m1 = imat;
                }
            }

            assert(indexc0m0 != -1);
            assert(indexc1m0 != -1);
            assert(indexc0m1 != -1);
            assert(indexc1m1 != -1);

            gmds::math::Vector v0m0 = (*varGrad)[cid0].m_grad[indexc0m0];
            gmds::math::Vector v1m0 = (*varGrad)[cid1].m_grad[indexc1m0];
            gmds::math::Vector v0m1 = (*varGrad)[cid0].m_grad[indexc0m1];
            gmds::math::Vector v1m1 = (*varGrad)[cid1].m_grad[indexc1m1];

            if(v0m0.isZero() || v1m0.isZero() || v0m1.isZero() || v1m1.isZero()) {
                (*varPCrossPt)[fid] = gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF);
                return;
            }
            if((v0m0 == gmds::math::Vector(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF))
            || (v1m0 == gmds::math::Vector(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF))
            || (v0m1 == gmds::math::Vector(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF))
            || (v1m1 == gmds::math::Vector(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF))) {
                (*varPCrossPt)[fid] = gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF);
                return;
            }


            gmds::math::Vector v0(v0m0 + tt * (v1m0 - v0m0));
            gmds::math::Vector v1(v0m1 + tt * (v1m1 - v0m1));
            v0.normalize();
            v1.normalize();

            // minus because in the two materials case they should be in opposite directions
            (*varPCrossVec)[fid] = ((1./2.) * (v0 - v1)).getNormalize();

        }
    };
/*----------------------------------------------------------------------------*/
    bool
    InterfaceNodesPos_computePcross_oneface_onemat_3D(const kmds::TCellID Acid0,
                                    const kmds::TCellID Acid1,
                                    const int AMat,
                                    const elg3d::FracPres* Afp,
                                    const kmds::Variable<gmds::math::Point>* AVarCenter,
                                    gmds::math::Point& APcross)
    {
        gmds::math::Point p0 = (*AVarCenter)[Acid0];
        gmds::math::Point p1 = (*AVarCenter)[Acid1];

        kmds::TFloat fp0 = Afp->getFracPres(AMat, Acid0);
        kmds::TFloat fp1 = Afp->getFracPres(AMat, Acid1);

        // special case : fp value is the same in both cells
        if(fp0 == fp1) {
            APcross = 1./2. * (p0 + p1);
            return false;
        }

        // special case : the expected iso value is outside the range [fp0 .. fp1] and is lower
        if((InterfaceNodesPos_ISOVALUEFP < fp0) && (InterfaceNodesPos_ISOVALUEFP < fp1)) {
            if(fp0 < fp1) {
                APcross = p0;
            } else {
                APcross = p1;
            }
            return false;
        }

        // special case : the expected iso value is outside the range [fp0 .. fp1] and is higher
        if((InterfaceNodesPos_ISOVALUEFP > fp0) && (InterfaceNodesPos_ISOVALUEFP > fp1)) {
            if(fp0 < fp1) {
                APcross = p1;
            } else {
                APcross = p0;
            }
            return false;
        }

        // this is the normal case : linear interpolation
        APcross = p0 + ((InterfaceNodesPos_ISOVALUEFP - fp0) / (fp1 - fp0)) * (p1 - p0);
        return true;
    }
/*----------------------------------------------------------------------------*/
    void
    InterfaceNodesPos_computeNodesNewPos_2D(const kmds::GrowingView<kmds::TCellID>* ASelectionInterfaceNodes,
                                            kmds::Mesh* AMesh,
                                            const kmds::Connectivity* c_N2F,
                                            const kmds::Connectivity* c_F2C,
                                            const MaterialAssignment* Ama,
                                            const elg3d::FracPres* Afp,
                                            const kmds::Variable<gmds::math::Point>* AVarMidpoints,
                                            const kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads,
                                            kmds::Variable<gmds::math::Point>* AVarNewPos)
    {

        // compute dual edge intersections near interfaces
        // first get the interface faces adjacent to interface nodes
//        kmds::TSize fcapacity = AMesh->getFaceCapacity();
//        Kokkos::UnorderedMap<kmds::TCellID, void> interfaceFacesSet(fcapacity);
//        Kokkos::parallel_for(ASelectionInterfaceNodes->getNbElems(),
//                             KOKKOS_LAMBDA(const int i) {
//                                 kmds::TCellID nid = ASelectionInterfaceNodes->get(i);
//                                 Kokkos::View<kmds::TCellID *> fids;
//                                 c_N2F->get(nid, fids);
//                                 for(int i_f=0; i_f<fids.size(); i_f++) {
//
//                                     Kokkos::View<kmds::TCellID *> cids;
//                                     c_F2C->get(fids(i_f), cids);
//
//                                     // check if it is an interface face
//                                     // we ignore boundary faces
//                                     if(cids.size() == 2) {
//
//                                        if(Ama->getMaterial(cids(0)) != Ama->getMaterial(cids(1))) {
//                                            Kokkos::UnorderedMapInsertResult res = interfaceFacesSet.insert(fids(i_f));
//                                            bool success = res.success();
//                                            bool exist = res.existing();
//                                            bool fail = res.failed();
//                                            assert(success || exist);
//                                            assert(!fail);
//                                        }
//                                     }
//                                 }
//                             }
//        );
//
//        std::cout<<"nbinterfacenodes "<<ASelectionInterfaceNodes->getNbElems()<<std::endl;
//        std::cout<<"facescapacity "<<fcapacity<<std::endl;
//        std::cout<<"nbfaces "<<AMesh->getNbFaces()<<std::endl;
//        std::cout<<"nbinterfacesfaces "<<interfaceFacesSet.size()<<std::endl;
        kmds::GrowingView<kmds::TCellID> edgesInterfaces("EDGES_ON_INTERFACES", AMesh->getNbEdges());
        elg3d::MaterialInterfaces_getEdgeOnInterfaces(AMesh, c_F2C, Ama, &edgesInterfaces);

        // then compute the pcross for each interface face
        kmds::Variable<gmds::math::Point>* AVarPCrossPt =
                AMesh->createVariable<gmds::math::Point>(gmds::math::Point (-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_EDGE, "tmp_pcrosspt");
        kmds::Variable<gmds::math::Vector>* AVarPCrossVec =
                AMesh->createVariable<gmds::math::Vector>(gmds::math::Vector (-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_EDGE, "tmp_pcrossvec");
        kmds::TSize nbInterfaceEdges =  edgesInterfaces.getNbElems();


        Kokkos::parallel_for(nbInterfaceEdges,
                             InterfaceNodesPos_computePcross_oneface_xD(&edgesInterfaces,
                                                                        AMesh,
                                                                        c_F2C,
                                                                        Afp,
                                                                        Ama,
                                                                        AVarMidpoints,
                                                                        AVarGrads,
                                                                        AVarPCrossPt,
                                                                        AVarPCrossVec
                             ));


        // for each interface node, compute its new position
        kmds::TSize nbInterfaceNodes = ASelectionInterfaceNodes->getNbElems();
        Kokkos::parallel_for(nbInterfaceNodes,
                             InterfaceNodesPos_computePos_onenode_2D(ASelectionInterfaceNodes,
                                                                     AMesh,
                                                                     c_N2F,
                                                                     AVarPCrossPt,
                                                                     AVarPCrossVec,
                                                                     AVarNewPos
                             )
        );


        // clean-up
        AMesh->deleteVariable(kmds::KMDS_EDGE, "tmp_pcrosspt");
        AMesh->deleteVariable(kmds::KMDS_EDGE, "tmp_pcrossvec");
    }
/*----------------------------------------------------------------------------*/
    void
    InterfaceNodesPos_computeNodesNewPos_3D(const kmds::GrowingView<kmds::TCellID>* ASelectionInterfaceNodes,
                                            kmds::Mesh* AMesh,
                                            const kmds::Connectivity* c_N2F,
                                            const kmds::Connectivity* c_F2C,
                                            const MaterialAssignment* Ama,
                                            const elg3d::FracPres* Afp,
                                            const kmds::Variable<gmds::math::Point>* AVarMidpoints,
                                            const kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads,
                                            kmds::Variable<gmds::math::Point>* AVarNewPos)
    {

        // compute dual edge intersections near interfaces
        // first get the interface faces adjacent to interface nodes
//        kmds::TSize fcapacity = AMesh->getFaceCapacity();
//        Kokkos::UnorderedMap<kmds::TCellID, void> interfaceFacesSet(fcapacity);
//        Kokkos::parallel_for(ASelectionInterfaceNodes->getNbElems(),
//                             KOKKOS_LAMBDA(const int i) {
//                                 kmds::TCellID nid = ASelectionInterfaceNodes->get(i);
//                                 Kokkos::View<kmds::TCellID *> fids;
//                                 c_N2F->get(nid, fids);
//                                 for(int i_f=0; i_f<fids.size(); i_f++) {
//
//                                     Kokkos::View<kmds::TCellID *> cids;
//                                     c_F2C->get(fids(i_f), cids);
//
//                                     // check if it is an interface face
//                                     // we ignore boundary faces
//                                     if(cids.size() == 2) {
//
//                                        if(Ama->getMaterial(cids(0)) != Ama->getMaterial(cids(1))) {
//                                            Kokkos::UnorderedMapInsertResult res = interfaceFacesSet.insert(fids(i_f));
//                                            bool success = res.success();
//                                            bool exist = res.existing();
//                                            bool fail = res.failed();
//                                            assert(success || exist);
//                                            assert(!fail);
//                                        }
//                                     }
//                                 }
//                             }
//        );
//
//        std::cout<<"nbinterfacenodes "<<ASelectionInterfaceNodes->getNbElems()<<std::endl;
//        std::cout<<"facescapacity "<<fcapacity<<std::endl;
//        std::cout<<"nbfaces "<<AMesh->getNbFaces()<<std::endl;
//        std::cout<<"nbinterfacesfaces "<<interfaceFacesSet.size()<<std::endl;
        kmds::GrowingView<kmds::TCellID> facesInterfaces("FACES_ON_INTERFACES", AMesh->getNbFaces());
        elg3d::MaterialInterfaces_getFaceOnInterfaces(AMesh, c_F2C, Ama, &facesInterfaces);

        // then compute the pcross for each interface face
        kmds::Variable<gmds::math::Point>* AVarPCrossPt =
                AMesh->createVariable<gmds::math::Point>(gmds::math::Point (-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_FACE, "tmp_pcrosspt");
        kmds::Variable<gmds::math::Vector>* AVarPCrossVec =
                AMesh->createVariable<gmds::math::Vector>(gmds::math::Vector (-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_FACE, "tmp_pcrossvec");


        kmds::TSize nbInterfaceFaces =  facesInterfaces.getNbElems();
        Kokkos::parallel_for(nbInterfaceFaces,
                             InterfaceNodesPos_computePcross_oneface_xD(&facesInterfaces,
                                                                        AMesh,
                                                                        c_F2C,
                                                                        Afp,
                                                                        Ama,
                                                                        AVarMidpoints,
                                                                        AVarGrads,
                                                                        AVarPCrossPt,
                                                                        AVarPCrossVec
                             ));


        // for each interface node, compute its new position
        kmds::TSize nbInterfaceNodes = ASelectionInterfaceNodes->getNbElems();
        Kokkos::parallel_for(nbInterfaceNodes,
                             InterfaceNodesPos_computePos_onenode_3D(ASelectionInterfaceNodes,
                                                                     AMesh,
                                                                     c_N2F,
                                                                     AVarPCrossPt,
                                                                     AVarPCrossVec,
                                                                     AVarNewPos
                             )
        );


        // clean-up
        AMesh->deleteVariable(kmds::KMDS_FACE, "tmp_pcrosspt");
        AMesh->deleteVariable(kmds::KMDS_FACE, "tmp_pcrossvec");
    }
/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
