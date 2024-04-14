/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    MaterialGradientComputation.cpp
 *  \author  legoff
 *  \date    05/01/2018
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/MaterialGradientComputation.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <set>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
// Usage of the Eigen Template library
#include <Eigen/Dense>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/math/Plane.h>
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
namespace elg3d {

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
    struct MaterialGradientComputation_computeGrad_leastsquare_onecell_2D {

        const kmds::GrowingView<kmds::TCellID> *cellIDsAccessor;

        kmds::Mesh *mesh;
        const kmds::Connectivity *c_C2C_byN;
        const elg3d::FracPres *fp;
        const elg3d::MaterialAssignment *ma;
        const kmds::Variable<gmds::math::Point>* varMidpoints;
        kmds::Variable<MaterialGradientComputation_midcellGradients>* varGrads;


        MaterialGradientComputation_computeGrad_leastsquare_onecell_2D(const kmds::GrowingView<kmds::TCellID>* ACellIDsAccessor_,
                                                                       const kmds::Connectivity* c_C2C_byN_,
                                                                       const elg3d::FracPres* fp_,
                                                                       const elg3d::MaterialAssignment *ma_,
                                                                       const kmds::Variable<gmds::math::Point>* AVarMidpoints_,
                                                                       kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads_
        )
                : cellIDsAccessor(ACellIDsAccessor_)
                , c_C2C_byN(c_C2C_byN_)
                , fp(fp_)
                , ma(ma_)
                , varMidpoints(AVarMidpoints_)
                , varGrads(AVarGrads_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int cid = cellIDsAccessor->get(i);

            // list the neighboring cells (neighbors by node)
            Kokkos::View<kmds::TCellID *> neighborCellsIDs;
            c_C2C_byN->get(cid, neighborCellsIDs);
            int nbCells = neighborCellsIDs.size();

            // list the materials we will compute the gradient of
            std::set<int> materials;
            materials.insert(ma->getMaterial(cid));

            for(int c = 0; c<nbCells; c++) {
                materials.insert(ma->getMaterial(neighborCellsIDs(c)));
            }

            int imat = 0;
            for(auto mat: materials) {

                gmds::math::Vector3d grad({0.,0.,0.});
                bool success = MaterialGradientComputation_computeGrad_leastsquare_onecell_onemat_2D(
                        cid,
                        mat,
                        fp,
                        c_C2C_byN,
                        varMidpoints,
                        grad
                );

                assert(success);
//                std::cout<<"cid "<<cid<<" "<<mat<<" "<<success<<std::endl;

                (*varGrads)[cid].m_matindex[imat] = mat;
                (*varGrads)[cid].m_grad[imat] = grad;

                imat++;
            }

        }
    };

/*----------------------------------------------------------------------------*/
    struct MaterialGradientComputation_computeGrad_leastsquare_allmat_onecell_2D {

        const kmds::GrowingView<kmds::TCellID> *cellIDsAccessor;

        kmds::Mesh *mesh;
        const kmds::Connectivity *c_C2C_byN;
        const elg3d::FracPres *fp;
        const int nbMat;
        const kmds::Variable<gmds::math::Point>* varMidpoints;
        kmds::Variable<MaterialGradientComputation_midcellGradients>* varGrads;


        MaterialGradientComputation_computeGrad_leastsquare_allmat_onecell_2D(const kmds::GrowingView<kmds::TCellID>* ACellIDsAccessor_,
                                                                              const kmds::Connectivity* c_C2C_byN_,
                                                                              const elg3d::FracPres* fp_,
                                                                              const int ANbMat_,
                                                                              const kmds::Variable<gmds::math::Point>* AVarMidpoints_,
                                                                              kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads_
        )
                : cellIDsAccessor(ACellIDsAccessor_)
                , c_C2C_byN(c_C2C_byN_)
                , fp(fp_)
                , nbMat(ANbMat_)
                , varMidpoints(AVarMidpoints_)
                , varGrads(AVarGrads_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int cid = cellIDsAccessor->get(i);

            // list the neighboring cells (neighbors by node)
            Kokkos::View<kmds::TCellID *> neighborCellsIDs;
            c_C2C_byN->get(cid, neighborCellsIDs);
            int nbCells = neighborCellsIDs.size();

            for(int imat=0; imat<nbMat; imat++) {

                gmds::math::Vector3d grad({0.,0.,0.});
                bool success = MaterialGradientComputation_computeGrad_leastsquare_onecell_onemat_2D(
                        cid,
                        imat,
                        fp,
                        c_C2C_byN,
                        varMidpoints,
                        grad
                );

                assert(success);
//                std::cout<<"cid "<<cid<<" "<<mat<<" "<<success<<std::endl;

                (*varGrads)[cid].m_matindex[imat] = imat;
                (*varGrads)[cid].m_grad[imat] = grad;
            }

        }
    };

    /*----------------------------------------------------------------------------*/
    struct MaterialGradientComputation_computeGrad_leastsquare_allmat_onecell_3D {

        const kmds::GrowingView<kmds::TCellID> *cellIDsAccessor;

        kmds::Mesh *mesh;
        const kmds::Connectivity *c_C2C_byN;
        const elg3d::FracPres *fp;
        const int nbMat;
        const kmds::Variable<gmds::math::Point>* varMidpoints;
        kmds::Variable<MaterialGradientComputation_midcellGradients>* varGrads;


        MaterialGradientComputation_computeGrad_leastsquare_allmat_onecell_3D(const kmds::GrowingView<kmds::TCellID>* ACellIDsAccessor_,
                                                                              const kmds::Connectivity* c_C2C_byN_,
                                                                              const elg3d::FracPres* fp_,
                                                                              const int ANbMat_,
                                                                              const kmds::Variable<gmds::math::Point>* AVarMidpoints_,
                                                                              kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads_
        )
                : cellIDsAccessor(ACellIDsAccessor_)
                , c_C2C_byN(c_C2C_byN_)
                , fp(fp_)
                , nbMat(ANbMat_)
                , varMidpoints(AVarMidpoints_)
                , varGrads(AVarGrads_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int cid = cellIDsAccessor->get(i);

            // list the neighboring cells (neighbors by node)
            Kokkos::View<kmds::TCellID *> neighborCellsIDs;
            c_C2C_byN->get(cid, neighborCellsIDs);
            int nbCells = neighborCellsIDs.size();

            for(int imat=0; imat<nbMat; imat++) {

                gmds::math::Vector3d grad({0.,0.,0.});
                bool success = MaterialGradientComputation_computeGrad_leastsquare_onecell_onemat_3D(
                        cid,
                        imat,
                        fp,
                        c_C2C_byN,
                        varMidpoints,
                        grad
                );

                assert(success);
//                std::cout<<"cid "<<cid<<" "<<mat<<" "<<success<<std::endl;

                (*varGrads)[cid].m_matindex[imat] = imat;
                (*varGrads)[cid].m_grad[imat] = grad;
            }

        }
    };

/*----------------------------------------------------------------------------*/
    bool
    MaterialGradientComputation_computeGrad_leastsquare_onecell_onemat_2D(const kmds::TCellID Acid,
                                                                          const int AMat,
                                                                          const elg3d::FracPres* Afp,
                                                                          const kmds::Connectivity* c_C2C_byN,
                                                                          const kmds::Variable<gmds::math::Point>* AVarMidpoints,
                                                                          gmds::math::Vector& AGrad)
    {
        gmds::math::Point pt = (*AVarMidpoints)[Acid];
        kmds::TFloat fp = Afp->getFracPres(AMat, Acid);

        Kokkos::View<kmds::TCellID *> cellIDs;
        c_C2C_byN->get(Acid, cellIDs);
        int nbCells = cellIDs.size();

        kmds::TCoord a00(0.);
        kmds::TCoord a01(0.);
        kmds::TCoord a02(0.);
        kmds::TCoord a10(0.);
        kmds::TCoord a11(0.);
        kmds::TCoord a12(0.);
        kmds::TCoord a20(0.);
        kmds::TCoord a21(0.);
        kmds::TCoord a22(0.);

        kmds::TFloat b0(0.);
        kmds::TFloat b1(0.);
        kmds::TFloat b2(0.);

        // compute weight for each cell
        // TODO investigate how the weight is defined
        kmds::TCoord weight_total = 0.;
        kmds::TCoord weights[nbCells];

        for(int icell=0; icell<nbCells; icell++) {
            kmds::TCellID cell = cellIDs[icell];

            weights[icell] = pt.distance2((*AVarMidpoints)[cell]);
            weight_total += weights[icell];
        }

        kmds::TCoord eps = weight_total / (100.0*nbCells);
        weight_total = 0.0;
        for (int ii = 0; ii< nbCells; ii++)
        {
            weights[ii] = 1.0 / (weights[ii] + eps);
            weight_total += weights[ii];
        }

        for(int icell=0; icell<nbCells; icell++) {
            kmds::TCellID cell = cellIDs[icell];

            kmds::TCoord dx = (*AVarMidpoints)[cell].X() - pt.X();
            kmds::TCoord dy = (*AVarMidpoints)[cell].Y() - pt.Y();

            kmds::TCoord w = weights[icell] / weight_total;


            a00 += w * dx * dx;
            a01 += w * dx * dy;

            a10 += w * dx * dy;
            a11 += w * dy * dy;


            kmds::TFloat dfpcell = w * (Afp->getFracPres(AMat, cell) - fp);

            b0 += dx * dfpcell;
            b1 += dy * dfpcell;
        }

        Eigen::Matrix2d Amatrix;
        Eigen::Vector2d Bvec;
        Eigen::Vector2d Xvec;

        Amatrix(0, 0) = a00;
        Amatrix(0, 1) = a01;

        Amatrix(1, 0) = a10;
        Amatrix(1, 1) = a11;


        // TODO check whether the system has no solution in which case return false
        kmds::TFloat det = a00*a11 - a10*a01;

//        std::cout<<"det "<<det<<std::endl;
//        std::cout<<Amatrix(0, 0)<<" ";
//        std::cout<<Amatrix(0, 1)<<" ";
//        std::cout<<Amatrix(0, 2)<<std::endl;
//
//        std::cout<<Amatrix(1, 0)<<" ";
//        std::cout<<Amatrix(1, 1)<<" ";
//        std::cout<<Amatrix(1, 2)<<std::endl;
//
//        std::cout<<Amatrix(2, 0)<<" ";
//        std::cout<<Amatrix(2, 1)<<" ";
//        std::cout<<Amatrix(2, 2)<<std::endl;

        if(det == 0.) {
            AGrad.setXYZ(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF);
            return false;
        }

        Bvec(0) = b0;
        Bvec(1) = b1;

        Eigen::LDLT<Eigen::Matrix2d> chol(Amatrix);
        Xvec = chol.solve(Bvec);

        AGrad.setXYZ(Xvec(0), Xvec(1),0.);
        AGrad.normalize();

//        if(Acid == 3 && AMat == 2) {
//            std::cout<<"Acid "<<Acid<<std::endl;
//            std::cout<<"line0 "<<a00 * Xvec(0) + a01 * Xvec(1) <<" b "<<b0<<std::endl;
//            std::cout<<"line0 "<<a10 * Xvec(0) + a11 * Xvec(1) <<" b "<<b1<<std::endl;
//
//            std::cout<<Xvec(0)<<" "<<Xvec(1)<<std::endl;
//            std::cout<< AGrad<<std::endl;
//        }
//
//        if(Acid == 4 && AMat == 2) {
//            std::cout<<"Acid "<<Acid<<std::endl;
//            std::cout<<"line0 "<<a00 * Xvec(0) + a01 * Xvec(1) <<" b "<<b0<<std::endl;
//            std::cout<<"line0 "<<a10 * Xvec(0) + a11 * Xvec(1) <<" b "<<b1<<std::endl;
//
//            std::cout<<Xvec(0)<<" "<<Xvec(1)<<std::endl;
//            std::cout<< AGrad<<std::endl;
//        }
//
//        if(Acid == 5 && AMat == 2) {
//            std::cout<<"Acid "<<Acid<<std::endl;
//            std::cout<<"line0 "<<a00 * Xvec(0) + a01 * Xvec(1) <<" b "<<b0<<std::endl;
//            std::cout<<"line0 "<<a10 * Xvec(0) + a11 * Xvec(1) <<" b "<<b1<<std::endl;
//
//            std::cout<<Xvec(0)<<" "<<Xvec(1)<<std::endl;
//            std::cout<< AGrad<<std::endl;
//        }

        return true;
    }
/*----------------------------------------------------------------------------*/
    void
    MaterialGradientComputation_leastsquare_2D(const kmds::GrowingView<kmds::TCellID>* ACellIDsAccessor,
                                               const kmds::Connectivity* c_C2C_byN,
                                               const elg3d::FracPres* Afp,
                                               const elg3d::MaterialAssignment *Ama,
                                               const kmds::Variable<gmds::math::Point>* AVarMidpoints,
                                               kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads)
    {
        kmds::TSize nbCells = ACellIDsAccessor->getNbElems();

        // compute the gradient at mid-cells
        // TODO only compute the gradient in cells near the interface
        Kokkos::parallel_for(nbCells, MaterialGradientComputation_computeGrad_leastsquare_onecell_2D(ACellIDsAccessor, c_C2C_byN, Afp, Ama, AVarMidpoints, AVarGrads));

    }

/*----------------------------------------------------------------------------*/
    void
    MaterialGradientComputation_leastsquare_allmat_2D(const kmds::GrowingView<kmds::TCellID>* ACellIDsAccessor,
                                                      const kmds::Connectivity* c_C2C_byN,
                                                      const elg3d::FracPres* Afp,
                                                      const kmds::Variable<gmds::math::Point>* AVarMidpoints,
                                                      kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads)
    {
        const kmds::TSize nbCells = ACellIDsAccessor->getNbElems();

        const kmds::TSize nbMat = Afp->getNbMaterials();
        if(nbMat > MaterialGradientComputation_MAXNBMATPERCELLBALL) {
            std::cout<<"MaterialGradientComputation_leastsquare_allmat_2D too many materials (" <<nbMat<<" compared to MAXNBMATPERCELLBALL ("<<MaterialGradientComputation_MAXNBMATPERCELLBALL<< ")"<<std::endl;
        }

        // compute the gradient at mid-cells
        Kokkos::parallel_for(nbCells, MaterialGradientComputation_computeGrad_leastsquare_allmat_onecell_2D(ACellIDsAccessor, c_C2C_byN, Afp, nbMat, AVarMidpoints, AVarGrads));

    }

/*----------------------------------------------------------------------------*/
    struct MaterialGradientComputation_verifInit_3D {

        const kmds::GrowingView<kmds::TCellID> *cellIDsAccessor;
        kmds::Variable<MaterialGradientComputation_midcellGradients>* varGrads;


        MaterialGradientComputation_verifInit_3D(const kmds::GrowingView<kmds::TCellID>* ACellIDsAccessor_,
                                                                       kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads_
        )
                : cellIDsAccessor(ACellIDsAccessor_)
                , varGrads(AVarGrads_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int cid = cellIDsAccessor->get(i);


            for(int imat=0; imat<MaterialGradientComputation_MAXNBMATPERCELLBALL; imat++) {

                assert((*varGrads)[cid].m_matindex[imat] == -1);
                assert((*varGrads)[cid].m_grad[imat] == gmds::math::Vector3d ({-HUGE_VALF, -HUGE_VALF, -HUGE_VALF}));
            }
        }
    };
/*----------------------------------------------------------------------------*/
    struct MaterialGradientComputation_computeGrad_leastsquare_onecell_3D {

        const kmds::GrowingView<kmds::TCellID> *cellIDsAccessor;

        kmds::Mesh *mesh;
        const kmds::Connectivity *c_C2C_byN;
        const elg3d::FracPres *fp;
        const elg3d::MaterialAssignment *ma;
        const kmds::Variable<gmds::math::Point>* varMidpoints;
        kmds::Variable<MaterialGradientComputation_midcellGradients>* varGrads;


        MaterialGradientComputation_computeGrad_leastsquare_onecell_3D(const kmds::GrowingView<kmds::TCellID>* ACellIDsAccessor_,
                                                             const kmds::Connectivity* c_C2C_byN_,
                                                             const elg3d::FracPres* fp_,
                                                             const elg3d::MaterialAssignment *ma_,
                                                             const kmds::Variable<gmds::math::Point>* AVarMidpoints_,
                                                             kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads_
        )
                : cellIDsAccessor(ACellIDsAccessor_)
                , c_C2C_byN(c_C2C_byN_)
                , fp(fp_)
                , ma(ma_)
                , varMidpoints(AVarMidpoints_)
                , varGrads(AVarGrads_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int cid = cellIDsAccessor->get(i);

            // list the neighboring cells (neighbors by node)
            Kokkos::View<kmds::TCellID *> neighborCellsIDs;
            c_C2C_byN->get(cid, neighborCellsIDs);
            int nbCells = neighborCellsIDs.size();


            // list the materials we will compute the gradient of
            std::set<int> materials;
            materials.insert(ma->getMaterial(cid));

            for(int c = 0; c<nbCells; c++) {
                kmds::TCellID id_tmp = neighborCellsIDs(c);
                assert(ma->getMaterial(id_tmp)>=0);
                materials.insert(ma->getMaterial(neighborCellsIDs(c)));
            }

            int imat = 0;
            for(auto mat: materials) {

                gmds::math::Vector3d grad({0.,0.,0.});
                bool success = MaterialGradientComputation_computeGrad_leastsquare_onecell_onemat_3D(
                        cid,
                        mat,
                        fp,
                        c_C2C_byN,
                        varMidpoints,
                        grad
                );

                assert(success);

                (*varGrads)[cid].m_matindex[imat] = mat;
                (*varGrads)[cid].m_grad[imat] = grad;

                imat++;
            }

        }
    };
/*----------------------------------------------------------------------------*/
    bool
    MaterialGradientComputation_computeGrad_leastsquare_onecell_onemat_3D(const kmds::TCellID Acid,
                                                                const int AMat,
                                                                const elg3d::FracPres* Afp,
                                                                const kmds::Connectivity* c_C2C_byN,
                                                                const kmds::Variable<gmds::math::Point>* AVarMidpoints,
                                                                gmds::math::Vector& AGrad)
    {
        gmds::math::Point pt = (*AVarMidpoints)[Acid];
        kmds::TFloat fp = Afp->getFracPres(AMat, Acid);

        Kokkos::View<kmds::TCellID *> cellIDs;
        c_C2C_byN->get(Acid, cellIDs);
        int nbCells = cellIDs.size();

        kmds::TCoord a00(0.);
        kmds::TCoord a01(0.);
        kmds::TCoord a02(0.);
        kmds::TCoord a10(0.);
        kmds::TCoord a11(0.);
        kmds::TCoord a12(0.);
        kmds::TCoord a20(0.);
        kmds::TCoord a21(0.);
        kmds::TCoord a22(0.);

        kmds::TFloat b0(0.);
        kmds::TFloat b1(0.);
        kmds::TFloat b2(0.);

        // compute weight for each cell
        // TODO investigate how the weight is defined
        kmds::TCoord weight_total = 0.;
        kmds::TCoord weights[nbCells];

        for(int icell=0; icell<nbCells; icell++) {
            kmds::TCellID cell = cellIDs[icell];

            weights[icell] = pt.distance2((*AVarMidpoints)[cell]);
            weight_total += weights[icell];
        }

        kmds::TCoord eps = weight_total / (100.0*nbCells);
        weight_total = 0.0;
        for (int ii = 0; ii< nbCells; ii++)
        {
            weights[ii] = 1.0 / (weights[ii] + eps);
            weight_total += weights[ii];
        }

        for(int icell=0; icell<nbCells; icell++) {
            kmds::TCellID cell = cellIDs[icell];

            kmds::TCoord dx = (*AVarMidpoints)[cell].X() - pt.X();
            kmds::TCoord dy = (*AVarMidpoints)[cell].Y() - pt.Y();
            kmds::TCoord dz = (*AVarMidpoints)[cell].Z() - pt.Z();

            kmds::TCoord w = weights[icell] / weight_total;


            a00 += w * dx * dx;
            a01 += w * dx * dy;
            a02 += w * dx * dz;

            a10 += w * dx * dy;
            a11 += w * dy * dy;
            a12 += w * dy * dz;

            a20 += w * dz * dx;
            a21 += w * dz * dy;
            a22 += w * dz * dz;

            kmds::TFloat dfpcell = w * (Afp->getFracPres(AMat, cell) - fp);

            b0 += dx * dfpcell;
            b1 += dy * dfpcell;
            b2 += dz * dfpcell;

        }

        Eigen::Matrix3d Amatrix;
        Eigen::Vector3d Bvec;
        Eigen::Vector3d Xvec;

        Amatrix(0, 0) = a00;
        Amatrix(0, 1) = a01;
        Amatrix(0, 2) = a02;

        Amatrix(1, 0) = a10;
        Amatrix(1, 1) = a11;
        Amatrix(1, 2) = a12;

        Amatrix(2, 0) = a20;
        Amatrix(2, 1) = a21;
        Amatrix(2, 2) = a22;

        // TODO check whether the system has no solution in which case return false
        kmds::TFloat det = a00*a11*a22 + a01*a12*a20 + a02*a10*a21 - (a20*a11*a02 + a21*a12*a00 + a22*a10*a01);


        if(det == 0.) {
            AGrad.setXYZ(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF);
            return false;
        }

        Bvec(0) = b0;
        Bvec(1) = b1;
        Bvec(2) = b2;

        Eigen::LDLT<Eigen::Matrix3d> chol(Amatrix);
        Xvec = chol.solve(Bvec);

        AGrad.setXYZ(Xvec(0), Xvec(1),Xvec(2));
        AGrad.normalize();

        if(std::isnan(Xvec(0)) ||
           std::isnan(Xvec(1)) ||
           std::isnan(Xvec(2)) ||
           std::isinf(Xvec(0)) ||
           std::isinf(Xvec(1)) ||
           std::isinf(Xvec(2))) {
            AGrad.setXYZ(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF);
            return false;
        }

        return true;
    }
/*----------------------------------------------------------------------------*/
    void
    MaterialGradientComputation_leastsquare_3D(const kmds::GrowingView<kmds::TCellID>* ACellIDsAccessor,
                                               const kmds::Connectivity* c_C2C_byN,
                                               const elg3d::FracPres* Afp,
                                               const elg3d::MaterialAssignment *Ama,
                                               const kmds::Variable<gmds::math::Point>* AVarMidpoints,
                                               kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads)
    {
        kmds::TSize nbCells = ACellIDsAccessor->getNbElems();


        std::cout<<"MaterialGradientComputation_leastsquare_3D begin verif"<<std::endl;

        Kokkos::parallel_for(nbCells, MaterialGradientComputation_verifInit_3D(ACellIDsAccessor, AVarGrads));


        std::cout<<"MaterialGradientComputation_leastsquare_3D begin computation"<<std::endl;

        // compute the gradient at mid-cells
        // TODO only compute the gradient in cells near the interface
        Kokkos::parallel_for(nbCells, MaterialGradientComputation_computeGrad_leastsquare_onecell_3D(ACellIDsAccessor, c_C2C_byN, Afp, Ama, AVarMidpoints, AVarGrads));

    }

    /*----------------------------------------------------------------------------*/
    void
    MaterialGradientComputation_leastsquare_allmat_3D(const kmds::GrowingView<kmds::TCellID>* ACellIDsAccessor,
                                                      const kmds::Connectivity* c_C2C_byN,
                                                      const elg3d::FracPres* Afp,
                                                      const kmds::Variable<gmds::math::Point>* AVarMidpoints,
                                                      kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads)
    {
        const kmds::TSize nbCells = ACellIDsAccessor->getNbElems();

        const kmds::TSize nbMat = Afp->getNbMaterials();
        if(nbMat > MaterialGradientComputation_MAXNBMATPERCELLBALL) {
            std::cout<<"MaterialGradientComputation_leastsquare_allmat_3D too many materials (" <<nbMat<<" compared to MAXNBMATPERCELLBALL ("<<MaterialGradientComputation_MAXNBMATPERCELLBALL<< ")"<<std::endl;
        }

        // compute the gradient at mid-cells
        Kokkos::parallel_for(nbCells, MaterialGradientComputation_computeGrad_leastsquare_allmat_onecell_3D(ACellIDsAccessor, c_C2C_byN, Afp, nbMat, AVarMidpoints, AVarGrads));

    }

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
