/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    SmartLaplacian.cpp
 *  \author  legoff
 *  \date    05/11/2018
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/SmartLaplacian.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
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
/*----------------------------------------------------------------------------*/
namespace elg3d {


    /*----------------------------------------------------------------------------*/
    struct SmartLaplacian_nonFixedNodesSelection
    {

        const kmds::GrowingView<kmds::TCellID>* nodesIDs;
        const kmds::Variable<bool>* var;

        kmds::GrowingView<kmds::TCellID>* selection;

        SmartLaplacian_nonFixedNodesSelection(

                const kmds::GrowingView<kmds::TCellID>* nodesIDs_,

                const kmds::Variable<bool>* var_,
                kmds::GrowingView<kmds::TCellID>* selection_
        )
                : nodesIDs(nodesIDs_)
                , var(var_)
                , selection(selection_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int nid = nodesIDs->get(i);

            if(!(*var)[nid]) {
                selection->push_back(nid);
            }
        }
    };


    /*----------------------------------------------------------------------------*/
    struct SmartLaplacian_smoothNode_2D
    {
        const kmds::GrowingView<kmds::TCellID>* selection;
        kmds::Mesh* mesh;
        const kmds::Connectivity* c_N2C;
        const kmds::Connectivity* c_N2N;
        const kmds::Variable<std::uintptr_t>* varNodeGeomAssociation;
        kmds::Variable<double>* varQuality;


        SmartLaplacian_smoothNode_2D(

                const kmds::GrowingView<kmds::TCellID>* selection_,
                kmds::Mesh* mesh_,
                const kmds::Connectivity* c_N2C_,
                const kmds::Connectivity* c_N2N_,
                const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation_,
                kmds::Variable<double>* varQuality_
        )
                : selection(selection_)
                , mesh(mesh_)
                , c_N2C(c_N2C_)
                , c_N2N(c_N2N_)
                , varNodeGeomAssociation(AVarNodeGeomAssociation_)
                , varQuality(varQuality_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int nid = selection->get(i);

            Kokkos::View<kmds::TCellID*> nids;
            c_N2N->get(nid, nids);

            gmds::math::Point pt(0.,0.,0.);

            for(int i_n=0; i_n<nids.size(); i_n++) {
                kmds::TCoord xyz[3];
                mesh->getNodeLocation(nids(i_n), xyz[0], xyz[1], xyz[2]);

                pt = pt + gmds::math::Point(xyz[0], xyz[1], xyz[2]);
            }

            pt = 1./nids.size() * pt;

            if((*varNodeGeomAssociation)[nid] != reinterpret_cast<std::uintptr_t>(nullptr)) {
                reinterpret_cast<gmds::cad::GeomEntity *>((*varNodeGeomAssociation)[nid])->project(pt);
            }

            kmds::TCoord xyz_old[3];
            mesh->getNodeLocation(nid, xyz_old[0], xyz_old[1], xyz_old[2]);

            mesh->setNodeLocation(nid, pt.X(), pt.Y(), pt.Z());

            // check whether minJ worsens
            kmds::TCoord minJ_old = HUGE_VALF;
            kmds::TCoord minJ_new = HUGE_VALF;

            Kokkos::View<kmds::TCellID *> cells;
            c_N2C->get(nid, cells);

            kmds::TCoord scaledJacobians[cells.size()];

            for (int c = 0; c < cells.size(); c++) {
                kmds::TCellID cid = cells[c];
                minJ_old = std::min(minJ_old, (*varQuality)[cid]);

                kmds::Face f = mesh->getFace(cid);
                scaledJacobians[c] = f.scaledJacobian();
                minJ_new = std::min(minJ_new, scaledJacobians[c]);
            }

            // revert location
            if(minJ_new < minJ_old) {
                mesh->setNodeLocation(nid, xyz_old[0], xyz_old[1], xyz_old[2]);
            } else {
                // update quality
                for(int c=0; c<cells.size(); c++) {
                    (*varQuality)[cells[c]] = scaledJacobians[c];
                }
            }
        }
    };

    /*----------------------------------------------------------------------------*/
    struct SmartLaplacian_smoothNode_optim_2D
    {
        const kmds::GrowingView<kmds::TCellID>* selection;
        kmds::Mesh* mesh;
        const kmds::Connectivity* c_N2C;
        const kmds::Connectivity* c_N2N;
        const kmds::Variable<std::uintptr_t>* varNodeGeomAssociation;
        kmds::Variable<double>* varQuality;


        SmartLaplacian_smoothNode_optim_2D(

                const kmds::GrowingView<kmds::TCellID>* selection_,
                kmds::Mesh* mesh_,
                const kmds::Connectivity* c_N2C_,
                const kmds::Connectivity* c_N2N_,
                const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation_,
                kmds::Variable<double>* varQuality_
        )
                : selection(selection_)
                , mesh(mesh_)
                , c_N2C(c_N2C_)
                , c_N2N(c_N2N_)
                , varNodeGeomAssociation(AVarNodeGeomAssociation_)
                , varQuality(varQuality_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int nid = selection->get(i);

            Kokkos::View<kmds::TCellID*> nids;
            c_N2N->get(nid, nids);

            gmds::math::Point pt(0.,0.,0.);

            for(int i_n=0; i_n<nids.size(); i_n++) {
                kmds::TCoord xyz[3];
                mesh->getNodeLocation(nids(i_n), xyz[0], xyz[1], xyz[2]);

                pt = pt + gmds::math::Point(xyz[0], xyz[1], xyz[2]);
            }

            pt = 1./nids.size() * pt;

            if((*varNodeGeomAssociation)[nid] != reinterpret_cast<std::uintptr_t>(nullptr)) {
                reinterpret_cast<gmds::cad::GeomEntity *>((*varNodeGeomAssociation)[nid])->project(pt);
            }

            kmds::TCoord xyz_old[3];
            mesh->getNodeLocation(nid, xyz_old[0], xyz_old[1], xyz_old[2]);

            gmds::math::Point pt_old(xyz_old[0], xyz_old[1], xyz_old[2]);
            gmds::math::Point positions[8];
            positions[0] = pt_old;
            positions[7] = pt;

            for(int ipt=1; ipt<7; ipt++) {
                gmds::math::Point pt_tmp = (1. - ipt/7.) * pt_old + ((ipt)/7.) * pt;
                if((*varNodeGeomAssociation)[nid] != reinterpret_cast<std::uintptr_t>(nullptr)) {
                    reinterpret_cast<gmds::cad::GeomEntity *>((*varNodeGeomAssociation)[nid])->project(pt_tmp);
                }
                positions[ipt] = pt_tmp;
            }

            kmds::TCoord maxminJ_all = -HUGE_VALF;
            int index = -1;

            Kokkos::View<kmds::TCellID *> cells;
            c_N2C->get(nid, cells);
            kmds::TCoord scaledJacobians[8][cells.size()];


            for(int ipt=0; ipt<8; ipt++) {

//                mesh->setNodeLocation(nid, pt.X(), pt.Y(), pt.Z());
                mesh->setNodeLocation(nid, positions[ipt]);

                // check whether minJ worsens
                kmds::TCoord minJ_new = HUGE_VALF;

                for (int c = 0; c < cells.size(); c++) {
                    kmds::TCellID cid = cells[c];

                    kmds::Face f = mesh->getFace(cid);
                    scaledJacobians[ipt][c] = f.scaledJacobian();
                    minJ_new = std::min(minJ_new, scaledJacobians[ipt][c]);
                }

                if(maxminJ_all < minJ_new) {
                    maxminJ_all = minJ_new;
                    index = ipt;
                }
            }

            mesh->setNodeLocation(nid, positions[index]);
            for(int c=0; c<cells.size(); c++) {
                (*varQuality)[cells[c]] = scaledJacobians[index][c];
            }
        }
    };

    /*----------------------------------------------------------------------------*/
    void
    SmartLaplacian_smoothNode_interface_optim_2D(const kmds::TCellID AID,
                                                 kmds::Mesh* mesh,
                                                 const kmds::Connectivity* c_N2C,
                                                 const kmds::Connectivity* c_N2N,
                                                 const gmds::cad::FACSurface * ASurf,
                                                 const elg3d::FacetedSurfaceGeomServices* AGeomServices,
                                                 kmds::Variable<double>* varQuality)
        {

            int nid = AID;

            Kokkos::View<kmds::TCellID*> nids;
            c_N2N->get(nid, nids);

            gmds::math::Point pt(0.,0.,0.);

            for(int i_n=0; i_n<nids.size(); i_n++) {
                kmds::TCoord xyz[3];
                mesh->getNodeLocation(nids(i_n), xyz[0], xyz[1], xyz[2]);

                pt = pt + gmds::math::Point(xyz[0], xyz[1], xyz[2]);
            }

            pt = 1./nids.size() * pt;


            AGeomServices->project(ASurf, pt);

            kmds::TCoord xyz_old[3];
            mesh->getNodeLocation(nid, xyz_old[0], xyz_old[1], xyz_old[2]);

            gmds::math::Point pt_old(xyz_old[0], xyz_old[1], xyz_old[2]);
            gmds::math::Point positions[8];
            for(int ipt=0; ipt<8; ipt++) {
                gmds::math::Point pt_tmp = (1. - ipt/7.) * pt_old + ((ipt)/7.) * pt;

                AGeomServices->project(ASurf, pt_tmp);

                positions[ipt] = pt_tmp;
            }

            kmds::TCoord maxminJ_all = -HUGE_VALF;
            int index = -1;

            for(int ipt=0; ipt<8; ipt++) {

//                mesh->setNodeLocation(nid, pt.X(), pt.Y(), pt.Z());
                mesh->setNodeLocation(nid, positions[ipt]);

                // check whether minJ worsens
                kmds::TCoord minJ_old = HUGE_VALF;
                kmds::TCoord minJ_new = HUGE_VALF;

                Kokkos::View<kmds::TCellID *> cells;
                c_N2C->get(nid, cells);

                kmds::TCoord scaledJacobians[cells.size()];

                for (int c = 0; c < cells.size(); c++) {
                    kmds::TCellID cid = cells[c];
                    minJ_old = std::min(minJ_old, (*varQuality)[cid]);

                    kmds::Face f = mesh->getFace(cid);
                    scaledJacobians[c] = f.scaledJacobian();
                    minJ_new = std::min(minJ_new, scaledJacobians[c]);
                }


                // TODO : we always choose the barycenter; why does it seem better ?
//                if(maxminJ_all < minJ_new) {
                    maxminJ_all = minJ_new;
                    index = ipt;

                    for(int c=0; c<cells.size(); c++) {
                        (*varQuality)[cells[c]] = scaledJacobians[c];
                    }
//                }
            }

            mesh->setNodeLocation(nid, positions[index]);

        }

    /*----------------------------------------------------------------------------*/
    void
    SmartLaplacian_smoothNode_interface_optim_doubleGeom_2D(const kmds::TCellID AID,
                                                            kmds::Mesh* mesh,
                                                            const kmds::Connectivity* c_N2C,
                                                            const kmds::Connectivity* c_N2N,
                                                            const kmds::Variable <std::uintptr_t> *AGeomassoc_interface_pixels,
                                                            const kmds::Variable <std::uintptr_t> *AGeomassoc_interface_boundingbox,
                                                            const elg3d::FacetedSurfaceGeomServices* AGeomServicesSurfaces,
                                                            const elg3d::FacetedCurveGeomServices* AGeomServicesCurves,
                                                            kmds::Variable<double>* varQuality)
    {

        int nid = AID;

        Kokkos::View<kmds::TCellID*> nids;
        c_N2N->get(nid, nids);

        gmds::math::Point pt(0.,0.,0.);

        if(nids.size() != 0) {
            for (int i_n = 0; i_n < nids.size(); i_n++) {
                kmds::TCoord xyz[3];
                mesh->getNodeLocation(nids(i_n), xyz[0], xyz[1], xyz[2]);

                pt = pt + gmds::math::Point(xyz[0], xyz[1], xyz[2]);
            }

            pt = 1. / nids.size() * pt;

        } else {
            pt = mesh->getNodeLocation(nid);
        }

        // projection on pixelated model
        const std::uintptr_t ge_pixels_ptr =  (*AGeomassoc_interface_pixels)[nid];
        if(ge_pixels_ptr == reinterpret_cast<std::uintptr_t>(nullptr)) {

            // default projection if no info on specific geom entity
            AGeomServicesSurfaces->project(pt);
        } else {

            gmds::cad::GeomEntity *ge_pixels = reinterpret_cast<gmds::cad::GeomEntity *> (ge_pixels_ptr);

            switch(ge_pixels->getDim()) {
                case 0 : {
                    gmds::cad::FACPoint *vert = reinterpret_cast<gmds::cad::FACPoint *> (ge_pixels);
                    vert->project(pt);
                    break;
                }
                case 1 : {
                    gmds::cad::FACCurve *curv = reinterpret_cast<gmds::cad::FACCurve *> (ge_pixels);
                    AGeomServicesCurves->project(curv, pt);
                    break;
                }
                case 2 : {
                    gmds::cad::FACSurface *surf = reinterpret_cast<gmds::cad::FACSurface *> (ge_pixels);
                    AGeomServicesSurfaces->project(surf, pt);
                    break;
                }
                default:
                    exit(-1);
            }
        }

        // computational domain geom association
        const std::uintptr_t ge_bb_ptr =  (*AGeomassoc_interface_boundingbox)[nid];
        if(ge_bb_ptr != reinterpret_cast<std::uintptr_t>(nullptr)) {
            gmds::cad::GeomEntity *ge_bb = reinterpret_cast<gmds::cad::GeomEntity *> (ge_bb_ptr);

            ge_bb->project(pt);
        }

        kmds::TCoord xyz_old[3];
        mesh->getNodeLocation(nid, xyz_old[0], xyz_old[1], xyz_old[2]);

        gmds::math::Point pt_old(xyz_old[0], xyz_old[1], xyz_old[2]);
        gmds::math::Point positions[8];
        for(int ipt=0; ipt<8; ipt++) {
            gmds::math::Point pt_tmp = (1. - ipt/7.) * pt_old + ((ipt)/7.) * pt;

//            AGeomServices->project(ASurf, pt_tmp);

            positions[ipt] = pt_tmp;


            if(ge_pixels_ptr == reinterpret_cast<std::uintptr_t>(nullptr)) {

                // default projection if no info on specific geom entity
                AGeomServicesSurfaces->project(pt_tmp);
            } else {

                gmds::cad::GeomEntity *ge_pixels = reinterpret_cast<gmds::cad::GeomEntity *> (ge_pixels_ptr);

                switch(ge_pixels->getDim()) {
                    case 0 : {
                        gmds::cad::FACPoint *vert = reinterpret_cast<gmds::cad::FACPoint *> (ge_pixels);
                        vert->project(pt_tmp);
                        break;
                    }
                    case 1 : {
                        gmds::cad::FACCurve *curv = reinterpret_cast<gmds::cad::FACCurve *> (ge_pixels);
                        AGeomServicesCurves->project(curv, pt_tmp);
                        break;
                    }
                    case 2 : {
                        gmds::cad::FACSurface *surf = reinterpret_cast<gmds::cad::FACSurface *> (ge_pixels);
                        AGeomServicesSurfaces->project(surf, pt_tmp);
                        break;
                    }
                    default:
                        exit(-1);
                }
            }

            // computational domain geom association
            const std::uintptr_t ge_bb_ptr =  (*AGeomassoc_interface_boundingbox)[nid];
            if(ge_bb_ptr != reinterpret_cast<std::uintptr_t>(nullptr)) {
                gmds::cad::GeomEntity *ge_bb = reinterpret_cast<gmds::cad::GeomEntity *> (ge_bb_ptr);

                ge_bb->project(pt_tmp);
            }

        }

        kmds::TCoord maxminJ_all = -HUGE_VALF;
        int index = -1;

        for(int ipt=0; ipt<8; ipt++) {

//                mesh->setNodeLocation(nid, pt.X(), pt.Y(), pt.Z());
            mesh->setNodeLocation(nid, positions[ipt]);

            // check whether minJ worsens
            kmds::TCoord minJ_old = HUGE_VALF;
            kmds::TCoord minJ_new = HUGE_VALF;

            Kokkos::View<kmds::TCellID *> cells;
            c_N2C->get(nid, cells);

            kmds::TCoord scaledJacobians[cells.size()];

            for (int c = 0; c < cells.size(); c++) {
                kmds::TCellID cid = cells[c];
                minJ_old = std::min(minJ_old, (*varQuality)[cid]);

                kmds::Face f = mesh->getFace(cid);
                scaledJacobians[c] = f.scaledJacobian();
                minJ_new = std::min(minJ_new, scaledJacobians[c]);
            }


            // TODO : we always choose the barycenter; why does it seem better ?
//                if(maxminJ_all < minJ_new) {
            maxminJ_all = minJ_new;
            index = ipt;

            for(int c=0; c<cells.size(); c++) {
                (*varQuality)[cells[c]] = scaledJacobians[c];
            }
//                }
        }

        mesh->setNodeLocation(nid, positions[index]);

    }

    /*----------------------------------------------------------------------------*/
    struct SmartLaplacian_smoothNode_3D
    {
        const kmds::GrowingView<kmds::TCellID>* selection;
        kmds::Mesh* mesh;
        const kmds::Connectivity* c_N2C;
        const kmds::Connectivity* c_N2N;
        const kmds::Variable<std::uintptr_t>* varNodeGeomAssociation;
        kmds::Variable<double>* varQuality;


        SmartLaplacian_smoothNode_3D(

                const kmds::GrowingView<kmds::TCellID>* selection_,
                kmds::Mesh* mesh_,
                const kmds::Connectivity* c_N2C_,
                const kmds::Connectivity* c_N2N_,
                const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation_,
                kmds::Variable<double>* varQuality_
        )
                : selection(selection_)
                , mesh(mesh_)
                , c_N2C(c_N2C_)
                , c_N2N(c_N2N_)
                , varNodeGeomAssociation(AVarNodeGeomAssociation_)
                , varQuality(varQuality_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int nid = selection->get(i);

            Kokkos::View<kmds::TCellID*> nids;
            c_N2N->get(nid, nids);

            gmds::math::Point pt(0.,0.,0.);

            for(int i_n=0; i_n<nids.size(); i_n++) {
                kmds::TCoord xyz[3];
                mesh->getNodeLocation(nids(i_n), xyz[0], xyz[1], xyz[2]);

                pt = pt + gmds::math::Point(xyz[0], xyz[1], xyz[2]);
            }

            pt = 1./nids.size() * pt;

            if((*varNodeGeomAssociation)[nid] != reinterpret_cast<std::uintptr_t>(nullptr)) {
                reinterpret_cast<gmds::cad::GeomEntity *>((*varNodeGeomAssociation)[nid])->project(pt);
            }

            kmds::TCoord xyz_old[3];
            mesh->getNodeLocation(nid, xyz_old[0], xyz_old[1], xyz_old[2]);

            mesh->setNodeLocation(nid, pt.X(), pt.Y(), pt.Z());

            // check whether minJ worsens
            kmds::TCoord minJ_old = HUGE_VALF;
            kmds::TCoord minJ_new = HUGE_VALF;

            Kokkos::View<kmds::TCellID*> cells;
            c_N2C->get(nid, cells);

            kmds::TCoord scaledJacobians[cells.size()];

            for(int c=0; c<cells.size(); c++) {
                kmds::TCellID cid = cells[c];
                minJ_old = std::min(minJ_old, (*varQuality)[cid]);

                kmds::Region r = mesh->getRegion(cid);
                scaledJacobians[c] = r.scaledJacobian();
                minJ_new = std::min(minJ_new, scaledJacobians[c]);
            }

            // revert location
            if(minJ_new < minJ_old) {
                mesh->setNodeLocation(nid, xyz_old[0], xyz_old[1], xyz_old[2]);
            } else {
                // update quality
                for(int c=0; c<cells.size(); c++) {
                    (*varQuality)[cells[c]] = scaledJacobians[c];
                }
            }
        }
    };

    /*----------------------------------------------------------------------------*/
    struct SmartLaplacian_smoothNode_optim_3D
    {
        const kmds::GrowingView<kmds::TCellID>* selection;
        kmds::Mesh* mesh;
        const kmds::Connectivity* c_N2C;
        const kmds::Connectivity* c_N2N;
        const kmds::Variable<std::uintptr_t>* varNodeGeomAssociation;
        kmds::Variable<double>* varQuality;


        SmartLaplacian_smoothNode_optim_3D(

                const kmds::GrowingView<kmds::TCellID>* selection_,
                kmds::Mesh* mesh_,
                const kmds::Connectivity* c_N2C_,
                const kmds::Connectivity* c_N2N_,
                const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation_,
                kmds::Variable<double>* varQuality_
        )
                : selection(selection_)
                , mesh(mesh_)
                , c_N2C(c_N2C_)
                , c_N2N(c_N2N_)
                , varNodeGeomAssociation(AVarNodeGeomAssociation_)
                , varQuality(varQuality_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int nid = selection->get(i);

            Kokkos::View<kmds::TCellID*> nids;
            c_N2N->get(nid, nids);

            gmds::math::Point pt(0.,0.,0.);

            for(int i_n=0; i_n<nids.size(); i_n++) {
                kmds::TCoord xyz[3];
                mesh->getNodeLocation(nids(i_n), xyz[0], xyz[1], xyz[2]);

                pt = pt + gmds::math::Point(xyz[0], xyz[1], xyz[2]);
            }

            pt = 1./nids.size() * pt;

            if((*varNodeGeomAssociation)[nid] != reinterpret_cast<std::uintptr_t>(nullptr)) {
                reinterpret_cast<gmds::cad::GeomEntity *>((*varNodeGeomAssociation)[nid])->project(pt);
            }

            kmds::TCoord xyz_old[3];
            mesh->getNodeLocation(nid, xyz_old[0], xyz_old[1], xyz_old[2]);

            gmds::math::Point pt_old(xyz_old[0], xyz_old[1], xyz_old[2]);
            gmds::math::Point positions[8];
            positions[0] = pt_old;
            positions[7] = pt;

            for(int ipt=1; ipt<7; ipt++) {
                gmds::math::Point pt_tmp = (1. - ipt/7.) * pt_old + ((ipt)/7.) * pt;
                if((*varNodeGeomAssociation)[nid] != reinterpret_cast<std::uintptr_t>(nullptr)) {
                    reinterpret_cast<gmds::cad::GeomEntity *>((*varNodeGeomAssociation)[nid])->project(pt_tmp);
                }
                positions[ipt] = pt_tmp;
            }

            kmds::TCoord maxminJ_all = -HUGE_VALF;
            int index = -1;

            Kokkos::View<kmds::TCellID *> cells;
            c_N2C->get(nid, cells);

            kmds::TCoord scaledJacobians[8][cells.size()];

            for(int ipt=0; ipt<8; ipt++) {

//                mesh->setNodeLocation(nid, pt.X(), pt.Y(), pt.Z());
                mesh->setNodeLocation(nid, positions[ipt]);

                // check whether minJ worsens
                kmds::TCoord minJ_new = HUGE_VALF;

                for (int c = 0; c < cells.size(); c++) {
                    kmds::TCellID cid = cells[c];

                    kmds::Region r = mesh->getRegion(cid);
                    scaledJacobians[ipt][c] = r.scaledJacobian();
                    minJ_new = std::min(minJ_new, scaledJacobians[ipt][c]);
                }

                if(maxminJ_all < minJ_new) {
                    maxminJ_all = minJ_new;
                    index = ipt;
                }
            }

            mesh->setNodeLocation(nid, positions[index]);
            for(int c=0; c<cells.size(); c++) {
                (*varQuality)[cells[c]] = scaledJacobians[index][c];
            }
        }
    };

    /*----------------------------------------------------------------------------*/
    struct SmartLaplacian_smoothNode_interface_3D
    {
        const kmds::GrowingView<kmds::TCellID>* selection;
        kmds::Mesh* mesh;
        const kmds::Connectivity* c_N2C;
        const kmds::Connectivity* c_N2N;
        const kmds::Variable<std::uintptr_t>* varNodeGeomAssociation;
        const kmds::Variable<std::uintptr_t>* varNodeGeomInterface;
        kmds::Variable<bool>* varIsInterface;
        kmds::Variable<double>* varQuality;


        SmartLaplacian_smoothNode_interface_3D(

                const kmds::GrowingView<kmds::TCellID>* selection_,
                kmds::Mesh* mesh_,
                const kmds::Connectivity* c_N2C_,
                const kmds::Connectivity* c_N2N_,
                const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation_,
                const kmds::Variable<std::uintptr_t>* AVarNodeGeomInterface_,
                kmds::Variable<bool>* varIsInterface_,
                kmds::Variable<double>* varQuality_
        )
                : selection(selection_)
                , mesh(mesh_)
                , c_N2C(c_N2C_)
                , c_N2N(c_N2N_)
                , varNodeGeomAssociation(AVarNodeGeomAssociation_)
                , varNodeGeomInterface(AVarNodeGeomInterface_)
                , varIsInterface(varIsInterface_)
                , varQuality(varQuality_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            int nid = selection->get(i);

            bool isInterface = (*varIsInterface)[nid];

            Kokkos::View<kmds::TCellID*> nids;
            c_N2N->get(nid, nids);

            gmds::math::Point pt(0.,0.,0.);

            int nbContrib = 0;

            for(int i_n=0; i_n<nids.size(); i_n++) {

                if(isInterface) {

                    if((*varIsInterface)[nids[i_n]]) {
                        kmds::TCoord xyz[3];
                        mesh->getNodeLocation(nids(i_n), xyz[0], xyz[1], xyz[2]);

                        pt = pt + gmds::math::Point(xyz[0], xyz[1], xyz[2]);

                        nbContrib++;
                    }
                } else {
                    kmds::TCoord xyz[3];
                    mesh->getNodeLocation(nids(i_n), xyz[0], xyz[1], xyz[2]);

                    pt = pt + gmds::math::Point(xyz[0], xyz[1], xyz[2]);

                    nbContrib++;
                }


            }

            pt = 1./nbContrib * pt;

            if((*varNodeGeomInterface)[nid] != reinterpret_cast<std::uintptr_t>(nullptr)) {
                reinterpret_cast<gmds::cad::GeomEntity *>((*varNodeGeomInterface)[nid])->project(pt);
            }

            if((*varNodeGeomAssociation)[nid] != reinterpret_cast<std::uintptr_t>(nullptr)) {
                reinterpret_cast<gmds::cad::GeomEntity *>((*varNodeGeomAssociation)[nid])->project(pt);
            }

            kmds::TCoord xyz_old[3];
            mesh->getNodeLocation(nid, xyz_old[0], xyz_old[1], xyz_old[2]);

            mesh->setNodeLocation(nid, pt.X(), pt.Y(), pt.Z());

            // check whether minJ worsens
            kmds::TCoord minJ_old = HUGE_VALF;
            kmds::TCoord minJ_new = HUGE_VALF;

            Kokkos::View<kmds::TCellID*> cells;
            c_N2C->get(nid, cells);

            kmds::TCoord scaledJacobians[cells.size()];

            for(int c=0; c<cells.size(); c++) {
                kmds::TCellID cid = cells[c];
                minJ_old = std::min(minJ_old, (*varQuality)[cid]);

                kmds::Region r = mesh->getRegion(cid);
                scaledJacobians[c] = r.scaledJacobian();
                minJ_new = std::min(minJ_new, scaledJacobians[c]);
            }

            // revert location
            if(minJ_new < minJ_old) {
                mesh->setNodeLocation(nid, xyz_old[0], xyz_old[1], xyz_old[2]);

            } else {
                 // update quality
                for(int c=0; c<cells.size(); c++) {
                    (*varQuality)[cells[c]] = scaledJacobians[c];
                }
            }
        }
    };

    /*----------------------------------------------------------------------------*/
    void
    smartLaplacian_2D(const int nbIter,
                      kmds::Mesh* AMesh,
                      const kmds::Connectivity* c_N2F,
                      const kmds::Connectivity* c_N2N,
                      const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                      const kmds::GrowingView<kmds::TCellID>* AFixedNodes)
    {
        Kokkos::Timer timer;
        timer.reset();

        // necessary init
        const int nbNodes = AMesh->getNbNodes();
        kmds::GrowingView<kmds::TCellID> nodesIDs("NODES", nbNodes);
        AMesh->getNodeIDs_dummy(&nodesIDs);

        const int nbCells = AMesh->getNbFaces();
        kmds::GrowingView<kmds::TCellID> cellsIDs("CELLS", nbCells);
        AMesh->getFaceIDs_dummy(&cellsIDs);

        kmds::Variable<bool>* varfixed = AMesh->createVariable<bool>(false, kmds::KMDS_NODE, "fixedNodes");
        Kokkos::parallel_for(AFixedNodes->getNbElems(), KOKKOS_LAMBDA(const int i) { (*varfixed)[AFixedNodes->get(i)] = true;});


        kmds::Variable<double>* varquality = AMesh->createVariable<double>(-HUGE_VALF, kmds::KMDS_FACE, "cellquality");
        Kokkos::parallel_for(nbCells,
                             KOKKOS_LAMBDA(const int i) {

                                 kmds::TCellID cid = cellsIDs.get(i);
                                 kmds::Face f = AMesh->getFace(cid);
                                 (*varquality)[cellsIDs.get(i)] = f.scaledJacobian() ;

                             });

//        std::cout<<"smartLaplacian_2D preptime "<<timer.seconds()<<std::endl;
        timer.reset();


        kmds::GrowingView<kmds::TCellID> nonFixedNodes("NODES_NONFIXED", nodesIDs.getNbElems());
        Kokkos::parallel_for(nbNodes, SmartLaplacian_nonFixedNodesSelection(&nodesIDs, varfixed, &nonFixedNodes));


//        std::cout<<"smartLaplacian_2D nonFixedNodes "<<timer.seconds()<<std::endl;
        timer.reset();

        const int nbNonFixedNodes = nonFixedNodes.getNbElems();

        // build graph
        // WARNING we hard-coded 20 as the max number of neighbors
        kmds::Graph graph("NODES_GRAPH", nbNonFixedNodes, 40);

        Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> kmapNode2Vertex(nbNonFixedNodes);

        ManifoldDetection_buildGraph_N_N2F2N(&nonFixedNodes, AMesh, c_N2F, &graph, &kmapNode2Vertex);
//        std::cout<<"smartLaplacian_2D buildGraph_N_N2F2N "<<timer.seconds()<<std::endl;
        timer.reset();

        // extract independent set
        graph.buildColoring(smartLaplacian_NBMAXCOLORS);


//        std::cout<<"smartLaplacian_2D buildColoring "<<timer.seconds()<<std::endl;
        timer.reset();


        // convert back from graph vertex to nodes IDs

        struct SmartLaplacian_vertex2Node
        {
            kmds::GrowingView<kmds::TCellID>* selection;
            const kmds::GrowingView<kmds::TCellID>* selection_map;

            SmartLaplacian_vertex2Node(kmds::GrowingView<kmds::TCellID>* selection_,
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

        for(int icolor=0; icolor<smartLaplacian_NBMAXCOLORS; icolor++) {

            kmds::GrowingView<kmds::TCellID> *nodesNonFixed_indepset = graph.getColoring(icolor);

            Kokkos::parallel_for(nodesNonFixed_indepset->getNbElems(),
                                 SmartLaplacian_vertex2Node(nodesNonFixed_indepset, &nonFixedNodes));

//            std::cout<<"smartLaplacian_2D colors "<<icolor<<" "<<nodesNonFixed_indepset->getNbElems()<<std::endl;

            if(nodesNonFixed_indepset->getNbElems() > 0 ) {
                colorIndexMax = icolor;
            }
        }

//        std::cout<<"smartLaplacian_2D vertex2nodes "<<timer.seconds()<<std::endl;
        timer.reset();


        int iter = 0;
        while(iter < nbIter) {




//            // convert back from graph vertex to nodes IDs
//
//            struct SmartLaplacian_vertex2Node
//            {
//                kmds::GrowingView<kmds::TCellID>* selection;
//                kmds::GrowingView<kmds::TCellID>* selection_map;
//
//                SmartLaplacian_vertex2Node(kmds::GrowingView<kmds::TCellID>* selection_, kmds::GrowingView<kmds::TCellID>* selection_map_)
//                        : selection(selection_)
//                        , selection_map(selection_map_)
//                {
//                }
//
//                KOKKOS_INLINE_FUNCTION
//                void
//                operator()(int i) const {
//                    int nid = selection->get(i);
//
//                    selection->set(i, selection_map->get(nid));
//                }
//            };

            for(int icolor=0; icolor<=colorIndexMax; icolor++) {

                kmds::GrowingView<kmds::TCellID>* nodesNonFixed_indepset = graph.getColoring(icolor);

//                Kokkos::parallel_for(nodesNonFixed_indepset->getNbElems(),
//                                     SmartLaplacian_vertex2Node(nodesNonFixed_indepset, &nonFixedNodes));

                Kokkos::parallel_for(nodesNonFixed_indepset->getNbElems(),
                                     SmartLaplacian_smoothNode_2D(nodesNonFixed_indepset,
                                                                  AMesh,
                                                                  c_N2F,
                                                                  c_N2N,
                                                                  AVarNodeGeomAssociation,
                                                                  varquality));

//                SmartLaplacian_smoothNode_optim_2D(nodesNonFixed_indepset,
//                                                   AMesh,
//                                                   c_N2F,
//                                                   c_N2N,
//                                                   AVarNodeGeomAssociation,
//                                                   varquality));
            }

            iter++;
        }

        std::cout<<"smartLaplacian_2D afteriter "<<timer.seconds()<<std::endl;

        // clean-up temporary data
        AMesh->deleteVariable(kmds::KMDS_NODE, "fixedNodes");
        AMesh->deleteVariable(kmds::KMDS_FACE, "cellquality");

//        std::cout<<"smartLaplacian_2D after deletevariables"<<std::endl;

    }
    /*----------------------------------------------------------------------------*/
    void
    smartLaplacian_3D(const int nbIter,
                      kmds::Mesh* AMesh,
                      const kmds::Connectivity* c_N2R,
                      const kmds::Connectivity* c_N2N,
                      const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                      const kmds::GrowingView<kmds::TCellID>* AFixedNodes)
    {
        Kokkos::Timer timer;
        timer.reset();

        // necessary init
        const int nbNodes = AMesh->getNbNodes();
        kmds::GrowingView<kmds::TCellID> nodesIDs("NODES", nbNodes);
        AMesh->getNodeIDs_dummy(&nodesIDs);

        const int nbCells = AMesh->getNbRegions();
        kmds::GrowingView<kmds::TCellID> cellsIDs("CELLS", nbCells);
        AMesh->getRegionIDs_dummy(&cellsIDs);

        kmds::Variable<bool>* varfixed = AMesh->createVariable<bool>(false, kmds::KMDS_NODE, "fixedNodes");
        Kokkos::parallel_for(AFixedNodes->getNbElems(), KOKKOS_LAMBDA(const int i) { (*varfixed)[AFixedNodes->get(i)] = true;});


        kmds::Variable<double>* varquality = AMesh->createVariable<double>(-HUGE_VALF, kmds::KMDS_REGION, "cellquality");
        Kokkos::parallel_for(nbCells,
                             KOKKOS_LAMBDA(const int i) {

                                 kmds::TCellID cid = cellsIDs.get(i);
                                 kmds::Region r = AMesh->getRegion(cid);
                                 (*varquality)[cid] = r.scaledJacobian() ;

                             });


        std::cout<<"smartLaplacian_3D preptime "<<timer.seconds()<<std::endl;
        timer.reset();


        kmds::GrowingView<kmds::TCellID> nonFixedNodes("NODES_NONFIXED", nodesIDs.getNbElems());
        Kokkos::parallel_for(nbNodes, SmartLaplacian_nonFixedNodesSelection(&nodesIDs, varfixed, &nonFixedNodes));


        std::cout<<"smartLaplacian_3D nonFixedNodes "<<timer.seconds()<<std::endl;
        timer.reset();


        const int nbNonFixedNodes = nonFixedNodes.getNbElems();

        // build graph
        // WARNING we hard-coded 20 as the max number of neighbors
        kmds::Graph graph("NODES_GRAPH", nbNonFixedNodes, 40);

        Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> kmapNode2Vertex(nbNonFixedNodes);

        ManifoldDetection_buildGraph_N_N2R2N(&nonFixedNodes, AMesh, c_N2R, &graph, &kmapNode2Vertex);

//        std::cout<<"smartLaplacian_3D buildGraph_N_N2R2N "<<timer.seconds()<<std::endl;
        timer.reset();

        // extract independent set
        graph.buildColoring(smartLaplacian_NBMAXCOLORS);

//        std::cout<<"smartLaplacian_3D buildColoring "<<timer.seconds()<<std::endl;
        timer.reset();

        // convert back from graph vertex to nodes IDs

        struct SmartLaplacian_vertex2Node
        {
            kmds::GrowingView<kmds::TCellID>* selection;
            const kmds::GrowingView<kmds::TCellID>* selection_map;

            SmartLaplacian_vertex2Node(kmds::GrowingView<kmds::TCellID>* selection_,
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

        for(int icolor=0; icolor<smartLaplacian_NBMAXCOLORS; icolor++) {

            kmds::GrowingView<kmds::TCellID> *nodesNonFixed_indepset = graph.getColoring(icolor);

            Kokkos::parallel_for(nodesNonFixed_indepset->getNbElems(),
                                 SmartLaplacian_vertex2Node(nodesNonFixed_indepset, &nonFixedNodes));

            if(nodesNonFixed_indepset->getNbElems() > 0 ) {
                colorIndexMax = icolor;
            }
//            std::cout<<"smartLaplacian_3D colors "<<icolor<<" "<<nodesNonFixed_indepset->getNbElems()<<std::endl;
        }

//        std::cout<<"smartLaplacian_3D vertex2nodes "<<timer.seconds()<<std::endl;
        timer.reset();


        int iter = 0;
        while(iter < nbIter) {




//            // convert back from graph vertex to nodes IDs
//
//            struct SmartLaplacian_vertex2Node
//            {
//                kmds::GrowingView<kmds::TCellID>* selection;
//                kmds::GrowingView<kmds::TCellID>* selection_map;
//
//                SmartLaplacian_vertex2Node(kmds::GrowingView<kmds::TCellID>* selection_, kmds::GrowingView<kmds::TCellID>* selection_map_)
//                        : selection(selection_)
//                        , selection_map(selection_map_)
//                {
//                }
//
//                KOKKOS_INLINE_FUNCTION
//                void
//                operator()(int i) const {
//                    int nid = selection->get(i);
//
//                    selection->set(i, selection_map->get(nid));
//                }
//            };

//            for(int icolor=0; icolor<smartLaplacian_NBMAXCOLORS; icolor++) {
            for(int icolor=0; icolor<=colorIndexMax; icolor++) {

                kmds::GrowingView<kmds::TCellID>* nodesNonFixed_indepset = graph.getColoring(icolor);

//                Kokkos::parallel_for(nodesNonFixed_indepset->getNbElems(),
//                                     SmartLaplacian_vertex2Node(nodesNonFixed_indepset, &nonFixedNodes));



                Kokkos::parallel_for(nodesNonFixed_indepset->getNbElems(),
//                                     SmartLaplacian_smoothNode_3D(nodesNonFixed_indepset,
//                                                                        AMesh,
//                                                                        c_N2R,
//                                                                        c_N2N,
//                                                                        AVarNodeGeomAssociation,
//                                                                        varquality));
                SmartLaplacian_smoothNode_optim_3D(nodesNonFixed_indepset,
                                                   AMesh,
                                                   c_N2R,
                                                   c_N2N,
                                                   AVarNodeGeomAssociation,
                                                   varquality));

//                std::cout<<"smartLaplacian_3D colors "<<icolor<<" "<<nodesNonFixed_indepset->getNbElems()<<" "<<timer.seconds()<<std::endl;
//                timer.reset();
            }


            iter++;
        }

        std::cout<<"smartLaplacian_3D afteriter "<<timer.seconds()<<std::endl;

//        std::cout<<"smartLaplacian_3D after nbiter"<<std::endl;

        // clean-up temporary data
        AMesh->deleteVariable(kmds::KMDS_NODE, "fixedNodes");
        AMesh->deleteVariable(kmds::KMDS_REGION, "cellquality");

//        std::cout<<"smartLaplacian_3D after deletevariables"<<std::endl;

    }
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
    void
    smartLaplacian_interface_3D(const int nbIter,
                                kmds::Mesh* AMesh,
                                const kmds::Connectivity* c_N2R,
                                const kmds::Connectivity* c_N2N,
                                const kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation,
                                const kmds::Variable<std::uintptr_t>* AVarNodeGeomInterface,
                                const kmds::GrowingView<kmds::TCellID>* AFixedNodes)
    {
        Kokkos::Timer timer;
        timer.reset();

        // necessary init
        const int nbNodes = AMesh->getNbNodes();
        kmds::GrowingView<kmds::TCellID> nodesIDs("NODES", nbNodes);
        AMesh->getNodeIDs_dummy(&nodesIDs);

        const int nbCells = AMesh->getNbRegions();
        kmds::GrowingView<kmds::TCellID> cellsIDs("CELLS", nbCells);
        AMesh->getRegionIDs_dummy(&cellsIDs);

        kmds::Variable<bool>* varfixed = AMesh->createVariable<bool>(false, kmds::KMDS_NODE, "fixedNodes");
//        Kokkos::parallel_for(AFixedNodes->getNbElems(), KOKKOS_LAMBDA(const int i) { (*varfixed)[AFixedNodes->get(i)] = true;});


        kmds::Variable<double>* varquality = AMesh->createVariable<double>(-HUGE_VALF, kmds::KMDS_REGION, "cellquality");
        Kokkos::parallel_for(nbCells,
                             KOKKOS_LAMBDA(const int i) {

                                 kmds::TCellID cid = cellsIDs.get(i);
                                 kmds::Region r = AMesh->getRegion(cid);
                                 (*varquality)[cid] = r.scaledJacobian() ;

                             });


        std::cout<<"smartLaplacian_3D preptime "<<timer.seconds()<<std::endl;
        timer.reset();


//        kmds::GrowingView<kmds::TCellID> nonFixedNodes("NODES_NONFIXED", nodesIDs.getNbElems());
//        Kokkos::parallel_for(nbNodes, SmartLaplacian_nonFixedNodesSelection(&nodesIDs, varfixed, &nonFixedNodes));
        kmds::GrowingView<kmds::TCellID> nonFixedNodes("NODES", nbNodes);
        AMesh->getNodeIDs_dummy(&nonFixedNodes);



        std::cout<<"smartLaplacian_3D nonFixedNodes "<<timer.seconds()<<std::endl;
        timer.reset();


        const int nbNonFixedNodes = nonFixedNodes.getNbElems();

        // build graph
        // WARNING we hard-coded 20 as the max number of neighbors
        kmds::Graph graph("NODES_GRAPH", nbNonFixedNodes, 40);

        Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> kmapNode2Vertex(nbNonFixedNodes);

        ManifoldDetection_buildGraph_N_N2R2N(&nonFixedNodes, AMesh, c_N2R, &graph, &kmapNode2Vertex);

//        std::cout<<"smartLaplacian_3D buildGraph_N_N2R2N "<<timer.seconds()<<std::endl;
        timer.reset();

        // extract independent set
        graph.buildColoring(smartLaplacian_NBMAXCOLORS);

//        std::cout<<"smartLaplacian_3D buildColoring "<<timer.seconds()<<std::endl;
        timer.reset();

        // convert back from graph vertex to nodes IDs

        struct SmartLaplacian_vertex2Node
        {
            kmds::GrowingView<kmds::TCellID>* selection;
            const kmds::GrowingView<kmds::TCellID>* selection_map;

            SmartLaplacian_vertex2Node(kmds::GrowingView<kmds::TCellID>* selection_,
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

        for(int icolor=0; icolor<smartLaplacian_NBMAXCOLORS; icolor++) {

            kmds::GrowingView<kmds::TCellID> *nodesNonFixed_indepset = graph.getColoring(icolor);

            Kokkos::parallel_for(nodesNonFixed_indepset->getNbElems(),
                                 SmartLaplacian_vertex2Node(nodesNonFixed_indepset, &nonFixedNodes));

            if(nodesNonFixed_indepset->getNbElems() > 0 ) {
                colorIndexMax = icolor;
            }
//            std::cout<<"smartLaplacian_3D colors "<<icolor<<" "<<nodesNonFixed_indepset->getNbElems()<<std::endl;
        }

//        std::cout<<"smartLaplacian_3D vertex2nodes "<<timer.seconds()<<std::endl;
        timer.reset();


        kmds::Variable<bool>* varIsInterface = AMesh->createVariable<bool>(false, kmds::KMDS_NODE, "tmp_varIsInterface");
        Kokkos::parallel_for(AFixedNodes->getNbElems(),
                             KOKKOS_LAMBDA(const int i) {

                                 kmds::TCellID nid = AFixedNodes->get(i);
                                 (*varIsInterface)[nid] = true;

                             });


        int iter = 0;
        while(iter < nbIter) {




//            // convert back from graph vertex to nodes IDs
//
//            struct SmartLaplacian_vertex2Node
//            {
//                kmds::GrowingView<kmds::TCellID>* selection;
//                kmds::GrowingView<kmds::TCellID>* selection_map;
//
//                SmartLaplacian_vertex2Node(kmds::GrowingView<kmds::TCellID>* selection_, kmds::GrowingView<kmds::TCellID>* selection_map_)
//                        : selection(selection_)
//                        , selection_map(selection_map_)
//                {
//                }
//
//                KOKKOS_INLINE_FUNCTION
//                void
//                operator()(int i) const {
//                    int nid = selection->get(i);
//
//                    selection->set(i, selection_map->get(nid));
//                }
//            };

//            for(int icolor=0; icolor<smartLaplacian_NBMAXCOLORS; icolor++) {
            for(int icolor=0; icolor<=colorIndexMax; icolor++) {

                kmds::GrowingView<kmds::TCellID>* nodesNonFixed_indepset = graph.getColoring(icolor);

//                Kokkos::parallel_for(nodesNonFixed_indepset->getNbElems(),
//                                     SmartLaplacian_vertex2Node(nodesNonFixed_indepset, &nonFixedNodes));



                Kokkos::parallel_for(nodesNonFixed_indepset->getNbElems(),
                                     SmartLaplacian_smoothNode_interface_3D(nodesNonFixed_indepset,
                                                                            AMesh,
                                                                            c_N2R,
                                                                            c_N2N,
                                                                            AVarNodeGeomAssociation,
                                                                            AVarNodeGeomInterface,
                                                                            varIsInterface,
                                                                            varquality));

//                std::cout<<"smartLaplacian_3D colors "<<icolor<<" "<<nodesNonFixed_indepset->getNbElems()<<" "<<timer.seconds()<<std::endl;
//                timer.reset();
            }


            iter++;
        }

        std::cout<<"smartLaplacian_3D afteriter "<<timer.seconds()<<std::endl;

//        std::cout<<"smartLaplacian_3D after nbiter"<<std::endl;

        // clean-up temporary data
        AMesh->deleteVariable(kmds::KMDS_NODE, "fixedNodes");
        AMesh->deleteVariable(kmds::KMDS_REGION, "cellquality");

        AMesh->deleteVariable(kmds::KMDS_NODE, "tmp_varIsInterface");

//        std::cout<<"smartLaplacian_3D after deletevariables"<<std::endl;

    }

    /*----------------------------------------------------------------------------*/
    void
    smartLaplacian_interface_fromF_3D(const int nbIter,
                                      kmds::Mesh* AMesh,
                                      const kmds::Connectivity* c_N2F,
                                      const kmds::Connectivity* c_N2N,
                                      const gmds::cad::FACSurface * ASurf,
                                      const elg3d::FacetedSurfaceGeomServices* AGeomServices)
    {
        Kokkos::Timer timer;
        timer.reset();

        // necessary init
        const int nbNodes = AMesh->getNbNodes();
        kmds::GrowingView<kmds::TCellID> nodesIDs("NODES", nbNodes);
        AMesh->getNodeIDs_dummy(&nodesIDs);

        const int nbCells = AMesh->getNbFaces();
        kmds::GrowingView<kmds::TCellID> cellsIDs("CELLS", nbCells);
        AMesh->getFaceIDs_dummy(&cellsIDs);

        kmds::Variable<bool>* varfixed = AMesh->createVariable<bool>(false, kmds::KMDS_NODE, "fixedNodes");


        kmds::Variable<double>* varquality = AMesh->createVariable<double>(-HUGE_VALF, kmds::KMDS_FACE, "cellquality");
        Kokkos::parallel_for(nbCells,
                             KOKKOS_LAMBDA(const int i) {

                                 kmds::TCellID cid = cellsIDs.get(i);
                                 kmds::Face f = AMesh->getFace(cid);
                                 (*varquality)[cellsIDs.get(i)] = f.scaledJacobian() ;

                             });

//        std::cout<<"smartLaplacian_2D preptime "<<timer.seconds()<<std::endl;
        timer.reset();


        kmds::GrowingView<kmds::TCellID> nonFixedNodes("NODES_NONFIXED", nodesIDs.getNbElems());
        Kokkos::parallel_for(nbNodes, SmartLaplacian_nonFixedNodesSelection(&nodesIDs, varfixed, &nonFixedNodes));


//        std::cout<<"smartLaplacian_2D nonFixedNodes "<<timer.seconds()<<std::endl;
        timer.reset();

        const int nbNonFixedNodes = nonFixedNodes.getNbElems();

        // build graph
        // WARNING we hard-coded 20 as the max number of neighbors
        kmds::Graph graph("NODES_GRAPH", nbNonFixedNodes, 40);

        Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> kmapNode2Vertex(nbNonFixedNodes);

        ManifoldDetection_buildGraph_N_N2F2N(&nonFixedNodes, AMesh, c_N2F, &graph, &kmapNode2Vertex);
//        std::cout<<"smartLaplacian_2D buildGraph_N_N2F2N "<<timer.seconds()<<std::endl;
        timer.reset();

        // extract independent set
        graph.buildColoring(smartLaplacian_NBMAXCOLORS);


//        std::cout<<"smartLaplacian_2D buildColoring "<<timer.seconds()<<std::endl;
        timer.reset();


        // convert back from graph vertex to nodes IDs

        struct SmartLaplacian_vertex2Node
        {
            kmds::GrowingView<kmds::TCellID>* selection;
            const kmds::GrowingView<kmds::TCellID>* selection_map;

            SmartLaplacian_vertex2Node(kmds::GrowingView<kmds::TCellID>* selection_,
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

        for(int icolor=0; icolor<smartLaplacian_NBMAXCOLORS; icolor++) {

            kmds::GrowingView<kmds::TCellID> *nodesNonFixed_indepset = graph.getColoring(icolor);

            Kokkos::parallel_for(nodesNonFixed_indepset->getNbElems(),
                                 SmartLaplacian_vertex2Node(nodesNonFixed_indepset, &nonFixedNodes));

//            std::cout<<"smartLaplacian_2D colors "<<icolor<<" "<<nodesNonFixed_indepset->getNbElems()<<std::endl;

            if(nodesNonFixed_indepset->getNbElems() > 0 ) {
                colorIndexMax = icolor;
            }
        }

//        std::cout<<"smartLaplacian_2D vertex2nodes "<<timer.seconds()<<std::endl;
        timer.reset();


        int iter = 0;
        while(iter < nbIter) {




//            // convert back from graph vertex to nodes IDs
//
//            struct SmartLaplacian_vertex2Node
//            {
//                kmds::GrowingView<kmds::TCellID>* selection;
//                kmds::GrowingView<kmds::TCellID>* selection_map;
//
//                SmartLaplacian_vertex2Node(kmds::GrowingView<kmds::TCellID>* selection_, kmds::GrowingView<kmds::TCellID>* selection_map_)
//                        : selection(selection_)
//                        , selection_map(selection_map_)
//                {
//                }
//
//                KOKKOS_INLINE_FUNCTION
//                void
//                operator()(int i) const {
//                    int nid = selection->get(i);
//
//                    selection->set(i, selection_map->get(nid));
//                }
//            };

            for(int icolor=0; icolor<=colorIndexMax; icolor++) {

                kmds::GrowingView<kmds::TCellID> *nodesNonFixed_indepset = graph.getColoring(icolor);

                // TODO : cannot be called in parallel because gts is not threadsafe
                for (int i = 0; i < nodesNonFixed_indepset->getNbElems(); i++) {

                    kmds::TCellID nid = nodesNonFixed_indepset->get(i);

                    SmartLaplacian_smoothNode_interface_optim_2D(nid,
                                                                 AMesh,
                                                                 c_N2F,
                                                                 c_N2N,
                                                                 ASurf,
                                                                 AGeomServices,
                                                                 varquality);
                }
            }

            iter++;
        }

        std::cout<<"smartLaplacian_2D afteriter "<<timer.seconds()<<std::endl;

        // clean-up temporary data
        AMesh->deleteVariable(kmds::KMDS_NODE, "fixedNodes");
        AMesh->deleteVariable(kmds::KMDS_FACE, "cellquality");

//        std::cout<<"smartLaplacian_2D after deletevariables"<<std::endl;

    }

    /*----------------------------------------------------------------------------*/
    void
    smartLaplacian_interface_fromF_doubleGeom_3D(const int nbIter,
                                                 kmds::Mesh* AMesh,
                                                 const kmds::Connectivity* c_N2F,
                                                 const kmds::Connectivity* c_N2N,
                                                 const kmds::Variable <std::uintptr_t> *AGeomassoc_interface_pixels,
                                                 const kmds::Variable <std::uintptr_t> *AGeomassoc_interface_boundingbox,
                                                 const elg3d::FacetedSurfaceGeomServices* AGeomServicesSurfaces,
                                                 const elg3d::FacetedCurveGeomServices* AGeomServicesCurves)
    {
        Kokkos::Timer timer;
        timer.reset();

        // necessary init
        const int nbNodes = AMesh->getNbNodes();
        kmds::GrowingView<kmds::TCellID> nodesIDs("NODES", nbNodes);
        AMesh->getNodeIDs_dummy(&nodesIDs);

        const int nbCells = AMesh->getNbFaces();
        kmds::GrowingView<kmds::TCellID> cellsIDs("CELLS", nbCells);
        AMesh->getFaceIDs_dummy(&cellsIDs);

        kmds::Variable<bool>* varfixed = AMesh->createVariable<bool>(false, kmds::KMDS_NODE, "fixedNodes");


        kmds::Variable<double>* varquality = AMesh->createVariable<double>(-HUGE_VALF, kmds::KMDS_FACE, "cellquality");
        Kokkos::parallel_for(nbCells,
                             KOKKOS_LAMBDA(const int i) {

                                 kmds::TCellID cid = cellsIDs.get(i);
                                 kmds::Face f = AMesh->getFace(cid);
                                 (*varquality)[cellsIDs.get(i)] = f.scaledJacobian() ;

                             });

//        std::cout<<"smartLaplacian_2D preptime "<<timer.seconds()<<std::endl;
        timer.reset();


        kmds::GrowingView<kmds::TCellID> nonFixedNodes("NODES_NONFIXED", nodesIDs.getNbElems());
        Kokkos::parallel_for(nbNodes, SmartLaplacian_nonFixedNodesSelection(&nodesIDs, varfixed, &nonFixedNodes));


//        std::cout<<"smartLaplacian_2D nonFixedNodes "<<timer.seconds()<<std::endl;
        timer.reset();

        const int nbNonFixedNodes = nonFixedNodes.getNbElems();

        // build graph
        // WARNING we hard-coded 20 as the max number of neighbors
        kmds::Graph graph("NODES_GRAPH", nbNonFixedNodes, 40);

        Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> kmapNode2Vertex(nbNonFixedNodes);

        ManifoldDetection_buildGraph_N_N2F2N(&nonFixedNodes, AMesh, c_N2F, &graph, &kmapNode2Vertex);
//        std::cout<<"smartLaplacian_2D buildGraph_N_N2F2N "<<timer.seconds()<<std::endl;
        timer.reset();

        // extract independent set
        graph.buildColoring(smartLaplacian_NBMAXCOLORS);


//        std::cout<<"smartLaplacian_2D buildColoring "<<timer.seconds()<<std::endl;
        timer.reset();


        // convert back from graph vertex to nodes IDs

        struct SmartLaplacian_vertex2Node
        {
            kmds::GrowingView<kmds::TCellID>* selection;
            const kmds::GrowingView<kmds::TCellID>* selection_map;

            SmartLaplacian_vertex2Node(kmds::GrowingView<kmds::TCellID>* selection_,
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

        for(int icolor=0; icolor<smartLaplacian_NBMAXCOLORS; icolor++) {

            kmds::GrowingView<kmds::TCellID> *nodesNonFixed_indepset = graph.getColoring(icolor);

            Kokkos::parallel_for(nodesNonFixed_indepset->getNbElems(),
                                 SmartLaplacian_vertex2Node(nodesNonFixed_indepset, &nonFixedNodes));

//            std::cout<<"smartLaplacian_2D colors "<<icolor<<" "<<nodesNonFixed_indepset->getNbElems()<<std::endl;

            if(nodesNonFixed_indepset->getNbElems() > 0 ) {
                colorIndexMax = icolor;
            }
        }

//        std::cout<<"smartLaplacian_2D vertex2nodes "<<timer.seconds()<<std::endl;
        timer.reset();


        int iter = 0;
        while(iter < nbIter) {




//            // convert back from graph vertex to nodes IDs
//
//            struct SmartLaplacian_vertex2Node
//            {
//                kmds::GrowingView<kmds::TCellID>* selection;
//                kmds::GrowingView<kmds::TCellID>* selection_map;
//
//                SmartLaplacian_vertex2Node(kmds::GrowingView<kmds::TCellID>* selection_, kmds::GrowingView<kmds::TCellID>* selection_map_)
//                        : selection(selection_)
//                        , selection_map(selection_map_)
//                {
//                }
//
//                KOKKOS_INLINE_FUNCTION
//                void
//                operator()(int i) const {
//                    int nid = selection->get(i);
//
//                    selection->set(i, selection_map->get(nid));
//                }
//            };

            for(int icolor=0; icolor<=colorIndexMax; icolor++) {

                kmds::GrowingView<kmds::TCellID> *nodesNonFixed_indepset = graph.getColoring(icolor);

                // TODO : cannot be called in parallel because gts is not threadsafe
                for (int i = 0; i < nodesNonFixed_indepset->getNbElems(); i++) {

                    kmds::TCellID nid = nodesNonFixed_indepset->get(i);

                    SmartLaplacian_smoothNode_interface_optim_doubleGeom_2D(nid,
                                                                            AMesh,
                                                                            c_N2F,
                                                                            c_N2N,
                                                                            AGeomassoc_interface_pixels,
                                                                            AGeomassoc_interface_boundingbox,
                                                                            AGeomServicesSurfaces,
                                                                            AGeomServicesCurves,
                                                                            varquality);
                }
            }

            iter++;
        }

        std::cout<<"smartLaplacian_2D afteriter "<<timer.seconds()<<std::endl;

        // clean-up temporary data
        AMesh->deleteVariable(kmds::KMDS_NODE, "fixedNodes");
        AMesh->deleteVariable(kmds::KMDS_FACE, "cellquality");

//        std::cout<<"smartLaplacian_2D after deletevariables"<<std::endl;

    }

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
