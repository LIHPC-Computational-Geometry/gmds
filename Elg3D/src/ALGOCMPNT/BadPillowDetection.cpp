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
/** \file    BadPillowDetection.cpp
 *  \author  legoff
 *  \date    02/01/201
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/BadPillowDetection.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <Kokkos_UnorderedMap.hpp>

extern "C" {
#include "ELG3D/ALGOCMPNT/r3d.h"
}

#include <glpk.h>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/math/Triangle.h>
#include <gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
namespace elg3d {


    const kmds::TCoord BADPILLOWINGDETECTION_OPPOSITE_DOTPRODUCT_THRESHOLD = - 0.5;

    const kmds::TCoord BADPILLOWINGDETECTION_GLPK_THRESHOLD = 0.1;
    const kmds::TCoord BADPILLOWINGDETECTION_R3D_THRESHOLD = 0.2;
    const kmds::TCoord BADPILLOWINGDETECTION_ORTHO_THRESHOLD = 0.1;



    /*----------------------------------------------------------------------------*/
//    struct BadPillowDetection_isNodeBad {
//        const kmds::GrowingView<kmds::TCellID> *selection;
//        kmds::Mesh *mesh;
//        const kmds::Connectivity* c_N2F;
//        const kmds::Connectivity* c_F2R;
//        const elg3d::MaterialAssignment* ma;
//        kmds::Variable<std::uintptr_t> *varNodeGeomAssociation;
//        std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> >* node2newNodes;
//
//
//        BadPillowDetection_isNodeBad(
//
//                const kmds::GrowingView<kmds::TCellID> *selection_,
//                kmds::Mesh *mesh_,
//                const kmds::Connectivity* c_N2C_,
//                const elg3d::MaterialAssignment* ma_,
//                kmds::Variable<std::uintptr_t> *AVarNodeGeomAssociation_,
//                std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> >* node2newNodes_
//        )
//                : selection(selection_)
//                , mesh(mesh_)
//                , c_N2C(c_N2C_)
//                , ma(ma_)
//                , varNodeGeomAssociation(AVarNodeGeomAssociation_)
//                , node2newNodes(node2newNodes_)
//        {
//        }
//
//        KOKKOS_INLINE_FUNCTION
//        void
//        operator()(int i) const {
//            int nid = selection->get(i);
//
//            Kokkos::View<kmds::TCellID*> cells;
//            c_N2C->get(nid, cells);
//
//            std::set<int> materials;
//            for(int i_c=0; i_c<cells.size(); i_c++) {
//                materials.insert(ma->getMaterial(cells(i_c)));
//            }
//
//
//            kmds::TCellID id0 = mesh->addNodes(materials.size());
//            kmds::TCoord xyz[3];
//            mesh->getNodeLocation(nid, xyz[0], xyz[1], xyz[2]);
//
//            int offset = 0;
//            for(auto mat: materials) {
//
//                // compute the new node position
//                gmds::math::Point pt_old(xyz[0], xyz[1], xyz[2]);
//                gmds::math::Point pt_new(0., 0., 0.);
//                int nbContrib = 0;
//                for(int i_c=0; i_c<cells.size(); i_c++) {
//                    if(ma->getMaterial(cells(i_c)) == mat) {
//                        gmds::math::Point pt_midpoint = mesh->getFace(cells(i_c)).midpoint();
//                        gmds::math::Point pt_tmp(pt_old + (1/5.) * gmds::math::Vector(pt_old, pt_midpoint));
//                        pt_new = pt_new + pt_tmp;
//                        nbContrib++;
//                    }
//                }
//
//                pt_new = (1./(double) nbContrib) * pt_new;
//                if((*varNodeGeomAssociation)[nid] != reinterpret_cast<std::uintptr_t>(nullptr)) {
////                    reinterpret_cast<gmds::cad::GeomEntity *>((*varNodeGeomAssociation)[nid])->project(pt_new);
//                }
//
//
//                mesh->setNodeLocation(id0 + offset, pt_new.X(), pt_new.Y(), pt_new.Z());
//
////                mesh->setNodeLocation(id0 + offset, xyz[0], xyz[1], xyz[2]);
//
//                (*node2newNodes)[mat].insert(nid, id0 + offset);
//                (*varNodeGeomAssociation)[id0 + offset] = (*varNodeGeomAssociation)[nid];
//                offset++;
//            }
//
//        }
//
//    };

    /*----------------------------------------------------------------------------*/
    bool BadPillowDetection_IsNodeBad_3D(kmds::Mesh* AMesh,
                                         const kmds::TCellID ANodeID,
                                         const std::vector<kmds::TCellID> ACellIDs,
                                         const std::vector<std::vector<kmds::FakeFace> >& Affs,
                                         const std::vector<std::vector<bool> >& AffIsOutward,
                                         const std::vector<std::vector<gmds::math::Vector> >& AffNormals)
    {
        // if there is only one cell, it is always ok
        // WARNING : We assume it is also the case with two cells (if it is not the case
        // it would mean that we have a flat cell between those two)
        if(ACellIDs.size() <= 2) {
            return false;
        }

        // we build this set of regions FakeFaces and we compute how many there are considering unicity
        std::map<kmds::FakeFace, int> nbFakeFaces;
        std::map<kmds::FakeFace, std::pair<bool, gmds::math::Vector> > pair_orientation_normals;

        for(int ic=0; ic<ACellIDs.size(); ic++) {

            kmds::TCellID cellid = ACellIDs[ic];

            for(int iff=0; iff<Affs[cellid].size(); iff++) {

                if(Affs[cellid][iff].hasNode(ANodeID)) {

                    nbFakeFaces[Affs[cellid][iff]]++;
                    pair_orientation_normals[Affs[cellid][iff]] = std::pair<bool, gmds::math::Vector> (AffIsOutward[cellid][iff], AffNormals[cellid][iff]);
                }
            }
        }


//        if(ANodeID == 13) {
//            std::cout<<"nbFakeFaces "<<nbFakeFaces.size()<<std::endl;
//        }

        // we keep the faces that appeared only once (non-internal faces)
        std::vector<gmds::math::Vector> normals_kept;
        for (auto ff: nbFakeFaces) {
            if (ff.second < 2) {

                //WARNING : remember that the normals are computed at ANodeID of the fakeface
                if(pair_orientation_normals[ff.first].first) {
                    normals_kept.push_back(pair_orientation_normals[ff.first].second);
                } else {
                    normals_kept.push_back((-1) * pair_orientation_normals[ff.first].second);
                }
            }
        }

//        if(ANodeID == 13) {
//            std::cout << "normals_kept " << normals_kept.size()<<std::endl;
//
//            for(auto n: normals_kept) {
//                std::cout<<n<<std::endl;
//            }
//        }

        // now check whether we have opposite normals
        for (auto n0: normals_kept) {
            for (auto n1: normals_kept) {
                kmds::TCoord dot = n0.dot(n1);

                if(dot < BADPILLOWINGDETECTION_OPPOSITE_DOTPRODUCT_THRESHOLD) {
                    return true;
                }
            }
        }

        return false;
    }

    /*----------------------------------------------------------------------------*/
    bool BadPillowDetection_IsNodeBad_glpk_3D(kmds::Mesh* AMesh,
                                              const kmds::TCellID ANodeID,
                                              const std::vector<kmds::TCellID> ACellIDs,
                                              const std::vector<std::vector<kmds::FakeFace> >& Affs,
                                              const std::vector<std::vector<bool> >& AffIsOutward,
                                              const std::vector<std::vector<gmds::math::Vector> >& AffNormals)
    {
        // if there is only one cell, it is always ok
        // WARNING : We assume it is also the case with two cells (if it is not the case
        // it would mean that we have a flat cell between those two)
        if(ACellIDs.size() <= 2) {
            return false;
        }

        // we build this set of regions FakeFaces and we compute how many there are considering unicity
        std::map<kmds::FakeFace, int> nbFakeFaces;
        std::map<kmds::FakeFace, std::pair<bool, gmds::math::Vector> > pair_orientation_normals;

        for(int ic=0; ic<ACellIDs.size(); ic++) {

            kmds::TCellID cellid = ACellIDs[ic];

            for(int iff=0; iff<Affs[cellid].size(); iff++) {

                if(Affs[cellid][iff].hasNode(ANodeID)) {

                    nbFakeFaces[Affs[cellid][iff]]++;
                    pair_orientation_normals[Affs[cellid][iff]] = std::pair<bool, gmds::math::Vector> (AffIsOutward[cellid][iff], AffNormals[cellid][iff]);
                }
            }
        }


        // we keep the faces that appeared only once (non-internal faces)
        std::vector<gmds::math::Vector> normals_kept;
        for (auto ff: nbFakeFaces) {
            if (ff.second < 2) {

                //WARNING : remember that the normals are computed at ANodeID of the fakeface
                if(pair_orientation_normals[ff.first].first) {
                    normals_kept.push_back(pair_orientation_normals[ff.first].second);
                } else {
                    normals_kept.push_back((-1) * pair_orientation_normals[ff.first].second);
                }
            }
        }


        // if there is zero or one normals kept, the problem has a solution
        if(normals_kept.size() < 2) {
            return false;
        }

        // =======================================
        // GLPK
        // =======================================

        glp_prob *lp;
        int *ia, *ja;
        double *ar;

        lp = glp_create_prob();
        glp_set_obj_dir(lp, GLP_MAX);

        // compute problem size
        int nbRows = 0;
        int nbCols = 0;
        int nnz = 0;

        // the three coordinates (x, y, z)
        nbCols += 3;

        // the half-spaces
        nbRows += normals_kept.size();
        nnz += nbCols * nbRows;

        // enforce constraints
        glp_add_rows(lp, nbRows);
        glp_add_cols(lp, nbCols);

        ia = (int *) calloc(1 + nnz, sizeof(int));
        ja = (int *) calloc(1 + nnz, sizeof(int));
        ar = (double *) calloc(1 + nnz, sizeof(double));

        int irow = 1; // start at 1, not 0 !
        int icol = 1; // start at 1, not 0 !
        int innz = 1; // start at 1, not 0 !


        double coefX = 0;
        double coefY = 0;
        double coefZ = 0;
        for(auto n: normals_kept) {
            coefX -= n.X();
            coefY -= n.Y();
            coefZ -= n.Z();
        }

        glp_set_col_name(lp, icol, "coordX");
        glp_set_col_kind(lp, icol, GLP_CV);
        glp_set_col_bnds(lp, icol, GLP_DB, -1, 1);
        glp_set_obj_coef(lp, icol, coefX);
        icol++;

        glp_set_col_name(lp, icol, "coordY");
        glp_set_col_kind(lp, icol, GLP_CV);
        glp_set_col_bnds(lp, icol, GLP_DB, -1, 1);
        glp_set_obj_coef(lp, icol, coefY);
        icol++;

        glp_set_col_name(lp, icol, "coordZ");
        glp_set_col_kind(lp, icol, GLP_CV);
        glp_set_col_bnds(lp, icol, GLP_DB, -1, 1);
        glp_set_obj_coef(lp, icol, coefZ);
        icol++;


        // one half-space per normal
        for(int n=0; n<normals_kept.size(); n++) {
            std::string row_name = "halfspace_" + std::to_string(n);
            glp_set_row_name(lp, irow, row_name.c_str());
            // WARNING: lower bound is cosinus(angle) with angle the angle between the (-1)*normal and the chosen direction ?
            glp_set_row_bnds(lp, irow, GLP_LO, BADPILLOWINGDETECTION_GLPK_THRESHOLD, 0.);

            ia[innz] = irow;
            ja[innz] = 1;
            ar[innz] = - normals_kept[n].X();
            innz++;

            ia[innz] = irow;
            ja[innz] = 2;
            ar[innz] = - normals_kept[n].Y();
            innz++;

            ia[innz] = irow;
            ja[innz] = 3;
            ar[innz] = - normals_kept[n].Z();
            innz++;

            irow++;
        }

        // add one because that is the indices starting index
//        for(int j=1; j<=nnz; j++) {
//            ja[j]++;
//        }

        glp_load_matrix(lp, nnz, ia, ja, ar);

         glp_term_out( GLP_OFF );
//        glp_term_out(GLP_ON);

        glp_smcp glpParams;
        glp_init_smcp(&glpParams);
        glpParams.presolve = GLP_ON;
        glpParams.meth = GLP_DUALP;


//        glp_write_lp(lp, NULL, "cplex.txt");

//        glp_write_mps(lp, GLP_MPS_FILE, NULL, "mps.txt");

        int glpErr = 0;
        glpErr = glp_simplex( lp, &glpParams);
//        glpErr = glp_simplex( lp, NULL);

        double value = glp_get_obj_val(lp);

//        std::cout<<"value "<<value<<std::endl;

        switch (glpErr) {
            case 0:
//                std::cout << "GLP OK" << std::endl;
                return false;
                break;
            default:
//                std::cout << "pb solving in GLP." << std::endl;
                return true;
                break;
        }

//        glpErr = glp_get_status(lp);
//        switch (glpErr) {
//            case GLP_UNDEF:
//                std::cout << " LP solution is undefined" << std::endl;
//                throw kmds::KException("SubMapping::boundaryDiscretization LP solution is undefined");
//                break;
//            case GLP_OPT:
//                std::cout << " LP solution is optimal" << std::endl;
//                break;
//            case GLP_FEAS:
//                std::cout << " LP solution is feasible" << std::endl;
//                break;
//            case GLP_INFEAS:
//                std::cout << "LP problem is infeasible" << std::endl;
//                throw kmds::KException("SubMapping::boundaryDiscretization LP problem is infeasible");
//                break;
//            case GLP_NOFEAS:
//                std::cout << "LP problem has no feasible solution" << std::endl;
//                throw kmds::KException("SubMapping::boundaryDiscretization LP problem has no feasible solution");
//                break;
//            case GLP_UNBND:
//                std::cout << "LP problem has unbounded solution" << std::endl;
//                throw kmds::KException("SubMapping::boundaryDiscretization LP problem has unbounded solution");
//                break;
//            default:
//                throw kmds::KException("SubMapping::boundaryDiscretization glp_simplex unknown return code.");
//        }
//
//
//        glp_print_sol(lp, "pixelAssignment.txt");


        return false;
    }

    /*----------------------------------------------------------------------------*/
    bool BadPillowDetection_IsNodeBad_r3d_3D(kmds::Mesh* AMesh,
                                             const kmds::TCellID ANodeID,
                                             const std::vector<kmds::TCellID> ACellIDs,
                                             const std::vector<std::vector<kmds::FakeFace> >& Affs,
                                             const std::vector<std::vector<bool> >& AffIsOutward,
                                             const std::vector<std::vector<gmds::math::Vector> >& AffNormals,
                                             const std::vector<std::vector<gmds::math::Triangle> >& AffTriangles)
    {
        // if there is only one cell, it is always ok
        // WARNING : We assume it is also the case with two cells (if it is not the case
        // it would mean that we have a flat cell between those two)
        if(ACellIDs.size() <= 2) {
            return false;
        }

        // we build this set of regions FakeFaces and we compute how many there are considering unicity
        std::map<kmds::FakeFace, int> nbFakeFaces;
        std::map<kmds::FakeFace, std::pair<bool, gmds::math::Vector> > pair_orientation_normals;
        std::map<kmds::FakeFace, gmds::math::Triangle> triangles;

        for(int ic=0; ic<ACellIDs.size(); ic++) {

            kmds::TCellID cellid = ACellIDs[ic];

            for(int iff=0; iff<Affs[cellid].size(); iff++) {

                if(Affs[cellid][iff].hasNode(ANodeID)) {

                    nbFakeFaces[Affs[cellid][iff]]++;
                    pair_orientation_normals[Affs[cellid][iff]] = std::pair<bool, gmds::math::Vector> (AffIsOutward[cellid][iff], AffNormals[cellid][iff]);
                    triangles[Affs[cellid][iff]] = AffTriangles[cellid][iff];
                }
            }
        }


        // we keep the faces that appeared only once (non-internal faces)
        std::vector<gmds::math::Vector> normals_kept;
        std::vector<gmds::math::Triangle> triangles_kept;
        for (auto ff: nbFakeFaces) {
            if (ff.second < 2) {

                //WARNING : remember that the normals are computed at ANodeID of the fakeface
                if(pair_orientation_normals[ff.first].first) {
                    normals_kept.push_back(pair_orientation_normals[ff.first].second);
                } else {
                    normals_kept.push_back((-1) * pair_orientation_normals[ff.first].second);
                }
                triangles_kept.push_back(triangles[ff.first]);
            }
        }


        // if there is zero or one normals kept, the problem has a solution
        if(normals_kept.size() < 2) {
            return false;
        }

        // =======================================
        // r3d
        // =======================================

        // build the hexahedron that represents the first half-space
        gmds::math::Vector p(triangles_kept[0].getPoint(0), triangles_kept[0].getPoint(1));
        p.normalize();
        const gmds::math::Vector q = p.cross(normals_kept[0]);


        r3d_poly src_cell;
        src_cell.nverts = 8;

        gmds::math::Point pt0 = ((-1.) * p - q).getPoint();
        gmds::math::Point pt1 = (        p - q).getPoint();
        gmds::math::Point pt2 = (        p + q).getPoint();
        gmds::math::Point pt3 = ((-1.) * p + q).getPoint();

        gmds::math::Vector r = normals_kept[0].opp();

        gmds::math::Point pt4 = ((-1.) * p - q + r).getPoint();
        gmds::math::Point pt5 = (        p - q + r).getPoint();
        gmds::math::Point pt6 = (        p + q + r).getPoint();
        gmds::math::Point pt7 = ((-1.) * p + q + r).getPoint();

        src_cell.verts[0].pos.xyz[0] = pt0.X();
        src_cell.verts[0].pos.xyz[1] = pt0.Y();
        src_cell.verts[0].pos.xyz[2] = pt0.Z();
        src_cell.verts[1].pos.xyz[0] = pt1.X();
        src_cell.verts[1].pos.xyz[1] = pt1.Y();
        src_cell.verts[1].pos.xyz[2] = pt1.Z();
        src_cell.verts[2].pos.xyz[0] = pt2.X();
        src_cell.verts[2].pos.xyz[1] = pt2.Y();
        src_cell.verts[2].pos.xyz[2] = pt2.Z();
        src_cell.verts[3].pos.xyz[0] = pt3.X();
        src_cell.verts[3].pos.xyz[1] = pt3.Y();
        src_cell.verts[3].pos.xyz[2] = pt3.Z();
        src_cell.verts[4].pos.xyz[0] = pt4.X();
        src_cell.verts[4].pos.xyz[1] = pt4.Y();
        src_cell.verts[4].pos.xyz[2] = pt4.Z();
        src_cell.verts[5].pos.xyz[0] = pt5.X();
        src_cell.verts[5].pos.xyz[1] = pt5.Y();
        src_cell.verts[5].pos.xyz[2] = pt5.Z();
        src_cell.verts[6].pos.xyz[0] = pt6.X();
        src_cell.verts[6].pos.xyz[1] = pt6.Y();
        src_cell.verts[6].pos.xyz[2] = pt6.Z();
        src_cell.verts[7].pos.xyz[0] = pt7.X();
        src_cell.verts[7].pos.xyz[1] = pt7.Y();
        src_cell.verts[7].pos.xyz[2] = pt7.Z();

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

        const int nbHalfspaces = normals_kept.size() - 1;
        r3d_plane planes[nbHalfspaces];

        // one half-space per normal
        for(int n=0; n<nbHalfspaces; n++) {

            planes[n].n.x = - normals_kept[n+1].X();
            planes[n].n.y = - normals_kept[n+1].Y();
            planes[n].n.z = - normals_kept[n+1].Z();
            planes[n].d = 0.;
        }

        const r3d_int POLY_ORDER = 2;
        r3d_real om[R3D_NUM_MOMENTS(POLY_ORDER)];

//        if(ANodeID == 13) {
//
//            for(int i=0; i<normals_kept.size(); i++) {
//                std::cout<<normals_kept[i]<<std::endl;
//            }
//
//            r3d_print(&src_cell);
//
//            r3d_reduce(&src_cell, om, POLY_ORDER);
//            std::cout<<"vol "<<om[0]<<std::endl;
//        }

        r3d_clip(&src_cell, planes, nbHalfspaces);

        r3d_reduce(&src_cell, om, POLY_ORDER);
        const double vol = om[0];

//        if(ANodeID == 13) {
//            std::cout<<"vol "<<vol<<std::endl;
//        }

        // TODO: what should be the volume threshold ?
        if(vol > BADPILLOWINGDETECTION_R3D_THRESHOLD) {
            return false;
        } else {
            return true;
        }
    }

    /*----------------------------------------------------------------------------*/
    bool BadPillowDetection_IsNodeBad_ortho_3D(kmds::Mesh* AMesh,
                                               const kmds::TCellID ANodeID,
                                               const std::vector<kmds::TCellID> ACellIDs,
                                               const std::vector<std::vector<kmds::FakeFace> >& Affs,
                                               const std::vector<std::vector<bool> >& AffIsOutward,
                                               const std::vector<std::vector<gmds::math::Vector> >& AffNormals,
                                               const std::vector<std::vector<gmds::math::Triangle> >& AffTriangles,
                                               const std::vector<std::vector<bool> >& ffBoundary)
    {
        // if there is only one cell, it is always ok
        // WARNING : We assume it is also the case with two cells (if it is not the case
        // it would mean that we have a flat cell between those two)
        if(ACellIDs.size() <= 2) {
            return false;
        }

//        if(ANodeID == 13 || ANodeID == 22){
//
//            std::cout<<"ANodeID "<<ANodeID<<" ACellIDs.size() "<<ACellIDs.size()<<std::endl;
//        }

        // we build this set of regions FakeFaces and we compute how many there are considering unicity
        std::map<kmds::FakeFace, int> nbFakeFaces;
        std::map<kmds::FakeFace, std::pair<bool, gmds::math::Vector> > pair_orientation_normals;
        std::map<kmds::FakeFace, gmds::math::Triangle> triangles;

        for(int ic=0; ic<ACellIDs.size(); ic++) {

            kmds::TCellID cellid = ACellIDs[ic];

            for(int iff=0; iff<Affs[cellid].size(); iff++) {

                if(Affs[cellid][iff].hasNode(ANodeID)) {

                    if(!ffBoundary[cellid][iff]) {

                        nbFakeFaces[Affs[cellid][iff]]++;
                        pair_orientation_normals[Affs[cellid][iff]] = std::pair<bool, gmds::math::Vector>(
                                AffIsOutward[cellid][iff], AffNormals[cellid][iff]);
                        triangles[Affs[cellid][iff]] = AffTriangles[cellid][iff];
                    }
                }
            }
        }


        // we keep the faces that appeared only once (non-internal faces)
        std::vector<gmds::math::Vector> normals_kept;
        std::vector<gmds::math::Triangle> triangles_kept;
        for (auto ff: nbFakeFaces) {
            if (ff.second < 2) {

                //WARNING : remember that the normals are computed at ANodeID of the fakeface
                if(pair_orientation_normals[ff.first].first) {
                    normals_kept.push_back(pair_orientation_normals[ff.first].second);
                } else {
                    normals_kept.push_back((-1) * pair_orientation_normals[ff.first].second);
                }
                triangles_kept.push_back(triangles[ff.first]);
            }
        }

//        if(ANodeID == 13 || ANodeID == 22){
//
//            std::cout<<"ANodeID "<<ANodeID<<" triangles_kept.size() "<<triangles_kept.size()<<std::endl;
//            for(int t=0; t<triangles_kept.size(); t++) {
//                std::cout<<triangles_kept[t]<<std::endl;
//                std::cout<<normals_kept[t]<<std::endl;
//                std::cout<<"angle "<<triangles_kept[t].angle()<<std::endl;
//            }
//        }


        // if there is zero or one normals kept, the problem has a solution
        if(normals_kept.size() < 2) {
            return false;
        }

        // =======================================
        // vector normal to the patch
        // =======================================
        gmds::math::Vector v;
        double area = 0.;

        for(int i=0; i<normals_kept.size(); i++) {
//            double weight = triangles_kept[i].area();
            double weight = triangles_kept[i].angle();

            v = v + weight * normals_kept[i];
            area += weight;
        }

        // TODO : we should ajust v so as to consider the geometric association of the nid node
        v = (1./area) * v;

//        if(ANodeID == 22 || ANodeID == 13) {
//            std::cout<<"ANodeID "<<ANodeID<<" v "<<v<<std::endl;
//        }

        // now check whether the normal to the patch is a good candidate direction for the duplicate node
        for(auto n: normals_kept) {

//            if(ANodeID == 22 || ANodeID == 13) {
//                std::cout<<"ANodeID "<<ANodeID<<" n "<<n<<" "<<v.dot(n)<<std::endl;
//            }

            if(v.dot(n) < BADPILLOWINGDETECTION_ORTHO_THRESHOLD) {
                return true;
            }
        }

        return false;
    }

//    /*----------------------------------------------------------------------------*/
//    struct Pillow_createTwoCells_3D {
//        const kmds::GrowingView<kmds::TCellID> *selection;
//        kmds::Mesh *mesh;
//        const kmds::Connectivity* c_F2C;
//        elg3d::MaterialAssignment* ma;
//        const std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> >* node2newNodes;
//
//
//        Pillow_createTwoCells_3D(
//                const kmds::GrowingView<kmds::TCellID> *selection_,
//                kmds::Mesh *mesh_,
//                const kmds::Connectivity* c_F2C_,
//                elg3d::MaterialAssignment* ma_,
//                std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> >* node2newNodes_
//        )
//                : selection(selection_)
//                , mesh(mesh_)
//                , c_F2C(c_F2C_)
//                , ma(ma_)
//                , node2newNodes(node2newNodes_)
//        {
//        }
//
//        KOKKOS_INLINE_FUNCTION
//        void
//        operator()(int i) const {
//            int fid = selection->get(i);
//
//            Kokkos::View<kmds::TCellID*> cells;
//            c_F2C->get(fid, cells);
//
//            kmds::Face f = mesh->getFace(fid);
//            Kokkos::View<kmds::TCellID*> nodes;
//            f.nodeIds(nodes);
//
//
//            // create two cells
//            kmds::TCellID cid0_new;
//            if(nodes.size() == 4) {
//                cid0_new = mesh->addHexahedra(2);
//            } else {
//                cid0_new = mesh->addPrism3s(2);
//            }
//
//
//            // determine which new nodes will be used for which cell
//            int mat0 = ma->getMaterial(cells(0));
//            int mat1 = ma->getMaterial(cells(1));
//
//            kmds::Region c0 = mesh->getRegion(cells(0));
//
//            bool isOutward0 = c0.isFaceOrientedOutward(nodes);
//
//            kmds::TCellID newIDs_0[8];
//            kmds::TCellID newIDs_1[8];
//
//            if(isOutward0) {
//                for(int i_n=0; i_n<nodes.size(); i_n++) {
//                    newIDs_0[i_n] = (*node2newNodes)[mat0].value_at((*node2newNodes)[mat0].find(nodes(i_n)));
//                    newIDs_0[i_n + nodes.size()] = nodes(i_n);
//                    newIDs_1[i_n] = nodes(i_n);
//                    newIDs_1[i_n + nodes.size()] = (*node2newNodes)[mat1].value_at((*node2newNodes)[mat1].find(nodes(i_n)));
//                }
//            } else {
//                for(int i_n=0; i_n<nodes.size(); i_n++) {
//                    newIDs_0[i_n] = nodes(i_n);
//                    newIDs_0[i_n + nodes.size()] = (*node2newNodes)[mat0].value_at((*node2newNodes)[mat0].find(nodes(i_n)));
//                    newIDs_1[i_n] = (*node2newNodes)[mat1].value_at((*node2newNodes)[mat1].find(nodes(i_n)));
//                    newIDs_1[i_n + nodes.size()] = nodes(i_n);
//                }
//            }
//
//            kmds::Region c0_new = mesh->getRegion(cid0_new);
//            kmds::Region c1_new = mesh->getRegion(cid0_new + 1);
//            c0_new.setNodes(newIDs_0, nodes.size() * 2);
//            c1_new.setNodes(newIDs_1, nodes.size() * 2);
//
//            ma->setMaterial(mat0, cid0_new);
//            ma->setMaterial(mat1, cid0_new + 1);
//        }
//
//    };
//
    /*----------------------------------------------------------------------------*/
    void
    BadPillowDetection_detect_3D(const kmds::GrowingView<kmds::TCellID>* AInterfaceFaces,
                                 const kmds::GrowingView<kmds::TCellID>* AInterfaceNodes,
                                 kmds::Mesh* AMesh,
                                 const kmds::Connectivity* c_N2F,
                                 const kmds::Connectivity* c_F2C,
                                 const elg3d::MaterialAssignment* Ama,
                                 kmds::GrowingView<kmds::TCellID>* ASelectionBadNodes)
    {

        // necessary init
//        kmds::Variable<bool>* varIsInterfaceNode = AMesh->createVariable<bool> (false, kmds::KMDS_NODE, "pillow_varIsInterfaceNode");
//        Kokkos::parallel_for(AInterfaceNodes->getNbElems(),
//                             KOKKOS_LAMBDA(const int i) {
//                                 (*varIsInterfaceNode)[AInterfaceNodes->get(i)] = true;
//                             });


//        std::vector<Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> > node2newNodes;
//        for(int imat=0; imat<Ama->getNbMaterials(); imat++) {
//            node2newNodes.push_back(Kokkos::UnorderedMap<kmds::TCellID, kmds::TCellID> (AInterfaceNodes->getNbElems()));
//        }


        // for each interface node check whether we have a bad configuration for each material
        // WARNING :  each material is manifold
//        Kokkos::parallel_for(AInterfaceNodes->getNbElems(), BadPillowDetection_isNodeBad(AInterfaceNodes,
//                                                                                         AMesh,
//                                                                                         c_N2F,
//                                                                                         c_F2C,
//                                                                                         Ama,
//                                                                                         ASelectionBadNodes));


//        // we need to get this cell container before creating the new cells
//        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", AMesh->getNbRegions());
//        AMesh->getRegionIDs_dummy(&cellIDs);
//
//
//        //
//        Kokkos::parallel_for(AInterfaceFaces->getNbElems(), Pillow_createTwoCells_3D(AInterfaceFaces,
//                                                                                     AMesh,
//                                                                                     c_F2C,
//                                                                                     Ama,
//                                                                                     &node2newNodes));
//
//
//        //
//// TODO get regions container, limit it to regions adjacent to interface nodes
//        Kokkos::parallel_for(cellIDs.getNbElems(), Pillow_replaceInCells_3D(&cellIDs,
//                                                                            AMesh,
//                                                                            Ama,
//                                                                            varIsInterfaceNode,
//                                                                            &node2newNodes));
//
//
//        // clean-up temporary data
//        AMesh->deleteVariable(kmds::KMDS_NODE, "pillow_varIsInterfaceNode");

    }
/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
