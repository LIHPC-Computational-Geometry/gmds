

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
/** \file    BoundingBoxGeomAssociation.cpp
 *  \author  legoff
 *  \date    04/23/2018
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/BoundingBoxGeomAssociation.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
namespace elg3d {

/*----------------------------------------------------------------------------*/
    void
    BoundingBoxGeomAssociation_init_2D(kmds::Mesh* AMesh,
                                       gmds::cad::FACManager* AFacGeomManager,
                                       kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation) {

        // first get the bounding box of the mesh
        kmds::TSize nbNodes = AMesh->getNbNodes();
        kmds::GrowingView<kmds::TCellID> nodeIDs("CELLS", nbNodes);
        AMesh->getNodeIDs_dummy(&nodeIDs);

        kmds::TCoord minXYZ[3];
        kmds::TCoord maxXYZ[3];

        minXYZ[0] = HUGE_VALF;
        minXYZ[1] = HUGE_VALF;
        minXYZ[2] = HUGE_VALF;
        maxXYZ[0] = -HUGE_VALF;
        maxXYZ[1] = -HUGE_VALF;
        maxXYZ[2] = -HUGE_VALF;

        for (int i = 0; i < nbNodes; i++) {

            kmds::TCellID nid = nodeIDs.get(i);
            kmds::Node n = AMesh->getNode(nid);
            kmds::TCoord xyz[3];
            n.getLocation(xyz[0], xyz[1], xyz[2]);

            if (xyz[0] < minXYZ[0]) minXYZ[0] = xyz[0];
            if (xyz[1] < minXYZ[1]) minXYZ[1] = xyz[1];
            if (xyz[2] < minXYZ[2]) minXYZ[2] = xyz[2];
            if (xyz[0] > maxXYZ[0]) maxXYZ[0] = xyz[0];
            if (xyz[1] > maxXYZ[1]) maxXYZ[1] = xyz[1];
            if (xyz[2] > maxXYZ[2]) maxXYZ[2] = xyz[2];
        }

        // then build the geom model
        gmds::Mesh &mesh = AFacGeomManager->getMeshView();
        gmds::Node corners[4];
        corners[0] = mesh.newNode(minXYZ[0], minXYZ[1], minXYZ[2]);
        corners[1] = mesh.newNode(maxXYZ[0], minXYZ[1], minXYZ[2]);
        corners[2] = mesh.newNode(maxXYZ[0], maxXYZ[1], minXYZ[2]);
        corners[3] = mesh.newNode(minXYZ[0], maxXYZ[1], minXYZ[2]);


        gmds::Face triangles[2];
        triangles[0] = mesh.newTriangle(corners[0], corners[1], corners[2]);
        triangles[1] = mesh.newTriangle(corners[0], corners[2], corners[3]);

        gmds::CellGroup<gmds::Node>* corner0 = mesh.newGroup<gmds::Node>("corner_0");
        corner0->add(corners[0]);
        gmds::CellGroup<gmds::Node>* corner1 = mesh.newGroup<gmds::Node>("corner_1");
        corner1->add(corners[1]);
        gmds::CellGroup<gmds::Node>* corner2 = mesh.newGroup<gmds::Node>("corner_2");
        corner2->add(corners[2]);
        gmds::CellGroup<gmds::Node>* corner3 = mesh.newGroup<gmds::Node>("corner_3");
        corner3->add(corners[3]);

        gmds::CellGroup<gmds::Node>* curve_0 = mesh.newGroup<gmds::Node>("curve_bottom");
        curve_0->add(corners[0]);
        curve_0->add(corners[1]);
        gmds::CellGroup<gmds::Node>* curve_1 = mesh.newGroup<gmds::Node>("curve_right");
        curve_1->add(corners[1]);
        curve_1->add(corners[2]);
        gmds::CellGroup<gmds::Node>* curve_2 = mesh.newGroup<gmds::Node>("curve_top");
        curve_2->add(corners[2]);
        curve_2->add(corners[3]);
        gmds::CellGroup<gmds::Node>* curve_3 = mesh.newGroup<gmds::Node>("curve_left");
        curve_3->add(corners[3]);
        curve_3->add(corners[0]);

        gmds::CellGroup<gmds::Face>* surf = mesh.newGroup<gmds::Face>("surf");
        surf->add(triangles[0]);
        surf->add(triangles[1]);

        AFacGeomManager->updateFromMesh();

        std::vector<gmds::cad::GeomPoint *> vertices;
        std::vector<gmds::cad::GeomCurve *> curves;
        std::vector<gmds::cad::GeomSurface *> surfaces;
        AFacGeomManager->getPoints(vertices);
        AFacGeomManager->getCurves(curves);
        AFacGeomManager->getSurfaces(surfaces);
        std::map<std::string, std::uintptr_t> name2Points;
        std::map<std::string, std::uintptr_t> name2Curves;
        std::map<std::string, std::uintptr_t> name2Surfaces;
        for (auto v: vertices) {
            name2Points[v->getName()] = reinterpret_cast<std::uintptr_t> (v);
        }
        for (auto c: curves) {
            name2Curves[c->getName()] = reinterpret_cast<std::uintptr_t> (c);
        }
        for (auto s: surfaces) {
            name2Surfaces[s->getName()] = reinterpret_cast<std::uintptr_t> (s);
        }

        double tol = gmds::math::Point(minXYZ[0], minXYZ[1], minXYZ[2]).distance(
                gmds::math::Point(maxXYZ[0], maxXYZ[1], maxXYZ[2])) * 10e-14;

        typedef enum {
            LEFT = 1, RIGHT = 1 << 1,
            BOTTOM = 1 << 2, TOP = 1 << 3
        } sides;

        for (int i = 0; i < nbNodes; i++) {

            kmds::TCellID nid = nodeIDs.get(i);
            kmds::Node n = AMesh->getNode(nid);
            kmds::TCoord xyz[3];
            n.getLocation(xyz[0], xyz[1], xyz[2]);

            gmds::math::Point pt(xyz[0], xyz[1], xyz[2]);
            int pos = 0;

            if (pt.X() < minXYZ[0] + tol) {
                pos = pos | LEFT;
            }
            if (pt.X() > maxXYZ[0] - tol) {
                pos = pos | RIGHT;
            }
            if (pt.Y() < minXYZ[1] + tol) {
                pos = pos | BOTTOM;
            }
            if (pt.Y() > maxXYZ[1] - tol) {
                pos = pos | TOP;
            }

            // associate to curves
            if ((pos | LEFT) == pos) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_left");
            }
            if ((pos | RIGHT) == pos) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_right");
            }
            if ((pos | BOTTOM) == pos) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_bottom");
            }
            if ((pos | TOP) == pos) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_top");
            }

            // associate to vertices
            if (((pos | BOTTOM) == pos) && ((pos | RIGHT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Points.at("corner_1");
            }
            if (((pos | BOTTOM) == pos) && ((pos | LEFT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Points.at("corner_0");
            }
            if (((pos | TOP) == pos) && ((pos | RIGHT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Points.at("corner_2");
            }
            if (((pos | TOP) == pos) && ((pos | LEFT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Points.at("corner_3");
            }

        }
    }
/*----------------------------------------------------------------------------*/
void
BoundingBoxGeomAssociation_init_3D(kmds::Mesh* AMesh,
                                   gmds::cad::FACManager* AFacGeomManager,
                                   kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation)
    {

        // first get the bounding box of the mesh
        kmds::TSize nbNodes = AMesh->getNbNodes();
        kmds::GrowingView<kmds::TCellID> nodeIDs("CELLS", nbNodes);
        AMesh->getNodeIDs_dummy(&nodeIDs);

        kmds::TCoord minXYZ[3];
        kmds::TCoord maxXYZ[3];

        minXYZ[0] = HUGE_VALF;
        minXYZ[1] = HUGE_VALF;
        minXYZ[2] = HUGE_VALF;
        maxXYZ[0] = -HUGE_VALF;
        maxXYZ[1] = -HUGE_VALF;
        maxXYZ[2] = -HUGE_VALF;

        for(int i=0; i<nbNodes; i++) {

            kmds::TCellID nid = nodeIDs.get(i);
            kmds::Node n = AMesh->getNode(nid);
            kmds::TCoord xyz[3];
            n.getLocation(xyz[0], xyz[1], xyz[2]);

            if(xyz[0] < minXYZ[0]) minXYZ[0] = xyz[0];
            if(xyz[1] < minXYZ[1]) minXYZ[1] = xyz[1];
            if(xyz[2] < minXYZ[2]) minXYZ[2] = xyz[2];
            if(xyz[0] > maxXYZ[0]) maxXYZ[0] = xyz[0];
            if(xyz[1] > maxXYZ[1]) maxXYZ[1] = xyz[1];
            if(xyz[2] > maxXYZ[2]) maxXYZ[2] = xyz[2];
        }

        // then build the geom model
        gmds::Mesh& mesh = AFacGeomManager->getMeshView();
        gmds::Node corners[8];
        corners[0] = mesh.newNode(minXYZ[0], minXYZ[1], minXYZ[2]);
        corners[1] = mesh.newNode(maxXYZ[0], minXYZ[1], minXYZ[2]);
        corners[2] = mesh.newNode(maxXYZ[0], maxXYZ[1], minXYZ[2]);
        corners[3] = mesh.newNode(minXYZ[0], maxXYZ[1], minXYZ[2]);
        corners[4] = mesh.newNode(minXYZ[0], minXYZ[1], maxXYZ[2]);
        corners[5] = mesh.newNode(maxXYZ[0], minXYZ[1], maxXYZ[2]);
        corners[6] = mesh.newNode(maxXYZ[0], maxXYZ[1], maxXYZ[2]);
        corners[7] = mesh.newNode(minXYZ[0], maxXYZ[1], maxXYZ[2]);

        gmds::Face triangles[12];
        triangles[ 0] = mesh.newTriangle(corners[0], corners[1], corners[5]); // front
        triangles[ 1] = mesh.newTriangle(corners[0], corners[5], corners[4]);
        triangles[ 2] = mesh.newTriangle(corners[2], corners[3], corners[7]); // back
        triangles[ 3] = mesh.newTriangle(corners[2], corners[7], corners[6]);
        triangles[ 4] = mesh.newTriangle(corners[3], corners[0], corners[4]); // left
        triangles[ 5] = mesh.newTriangle(corners[3], corners[4], corners[7]);
        triangles[ 6] = mesh.newTriangle(corners[1], corners[2], corners[6]); // right
        triangles[ 7] = mesh.newTriangle(corners[1], corners[6], corners[5]);
        triangles[ 8] = mesh.newTriangle(corners[0], corners[3], corners[2]); // bottom
        triangles[ 9] = mesh.newTriangle(corners[0], corners[2], corners[1]);
        triangles[10] = mesh.newTriangle(corners[4], corners[5], corners[6]); // top
        triangles[11] = mesh.newTriangle(corners[4], corners[6], corners[7]);

        gmds::CellGroup<gmds::Node>* corner0 = mesh.newGroup<gmds::Node>("corner_0");
        corner0->add(corners[0]);
        gmds::CellGroup<gmds::Node>* corner1 = mesh.newGroup<gmds::Node>("corner_1");
        corner1->add(corners[1]);
        gmds::CellGroup<gmds::Node>* corner2 = mesh.newGroup<gmds::Node>("corner_2");
        corner2->add(corners[2]);
        gmds::CellGroup<gmds::Node>* corner3 = mesh.newGroup<gmds::Node>("corner_3");
        corner3->add(corners[3]);
        gmds::CellGroup<gmds::Node>* corner4 = mesh.newGroup<gmds::Node>("corner_4");
        corner4->add(corners[4]);
        gmds::CellGroup<gmds::Node>* corner5 = mesh.newGroup<gmds::Node>("corner_5");
        corner5->add(corners[5]);
        gmds::CellGroup<gmds::Node>* corner6 = mesh.newGroup<gmds::Node>("corner_6");
        corner6->add(corners[6]);
        gmds::CellGroup<gmds::Node>* corner7 = mesh.newGroup<gmds::Node>("corner_7");
        corner7->add(corners[7]);

        gmds::CellGroup<gmds::Node>* curve_0 = mesh.newGroup<gmds::Node>("curve_0");
        curve_0->add(corners[0]);
        curve_0->add(corners[1]);
        gmds::CellGroup<gmds::Node>* curve_1 = mesh.newGroup<gmds::Node>("curve_1");
        curve_1->add(corners[1]);
        curve_1->add(corners[2]);
        gmds::CellGroup<gmds::Node>* curve_2 = mesh.newGroup<gmds::Node>("curve_2");
        curve_2->add(corners[2]);
        curve_2->add(corners[3]);
        gmds::CellGroup<gmds::Node>* curve_3 = mesh.newGroup<gmds::Node>("curve_3");
        curve_3->add(corners[3]);
        curve_3->add(corners[0]);
        gmds::CellGroup<gmds::Node>* curve_4 = mesh.newGroup<gmds::Node>("curve_4");
        curve_4->add(corners[4]);
        curve_4->add(corners[5]);
        gmds::CellGroup<gmds::Node>* curve_5 = mesh.newGroup<gmds::Node>("curve_5");
        curve_5->add(corners[5]);
        curve_5->add(corners[6]);
        gmds::CellGroup<gmds::Node>* curve_6 = mesh.newGroup<gmds::Node>("curve_6");
        curve_6->add(corners[6]);
        curve_6->add(corners[7]);
        gmds::CellGroup<gmds::Node>* curve_7 = mesh.newGroup<gmds::Node>("curve_7");
        curve_7->add(corners[7]);
        curve_7->add(corners[4]);
        gmds::CellGroup<gmds::Node>* curve_8 = mesh.newGroup<gmds::Node>("curve_8");
        curve_8->add(corners[0]);
        curve_8->add(corners[4]);
        gmds::CellGroup<gmds::Node>* curve_9 = mesh.newGroup<gmds::Node>("curve_9");
        curve_9->add(corners[1]);
        curve_9->add(corners[5]);
        gmds::CellGroup<gmds::Node>* curve_10 = mesh.newGroup<gmds::Node>("curve_10");
        curve_10->add(corners[2]);
        curve_10->add(corners[6]);
        gmds::CellGroup<gmds::Node>* curve_11 = mesh.newGroup<gmds::Node>("curve_11");
        curve_11->add(corners[3]);
        curve_11->add(corners[7]);

        gmds::CellGroup<gmds::Face>* surf_front = mesh.newGroup<gmds::Face>("surf_front");
        surf_front->add(triangles[0]);
        surf_front->add(triangles[1]);
        gmds::CellGroup<gmds::Face>* surf_back = mesh.newGroup<gmds::Face>("surf_back");
        surf_back->add(triangles[2]);
        surf_back->add(triangles[3]);
        gmds::CellGroup<gmds::Face>* surf_left = mesh.newGroup<gmds::Face>("surf_left");
        surf_left->add(triangles[4]);
        surf_left->add(triangles[5]);
        gmds::CellGroup<gmds::Face>* surf_right = mesh.newGroup<gmds::Face>("surf_right");
        surf_right->add(triangles[6]);
        surf_right->add(triangles[7]);
        gmds::CellGroup<gmds::Face>* surf_bottom = mesh.newGroup<gmds::Face>("surf_bottom");
        surf_bottom->add(triangles[8]);
        surf_bottom->add(triangles[9]);
        gmds::CellGroup<gmds::Face>* surf_top = mesh.newGroup<gmds::Face>("surf_top");
        surf_top->add(triangles[10]);
        surf_top->add(triangles[11]);

        AFacGeomManager->updateFromMesh();

        std::vector<gmds::cad::GeomPoint*> vertices;
        std::vector<gmds::cad::GeomCurve*> curves;
        std::vector<gmds::cad::GeomSurface*> surfaces;
        AFacGeomManager->getPoints(vertices);
        AFacGeomManager->getCurves(curves);
        AFacGeomManager->getSurfaces(surfaces);
        std::map<std::string, std::uintptr_t> name2Points;
        std::map<std::string, std::uintptr_t> name2Curves;
        std::map<std::string, std::uintptr_t> name2Surfaces;
        for(auto v: vertices) {
            name2Points[v->getName()] = reinterpret_cast<std::uintptr_t> (v);
        }
        for(auto c: curves) {
            name2Curves[c->getName()] = reinterpret_cast<std::uintptr_t> (c);
        }
        for(auto s: surfaces) {
            name2Surfaces[s->getName()] = reinterpret_cast<std::uintptr_t> (s);
        }

        double tol = gmds::math::Point(minXYZ[0],minXYZ[1],minXYZ[2]).distance(gmds::math::Point(maxXYZ[0],maxXYZ[1],maxXYZ[2])) * 10e-14;

        typedef enum {
            LEFT = 1		  , RIGHT = 1 << 1  , FRONT = 1 << 2  , BACK = 1 << 3  ,
            BOTTOM = 1 << 4  , TOP = 1 << 5
        } sides;

        for(int i=0; i<nbNodes; i++) {

            kmds::TCellID nid = nodeIDs.get(i);
            kmds::Node n = AMesh->getNode(nid);
            kmds::TCoord xyz[3];
            n.getLocation(xyz[0], xyz[1], xyz[2]);

            gmds::math::Point pt(xyz[0], xyz[1], xyz[2]);
            int pos = 0;

            if(pt.X() < minXYZ[0]+tol ) {
                pos = pos|LEFT;
            }
            if(pt.X() > maxXYZ[0]-tol ) {
                pos = pos|RIGHT;
            }
            if(pt.Y() < minXYZ[1]+tol ) {
                pos = pos|FRONT;
            }
            if(pt.Y() > maxXYZ[1]-tol ) {
                pos = pos|BACK;
            }
            if(pt.Z() < minXYZ[2]+tol ) {
                pos = pos|BOTTOM;
            }
            if(pt.Z() > maxXYZ[2]-tol ) {
                pos = pos|TOP;
            }

            // associate to surfaces
            if((pos|LEFT) == pos) {
                (*AVarNodeGeomAssociation)[nid] = name2Surfaces.at("surf_left");
            }
            if((pos|RIGHT) == pos) {
                (*AVarNodeGeomAssociation)[nid] = name2Surfaces.at("surf_right");
            }
            if((pos|FRONT) == pos) {
                (*AVarNodeGeomAssociation)[nid] = name2Surfaces.at("surf_front");
            }
            if((pos|BACK) == pos) {
                (*AVarNodeGeomAssociation)[nid] = name2Surfaces.at("surf_back");
            }
            if((pos|BOTTOM) == pos) {
                (*AVarNodeGeomAssociation)[nid] = name2Surfaces.at("surf_bottom");
            }
            if((pos|TOP) == pos) {
                (*AVarNodeGeomAssociation)[nid] = name2Surfaces.at("surf_top");
            }

            // associate to curves
            if(((pos|BOTTOM) == pos) && ((pos|FRONT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_0");
            }
            if(((pos|BOTTOM) == pos) && ((pos|RIGHT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_1");
            }
            if(((pos|BOTTOM) == pos) && ((pos|BACK) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_2");
            }
            if(((pos|BOTTOM) == pos) && ((pos|LEFT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_3");
            }
            if(((pos|TOP) == pos) && ((pos|FRONT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_4");
            }
            if(((pos|TOP) == pos) && ((pos|RIGHT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_5");
            }
            if(((pos|TOP) == pos) && ((pos|BACK) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_6");
            }
            if(((pos|TOP) == pos) && ((pos|LEFT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_7");
            }
            if(((pos|LEFT) == pos) && ((pos|FRONT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_8");
            }
            if(((pos|FRONT) == pos) && ((pos|RIGHT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_9");
            }
            if(((pos|RIGHT) == pos) && ((pos|BACK) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_10");
            }
            if(((pos|BACK) == pos) && ((pos|LEFT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_11");
            }

            // associate to vertices
            if(((pos|BOTTOM) == pos) && ((pos|FRONT) == pos) && ((pos|LEFT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Points.at("corner_0");
            }
            if(((pos|BOTTOM) == pos) && ((pos|FRONT) == pos) && ((pos|RIGHT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Points.at("corner_1");
            }
            if(((pos|BOTTOM) == pos) && ((pos|BACK) == pos) && ((pos|RIGHT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Points.at("corner_2");
            }
            if(((pos|BOTTOM) == pos) && ((pos|BACK) == pos) && ((pos|LEFT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Points.at("corner_3");
            }
            if(((pos|TOP) == pos) && ((pos|FRONT) == pos) && ((pos|LEFT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Points.at("corner_4");
            }
            if(((pos|TOP) == pos) && ((pos|FRONT) == pos) && ((pos|RIGHT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Points.at("corner_5");
            }
            if(((pos|TOP) == pos) && ((pos|BACK) == pos) && ((pos|RIGHT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Points.at("corner_6");
            }
            if(((pos|TOP) == pos) && ((pos|BACK) == pos) && ((pos|LEFT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Points.at("corner_7");
            }
        }
    }

    /*----------------------------------------------------------------------------*/
    void
    BoundingBoxGeomAssociation_assoc_3D(kmds::Mesh* AMesh,
                                        const double minXYZ[3],
                                        const double maxXYZ[3],
                                        const gmds::cad::FACManager* AFacGeomManager,
                                        kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation)
    {

        // first get the bounding box of the mesh
        kmds::TSize nbNodes = AMesh->getNbNodes();
        kmds::GrowingView<kmds::TCellID> nodeIDs("CELLS", nbNodes);
        AMesh->getNodeIDs(&nodeIDs);

//        kmds::TCoord minXYZ[3];
//        kmds::TCoord maxXYZ[3];
//
//        minXYZ[0] = HUGE_VALF;
//        minXYZ[1] = HUGE_VALF;
//        minXYZ[2] = HUGE_VALF;
//        maxXYZ[0] = -HUGE_VALF;
//        maxXYZ[1] = -HUGE_VALF;
//        maxXYZ[2] = -HUGE_VALF;
//
//        for(int i=0; i<nbNodes; i++) {
//
//            kmds::TCellID nid = nodeIDs.get(i);
//            kmds::Node n = AMesh->getNode(nid);
//            kmds::TCoord xyz[3];
//            n.getLocation(xyz[0], xyz[1], xyz[2]);
//
//            if(xyz[0] < minXYZ[0]) minXYZ[0] = xyz[0];
//            if(xyz[1] < minXYZ[1]) minXYZ[1] = xyz[1];
//            if(xyz[2] < minXYZ[2]) minXYZ[2] = xyz[2];
//            if(xyz[0] > maxXYZ[0]) maxXYZ[0] = xyz[0];
//            if(xyz[1] > maxXYZ[1]) maxXYZ[1] = xyz[1];
//            if(xyz[2] > maxXYZ[2]) maxXYZ[2] = xyz[2];
//        }
//
//        // then build the geom model
//        gmds::IGMesh& mesh = AFacGeomManager->getMeshView();
//        gmds::Node corners[8];
//        corners[0] = mesh.newNode(minXYZ[0], minXYZ[1], minXYZ[2]);
//        corners[1] = mesh.newNode(maxXYZ[0], minXYZ[1], minXYZ[2]);
//        corners[2] = mesh.newNode(maxXYZ[0], maxXYZ[1], minXYZ[2]);
//        corners[3] = mesh.newNode(minXYZ[0], maxXYZ[1], minXYZ[2]);
//        corners[4] = mesh.newNode(minXYZ[0], minXYZ[1], maxXYZ[2]);
//        corners[5] = mesh.newNode(maxXYZ[0], minXYZ[1], maxXYZ[2]);
//        corners[6] = mesh.newNode(maxXYZ[0], maxXYZ[1], maxXYZ[2]);
//        corners[7] = mesh.newNode(minXYZ[0], maxXYZ[1], maxXYZ[2]);
//
//        gmds::Face triangles[12];
//        triangles[ 0] = mesh.newTriangle(corners[0], corners[1], corners[5]); // front
//        triangles[ 1] = mesh.newTriangle(corners[0], corners[5], corners[4]);
//        triangles[ 2] = mesh.newTriangle(corners[2], corners[3], corners[7]); // back
//        triangles[ 3] = mesh.newTriangle(corners[2], corners[7], corners[6]);
//        triangles[ 4] = mesh.newTriangle(corners[3], corners[0], corners[4]); // left
//        triangles[ 5] = mesh.newTriangle(corners[3], corners[4], corners[7]);
//        triangles[ 6] = mesh.newTriangle(corners[1], corners[2], corners[6]); // right
//        triangles[ 7] = mesh.newTriangle(corners[1], corners[6], corners[5]);
//        triangles[ 8] = mesh.newTriangle(corners[0], corners[3], corners[2]); // bottom
//        triangles[ 9] = mesh.newTriangle(corners[0], corners[2], corners[1]);
//        triangles[10] = mesh.newTriangle(corners[4], corners[5], corners[6]); // top
//        triangles[11] = mesh.newTriangle(corners[4], corners[6], corners[7]);
//
//        gmds::CellGroup<gmds::Node>* corner0 = mesh.newGroup<gmds::Node>("corner_0");
//        corner0->add(corners[0]);
//        gmds::CellGroup<gmds::Node>* corner1 = mesh.newGroup<gmds::Node>("corner_1");
//        corner1->add(corners[1]);
//        gmds::CellGroup<gmds::Node>* corner2 = mesh.newGroup<gmds::Node>("corner_2");
//        corner2->add(corners[2]);
//        gmds::CellGroup<gmds::Node>* corner3 = mesh.newGroup<gmds::Node>("corner_3");
//        corner3->add(corners[3]);
//        gmds::CellGroup<gmds::Node>* corner4 = mesh.newGroup<gmds::Node>("corner_4");
//        corner4->add(corners[4]);
//        gmds::CellGroup<gmds::Node>* corner5 = mesh.newGroup<gmds::Node>("corner_5");
//        corner5->add(corners[5]);
//        gmds::CellGroup<gmds::Node>* corner6 = mesh.newGroup<gmds::Node>("corner_6");
//        corner6->add(corners[6]);
//        gmds::CellGroup<gmds::Node>* corner7 = mesh.newGroup<gmds::Node>("corner_7");
//        corner7->add(corners[7]);
//
//        gmds::CellGroup<gmds::Node>* curve_0 = mesh.newGroup<gmds::Node>("curve_0");
//        curve_0->add(corners[0]);
//        curve_0->add(corners[1]);
//        gmds::CellGroup<gmds::Node>* curve_1 = mesh.newGroup<gmds::Node>("curve_1");
//        curve_1->add(corners[1]);
//        curve_1->add(corners[2]);
//        gmds::CellGroup<gmds::Node>* curve_2 = mesh.newGroup<gmds::Node>("curve_2");
//        curve_2->add(corners[2]);
//        curve_2->add(corners[3]);
//        gmds::CellGroup<gmds::Node>* curve_3 = mesh.newGroup<gmds::Node>("curve_3");
//        curve_3->add(corners[3]);
//        curve_3->add(corners[0]);
//        gmds::CellGroup<gmds::Node>* curve_4 = mesh.newGroup<gmds::Node>("curve_4");
//        curve_4->add(corners[4]);
//        curve_4->add(corners[5]);
//        gmds::CellGroup<gmds::Node>* curve_5 = mesh.newGroup<gmds::Node>("curve_5");
//        curve_5->add(corners[5]);
//        curve_5->add(corners[6]);
//        gmds::CellGroup<gmds::Node>* curve_6 = mesh.newGroup<gmds::Node>("curve_6");
//        curve_6->add(corners[6]);
//        curve_6->add(corners[7]);
//        gmds::CellGroup<gmds::Node>* curve_7 = mesh.newGroup<gmds::Node>("curve_7");
//        curve_7->add(corners[7]);
//        curve_7->add(corners[4]);
//        gmds::CellGroup<gmds::Node>* curve_8 = mesh.newGroup<gmds::Node>("curve_8");
//        curve_8->add(corners[0]);
//        curve_8->add(corners[4]);
//        gmds::CellGroup<gmds::Node>* curve_9 = mesh.newGroup<gmds::Node>("curve_9");
//        curve_9->add(corners[1]);
//        curve_9->add(corners[5]);
//        gmds::CellGroup<gmds::Node>* curve_10 = mesh.newGroup<gmds::Node>("curve_10");
//        curve_10->add(corners[2]);
//        curve_10->add(corners[6]);
//        gmds::CellGroup<gmds::Node>* curve_11 = mesh.newGroup<gmds::Node>("curve_11");
//        curve_11->add(corners[3]);
//        curve_11->add(corners[7]);
//
//        gmds::CellGroup<gmds::Face>* surf_front = mesh.newGroup<gmds::Face>("surf_front");
//        surf_front->add(triangles[0]);
//        surf_front->add(triangles[1]);
//        gmds::CellGroup<gmds::Face>* surf_back = mesh.newGroup<gmds::Face>("surf_back");
//        surf_back->add(triangles[2]);
//        surf_back->add(triangles[3]);
//        gmds::CellGroup<gmds::Face>* surf_left = mesh.newGroup<gmds::Face>("surf_left");
//        surf_left->add(triangles[4]);
//        surf_left->add(triangles[5]);
//        gmds::CellGroup<gmds::Face>* surf_right = mesh.newGroup<gmds::Face>("surf_right");
//        surf_right->add(triangles[6]);
//        surf_right->add(triangles[7]);
//        gmds::CellGroup<gmds::Face>* surf_bottom = mesh.newGroup<gmds::Face>("surf_bottom");
//        surf_bottom->add(triangles[8]);
//        surf_bottom->add(triangles[9]);
//        gmds::CellGroup<gmds::Face>* surf_top = mesh.newGroup<gmds::Face>("surf_top");
//        surf_top->add(triangles[10]);
//        surf_top->add(triangles[11]);
//
//        AFacGeomManager->updateFromMesh();

        std::vector<gmds::cad::GeomPoint*> vertices;
        std::vector<gmds::cad::GeomCurve*> curves;
        std::vector<gmds::cad::GeomSurface*> surfaces;
        AFacGeomManager->getPoints(vertices);
        AFacGeomManager->getCurves(curves);
        AFacGeomManager->getSurfaces(surfaces);
        std::map<std::string, std::uintptr_t> name2Points;
        std::map<std::string, std::uintptr_t> name2Curves;
        std::map<std::string, std::uintptr_t> name2Surfaces;
        for(auto v: vertices) {
            name2Points[v->getName()] = reinterpret_cast<std::uintptr_t> (v);
        }
        for(auto c: curves) {
            name2Curves[c->getName()] = reinterpret_cast<std::uintptr_t> (c);
        }
        for(auto s: surfaces) {
            name2Surfaces[s->getName()] = reinterpret_cast<std::uintptr_t> (s);
        }

        double tol = gmds::math::Point(minXYZ[0],minXYZ[1],minXYZ[2]).distance(gmds::math::Point(maxXYZ[0],maxXYZ[1],maxXYZ[2])) * 10e-14;

        typedef enum {
            LEFT = 1		  , RIGHT = 1 << 1  , FRONT = 1 << 2  , BACK = 1 << 3  ,
            BOTTOM = 1 << 4  , TOP = 1 << 5
        } sides;

        for(int i=0; i<nbNodes; i++) {

            kmds::TCellID nid = nodeIDs.get(i);
            kmds::Node n = AMesh->getNode(nid);
            kmds::TCoord xyz[3];
            n.getLocation(xyz[0], xyz[1], xyz[2]);

            gmds::math::Point pt(xyz[0], xyz[1], xyz[2]);
            int pos = 0;

            if(pt.X() < minXYZ[0]+tol ) {
                pos = pos|LEFT;
            }
            if(pt.X() > maxXYZ[0]-tol ) {
                pos = pos|RIGHT;
            }
            if(pt.Y() < minXYZ[1]+tol ) {
                pos = pos|FRONT;
            }
            if(pt.Y() > maxXYZ[1]-tol ) {
                pos = pos|BACK;
            }
            if(pt.Z() < minXYZ[2]+tol ) {
                pos = pos|BOTTOM;
            }
            if(pt.Z() > maxXYZ[2]-tol ) {
                pos = pos|TOP;
            }

            // associate to surfaces
            if((pos|LEFT) == pos) {
                (*AVarNodeGeomAssociation)[nid] = name2Surfaces.at("surf_left");
            }
            if((pos|RIGHT) == pos) {
                (*AVarNodeGeomAssociation)[nid] = name2Surfaces.at("surf_right");
            }
            if((pos|FRONT) == pos) {
                (*AVarNodeGeomAssociation)[nid] = name2Surfaces.at("surf_front");
            }
            if((pos|BACK) == pos) {
                (*AVarNodeGeomAssociation)[nid] = name2Surfaces.at("surf_back");
            }
            if((pos|BOTTOM) == pos) {
                (*AVarNodeGeomAssociation)[nid] = name2Surfaces.at("surf_bottom");
            }
            if((pos|TOP) == pos) {
                (*AVarNodeGeomAssociation)[nid] = name2Surfaces.at("surf_top");
            }

            // associate to curves
            if(((pos|BOTTOM) == pos) && ((pos|FRONT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_0");
            }
            if(((pos|BOTTOM) == pos) && ((pos|RIGHT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_1");
            }
            if(((pos|BOTTOM) == pos) && ((pos|BACK) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_2");
            }
            if(((pos|BOTTOM) == pos) && ((pos|LEFT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_3");
            }
            if(((pos|TOP) == pos) && ((pos|FRONT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_4");
            }
            if(((pos|TOP) == pos) && ((pos|RIGHT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_5");
            }
            if(((pos|TOP) == pos) && ((pos|BACK) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_6");
            }
            if(((pos|TOP) == pos) && ((pos|LEFT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_7");
            }
            if(((pos|LEFT) == pos) && ((pos|FRONT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_8");
            }
            if(((pos|FRONT) == pos) && ((pos|RIGHT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_9");
            }
            if(((pos|RIGHT) == pos) && ((pos|BACK) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_10");
            }
            if(((pos|BACK) == pos) && ((pos|LEFT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Curves.at("curve_11");
            }

            // associate to vertices
            if(((pos|BOTTOM) == pos) && ((pos|FRONT) == pos) && ((pos|LEFT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Points.at("corner_0");
            }
            if(((pos|BOTTOM) == pos) && ((pos|FRONT) == pos) && ((pos|RIGHT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Points.at("corner_1");
            }
            if(((pos|BOTTOM) == pos) && ((pos|BACK) == pos) && ((pos|RIGHT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Points.at("corner_2");
            }
            if(((pos|BOTTOM) == pos) && ((pos|BACK) == pos) && ((pos|LEFT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Points.at("corner_3");
            }
            if(((pos|TOP) == pos) && ((pos|FRONT) == pos) && ((pos|LEFT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Points.at("corner_4");
            }
            if(((pos|TOP) == pos) && ((pos|FRONT) == pos) && ((pos|RIGHT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Points.at("corner_5");
            }
            if(((pos|TOP) == pos) && ((pos|BACK) == pos) && ((pos|RIGHT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Points.at("corner_6");
            }
            if(((pos|TOP) == pos) && ((pos|BACK) == pos) && ((pos|LEFT) == pos)) {
                (*AVarNodeGeomAssociation)[nid] = name2Points.at("corner_7");
            }
        }
    }
/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
