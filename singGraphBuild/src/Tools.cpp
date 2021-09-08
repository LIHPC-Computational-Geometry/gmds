/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux (2015)
 *
 * franck.ledoux@cea.fr
 *
 * The FRAME software is a computer program whose purpose is to provide a set
 * of algorithms to build 2D and 3D meshes using frame field concept. The main
 * focus of these algorithms is quadrilateral and hexahedral meshing.
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
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*----------------------------------------------------------------------------*/
/*
 * Tools.cpp
 *
 *  Created on: Sept. 18, 2015
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/math/Constants.h>
#include <gmds/math/Line.h>
#include <gmds/math/Numerics.h>
#include <gmds/math/Plane.h>
#include <gmds/math/Quaternion.h>
#include <gmds/math/Ray.h>
#include <gmds/math/Segment.h>
#include <gmds/singGraphBuild/Tools.h>
/*----------------------------------------------------------------------------*/
#include <set>
#include <sstream>
#include <utility>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace std;
/*----------------------------------------------------------------------------*/
int
Tools::removeBoundarySlivers(Mesh *AMesh)
{
	/* If the normals of the two facets that are on the boundary  is larger
	 * than this threshold, then the tet is considered as a sliver */
	//
	//
	//        vector<index_t> to_remove(mesh.cells.nb(), 0);
	//        index_t nb_slivers_on_border = 0;
	//
	//
	//        for(index_t t=0; t<mesh.cells.nb(); ++t) {
	//            int nb_f_on_border =
	//            (mesh.cells.tet_adjacent(t,0) == NO_CELL) +
	//            (mesh.cells.tet_adjacent(t,1) == NO_CELL) +
	//            (mesh.cells.tet_adjacent(t,2) == NO_CELL) +
	//            (mesh.cells.tet_adjacent(t,3) == NO_CELL) ;
	//            if(nb_f_on_border == 2) {
	//                vec3 normals[2];
	//                index_t cur_normal=0;
	//                for(index_t f=0; f<4; ++f) {
	//                    if(mesh.cells.tet_adjacent(t,f) == NO_CELL) {
	//                        const vec3& p1 = Geom::mesh_vertex(
	//                                                           mesh, mesh.cells.tet_facet_vertex(t,f,0)
	//                                                           );
	//                        const vec3& p2 = Geom::mesh_vertex(
	//                                                           mesh, mesh.cells.tet_facet_vertex(t,f,1)
	//                                                           );
	//                        const vec3& p3 = Geom::mesh_vertex(
	//                                                           mesh, mesh.cells.tet_facet_vertex(t,f,2)
	//                                                           );
	//                        normals[cur_normal]=normalize(cross(p2 - p1, p3 - p1));
	//                        ++cur_normal;
	//                    }
	//                }
	//                if(dot(normals[0],normals[1]) >= sliver_cosangle) {
	//                    to_remove[t] = 1;
	//                    ++nb_slivers_on_border;
	//                }
	//            }
	//        }
	//    }
	//
	//
	gmds::MeshModel new_mod(DIM3 | R | N | F | F2N | R2N | R2F | F2R);
	AMesh->changeModel(new_mod);

	MeshDoctor doc(AMesh);
	doc.buildFacesAndR2F();
	doc.updateUpwardConnectivity();

	std::set<Region> slivers;
	for (auto f_id : AMesh->faces()) {
		Face f = AMesh->get<Face>(f_id);
		std::vector<TCellID> reg_f = f.getIDs<Region>();
		if (reg_f.size() != 1) continue;

		// so f is a boundary face.
		std::vector<Node> n = f.get<Node>();
		// We get is only adjacent tetrahedral element
		Region r = f.get<Region>()[0];

		Node n_opp;
		std::vector<Node> r_nodes = r.get<Node>();
		auto found_opp = false;
		for (int i = 0; i < 4 && !found_opp; i++) {
			if (r_nodes[i].id() != n[0].id() && r_nodes[i].id() != n[1].id() && r_nodes[i].id() != n[2].id()) {
				n_opp = r_nodes[i];
				found_opp = true;
			}
		}
		if (!found_opp) throw GMDSException("Topological error in bnd sliver deletion");

		math::Point p = n_opp.getPoint();
		math::Plane pl(n[0].getPoint(), n[1].getPoint(), n[2].getPoint());
		math::Point pr = pl.project(p);
		math::Vector3d v01(n[0].getPoint(), n[1].getPoint());
		math::Vector3d v02(n[0].getPoint(), n[2].getPoint());
		math::Vector3d v12(n[1].getPoint(), n[2].getPoint());

		double min_dist = math::min3(v01.norm(), v02.norm(), v12.norm());

		if (p.distance(pr) < 0.1 * min_dist) slivers.insert(r);
	}
	std::set<TCellID> faces_to_remove;
	for (auto s : slivers) {
		std::vector<Node> n = s.get<Node>();
		for (auto ni : n) {
			ni.remove(s);
		}

		std::vector<Face> f = s.get<Face>();
		for (auto fi : f) {
			std::vector<TCellID> reg_fi = fi.getIDs<Region>();
			if (reg_fi.size() < 2) {
				faces_to_remove.insert(fi.id());
			}

			fi.remove(s);
		}

		AMesh->deleteRegion(s);
	}

	for (auto fi : faces_to_remove) {
		Face fdi = AMesh->get<Face>(fi);

		AMesh->deleteFace(fi);
	}

	return (!slivers.empty());
}
/*----------------------------------------------------------------------------*/
Tools::Tools(Mesh *AMesh, Variable<math::Cross2D> *AField, Variable<math::AxisAngleRotation> *ARotField) :
  m_mesh(AMesh), m_field(AField), m_rot_field(ARotField)
{
}
/*----------------------------------------------------------------------------*/
math::Chart
Tools::computeChartIn(const math::Point &APnt, const Face &AFace)
{
	std::vector<Node> n = AFace.get<Node>();
	//===================================================================
	// STEP 1 - We compute the location of APnt into AFace
	//===================================================================
	double coeff[3] = {0, 0, 0};

	math::Point::computeBarycentric(n[0].getPoint(), n[1].getPoint(), n[2].getPoint(), APnt, coeff[0], coeff[1], coeff[2]);

	//===================================================================
	// STEP 2 - We extract the quaternion representation of each frame  at
	//         the face corners
	//===================================================================
	std::vector<math::AxisAngleRotation> r;
	r.resize(3);
	for (int i = 0; i < 3; i++) {
		r[i] = (*m_rot_field)[n[i].id()];
	}

	std::vector<math::Quaternion> qs;
	qs.resize(3);
	for (int i = 0; i < 3; i++) {
		qs[i] = math::Quaternion(r[i].toChart());
	}

	std::vector<TCoord> ws;
	ws.resize(3);
	ws[0] = coeff[0];
	ws[1] = coeff[1];
	ws[2] = coeff[2];

	//===================================================================
	// STEP 3 - We compute the mean quaternion and return the corresponding
	//         chart
	//===================================================================
	math::Quaternion q = math::Quaternion::mean(qs, ws);

	return math::Chart(q);
}
/*-----------------------------------------------------------------------------*/

bool
Tools::isPntInTri(const math::Point &AP, const Face &ATri, bool &AOnEdge0, bool &AOnEdge1, bool &AOnEdge2, double &lambda1, double &lambda2)
{
	std::vector<Node> n = ATri.get<Node>();

	math::Point pnt[3] = {n[0].getPoint(), n[1].getPoint(), n[2].getPoint()};

	math::Plane plane(pnt[0], pnt[1], pnt[2]);
	math::Point p = plane.project(AP);

	// We look if p is insie or outside of the triangle defined by ATri
	double det;
	det = ((pnt[1].Y() - pnt[2].Y()) * (pnt[0].X() - pnt[2].X()) + (pnt[2].X() - pnt[1].X()) * (pnt[0].Y() - pnt[2].Y()));
	lambda1 = ((pnt[1].Y() - pnt[2].Y()) * (p.X() - pnt[2].X()) + (pnt[2].X() - pnt[1].X()) * (p.Y() - pnt[2].Y()));
	lambda1 = lambda1 / det;

	lambda2 = (double) ((pnt[2].Y() - pnt[0].Y()) * (p.X() - pnt[2].X()) + (pnt[0].X() - pnt[2].X()) * (p.Y() - pnt[2].Y()));
	lambda2 = (double) lambda2 / det;

	if (lambda1 >= 0 && lambda2 >= 0 && (lambda1 + lambda2) <= 1 && lambda1 <= 1 && lambda2 <= 1) {
		if (lambda1 == 0) AOnEdge0 = true;
		if (lambda2 == 0) AOnEdge1 = true;
		if ((lambda1 + lambda2) == 0) AOnEdge2 = true;
		return true;
	}
	return false;
}

bool
Tools::isPntInTri(const math::Point &AP, const Face &ATri)
{
	std::vector<Node> n = ATri.get<Node>();

	math::Point pnt[3] = {n[0].getPoint(), n[1].getPoint(), n[2].getPoint()};

	math::Plane plane(pnt[0], pnt[1], pnt[2]);
	math::Point p = plane.project(AP);

	// We look if p is insie or outside of the triangle defined by ATri
	double det;
	det = ((pnt[1].Y() - pnt[2].Y()) * (pnt[0].X() - pnt[2].X()) + (pnt[2].X() - pnt[1].X()) * (pnt[0].Y() - pnt[2].Y()));

	double lambda1 = ((pnt[1].Y() - pnt[2].Y()) * (p.X() - pnt[2].X()) + (pnt[2].X() - pnt[1].X()) * (p.Y() - pnt[2].Y()));
	lambda1 = lambda1 / det;

	double lambda2 = (double) ((pnt[2].Y() - pnt[0].Y()) * (p.X() - pnt[2].X()) + (pnt[0].X() - pnt[2].X()) * (p.Y() - pnt[2].Y()));
	lambda2 = lambda2 / det;

	if (lambda1 >= 0 && lambda2 >= 0 && (lambda1 + lambda2) <= 1 && lambda1 <= 1 && lambda2 <= 1) {

		return true;
	}
	return false;
}
/*----------------------------------------------------------------------------*/
void
Tools::computeFuzzyHeuns(const math::Point &AFromPnt,
                         const math::Vector3d &AFromDir,
                         const std::vector<Face> &AFaces,
                         const std::vector<math::Triangle> &ATri,
                         math::Point &AToPnt,
                         math::Vector3d &AToDir,
                         int &AToFaceId)
{
	// check whether the line intersects the triangle
	math::Line ray(AFromPnt, math::Vector3d(AFromDir.X(), AFromDir.Y(), AFromDir.Z()));

	//===================================================================
	double param[4] = {-1, -1, -1, -1};
	double temp;
	math::Point p[4];
	for (auto i = 0; i < 4; i++) {
		math::Plane pli = ATri[i].getPlaneIncluding();
		if (ray.intersect3D(pli, p[i], param[i])) {
			// intersect the plane, but the triangle??
			// before computing the bar coordinate, we eliminate intersection
			// at the infinity almost due to almost parallel ray and plane
			std::vector<Node> n = AFaces[i].get<Node>();
			/* math::Point pn[3] ={
			    n[0].getPoint(),
			    n[1].getPoint(),
			    n[2].getPoint()};*/
			// We compute the barycentric coords.
			bool on_edge[3] = {false, false, false};
			if (!isPntInTri(p[i], AFaces[i], on_edge[0], on_edge[1], on_edge[2], temp, temp)) {
				param[i] = -1;
			};
			/* if(!isIn(p[i],AFaces[i], on_edge[0],on_edge[1],on_edge[2])){
			         param[i]=-1;
			     };*/
			//            double coeff[3]={0, 0, 0};

			//            p[i]=pli.project(p[i]);
			//            math::Point::computeBarycentric(pn[0], pn[1], pn[2], p[i],
			//                                            coeff[0],coeff[1],
			//                                            coeff[2]);
			//
			//            if(coeff[0]<0 || coeff[1]<0 || coeff[2]<0 ){
			//                param[i]=-1;
			//            }
		}
		else {
			param[i] = -1;
		}
	}
	double best_param = param[0];
	auto out_index = 0;

	for (auto i = 1; i < 4; i++) {
		if (param[i] > best_param) {
			out_index = i;
			best_param = param[i];
		}
	}
	if (best_param <= 1e-8) throw GMDSException("Tools::computeFuzzyHeuns: No out face (1)");
	AToPnt = p[out_index];

	//===================================================================
	// compute the frame in out_pnt
	math::Chart ci = computeChartIn(AToPnt, AFaces[out_index]);

	//===================================================================
	// among the 6 vectors of ci, we take the one which is the
	// best aligned with dir and we start the process a second
	// time
	math::Vector3d ci_vectors[6] = {ci.X(), -ci.X(), ci.Y(), -ci.Y(), ci.Z(), -ci.Z()};
	math::Vector3d heuns_corr = ci_vectors[0];
	double best_align_dot = AFromDir.dot(ci_vectors[0]);

	for (int i = 0; i < 6; i++) {
		if (best_align_dot < AFromDir.dot(ci_vectors[i])) {
			heuns_corr = ci_vectors[i];
			best_align_dot = AFromDir.dot(ci_vectors[i]);
		}
	}

	AToDir = heuns_corr;
	AToFaceId = out_index;
}

/*----------------------------------------------------------------------------*/
void
Tools::computeFuzzyHeuns(const math::Point &AFromPnt,
                         const math::Vector3d &AFromDir,
                         const std::vector<Edge> &AEdges,
                         const std::vector<math::Segment> &ASeg,
                         math::Point &AToPnt,
                         math::Vector3d &AToDir,
                         int &AToFaceId)
{
	//===================================================================
	math::Point ray_from = AFromPnt;
	math::Point ray_to = AFromPnt + math::Vector3d(AFromDir.X(), AFromDir.Y(), AFromDir.Z());
	double param[3] = {-1, -1, -1};
	double param_seg[3] = {-1, -1, -1};
	math::Point p[3];
	for (auto i = 0; i < 3; i++) {
		math::Point s0 = ASeg[i].getPoint(0);
		math::Point s1 = ASeg[i].getPoint(1);
		math::Point ray_to_i = math::Plane(ray_from, s0, s1).project(ray_to);
		math::Ray ray(ray_from, ray_to_i);

		ray.intersect3D(ASeg[i], p[i], param_seg[i], param[i]);
		//        if(ok){
		//            std::cout<<"INTERSECT Edge "<<AEdges[i].getIDs<Node>()[0]<<" - "
		//            <<AEdges[i].getIDs<Node>()[1]<<": "<<param_seg[i]<<"(s), "
		//            <<param[i]<<"(r)"<<std::endl;
		//        }
		//        else{
		//            std::cout<<"NOT INTERSECT Edge "<<AEdges[i].getIDs<Node>()[0]<<" - "
		//            <<AEdges[i].getIDs<Node>()[1]<<": "<<param_seg[i]<<"(s), "
		//            <<param[i]<<"(r)"<<std::endl;
		//
		//        }
	}
	double best_param = param[0];
	auto out_index = 0;

	for (auto i = 1; i < 3; i++) {
		if (param[i] > best_param) {
			out_index = i;
			best_param = param[i];
		}
	}
	if (best_param <= 1e-8) {

		for (int i = 0; i < 3; i++) {
			std::cout << "Edge " << AEdges[i].getIDs<Node>()[0] << " - " << AEdges[i].getIDs<Node>()[1] << ": " << param_seg[i] << "(s), " << param[i] << "(r)"
			          << std::endl;
		}
		throw GMDSException("Tools::computeFuzzyHeuns: No out face (2)");
	}
	double out_param = param_seg[out_index];
	AToPnt = p[out_index];
	//===================================================================
	// compute the frame in out_pnt
	Edge e = AEdges[out_index];
	std::vector<TCellID> ne = e.getIDs<Node>();
	math::AxisAngleRotation r[2] = {(*m_rot_field)[ne[0]], (*m_rot_field)[ne[1]]};

	std::vector<math::Quaternion> qs;
	qs.resize(2);
	for (int i = 0; i < 2; i++) {
		qs[i] = math::Quaternion(r[i].toChart());
	}

	std::vector<TCoord> ws;
	ws.resize(2);
	ws[0] = out_param;
	ws[1] = 1 - out_param;

	math::Quaternion q = math::Quaternion::mean(qs, ws);

	math::Chart ci = math::Chart(q);
	//===================================================================
	// among the 6 vectors of ci, we take the one which is the
	// best aligned with dir and we start the process a second
	// time
	math::Vector3d ci_vectors[6] = {ci.X(), -ci.X(), ci.Y(), -ci.Y(), ci.Z(), -ci.Z()};
	math::Vector3d heuns_corr = ci_vectors[0];
	double best_align_dot = AFromDir.dot(ci_vectors[0]);

	for (int i = 0; i < 6; i++) {
		if (best_align_dot < AFromDir.dot(ci_vectors[i])) {
			heuns_corr = ci_vectors[i];
			best_align_dot = AFromDir.dot(ci_vectors[i]);
		}
	}

	AToDir = heuns_corr;
	AToFaceId = out_index;
}
/*----------------------------------------------------------------------------*/
bool
Tools::followFlow(const PointVolumetricData &AData, const double AMaxDist, math::Point &APnt)
{
	//   std::cout<<"================FLOW VOLUME ================"<<std::endl;
	double current_dist = 0;
	math::Vector3d dir(AData.dir);
	math::Point from = AData.pnt;
	gmds::Region current_cell = AData.tet;

	while (current_dist < AMaxDist) {
		//        if(isFFSingular(current_cell))
		//            return false;

		// to avoid starting in the wrong being on an edge
		from = 0.95 * from + 0.05 * current_cell.center();

		//================================================
		// We build the triangular shapes corresponding to each
		// tet face
		std::vector<Face> cf = current_cell.get<Face>();
		//        std::cout<<"------------------"<<std::endl;
		//        std::cout<<"From point "<<from<<std::endl;
		//        math::Vector3d dir_v(dir.X(),dir.Y(), dir.Z());
		//        std::cout<<"        to "<<from+dir_v<<std::endl;
		//        std::cout<<"Region "<<current_cell.id()<<std::endl;

		std::vector<math::Triangle> t;
		t.resize(4);
		for (unsigned int i = 0; i < 4; i++) {
			std::vector<Node> cn = cf[i].get<Node>();
			t[i] = math::Triangle(cn[0].getPoint(), cn[1].getPoint(), cn[2].getPoint());
		}
		//================================================
		// PHASE 1
		//================================================
		math::Point out_pnt;
		math::Vector3d out_dir;
		int out_index = -1;
		computeFuzzyHeuns(from, dir, cf, t, out_pnt, out_dir, out_index);
		//================================================
		// PHASE 2
		//================================================
		dir = 0.5 * dir + 0.5 * out_dir;
		computeFuzzyHeuns(from, dir, cf, t, out_pnt, out_dir, out_index);

		//================================================
		// FINALIZATION
		//================================================
		// Now we have our point to get out of current_cell
		math::Vector3d v(from, out_pnt);
		if (current_dist + v.norm() >= AMaxDist) {
			v.normalize();
			double remaining_dist = AMaxDist - current_dist;

			// we stop in this tet
			APnt = from + remaining_dist * v;
			return true;
		}
		else {
			Face out_face = cf[out_index];
			// std::cout<<"Out face: "<<out_face.id()<<std::endl;
			// Fuzzy approach to avoid topological cases
			math::Point new_from = 0.95 * out_pnt + 0.05 * out_face.center();
			dir = out_dir;
			current_dist += math::Vector3d(from, new_from).norm();
			from = new_from;
			// We compute the next current cell using the face we go
			// from
			std::vector<Region> adj_out_face = out_face.get<Region>();
			if (adj_out_face.size() == 1) {
				// we are on the boudary even if we have not reach the
				// expected distance
				APnt = from;
				return true;
			}
			else {
				if (adj_out_face[0].id() == current_cell.id())
					current_cell = adj_out_face[1];
				else
					current_cell = adj_out_face[0];
			}
			// std::cout<<"\t -->next region: "<<current_cell.id()<<std::endl;
			// not the same current dist due to the fuzzy displacement
			if (current_dist >= AMaxDist) {
				// we stop here so
				APnt = from;
				return true;
			}
		}     // else
	};

	throw GMDSException("Error in Tools::followFlow");
}
/*---------------------------------------------------------------------------*/
math::Chart::Mapping
Tools::getRij(const TCellID AFrom, const TCellID ATo) const
{
	math::AxisAngleRotation rotation_from = (*m_rot_field)[AFrom];
	math::AxisAngleRotation rotation_to = (*m_rot_field)[ATo];

	math::Chart chart_from = rotation_from.toChart();
	math::Chart chart_to = rotation_to.toChart();

	return math::Chart::Mapping(chart_from, chart_to);
}
/*---------------------------------------------------------------------------*/
bool
Tools::isFFSingular(const Face &AF)
{
	std::vector<TCellID> n = AF.getIDs<Node>();
	math::Chart::Mapping m01 = getRij(n[0], n[1]);
	math::Chart::Mapping m12 = getRij(n[1], n[2]);
	math::Chart::Mapping m20 = getRij(n[2], n[0]);

	return !(m20 * m12 * m01).isIdentity();
}
/*---------------------------------------------------------------------------*/
bool
Tools::isFFSingular(const Region &AR)
{
	std::vector<Face> f = AR.get<Face>();
	return (isFFSingular(f[0]) || isFFSingular(f[1]) || isFFSingular(f[2]) || isFFSingular(f[3]));
}
/*----------------------------------------------------------------------------*/
bool
Tools::followFlow(const PointSurfacicData &AData, const double AMaxDist, const int AMarkEdgeOnCurve, math::Point &APnt)
{
	//    std::cout<<"==================== FLOW SURF ===================="<<std::endl;

	double current_dist = 0;
	math::Vector3d dir(AData.dir);
	math::Point from = AData.pnt;
	gmds::Face current_cell = AData.tri;

	while (current_dist < AMaxDist) {
		// std::cout<<"************************"<<std::endl;
		// if(isFFSingular(current_cell))
		// return false;
		// to avoid starting in the wrong being on an edge
		from = 0.95 * from + 0.05 * current_cell.center();

		// out_dir lives in the plan of the face we come from, while
		// dir must live in the plan of the next face, now current_cell
		math::Point p(from.X() + dir.X(), from.Y() + dir.Y(), from.Z() + dir.Z());
		std::vector<Node> cur_nodes = current_cell.get<Node>();
		math::Plane pl(cur_nodes[0].getPoint(), cur_nodes[1].getPoint(), cur_nodes[2].getPoint());
		math::Point proj = pl.project(p);
		dir = math::Vector3d(from, proj);

		// std::cout<<"From point "<<from<<std::endl;
		// math::Vector3d dir_v(dir.X(),dir.Y(), dir.Z());
		// std::cout<<"        to "<<from+dir_v<<std::endl;
		// std::cout<<"Face "<<current_cell.id()<<std::endl;
		//================================================
		// We build the triangular shapes corresponding to each
		// tet face
		std::vector<Edge> cf = current_cell.get<Edge>();
		std::vector<math::Segment> s;
		s.resize(3);
		for (unsigned int i = 0; i < 3; i++) {
			std::vector<Node> cn = cf[i].get<Node>();
			s[i] = math::Segment(cn[0].getPoint(), cn[1].getPoint());
		}
		//================================================
		// PHASE 1
		//================================================
		math::Point out_pnt;
		math::Vector3d out_dir;
		int out_index = -1;
		computeFuzzyHeuns(from, dir, cf, s, out_pnt, out_dir, out_index);
		//================================================
		// PHASE 2
		//================================================
		dir = 0.5 * dir + 0.5 * out_dir;
		math::Point p2(from.X() + dir.X(), from.Y() + dir.Y(), from.Z() + dir.Z());
		dir = math::Vector3d(from, pl.project(p2));

		computeFuzzyHeuns(from, dir, cf, s, out_pnt, out_dir, out_index);

		//================================================
		// FINALIZATION
		//================================================
		// Now we have our point to get out of current_cell
		math::Vector3d v(from, out_pnt);
		// std::cout<<"OUT PNT: "<<out_pnt<<std::endl;
		if (current_dist + v.norm() >= AMaxDist) {
			v.normalize();
			double remaining_dist = AMaxDist - current_dist;

			// we stop in this tet
			APnt = from + remaining_dist * v;
			return true;
		}
		else {
			Edge out_edge = cf[out_index];
			// Fuzzy approach to avoid topological cases
			from = 0.95 * out_pnt + 0.05 * out_edge.center();
			current_dist += v.norm();
			// We compute the next current cell using the face we go
			// from
			std::vector<Face> adj_out_edge = out_edge.get<Face>();
			std::vector<Face> adj_bnd;
			for (auto f : adj_out_edge) {
				std::vector<Region> adj_f = f.get<Region>();
				if (adj_f.size() == 1) adj_bnd.push_back(f);
			}
			if (adj_bnd.size() == 1 || m_mesh->isMarked(out_edge, AMarkEdgeOnCurve)) {
				// we are on the boudary even if we have not reach the
				// expected distance
				APnt = from;
				return true;
			}
			else {
				if (adj_bnd[0].id() == current_cell.id())
					current_cell = adj_bnd[1];
				else
					current_cell = adj_bnd[0];
			}
			dir = out_dir;
		}
	};

	throw GMDSException("Error in Tools::followFlow");
}
/*----------------------------------------------------------------------------*/
void
Tools::traverseTriangle(const Face &AFace,
                        const math::Point &AInPnt,
                        const math::Vector3d &AInVec,
                        const int AInCellDim,
                        const TCellID AInCellID,
                        math::Point &AOutPnt,
                        math::Vector3d &AOutVec,
                        int &AOutCellDim,
                        TCellID &AOutCellID,
                        double &streamlineDeviation)
{
	if (AInCellDim == 0) {
		Node from_node = m_mesh->get<Node>(AInCellID);
		traverseTriangle(AFace, from_node, AInPnt, AInVec, AOutPnt, AOutVec, AOutCellDim, AOutCellID, streamlineDeviation);
	}
	else if (AInCellDim == 1) {
		Edge from_edge = m_mesh->get<Edge>(AInCellID);
		traverseTriangle(AFace, from_edge, AInPnt, AInVec, AOutPnt, AOutVec, AOutCellDim, AOutCellID, streamlineDeviation);
	}
	else
		throw GMDSException("findNextCell: we can only come from a node or an edge");
}

/*----------------------------------------------------------------------------*/
void
Tools::traverseTriangle(const Face &AFace,
                        const Node &ANode,
                        const math::Point &AInPnt,
                        const math::Vector3d &AInVec,
                        math::Point &AOutPnt,
                        math::Vector3d &AOutVec,
                        int &AOutCellDim,
                        TCellID &AOutCellID,
                        double &streamlineDeviation)
{
	Node node_from = ANode;
	//=====================================================================
	// We look for the opposite edge
	//=====================================================================
	Edge opposite_edge;
	std::vector<Edge> current_edges = AFace.get<Edge>();
	for (unsigned int i = 0; i < current_edges.size(); i++) {
		Edge ei = current_edges[i];
		std::vector<Node> ei_nodes = ei.get<Node>();
		if (node_from != ei_nodes[0] && node_from != ei_nodes[1]) opposite_edge = ei;
	}

	//=====================================================================
	// Other nodes are the nodes of the opposite edge
	//=====================================================================
	std::vector<Node> other_nodes = opposite_edge.get<Node>();

	//=====================================================================
	// STEP 2: We compute a first out pnt and vector
	//=====================================================================
	//  std::cout << "========================= Heun's 1 ===================="
	//	    << std::endl;
	math::Point out_pnt;
	math::Vector3d out_vec;
	math::Vector3d in_vec = AInVec;
	in_vec.normalize();
	double deviation = 0.0;
	int out_id;

	out_id = heunsComputation(     // node_from,
	   AInPnt, in_vec, other_nodes[0], other_nodes[1], opposite_edge, out_pnt, out_vec, deviation);

	// cout<<"out_id "<<out_id<<endl;
	/*
  int out_id = RK4Computation(node_from,
	                            AInPnt, in_vec,
	                            other_nodes[0],
	                            other_nodes[1],
	                            opposite_edge,
	                            out_pnt, out_vec, deviation, AFace, stepSize);*/

	//  std::cout << "First computed out point: "<< out_pnt << std::endl;
	//  std::cout << "First computed out vect: " << out_vec << std::endl;
	//=====================================================================
	// STEP 3: We recompute out pnt and vector using Heun's method
	//=====================================================================
	//  std::cout << "========================= Heun's 2 ====================" << std::endl;
	//  std::cout << "Recomputation via Heun's method" << std::endl;
	// math::Vector3d out_vec = cross(out_q, ACurrentFace.normal()).closestVector(AINVec);

	// We compute the approximation only if we are always in the face
	if (out_id == 3) {
		math::Vector3d v_Heun = 0.5 * AInVec + 0.5 * out_vec;
		v_Heun.normalize();

		out_id = heunsComputation(     // node_from,
		   AInPnt, v_Heun, other_nodes[0], other_nodes[1], opposite_edge, out_pnt, out_vec, deviation);
		streamlineDeviation = streamlineDeviation + fabs(1.0 - in_vec.dot(v_Heun)) + fabs(1.0 - out_vec.dot(v_Heun));
		// cout<<"streamlineDeviation for edge "<<AInPnt<<" , "<<out_pnt<<endl;
	}
	else {
		streamlineDeviation = streamlineDeviation + fabs(1.0 - in_vec.dot(out_vec));
		// cout<<"streamlineDeviation for edge "<<AInPnt<<" , "<<out_pnt<<endl;
	}
	//  std::cout << "Second computed out point: " << out_pnt << std::endl;

	AOutPnt = out_pnt;
	AOutVec = out_vec;

	if (out_id == 1 || out_id == 2)
		AOutCellDim = 0;
	else if (out_id == 3)
		AOutCellDim = 1;
	else {
		std::cout << "OUT ID: " << out_id << std::endl;
		throw GMDSException("Wrong output cell dimension in Heun's method");
	}

	if (out_id == 1) {
		AOutCellID = other_nodes[0].id();
	}
	else if (out_id == 2) {
		AOutCellID = other_nodes[1].id();
	}
	else if (out_id == 3) {
		AOutCellID = opposite_edge.id();
	}
}
/*----------------------------------------------------------------------------*/
void
Tools::traverseTriangle(const Face &AFace,
                        const Edge &AEdge,
                        const math::Point &AInPnt,
                        const math::Vector3d &AInVec,
                        math::Point &AOutPnt,
                        math::Vector3d &AOutVec,
                        int &AOutCellDim,
                        TCellID &AOutCellID,
                        double &streamlineDeviation)
{
	Edge edge_in = AEdge;
	Edge other_edges[2];
	int nb_edges = 0;
	std::vector<Edge> current_edges = AFace.get<Edge>();
	for (unsigned int i = 0; i < current_edges.size(); i++) {
		if (current_edges[i] != edge_in) other_edges[nb_edges++] = current_edges[i];
	}

	Node other_node;
	std::vector<Node> current_nodes = AFace.get<Node>();
	std::vector<Node> nodes_in = edge_in.get<Node>();
	for (unsigned int i = 0; i < current_nodes.size(); i++) {
		Node current_node = current_nodes[i];
		if (current_node != nodes_in[0] && current_node != nodes_in[1]) other_node = current_node;
	}

	// std::cout << "EDGE IN  : " << nodes_in[0].id() << ", " << nodes_in[1].id() << std::endl;

	std::vector<Node> nodes_out0 = other_edges[0].get<Node>();
	std::vector<Node> nodes_out1 = other_edges[1].get<Node>();
	//  std::cout << "OTHER TRIANGLE NODE: " << other_node.id() << std::endl;

	//=====================================================================
	// STEP 2: We compute a first out pnt and vector
	//=====================================================================
	//  std::cout << "========================= Heun's 1 ===================="
	//	    << std::endl;
	math::Point out_pnt;
	math::Vector3d out_vec;
	math::Vector3d in_vec = AInVec;
	in_vec.normalize();
	double deviation = 0.0;
	int out_id;

	out_id = heunsComputation(     // edge_in,
	   AInPnt, in_vec, other_node, edge_in.get<Node>()[0], edge_in.get<Node>()[1], other_edges[0], other_edges[1], out_pnt, out_vec, deviation);

	// cout<<"out_id "<<out_id<<endl;

	/*
  int out_id = RK4Computation(edge_in,
	                            AInPnt, in_vec,
	                            edge_in.get<Node>()[0],
	                            edge_in.get<Node>()[1],
	                            other_edges[0],
	                            other_edges[1],
	                            out_pnt, out_vec, deviation, AFace, stepSize);*/

	// std::cout << "First computed out point: " << out_pnt << std::endl;
	// std::cout << "First computed out vect: " << out_vec << std::endl;
	//=====================================================================
	// STEP 3: We recompute out pnt and vector using Heun's method
	//=====================================================================
	// std::cout << "========================= Heun's 2 ====================" << std::endl;
	// std::cout << "Recomputation via Heun's method" << std::endl;

	if ((out_id == 1 || out_id == 4 || out_id == 5)) {
		math::Vector3d v_Heun = 0.5 * AInVec + 0.5 * out_vec;
		v_Heun.normalize();

		out_id = heunsComputation(     // edge_in,
		   AInPnt, v_Heun, other_node, edge_in.get<Node>()[0], edge_in.get<Node>()[1], other_edges[0], other_edges[1], out_pnt, out_vec, deviation);
		streamlineDeviation = streamlineDeviation + fabs(1.0 - in_vec.dot(v_Heun)) + fabs(1.0 - out_vec.dot(v_Heun));
	}
	else {
		streamlineDeviation = streamlineDeviation + fabs(1.0 - in_vec.dot(out_vec));
	}

	AOutPnt = out_pnt;
	AOutVec = out_vec;

	if (out_id == 1 || out_id == 2 || out_id == 3)
		AOutCellDim = 0;
	else if (out_id == 4 || out_id == 5)
		AOutCellDim = 1;
	else
		throw GMDSException("Wrong output cell dimension in Heun's method");

	if (out_id == 1) {
		AOutCellID = other_node.id();
	}
	else if (out_id == 2) {
		AOutCellID = edge_in.get<Node>()[0].id();
	}
	else if (out_id == 3) {
		AOutCellID = edge_in.get<Node>()[1].id();
	}
	else if (out_id == 4) {
		AOutCellID = other_edges[0].id();
	}
	else {
		AOutCellID = other_edges[1].id();
	}
	// std::cout<<"OUT CELL: "<<AOutCellDim<<" - "<<AOutCellIconst gmds::Edge&         AInEdge,D<<std::endl;
}
/*----------------------------------------------------------------------------*/
int
Tools::heunsComputation(     // const Edge&         AINEdge,
   const math::Point &AINPnt,
   const math::Vector3d &AINVec,
   const Node &AOPPNode,
   const Node &AINNode1,
   const Node &AINNode2,
   const Edge &AOUTEdge1,
   const Edge &AOUTEdge2,
   math::Point &AOUTPnt,
   math::Vector3d &AOUTVec,
   double &deviation)
{
	math::Point pout = AINPnt + AINVec;

	// std::cout << "===================================" << std::endl;
	// std::cout << "compute OUT cross from edge "<<AINEdge.id() << std::endl;
	// std::cout << "\t AIN Pnt: " << AINPnt << std::endl;
	// std::cout << "\t AIN Vec: " << AINVec << std::endl;
	// std::cout << "\t \t out pnt: " << pout << std::endl;
	// std::cout << "\t Node candidate : " << AOPPNode.id() << std::endl;
	// std::cout << "\t Edge candidate1: "
	// 	    << AOUTEdge1.get<Node>()[0].id() << ",  "
	// 	    << AOUTEdge1.get<Node>()[1].id() << std::endl;
	// std::cout << "\t Edge candidate2: "
	// 	    << AOUTEdge2.get<Node>()[0].id() << ",  "
	// 	    << AOUTEdge2.get<Node>()[1].id() << std::endl;
	// std::cout << "Node " << AOPPNode.id() << " - " << AOPPNode.getPoint() << std::endl;
	// std::cout << "Node " << AINNode1.id() << " - " << AINNode1.getPoint() << std::endl;
	// std::cout << "Node " << AINNode2.id() << " - " << AINNode2.getPoint() << std::endl;

	math::Vector3d v_in = AINVec;
	v_in.normalize();

	//================================================
	// Go through the opposite node AOPPNode???
	//================================================
	math::Point opp_node_loc = AOPPNode.getPoint();
	math::Vector3d v_opp(AINPnt, opp_node_loc);
	v_opp.normalize();

	if (math::near(v_opp.dot(v_in) - 1, 0)) {
		// std::cout << "INTERSECT NODE --> out in a node" << std::endl;
		// intersected point = node
		AOUTPnt = opp_node_loc;
		computeOutVectorAtPoint(AOPPNode, AINVec, AOUTVec);
		deviation = deviation + fabs(1.0 - v_in.dot(v_opp)) + fabs(1.0 - AOUTVec.dot(v_opp));
		return 1;
	}
	//================================================
	// Go through AINNode1???
	//================================================
	math::Point node1_loc = AINNode1.getPoint();
	math::Vector3d v1(AINPnt, node1_loc);
	v1.normalize();

	if (math::near(v1.dot(v_in) - 1, 0)) {
		// std::cout << "INTERSECT NODE --> out in a node" << std::endl;
		// intersected point = node
		AOUTPnt = node1_loc;
		computeOutVectorAtPoint(AINNode1, AINVec, AOUTVec);
		deviation = deviation + fabs(1.0 - v_in.dot(v1)) + fabs(1.0 - AOUTVec.dot(v1));
		return 2;
	}
	//================================================
	// Go through AINNode2???
	//================================================
	math::Point node2_loc = AINNode2.getPoint();
	math::Vector3d v2(AINPnt, node2_loc);
	v2.normalize();

	if (math::near(v2.dot(v_in) - 1, 0)) {
		// std::cout << "INTERSECT NODE --> out in a node" << std::endl;
		// intersected point = node
		AOUTPnt = node2_loc;
		computeOutVectorAtPoint(AINNode2, AINVec, AOUTVec);
		deviation = deviation + fabs(1.0 - v_in.dot(v2)) + fabs(1.0 - AOUTVec.dot(v2));
		return 3;
	}

	//================================================
	// Go through the first edge ???
	// And not through one of its end points due to
	// previous tests.
	//================================================
	bool intersectEdge1 = computeOutVectorFromRayAndEdge(AOUTEdge1, AINPnt, AINVec, AOUTPnt, AOUTVec, deviation);

	// std::cout<<"First intersection: "<<intersectEdge1<<std::endl;
	if (intersectEdge1) return 4;
	//================================================
	// Go through the second edge ???
	// And not through one of its end points due to
	// previous tests.
	//================================================
	bool intersectEdge2 = computeOutVectorFromRayAndEdge(AOUTEdge2, AINPnt, AINVec, AOUTPnt, AOUTVec, deviation);
	// std::cout<<"Second intersection: "<<intersectEdge2<<std::endl;
	if (intersectEdge2) return 5;

	//================================================
	// If we arrive here, it means that the intersection
	// of our ray is out of the triangle
	//
	// Our choice is to take the closest end point
	//================================================
	double value_opp = v_opp.dot(v_in);
	double value_1 = v1.dot(v_in);
	double value_2 = v2.dot(v_in);

	if (value_opp >= value_1 && value_opp >= value_2) {
		// the outpoint is the first one
		AOUTPnt = opp_node_loc;
		computeOutVectorAtPoint(AOPPNode, AINVec, AOUTVec);
		deviation = deviation + fabs(1.0 - value_opp) + fabs(1.0 - AOUTVec.dot(v_opp));
		return 1;
	}
	else if (value_1 >= value_opp && value_1 >= value_2) {
		// the outpoint is the second one
		AOUTPnt = node1_loc;
		computeOutVectorAtPoint(AINNode1, AINVec, AOUTVec);
		deviation = deviation + fabs(1.0 - value_1) + fabs(1.0 - AOUTVec.dot(v1));
		return 2;
	}
	else {
		AOUTPnt = node2_loc;
		computeOutVectorAtPoint(AINNode2, AINVec, AOUTVec);
		deviation = deviation + fabs(1.0 - value_2) + fabs(1.0 - AOUTVec.dot(v2));
		return 3;
	}
	throw GMDSException("ERROR: No out point in Tools::heunsComputation");
}
/*----------------------------------------------------------------------------*/
int
Tools::heunsComputation(     // const Node&         AFROMNode,
   const math::Point &AINPnt,
   const math::Vector3d &AINVec,
   const Node &AOPPNode1,
   const Node &AOPPNode2,
   const Edge &AOPPEdge,
   math::Point &AOUTPnt,
   math::Vector3d &AOUTVec,
   double &deviation)
{
	// std::cout << "===================================" << std::endl;
	// std::cout << "compute OUT cross from node"<<AFROMNode.id() << std::endl;
	// std::cout << "\t AIN Pnt: " << AINPnt << std::endl;
	// std::cout << "\t AIN Vec: " << AINVec << std::endl;
	// std::cout << "\t \t out pnt: " << pout << std::endl;
	// std::cout << "\t Node candidate 1: " << AOPPNode1.id() << std::endl;
	// std::cout << "\t Node candidate 2: " << AOPPNode2.id() << std::endl;
	// std::cout << "\t Edge candidate: "
	// 	    << AOPPEdge.get<Node>()[0].id() << ",  "
	// 	    << AOPPEdge.get<Node>()[1].id() << std::endl;

	math::Vector3d v_in = AINVec;
	v_in.normalize();

	//================================================
	// Go through the opposite node AOPPNode1???
	//================================================
	math::Point opp_node_loc1 = AOPPNode1.getPoint();
	math::Vector3d v_opp1(AINPnt, opp_node_loc1);
	v_opp1.normalize();

	if (math::near(v_opp1.dot(v_in) - 1, 0)) {
		// std::cout << "INTERSECT NODE --> out in a node" << std::endl;
		// intersected point = node
		AOUTPnt = opp_node_loc1;
		computeOutVectorAtPoint(AOPPNode1, AINVec, AOUTVec);
		deviation = deviation + fabs(1.0 - v_in.dot(v_opp1)) + fabs(1.0 - AOUTVec.dot(v_opp1));
		return 1;
	}
	//================================================
	// Go through the opposite node AOPPNode1???
	//================================================
	math::Point opp_node_loc2 = AOPPNode2.getPoint();
	math::Vector3d v_opp2(AINPnt, opp_node_loc2);
	v_opp2.normalize();

	if (math::near(v_opp2.dot(v_in) - 1, 0)) {
		// std::cout << "INTERSECT NODE --> out in a node" << std::endl;
		// intersected point = node
		AOUTPnt = opp_node_loc2;
		computeOutVectorAtPoint(AOPPNode2, AINVec, AOUTVec);
		deviation = deviation + fabs(1.0 - v_in.dot(v_opp2)) + fabs(1.0 - AOUTVec.dot(v_opp2));
		return 2;
	}

	//================================================
	// Go through the opposite edge
	// And not through one of its end points due to
	// previous tests.
	//================================================
	bool intersectEdge = computeOutVectorFromRayAndEdge(AOPPEdge, AINPnt, AINVec, AOUTPnt, AOUTVec, deviation);

	if (intersectEdge) return 3;

	//================================================
	// If we arrive here, it means that the intersection
	// of our ray with the opposite edge is out of the
	// segment defined by AOppEdge
	// Our choice is to take the closest end point
	//================================================
	double value1 = v_opp1.dot(v_in);
	double value2 = v_opp2.dot(v_in);
	if (value1 > value2) {
		// the outpoint is the first one
		AOUTPnt = opp_node_loc1;
		computeOutVectorAtPoint(AOPPNode1, AINVec, AOUTVec);
		deviation = deviation + fabs(1.0 - v_in.dot(v_opp1)) + fabs(1.0 - AOUTVec.dot(v_opp1));
		return 1;
	}
	else {
		// the outpoint is the second one
		AOUTPnt = opp_node_loc2;
		computeOutVectorAtPoint(AOPPNode2, AINVec, AOUTVec);
		deviation = deviation + fabs(1.0 - v_in.dot(v_opp2)) + fabs(1.0 - AOUTVec.dot(v_opp2));
		return 2;
	}

	return 0;
}

/*----------------------------------------------------------------------------*/
void
Tools::RK4Computation(const math::Point &AINPnt,
                      const math::Vector3d &AINVec,
                      math::Point &AOUTPnt,
                      math::Vector3d &AOUTVec,
                      double &deviation,
                      const double &stepSize,
                      bool &AEndOnBdry,
                      std::vector<gmds::TCellID> &ATriangles,
                      int &AToCellDim,
                      gmds::TCellID &AToCellID)
{
	/* ∆t - stepSize
	  * v_in - vector at the intial point p
	  * • calculate the vector k_1 = ∆t * v_in
• determine the vector v_next1, which is the vector corresponding to the point p + k_1/2
• calculate the vector k_2 = ∆t * v_next1
• determine the vector v_next2, which is the vector corresponding to the point p + k_2/2
• calculate the vector k_3 = ∆t * v_next2;
• determine the vector v_next3, which is the vector corresponding to the point p + k_3
• calculate the vector k_4 = ∆t * v_next3
• calculate next_point = p + (k_1 + 2*k_2 + 2*k_3 + k_4 )/6
however, at the last point our directional vector should be taken from the field (given its location)
*/
	std::vector<gmds::TCellID> tempTriangles;
	math::Vector3d v_in = AINVec;
	v_in.normalize();

	math::Vector3d k_1, k_2, k_3, k_4, v_next1, v_next2, v_next3, v_next4;
	math::Point point_1, point_2, point_3;
	AEndOnBdry = false;
	k_1 = stepSize * v_in;
	point_1 = AINPnt + (k_1 / 2.0);
	// it can happen that our start point is already on the bdry - check outside

	if (!AEndOnBdry) {
		bool isFinal = false;
		/*unsigned int previousVisitedFace;
		if(ATriangles.size()==1)
		previousVisitedFace = ATriangles[0];
		else
		previousVisitedFace = ATriangles[ATriangles.size()-2];
		*/
		std::vector<gmds::TCellID> tempTrianglesNew;
		// gmds::TCellID AToCellIDNew = AToCellID;
		// int AToCellDimNew = AToCellDim;
		findTriangleAndNextVectorRK4(AINPnt, point_1, v_in, v_next1, AEndOnBdry, tempTrianglesNew, AToCellDim, AToCellID, isFinal, ATriangles);

		if (!tempTrianglesNew.empty()) {
			tempTriangles.clear();
			tempTriangles.insert(tempTriangles.end(), tempTrianglesNew.begin(), tempTrianglesNew.end());
			// AToCellID = AToCellIDNew;
			// AToCellDim = AToCellDimNew;
		}
		AOUTPnt = point_1;
		AOUTVec = v_next1;
		if (!AEndOnBdry) {
			k_2 = stepSize * v_next1;
			point_2 = AINPnt + (k_2 / 2.0);
			// gmds::TCellID AToCellIDNew = AToCellID;
			// int AToCellDimNew = AToCellDim;
			std::vector<gmds::TCellID> tempTrianglesNew;
			findTriangleAndNextVectorRK4(AINPnt, point_2, v_in, v_next2, AEndOnBdry, tempTrianglesNew, AToCellDim, AToCellID, isFinal, ATriangles);
			if (!tempTrianglesNew.empty()) {
				tempTriangles.clear();
				tempTriangles.insert(tempTriangles.end(), tempTrianglesNew.begin(), tempTrianglesNew.end());
				// AToCellID = AToCellIDNew;
				// AToCellDim = AToCellDimNew;
			}
			AOUTPnt = point_2;
			AOUTVec = v_next2;
			if (!AEndOnBdry) {
				k_3 = stepSize * v_next2;
				point_3 = AINPnt + k_3;
				// gmds::TCellID AToCellIDNew = AToCellID;
				// int AToCellDimNew = AToCellDim;
				std::vector<gmds::TCellID> tempTrianglesNew;

				findTriangleAndNextVectorRK4(AINPnt, point_3, v_in, v_next3, AEndOnBdry, tempTrianglesNew, AToCellDim, AToCellID, isFinal, ATriangles);
				if (!tempTrianglesNew.empty()) {
					tempTriangles.clear();
					tempTriangles.insert(tempTriangles.end(), tempTrianglesNew.begin(), tempTrianglesNew.end());
					// AToCellID = AToCellIDNew;
					// AToCellDim = AToCellDimNew;
				}
				AOUTPnt = point_3;
				AOUTVec = v_next3;
				if (!AEndOnBdry) {

					k_4 = stepSize * v_next3;

					AOUTPnt = AINPnt + k_1 / 6.0 + k_2 / 3.0 + k_3 / 3.0 + k_4 / 6.0;

					isFinal = true;
					// gmds::TCellID AToCellIDNew = AToCellID;
					/// int AToCellDimNew = AToCellDim;
					std::vector<gmds::TCellID> tempTrianglesNew;
					findTriangleAndNextVectorRK4(AINPnt, AOUTPnt, v_in, AOUTVec, AEndOnBdry, tempTrianglesNew, AToCellDim, AToCellID, isFinal, ATriangles);
					if (!tempTrianglesNew.empty()) {
						tempTriangles.clear();
						tempTriangles.insert(tempTriangles.end(), tempTrianglesNew.begin(), tempTrianglesNew.end());
						// AToCellID = AToCellIDNew;
						// AToCellDim = AToCellDimNew;
					}
				}
			}
		}

		if (AOUTVec.isZero()) AOUTVec = v_in;
		ATriangles.insert(ATriangles.end(), tempTriangles.begin(), tempTriangles.end());
		if ((ATriangles.back() != AToCellID) && (AToCellDim == 2)) {
			// cout<<"!!! ATriangles.back()= "<<ATriangles.back()<<" , AToCellID= "<<AToCellID<<endl;
			AToCellID = ATriangles.back();
		}
	}
}

/*----------------------------------------------------------------------------*/
void
Tools::computeOutVectorAtPoint(const Node &ANode, const math::Vector3d &AInVec, math::Vector3d &AOutVec)
{
	math::Cross2D c = (*m_field)[ANode.id()];
	std::vector<math::Vector3d> c_vectors = c.componentVectors();

	int ref_index = 0;
	double ref_dot_product = AInVec.dot(c_vectors[0]);
	for (int i = 1; i < 4; i++) {
		double dot_product_i = AInVec.dot(c_vectors[i]);
		if (dot_product_i > ref_dot_product) {
			ref_dot_product = dot_product_i;
			ref_index = i;
		}
	}
	AOutVec = c_vectors[ref_index];
}

/*----------------------------------------------------------------------------*/
bool
Tools::computeOutVectorFromRayAndEdge(
   const Edge &AEdge, const math::Point &AINPnt, const math::Vector3d &AINVec, math::Point &AOUTPnt, math::Vector3d &AOUTVec, double &deviation)
{
	double epsilon = (double) math::Constants::EPSILON;

	Node n0 = AEdge.get<Node>()[0];
	Node n1 = AEdge.get<Node>()[1];
	math::Segment seg(n0.getPoint(), n1.getPoint());

	math::Point pnt_intersection;
	double param_intersection = 0.0;
	math::Ray ray(AINPnt, AINVec);
	bool intersectEdge = ray.SecondMetIntersect2D(seg, pnt_intersection, param_intersection, epsilon);
	if (!intersectEdge) return false;

	math::Cross2D c0 = (*m_field)[n0.id()];
	math::Cross2D c1 = (*m_field)[n1.id()];
	math::Cross2D meanCross = math::Cross2D::mean(c0, 1 - param_intersection, c1, param_intersection);

	AOUTVec = meanCross.closestComponentVector(AINVec);
	AOUTPnt = pnt_intersection;

	return true;

}

/*----------------------------------------------------------------------------*/
void
Tools::findNextCell(
   const math::Point &AFromPnt, const math::Vector3d &AFromVec, const int AFromCellDim, const TCellID AFromCellID, int &AToCellDim, TCellID &AToCellID)
{
	if (AFromCellDim == 0) {
		Node from_node = m_mesh->get<Node>(AFromCellID);
		findNextCell(AFromPnt, AFromVec, from_node, AToCellDim, AToCellID);
	}
	else if (AFromCellDim == 1) {
		Edge from_edge = m_mesh->get<Edge>(AFromCellID);
		findNextCell(AFromPnt, AFromVec, from_edge, AToCellDim, AToCellID);
	}
	else
		throw GMDSException("findNextCell: we can only arrive from a node or an edge");
}
/*----------------------------------------------------------------------------*/
void
Tools::findNextCell(const math::Point &AFromPnt, const math::Vector3d &AFromVec, const Node &AFromNode, int &AToCellDim, TCellID &AToCellID)
{
	std::vector<Face> adj_faces = AFromNode.get<Face>();

	//==============================================================
	// LOOK FOR A FACE
	//==============================================================
	bool find_next_cell = false;
	Face next_face;
	for (unsigned int i = 0; i < adj_faces.size() && !find_next_cell; i++) {
		Face current_face = adj_faces[i];

		if (isGoingInto(AFromPnt, AFromVec, AFromNode, current_face)) {
			next_face = current_face;
			find_next_cell = true;
		}
	}     // for(unsigned int i=0; i<adj_faces.size(); i++){

	if (!find_next_cell) throw GMDSException("findNextCell: no right face for a node");

	//==============================================================
	// LOOK FOR AN EDGE
	//==============================================================
	// We have found a candidate face, now we check if edges are not
	// a better answer
	bool find_better_edge = false;
	Edge next_edge;
	std::vector<Edge> next_edges = next_face.get<Edge>();
	for (unsigned int i = 0; i < next_edges.size() && !find_better_edge; i++) {
		Edge current_edge = next_edges[i];
		if (isAlong(AFromVec, AFromNode, current_edge)) {
			next_edge = current_edge;
			find_better_edge = true;
		}
	}
	if (find_better_edge) {
		AToCellDim = 1;
		AToCellID = next_edge.id();
	}
	else if (find_next_cell) {
		AToCellDim = 2;
		AToCellID = next_face.id();
	}
}
/*----------------------------------------------------------------------------*/
void
Tools::findNextCell(const math::Point &AFromPnt, const math::Vector3d &AFromVec, const Edge &AFromEdge, int &AToCellDim, TCellID &AToCellID)
{
	std::vector<Face> adj_faces = AFromEdge.get<Face>();
	bool find_next_cell = false;

	for (unsigned int i = 0; i < adj_faces.size() && !find_next_cell; i++) {
		Face current_face = adj_faces[i];

		if (isGoingInto(AFromPnt, AFromVec, AFromEdge, current_face)) {
			find_next_cell = true;
			AToCellDim = 2;
			AToCellID = current_face.id();
		}
	}     // for(unsigned int i=0; i<adj_faces.size(); i++){
}
/*----------------------------------------------------------------------------*/
bool
Tools::isGoingInto(const math::Point &APnt, const math::Vector3d &AVec, const Node &AFromNode, const Face &AFace)
{
	// APnt is located in AFromNode
	double temp_epsilon = (double) math::Constants::EPSILON;
	// We take the nodes of AFace opposite to AFromNode
	std::vector<Node> current_nodes = AFace.get<Node>();
	std::vector<Node> other_nodes;
	for (unsigned int j = 0; j < current_nodes.size(); j++) {
		if (current_nodes[j].id() != AFromNode.id()) other_nodes.push_back(current_nodes[j]);
	}
	if (other_nodes.size() != 2) throw GMDSException("isGoingInto(): Error, I don't find 2 opposite nodes");

	math::Ray from_ray(APnt, AVec);
	math::Segment opp_seg(other_nodes[0].getPoint(), other_nodes[1].getPoint());
	math::Point intersection_pnt;
	double intersection_param;
	return from_ray.SecondMetIntersect2D(opp_seg, intersection_pnt, intersection_param, temp_epsilon);
}
/*----------------------------------------------------------------------------*/
bool
Tools::isGoingInto(const math::Point &APnt, const math::Vector3d &AVec, const Edge &AFromEdge, const Face &AFace)
{
	// APnt is located on AFromEdge

	// We take the node of AFace that is not incident to AFromEdge
	std::vector<Node> current_nodes = AFace.get<Node>();
	std::vector<Node> from_nodes = AFromEdge.get<Node>();
	Node opp_node;
	for (unsigned int j = 0; j < current_nodes.size(); j++) {
		if (current_nodes[j] != from_nodes[0] && current_nodes[j] != from_nodes[1]) opp_node = current_nodes[j];
	}
	if (opp_node.id() == NullID) throw GMDSException("isGoingInto(): Error, no opposite node");

	math::Point opp_pnt = opp_node.getPoint();
	math::Vector3d edge_vector(from_nodes[0].getPoint(), from_nodes[1].getPoint());
	math::Vector3d edge_witness(from_nodes[0].getPoint(), opp_pnt);
	math::Vector3d edge_ortho = edge_vector.cross(math::Vector3d(0, 0, 1));

	if (edge_ortho.dot(edge_witness) < 0) edge_ortho = edge_ortho.opp();

	return (AVec.dot(edge_ortho) >= 0.0);
}
/*----------------------------------------------------------------------------*/
bool
Tools::isAlong(const math::Vector3d &AVec, const Node &AFromNode, Edge &AEdge)
{
	std::vector<Node> n = AEdge.get<Node>();

	Node other_node;
	bool found_origin = false;
	for (unsigned int i = 0; i < n.size(); i++) {
		Node current_node = n[i];
		if (current_node.id() == AFromNode.id())
			found_origin = true;
		else
			other_node = current_node;
	}
	if (!found_origin) return false;

	math::Vector3d v_edge(AFromNode.getPoint(), other_node.getPoint());
	v_edge.normalize();

	math::Vector3d v_dir = AVec;
	v_dir.normalize();

	return math::near(v_dir.dot(v_edge) - 1, 0);
}

/*----------------------------------------------------------------------------*/
void
Tools::getNextVectorRK4(     // const gmds::math::Point& point_1,
   const gmds::Face &AFace,
   vector<double> &lambdas,
   const gmds::math::Vector3d &v_in,
   gmds::math::Vector3d &v_next)
{
	vector<math::Cross2D> faceCrosses(3);
	vector<math::Vector3d> temp_comp_vectors(4);
	vector<gmds::Node> current_nodes = AFace.get<gmds::Node>();

	faceCrosses[0] = (*m_field)[current_nodes[0].id()];
	std::vector<math::Vector3d> c_vectors0 = faceCrosses[0].componentVectors();
	faceCrosses[1] = (*m_field)[current_nodes[1].id()];
	std::vector<math::Vector3d> c_vectors1 = faceCrosses[1].componentVectors();
	faceCrosses[2] = (*m_field)[current_nodes[2].id()];
	std::vector<math::Vector3d> c_vectors2 = faceCrosses[2].componentVectors();
	for (unsigned int i = 0; i < 4; i++) {
		temp_comp_vectors[i] = c_vectors0[i] * lambdas[0] + c_vectors1[i] * lambdas[1] + c_vectors2[i] * lambdas[2];
		// get closest comp vector to v_in
	}

	math::Cross2D cross_next = math::Cross2D::meanNotMedian(faceCrosses, lambdas, 1);
	v_next = cross_next.closestComponentVector(v_in);
}
/*----------------------------------------------------------------------------*/

void
Tools::findTriangleAndNextVectorRK4(const gmds::math::Point &AInPnt,
                                    gmds::math::Point &point_x,
                                    const gmds::math::Vector3d &v_in,
                                    gmds::math::Vector3d &v_next,
                                    bool &AEndOnBdry,
                                    std::vector<gmds::TCellID> &ATriangles,
                                    int &AToCellDim,
                                    gmds::TCellID &AToCellID,
                                    bool &isFinal,
                                    vector<gmds::TCellID> &previousVisitedFaces)
{
	// WARNING if i end up on a node (or very close to one), since i treat only edges -> i limit the possibilities for next cell

	double temp_epsilon = (double) math::Constants::EPSILON;
	bool AOnEdge0, AOnEdge1, AOnEdge2;
	vector<double> lambdas(3);
	AEndOnBdry = false;
	bool withComments = false;
	if (withComments) cout << "findTriangleAndNextVectorRK4" << endl;

	vector<bool> visited_faces(m_mesh->getNbFaces(), false), visited_edges(m_mesh->getNbEdges(), false);

	vector<gmds::Face> current_faces;

	if (AToCellDim == 0) {
		current_faces = m_mesh->get<gmds::Node>(AToCellID).get<gmds::Face>();
	}
	else {
		if (AToCellDim == 1) {
			current_faces = m_mesh->get<gmds::Edge>(AToCellID).get<gmds::Face>();
		}
		else
			current_faces.push_back(m_mesh->get<gmds::Face>(AToCellID));
	}
	if (withComments) cout << "AToCellDim " << AToCellDim << endl;
	bool found_v = false;

	if (previousVisitedFaces.size() == 1)
		visited_faces[previousVisitedFaces[0]] = true;
	else {

		for (unsigned int i = 0; i < previousVisitedFaces.size() - 1; i++)
			visited_faces[previousVisitedFaces[i]] = true;
	}

	if (withComments) {
		cout << "current_faces" << endl;
		for (unsigned int i = 0; i < current_faces.size(); i++)
			cout << current_faces[i].id() << endl;
	}

	while ((!found_v) && (!AEndOnBdry) && (!current_faces.empty())) {     //&&(!current_faces.empty())

		gmds::Face current_face = current_faces.back();

		current_faces.pop_back();

		if (!visited_faces[current_face.id()]) {
			if (withComments) {
				cout << "point_x " << point_x << " InTri " << current_face.id() << " -> "
				     << isPntInTri(point_x, current_face, AOnEdge0, AOnEdge1, AOnEdge2, lambdas[0], lambdas[1]) << endl;
			}
			visited_faces[current_face.id()] = true;
			if (isPntInTri(point_x, current_face, AOnEdge0, AOnEdge1, AOnEdge2, lambdas[0], lambdas[1])) {
				lambdas[2] = 1.0 - lambdas[0] - lambdas[1];
				if (withComments) {
					cout << "point_x " << point_x << " is inside tri " << current_face.id() << endl;
				}
				getNextVectorRK4(/*point_x,*/ current_face, lambdas, v_in, v_next);

				found_v = true;
				if (isFinal) {
					AToCellDim = 2;
					AToCellID = current_face.id();
					if (previousVisitedFaces.back() != current_face.id()) ATriangles.push_back(current_face.id());
				}
				// WARNING here it could be on a bdry edge or node; check outside
			}
			else {     // we got out of the current triangle, check where we are
				if (withComments) {
					cout << "point_x " << point_x << " is NOT!!! inside tri " << current_face.id() << endl;
				}
				math::Vector3d current_vector(AInPnt, point_x);

				math::Ray ray(AInPnt, current_vector);
				// it might be that the next point is situated in a face we cannot reach by edge adjacency
				// while((!found_v)&&(!AEndOnBdry)){

				vector<gmds::Edge> current_edges = current_face.get<gmds::Edge>();
				math::Point pnt_intersection;
				double param_intersection = 0.0;

				for (unsigned int i = 0; i < 3; i++) {
					if (!visited_edges[current_edges[i].id()]) {
						math::Segment seg(current_edges[i].get<Node>()[0].getPoint(), current_edges[i].get<Node>()[1].getPoint());
						if (withComments) {
							cout << "?is ray " << AInPnt << " -> " << point_x << " intersecting any edge of " << current_face.id() << endl;
						}
						// ray.SecondMetIntersect2D(seg, pnt_intersection, param_intersection, m_temp_epsilon);
						if (ray.SecondMetIntersect2D(seg, pnt_intersection, param_intersection, temp_epsilon)) {
							if (withComments) cout << "yes" << endl;
							vector<gmds::Face> adj_faces = current_edges[i].get<gmds::Face>();

							if (adj_faces.size() == 1) {

								AEndOnBdry = true;

								point_x = pnt_intersection;

								current_faces.push_back(adj_faces[0]);
								if (!visited_faces[adj_faces[0].id()]) ATriangles.push_back(adj_faces[0].id());

								if (param_intersection < math::Constants::EPSILON) {
									AToCellDim = 0;     // point
									AToCellID = current_edges[i].get<Node>()[0].id();
									point_x = current_edges[i].get<Node>()[0].getPoint();
								}
								else {
									if (param_intersection > (1.0 - math::Constants::EPSILON)) {
										AToCellDim = 0;     // point
										AToCellID = current_edges[i].get<Node>()[1].id();
										point_x = current_edges[i].get<Node>()[1].getPoint();
									}
									else {
										AToCellDim = 1;
										AToCellID = current_edges[i].id();
									}
								}

								break;
							}
							else {
								if (!visited_faces[adj_faces[0].id()]) ATriangles.push_back(adj_faces[0].id());
								if (!visited_faces[adj_faces[1].id()]) ATriangles.push_back(adj_faces[1].id());
								if ((adj_faces[0] != current_face) && (!visited_faces[adj_faces[0].id()])) {
									// AToCellDim = 2;
									// AToCellID = adj_faces[0].id();
									current_faces.push_back(adj_faces[0]);
									if (withComments) {
										cout << "? is point_x " << point_x << " is inside adj_faces0 " << adj_faces[0].id() << endl;
									}
									if (isPntInTri(point_x, adj_faces[0], AOnEdge0, AOnEdge1, AOnEdge2, lambdas[0], lambdas[1])) {
										lambdas[2] = 1.0 - lambdas[0] - lambdas[1];
										if (withComments) cout << "yes" << endl;
										getNextVectorRK4(/*point_x,*/ adj_faces[0], lambdas, v_in, v_next);
										found_v = true;

										if (isFinal) {
											AToCellDim = 2;
											AToCellID = adj_faces[0].id();
											// ATriangles.push_back(adj_faces[0].id());
										}
										break;
									}
									if (!visited_faces[adj_faces[0].id()]) {
										current_faces.push_back(adj_faces[0]);
									}
									visited_edges[current_edges[i].id()] = true;
								}
								else {
									// AToCellDim = 2;
									// AToCellID = adj_faces[1].id();
									current_faces.push_back(adj_faces[1]);
									if (withComments) {
										cout << "? is point_x " << point_x << " is inside adj_faces1 " << adj_faces[1].id() << endl;
									}
									if (isPntInTri(point_x, adj_faces[1], AOnEdge0, AOnEdge1, AOnEdge2, lambdas[0], lambdas[1])) {
										if (withComments) cout << "yes" << endl;
										lambdas[2] = 1.0 - lambdas[0] - lambdas[1];
										getNextVectorRK4(/*point_x,*/ adj_faces[1], lambdas, v_in, v_next);
										found_v = true;
										if (isFinal) {
											// ATriangles.push_back(adj_faces[1].id());
											AToCellDim = 2;
											AToCellID = adj_faces[1].id();
										}
										break;
									}
									if (!visited_faces[adj_faces[1].id()]) {
										current_faces.push_back(adj_faces[1]);
									}
									visited_edges[current_edges[i].id()] = true;
								}
							}
						}
						visited_edges[current_edges[i].id()] = true;
					}
				}
				//}
			}
		}
	}
}
gmds::Node
Tools::getOpposedNodeOnEdge(const gmds::Node &ANode, const gmds::Edge &AEdge)
{
	const auto &nodes = AEdge.get<Node>();

	return nodes.at(0).id() != ANode.id() ? nodes.at(0) : nodes.at(1);
}

gmds::Face
Tools::getNextFace(const gmds::TCellID &ANodeID, const gmds::TCellID &AFaceID, gmds::Edge &AEdge)
{
	// find opposed Face
	const auto &faces = AEdge.get<Face>();
	gmds::Face opposedFace = faces[0].id() == AFaceID ? faces[1] : faces[0];

	// set next opposed edge
	for (const auto &ei : opposedFace.get<Edge>()) {
		if (ei.id() != AEdge.id()) {
			const auto &nodeIDs = ei.getIDs<Node>();
			if (nodeIDs[0] == ANodeID || nodeIDs[1] == ANodeID) {
				AEdge = ei;
				break;
			}
		}
	}
	return opposedFace;
}

gmds::Node
Tools::getCommonNode(const gmds::Edge &AEdge1, const gmds::Edge &AEdge2)
{
	for (const auto &ni : AEdge1.get<Node>()) {
		for (const auto &nj : AEdge2.get<Node>()) {
			if (ni == nj) {
				return ni;
			}
		}
	}
	throw GMDSException("No common Node found");
}
gmds::Node
Tools::getCommonNode(const gmds::Face &AFace1, const gmds::Face &AFace2)
{
	for (const auto &ni : AFace1.get<Node>()) {
		for (const auto &nj : AFace2.get<Node>()) {
			if (ni == nj) {
				return ni;
			}
		}
	}
	throw GMDSException("No common Node found");
}

gmds::Edge
Tools::getCommonEdge(const gmds::Face &AFace1, const gmds::Face &AFace2)
{
	for (const auto &ei : AFace1.get<Edge>()) {
		for (const auto &ej : AFace2.get<Edge>()) {
			if (ei.id() == ej.id()) {
				return ei;
			}
		}
	}
	throw GMDSException("No common Edge found");
}

bool
Tools::isAdjacency(const gmds::TCellID &AFace1, const gmds::TCellID &AFace2)
{
	return isAdjacency(m_mesh->get<Face>(AFace1), m_mesh->get<Face>(AFace2));
}

bool
Tools::isAdjacency(const gmds::Face &AFace1, const gmds::Face &AFace2)
{
	for (const auto &ei : AFace1.get<Edge>()) {
		for (const auto &ej : AFace2.get<Edge>()) {
			if (ei.id() == ej.id()) {
				return true;
			}
		}
	}
	return false;
}
