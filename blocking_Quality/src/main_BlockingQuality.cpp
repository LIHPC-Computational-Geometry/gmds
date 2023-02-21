//
// Created by bourmaudp on 02/12/22.
//
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <string>


#include <gmds/ig/Mesh.h>
#include <gmds/ig/Node.h>
#include <gmds/quality/QuadQuality.h>
#include <gmds/quality/HexQuality.h>
#include <gmds/blocking_Quality/BlockingQuality.h>
#include <gmds/igalgo/GridBuilder.h>

#include <gmds/blockMesher/BlockMesher.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/cad/GeomPoint.h>
#include <gmds/cad/GeomCurve.h>
#include <gmds/cad/GeomSurface.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/smoothy/LaplacianSmoother.h>
#include <gmds/blocking_Quality/LinkerBlockingGeom.h>


#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gmds/ig/MeshDoctor.h>

#include <gmds/blocking_Quality/mainBis1.h>


//================================================================================
using namespace gmds;
int main(int argc, char* argv[]){

// LES FORMES

/*=========================================================*/
// CUBE SIMPLE POUR LE RUBIKSCUBE DE TAILLE 2x2x2
Mesh cubeRubiksCube(MeshModel(DIM3 | R | F | E | N |
	                           R2N | R2F | R2E |
	                           F2N | F2R | F2E |
	                           E2F | E2N | E2R| N2E | N2R | N2F));


Node crn0 = cubeRubiksCube.newNode(math::Point(0,0,0));
Node crn1 = cubeRubiksCube.newNode(math::Point(4,0,0));
Node crn2 = cubeRubiksCube.newNode(math::Point(4,0,4));
Node crn3 = cubeRubiksCube.newNode(math::Point(0,0,4));
Node crn4 = cubeRubiksCube.newNode(math::Point(0,4,0));
Node crn5 = cubeRubiksCube.newNode(math::Point(4,4,0));
Node crn6 = cubeRubiksCube.newNode(math::Point(4,4,4));
Node crn7 = cubeRubiksCube.newNode(math::Point(0,4,4));


Region crr1 = cubeRubiksCube.newHex(crn3,crn2,crn1,crn0,crn7,crn6,crn5,crn4);

MeshDoctor cubeDoctor(&cubeRubiksCube);
cubeDoctor.buildFacesAndR2F();
cubeDoctor.buildEdgesAndX2E();
cubeDoctor.updateUpwardConnectivity();




Variable<int> *boundaryEdgeCubeRubiksCube = cubeRubiksCube.newVariable<int,gmds::GMDS_EDGE>("boundaryEdge");

for (auto e : cubeRubiksCube.edges()){
	//std::cout<<"in edge : "<<e<<std::endl;
	//std::cout<<"Nb regions : "<<rubiksCube.get<Edge>(e).nbRegions()<<std::endl;
	if(cubeRubiksCube.get<Edge>(e).nbRegions()== 1){
		//std::cout<<"dans if "<<std::endl;
		boundaryEdgeCubeRubiksCube->set(e,1);
	}
	else{
		//std::cout<<"Dans else"<<std::endl;
		boundaryEdgeCubeRubiksCube->set(e,0);
	}
}


/*=========================================================*/
// FORME GRILLE 3D

Mesh rubiksCube(MeshModel(DIM3 | R | F | E | N |
	                       R2N | R2F | R2E | F2N | F2R | F2E | E2F | E2N |E2R| N2E | N2R | N2F));


GridBuilder gridBuilder(&rubiksCube,3);
double step = 2;
int nb = 3;
int AXnb = nb;
int AYnb = nb;
int AZnb = nb;

gridBuilder.execute(AXnb,step,AYnb,step,AZnb,step);
MeshDoctor doctor(&rubiksCube);
doctor.buildFacesAndR2F();
doctor.buildEdgesAndX2E();
doctor.updateUpwardConnectivity();


Variable<int> *boundaryEdgeRubiksCube = rubiksCube.newVariable<int,gmds::GMDS_EDGE>("boundaryEdge");

	for(auto r : rubiksCube.regions()){
	   Region region = rubiksCube.get<Region>(r);
	   std::cout<<"region "<< r <<std::endl;
	   std::vector<TCellID> edges;
	   region.delegateGetEdgeIDs(edges);
	   for(auto e : edges){
	      rubiksCube.get<Edge>(e).add(region);
	   }
	}

for (auto e : rubiksCube.edges()){
	//std::cout<<"in edge : "<<e<<std::endl;
	//std::cout<<"Nb regions : "<<rubiksCube.get<Edge>(e).nbRegions()<<std::endl;
	if(rubiksCube.get<Edge>(e).nbRegions()== 1){
		//std::cout<<"dans if "<<std::endl;
		boundaryEdgeRubiksCube->set(e,1);
	}
	else{
		//std::cout<<"Dans else"<<std::endl;
		boundaryEdgeRubiksCube->set(e,0);
	}
}

//=================================================================================================
//=================================FORME BLOCKING NON VALIDE ======================================
//=================================================================================================
/*
Mesh badRubiksCube(MeshModel(DIM3 | R | F | E | N |
	                       R2N | R2F | R2E | F2N | F2R | F2E | E2F | E2N |E2R| N2E | N2R | N2F));


Node brcn0 = badRubiksCube.newNode(math::Point(0,0,0));
Node brcn1 = badRubiksCube.newNode(math::Point(4,0,0));
Node brcn2 = badRubiksCube.newNode(math::Point(2,2,0));
Node brcn3 = badRubiksCube.newNode(math::Point(0,4,0));
Node brcn4 = badRubiksCube.newNode(math::Point(0,0,4));
Node brcn5 = badRubiksCube.newNode(math::Point(4,0,4));
Node brcn6 = badRubiksCube.newNode(math::Point(2,2,4));
Node brcn7 = badRubiksCube.newNode(math::Point(0,4,4));


Region brcr0 = badRubiksCube.newHex(crn3,crn2,crn1,crn0,crn7,crn6,crn5,crn4);


MeshDoctor badDoc(&badRubiksCube);
badDoc.buildFacesAndR2F();
badDoc.buildEdgesAndX2E();
badDoc.updateUpwardConnectivity();

*/
//=================================================================================================
//=================================================================================================
//=================================================================================================
//=================================================================================================
//=================================================================================================
//=================================================================================================


cad::FACManager geom_manager;
geom_manager.initFrom3DMesh(&cubeRubiksCube);

cad::GeomMeshLinker linker(&rubiksCube,&geom_manager);
std::cout<<"Geom info ("<<geom_manager.getNbVolumes()<<", "
	       <<geom_manager.getNbSurfaces()<<", "
	       <<geom_manager.getNbCurves()<<", "
	       <<geom_manager.getNbPoints()<<")"<<std::endl;
//=================================================================================================
//=================================================================================================


std::vector<cad::GeomSurface*> surfaces;
std::vector<cad::GeomCurve*> curves;
std::vector<cad::GeomPoint*> points;

geom_manager.getSurfaces(surfaces);
geom_manager.getCurves(curves);
geom_manager.getPoints(points);


LinkerBlockingGeom validBlocking(&rubiksCube,&geom_manager);
validBlocking.execute();


}


/*
 *
 *
 *
cad::FACManager geom_manager;
geom_manager.initFrom3DMesh(&cubeRubiksCube);

cad::GeomMeshLinker linker(&rubiksCube,&geom_manager);
std::cout<<"Geom info ("<<geom_manager.getNbVolumes()<<", "
          <<geom_manager.getNbSurfaces()<<", "
          <<geom_manager.getNbCurves()<<", "
          <<geom_manager.getNbPoints()<<")"<<std::endl;
//=================================================================================================
//=================================================================================================


std::vector<cad::GeomSurface*> surfaces;
std::vector<cad::GeomCurve*> curves;
std::vector<cad::GeomPoint*> points;

geom_manager.getSurfaces(surfaces);
geom_manager.getCurves(curves);
geom_manager.getPoints(points);

//==================================================================
//First, we classify each face

//==================================================================
// MARK ALL THE BOUNDARY CELL OF THE INIT MESH
//==================================================================

std::cout<<"> Start block boundary retrieval"<<std::endl;
//we get all the nodes that are on the mesh boundary
BoundaryOperator op(&rubiksCube);
auto mark_node_NAN = rubiksCube.newMark<Node>();
auto mark_node_on_pnt = rubiksCube.newMark<Node>();
auto mark_node_on_crv = rubiksCube.newMark<Node>();
auto mark_node_on_srf = rubiksCube.newMark<Node>();
auto mark_edge_on_crv = rubiksCube.newMark<Edge>();
auto mark_edge_on_srf = rubiksCube.newMark<Edge>();
auto mark_face_on_srf = rubiksCube.newMark<Face>();

op.markCellOnGeometry(mark_face_on_srf,
                      mark_edge_on_srf,
                      mark_node_on_srf,
                      mark_edge_on_crv,
                      mark_node_on_crv,
                      mark_node_on_pnt,
                      mark_node_NAN);

for(auto f_id: rubiksCube.faces()){
Face f= rubiksCube.get<Face>(f_id);
if(rubiksCube.isMarked(f, mark_face_on_srf) ){
	//we've got a boundary face
	math::Point p = f.center();
	math::Vector3d v = f.normal();
	//we ensure to get the outward normal that goes from the inside block
	//towards the outside
	//as f is a boundary face, it is connected to only one block
	Region r = f.get<Region>()[0];
	math::Vector3d v2inside= r.center()-p;
	if(v.dot(v2inside)>0){
		//v points inside
		v = v.opp();
	}

	v.normalize();
	std::cout<<"FACE TO PROJECT "<<f_id<<" with center point  "<<p<<" and normal direction "<<v<<std::endl;
	double min_dist = 100000;
	int min_entity_dim=-1;
	int min_entity_id = -1;
	std::map<TCellID , double> surf_dist;
	std::map<TCellID , double> surf_dot;
	bool found_surf=false;
	for(auto s:surfaces){
		math::Point closest_pnt=p;
		s->project(closest_pnt);
		double dist = p.distance2(closest_pnt);
		math::Vector v_proj=closest_pnt-p;
		v_proj.normalize();
		surf_dot[s->id()]=(v.dot(v_proj));
		surf_dist[s->id()]=dist;

		if(dist<min_dist) {
			//either the points are the same or we check the vector orientation
			if (dist<1e-4 || fabs(v.dot(v_proj)) > 0.2 ) {
				bool keep_projection = true;

				if (v.dot(v_proj)<0){
					//we have to check that we do not cross the volume. It can occur for thin curved
					//shape (like cylinder with a hole)

					//We take other points closed to the face corners and we project all of them, if they go to
					// the same surf, we keep it. Otherwise, we will put another surface later.
					keep_projection=false;
					std::vector<Node> corners = f.get<Node>();
					std::vector<math::Point> corner_pnts;
					for(auto c:corners){
						corner_pnts.push_back(0.1*p+0.9* c.point());
					}
					std::set<int> corner_surf;
					for(auto cp:corner_pnts) {
						double cp_min_dist = 100000;
						auto cp_min_surf_id = -1;
						for (auto s: surfaces) {
							math::Point closest_pnt = cp;
							s->project(closest_pnt);

							double local_dist = cp.distance2(closest_pnt);
							if (local_dist < cp_min_dist) {
								cp_min_dist = local_dist;
								cp_min_surf_id = s->id();
							}
						}
						corner_surf.insert(cp_min_surf_id);
					}
					if(corner_surf.size()==1 && (*corner_surf.begin()==s->id())){
						//ok we project on it
						keep_projection= true;
					}

				}
				if(keep_projection) {
					min_dist = dist;
					min_entity_dim = 2;
					min_entity_id = s->id();
					found_surf = true;
				}
			}
		}
	}
	if(!found_surf){
		min_dist = 100000;
		min_entity_dim=-1;
		min_entity_id = -1;
		//means we have a block face in a concave area
		for(auto s:surfaces){
			math::Point closest_pnt=p;
			s->project(closest_pnt);

			double dist = p.distance2(closest_pnt);
			math::Vector v_proj=closest_pnt-p;
			v_proj.normalize();
			surf_dot[s->id()]=(v.dot(v_proj));
			surf_dist[s->id()]=dist;
			if(dist<min_dist) {
				min_dist = dist;
				min_entity_dim = 2;
				min_entity_id = s->id();
				found_surf=true;

			}
		}
	}
	if (min_entity_dim==2){
		std::cout<<"==> Link face "<<f_id<<" on surf "<<min_entity_id<<" with distance:" <<min_dist<<std::endl;
		linker.linkFaceToSurface(f_id,min_entity_id);
	}
	else{
		std::cout<<"\t ====> Link error for classifying a face"<<std::endl;
		throw GMDSException("Link error for classifying a face");
	}

}
}
//==================================================================
//Second, we classify each edge
for(auto e_id: rubiksCube.edges()){
Edge e= rubiksCube.get<Edge>(e_id);
if(rubiksCube.isMarked(e, mark_edge_on_crv) ||
	 rubiksCube.isMarked(e, mark_edge_on_srf) ){
	//we've got a boundary edge, now we get the 2 boundary faces
	// around e
	std::vector<Node> e_nodes = e.get<Node>();
	std::vector<TCellID> adj_face_id = e.getIDs<Face>();
	std::vector<TCellID> adj_bnd_face_ids;
	for(auto f_id:adj_face_id){
		if(rubiksCube.isMarked<Face>(f_id, mark_face_on_srf)){
			adj_bnd_face_ids.push_back(f_id);
		}
	}
	if(adj_bnd_face_ids.size()!=2){
		std::cout<<"ERROR: One boundary edge is adjacent to more "
			       <<"than 2 boundary faces ("
			       <<"face "<<e_id<<")"<<std::endl;
	}
	auto surf0_id = linker.getGeomId<Face>(adj_bnd_face_ids[0]);
	auto surf1_id = linker.getGeomId<Face>(adj_bnd_face_ids[1]);

	if(surf0_id==surf1_id){
		//edge embedded in the surface
		linker.linkEdgeToSurface(e_id,surf1_id);
	}
	else{
		//it must be an edge classified on curves
		// Warning, like for the face, it is not necessary the geometrically closest to be connected to.
		//We check which curve is bounded by surf0_id and surf1_id
		bool found = false;

		std::vector<cad::GeomCurve*> candidate_curves;
		for(auto ci:curves) {
			std::vector<cad::GeomSurface *> c_surfs = ci->surfaces();
			if (c_surfs.size() == 2) {
				if ((c_surfs[0]->id() == surf0_id && c_surfs[1]->id() == surf1_id) ||
					 (c_surfs[1]->id() == surf0_id && c_surfs[0]->id() == surf1_id)) {
					std::cout<<"Curve "<<ci->id()<<" adj to surfaces "<<c_surfs[0]->id()<<" and "<<c_surfs[1]->id()<<std::endl;
					found = true;
					candidate_curves.push_back(ci);
				}
			}
		}
		if(found) {
			if(candidate_curves.size()==1) {
				linker.linkEdgeToCurve(e_id, candidate_curves[0]->id());
			}
			else{
				//we project on each of curve and keep the closest one
				double min_dist = 1e6;
				cad::GeomCurve* selected_curve = NULL;
				math::Point center_edge = e.center();
				for(auto ci:candidate_curves) {
					math::Point pi = center_edge;
					ci->project(pi);
					auto di = pi.distance2(center_edge);
					if(di<min_dist){
						min_dist=di;
						selected_curve = ci;
					}
				}

				linker.linkEdgeToCurve(e_id, selected_curve->id());

			}
		}
		else{
			throw GMDSException("BlockMesher error: impossible to link a block edge onto the geometry");
		}
	}
}
}
//==================================================================
//we classify each node
auto on_pnt=0, on_curve=0, on_surf=0;
for(auto n_id: rubiksCube.nodes()){
Node n= rubiksCube.get<Node>(n_id);
if(rubiksCube.isMarked(n, mark_node_on_pnt) ||
	 rubiksCube.isMarked(n, mark_node_on_crv) ||
	 rubiksCube.isMarked(n, mark_node_on_srf) ){
	//we've got a boundary node
	//As face and edges are already classified, we used it
	// A node is on a surface if all its adjacent boundary faces are on the same
	// surface. If they are on two it is on a curve potentially
	std::vector<TCellID> adj_faces = n.getIDs<Face>();
	std::vector<TCellID> adj_bnd_faces;
	for(auto id_face:adj_faces){
		if(rubiksCube.isMarked<Face>(id_face, mark_face_on_srf)) {
			adj_bnd_faces.push_back(id_face);
		}
	}
	std::set<int> bnd_surfaces;

	for (auto  id_face:adj_bnd_faces){
		bnd_surfaces.insert(linker.getGeomId<Face>(id_face));
	}
	int min_entity_dim=-1;
	int min_entity_id = -1;
	double min_dist = 100000;

	if(bnd_surfaces.size()==1){
		//on a single surface
		min_entity_dim=2;
		min_entity_id=*(bnd_surfaces.begin());
	}
	else if(bnd_surfaces.size()==2){
		//on a curve or a point that is connected to a single curve
		math::Point node_loc = n.point();
		for(auto p:points){
			math::Point closest_pnt=node_loc;
			double dist = node_loc.distance(p->point());
			if(dist<min_dist){
				min_dist =dist;
				min_entity_dim=0;
				min_entity_id=p->id();
			}
		}
		for(auto c:curves){
			math::Point closest_pnt=node_loc;
			c->project(closest_pnt);
			double dist = node_loc.distance(closest_pnt);
			if(dist<min_dist){
				//WARNING: Take care of this trick that is not good at all but mandatory to be staying on the
				// curve and not on the surface
				if(dist<1e-4)
					dist=0;
				min_dist =dist;
				min_entity_dim=1;
				min_entity_id=c->id();
			}
		}

	} else{
		//On a point!!
		math::Point node_loc = n.point();
		for(auto p:points){
			math::Point closest_pnt=node_loc;
			double dist = node_loc.distance(p->point());
			if(dist<min_dist){
				min_dist =dist;
				min_entity_dim=0;
				min_entity_id=p->id();
			}
		}
	}

	if(min_entity_dim==0){
		on_pnt++;
		linker.linkNodeToPoint(n_id,min_entity_id);
		std::cout<<"Node "<<n_id<<" is on point "<<min_entity_id<<std::endl;
	}
	else if (min_entity_dim==1){
		on_curve++;
		linker.linkNodeToCurve(n_id,min_entity_id);
	}
	else if (min_entity_dim==2){
		on_surf++;
		linker.linkNodeToSurface(n_id,min_entity_id);
	}
	else{
		throw GMDSException("Link error for classifying a node");
	}

}
}
std::cout<<"  info [node classified on points "<<on_pnt;
std::cout<<", on curves "<<on_curve;
std::cout<<", on surfs "<<on_surf<<"]"<<std::endl;

linker.writeVTKDebugMesh("linker_debug_1.vtk");

//=================================================================================================
//=================================================================================================
BlockMesher bm(&rubiksCube,&linker);
int edge_discretization = 4;

auto val = bm.execute(edge_discretization);

std::map<TCellID, Cell::Data> mesh_node_info = bm.mesh_node_classification();
std::map<TCellID, Cell::Data> mesh_edge_info = bm.mesh_edge_classification();
std::map<TCellID, Cell::Data> mesh_face_info = bm.mesh_face_classification();

std::cout<<"> Start block smoothing"<<std::endl;
//bm.mesh()->changeModel(MeshModel(DIM3|R|N|E|F|R2N|F2N|E2N|N2F|N2E|N2R));
MeshDoctor doc3(bm.mesh());
doc3.updateUpwardConnectivity();
cad::GeomMeshLinker linker_final(bm.mesh(),&geom_manager);
for(auto info:mesh_node_info){
Node n = bm.mesh()->get<Node>(info.first);
if(info.second.dim==0) {
	linker_final.linkNodeToPoint(n.id(), info.second.id);
}
else if(info.second.dim==1) {
	linker_final.linkNodeToCurve(n.id(), info.second.id);
}
else if(info.second.dim==2)
	linker_final.linkNodeToSurface(n.id(),info.second.id);
}
for(auto info:mesh_edge_info){
std::cout<<"Dans le for"<<std::endl;
Edge e = bm.mesh()->get<Edge>(info.first);
if(info.second.dim==1) {
	linker_final.linkEdgeToCurve(e.id(), info.second.id);
	std::cout<<"Edge "<<e.getIDs<Node>()[0]<<"-"<<e.getIDs<Node>()[1]<<" on curve "<<info.second.id<<std::endl;
}
}
for(auto info:mesh_face_info){
Face f = bm.mesh()->get<Face>(info.first);
if(info.second.dim==2)
	linker_final.linkFaceToSurface(f.id(),info.second.id);
}

 */



//======================= BOU de CODE obsolete =======================================


/*
//==================================================================
//First, we classify each face

//==================================================================
// MARK ALL THE BOUNDARY CELL OF THE INIT MESH
//==================================================================

std::cout<<"> Start block boundary retrieval"<<std::endl;
//we get all the nodes that are on the mesh boundary
BoundaryOperator op(&badRubiksCube);
auto mark_node_NAN = badRubiksCube.newMark<Node>();
auto mark_node_on_pnt = badRubiksCube.newMark<Node>();
auto mark_node_on_crv = badRubiksCube.newMark<Node>();
auto mark_node_on_srf = badRubiksCube.newMark<Node>();
auto mark_edge_on_crv = badRubiksCube.newMark<Edge>();
auto mark_edge_on_srf = badRubiksCube.newMark<Edge>();
auto mark_face_on_srf = badRubiksCube.newMark<Face>();

op.markCellOnGeometry(mark_face_on_srf,
                      mark_edge_on_srf,
                      mark_node_on_srf,
                      mark_edge_on_crv,
                      mark_node_on_crv,
                      mark_node_on_pnt,
                      mark_node_NAN);

for(auto f_id: badRubiksCube.faces()){
Face f= badRubiksCube.get<Face>(f_id);
if(badRubiksCube.isMarked(f, mark_face_on_srf) ){
	//we've got a boundary face
	math::Point p = f.center();
	math::Vector3d v = f.normal();
	//we ensure to get the outward normal that goes from the inside block
	//towards the outside
	//as f is a boundary face, it is connected to only one block
	Region r = f.get<Region>()[0];
	math::Vector3d v2inside= r.center()-p;
	if(v.dot(v2inside)>0){
		//v points inside
		v = v.opp();
	}

	v.normalize();
	std::cout<<"FACE TO PROJECT "<<f_id<<" with center point  "<<p<<" and normal direction "<<v<<std::endl;
	double min_dist = 100000;
	int min_entity_dim=-1;
	int min_entity_id = -1;
	std::map<TCellID , double> surf_dist;
	std::map<TCellID , double> surf_dot;
	bool found_surf=false;
	for(auto s:surfaces){
		math::Point closest_pnt=p;
		s->project(closest_pnt);
		double dist = p.distance2(closest_pnt);
		math::Vector v_proj=closest_pnt-p;
		v_proj.normalize();
		surf_dot[s->id()]=(v.dot(v_proj));
		surf_dist[s->id()]=dist;

		if(dist<min_dist) {
			//either the points are the same or we check the vector orientation
			if (dist<1e-4 || fabs(v.dot(v_proj)) > 0.2 ) {
				bool keep_projection = true;

				if (v.dot(v_proj)<0){
					//we have to check that we do not cross the volume. It can occur for thin curved
					//shape (like cylinder with a hole)

					//We take other points closed to the face corners and we project all of them, if they go to
					// the same surf, we keep it. Otherwise, we will put another surface later.
					keep_projection=false;
					std::vector<Node> corners = f.get<Node>();
					std::vector<math::Point> corner_pnts;
					for(auto c:corners){
						corner_pnts.push_back(0.1*p+0.9* c.point());
					}
					std::set<int> corner_surf;
					for(auto cp:corner_pnts) {
						double cp_min_dist = 100000;
						auto cp_min_surf_id = -1;
						for (auto s: surfaces) {
							math::Point closest_pnt = cp;
							s->project(closest_pnt);

							double local_dist = cp.distance2(closest_pnt);
							if (local_dist < cp_min_dist) {
								cp_min_dist = local_dist;
								cp_min_surf_id = s->id();
							}
						}
						corner_surf.insert(cp_min_surf_id);
					}
					if(corner_surf.size()==1 && (*corner_surf.begin()==s->id())){
						//ok we project on it
						keep_projection= true;
					}

				}
				if(keep_projection) {
					min_dist = dist;
					min_entity_dim = 2;
					min_entity_id = s->id();
					found_surf = true;
				}
			}
		}
	}
	if(!found_surf){
		min_dist = 100000;
		min_entity_dim=-1;
		min_entity_id = -1;
		//means we have a block face in a concave area
		for(auto s:surfaces){
			math::Point closest_pnt=p;
			s->project(closest_pnt);

			double dist = p.distance2(closest_pnt);
			math::Vector v_proj=closest_pnt-p;
			v_proj.normalize();
			surf_dot[s->id()]=(v.dot(v_proj));
			surf_dist[s->id()]=dist;
			if(dist<min_dist) {
				min_dist = dist;
				min_entity_dim = 2;
				min_entity_id = s->id();
				found_surf=true;

			}
		}
	}
	if (min_entity_dim==2){
		std::cout<<"==> Link face "<<f_id<<" on surf "<<min_entity_id<<" with distance:" <<min_dist<<std::endl;
		linker.linkFaceToSurface(f_id,min_entity_id);
	}
	else{
		std::cout<<"\t ====> Link error for classifying a face"<<std::endl;
		throw GMDSException("Link error for classifying a face");
	}

}
}
//==================================================================
//Second, we classify each edge
for(auto e_id: badRubiksCube.edges()){
Edge e= badRubiksCube.get<Edge>(e_id);
if(badRubiksCube.isMarked(e, mark_edge_on_crv) ||
	 badRubiksCube.isMarked(e, mark_edge_on_srf) ){
	//we've got a boundary edge, now we get the 2 boundary faces
	// around e
	std::vector<Node> e_nodes = e.get<Node>();
	std::vector<TCellID> adj_face_id = e.getIDs<Face>();
	std::vector<TCellID> adj_bnd_face_ids;
	for(auto f_id:adj_face_id){
		if(badRubiksCube.isMarked<Face>(f_id, mark_face_on_srf)){
			adj_bnd_face_ids.push_back(f_id);
		}
	}
	if(adj_bnd_face_ids.size()!=2){
		std::cout<<"ERROR: One boundary edge is adjacent to more "
			       <<"than 2 boundary faces ("
			       <<"face "<<e_id<<")"<<std::endl;
	}
	auto surf0_id = linker.getGeomId<Face>(adj_bnd_face_ids[0]);
	auto surf1_id = linker.getGeomId<Face>(adj_bnd_face_ids[1]);

	if(surf0_id==surf1_id){
		//edge embedded in the surface
		linker.linkEdgeToSurface(e_id,surf1_id);
	}
	else{
		//it must be an edge classified on curves
		// Warning, like for the face, it is not necessary the geometrically closest to be connected to.
		//We check which curve is bounded by surf0_id and surf1_id
		bool found = false;

		std::vector<cad::GeomCurve*> candidate_curves;
		for(auto ci:curves) {
			std::vector<cad::GeomSurface *> c_surfs = ci->surfaces();
			if (c_surfs.size() == 2) {
				if ((c_surfs[0]->id() == surf0_id && c_surfs[1]->id() == surf1_id) ||
					 (c_surfs[1]->id() == surf0_id && c_surfs[0]->id() == surf1_id)) {
					std::cout<<"Curve "<<ci->id()<<" adj to surfaces "<<c_surfs[0]->id()<<" and "<<c_surfs[1]->id()<<std::endl;
					found = true;
					candidate_curves.push_back(ci);
				}
			}
		}
		if(found) {
			if(candidate_curves.size()==1) {
				linker.linkEdgeToCurve(e_id, candidate_curves[0]->id());
			}
			else{
				//we project on each of curve and keep the closest one
				double min_dist = 1e6;
				cad::GeomCurve* selected_curve = NULL;
				math::Point center_edge = e.center();
				for(auto ci:candidate_curves) {
					math::Point pi = center_edge;
					ci->project(pi);
					auto di = pi.distance2(center_edge);
					if(di<min_dist){
						min_dist=di;
						selected_curve = ci;
					}
				}

				linker.linkEdgeToCurve(e_id, selected_curve->id());

			}
		}
		else{
			throw GMDSException("BlockMesher error: impossible to link a block edge onto the geometry");
		}
	}
}
}
//==================================================================
//we classify each node
auto on_pnt=0, on_curve=0, on_surf=0;
for(auto n_id: badRubiksCube.nodes()){
Node n= badRubiksCube.get<Node>(n_id);
if(badRubiksCube.isMarked(n, mark_node_on_pnt) ||
	 badRubiksCube.isMarked(n, mark_node_on_crv) ||
	 badRubiksCube.isMarked(n, mark_node_on_srf) ){
	//we've got a boundary node
	//As face and edges are already classified, we used it
	// A node is on a surface if all its adjacent boundary faces are on the same
	// surface. If they are on two it is on a curve potentially
	std::vector<TCellID> adj_faces = n.getIDs<Face>();
	std::vector<TCellID> adj_bnd_faces;
	for(auto id_face:adj_faces){
		if(badRubiksCube.isMarked<Face>(id_face, mark_face_on_srf)) {
			adj_bnd_faces.push_back(id_face);
		}
	}
	std::set<int> bnd_surfaces;

	for (auto  id_face:adj_bnd_faces){
		bnd_surfaces.insert(linker.getGeomId<Face>(id_face));
	}
	int min_entity_dim=-1;
	int min_entity_id = -1;
	double min_dist = 100000;

	if(bnd_surfaces.size()==1){
		//on a single surface
		min_entity_dim=2;
		min_entity_id=*(bnd_surfaces.begin());
	}
	else if(bnd_surfaces.size()==2){
		//on a curve or a point that is connected to a single curve
		math::Point node_loc = n.point();
		for(auto p:points){
			math::Point closest_pnt=node_loc;
			double dist = node_loc.distance(p->point());
			if(dist<min_dist){
				min_dist =dist;
				min_entity_dim=0;
				min_entity_id=p->id();
			}
		}
		for(auto c:curves){
			math::Point closest_pnt=node_loc;
			c->project(closest_pnt);
			double dist = node_loc.distance(closest_pnt);
			if(dist<min_dist){
				//WARNING: Take care of this trick that is not good at all but mandatory to be staying on the
				// curve and not on the surface
				if(dist<1e-4)
					dist=0;
				min_dist =dist;
				min_entity_dim=1;
				min_entity_id=c->id();
			}
		}

	} else{
		//On a point!!
		math::Point node_loc = n.point();
		for(auto p:points){
			math::Point closest_pnt=node_loc;
			double dist = node_loc.distance(p->point());
			if(dist<min_dist){
				min_dist =dist;
				min_entity_dim=0;
				min_entity_id=p->id();
			}
		}
	}

	if(min_entity_dim==0){
		on_pnt++;
		linker.linkNodeToPoint(n_id,min_entity_id);
		std::cout<<"Node "<<n_id<<" is on point "<<min_entity_id<<std::endl;
	}
	else if (min_entity_dim==1){
		on_curve++;
		linker.linkNodeToCurve(n_id,min_entity_id);
	}
	else if (min_entity_dim==2){
		on_surf++;
		linker.linkNodeToSurface(n_id,min_entity_id);
	}
	else{
		throw GMDSException("Link error for classifying a node");
	}

}
}
std::cout<<"  info [node classified on points "<<on_pnt;
std::cout<<", on curves "<<on_curve;
std::cout<<", on surfs "<<on_surf<<"]"<<std::endl;

linker.writeVTKDebugMesh("linker_debug.vtk");

//=================================================================================================
//=================================================================================================
BlockMesher bm(&badRubiksCube,&linker);
int edge_discretization = 1;

auto val = bm.execute(edge_discretization);

std::map<TCellID, Cell::Data> mesh_node_info = bm.mesh_node_classification();
std::map<TCellID, Cell::Data> mesh_edge_info = bm.mesh_edge_classification();
std::map<TCellID, Cell::Data> mesh_face_info = bm.mesh_face_classification();

std::cout<<"> Start block smoothing"<<std::endl;
//bm.mesh()->changeModel(MeshModel(DIM3|R|N|E|F|R2N|F2N|E2N|N2F|N2E|N2R));
MeshDoctor doc3(bm.mesh());
doc3.updateUpwardConnectivity();
cad::GeomMeshLinker linker_final(bm.mesh(),&geom_manager);
for(auto info:mesh_node_info){
Node n = bm.mesh()->get<Node>(info.first);
if(info.second.dim==0) {
	linker_final.linkNodeToPoint(n.id(), info.second.id);
}
else if(info.second.dim==1) {
	linker_final.linkNodeToCurve(n.id(), info.second.id);
}
else if(info.second.dim==2)
	linker_final.linkNodeToSurface(n.id(),info.second.id);
}
for(auto info:mesh_edge_info){
Edge e = bm.mesh()->get<Edge>(info.first);
if(info.second.dim==1) {
	linker_final.linkEdgeToCurve(e.id(), info.second.id);
	std::cout<<"Edge "<<e.getIDs<Node>()[0]<<"-"<<e.getIDs<Node>()[1]<<" on curve "<<info.second.id<<std::endl;
}
}
for(auto info:mesh_face_info){
Face f = bm.mesh()->get<Face>(info.first);
if(info.second.dim==2)
	linker_final.linkFaceToSurface(f.id(),info.second.id);
}

std::cout<<"> Start final mesh writing"<<std::endl;
IGMeshIOService ioService_mesh(bm.mesh());
VTKWriter vtkWriter_mesh(&ioService_mesh);
vtkWriter_mesh.setCellOptions(N|F);
vtkWriter_mesh.setDataOptions(N|F);
vtkWriter_mesh.write("mBlockMesher.vtk");

VTKWriter vtkWriter_curves(&ioService_mesh);
vtkWriter_curves.setCellOptions(N|E);
vtkWriter_curves.setDataOptions(N|E);
vtkWriter_curves.write("mesh_curves1.vtk");

//SAVE BLOCKING VTK
/*
IGMeshIOService ioService_blocking(&rubiksCube);
VTKWriter vtkWriter_blocking(&ioService_blocking);
vtkWriter_blocking.setCellOptions(N|F);
vtkWriter_blocking.setDataOptions(N|F);
vtkWriter_blocking.write("rubiksCube.vtk");
*/
//SAVE GEOMETRY VTK
   /*
IGMeshIOService ioService_geom(&cubeRubiksCube);
VTKWriter vtkWriter_geom(&ioService_geom);
vtkWriter_geom.setCellOptions(N|F);
vtkWriter_geom.setDataOptions(N|F);
vtkWriter_geom.write("cubeRubiksCube.vtk");



//SAVE BAD RUBIKSCUBE
IGMeshIOService ioService_bad(&badRubiksCube);
VTKWriter vtkWriter_bad(&ioService_bad);
vtkWriter_bad.setCellOptions(N|F);
vtkWriter_bad.setDataOptions(N|F);
vtkWriter_bad.write("badRubiksCube.vtk");

//=================================================================================================
//======================================= END LINK ================================================
//=================================================================================================

std::cout << "======== Task done by blocker =========" << std::endl;

int idChecked = 0;
std::cout<<"info du linker Node 0, info 1 : "<<linker_final.getGeomInfo(bm.mesh()->get<Node>(idChecked)).first
          <<" et info 2 : "<<
   linker_final.getGeomInfo(bm.mesh()->get<Node>(idChecked)).second<<std::endl;

std::cout<<"Check info Arete : "<<std::endl;
Edge theEdgeSelected = bm.mesh()->get<Edge>(idChecked);
std::cout<<"Noeud 1 : "<<theEdgeSelected.get<Node>()[0]<<std::endl;
std::cout<<"Noeud 2 : "<<theEdgeSelected.get<Node>()[1]<<std::endl;
std::cout<<"========================================================="<<std::endl;
std::cout<<"Check info Curves : "<<std::endl;
Edge theEdgeSelectedCurve = cubeRubiksCube.get<Edge>(11);
std::cout<<"Noeud 1 : "<<theEdgeSelectedCurve.get<Node>()[0]<<std::endl;
std::cout<<"Noeud 2 : "<<theEdgeSelectedCurve.get<Node>()[1]<<std::endl;
/*
std::cout<<"Lien Nodes2Points valide ? " <<linkNodes2Points(&badRubiksCube,&cubeRubiksCube,&linker_final)<<std::endl;
std::cout<<"Lien Edges2Curves valide ? " <<linkEdges2Curves(&badRubiksCube,&cubeRubiksCube,&linker_final)<<std::endl;
std::cout<<"Lien Faces2Surfaces valide ? " <<linkFaces2Surfaces(&badRubiksCube,&cubeRubiksCube,&linker_final)<<std::endl;
*/
/*
noLinkNodes2Points(&badRubiksCube,&cubeRubiksCube,&linker_final);

std::cout<<"Element dans nodes : "<<cubeRubiksCube.nodes().getNbElements()	<<std::endl;

//=================================================================================================
//=================================================================================================
//=================================================================================================
//=================================================================================================
//=================================================================================================
//=================================================================================================
*/