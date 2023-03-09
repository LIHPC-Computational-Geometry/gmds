//
// Created by bourmaudp on 04/01/23.
//

#include <gmds/blocking_Quality/LinkerBlockingGeom.h>

using namespace gmds;

LinkerBlockingGeom::LinkerBlockingGeom(Mesh *ABlocks, cad::FACManager *AGeom):
  m_blocks(ABlocks),
  m_geom(AGeom){;}

void
LinkerBlockingGeom::execute(cad::GeomMeshLinker *ALinker) 
{
	std::vector<cad::GeomSurface*> surfaces;
	std::vector<cad::GeomCurve*> curves;
	std::vector<cad::GeomPoint*> points;

	m_geom->getSurfaces(surfaces);
	m_geom->getCurves(curves);
	m_geom->getPoints(points);


	//==================================================================
	// MARK ALL THE BOUNDARY CELL OF THE INIT MESH
	//==================================================================

	std::cout<<"> Start block boundary retrieval"<<std::endl;
	//we get all the nodes that are on the mesh boundary
	BoundaryOperator op(m_blocks);
	auto mark_node_NAN = m_blocks->newMark<Node>();
	auto mark_node_on_pnt = m_blocks->newMark<Node>();
	auto mark_node_on_crv = m_blocks->newMark<Node>();
	auto mark_node_on_srf = m_blocks->newMark<Node>();
	auto mark_edge_on_crv = m_blocks->newMark<Edge>();
	auto mark_edge_on_srf = m_blocks->newMark<Edge>();
	auto mark_face_on_srf = m_blocks->newMark<Face>();

	op.markCellOnGeometry(mark_face_on_srf,
	                      mark_edge_on_srf,
	                      mark_node_on_srf,
	                      mark_edge_on_crv,
	                      mark_node_on_crv,
	                      mark_node_on_pnt,
	                      mark_node_NAN);
	//==================================================================
	//First, we classify each face


	for(auto f_id: m_blocks->faces()){
		Face f= m_blocks->get<Face>(f_id);
		if(m_blocks->isMarked(f, mark_face_on_srf) ){
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
				ALinker->linkFaceToSurface(f_id,min_entity_id);
			}
			else{
				std::cout<<"\t ====> Link error for classifying a face"<<std::endl;
				throw GMDSException("Link error for classifying a face");
			}

		}
	}


	//==================================================================
	//Second, we classify each edge
	for(auto e_id: m_blocks->edges()){
		Edge e= m_blocks->get<Edge>(e_id);
		if(m_blocks->isMarked(e, mark_edge_on_crv) ||
		    m_blocks->isMarked(e, mark_edge_on_srf) ){
			//we've got a boundary edge, now we get the 2 boundary faces
			// around e
			std::vector<Node> e_nodes = e.get<Node>();
			std::vector<TCellID> adj_face_id = e.getIDs<Face>();
			std::vector<TCellID> adj_bnd_face_ids;
			for(auto f_id:adj_face_id){
				if(m_blocks->isMarked<Face>(f_id, mark_face_on_srf)){
					adj_bnd_face_ids.push_back(f_id);
				}
			}
			if(adj_bnd_face_ids.size()!=2){
				std::cout<<"ERROR: One boundary edge is adjacent to more "
				          <<"than 2 boundary faces ("
				          <<"face "<<e_id<<")"<<std::endl;
			}
			auto surf0_id = ALinker->getGeomId<Face>(adj_bnd_face_ids[0]);
			auto surf1_id = ALinker->getGeomId<Face>(adj_bnd_face_ids[1]);

			if(surf0_id==surf1_id){
				//edge embedded in the surface
				ALinker->linkEdgeToSurface(e_id,surf1_id);
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
						ALinker->linkEdgeToCurve(e_id, candidate_curves[0]->id());
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

						ALinker->linkEdgeToCurve(e_id, selected_curve->id());

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
	for(auto n_id: m_blocks->nodes()){
		Node n= m_blocks->get<Node>(n_id);
		if(m_blocks->isMarked(n, mark_node_on_pnt) ||
		    m_blocks->isMarked(n, mark_node_on_crv) ||
		    m_blocks->isMarked(n, mark_node_on_srf) ){
			//we've got a boundary node
			//As face and edges are already classified, we used it
			// A node is on a surface if all its adjacent boundary faces are on the same
			// surface. If they are on two it is on a curve potentially
			std::vector<TCellID> adj_faces = n.getIDs<Face>();
			std::vector<TCellID> adj_bnd_faces;
			for(auto id_face:adj_faces){
				if(m_blocks->isMarked<Face>(id_face, mark_face_on_srf)) {
					adj_bnd_faces.push_back(id_face);
				}
			}
			std::set<int> bnd_surfaces;

			for (auto  id_face:adj_bnd_faces){
				bnd_surfaces.insert(ALinker->getGeomId<Face>(id_face));
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
				ALinker->linkNodeToPoint(n_id,min_entity_id);
				std::cout<<"Node "<<n_id<<" is on point "<<min_entity_id<<std::endl;
			}
			else if (min_entity_dim==1){
				on_curve++;
				ALinker->linkNodeToCurve(n_id,min_entity_id);
			}
			else if (min_entity_dim==2){
				on_surf++;
				ALinker->linkNodeToSurface(n_id,min_entity_id);
			}
			else{
				throw GMDSException("Link error for classifying a node");
			}

		}
	}
	std::cout<<"  info [node classified on points "<<on_pnt;
	std::cout<<", on curves "<<on_curve;
	std::cout<<", on surfs "<<on_surf<<"]"<<std::endl;

	ALinker->writeVTKDebugMesh("linker_debug_1.vtk");

}
LinkerBlockingGeom::~LinkerBlockingGeom() = default;
