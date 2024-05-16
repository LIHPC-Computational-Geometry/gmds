/*----------------------------------------------------------------------------*/
#include <gmds/blocking/CurvedBlockingClassifier.h>
#include <gmds/utils/Exception.h>
/*----------------------------------------------------------------------------*/
#include <limits>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::blocking;
/*----------------------------------------------------------------------------*/
CurvedBlockingClassifier::CurvedBlockingClassifier(gmds::blocking::CurvedBlocking *ABlocking) : m_blocking(ABlocking), m_geom_model(ABlocking->geom_model()) {}
/*----------------------------------------------------------------------------*/
CurvedBlockingClassifier::~CurvedBlockingClassifier() {}
/*----------------------------------------------------------------------------*/
ClassificationErrors
CurvedBlockingClassifier::detect_classification_errors(ClassificationErrors &AErrors)
{
	// 1) We check the geometric issues first
	std::vector<cad::GeomPoint *> geom_points;
	m_geom_model->getPoints(geom_points);
	for (auto p : geom_points) {
		auto [found, n] = find_node_classified_on(p);
		if (!found) {
			AErrors.non_captured_points.push_back(p->id());
		}
	}
	// 2) We check the geometric curves
	std::vector<cad::GeomCurve *> geom_curves;
	m_geom_model->getCurves(geom_curves);
	for (auto c : geom_curves) {
		auto [found, n] = find_edge_classified_on(c);
		if (!found) {
			AErrors.non_captured_curves.push_back(c->id());
		}
	}
	// 3) We check the geometric surfaces
	std::vector<cad::GeomSurface *> geom_surfaces;
	m_geom_model->getSurfaces(geom_surfaces);
	for (auto s : geom_surfaces) {
		auto [found, n] = find_face_classified_on(s);
		//if all points and all curves are not captured, we dont class the faces so, the surfaces are not captured
		if(AErrors.non_captured_points.size()!=0 && AErrors.non_captured_curves.size()!=0){
			AErrors.non_captured_surfaces.push_back(s->id());
			continue;
		}
		else if (!found) {
			AErrors.non_captured_surfaces.push_back(s->id());
		}
	}
	return AErrors;
}
/*----------------------------------------------------------------------------*/
void
CurvedBlockingClassifier::clear_classification()
{
	GMap3 *gm = m_blocking->gmap();
	for (auto a : gm->attributes<0>()) {
		a.info().geom_id = NullID;
		a.info().geom_dim = 4;
	}
	for (auto a : gm->attributes<1>()) {
		a.info().geom_id = NullID;
		a.info().geom_dim = 4;
	}
	for (auto a : gm->attributes<2>()) {
		a.info().geom_id = NullID;
		a.info().geom_dim = 4;
	}
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void
CurvedBlockingClassifier::classify_nodes(ClassificationErrors &AErrors, const double AMaxDistance, const double APointSnapDistance)
{
	GMap3 *gm = m_blocking->gmap();
	std::vector<cad::GeomPoint *> geom_points;
	std::vector<cad::GeomCurve *> geom_curves;
	std::vector<cad::GeomSurface *> geom_surfaces;
	m_geom_model->getPoints(geom_points);
	m_geom_model->getCurves(geom_curves);
	m_geom_model->getSurfaces(geom_surfaces);

	// initial projection stage
	for (auto &current_node : gm->attributes<0>()) {
		math::Point p = current_node.info().point;

		std::vector<cad::GeomEntity *> cells;
		cells.insert(cells.end(), geom_points.begin(), geom_points.end());
		auto [closest_pnt_dist, closest_pnt_id, closest_pnt_loc] = get_closest_cell(p, cells);
		cells.clear();
		cells.insert(cells.end(), geom_curves.begin(), geom_curves.end());
		auto [closest_curv_dist, closest_curv_id, closest_curv_loc] = get_closest_cell(p, cells);
		cells.clear();
		cells.insert(cells.end(), geom_surfaces.begin(), geom_surfaces.end());
		auto [closest_surf_dist, closest_surf_id, closest_surf_loc] = get_closest_cell(p, cells);

		if (closest_pnt_dist <= closest_curv_dist && closest_pnt_dist <= closest_surf_dist && closest_pnt_dist <= AMaxDistance) {     // On point
			current_node.info().geom_dim = 0;
			current_node.info().geom_id = closest_pnt_id;
			// projection is done next line
			current_node.info().point = closest_pnt_loc;
		}
		else if (closest_curv_dist < closest_pnt_dist && closest_curv_dist <= closest_surf_dist && closest_curv_dist <= AMaxDistance) {     // On curve
			current_node.info().geom_dim = 1;
			current_node.info().geom_id = closest_curv_id;
			// projection is done next line
			current_node.info().point = closest_curv_loc;
		}
		else if (closest_surf_dist < closest_pnt_dist && closest_surf_dist < closest_curv_dist && closest_surf_dist <= AMaxDistance) {     // On surface
			current_node.info().geom_dim = 2;
			current_node.info().geom_id = closest_surf_id;
			// projection is done next line
			current_node.info().point = closest_surf_loc;
		}
		else {
			// the node is not classified and keep is location
			current_node.info().geom_dim = 4;
			current_node.info().geom_id = NullID;
			AErrors.non_classified_nodes.push_back(current_node.info().topo_id);
		}
	}
	// snapping to the closest point if mandatory
	for (auto &current_node : gm->attributes<0>()) {
		// we only consider nodes that are not classified on points
		if (current_node.info().geom_dim != 0) {
			math::Point p = current_node.info().point;
			std::vector<cad::GeomEntity *> cells;
			cells.insert(cells.end(), geom_points.begin(), geom_points.end());
			auto [closest_pnt_dist, closest_pnt_id, closest_pnt_loc] = get_closest_cell(p, cells);
			if (closest_pnt_dist < APointSnapDistance) {     // On point
				current_node.info().geom_dim = 0;
				current_node.info().geom_id = closest_pnt_id;
				// projection is done next line
				current_node.info().point = closest_pnt_loc;
			}
		}
	}
}
/*----------------------------------------------------------------------------*/
void
CurvedBlockingClassifier::classify_edges(gmds::blocking::ClassificationErrors &AErrors)
{
	GMap3 *gm = m_blocking->gmap();
	std::vector<cad::GeomPoint *> geom_points;
	std::vector<cad::GeomCurve *> geom_curves;
	std::vector<cad::GeomSurface *> geom_surfaces;
	m_geom_model->getPoints(geom_points);
	m_geom_model->getCurves(geom_curves);
	m_geom_model->getSurfaces(geom_surfaces);

	// initial projection stage
	for (auto it = gm->attributes<1>().begin(), itend = gm->attributes<1>().end(); it != itend; ++it) {
		std::vector<CurvedBlocking::Node> ending_nodes = m_blocking->get_nodes_of_edge(it);
		auto geo_d0 = ending_nodes[0]->info().geom_dim;
		auto geo_d1 = ending_nodes[1]->info().geom_dim;
		auto geo_i0 = ending_nodes[0]->info().geom_id;
		auto geo_i1 = ending_nodes[1]->info().geom_id;

		/* We list possible configuration of ending nodes classification:
		 * We decide to not classify an edge on a surface during this process.
		 * IMPORTANT: We'll classify the edges in this case, when we'll classify the faces
		 * 1) Nodes are on different geom points. If those points have a common curve, then the
		 *    edge is on this curve. If they have several common curves,we don't know.
		 * 1bis) We decide to don't
		 * 2) Nodes are on different geom points. If those points have no common curve, but a
		 *    common surface, then the edge is on the surface.
		 * 3) Nodes are on the same curve, then is the edge.
		 * 4) One node is on a point P and the other on a curve, then the edge is
       *    on this curve if it is adjacent to point P
		 * 5) Otherwise we don't know how to classify the edge
		 * 6) One node n0 is on a point P and an edge e0 part of the node is almost tangent of a curve, \\TODO
		 * we check if the other node n1 of the edge is unclassified. If that is the case, we classify
		 * the edge e0 and the node n1 on the curve
		 * 7) Close to the configuration 6, but the first node is on a curve. We can use the same principle \\TODO
		 */
		if (geo_d0==0 && geo_d1==0 && geo_i0!=geo_i1){
			//We look for a common curve
			cad::GeomPoint* p0 = m_geom_model->getPoint(geo_i0);
			cad::GeomPoint* p1 = m_geom_model->getPoint(geo_i1);
			auto curve_id = m_geom_model->getCommonCurve(p0,p1);
			if(curve_id!=-1){
				//We have a common curve (CONFIGURATION 1)
				it->info().geom_dim = 1;
				it->info().geom_id = curve_id;
			}
			else{
				// Nothing (CONFIGURATION 5)
				it->info().geom_dim = 4;
				it->info().geom_id = NullID;
				AErrors.non_classified_edges.push_back(it->info().topo_id);
			}
		}
		else if (geo_d0==1 && geo_d1==1 && geo_i0==geo_i1){
			// On the same curve (CONFIGURATION 3)
			it->info().geom_dim = 1;
			it->info().geom_id = geo_i0;
		}
		else if ((geo_d0==0 && geo_d1==1) || (geo_d0==1 && geo_d1==0)){
			//we check if the point is adjacent to the curve
			auto p_id = (geo_d0<geo_d1)?geo_i0:geo_i1;
			auto c_id = (geo_d0>geo_d1)?geo_i0:geo_i1;
			cad::GeomPoint* p = m_geom_model->getPoint(p_id);
			cad::GeomCurve* c = m_geom_model->getCurve(c_id);
			std::vector<cad::GeomPoint*> c_points = c->points();
			if(c_points[0]==p ||c_points[1]==p){
				// On a curve (CONFIGURATION 4)
				it->info().geom_dim = 1;
				it->info().geom_id = c_id;
			}
		}
		else {
			// Nothing (CONFIGURATION 5)
			it->info().geom_dim = 4;
			it->info().geom_id = NullID;
			AErrors.non_classified_edges.push_back(it->info().topo_id);
		}
	}

	//Now that edges are clasified on curves
}
/*----------------------------------------------------------------------------*/
void
CurvedBlockingClassifier::classify_faces(gmds::blocking::ClassificationErrors &AErrors)
{
	GMap3 *gm = m_blocking->gmap();
	std::vector<cad::GeomPoint *> geom_points;
	std::vector<cad::GeomCurve *> geom_curves;
	std::vector<cad::GeomSurface *> geom_surfaces;
	m_geom_model->getPoints(geom_points);
	m_geom_model->getCurves(geom_curves);
	m_geom_model->getSurfaces(geom_surfaces);

	//check if all points and all curves are captured, we can classify the faces
	if(AErrors.non_captured_curves.empty() && AErrors.non_captured_points.empty()) {
		auto faces = m_blocking->get_all_faces();
		auto map_faces_colored = blocking_color_faces();

		//We class the faces with a surface
		for(auto s : geom_surfaces){

			std::vector<cad::GeomPoint *> s_points = s->points();
			std::vector<cad::GeomCurve *> s_curves = s->curves();
			int color_of_this_surface=-1;

			for(auto f : map_faces_colored) {
				int nb_edges_on_curve = 0;
				int nb_nodes_on_point = 0;

				auto nodes_f = m_blocking->get_nodes_of_face(f.first);
				auto edges_f = m_blocking->get_edges_of_face(f.first);
				// We check if the face have 1 node class on a point and 2 edges class on curves
				for (auto n : nodes_f) {
					if (n->info().geom_dim == 0) {
						for (auto p : s_points) {
							if (p->id() == n->info().geom_id) {
								nb_nodes_on_point++;
							}
						}
					}
				}
				for (auto e : edges_f) {
					if (e->info().geom_dim == 1 && m_blocking->get_faces_of_edge(e).size() ==2) {
						for (auto c : s_curves) {
							if (c->id() == e->info().geom_id) {
								nb_edges_on_curve++;
							}
						}
					}
				}

				if (nb_nodes_on_point >= 1 && nb_edges_on_curve >= 2) {
					color_of_this_surface = f.second;
				}
			}
			for(auto f : map_faces_colored){
				//We do nothing if the color of the face is 0, and any edges is class on a curve, because the face is not on the boundary, but its poss
				if(f.second == 0){
						f.first->info().geom_id = -1;
						f.first->info().geom_dim= 4;
				}

				else if(f.second == color_of_this_surface){
					if(f.first->info().geom_dim==4) {
						f.first->info().geom_id = s->id();
						f.first->info().geom_dim = s->dim();
					}
					//We classify the element of the face not class on the surface
					auto nodes_f = m_blocking->get_nodes_of_face(f.first);
					auto edges_f = m_blocking->get_edges_of_face(f.first);
					//First, the nodes
					for(auto n : nodes_f){
						if(n->info().geom_dim ==4){
							n->info().geom_dim =s->dim();
							n->info().geom_id =s->id();
						}
					}
					//And the edges
					for(auto e : edges_f){
						if(e->info().geom_dim ==4){
							e->info().geom_dim =s->dim();
							e->info().geom_id =s->id();
						}
					}
				}
			}
		}
	}
	auto nodes = m_blocking->get_all_nodes();
	AErrors.non_classified_nodes.clear();
	for(auto n : nodes) {
		if (n->info().geom_dim == 4) {
			AErrors.non_classified_nodes.push_back(n->info().topo_id);
		}
	}

	AErrors.non_classified_edges.clear();

	auto edges = m_blocking->get_all_edges();
	for(auto e : edges){
		if(e->info().geom_dim == 4){
			AErrors.non_classified_edges.push_back(e->info().topo_id);
		}
	}
	AErrors.non_classified_faces.clear();

	auto faces = m_blocking->get_all_faces();
	for(auto f : faces){
		if(f->info().geom_dim == 4){
			AErrors.non_classified_faces.push_back(f->info().topo_id);
		}
	}
}
/*----------------------------------------------------------------------------*/
ClassificationErrors
CurvedBlockingClassifier::classify(const double AMaxDistance, const double APointSnapDistance)
{
	ClassificationErrors errors;
	//============ (1) We classify nodes =================
	classify_nodes(errors, AMaxDistance, APointSnapDistance);

	//============ (2) We classify edges =================
	classify_edges(errors);


	detect_classification_errors(errors);

	if(errors.non_captured_points.empty() && errors.non_captured_curves.empty()){


		//============ (2) We classify faces =================
		classify_faces(errors);

		errors.non_captured_surfaces.clear();

		detect_classification_errors(errors);
	}
	else{
		errors.non_classified_nodes.clear();
		for(auto n : m_blocking->get_all_nodes()){
			if(n->info().geom_dim==4) {
				errors.non_classified_nodes.push_back(n->info().topo_id);
			}
		}
		errors.non_classified_edges.clear();
		for(auto e : m_blocking->get_all_edges()){
			if(e->info().geom_dim==4) {
				errors.non_classified_edges.push_back(e->info().topo_id);
			}
		}
		for(auto f : m_blocking->get_all_faces()){
			errors.non_classified_faces.push_back(f->info().topo_id);
		}
	}


	return errors;
}

/*----------------------------------------------------------------------------*/
std::pair<bool, CurvedBlocking::Node>
CurvedBlockingClassifier::find_node_classified_on(cad::GeomPoint *AP)
{
	GMap3 *gm = m_blocking->gmap();
	for (auto it = gm->attributes<0>().begin(), itend = gm->attributes<0>().end(); it != itend; ++it) {
		if (it->info().geom_dim == 0 && it->info().geom_id == AP->id()) return {true, it};
	}
	return {false, gm->attributes<0>().begin()};
}

/*----------------------------------------------------------------------------*/
std::pair<bool, std::vector<CurvedBlocking::Edge>>
CurvedBlockingClassifier::find_edge_classified_on(cad::GeomCurve *AC)
{

	GMap3 *gm = m_blocking->gmap();

	std::vector<CurvedBlocking::Edge> edgesOnCurve;
	std::vector<GMap3::Attribute_handle<0>::type> nodesOnCurve;

	for (auto it = gm->attributes<0>().begin(), itend = gm->attributes<0>().end(); it != itend; ++it) {
		if (it->info().geom_dim == 1 && it->info().geom_id == AC->id()){
			nodesOnCurve.push_back(it);
		}
	}
	for (auto it = gm->attributes<1>().begin(), itend = gm->attributes<1>().end(); it != itend; ++it) {
		if (it->info().geom_dim == 1 && it->info().geom_id == AC->id()){
			edgesOnCurve.push_back(it);
		}
	}
	for(auto edge : edgesOnCurve){
		std::vector<CurvedBlocking::Node> nodes;
		nodes.reserve(2);
		Dart3 d = edge->dart();
		nodes[0] = gm->attribute<0>(d);
		nodes[1] = gm->attribute<0>(gm->alpha<0>(d));

		auto point0 = AC->points()[0];
		auto point1 = AC->points()[1];

		if(nodes[0]->info().geom_dim == 0 && (nodes[1]->info().geom_dim == 0 )
			&& ((nodes[0]->info().geom_id == point0->id() ) ||(nodes[0]->info().geom_id == point1->id() ))
			&& ((nodes[1]->info().geom_id == point0->id() ) ||(nodes[1]->info().geom_id == point1->id() )) ){

			return {true,edgesOnCurve};
		}

		else if(nodes[0]->info().geom_dim == 0 &&
				  (nodes[0]->info().geom_id==point0->id() || nodes[0]->info().geom_id == point1->id()) &&
				  (nodes[1]->info().geom_dim == 1  && nodes[1]->info().geom_id == AC->id())){

			auto edgesOfN1 = m_blocking->get_edges_of_node(nodes[1]);
			unsigned int nbEdgesOnCurve=0;
			for(auto e : edgesOfN1){
				if(e->info().geom_dim == 1 && e->info().geom_id==AC->id()){
					nbEdgesOnCurve++;
				}
			}
			if(nbEdgesOnCurve==2){
				return{ true,edgesOnCurve};
			}
		}

		else if(nodes[1]->info().geom_dim == 0 &&
				  (nodes[1]->info().geom_id == point0->id() || nodes[1]->info().geom_id==point1->id()) ){

			auto edgesOfN0 = m_blocking->get_edges_of_node(nodes[0]);
			unsigned int nbEdgesOnCurve=0;
			for(auto e : edgesOfN0){
				if(e->info().geom_dim == 1 && e->info().geom_id==AC->id()){
					nbEdgesOnCurve++;
				}
			}
			if(nbEdgesOnCurve==2){
				return{ true,edgesOnCurve};
			}

		}
		else if((nodes[0]->info().geom_dim=1) && (nodes[1]->info().geom_dim=1)){
			unsigned int nbEdgesOnCurve=0;
			auto edgesOfN0 = m_blocking->get_edges_of_node(nodes[0]);
			for(auto e : edgesOfN0){
				if(e->info().geom_dim == 1 && e->info().geom_id==AC->id()){
					nbEdgesOnCurve++;
				}
			}
			auto edgesOfN1 = m_blocking->get_edges_of_node(nodes[1]);
			for(auto e : edgesOfN1){
				if(e->info().geom_dim == 1 && e->info().geom_id==AC->id()){
					nbEdgesOnCurve++;
				}
			}

			if(nbEdgesOnCurve==4){
				return{ true,edgesOnCurve};
			}
		}
	}
	return {false, edgesOnCurve};
}
/*----------------------------------------------------------------------------*/
std::pair<bool, std::vector<CurvedBlocking::Face>>
CurvedBlockingClassifier::find_face_classified_on(cad::GeomSurface *AS)
{
	std::vector<CurvedBlocking::Face> facesOnSurface;
	auto allFaces = m_blocking->get_all_faces();

	for(auto f : allFaces){
				if(f->info().geom_dim == AS->dim() && f->info().geom_id == AS->id()){
			facesOnSurface.push_back(f);
		}
	}
	auto loopOfS = AS->loops();
	auto pointsOfS = AS->points();
	auto curvesOfS = AS->curves();
	std::set<int> listCurvesCapt;


	if(facesOnSurface.size()==0){
		return {false, facesOnSurface};
	}

	//We get all id curves capt by the edges of all the faces class on the surface
	for(auto f_S : facesOnSurface){
		auto edges_F = m_blocking->get_edges_of_face(f_S);
		for (auto e : edges_F){
			if(e->info().geom_dim == 1){
				listCurvesCapt.insert(e->info().geom_id);
			}
		}
	}
	for(auto l : loopOfS){
		for(auto c : l){
			if(listCurvesCapt.count(c->id()) ==0){
				return {false, facesOnSurface};
			}
		}
	}
	return {true, facesOnSurface};



}
/*----------------------------------------------------------------------------*/
std::tuple<double, int, math::Point>
CurvedBlockingClassifier::get_closest_cell(const math::Point &AP, const std::vector<cad::GeomEntity *> &AGeomCells)
{
	math::Point closest_point = AGeomCells[0]->closestPoint(AP);
	double closest_distance = AP.distance(closest_point);
	int closest_id = AGeomCells[0]->id();
	for (auto geom_crv : AGeomCells) {
		math::Point current_point = geom_crv->closestPoint(AP);
		double current_dist = AP.distance(current_point);
		if (current_dist < closest_distance) {
			closest_distance = current_dist;
			closest_point = current_point;
			closest_id = geom_crv->id();
		}
	}
	return {closest_distance, closest_id, closest_point};
}

/*----------------------------------------------------------------------------*/
bool CurvedBlockingClassifier::boundary_surface_captured(cad::GeomSurface *AS){
	auto curves = AS->curves();
	for(auto c : curves){
		auto [found, n] = find_edge_classified_on(c);
		if (!found) {
			return false ;
		}
	}
	return true;
}

/*----------------------------------------------------------------------------*/
std::map<CurvedBlocking::Face,int>
CurvedBlockingClassifier::blocking_color_faces()
{
	std::map<CurvedBlocking::Face,int> faces_colored;
	auto allFaces = m_blocking->get_all_faces();
	auto allEdges = m_blocking->get_all_edges();

	for(auto aF : allFaces){
		faces_colored.insert(std::pair<CurvedBlocking::Face,int>(aF,0));
	}
	int current_color = 0;

	bool finish = false;

	while(!finish){
		current_color ++;
		CurvedBlocking::Face currentFace = NULL;
		for(auto aFC : faces_colored){
			if(aFC.second == 0 && (m_blocking->get_blocks_of_face(aFC.first).size() == 1)){
				currentFace = aFC.first;
				break;
			}
		}
		if(currentFace==NULL){
			finish=true;
		}
		else{
			std::set<TCellID> front;
			front.insert(currentFace->info().topo_id);
			while (!front.empty()){
				TCellID current_id = *front.begin();
				front.erase(front.begin());
				CurvedBlocking::Face aFace;
				for(auto f : allFaces) {
					if (f->info().topo_id == current_id) {
						aFace = f;
						break;
					}
				}
				faces_colored[aFace]=current_color;
				auto edges_f = m_blocking->get_edges_of_face(aFace);
				for(auto e : edges_f){
					//if edge not classified on a curve
					if(e->info().geom_dim!=1){
						//We get her faces
						auto faces_e = m_blocking->get_faces_of_edge(e);
						for(auto f : faces_e){
							if(faces_colored[f]==0 && m_blocking->get_blocks_of_face(f).size()==1){
								front.insert(f->info().topo_id);
							}
						}
					}
				}
			}
		}
	}
	return faces_colored;
}
/*----------------------------------------------------------------------------*/

std::pair<CurvedBlocking::Edge, double>
CurvedBlockingClassifier::get_cut_info(int pointId, std::vector<std::vector<CurvedBlocking::Edge>>& AllEdges)
{
	std::pair<CurvedBlocking::Edge,double> paramCut;

	//============================================
	auto noCaptPoint0 = m_geom_model->getPoint(pointId);
	gmds::math::Point p(noCaptPoint0->X(),noCaptPoint0->Y(),noCaptPoint0->Z());

	auto listEdgesPara = m_blocking->get_all_sheet_edge_sets();
	std::vector<gmds::blocking::CurvedBlocking::Edge > listEdgesSplitable;
	double distMini = std::numeric_limits<double>::max();
	for(auto edges : listEdgesPara){
		auto projInfo = m_blocking->get_projection_info(p,edges);
		for(int i =0; i< projInfo.size();i++){
			if(projInfo[i].second<1 && projInfo[i].second>0 && projInfo[i].first <distMini){
				paramCut.first = edges.at(i);
				paramCut.second = projInfo[i].second;
				distMini=projInfo[i].first;
			}
		}
	}

	return paramCut;
}
/*----------------------------------------------------------------------------*/
std::vector<std::pair<TCellID ,double>>
CurvedBlockingClassifier::list_Possible_Cuts()
{
    std::vector<std::pair<TCellID ,double>> list_actions;
    ClassificationErrors no_capt_elements;
	 this->detect_classification_errors(no_capt_elements);
    auto no_points_capt = no_capt_elements.non_captured_points;
    auto no_curves_capt = no_capt_elements.non_captured_curves;

    auto allEdges = m_blocking->get_all_edges();


    for(auto p : no_points_capt){
        if(m_blocking->check_capt_element(p,0)){
            auto point = m_geom_model->getPoint(p);
            math::Point mathPoint(point->X(),point->Y(),point->Z());
            auto projInfo = m_blocking->get_projection_info(mathPoint,allEdges);
            auto action = m_blocking->get_cut_info(p);
			   bool pairOnList = false;
			   for(auto it : list_actions) {
				   if (it.first == action.first->info().topo_id && it.second == action.second) {
					   pairOnList = true;
				   }
			   }
			   if(pairOnList == false && action.second > 0.01 && action.second<0.99){
				   std::pair<TCellID ,double> actionG (action.first->info().topo_id,action.second);
					   list_actions.push_back(actionG);
			   }
        }
    }

    for(auto c : no_curves_capt){
        if(m_blocking->check_capt_element(c,1)){
            std::vector<std::pair<CurvedBlocking::Edge ,double>> list_actions_curve;
            auto theCurve = m_geom_model->getCurve(c);
            gmds::TCoord minXYX[3];
            gmds::TCoord maxXYX[3];

            theCurve->computeBoundingBox(minXYX,maxXYX);

            gmds::math::Point minPoint(minXYX[0],minXYX[1],minXYX[2]);
            if(m_blocking->get_cut_info(minPoint).second!=2){
                list_actions_curve.push_back(m_blocking->get_cut_info(minPoint));
            }

            gmds::math::Point maxPoint(maxXYX[0],maxXYX[1],maxXYX[2]);
            if(m_blocking->get_cut_info(maxPoint).second!=2){
                list_actions_curve.push_back(m_blocking->get_cut_info(maxPoint));
            }

            gmds::math::Point p0 (minXYX[0],minXYX[1],maxXYX[2]);
            if(m_blocking->get_cut_info(p0).second!=2){
                list_actions_curve.push_back(m_blocking->get_cut_info(p0));
            }

            gmds::math::Point p1 (maxXYX[0],minXYX[1],maxXYX[2]);
            if(m_blocking->get_cut_info(p1).second!=2){
                list_actions_curve.push_back(m_blocking->get_cut_info(p1));
            }

            gmds::math::Point p2 (maxXYX[0],minXYX[1],minXYX[2]);
            if(m_blocking->get_cut_info(p2).second!=2){
                list_actions_curve.push_back(m_blocking->get_cut_info(p2));
            }

            gmds::math::Point p3 (minXYX[0],maxXYX[1],maxXYX[2]);
            if(m_blocking->get_cut_info(p3).second!=2){
                list_actions_curve.push_back(m_blocking->get_cut_info(p3));
            }

            gmds::math::Point p4 (minXYX[0],maxXYX[1],minXYX[2]);
            if(m_blocking->get_cut_info(p4).second!=2){
                list_actions_curve.push_back(m_blocking->get_cut_info(p4));
            }

            gmds::math::Point p5 (maxXYX[0],maxXYX[1],minXYX[2]);
            if(m_blocking->get_cut_info(p5).second!=2){
                list_actions_curve.push_back(m_blocking->get_cut_info(p5));
            }


            for(auto i : list_actions_curve){
                bool pairOnList = false;
                for(auto it : list_actions) {
                    if (it.first == i.first->info().topo_id && it.second == i.second) {
                        pairOnList = true;
                    }
                }
                if(pairOnList == false){
                    std::pair<TCellID ,double> actionG (i.first->info().topo_id,i.second);
                    list_actions.push_back(actionG);
                }
            }
        }
    }

    return list_actions;
}
/*----------------------------------------------------------------------------*/
bool
CurvedBlockingClassifier::checkValidity(ClassificationErrors &AErrors) {
    if(m_blocking->get_all_nodes().size() < m_geom_model->getNbPoints() ||
        m_blocking->get_all_edges().size() < m_geom_model->getNbCurves() ||
        m_blocking->get_all_faces().size() < m_geom_model->getNbSurfaces()){

        return false;
    }
    else if(AErrors.non_captured_points.size() != 0 ||
            AErrors.non_captured_curves.size() != 0 ||
            AErrors.non_captured_surfaces.size() != 0){
        return false;
    }
    return true;
}
/*----------------------------------------------------------------------------*/
