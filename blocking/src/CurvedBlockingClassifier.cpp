/*----------------------------------------------------------------------------*/
#include <gmds/blocking/CurvedBlockingClassifier.h>
#include <gmds/utils/Exception.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::blocking;
/*----------------------------------------------------------------------------*/
CurvedBlockingClassifier::CurvedBlockingClassifier(gmds::blocking::CurvedBlocking *ABlocking) : m_blocking(ABlocking), m_geom_model(ABlocking->geom_model()) {}
/*----------------------------------------------------------------------------*/
CurvedBlockingClassifier::~CurvedBlockingClassifier() {}
/*----------------------------------------------------------------------------*/
ClassificationErrors
CurvedBlockingClassifier::detect_classification_errors()
{
	ClassificationErrors errors;
	// 1) We check the geometric issues first
	std::vector<cad::GeomPoint *> geom_points;
	m_geom_model->getPoints(geom_points);
	for (auto p : geom_points) {
		auto [found, n] = find_node_classified_on(p);
		if (!found) {
			errors.non_captured_points.push_back(p->id());
		}
	}
	// 2) We check the geometric curves
	std::vector<cad::GeomCurve *> geom_curves;
	m_geom_model->getCurves(geom_curves);
	for (auto c : geom_curves) {
		auto [found, n] = find_edge_classified_on(c);
		if (!found) {
			errors.non_captured_curves.push_back(c->id());
		}
	}
	return errors;
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
		 * 1) Nodes are on different geom points. If those points have a common curve, then the
		 *    edge is on this curve. If they have several common curves,we don't know.
		 * 2) Nodes are on different geom points. If those points have no common curve, but a
		 *    common surface, then the edge is on this surface.
		 * 3) Nodes are on the same curve or surface, then is the edge.
		 * 4) One node is on a point P and the other on a curve or a surface, then the edge is
       *    on this curve or surface if it is adjacent to point P
		 * 5) One node is on a curve C and the other on a surface S, then the edge is
		 *    on S if S is adjacent to point C
		 * 6) Otherwise we don't know how to classify the edge
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
				// Nothing (CONFIGURATION 6)
				it->info().geom_dim = 4;
				it->info().geom_id = NullID;
				AErrors.non_classified_edges.push_back(it->info().topo_id);
				//We look for a common surface
				//auto surf_ids = m_geom_model->getCommonSurfaces(p0,p1);
				/*if(surf_ids.size()==1){
					//We have a common surface (CONFIGURATION 2)
					std::cout<<"DANS CONFIG 2"<<std::endl;
					it->info().geom_dim = 2;
					it->info().geom_id = surf_ids[0];
				}*/
			}
		}
		else if (geo_d0==1 && geo_d1==1 && geo_i0==geo_i1){
			// On the same curve (CONFIGURATION 3.a)
			it->info().geom_dim = 1;
			it->info().geom_id = geo_i0;
		}
		/*else if (geo_d0==2 && geo_d1==2 && geo_i0==geo_i1){
			// On the same surface (CONFIGURATION 3.b), we decide to don't classify right now, if the surface is not captured
			it->info().geom_dim = 2;
			it->info().geom_id = geo_i0;
		}*/

		else if ((geo_d0==0 && geo_d1==1) || (geo_d0==1 && geo_d1==0)){
			//we check if the point is adjacent to the curve
			auto p_id = (geo_d0<geo_d1)?geo_d0:geo_d1;
			auto c_id = (geo_d0>geo_d1)?geo_d0:geo_d1;
			cad::GeomPoint* p = m_geom_model->getPoint(p_id);
			cad::GeomCurve* c = m_geom_model->getCurve(c_id);
			std::vector<cad::GeomPoint*> c_points = c->points();
			if(c_points[0]==p ||c_points[1]==p){
				// On a curve (CONFIGURATION 4.a)
				it->info().geom_dim = 1;
				it->info().geom_id = c_id;
			}
		}
		else if ((geo_d0==0 && geo_d1==2) || (geo_d0==2 && geo_d1==0)){
			//we check if the point is adjacent to the surface
			auto p_id = (geo_d0<geo_d1)?geo_d0:geo_d1;
			auto s_id = (geo_d0>geo_d1)?geo_d0:geo_d1;
			cad::GeomPoint* p = m_geom_model->getPoint(p_id);
			cad::GeomSurface* s = m_geom_model->getSurface(s_id);
			std::vector<cad::GeomPoint*> s_points = s->points();
			for(auto s_p : s_points){
				if(s_p==p){
					// On a surface (CONFIGURATION 4.b)
					it->info().geom_dim = 2;
					it->info().geom_id = s_id;
				}
			}
		}
		else if ((geo_d0==1 && geo_d1==2) || (geo_d0==2 && geo_d1==1)){
			//we check if the point is adjacent to the surface
			auto c_id = (geo_d0<geo_d1)?geo_d0:geo_d1;
			auto s_id = (geo_d0>geo_d1)?geo_d0:geo_d1;
			cad::GeomCurve* c = m_geom_model->getCurve(c_id);
			cad::GeomSurface* s = m_geom_model->getSurface(s_id);
			std::vector<cad::GeomCurve*> s_curves = s->curves();
			for(auto s_c : s_curves){
				if(s_c==c){
					// On a surface (CONFIGURATION 5)
					it->info().geom_dim = 2;
					it->info().geom_id = s_id;
				}
			}
		}
		else {
			// Nothing (CONFIGURATION 6)
			it->info().geom_dim = 4;
			it->info().geom_id = NullID;
			AErrors.non_classified_edges.push_back(it->info().topo_id);
		}
	}
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


	auto errors_Capt = detect_classification_errors();
	//check if all points and all curves are captured
	if(detect_classification_errors().non_captured_surfaces.empty() && detect_classification_errors().non_captured_points.empty()){

	}


	// initial projection stage
	for (auto it = gm->attributes<2>().begin(), itend = gm->attributes<2>().end(); it != itend; ++it) {

		std::vector<CurvedBlocking::Node> face_nodes = m_blocking->get_nodes_of_face(it);
		std::vector<CurvedBlocking::Edge> face_edges = m_blocking->get_edges_of_face(it);
		auto geo_n_d0 = face_nodes[0]->info().geom_dim;
		auto geo_n_d1 = face_nodes[1]->info().geom_dim;
		auto geo_n_d2 = face_nodes[2]->info().geom_dim;
		auto geo_n_d3 = face_nodes[3]->info().geom_dim;
		auto geo_n_i0 = face_nodes[0]->info().geom_id;
		auto geo_n_i1 = face_nodes[1]->info().geom_id;
		auto geo_n_i2 = face_nodes[2]->info().geom_id;
		auto geo_n_i3 = face_nodes[3]->info().geom_id;


		auto geo_e_d0 = face_edges[0]->info().geom_dim;
		auto geo_e_d1 = face_edges[1]->info().geom_dim;
		auto geo_e_d2 = face_edges[2]->info().geom_dim;
		auto geo_e_d3 = face_edges[3]->info().geom_dim;
		auto geo_e_i0 = face_edges[0]->info().geom_id;
		auto geo_e_i1 = face_edges[1]->info().geom_id;
		auto geo_e_i2 = face_edges[2]->info().geom_id;
		auto geo_e_i3 = face_edges[3]->info().geom_id;

		std::cout<<"Noeud"<<std::endl;
		for(int i=0;i<4;i++){
			std::cout<<face_nodes[i]->info().geom_id<<std::endl;
		}
		std::cout<<"Arete"<<std::endl;
		for(int i=0;i<4;i++){
			std::cout<<face_edges[i]->info().geom_id<<std::endl;
		}

		/* We list possible configuration of elements classification:
		 * 1) Nodes are on different geom points. If those points have a common curve, then the
		 *    edge is on this curve. If they have several common curves,we don't know.
		 * 2) Nodes are on different geom points. If those points have no common curve, but a
		 *    common surface, then the edge is on this surface.
		 * 3) Nodes are on the same curve or surface, then is the edge.
		 * 4) One node is on a point P and the other on a curve or a surface, then the edge is
       *    on this curve or surface if it is adjacent to point P
		 * 5) One node is on a curve C and the other on a surface S, then the edge is
		 *    on S if S is adjacent to point C
		 * 6) Otherwise we don't know how to classify the edge
		 */
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

	//============ (2) We classify faces =================
	classify_faces(errors);

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
