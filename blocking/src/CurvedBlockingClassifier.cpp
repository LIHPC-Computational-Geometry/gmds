/*----------------------------------------------------------------------------*/
#include <gmds/blocking/CurvedBlockingClassifier.h>
#include <gmds/utils/Exception.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::blocking;
/*----------------------------------------------------------------------------*/
CurvedBlockingClassifier::CurvedBlockingClassifier(gmds::blocking::CurvedBlocking *ABlocking)
: m_blocking(ABlocking), m_geom_model(ABlocking->geom_model()){
}
/*----------------------------------------------------------------------------*/
CurvedBlockingClassifier::~CurvedBlockingClassifier()
{}
/*----------------------------------------------------------------------------*/
ClassificationErrors
CurvedBlockingClassifier::detect_classification_errors()
{
	ClassificationErrors errors;
	// 1) We check the geometric issues first
	std::vector<cad::GeomPoint*> geom_points;
	m_geom_model->getPoints(geom_points);
	for (auto p : geom_points) {
		auto [found, n] = find_node_classified_on(p);
		if (!found) {
			errors.non_captured_points.push_back(p->id());
		}
	}
}
/*----------------------------------------------------------------------------*/
void CurvedBlockingClassifier::clear_classification()
{
	GMap3* gm = m_blocking->gmap();
	for (auto a : gm->attributes<0>()){
		a.info().geom_id=NullID;
		a.info().geom_dim=4;
	}
	for (auto a : gm->attributes<1>()){
		a.info().geom_id=NullID;
		a.info().geom_dim=4;
	}
	for (auto a : gm->attributes<2>()){
		a.info().geom_id=NullID;
		a.info().geom_dim=4;
	}
}
/*----------------------------------------------------------------------------*/
void CurvedBlockingClassifier::classify_nodes(ClassificationErrors& AErrors, const double AMaxDistance,const double APointSnapDistance)
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
		//we only consider nodes that are not classified on points
		if(current_node.info().geom_dim!=0) {
			math::Point p = current_node.info().point;
			std::vector<cad::GeomEntity *> cells;
			cells.insert(cells.end(), geom_points.begin(), geom_points.end());
			auto [closest_pnt_dist, closest_pnt_id, closest_pnt_loc] = get_closest_cell(p, cells);
			std::cout << "Distance = " << closest_pnt_dist << std::endl;
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
ClassificationErrors CurvedBlockingClassifier::classify(const double AMaxDistance,const double APointSnapDistance){
	ClassificationErrors errors;

	//============ (1) We classify nodes =================
	classify_nodes(errors,AMaxDistance, APointSnapDistance);

	return errors;
}

/*----------------------------------------------------------------------------*/
std::pair<bool, CurvedBlocking::Node>
CurvedBlockingClassifier::find_node_classified_on(cad::GeomPoint* AP){
	GMap3* gm = m_blocking->gmap();
	for (auto it = gm->attributes<0>().begin(), itend = gm->attributes<0>().end(); it != itend; ++it) {
		if (it->info().geom_dim==0 && it->info().geom_id==AP->id())
			return {true, it};
	}
	return {false, gm->attributes<0>().begin()};
}

/*----------------------------------------------------------------------------*/
std::tuple<double, int, math::Point>
CurvedBlockingClassifier::get_closest_cell(const math::Point& AP, const std::vector<cad::GeomEntity*>& AGeomCells){
	math::Point closest_point = AGeomCells[0]->closestPoint(AP);
	double closest_distance = AP.distance(closest_point);
	int closest_id = AGeomCells[0]->id();
	for(auto geom_crv:AGeomCells){
		math::Point current_point = geom_crv->closestPoint(AP);
		double current_dist = AP.distance(current_point);
		if(current_dist<closest_distance){
			closest_distance = current_dist;
			closest_point = current_point;
			closest_id = geom_crv->id();
		}
	}
	return {closest_distance,closest_id,closest_point};
}
/*----------------------------------------------------------------------------*/
