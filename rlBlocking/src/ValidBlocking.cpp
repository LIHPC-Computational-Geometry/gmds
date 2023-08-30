//
// Created by bourmaudp on 21/02/23.
//

#include <gmds/rlBlocking/ValidBlocking.h>

using namespace gmds;

ValidBlocking::ValidBlocking(Mesh *ABlocks, cad::FACManager *AGeom, cad::GeomMeshLinker *ALinker) : m_blocks(ABlocks), m_geom(AGeom), a_linker(ALinker)
{
	;
}

ValidBlocking::~ValidBlocking() = default;

bool
ValidBlocking::execute()
{
	// check node validity
	if (!checkValidNodes()) {
		return false;
	}

	// check edges validity
	if (!checkValidEdges()) {
		return false;
	}

	// check faces validity
	if (!checkValidFaces()) {
		return false;
	}

	return true;
}

bool
ValidBlocking::checkValidNodes()
{
	auto nb_points = m_geom->getNbPoints();
	TInt check_nb_points = 0;
	for (auto n_id : m_blocks->nodes()) {
		if (checkValidityNode(n_id)) {
			check_nb_points++;
		}
	}
	if (check_nb_points != nb_points) {
		return false;
	}
	else {
		return true;
	}
}
bool
ValidBlocking::checkValidityNode(TInt ANodeId)
{
	Node n = m_blocks->get<Node>(ANodeId);
	auto dim_element = a_linker->getGeomDim(n);
	auto id_element = a_linker->getGeomId(n);
	auto nb_points = m_geom->getNbPoints();
	TInt check_nb_points = 0;
	std::vector<int> list_point_check;
	if (dim_element == 1 && (std::count(list_point_check.begin(), list_point_check.end(), id_element) || list_point_check.empty())) {
		check_nb_points++;
		list_point_check.push_back(id_element);
		return true;
	}
	else {
		return false;
	}
}

bool
ValidBlocking::checkValidEdges()
{
	auto nb_curves = m_geom->getNbCurves();
	TInt check_nb_curves = 0;
	std::vector<cad::GeomCurve *> all_curves;
	m_geom->getCurves(all_curves);
	std::map<std::vector<std::pair<int, TCellID>>, TCellID> mapCurvesElements;
	std::vector<TCellID> vectorNodes;
	std::vector<TCellID> vectorEdges;
	bool checkNode0Extern = false;
	bool checkNode1Extern = false;

	for (auto c_pointeur : all_curves) {
		mapCurvesElements.insert(std::make_pair(gmds::ValidBlocking::elementsOnCurve(c_pointeur), c_pointeur->id()));
	}
	for (auto item : mapCurvesElements) {

		for (auto itemBis : item.first) {
			if (itemBis.first == 0) {
				vectorNodes.push_back(itemBis.second);
			}
			if (itemBis.first == 1) {
				vectorEdges.push_back(itemBis.second);
			}
		}
	}
	for (auto item : vectorEdges) {
		Edge e = m_blocks->get<Edge>(item);
		Node n0 = e.get<Node>()[0];
		Node n1 = e.get<Node>()[1];
		if (a_linker->getGeomDim(n0) == 1) {
			auto curve = m_geom->getCurve(a_linker->getGeomId(e));
			auto curve_id = curve->id();
			auto p0 = m_geom->getPoint(curve->points()[0]->id());
			auto p1 = m_geom->getPoint(curve->points()[1]->id());
			;
			if (!(a_linker->getGeomId(n0) == p0->id() || a_linker->getGeomId(n0) == p1->id())) {
				return false;
			}
		}
		else if (a_linker->getGeomDim(n1) == 1) {
			auto curve = m_geom->getCurve(a_linker->getGeomId(e));
			auto p0 = m_geom->getPoint(curve->points()[0]->id());
			auto p1 = m_geom->getPoint(curve->points()[1]->id());
			;
			if (!(a_linker->getGeomId(n1) == p0->id() || a_linker->getGeomId(n1) == p1->id())) {
				return false;
			}
		}
		else {
			if (!(std::count(vectorNodes.begin(), vectorNodes.end(), n0.id()) && std::count(vectorNodes.begin(), vectorNodes.end(), n1.id()))) {
				return false;
			}
		}
	}

	return true;
}

bool
ValidBlocking::checkValidityEdge(TInt AEdgeId)
{
	Edge e = m_blocks->get<Edge>(AEdgeId);
	return true;
}

std::vector<std::pair<int, TCellID>>
ValidBlocking::elementsOnCurve(cad::GeomCurve *ACurve)
{
	std::vector<std::pair<int, TCellID>> vectorElements;

	auto id_curve = ACurve->id();
	// get nodes link to the curve
	for (auto n_id : m_blocks->nodes()) {
		auto dim_link = a_linker->getGeomDim(m_blocks->get<Node>(n_id));
		auto id_link = a_linker->getGeomId(m_blocks->get<Node>(n_id));
		if (dim_link == 2 && id_link == id_curve) {
			vectorElements.push_back(std::make_pair(0, n_id));
		}
	}

	// get edges link to the curve
	for (auto e_id : m_blocks->edges()) {
		auto dim_link = a_linker->getGeomDim(m_blocks->get<Edge>(e_id));
		auto id_link = a_linker->getGeomId(m_blocks->get<Edge>(e_id));
		if (dim_link == 2 && id_link == id_curve) {
			vectorElements.push_back(std::make_pair(1, e_id));
		}
	}
	return vectorElements;
}

bool
ValidBlocking::checkValidFaces()
{
	for (auto f_id : m_blocks->faces()) {
		if (!checkValidityFace(f_id)) {
			return false;
		}
	}
	return true;
}

bool
ValidBlocking::checkValidityFace(TInt AFaceId)
{
	Face f = m_blocks->get<Face>(AFaceId);
	auto link_f_id = a_linker->getGeomId(f);
	auto link_f_dim = a_linker->getGeomDim(f);
	// Get edges of f and link
	auto e0 = f.get<Edge>()[0];
	auto link_e0_id = a_linker->getGeomId(e0);
	auto link_e0_dim = a_linker->getGeomDim(e0);

	auto e1 = f.get<Edge>()[1];
	auto link_e1_id = a_linker->getGeomId(e1);
	auto link_e1_dim = a_linker->getGeomDim(e1);

	auto e2 = f.get<Edge>()[2];
	auto link_e2_id = a_linker->getGeomId(e2);
	auto link_e2_dim = a_linker->getGeomDim(e2);

	auto e3 = f.get<Edge>()[3];
	auto link_e3_id = a_linker->getGeomId(e3);
	auto link_e3_dim = a_linker->getGeomDim(e3);

	// f linked false
	if (link_f_dim < 3) {
		return false;
	}

	// f linked on surface
	else if (link_f_dim == 3) {
		return true;
	}
	// f linked on volume
	else {
		return true;
	}
}

bool
checkValidityElementsFace(TInt AFaceId)
{

}
