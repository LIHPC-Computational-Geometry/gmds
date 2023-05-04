//
// Created by bourmaudp on 21/02/23.
//

#include <gmds/blocking_Quality/ValidBlocking.h>

using namespace gmds;

ValidBlocking::ValidBlocking(Mesh *ABlocks, cad::FACManager *AGeom, cad::GeomMeshLinker *ALinker,std::map<std::vector<TCellID>,int> *elementsNoClassified) :
  m_blocks(ABlocks), m_geom(AGeom), a_linker(ALinker),elements_No_Classified(elementsNoClassified)
{
	;
}

ValidBlocking::~ValidBlocking() = default;

bool
ValidBlocking::execute()
{

	// check node validity
	/*if (!checkValidNodes()) {
		return false;
	}*/
	// check points have a node linked
	if (!checkPoints()) {
		return false;
	}

	// check edges validity
	/*if (!checkValidEdges()) {
		return false;
	}*/
	// check curves have path edges linked
	if (!checkCurves()) {
		return false;
	}

	// check faces validity
	/*if (!checkValidFaces()) {
		return false;
	}*/

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
	if (dim_element == cad::GeomMeshLinker::LINK_POINT && (std::count(list_point_check.begin(), list_point_check.end(), id_element) || list_point_check.empty())) {
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
		if (a_linker->getGeomDim(n0) == cad::GeomMeshLinker::LINK_POINT) {
			auto curve = m_geom->getCurve(a_linker->getGeomId(e));
			auto curve_id = curve->id();
			auto p0 = m_geom->getPoint(curve->points()[0]->id());
			auto p1 = m_geom->getPoint(curve->points()[1]->id());
			;
			if (!(a_linker->getGeomId(n0) == p0->id() || a_linker->getGeomId(n0) == p1->id())) {
				return false;
			}
		}
		else if (a_linker->getGeomDim(n1) == cad::GeomMeshLinker::LINK_POINT) {
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
		if (dim_link == cad::GeomMeshLinker::LINK_CURVE && id_link == id_curve) {
			vectorElements.push_back(std::make_pair(0, n_id));
		}
	}

	// get edges link to the curve
	for (auto e_id : m_blocks->edges()) {
		auto dim_link = a_linker->getGeomDim(m_blocks->get<Edge>(e_id));
		auto id_link = a_linker->getGeomId(m_blocks->get<Edge>(e_id));
		if (dim_link == cad::GeomMeshLinker::LINK_CURVE  && id_link == id_curve) {
			vectorElements.push_back(std::make_pair(1, e_id));
		}
	}
	return vectorElements;
}

bool
ValidBlocking::checkValidFaces()
{

	for (auto f_id : m_blocks->faces()) {
		if ((!checkValidityFace(f_id)) && (!checkValidityElementsFace(f_id))) {
			return false;
		}
	}
	if (!checkFacesLinked()) {
		return false;
	}
	return true;
}

bool
ValidBlocking::checkFacesLinked()
{
	std::map<std::vector<Face>, cad::GeomSurface *> mapFaces;
	std::vector<cad::GeomSurface *> listSurfaces;
	m_geom->getSurfaces(listSurfaces);

	for (auto s : listSurfaces) {
		std::vector<Face> listFaces;
		std::vector<Edge> listEdges;
		for (auto f : m_blocks->faces()) {
			Face AFace = m_blocks->get<Face>(f);
			auto link_id = a_linker->getGeomId(AFace);
			auto link_dim = a_linker->getGeomDim(AFace);

			if ((link_dim == cad::GeomMeshLinker::LINK_SURFACE) && (link_id == s->id())) {
				auto edgesVect = AFace.get<Edge>();
				for (auto AEdge : edgesVect) {
					if ((AEdge.nbFaces() != 3) && (a_linker->getGeomDim(AEdge) == cad::GeomMeshLinker::LINK_SURFACE)) {
						return false;
					}
				}
			}
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

	// si pas de link avec volume de la geometry, le cas pour le moment
	if (link_f_dim == cad::GeomMeshLinker::NO_LINK) {
		return true;
	}
	// f linked false
	else if (link_f_dim < 3) {
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
ValidBlocking::checkValidityElementsFace(TInt AFaceId)
{
	Face f = m_blocks->get<Face>(AFaceId);
	auto link_f_id = a_linker->getGeomId(f);
	auto link_f_dim = a_linker->getGeomDim(f);

	std::vector<std::pair<TCellID, int>> vectorElementsFace;

	// Get edges of f and link
	auto e0 = f.get<Edge>()[0];
	auto link_e0_id = a_linker->getGeomId(e0);
	auto link_e0_dim = a_linker->getGeomDim(e0);
	vectorElementsFace.push_back(std::make_pair(link_e0_id, link_e0_dim));

	auto e1 = f.get<Edge>()[1];
	auto link_e1_id = a_linker->getGeomId(e1);
	auto link_e1_dim = a_linker->getGeomDim(e1);
	vectorElementsFace.push_back(std::make_pair(link_e1_id, link_e1_dim));

	auto e2 = f.get<Edge>()[2];
	auto link_e2_id = a_linker->getGeomId(e2);
	auto link_e2_dim = a_linker->getGeomDim(e2);
	vectorElementsFace.push_back(std::make_pair(link_e2_id, link_e2_dim));

	auto e3 = f.get<Edge>()[3];
	auto link_e3_id = a_linker->getGeomId(e3);
	auto link_e3_dim = a_linker->getGeomDim(e3);
	vectorElementsFace.push_back(std::make_pair(link_e3_id, link_e3_dim));

	// linked on surface
	if (link_f_dim == cad::GeomMeshLinker::LINK_SURFACE) {
		auto surface = m_geom->getSurface(link_f_id);
		auto curves = surface->curves();
		std::vector<int> id_curves;
		for (auto c : curves) {
			id_curves.push_back(c->id());
		}
		for (auto e : vectorElementsFace) {
			// Edges link on curves
			if (e.second == 2) {
				if (!(std::count(id_curves.begin(), id_curves.end(), e.first))) {
					return false;
				}
			}
			// Edges link on surfaces
			else if (e.second == 3) {
				if (e.first != link_f_id) {
					return false;
				}
			}
		}
	}
	// linked on volume
	else {
		auto volume = m_geom->getVolume(link_f_id);
	}
	return true;
}

bool
ValidBlocking::checkGeomValidity()
{
	if(!checkPoints() ) {
		returnNoValidsPoints();
		return false;
	}
	else if(!checkCurves()){
		returnNoValidsCurves();
		return false;
	}
	return true;

}

bool
ValidBlocking::checkPoints()
{
	std::vector<cad::GeomPoint *> pointsList;
	m_geom->getPoints(pointsList);

	std::vector<TCellID> listPointsNoClassified;

	for (auto p : pointsList) {
		if(!checkAPoint(p)){
			listPointsNoClassified.push_back(p->id());
		}
	}
	if(!listPointsNoClassified.empty()){
		return false;
	}
	else {
		return true;
	}
}

bool ValidBlocking::checkAPoint(cad::GeomPoint *APoint){
	unsigned int nbNodesLinked = 0;
	int APointId = APoint->id();
	for (auto n : m_blocks->nodes()) {
		Node ANode = m_blocks->get<Node>(n);
		auto link_id = a_linker->getGeomId(ANode);
		auto link_dim = a_linker->getGeomDim(ANode);
		if ((link_dim == cad::GeomMeshLinker::LINK_POINT) && (link_id == APointId)) {
			nbNodesLinked++;
		}
	}
	if (nbNodesLinked != 1) {
		return false;
	}
	else{return true;}
}



bool
ValidBlocking::checkSuiteOfElements(std::vector<std::pair<int, TCellID>> AElementsOnCurve)
{
	std::vector<Node> vectorNodes;
	std::vector<Edge> vectorEdges;
	for (auto item : AElementsOnCurve) {
		if (item.first == 0) {
			Node ANode = m_blocks->get<Node>(item.second);
			vectorNodes.push_back(m_blocks->get<Node>(item.second));
		}
		else if (item.first == 1) {
			Edge AEdge = m_blocks->get<Edge>(item.second);
			vectorEdges.push_back(AEdge);
		}
	}

	if (vectorNodes.size() + 1 != vectorEdges.size()) {
		return false;
	}

	for (auto ANode : vectorNodes) {
		int nbEdgesFromANode = 0;
		auto listAllEdgesNodes = ANode.get<Edge>();
		std::vector<Edge> listEdgesNodes;
		for (auto Edge : listAllEdgesNodes) {
			if (std::count(vectorEdges.begin(), vectorEdges.end(), Edge)) {
				listEdgesNodes.push_back(Edge);
			}
		}
		for (auto AEdge : listEdgesNodes) {
			Node n0 = AEdge.get<Node>()[0];
			Node n1 = AEdge.get<Node>()[1];

			if ((n0.id() == ANode.id() || n1.id() == ANode.id())) {
				nbEdgesFromANode++;
			}
		}
		if (nbEdgesFromANode != 2) {
			return false;
		}
	}
	return true;
}

bool
ValidBlocking::checkExtremityCurve(std::vector<Edge> AEdgesVector, std::vector<Node> ANodesVector, cad::GeomCurve *ACurve)
{
	for(auto AEdge : AEdgesVector){
		Node n0 = AEdge.get<Node>()[0];
		Node n1 = AEdge.get<Node>()[1];
		auto pointsVector = ACurve->points();
		std::vector<int> idPointVector;
		for(auto p : pointsVector){
			idPointVector.push_back(p->id());
		}
		if(!(std::count(ANodesVector.begin(), ANodesVector.end(), n0) && std::count(ANodesVector.begin(), ANodesVector.end(), n1))){
			auto n0LinkId = a_linker->getGeomId(n0);

			auto n1LinkId = a_linker->getGeomId(n1);
			if(!( std::count(idPointVector.begin(), idPointVector.end(), n0LinkId) ||  std::count(idPointVector.begin(), idPointVector.end(), n1LinkId))){
				return false;
			}
		}
	}
	return true;
}

bool
ValidBlocking::checkCurves()
{
	std::vector<cad::GeomCurve *> curvesList;
	m_geom->getCurves(curvesList);
	for (auto ACurve : curvesList) {
		if(!checkACurve(ACurve)){
			return false;
		}
	}
	return true;
}

bool
ValidBlocking::checkACurve(cad::GeomCurve *ACurve)
{
	int ACurveId = ACurve->id();
	std::vector<std::pair<int, TCellID>> vectElmntCurve = elementsOnCurve(ACurve);
	std::vector<Node> vectorNodes;
	std::vector<Edge> vectorEdges;
	for (auto item : vectElmntCurve) {
		if (item.first == 0) {
			Node ANode = m_blocks->get<Node>(item.second);
			vectorNodes.push_back(m_blocks->get<Node>(item.second));
		}
		else if (item.first == 1) {
			Edge AEdge = m_blocks->get<Edge>(item.second);
			Node n0 = AEdge.get<Node>()[0];
			Node n1 = AEdge.get<Node>()[1];
			auto link_AEdge_ID = a_linker->getGeomId(AEdge);
			auto link_AEdge_Dim = a_linker->getGeomDim(AEdge);

			if (link_AEdge_ID != ACurveId && link_AEdge_Dim != cad::GeomMeshLinker::LINK_CURVE) {
				return false;
			}
			vectorEdges.push_back(AEdge);
		}
	}
	// check suite of elements
	if (!checkSuiteOfElements(vectElmntCurve)) {
		return false;
	}
	// check if extremity of suite linked of extremity of curve
	if (!checkExtremityCurve(vectorEdges, vectorNodes, ACurve)) {
		return false;
	}
}

std::map<int,TCellID>
ValidBlocking::returnNoValidsPoints(){
	std::map<int,TCellID> noValidsPoints;
	return noValidsPoints;
}

std::map<int,TCellID>
ValidBlocking::returnNoValidsCurves(){
	std::map<int,TCellID> noValidsCurves;
	return noValidsCurves;
}