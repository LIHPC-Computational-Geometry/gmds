/*---------------------------------------------------------------------------*/
/*
 * SingularityGraph.cpp
 *
 *  Created on: 13 juil. 2014
 *      Author: F. Ledoux
 */
/*---------------------------------------------------------------------------*/
#include <gmds/singGraphBuild/SingularityGraph.h>
#include <gmds/singGraphBuild/Tools.h>
/*---------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/math/Numerics.h>
#include <gmds/math/Segment.h>
#include <gmds/math/Triangle.h>
#include <gmds/utils/Exception.h>
/*---------------------------------------------------------------------------*/
using namespace gmds;
using namespace std;
/*---------------------------------------------------------------------------*/
SingularityGraph::SingularityGraph(Mesh *AMesh) : m_mesh(AMesh) {}
/*---------------------------------------------------------------------------*/
SingularityGraph::~SingularityGraph()
{
	for (unsigned int i = 0; i < m_points.size(); i++)
		delete m_points[i];
	for (unsigned int i = 0; i < m_lines.size(); i++)
		delete m_lines[i];
	for (unsigned int i = 0; i < m_patchs.size(); i++)
		delete m_patchs[i];
}
/*---------------------------------------------------------------------------*/
std::vector<SingularityPoint *> &
SingularityGraph::getPoints()
{
	return m_points;
}
/*---------------------------------------------------------------------------*/
std::vector<VolumeSingularityLine *>
SingularityGraph::getVolumeLines()
{
	std::vector<VolumeSingularityLine *> l;
	for (unsigned int i = 0; i < m_lines.size(); i++) {
		if (m_lines[i]->getType() == SingularityLine::VOLUME) l.push_back(dynamic_cast<VolumeSingularityLine *>(m_lines[i]));
	}
	return l;
}
/*---------------------------------------------------------------------------*/
std::vector<SurfaceSingularityLine *>
SingularityGraph::getSurfaceLines()
{
	std::vector<SurfaceSingularityLine *> l;
	for (unsigned int i = 0; i < m_lines.size(); i++) {
		if (m_lines[i]->getType() == SingularityLine::SURFACE) l.push_back(dynamic_cast<SurfaceSingularityLine *>(m_lines[i]));
	}
	return l;
}
/*---------------------------------------------------------------------------*/
std::vector<CurveSingularityLine *>
SingularityGraph::getCurveLines()
{
	std::vector<CurveSingularityLine *> l;
	for (unsigned int i = 0; i < m_lines.size(); i++) {
		if (m_lines[i]->getType() == SingularityLine::CURVE) l.push_back(dynamic_cast<CurveSingularityLine *>(m_lines[i]));
	}
	return l;
}

/*---------------------------------------------------------------------------*/
std::vector<VolumeSingularityPoint *>
SingularityGraph::getVolumePoints()
{
	std::vector<VolumeSingularityPoint *> vol_points;
	for (unsigned int i = 0; i < m_points.size(); i++) {
		if (m_points[i]->getGeomType() == SingularityPoint::VOLUME) {
			VolumeSingularityPoint *vol_pnt = dynamic_cast<VolumeSingularityPoint *>(m_points[i]);
			vol_points.push_back(vol_pnt);
		}
	}
	return vol_points;
}
/*---------------------------------------------------------------------------*/
std::vector<VertexSingularityPoint *>
SingularityGraph::getVertexPoints()
{
	std::vector<VertexSingularityPoint *> points;
	for (unsigned int i = 0; i < m_points.size(); i++) {
		if (m_points[i]->getGeomType() == SingularityPoint::VERTEX) points.push_back(dynamic_cast<VertexSingularityPoint *>(m_points[i]));
	}
	return points;
}
/*---------------------------------------------------------------------------*/
std::vector<CurveSingularityPoint *>
SingularityGraph::getCurvePoints()
{
	std::vector<CurveSingularityPoint *> points;
	for (unsigned int i = 0; i < m_points.size(); i++) {
		if (m_points[i]->getGeomType() == SingularityPoint::CURVE) points.push_back(dynamic_cast<CurveSingularityPoint *>(m_points[i]));
	}
	return points;
}
/*---------------------------------------------------------------------------*/
std::vector<SurfaceSingularityPoint *>
SingularityGraph::getSurfacePoints()
{
	std::vector<SurfaceSingularityPoint *> points;
	for (unsigned int i = 0; i < m_points.size(); i++) {
		if (m_points[i]->getGeomType() == SingularityPoint::SURFACE) points.push_back(dynamic_cast<SurfaceSingularityPoint *>(m_points[i]));
	}
	return points;
}
/*---------------------------------------------------------------------------*/
std::vector<SingularityPatch *>
SingularityGraph::getSurfacePatchs()
{
	return m_patchs;
}
/*---------------------------------------------------------------------------*/
std::vector<SingularityLine *> &
SingularityGraph::getLines()
{
	return m_lines;
}
/*---------------------------------------------------------------------------*/
SurfaceSingularityLine *
SingularityGraph::newSurfaceLine()
{
	SurfaceSingularityLine *l = new SurfaceSingularityLine();
	m_lines.push_back(l);
	l->setNumber(getNbLines());
	return l;
}
/*---------------------------------------------------------------------------*/
SingularityPatch *
SingularityGraph::newSurfacePatch()
{
	SingularityPatch *p = new SingularityPatch();
	m_patchs.push_back(p);
	return p;
}
/*---------------------------------------------------------------------------*/
SurfaceSingularityPoint *
SingularityGraph::newSurfacePoint()
{
	SurfaceSingularityPoint *p = new SurfaceSingularityPoint(m_mesh);
	p->setNumber(getNbPoints());
	m_points.push_back(p);
	return p;
}
/*---------------------------------------------------------------------------*/
VolumeSingularityLine *
SingularityGraph::newVolumeLine()
{
	VolumeSingularityLine *l = new VolumeSingularityLine();
	m_lines.push_back(l);
	l->setNumber(getNbLines());
	return l;
}
/*---------------------------------------------------------------------------*/
VolumeSingularityPoint *
SingularityGraph::newVolumePoint()
{
	VolumeSingularityPoint *p = new VolumeSingularityPoint(m_mesh);
	p->setNumber(getNbPoints());
	m_points.push_back(p);
	return p;
}
/*---------------------------------------------------------------------------*/
CurveSingularityPoint *
SingularityGraph::newCurvePoint()
{
	CurveSingularityPoint *p = new CurveSingularityPoint(m_mesh);
	p->setNumber(getNbPoints());
	m_points.push_back(p);
	return p;
}
/*---------------------------------------------------------------------------*/
CurveSingularityLine *
SingularityGraph::newCurveLine()
{
	CurveSingularityLine *l = new CurveSingularityLine();
	m_lines.push_back(l);
	l->setNumber(getNbLines());
	return l;
}
/*---------------------------------------------------------------------------*/
VertexSingularityPoint *
SingularityGraph::newVertexPoint()
{
	VertexSingularityPoint *p = new VertexSingularityPoint(m_mesh);
	p->setNumber(getNbPoints());
	m_points.push_back(p);
	return p;
}
/*---------------------------------------------------------------------------*/
int
SingularityGraph::getNbPoints() const
{
	return m_points.size();
}
/*---------------------------------------------------------------------------*/
int
SingularityGraph::getNbSurfacePatches() const
{
	return m_patchs.size();
}
/*---------------------------------------------------------------------------*/
int
SingularityGraph::getNbLines() const
{
	return m_lines.size();
}
/*---------------------------------------------------------------------------*/
int
SingularityGraph::getNbVolumePoints() const
{
	int nbPnts = 0;

	for (unsigned int i = 0; i < m_points.size(); i++) {
		SingularityPoint *p = m_points[i];
		if (p->getGeomType() == SingularityPoint::VOLUME) {
			nbPnts++;
		}
	}
	return nbPnts;
}
/*---------------------------------------------------------------------------*/
int
SingularityGraph::getNbSurfacePoints() const
{
	int nbPnts = 0;
	for (unsigned int i = 0; i < m_points.size(); i++) {
		SingularityPoint *p = m_points[i];
		if (p->getGeomType() == SingularityPoint::SURFACE) nbPnts++;
	}
	return nbPnts;
}
/*---------------------------------------------------------------------------*/
int
SingularityGraph::getNbCurvePoints() const
{
	int nbPnts = 0;
	for (unsigned int i = 0; i < m_points.size(); i++) {
		SingularityPoint *p = m_points[i];
		if (p->getGeomType() == SingularityPoint::CURVE) nbPnts++;
	}
	return nbPnts;
}
/*---------------------------------------------------------------------------*/
int
SingularityGraph::getNbVertexPoints() const
{
	int nbPnts = 0;
	for (unsigned int i = 0; i < m_points.size(); i++) {
		SingularityPoint *p = m_points[i];
		if (p->getGeomType() == SingularityPoint::VERTEX) nbPnts++;
	}
	return nbPnts;
}
/*---------------------------------------------------------------------------*/
int
SingularityGraph::getNbVolumeLines() const
{
	int nb = 0;
	for (unsigned int i = 0; i < m_lines.size(); i++) {
		SingularityLine *l = m_lines[i];
		if (l->getType() == SingularityLine::VOLUME) nb++;
	}
	return nb;
}
/*---------------------------------------------------------------------------*/
int
SingularityGraph::getNbSurfaceLines() const
{
	int nb = 0;
	for (unsigned int i = 0; i < m_lines.size(); i++) {
		SingularityLine *l = m_lines[i];
		if (l->getType() == SingularityLine::SURFACE) nb++;
	}
	return nb;
}
/*---------------------------------------------------------------------------*/
int
SingularityGraph::getNbCurveLines() const
{
	int nb = 0;
	for (unsigned int i = 0; i < m_lines.size(); i++) {
		SingularityLine *l = m_lines[i];
		if (l->getType() == SingularityLine::CURVE) nb++;
	}
	return nb;
}
/*---------------------------------------------------------------------------*/
void
findPotentialCurveToSPlit(const std::vector<CurveSingularityLine *> &Alines, std::vector<std::pair<int, int>> &APotentialCurveToSplit, const gmds::Edge &AEdge)
{

	for (unsigned int i = 0; i < Alines.size(); i++) {
		CurveSingularityLine *l = Alines[i];
		l->healOrientation();
		std::vector<gmds::Edge> l_edges = l->getMeshEdges();
		for (unsigned int j = 0; j < l_edges.size(); j++) {
			if (l_edges[j].id() == AEdge.id()) {
				APotentialCurveToSplit.push_back(std::make_pair(i, j));
				// if the edge is not one one of the end of then line, then no other line could be candidate:
				if (j != 0 && j != l_edges.size() - 1) {
					return;
				}
			}
		}
	}
}

void
SingularityGraph::splitCurveLine(
   const gmds::math::Point &APnt, const gmds::math::Vector3d &AVec, const gmds::Edge &AEdge, SingularityPoint *&ASing, SingularityPoint::Slot *&ASlot)
{
	//==============================================================
	// STEP 1 - We find the line that must be splitted
	//==============================================================
	std::vector<CurveSingularityLine *> lines = getCurveLines();
	CurveSingularityLine *split_curve = nullptr;
	int curve_index = -1;

	std::vector<std::pair<int, int>> potentialCurveToSplit;
	findPotentialCurveToSPlit(lines, potentialCurveToSplit, AEdge);

	if (potentialCurveToSplit.empty()) {
		throw gmds::GMDSException("ERROR: No curve line found");
	}

	curve_index = potentialCurveToSplit[0].second;
	split_curve = lines[potentialCurveToSplit[0].first];

	for (int i = 0; i < potentialCurveToSplit.size() - 1; i++) {

		const auto &curve_points = split_curve->getDiscretizationPoints();
		const auto lastSegment = gmds::math::Segment(curve_points.back(), curve_points.end()[-2]);
		if (lastSegment.isIn(APnt)) break;
		const auto firstSegment = gmds::math::Segment(curve_points.front(), curve_points[1]);
		if (firstSegment.isIn(APnt)) break;

		curve_index = potentialCurveToSplit[i + 1].second;
		split_curve = lines[potentialCurveToSplit[i + 1].first];
	}

	CurveSingularityLine *newLine = newCurveLine();
	const auto &oldDiscretization = split_curve->getDiscretizationPoints();

	//==============================================================
	// STEP 2 - CREATION OF THE SINGULARITY POINT
	//==============================================================

	ASing = newCurvePoint();
	ASing->setLocation(APnt);
	ASing->addMeshEdge(AEdge);

	// create its 3 slots :
	const auto slotDir = math::Vector3d(oldDiscretization[curve_index], oldDiscretization[curve_index + 1]).normalize();

	SingularityPoint::Slot *newSlotOnSplittedCurve = ASing->newSlot(APnt, slotDir.opp(), AEdge.id(), 1, true, split_curve, slotDir);
	SingularityPoint::Slot *newSlotOnNewLine = ASing->newSlot(APnt, slotDir, AEdge.id(), 1, true, newLine, slotDir.opp());
	ASlot = ASing->newSlot(APnt, AVec.opp(), AEdge.id(), 1, true);

	//==============================================================
	// STEP 3 - CONNECTION OF LINES TO POINT
	//==============================================================
	SingularityPoint::Slot *oldSlotOnSplittedCurve = split_curve->removeSlot();
	split_curve->addSlot(newSlotOnSplittedCurve);

	// new line from last_end_point to ASing (inverted direction)
	oldSlotOnSplittedCurve->line = newLine;
	newLine->addSlot(oldSlotOnSplittedCurve);
	newLine->addSlot(newSlotOnNewLine);

	const std::vector<Edge> &init_edges = split_curve->getMeshEdges();

	std::vector<Edge> first_part, second_part;
	std::vector<gmds::math::Point> first_disc, second_disc;

	// The edge were both curves met is shared by both
	for (int i = 0; i <= curve_index; ++i) {
		first_part.push_back(init_edges[i]);
		first_disc.push_back(oldDiscretization[i]);
	}

	for (int i = init_edges.size() - 1; i >= curve_index; --i) {
		second_part.push_back(init_edges[i]);
		second_disc.push_back(oldDiscretization[i + 1]);
	}

	split_curve->setMeshEdges(first_part);
	newLine->setMeshEdges(second_part);

	first_disc.push_back(APnt);
	second_disc.push_back(APnt);
	split_curve->setDiscretizationPoints(first_disc);
	newLine->setDiscretizationPoints(second_disc);
}
/*---------------------------------------------------------------------------*/
void
SingularityGraph::splitCurveLine(
   const gmds::math::Point &APnt, const gmds::math::Vector3d &AVec, const gmds::Node &ANode, SingularityPoint *&ASing, SingularityPoint::Slot *&ASlot)
{
	//==============================================================
	// STEP 1 - We find the line that must be splitted
	//==============================================================
	std::vector<CurveSingularityLine *> lines = getCurveLines();

	bool line_found = false;
	CurveSingularityLine *split_curve = 0;

	int splitted_edge_id = -1;

	for (unsigned int i = 0; !line_found && i < lines.size(); i++) {
		CurveSingularityLine *l = lines[i];
		l->healOrientation();
		std::vector<gmds::Edge> l_edges = l->getMeshEdges();

		for (unsigned int j = 0; !line_found && j < l_edges.size(); j++) {

			std::vector<gmds::Node> nodes_j = l_edges[j].get<Node>();
			if (nodes_j[0].id() == ANode.id()) {
				line_found = true;
				splitted_edge_id = j;
				split_curve = l;
			}
			else if (nodes_j[1].id() == ANode.id()) {
				line_found = true;
				splitted_edge_id = j;
				split_curve = l;
			}
		}
	}
	if (!line_found) throw gmds::GMDSException("ERROR: No curve line found");

	//==============================================================
	// STEP 2 - CREATION OF THE SINGULARITY POINT
	//==============================================================

	ASing = newCurvePoint();
	ASing->setLocation(APnt);
	ASing->addMeshNode(ANode);

	//==============================================================
	// STEP 3 - CONNECTION OF LINES TO EDGES
	//==============================================================
	// old line keeps its first part

	SingularityPoint::Slot *last_end_slot = split_curve->removeSlot();
	split_curve->addSlot(ASlot);

	// new line from last_end_point to ASing (inverted direction)
	CurveSingularityLine *new_half_line = newCurveLine();
	new_half_line->addSlot(last_end_slot);
	new_half_line->addSlot(ASlot);

	std::vector<Edge> first_part, second_part;
	std::vector<Edge> init_edges = split_curve->getMeshEdges();

	for (unsigned int j = 0; j <= splitted_edge_id; j++) {
		first_part.push_back(init_edges[j]);
	}

	for (unsigned int j = init_edges.size() - 1; j >= splitted_edge_id; j--) {
		second_part.push_back(init_edges[j]);
	}

	split_curve->setMeshEdges(first_part);
	new_half_line->setMeshEdges(second_part);

	//==============================================================
	// STEP 4 - CONNECTION OF LINES TO POINTS
	//==============================================================
	std::vector<gmds::math::Point> new_discretization;

	// 1st curve, i.e. split_curve
	SingularityPoint *first_end = split_curve->getEndPoints()[0];
	gmds::math::Point prev_pnt = first_end->getLocation();

	bool first_end_on_vertex = (first_end->getGeomType() == SingularityPoint::VERTEX);
	new_discretization.push_back(prev_pnt);

	if (first_part.size() > 1) {
		if (!first_end_on_vertex) {
			Edge current_edge = first_part[0];
			Edge next_edge = first_part[1];

			std::vector<Node> current_nodes = current_edge.get<Node>();
			std::vector<Node> next_nodes = next_edge.get<Node>();

			Node common_node = Tools::getCommonNode(current_edge, next_edge);

			prev_pnt = common_node.point();
			new_discretization.push_back(prev_pnt);
		}
		else {     // initial sing. point on a mesh edge
			gmds::Edge current_edge = first_part[0];
			std::vector<Node> current_nodes = current_edge.get<Node>();

			gmds::math::Point p0 = current_nodes[0].point();
			gmds::math::Point p1 = current_nodes[1].point();

			double d0 = prev_pnt.distance2(p0);
			double d1 = prev_pnt.distance2(p1);

			if (d0 < d1) {
				prev_pnt = p1;
			}
			else {
				prev_pnt = p0;
			}

			new_discretization.push_back(prev_pnt);
		}

		for (unsigned int i = 1; i < first_part.size() - 1; i++) {
			Edge current_edge = first_part[i];
			std::vector<Node> current_nodes = current_edge.get<Node>();

			gmds::math::Point p0 = current_nodes[0].point();
			gmds::math::Point p1 = current_nodes[1].point();

			double d0 = prev_pnt.distance2(p0);
			double d1 = prev_pnt.distance2(p1);

			if (d0 < d1) {
				new_discretization.push_back(p1);
				prev_pnt = p1;
			}
			else {
				new_discretization.push_back(p0);
				prev_pnt = p0;
			}
		}
	}

	// we have reached the end of the firs line
	// we add the last point
	new_discretization.push_back(APnt);

	split_curve->setDiscretizationPoints(new_discretization);

	//=========================================================
	// 2nd curve, i.e. new_half_line

	new_discretization.clear();
	first_end = new_half_line->getEndPoints()[0];
	prev_pnt = first_end->getLocation();

	first_end_on_vertex = (first_end->getGeomType() == SingularityPoint::VERTEX);
	new_discretization.push_back(prev_pnt);

	if (second_part.size() > 1) {
		if (!first_end_on_vertex) {
			Edge current_edge = second_part[0];
			Edge next_edge = second_part[1];

			std::vector<Node> current_nodes = current_edge.get<Node>();
			std::vector<Node> next_nodes = next_edge.get<Node>();

			Node common_node = Tools::getCommonNode(current_edge, next_edge);

			prev_pnt = common_node.point();
			new_discretization.push_back(prev_pnt);
		}
		else {     // initial sing. point on a mesh edge
			Edge current_edge = second_part[0];
			std::vector<Node> current_nodes = current_edge.get<Node>();

			gmds::math::Point p0 = current_nodes[0].point();
			gmds::math::Point p1 = current_nodes[1].point();

			double d0 = prev_pnt.distance2(p0);
			double d1 = prev_pnt.distance2(p1);

			if (d0 < d1) {
				prev_pnt = p1;
			}
			else {
				prev_pnt = p0;
			}
			new_discretization.push_back(prev_pnt);
		}

		for (unsigned int i = 1; i < second_part.size() - 1; i++) {
			Edge current_edge = second_part[i];
			std::vector<Node> current_nodes = current_edge.get<Node>();

			gmds::math::Point p0 = current_nodes[0].point();
			gmds::math::Point p1 = current_nodes[1].point();

			double d0 = prev_pnt.distance2(p0);
			double d1 = prev_pnt.distance2(p1);

			if (d0 < d1) {
				new_discretization.push_back(p1);
				prev_pnt = p1;
			}
			else {
				new_discretization.push_back(p0);
				prev_pnt = p0;
			}
		}
	}

	// we have reached the end of the second line
	// we add the last point
	new_discretization.push_back(APnt);
	new_half_line->setDiscretizationPoints(new_discretization);

	// INVERSE CONNECTIVITY
	split_curve->healOrientation();
	new_half_line->healOrientation();

	ASing->addLine(split_curve);
	ASing->addLine(new_half_line);
	// OPEN SLOTS OF ASing
	ASing->newSlot(APnt, AVec,        // slot geometrical informaiton
	               ANode.id(), 0,     // cell owner info (id+dim)
	               true);             // it is on the surface
}
/*---------------------------------------------------------------------------*/
void
SingularityGraph::splitSurfaceLineSimple(SurfaceSingularityPoint *APnt, SurfaceSingularityLine *ALine, unsigned int APrevPointID)
{
	gmds::math::Point pnt_loc = APnt->getLocation();
	//==============================================================
	// check the intersected segment
	//==============================================================
	const std::vector<gmds::math::Point> &pnts = ALine->getDiscretizationPoints();
	const gmds::math::Segment segmentContainingSingularity(pnts[APrevPointID], pnts[APrevPointID + 1]);

	if (!segmentContainingSingularity.isIn2ndMethod(pnt_loc)) {
		cout << "problem: point " << pnt_loc << " is not found as belonging to the line " << ALine->getNumber() << endl;
		throw gmds::GMDSException("splitSurfaceLine:Error during surface line splitting");
	}

	SurfaceSingularityLine *newLine = newSurfaceLine();

	//==============================================================
	// compute discretization
	//==============================================================
	if (pnt_loc == pnts[APrevPointID]) APrevPointID--;

	std::vector<gmds::math::Point> first_part, second_part;

	first_part.insert(first_part.end(), pnts.begin(), pnts.begin() + APrevPointID + 1);
	first_part.push_back(pnt_loc);
	second_part.push_back(pnt_loc);
	second_part.insert(second_part.end(), pnts.begin() + APrevPointID + 1, pnts.end());

	ALine->setDiscretizationPoints(first_part);
	newLine->setDiscretizationPoints(second_part);

	//==============================================================
	// link lines with slots
	//==============================================================
	// old line keeps its first part
	SingularityPoint::Slot *oldLine_oldSlot = ALine->removeSlot();
	SingularityPoint::Slot *oldLine_newSlot = APnt->addLine(ALine);
	ALine->addSlot(oldLine_newSlot);

	// new line from old_end_point to ASing
	SingularityPoint::Slot *newLine_newSlot = APnt->addLine(newLine);
	newLine->addSlot(newLine_newSlot);
	newLine->addSlot(oldLine_oldSlot);
}
/*---------------------------------------------------------------------------*/
void
SingularityGraph::splitSurfaceLine(SurfaceSingularityPoint *APnt, SurfaceSingularityLine *ALine)
{
	gmds::math::Point pnt_loc = APnt->getLocation();
	//==============================================================
	// We found the interesected segment
	//==============================================================
	std::vector<gmds::math::Point> &pnts = ALine->getDiscretizationPoints();
	int index_pnt = 0;
	bool found_pnt = false;

	for (unsigned int i = 0; !found_pnt && i < pnts.size() - 1; i++) {     // 1;-2
		gmds::math::Segment sij(pnts[i], pnts[i + 1]);

		if (sij.isIn2ndMethod(pnt_loc)) {
			found_pnt = true;
			index_pnt = i;
		}
	}

	if (!found_pnt) {
		cout << "problem: point " << pnt_loc << " is not found as belonging to the line " << ALine->getNumber() << endl;
		throw gmds::GMDSException("splitSurfaceLine:Error during surface line splitting");
	}

	//==============================================================
	// old line keeps its first part
	SingularityPoint::Slot *snd_end_slot = ALine->removeSlot();
	SingularityPoint::Slot *new_snd_slot = APnt->addLine(ALine);
	ALine->addSlot(new_snd_slot);

	//==============================================================
	// new line from last_end_point to ASing (inverted direction)
	SurfaceSingularityLine *new_half_line = newSurfaceLine();
	SingularityPoint::Slot *new_slot = APnt->addLine(new_half_line);
	new_half_line->addSlot(snd_end_slot);
	new_half_line->addSlot(new_slot);

	//==============================================================
	// discretization are computed
	//==============================================================

	std::vector<gmds::math::Point> first_part, second_part;

	for (unsigned int i = 0; i < index_pnt + 1; i++) {
		first_part.push_back(pnts[i]);
	}
	first_part.push_back(pnt_loc);

	for (unsigned int i = pnts.size() - 1; i > index_pnt; i--) {
		second_part.push_back(pnts[i]);
	}
	second_part.push_back(pnt_loc);

	ALine->setDiscretizationPoints(first_part);
	new_half_line->setDiscretizationPoints(second_part);

	ALine->healOrientation();
	new_half_line->healOrientation();

	//==============================================================
	// We recompute traversed faces
	//==============================================================
	std::vector<gmds::TCellID> prev_faces = ALine->getTraversedFaces();
	std::vector<gmds::TCellID> new_faces_1, new_faces_2;

	for (unsigned int i = 0; i < prev_faces.size(); i++) {
		Face f = m_mesh->get<Face>(prev_faces[i]);
		std::vector<Node> f_nodes = f.get<Node>();
		math::Triangle t(f_nodes[0].point(), f_nodes[1].point(), f_nodes[2].point());
		bool found = false;
		for (unsigned int j = 1; !found && j < first_part.size(); j++) {
			math::Point pj = first_part[j];
			if (t.isIn(pj)) {
				found = true;
				new_faces_1.push_back(f.id());
			}
			if (j != 0) {
				math::Point pk = first_part[j - 1];
				if (t.intersect(math::Segment(pj, pk))) {
					found_pnt = true;
					new_faces_1.push_back(f.id());
				}
			}
		}
		found = false;
		for (unsigned int j = 1; !found && j < second_part.size(); j++) {
			math::Point pj = second_part[j];
			if (t.isIn(pj)) {
				found = true;
				new_faces_2.push_back(f.id());
			}
			if (j != 0) {
				math::Point pk = second_part[j - 1];
				if (t.intersect(math::Segment(pj, pk))) {
					found_pnt = true;
					new_faces_2.push_back(f.id());
				}
			}
		}
	}

	ALine->setTraversedFaces(new_faces_1);
	new_half_line->setTraversedFaces(new_faces_2);
}
/*---------------------------------------------------------------------------*/
void
SingularityGraph::removeLine(SingularityLine *ALine)
{
	if (ALine == 0) return;

	bool found_line = false;
	for (unsigned int i = 0; i < m_lines.size() && !found_line; i++) {
		SingularityLine *current_line = m_lines[i];
		if (current_line == ALine) {
			found_line = true;
			m_lines[i] = m_lines.back();
			m_lines.pop_back();
			delete ALine;
		}
	}
}
/*---------------------------------------------------------------------------*/
void
SingularityGraph::removePoint(SingularityPoint *APoint)
{
	if (APoint == 0) return;

	bool found_point = false;
	for (unsigned int i = 0; i < m_points.size() && !found_point; i++) {
		SingularityPoint *current_point = m_points[i];
		if (current_point == APoint) {
			found_point = true;
			m_points.erase(m_points.begin() + i);
		}
	}

	unsigned int i = 0;
	while (i < m_lines.size()) {
		vector<SingularityPoint *> current_points = m_lines[i]->getEndPoints();
		for (unsigned int j = 0; j < 2; j++) {
			if (current_points[j] == APoint) {
				m_lines.erase(m_lines.begin() + i);
				i--;
			}
		}
		i++;
	}
	if (found_point) delete APoint;
}
/*---------------------------------------------------------------------------*/
void
SingularityGraph::removeCurveLine(SingularityLine *ALine)
{
	// it can't be the last added line; that is the surfaceLine
	if (ALine == 0) return;

	bool found_line = false;
	for (unsigned int i = 0; i < m_lines.size() && !found_line; i++) {
		SingularityLine *current_line = m_lines[i];
		if (current_line == ALine) {
			found_line = true;
			for (unsigned int j = i; j < m_lines.size() - 1; j++)
				m_lines[j] = m_lines[j + 1];
			m_lines.pop_back();
			delete ALine;
		}
	}
}
/*---------------------------------------------------------------------------*/
void
SingularityGraph::buildSurfacePatchs()
{
	//========================================================================
	// LOOP to build patches from lines
	//========================================================================
	// First, create all patch from surface singularity lines.
	// From each surface singularity lines, two patchs are created, one following the line direction, the other in reverse.

	std::map<SingularityLine *, bool> direct, reverse;
	for (const auto line : m_lines) {
		direct[line] = false;
		reverse[line] = false;
	}

	std::vector<SingularityLine *> checkingLines;
	for (const auto initialLine : m_lines) {

		if (initialLine->getType() == SingularityLine::CURVE) continue;

		if (!direct[initialLine]) {     // First traversal
			direct[initialLine] = true;

			// A patch is created
			SingularityPatch *patch = newSurfacePatch();
			initialLine->addPatch(patch);
			std::vector<SingularityPoint *> points = initialLine->getEndPoints();

			SingularityPoint *first_point = points[0];

			patch->addPoint(first_point);
			patch->addLine(initialLine);
			SingularityPoint *current_point = points[1];
			SingularityLine *current_line = current_point->nextLine(initialLine);

			while (current_line->getNumber() != initialLine->getNumber()) {

				current_line->addPatch(patch);
				patch->addPoint(current_point);
				patch->addLine(current_line);

				points = current_line->getEndPoints();
				if (current_point == points[0]) {
					current_point = points[1];
					direct[current_line] = true;
				}
				else {
					current_point = points[0];
					reverse[current_line] = true;
				}

				current_line = current_point->nextLine(current_line);
			}

			patch->getLines(checkingLines);
			if (checkingLines.size() != 4) {
				throw GMDSException("Error during patch building, patch does not have 4 lines, line ID : " + std::to_string(current_line->getNumber()));
			}
		}
		if (!reverse[initialLine]) {     // Second traversal

			SingularityPatch *patch = newSurfacePatch();
			initialLine->addPatch(patch);
			std::vector<SingularityPoint *> points = initialLine->getEndPoints();

			SingularityPoint *first_point = points[1];

			patch->addPoint(first_point);
			patch->addLine(initialLine);
			SingularityPoint *current_point = points[0];
			SingularityLine *current_line = current_point->nextLine(initialLine);

			while (current_line->getNumber() != initialLine->getNumber()) {

				current_line->addPatch(patch);
				patch->addPoint(current_point);
				patch->addLine(current_line);

				points = current_line->getEndPoints();
				if (current_point == points[0]) {
					current_point = points[1];
					direct[current_line] = true;
				}
				else {
					current_point = points[0];
					reverse[current_line] = true;
				}

				current_line = current_point->nextLine(current_line);
			}
			patch->getLines(checkingLines);
			if (checkingLines.size() != 4) {
				throw GMDSException("Error during patch building, patch does not have 4 lines, line ID : " + std::to_string(current_line->getNumber()));
			}
		}
	}
	// At this point, isolated parts of the mesh containing no surface singularity lines may remain.
	// This time, the corresponding patch are build from the boundary lines.
	for (const auto initialLine : m_lines) {
		if (!direct[initialLine] && !reverse[initialLine]) {
			direct[initialLine] = true;

			SingularityPatch *patch = newSurfacePatch();
			initialLine->addPatch(patch);

			auto singPoints = initialLine->getEndPoints();
			const auto firstPoint = singPoints.front();
			auto nextPoint = singPoints.back();

			patch->addLine(initialLine);
			patch->addPoint(nextPoint);

			auto nextLine = initialLine;
			for (size_t i = 0; i < 3; ++i) {
				const auto slots = nextPoint->getSlots();
				if (slots.size() != 2) {
					throw GMDSException("Error during patch building, bdry Singularity have more than 2 slots, sing ID : " + std::to_string(nextPoint->getNumber()));
				}

				nextLine = slots.front()->line == nextLine ? slots.back()->line : slots.front()->line;
				singPoints = nextLine->getEndPoints();
				nextPoint = singPoints.front() == nextPoint ? singPoints.back() : singPoints.front();

				nextLine->addPatch(patch);
				patch->addLine(nextLine);
				patch->addPoint(nextPoint);
				direct[nextLine] = true;
			}
			if (nextPoint != firstPoint) {
				throw GMDSException("Error during patch building, patch does not have 4 lines, line ID : " + std::to_string(nextLine->getNumber()));
			}
		}
	}
}

/*---------------------------------------------------------------------------*/
void
SingularityGraph::buildCurveSurfacePatchs(unsigned int &number_of_control_points)
{
	//========================================================================
	// We start from a set of lines and points that must be connected
	// in an acceptable configuration. Wrong configurations will not
	// be detected leading to a wrong patch decomposition.
	//========================================================================

	cout << "buildCurveSurfacePatchs" << endl;
	// For each point, we remove the connected lines
	vector<SingularityPoint *> pi = m_points;
	for (unsigned int i = 0; i < m_points.size(); i++) {
		// SingularityPoint *pi = m_points[i];
		pi[i]->clearSlots();
	}

	// For each line, we remove connected points.
	const vector<SingularityLine *> &li = m_lines;
	for (unsigned int i = 0; i < li.size(); i++) {
		if (li[i]->getEndPoints().size() != 2) {
			throw GMDSException("Error during topoloy rebuilding: a line is not adjacent to 2 points");
		}
	}
	//========================================================================
	// TOPOLOGY REBUILDING
	//========================================================================
	// Now we reconnnect to have a clean configuration. Each curve must have
	// two end points.
	// WARNING if a line is close to an interior boundary(problematic concave bdry), simple Bezier could compute the line as crossing that bdry;
	// solution: if any of the lines points are on the boundary, choose a number of control points higher than the number of points the line originally has
	// also: start with patches on the boundary
	gmds::math::Point zeroPoint(0.0, 0.0, 0.0);
	int maxLineNo = 0;
	for (unsigned int i = 0; i < li.size(); i++) {
		maxLineNo = std::max(li[i]->getNumber(), maxLineNo);
	}

	vector<bool> refinedLines(maxLineNo, false);
	for (unsigned int i = 0; i < li.size(); i++) {

		std::vector<SingularityPoint *> points = li[i]->getEndPoints();

		SingularityPoint *current_point = points[1];
		SingularityLine *current_line = current_point->nextLine(li[i]);
		bool boundaryPatch = false;

		while (current_line->getNumber() != li[i]->getNumber()) {

			points = current_line->getEndPoints();

			if (current_line->getType() == SingularityLine::CURVE) {
				boundaryPatch = true;
				break;
			}

			if (current_point == points[0]) {
				current_point = points[1];
			}
			else {
				current_point = points[0];
			}

			current_line = current_point->nextLine(current_line);
		}
		if (boundaryPatch) {
			// lines2Refine[cont] = i;
			std::vector<SingularityPoint *> points = li[i]->getEndPoints();
			SingularityPoint *current_point = points[1];
			SingularityLine *current_line = current_point->nextLine(li[i]);

			while (current_line->getNumber() != li[i]->getNumber()) {

				points = current_line->getEndPoints();
				std::vector<math::Point> li_points = li[i]->getDiscretizationPoints();
				// technically here we have to create the Bezier curve;
				// either define manually level of detail and proceed as before, either need a different format in SingularityLine
				math::Point p1 = li_points[0];
				math::Point p2 = li_points.back();

				if (!refinedLines[li[i]->getNumber() - 1]) {
					if (li[i]->getType() != 0) {
						// interior line
						unsigned int number_of_original_points = li_points.size();
						unsigned int local_number_of_control_points = number_of_original_points;     // number_of_control_points;
						local_number_of_control_points = std::min(local_number_of_control_points, number_of_original_points);
						unsigned int number_of_steps = 2 * local_number_of_control_points;        //(precision) this should be set
						std::vector<math::Point> new_line_points(number_of_steps, zeroPoint);     // new_line_points(li_points);//
						// new_line_points[0] = li_points[0];
						for (unsigned int j = 0; j < number_of_steps; j++) {
							double u = (double) j / (double) (number_of_steps - 1);
							double ratioPoints = (double) (number_of_original_points - 1) / (local_number_of_control_points - 1);
							double bm = BernsteinBasis(0, local_number_of_control_points - 1, u);
							new_line_points[j].X() = new_line_points[j].X() + bm * li_points[0].X();
							new_line_points[j].Y() = new_line_points[j].Y() + bm * li_points[0].Y();
							new_line_points[j].Z() = new_line_points[j].Z() + bm * li_points[0].Z();
							bm = BernsteinBasis(local_number_of_control_points - 1, local_number_of_control_points - 1, u);

							new_line_points[j].X() = new_line_points[j].X() + bm * li_points.back().X();
							new_line_points[j].Y() = new_line_points[j].Y() + bm * li_points.back().Y();
							new_line_points[j].Z() = new_line_points[j].Z() + bm * li_points.back().Z();

							for (unsigned int k = 1; k < local_number_of_control_points - 1; k++) {

								double bm = BernsteinBasis(k, local_number_of_control_points - 1, u);

								math::Point tempPoint;     // = li_points[k];
								if (k * ratioPoints == (unsigned int) (k * ratioPoints)) {
									tempPoint = li_points[(unsigned int) (k * ratioPoints)];
								}
								else {

									double interpFact = (double) k * ratioPoints - (unsigned int) (k * ratioPoints);

									tempPoint =
									   interpolateBtwPoints(li_points[(unsigned int) (k * ratioPoints)], li_points[(unsigned int) (k * ratioPoints) + 1], interpFact);
								}

								new_line_points[j].X() = new_line_points[j].X() + bm * tempPoint.X();
								new_line_points[j].Y() = new_line_points[j].Y() + bm * tempPoint.Y();
								new_line_points[j].Z() = new_line_points[j].Z() + bm * tempPoint.Z();
							}
						}

						li[i]->setDiscretizationPoints(new_line_points);
						refinedLines[li[i]->getNumber() - 1] = true;
					}
					// else bdry line is left as is
					/*double nb_connections=0;
					li[i]->removeAllSingularityPoints();
					for(unsigned int j=0; j<pi.size(); j++) {

					   if(math::near(p1.distance(pi[j]->getLocation()),0.0) ||
					      math::near(p2.distance(pi[j]->getLocation()),0.0) ) {

					         pi[j]->addLine(li[i]);
					      li[i]->addSingularityPoint(pi[j]);
					      nb_connections++;
					      }
					}
					if(nb_connections!=2){
					   throw GMDSException("Error during topoloy rebuilding: a line is not adjacent to 2 points");
					}*/
				}
				if (current_point == points[0]) {
					current_point = points[1];
				}
				else {
					current_point = points[0];
				}

				current_line = current_point->nextLine(current_line);
			}
		}
	}
	for (unsigned int i = 0; i < li.size(); i++) {

		if (!refinedLines[li[i]->getNumber() - 1]) {

			std::vector<math::Point> li_points = li[i]->getDiscretizationPoints();
			// technically here we have to create the Bezier curve;
			// either define manually level of detail and proceed as before, either need a different format in SingularityLine
			math::Point p1 = li_points[0];
			math::Point p2 = li_points.back();

			if (li[i]->getType() != 0) {
				// interior line
				unsigned int number_of_original_points = li_points.size();
				unsigned int local_number_of_control_points = number_of_control_points;
				local_number_of_control_points = std::min(local_number_of_control_points, number_of_original_points);
				unsigned int number_of_steps = 2 * local_number_of_control_points;        //(precision) this should be set
				std::vector<math::Point> new_line_points(number_of_steps, zeroPoint);     // new_line_points(li_points);//
				// new_line_points[0] = li_points[0];
				for (unsigned int j = 0; j < number_of_steps; j++) {
					double u = (double) j / (double) (number_of_steps - 1);
					double ratioPoints = (double) (number_of_original_points - 1) / (local_number_of_control_points - 1);
					double bm = BernsteinBasis(0, local_number_of_control_points - 1, u);
					new_line_points[j].X() = new_line_points[j].X() + bm * li_points[0].X();
					new_line_points[j].Y() = new_line_points[j].Y() + bm * li_points[0].Y();
					new_line_points[j].Z() = new_line_points[j].Z() + bm * li_points[0].Z();
					bm = BernsteinBasis(local_number_of_control_points - 1, local_number_of_control_points - 1, u);

					new_line_points[j].X() = new_line_points[j].X() + bm * li_points.back().X();
					new_line_points[j].Y() = new_line_points[j].Y() + bm * li_points.back().Y();
					new_line_points[j].Z() = new_line_points[j].Z() + bm * li_points.back().Z();

					for (unsigned int k = 1; k < local_number_of_control_points - 1; k++) {

						double bm = BernsteinBasis(k, local_number_of_control_points - 1, u);

						math::Point tempPoint;     // = li_points[k];
						if (k * ratioPoints == (unsigned int) (k * ratioPoints)) {
							tempPoint = li_points[(unsigned int) (k * ratioPoints)];
						}
						else {

							double interpFact = (double) k * ratioPoints - (unsigned int) (k * ratioPoints);

							tempPoint = interpolateBtwPoints(li_points[(unsigned int) (k * ratioPoints)], li_points[(unsigned int) (k * ratioPoints) + 1], interpFact);
						}

						new_line_points[j].X() = new_line_points[j].X() + bm * tempPoint.X();
						new_line_points[j].Y() = new_line_points[j].Y() + bm * tempPoint.Y();
						new_line_points[j].Z() = new_line_points[j].Z() + bm * tempPoint.Z();
					}
					// cout<<"new_line_points["<<j<<"] "<<new_line_points[j]<<endl;
				}

				li[i]->setDiscretizationPoints(new_line_points);
			}
			else {
				// simple interpolation for bdry or leave as it is
			}
			/*double nb_connections=0;

			for(unsigned int j=0; j<pi.size(); j++) {

			      if(math::near(p1.distance(pi[j]->getLocation()),0.0) ||
			         math::near(p2.distance(pi[j]->getLocation()),0.0) ) {

			         pi[j]->addLine(li[i]);
			         li[i]->addSingularityPoint(pi[j]);
			         //cout<<"ads sing point "<<pi[j]->getLocation()<<endl;
			         nb_connections++;
			      }
			}
			if(nb_connections!=2){
			   throw GMDSException("Error during topoloy rebuilding: a line is not adjacent to 2 points");
			}*/

			refinedLines[li[i]->getNumber() - 1] = true;
		}
	}

	// We store an integer for each curve to indicate how many times it has
	// been traversed. It must be 1 (on the boundary) or 2 (inside). To detect
	// if a curve is on the boundary, we check if it classified onto a curve.
	std::map<SingularityLine *, bool> direct, reverse;
	for (unsigned int i = 0; i < li.size(); i++) {
		direct[li[i]] = false;
		reverse[li[i]] = false;
	}
	//========================================================================
	// LOOP to build patches from lines
	//========================================================================
	// for(unsigned int i=0;i<li.size(); i++) {
	unsigned int i = 0;
	while (i < li.size()) {
		if (li[i]->getType() == SingularityLine::CURVE) {
			i++;
			continue;
		}
		if (!direct[li[i]]) {     // first time at least one patch to build
			// We never work with curve line as first
			//=====================================================
			// First traversal
			//=====================================================
			direct[li[i]] = true;

			// A patch is created
			SingularityPatch *patch = newSurfacePatch();

			std::vector<SingularityPoint *> points = li[i]->getEndPoints();

			SingularityPoint *first_point = points[0];

			patch->addPoint(first_point);
			// cout<<"patch->addPoint "<<first_point->getLocation()<<endl;
			patch->addLine(li[i]);

			SingularityPoint *current_point = points[1];

			SingularityLine *current_line = current_point->nextLine(li[i]);

			while (current_line->getNumber() != li[i]->getNumber()) {
				patch->addPoint(current_point);
				patch->addLine(current_line);
				// cout<<"patch->addPoint "<<current_point->getLocation()<<endl;

				points = current_line->getEndPoints();

				if (current_point == points[0]) {
					current_point = points[1];
					direct[current_line] = true;
				}
				else {
					current_point = points[0];
					reverse[current_line] = true;
				}

				current_line = current_point->nextLine(current_line);
			}

			// we decrement the index to traverse it again
			i--;
		}
		else if (!reverse[li[i]]) {
			//=====================================================
			// Second traversal
			//=====================================================
			// second time, a patch can be built if the curve is not on the boundary.
			// Second time we go from the 2nd end point towards the 1st one A patch is
			// created

			// cout<<"!reverse[li[i]] "<<endl;
			// cout<<"A patch is created "<<endl;
			SingularityPatch *patch = newSurfacePatch();

			std::vector<SingularityPoint *> points = li[i]->getEndPoints();

			SingularityPoint *first_point = points[1];

			patch->addPoint(first_point);
			patch->addLine(li[i]);
			// cout<<"patch->addPoint "<<first_point->getLocation()<<endl;

			SingularityPoint *current_point = points[0];
			SingularityLine *current_line = current_point->nextLine(li[i]);

			while (current_line->getNumber() != li[i]->getNumber()) {

				patch->addPoint(current_point);
				patch->addLine(current_line);
				// cout<<"patch->addPoint "<<current_point->getLocation()<<endl;
				points = current_line->getEndPoints();
				// cout<<"points: "<<points[0]->getLocation()<<" , "<<points[1]->getLocation()<<endl;
				if (current_point == points[0]) {
					current_point = points[1];
					direct[current_line] = true;
				}
				else {
					current_point = points[0];
					reverse[current_line] = true;
				}

				current_line = current_point->nextLine(current_line);
			}
		}
		i++;
	}
}

/*---------------------------------------------------------------------------*/
double
SingularityGraph::BernsteinBasis(unsigned int i, unsigned int n, double u)
{
	unsigned int ni = n - i;
	double r = 1.0;
	if ((u == 0.0) && (i == 0.0))     // pow(u,i)
		return r;

	if ((u == 1.0) && (i == n))     // pow(1-u,n-i);
		return r;

	r = binomialCoeff(n, i);
	// cout<<"r= "<<r<<" for (n="<<n<<" ,i="<<i<<") and u= "<<u<<endl;

	r = r * pow(u, i);
	// cout<<"r= "<<r<<endl;
	r = r * pow(1.0 - u, n - i);
	// cout<<"r= "<<r<<endl;
	return r;
}
/*---------------------------------------------------------------------------*/
double
SingularityGraph::binomialCoeff(unsigned int &n, unsigned int &k)
{
	double res = 1.0;
	unsigned int copyk = k;

	// Since C(n, k) = C(n, n-k)
	if (copyk > n - copyk) copyk = n - copyk;

	for (unsigned int i = 0; i < copyk; i++) {
		res = res * (double) (n - i);
		res = res / (double) (i + 1);
	}

	return res;
}

/*---------------------------------------------------------------------------*/
void
SingularityGraph::createVTKOutputFile(const std::string &AFileName) const
{
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
	if (m_patchs.empty()) {

		gmds::Variable<int> *var = nullptr;
		if (m.hasVariable(GMDS_FACE, "lines"))
			var = m.getVariable<int, GMDS_FACE>("lines");
		else
			var = m.newVariable<int, GMDS_FACE>("lines");

		for (unsigned int i = 0; i < m_lines.size(); i++) {
			SingularityLine *current_line = m_lines[i];
			std::vector<gmds::math::Point> &points = current_line->getDiscretizationPoints();
			for (unsigned int j = 0; j < points.size() - 1; j++) {
				gmds::Node n1 = m.newNode(points[j].X(), points[j].Y(), points[j].Z());
				gmds::Node n2 = m.newNode(points[j + 1].X(), points[j + 1].Y(), points[j + 1].Z());
				gmds::Face f = m.newTriangle(n1, n1, n2);
				(*var)[f.id()] = current_line->getNumber();
			}
		}
	}
	else {

		for (unsigned int i = 0; i < m_patchs.size(); i++) {

			SingularityPatch *current_patch = m_patchs[i];

			std::vector<SingularityPoint *> current_points;
			current_patch->getPoints(current_points);

			std::vector<Node> current_nodes;
			for (unsigned int j = 0; j < current_points.size(); j++) {
				gmds::math::Point p = current_points[j]->getLocation();
				current_nodes.push_back(m.newNode(p.X(), p.Y(), p.Z()));
			}
			m.newFace(current_nodes);
		}
	}

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N | gmds::F);
	vtkWriter.setDataOptions(gmds::N | gmds::F);
	vtkWriter.write(AFileName);
}
/*---------------------------------------------------------------------------*/
void
SingularityGraph::createVTKOutputFile(const std::string &AFileName, bool &curvePatches) const
{
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
	if (m_patchs.empty()) {
		gmds::Variable<int> *var = 0;
		try {
			var = m.getVariable<int, GMDS_FACE>("lines");
		}
		catch (gmds::GMDSException &) {
			var = m.newVariable<int, GMDS_FACE>("lines");
		}

		for (unsigned int i = 0; i < m_lines.size(); i++) {
			cout << "i= " << i << endl;
			SingularityLine *current_line = m_lines[i];
			std::vector<gmds::math::Point> &points = current_line->getDiscretizationPoints();
			for (unsigned int j = 0; j < points.size() - 1; j++) {
				gmds::Node n1 = m.newNode(points[j].X(), points[j].Y(), points[j].Z());
				gmds::Node n2 = m.newNode(points[j + 1].X(), points[j + 1].Y(), points[j + 1].Z());
				gmds::Face f = m.newTriangle(n1, n1, n2);
				(*var)[f.id()] = current_line->getNumber();
			}
		}
	}
	else {
		if (!curvePatches) {
			for (unsigned int i = 0; i < m_patchs.size(); i++) {

				SingularityPatch *current_patch = m_patchs[i];

				std::vector<SingularityPoint *> current_points;
				current_patch->getPoints(current_points);

				std::vector<Node> current_nodes;
				for (unsigned int j = 0; j < current_points.size(); j++) {
					gmds::math::Point p = current_points[j]->getLocation();
					current_nodes.push_back(m.newNode(p.X(), p.Y(), p.Z()));
				}
				m.newFace(current_nodes);
			}
		}
		else {
			// curvePatches
			for (unsigned int i = 0; i < m_patchs.size(); i++) {
				gmds::Mesh patchMesh(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
				SingularityPatch *current_patch = m_patchs[i];
				std::vector<SingularityPoint *> current_points;
				current_patch->getPoints(current_points);

				std::vector<Node> current_nodes;
				for (unsigned int j = 0; j < current_points.size(); j++) {
					gmds::math::Point p = current_points[j]->getLocation();
					current_nodes.push_back(m.newNode(p.X(), p.Y(), p.Z()));
				}
				m.newFace(current_nodes);

				std::vector<SingularityLine *> current_lines;
				current_patch->getLines(current_lines);

				for (unsigned int j = 0; j < current_lines.size(); j++) {
					std::vector<gmds::math::Point> line_discretization_points = current_lines[j]->getDiscretizationPoints();
					math::Point APoint = line_discretization_points[0];
					for (unsigned int k = 1; k < line_discretization_points.size(); k++) {
						math::Point prevPoint = APoint;
						APoint = line_discretization_points[k];
						gmds::Node node0 = patchMesh.newNode(prevPoint.X(), prevPoint.Y(), prevPoint.Z());
						gmds::Node node1 = patchMesh.newNode(APoint.X(), APoint.Y(), APoint.Z());
						gmds::Face f = patchMesh.newTriangle(node0, node1, node1);
					}
				}
				gmds::IGMeshIOService ioService(&patchMesh);
				gmds::VTKWriter vtkWriter(&ioService);
				vtkWriter.setCellOptions(gmds::N | gmds::F);
				vtkWriter.setDataOptions(gmds::N | gmds::F);
				std::string file_name = AFileName + "_patch_" + std::to_string(i) + ".vtk";
				vtkWriter.write(file_name);
			}
		}
	}

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N | gmds::F);
	vtkWriter.setDataOptions(gmds::N | gmds::F);
	vtkWriter.write(AFileName);
}

void
SingularityGraph::computeLinkedEdges()
{
	std::list<SingularityLine *> boundaryLines;
	std::list<SingularityLine *> surfaceLines;
	for (const auto line : getLines()) {
		if (line->getType() == SingularityLine::ESingularityGeomLineType::CURVE)
			boundaryLines.push_back(line);
		else
			surfaceLines.push_back(line);
	}
	int groupID = 0;
	while (!boundaryLines.empty()) {
		SingularityLine *line = boundaryLines.front();
		boundaryLines.pop_front();

		line->setLinkedEdgeID(groupID);
		auto prevPatch = line->getPatchs()[0];
		SingularityLine *nextline = prevPatch->computeOpposedLine(line);
		nextline->setLinkedEdgeID(groupID);
		while (nextline->getType() != SingularityLine::CURVE) {
			surfaceLines.remove(nextline);
			for (const auto patch : nextline->getPatchs()) {
				if (patch != prevPatch) {
					nextline = patch->computeOpposedLine(nextline);
					nextline->setLinkedEdgeID(groupID);
					prevPatch = patch;
					break;
				}
			}
		}
		boundaryLines.remove(nextline);
		groupID++;
	}
	while (!surfaceLines.empty()) {     // find cycling edge group now
		SingularityLine *line = surfaceLines.front();
		surfaceLines.pop_front();

		line->setLinkedEdgeID(groupID);
		auto prevPatch = line->getPatchs()[0];
		SingularityLine *nextline = prevPatch->computeOpposedLine(line);
		nextline->setLinkedEdgeID(groupID);
		while (nextline != line) {
			surfaceLines.remove(nextline);
			for (const auto patch : nextline->getPatchs()) {
				if (patch != prevPatch) {
					nextline = patch->computeOpposedLine(nextline);
					nextline->setLinkedEdgeID(groupID);
					prevPatch = patch;
					break;
				}
			}
		}
		groupID++;
	}
	m_nGroupLinkedEdges = groupID;
}

/*---------------------------------------------------------------------------*/
gmds::math::Point
SingularityGraph::interpolateBtwPoints(gmds::math::Point &point1, gmds::math::Point &point2, double &interpFact)
{
	gmds::math::Point interpolatedPoint;

	interpolatedPoint.X() = point1.X() * interpFact + point2.X() * (1.0 - interpFact);
	interpolatedPoint.Y() = point1.Y() * interpFact + point2.Y() * (1.0 - interpFact);
	interpolatedPoint.Z() = point1.Z() * interpFact + point2.Z() * (1.0 - interpFact);

	return interpolatedPoint;
}
