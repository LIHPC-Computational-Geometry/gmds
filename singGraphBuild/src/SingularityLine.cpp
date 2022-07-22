/*----------------------------------------------------------------------------*/
/*
 * SingularityLine.cpp
 *
 *  Created on: 13 juil. 2014
 *      Author: bibi
 */
/*----------------------------------------------------------------------------*/
#include <gmds/math/Numerics.h>
#include <gmds/math/Segment.h>
#include <gmds/singGraphBuild/SingularityLine.h>
using namespace std;
/*----------------------------------------------------------------------------*/
SingularityLine::SingularityLine(const ESingularityGeomLineType AType) : m_type(AType) {}
/*----------------------------------------------------------------------------*/
SingularityLine::~SingularityLine() {}
/*----------------------------------------------------------------------------*/
void
SingularityLine::addSlot(SingularityPoint::Slot *ASlot)
{
	if (m_slots.size() >= 2) {
		throw gmds::GMDSException("SingularityLine::A singularity line can only have two end points");
	}
	ASlot->line = this;
	m_slots.push_back(ASlot);
}
void
SingularityLine::addSingularityPoint(SingularityPoint *ASing)
{
	throw gmds::GMDSException("Unimplemented method");
}
/*----------------------------------------------------------------------------*/
void
SingularityLine::addDiscretizationPoint(const gmds::math::Point &APnt)
{
	m_discretization_points.push_back(APnt);
}
/*----------------------------------------------------------------------------*/
void
SingularityLine::setDiscretizationPoints(const std::vector<gmds::math::Point> &APoints)
{
	// m_discretization_points.clear();
	m_discretization_points = APoints;
}
/*----------------------------------------------------------------------------*/
std::vector<gmds::math::Point> &
SingularityLine::getDiscretizationPoints()
{
	return m_discretization_points;
}
/*----------------------------------------------------------------------------*/
std::vector<SingularityPoint *>
SingularityLine::getEndPoints() const
{
	return {m_slots[0]->from_point, m_slots[1]->from_point};
}
std::vector<SingularityPoint::Slot *> &
SingularityLine::getSlots()
{
	return m_slots;
}
/*----------------------------------------------------------------------------*/
SingularityPoint::Slot *
SingularityLine::getSlotAtLocation(const unsigned int &APointID)
{
	const gmds::math::Point &APnt = m_discretization_points[APointID];
	for (const auto slot : getSlots())
		if (gmds::math::near(APnt.distance2(slot->location), 0.0)) return slot;
	return nullptr;
}
/*----------------------------------------------------------------------------*/
bool
SingularityLine::healOrientation()
{ /*WARNING  this should work in general, for the examples given;
	however special cases can appear(especially in 3d); I should make a healOrientatinByConnectivity   */
	// We need two end points
	if (m_slots.size() != 2) {
		return false;
	}
	return false;
	gmds::math::Point first_end = m_slots[0]->from_point->getLocation();
	gmds::math::Point second_end = m_slots[1]->from_point->getLocation();

	gmds::math::Point first_discretization_end = m_discretization_points[0];
	gmds::math::Point second_discretization_end = m_discretization_points[1];

	double d1 = first_end.distance(first_discretization_end);
	double d2 = first_end.distance(second_discretization_end);
	double d3 = second_end.distance(first_discretization_end);
	double d4 = second_end.distance(second_discretization_end);

	// We need that discretization points link end points
	if (!gmds::math::near(d1, 0.0) && !gmds::math::near(d2, 0.0) && !gmds::math::near(d3, 0.0) && !gmds::math::near(d4, 0.0)) return false;

	if (!gmds::math::near(d1, 0.0)) {
		// we have to inverse the line discretization
		std::vector<gmds::math::Point> tmp_points = m_discretization_points;

		unsigned int nb_elts = m_discretization_points.size();
		for (unsigned int i = 0; i < nb_elts; i++)
			m_discretization_points[i] = tmp_points[nb_elts - 1 - i];
	}

	return true;
}
/*----------------------------------------------------------------------------*/
bool
CurveSingularityLine::healOrientation()
{
	return true;
	bool parent_result = SingularityLine::healOrientation();

	gmds::math::Point p1 = m_ordered_mesh_edges[0].get<gmds::Node>()[0].point();
	gmds::math::Point p2 = m_ordered_mesh_edges[m_ordered_mesh_edges.size() - 1].get<gmds::Node>()[0].point();
	gmds::math::Point first_end = m_slots[0]->from_point->getLocation();
	double d1 = first_end.distance(p1);
	double d2 = first_end.distance(p2);
	if (d1 > d2) {
		// we have to inverse the line discretization
		std::vector<gmds::Edge> tmp_edges = m_ordered_mesh_edges;

		unsigned int nb_elts = m_ordered_mesh_edges.size();
		for (unsigned int i = 0; i < nb_elts; i++) {
			m_ordered_mesh_edges[i] = tmp_edges[nb_elts - 1 - i];
		}
	}

	return parent_result;
}
/*----------------------------------------------------------------------------*/
void
SurfaceSingularityLine::smooth(const int ANbStep)
{
	std::cout << "Smoothing is not implemented for surface singularity lines!!" << std::endl;
}
/*----------------------------------------------------------------------------*/
void
SurfaceSingularityLine::addTraversedFace(const gmds::TCellID AID)
{
	m_traversed_faces_id.push_back(AID);
}
/*----------------------------------------------------------------------------*/
std::vector<gmds::TCellID>
SurfaceSingularityLine::getTraversedFaces() const
{
	return m_traversed_faces_id;
}
/*----------------------------------------------------------------------------*/
bool
SurfaceSingularityLine::isTraversed(const gmds::TCellID AFaceID) const
{
	for (unsigned int i = 0; i < m_traversed_faces_id.size(); i++) {
		if (AFaceID == m_traversed_faces_id[i]) return true;
	}

	return false;
}
/*----------------------------------------------------------------------------*/
gmds::math::Vector3d
SingularityLine::getTangent(const double AParam) const
{
	if (AParam < 0 || AParam > 1) throw gmds::GMDSException("The line parameter must be in [0,1]");

	if (m_discretization_points.empty()) throw gmds::GMDSException("Line discretization must be known to compute tangent");

	gmds::math::Vector3d t;
	if (AParam == 0) {
		gmds::math::Vector3d v1=m_discretization_points[1]- m_discretization_points[0];
		gmds::math::Vector3d v2=m_discretization_points[2]- m_discretization_points[1];
		t = 0.5 * (v1 + v2);
	}
	else if (AParam == 1) {
		gmds::math::Vector3d v1=m_discretization_points[m_discretization_points.size() - 1]- m_discretization_points[m_discretization_points.size() - 2];
		gmds::math::Vector3d v2=m_discretization_points[m_discretization_points.size() - 2]- m_discretization_points[m_discretization_points.size() - 3];
		t = 0.5 * (v1 + v2);
	}
	else {
		double l = AParam * length();
		double c = 0;
		int i = 0;
		while (c < l) {
			gmds::math::Point pi = m_discretization_points[i];
			gmds::math::Point pj = m_discretization_points[i + 1];
			double dij = pi.distance(pj);
			if (c + dij >= l) {
				// stop on this edge or on pj
				gmds::math::Vector3d vij=pj-pi;
				// we smooth the tangent with previous and next segment if we can
				if (i == 0) {
					// only next segment
					gmds::math::Point pk = m_discretization_points[i + 2];
					gmds::math::Vector3d vjk=pj-pi;
					return 0.5 * (vij + vjk);
				}
				else if (i == m_discretization_points.size() - 2) {
					// only prev segment
					gmds::math::Point ph = m_discretization_points[i - 1];
					gmds::math::Vector3d vhi=ph-pi;
					return 0.5 * (vij + vhi);
				}
				else {
					// prev and next segments
					gmds::math::Point ph = m_discretization_points[i - 1];
					gmds::math::Point pk = m_discretization_points[i + 2];
					gmds::math::Vector3d vhi=pi-ph;
					gmds::math::Vector3d vjk=pi-pj;
					double r = 1.0 / 3.0;
					gmds::math::Vector3d v = r * (vhi + vij + vjk);
					v.normalize();
					return v;
				}
			}
			else {
				// execute to the next edge
				c += dij;
				i++;
			}

		}     // while (c<l)

	}     // else

	t.normalize();
	return t;
}

/*----------------------------------------------------------------------------*/
gmds::math::Point
SingularityLine::getPoint(const double AParam) const
{
	if (AParam < 0 || AParam > 1) throw gmds::GMDSException("The line parameter must be in [0,1]");

	if (m_discretization_points.empty()) throw gmds::GMDSException("Line discretization must be known to compute tangent");

	gmds::math::Point t;
	if (AParam == 0)
		t = m_discretization_points[0];
	else if (AParam == 1)
		t = m_discretization_points[m_discretization_points.size() - 1];
	else {
		double l = AParam * length();
		double c = 0;
		int i = 0;
		while (c < l) {
			gmds::math::Point pi = m_discretization_points[i];
			gmds::math::Point pj = m_discretization_points[i + 1];
			double dij = pi.distance(pj);
			if (c + dij > l) {
				// stop on this edge
				auto p = (l - c) / dij;
				return (1 - p) * pi + p * pj;
			}
			else if (c + dij == l) {
				// stop on pj
				return pj;
			}
			else {
				// execute to the next edge
				c += dij;
				i++;
			}

		}     // while (c<l)

	}     // else

	return t;
}
/*----------------------------------------------------------------------------*/
void
VolumeSingularityLine::smooth(const int ANbStep)
{
	std::vector<gmds::math::Point> new_discretization;
	new_discretization.resize(m_discretization_points.size());
	if (m_discretization_points.size() < 5) return;

	for (int k = 0; k < ANbStep; k++) {
		for (unsigned int i = 0; i < m_discretization_points.size(); i++) {
			if (i == 0 || i == m_discretization_points.size() - 1) {
				new_discretization[i] = m_discretization_points[i];
			}
			else if (i == 1) {
				gmds::math::Point p0 = m_discretization_points[0];
				gmds::math::Point p1 = m_discretization_points[1];
				gmds::math::Point p2 = m_discretization_points[2];
				gmds::math::Point p3 = m_discretization_points[3];
				gmds::math::Point p4 = m_discretization_points[4];
				new_discretization[i] = gmds::math::Point(0.309 * p0.X() + 0.382 * p1.X() + 0.242 * p2.X() + 0.061 * p3.X() + 0.006 * p4.X(),
				                                          0.309 * p0.Y() + 0.382 * p1.Y() + 0.242 * p2.Y() + 0.061 * p3.Y() + 0.006 * p4.Y(),
				                                          0.309 * p0.Z() + 0.382 * p1.Z() + 0.242 * p2.Z() + 0.061 * p3.Z() + 0.006 * p4.Z());
			}
			else if (i == m_discretization_points.size() - 2) {
				gmds::math::Point p0 = m_discretization_points[m_discretization_points.size() - 1];
				gmds::math::Point p1 = m_discretization_points[m_discretization_points.size() - 2];
				gmds::math::Point p2 = m_discretization_points[m_discretization_points.size() - 3];
				gmds::math::Point p3 = m_discretization_points[m_discretization_points.size() - 4];
				gmds::math::Point p4 = m_discretization_points[m_discretization_points.size() - 5];
				new_discretization[i] = gmds::math::Point(0.309 * p0.X() + 0.382 * p1.X() + 0.242 * p2.X() + 0.061 * p3.X() + 0.006 * p4.X(),
				                                          0.309 * p0.Y() + 0.382 * p1.Y() + 0.242 * p2.Y() + 0.061 * p3.Y() + 0.006 * p4.Y(),
				                                          0.309 * p0.Z() + 0.382 * p1.Z() + 0.242 * p2.Z() + 0.061 * p3.Z() + 0.006 * p4.Z());
			}
			else if (i == 2) {
				gmds::math::Point p0 = m_discretization_points[0];
				gmds::math::Point p1 = m_discretization_points[1];
				gmds::math::Point p2 = m_discretization_points[2];
				gmds::math::Point p3 = m_discretization_points[3];
				gmds::math::Point p4 = m_discretization_points[4];
				gmds::math::Point p5 = m_discretization_points[5];
				new_discretization[i] = gmds::math::Point(0.067 * p0.X() + 0.242 * p1.X() + 0.382 * p2.X() + 0.242 * p3.X() + 0.061 * p4.X() + 0.006 * p5.X(),
				                                          0.067 * p0.Y() + 0.242 * p1.Y() + 0.382 * p2.Y() + 0.242 * p3.Y() + 0.061 * p4.Y() + 0.006 * p5.Y(),
				                                          0.067 * p0.Z() + 0.242 * p1.Z() + 0.382 * p2.Z() + 0.242 * p3.Z() + 0.061 * p4.Z() + 0.006 * p5.Z());
			}
			else if (i == m_discretization_points.size() - 3) {
				gmds::math::Point p0 = m_discretization_points[m_discretization_points.size() - 1];
				gmds::math::Point p1 = m_discretization_points[m_discretization_points.size() - 2];
				gmds::math::Point p2 = m_discretization_points[m_discretization_points.size() - 3];
				gmds::math::Point p3 = m_discretization_points[m_discretization_points.size() - 4];
				gmds::math::Point p4 = m_discretization_points[m_discretization_points.size() - 5];
				gmds::math::Point p5 = m_discretization_points[m_discretization_points.size() - 6];
				new_discretization[i] = gmds::math::Point(0.067 * p0.X() + 0.242 * p1.X() + 0.382 * p2.X() + 0.242 * p3.X() + 0.061 * p4.X() + 0.006 * p5.X(),
				                                          0.067 * p0.Y() + 0.242 * p1.Y() + 0.382 * p2.Y() + 0.242 * p3.Y() + 0.061 * p4.Y() + 0.006 * p5.Y(),
				                                          0.067 * p0.Z() + 0.242 * p1.Z() + 0.382 * p2.Z() + 0.242 * p3.Z() + 0.061 * p4.Z() + 0.006 * p5.Z());
			}
			else {
				gmds::math::Point p0 = m_discretization_points[i - 3];
				gmds::math::Point p1 = m_discretization_points[i - 2];
				gmds::math::Point p2 = m_discretization_points[i - 1];
				gmds::math::Point p3 = m_discretization_points[i];
				gmds::math::Point p4 = m_discretization_points[i + 1];
				gmds::math::Point p5 = m_discretization_points[i + 2];
				gmds::math::Point p6 = m_discretization_points[i + 3];

				new_discretization[i] =
				   gmds::math::Point(0.006 * p0.X() + 0.061 * p1.X() + 0.242 * p2.X() + 0.382 * p3.X() + 0.242 * p4.X() + 0.061 * p5.X() + 0.006 * p6.X(),
				                     0.006 * p0.Y() + 0.061 * p1.Y() + 0.242 * p2.Y() + 0.382 * p3.Y() + 0.242 * p4.Y() + 0.061 * p5.Y() + 0.006 * p6.Y(),
				                     0.006 * p0.Z() + 0.061 * p1.Z() + 0.242 * p2.Z() + 0.382 * p3.Z() + 0.242 * p4.Z() + 0.061 * p5.Z() + 0.006 * p6.Z());
			}
		}

		m_discretization_points.clear();
		m_discretization_points.insert(m_discretization_points.end(), new_discretization.begin(), new_discretization.end());
	}
}
/*----------------------------------------------------------------------------*/
SingularityPoint::Slot *
SingularityLine::removeSlot()
{
	SingularityPoint::Slot *slot = m_slots.back();
	m_slots.pop_back();
	return slot;
}
/*----------------------------------------------------------------------------*/
void
SingularityLine::removeAllSingularityPoints()
{
	m_slots.clear();
}
/*----------------------------------------------------------------------------*/
double
SingularityLine::length() const
{
	double l = 0;
	for (unsigned int i = 0; i < m_discretization_points.size() - 1; i++) {
		gmds::math::Point pi = m_discretization_points[i];
		gmds::math::Point pj = m_discretization_points[i + 1];
		l += pi.distance(pj);
	}
	return l;
}
/*----------------------------------------------------------------------------*/
void
SurfaceSingularityLine::setTraversedFaces(const std::vector<gmds::TCellID> &AIDs)
{
	m_traversed_faces_id = AIDs;
}
/*----------------------------------------------------------------------------*/
bool
SurfaceSingularityLine::getIntersectionPoint(SurfaceSingularityLine *ALine, const gmds::math::Point &ARefPnt, const double ARefRadius, gmds::math::Point &APnt)
{
	std::vector<gmds::math::Point> &other_pnts = ALine->getDiscretizationPoints();
	std::vector<gmds::math::Segment> segments_1, segments_2;

	for (unsigned int i = 0; i < m_discretization_points.size() - 1; i++) {
		gmds::math::Point pi = m_discretization_points[i];
		gmds::math::Point pj = m_discretization_points[i + 1];
		if (ARefPnt.distance2(pi) <= ARefRadius) segments_1.push_back(gmds::math::Segment(pi, pj));
	}

	for (unsigned int i = 0; i < other_pnts.size() - 1; i++) {
		gmds::math::Point pi = other_pnts[i];
		gmds::math::Point pj = other_pnts[i + 1];
		if (ARefPnt.distance2(pi) <= ARefRadius) segments_2.push_back(gmds::math::Segment(pi, pj));
	}

	bool found_intersection = false;

	for (unsigned int i = 0; !found_intersection && i < segments_1.size(); i++) {

		gmds::math::Segment si = segments_1[i];

		for (unsigned int j = 0; !found_intersection && j < segments_2.size(); j++) {
			gmds::math::Segment sj = segments_2[j];

			gmds::math::Point intersection;
			double intersection_param;
			if (si.intersect2D(sj, intersection)) {
				found_intersection = true;
				APnt = intersection;
			}
		}
	}
	return found_intersection;
}
/*----------------------------------------------------------------------------*/
