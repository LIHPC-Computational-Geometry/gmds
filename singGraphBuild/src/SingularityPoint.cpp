/*----------------------------------------------------------------------------*/
/*
 * SingularityPoint.cpp
 *
 *  Created on: 13 juil. 2014
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
#include <gmds/math/Numerics.h>
/*----------------------------------------------------------------------------*/
#include <gmds/singGraphBuild/SingularityLine.h>
#include <gmds/singGraphBuild/SingularityPoint.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace std;
/*----------------------------------------------------------------------------*/
SingularityPoint::SingularityPoint(gmds::Mesh *AOwner, const ESingularityType AType, const ESingularityGeomType AGType) :
  m_mesh(AOwner), m_type(AType), m_geom_type(AGType)
{
}
/*----------------------------------------------------------------------------*/
SingularityPoint::SingularityPoint(
   gmds::Mesh *AOwner, const double AX, const double AY, const double AZ, const ESingularityType AType, const ESingularityGeomType AGType) :
  m_mesh(AOwner), m_type(AType), m_geom_type(AGType)
{
	setXYZ(AX, AY, AZ);
}
/*----------------------------------------------------------------------------*/
SingularityPoint::~SingularityPoint()
{
	for (unsigned int i = 0; i < m_slots.size(); i++) {
		if (m_slots[i]) delete m_slots[i];
	}
}
/*----------------------------------------------------------------------------*/
void
SingularityPoint::setXYZ(const double AX, const double AY, const double AZ)
{
	m_location.setXYZ(AX, AY, AZ);
}
/*----------------------------------------------------------------------------*/
void
SingularityPoint::setLocation(const gmds::math::Point &APnt)
{
	m_location.setXYZ(APnt.X(), APnt.Y(), APnt.Z());
}
/*----------------------------------------------------------------------------*/
gmds::math::Point
SingularityPoint::getLocation() const
{
	return m_location;
}
/*----------------------------------------------------------------------------*/
SingularityPoint::Slot *
SingularityPoint::newSlot(const gmds::math::Point &APnt,
                          const gmds::math::Vector3d &AVec,
                          const gmds::TCellID &ACellID,
                          const int ACellDim,
                          const bool AIsOnSurf,
                          SingularityLine *ALine,
                          const gmds::math::Vector3d ALineDirection)
{
	Slot *s = new Slot();
	s->location = APnt;
	s->direction = AVec;
	s->from_point = this;
	s->starting_cell_id = ACellID;
	s->starting_cell_dim = ACellDim;
	s->isOnSurface = AIsOnSurf;

	if (ALine) {
		s->line = ALine;
		s->isLaunched = true;
		s->line_direction = ALineDirection;
		// cout<<"s->line_direction "<<s->line_direction[0]<<" "<<s->line_direction[1]<<" "<<s->line_direction[2]<<endl;
	}
	else {
		s->line = 0;
		s->isLaunched = false;
	}

	m_slots.push_back(s);
	return m_slots.back();
}
/*----------------------------------------------------------------------------*/
SingularityPoint::Slot *
SingularityPoint::newGeomSlot(const gmds::math::Point &APnt, SingularityLine *ALine)
{

	// We have to compute the slot direction. To do it, we rely on the line
	// discretization
	std::vector<math::Point> line_disc = ALine->getDiscretizationPoints();
	if (line_disc.empty()) throw GMDSException("A line must be discretized before being connected to a sing. point");

	math::Point p_begin = line_disc[0];
	math::Point p_end = line_disc[line_disc.size() - 1];

	double d_begin = APnt.distance(p_begin);
	double d_end = APnt.distance(p_end);

	math::Vector3d slot_dir;
	if (d_begin <= d_end) {

		// start from the first point, we look for a point quite far enough
		int snd_index = 1;
		math::Point snd_pnt = line_disc[snd_index];
		while (math::near(p_begin.distance(snd_pnt), 0.0))
			snd_pnt = line_disc[++snd_index];

		slot_dir = math::Vector3d(p_begin, snd_pnt);
	}
	else {
		// start from the last point, we look for a point quite far enough
		int snd_index = line_disc.size() - 2;
		math::Point snd_pnt = line_disc[snd_index];
		while (math::near(p_end.distance(snd_pnt), 0.0))
			snd_pnt = line_disc[--snd_index];

		slot_dir = math::Vector3d(p_end, snd_pnt);
	}
	return newSlot(APnt,                // slot location
	               slot_dir,            // slot direction
	               0,                   // No linked cell (id)
	               0,                   // No linked cell (dim)
	               true,                // Always on surface (Maybe false in the future)
	               ALine,               // Connected line
	               slot_dir.opp());     // Line direction is the slot direction
}
/*----------------------------------------------------------------------------*/
bool
SingularityPoint::connectLine(SingularityLine *ALine, const math::Vector3d &AVec)
{
	// We keep the normalized vector equivalent to AVec

	bool found_slot = false;
	for (unsigned int i = 0; !found_slot && i < m_slots.size(); i++) {
		Slot *si = m_slots[i];
		gmds::math::Vector3d dir = si->direction;
		if (fabs(dir.angle(AVec)) < 15) {
			if (si->isLaunched) {
				throw GMDSException("Cannot add a line in an already used slot!!");
			}
			else {

				si->isLaunched = true;
				si->line = ALine;
				si->line_direction = AVec;
				found_slot = true;
			}
		}
	}

	if (!found_slot) return false;

	return true;
}
/*----------------------------------------------------------------------------*/
std::vector<SingularityLine *>
SingularityPoint::getLines()
{
	std::vector<SingularityLine *> adj_lines;
	for (unsigned int i = 0; i < m_slots.size(); i++) {
		Slot *si = m_slots[i];
		if (si->isLaunched) adj_lines.push_back(si->line);
	}

	return adj_lines;
}
/*----------------------------------------------------------------------------*/
SingularityLine *
SingularityPoint::nextLine(SingularityLine *ALine)
{
	// the singularity point must be classified onto a point, a curve or a
	// surface
	if (getGeomType() == VOLUME || getGeomType() == GEOM_UNDEF) return 0;

	// Entry line must be on the boundary
	if (!ALine->isOnBoundary()) return 0;

	// ALine and (*this) must be connected
	std::vector<SingularityPoint *> line_points = ALine->getEndPoints();
	if (line_points[0] != this && line_points[1] != this) return 0;

	// We get all the boundary lines of this
	// TODO: why is it searching among all slots???
	std::vector<SingularityLine *> candidates;
	std::vector<Slot *> candidates_slot;
	Slot *ref_slot = nullptr;

	for (Slot *slot : m_slots) {
		if (slot->isLaunched) {
			SingularityLine *line = slot->line;
			if (line == ALine) {
				ref_slot = slot;
			}
			else {
				candidates.push_back(line);
				candidates_slot.push_back(slot);
			}
		}
	}
	if (ref_slot == 0) throw GMDSException(" SingularityPoint::nextLine - Error, input line and point are not connected");
	// Now among all the candidates, we must find the next line in a counter-clockwise direction.

	math::Vector3d d_ref = ref_slot->direction;
	double angle = d_ref.angleIn02PI(candidates_slot[0]->direction);
	double next_id = 0;

	for (unsigned int i = 1; i < candidates_slot.size(); i++) {
		Slot *candidate_i = candidates_slot[i];
		math::Vector3d di = candidate_i->direction;
		double angle_i = d_ref.angleIn02PI(di);

		if (angle_i > angle) {
			angle = angle_i;
			next_id = i;
		}
	}
	SingularityLine *next_line = candidates[next_id];
	return next_line;
}
/*----------------------------------------------------------------------------*/
SingularityPoint::Slot *
SingularityPoint::addLine(SingularityLine *ALine)
{
	return newGeomSlot(this->getLocation(), ALine);
}
/*----------------------------------------------------------------------------*/
void
SingularityPoint::clearSlots()
{
	for (unsigned int i = 0; i < m_slots.size(); i++) {
		Slot *si = m_slots[i];
		delete si;
	}
	m_slots.clear();
}
/*----------------------------------------------------------------------------*/
void
SingularityPoint::sortSlots()
{
	std::sort(m_slots.begin(), m_slots.end(), [](Slot *a, Slot *b) -> bool {
		if (a->direction.Y() >= 0) {     // a between 0 and 180
			if (b->direction.Y() < 0)     // b between 180 and 360
				return false;
			return a->direction.X() < b->direction.X();
		}
		else {                           // a between 180 and 360
			if (b->direction.Y() > 0)     // b between 0 and 180
				return true;
			return a->direction.X() > b->direction.X();
		}
	});
}
/*----------------------------------------------------------------------------*/
template<> std::vector<Node>
SingularityPoint::getMesh<Node>()
{
	std::vector<Node> nodes;
	for (unsigned int i = 0; i < m_dim_owner.size(); i++) {
		if (m_dim_owner[i] == 0) {     // we have a node
			nodes.push_back(m_mesh->get<Node>(m_cell_id_owner[i]));
		}
	}
	return nodes;
}
/*----------------------------------------------------------------------------*/
template<> std::vector<Edge>
SingularityPoint::getMesh<Edge>()
{
	std::vector<Edge> edges;
	for (unsigned int i = 0; i < m_dim_owner.size(); i++) {
		if (m_dim_owner[i] == 1) {     // we have a node
			edges.push_back(m_mesh->get<Edge>(m_cell_id_owner[i]));
		}
	}
	return edges;
}
/*----------------------------------------------------------------------------*/
template<> std::vector<Face>
SingularityPoint::getMesh<Face>()
{
	std::vector<Face> faces;
	for (unsigned int i = 0; i < m_dim_owner.size(); i++) {
		if (m_dim_owner[i] == 2) {     // we have a node
			faces.push_back(m_mesh->get<Face>(m_cell_id_owner[i]));
		}
	}
	return faces;
}
/*----------------------------------------------------------------------------*/
template<> std::vector<Region>
SingularityPoint::getMesh<Region>()
{
	std::vector<Region> regions;
	for (unsigned int i = 0; i < m_dim_owner.size(); i++) {
		if (m_dim_owner[i] == 3) {     // we have a node
			regions.push_back(m_mesh->get<Region>(m_cell_id_owner[i]));
		}
	}
	return regions;
}
/*----------------------------------------------------------------------------*/
SingularityPoint::Slot::Slot() :
  isLaunched(false),
  location(math::Point(0, 0, 0)),
  direction(math::Vector3d(0, 0, 0)),
  starting_cell_id(NullID),
  starting_cell_dim(0),
  line(0),
  isOnSurface(false),
  isFreeze(false)
{
	;
}
/*----------------------------------------------------------------------------*/
bool
VolumeSingularityPoint::addLine(SingularityLine *ALine, const gmds::Face &AIncomingFace)
{
	// We keep the normalized vector equivalent to AVec

	bool found_slot = false;
	for (unsigned int i = 0; !found_slot && i < m_slots.size(); i++) {
		Slot *si = m_slots[i];
		if (AIncomingFace.id() == si->starting_cell_id && si->starting_cell_dim == 2) {

			if (si->isLaunched) {
				throw GMDSException("Cannot add a line in an already used slot!!");
			}
			else {
				m_slots[i]->isLaunched = true;
				m_slots[i]->line = ALine;
				found_slot = true;
			}
		}

	}     // for (unsigned int i = 0; !found_slot && i < m_slots.size(); i++)

	if (!found_slot) return false;

	return true;
}
/*----------------------------------------------------------------------------*/
gmds::Node
VertexSingularityPoint::getMeshNode()
{
	if (m_cell_id_owner.size() != 1) throw GMDSException("A geom. singularity point must be associated to EXACTLY one mesh cell");

	if (m_dim_owner[0] != 0) throw GMDSException("A surf. singularity point must be associated to one 0-cell");

	return m_mesh->get<Node>(m_cell_id_owner[0]);
}
/*----------------------------------------------------------------------------*/
void
SingularityPoint::addMeshNode(const gmds::Node &ANode)
{
	m_cell_id_owner.push_back(ANode.id());
	m_dim_owner.push_back(0);
}

/*----------------------------------------------------------------------------*/
void
SingularityPoint::addMeshEdge(const gmds::Edge &AEdge)
{
	m_cell_id_owner.push_back(AEdge.id());
	m_dim_owner.push_back(1);
}
/*----------------------------------------------------------------------------*/
gmds::Face
SurfaceSingularityPoint::getMeshFace()
{
	if (m_cell_id_owner.size() != 1) throw GMDSException("A surf. singularity point must be associated to EXACTLY one mesh cell");

	if (m_dim_owner[0] != 2) throw GMDSException("A surf. singularity point must be associated to one 2-cell");

	return m_mesh->get<Face>(m_cell_id_owner[0]);
}
/*----------------------------------------------------------------------------*/
void
SingularityPoint::addMeshFace(const gmds::Face &AFace)
{
	m_cell_id_owner.push_back(AFace.id());
	m_dim_owner.push_back(2);
}
/*----------------------------------------------------------------------------*/
bool
SurfaceSingularityPoint::addLineFromVolume(SingularityLine *ALine)
{
	bool found_slot = false;
	for (unsigned int i = 0; !found_slot && i < m_slots.size(); i++) {
		Slot *si = m_slots[i];
		if (!si->isOnSurface) {
			if (si->isLaunched) {
				//	throw GMDSException("Cannot add a line in an already used slot!!");
			}
			else {

				si->isLaunched = true;
				si->line = ALine;
				found_slot = true;
			}
		}
	}

	if (!found_slot) return false;

	return true;
}
/*----------------------------------------------------------------------------*/
gmds::Region
VolumeSingularityPoint::getMeshRegion()
{
	if (m_cell_id_owner.size() != 1) throw GMDSException("A surf. singularity point must be associated to EXACTLY one mesh cell");

	if (m_dim_owner[0] != 3) throw GMDSException("A surf. singularity point must be associated to one 3-cell");

	return m_mesh->get<Region>(m_cell_id_owner[0]);
}
/*----------------------------------------------------------------------------*/
std::vector<gmds::Region>
VolumeSingularityPoint::getMeshRegions()
{

	for (unsigned int i = 0; i < m_dim_owner.size(); i++) {
		if (m_dim_owner[i] != 3) {
			std::cout << "A vol. singularity point must only be associated to 3-cells" << std::endl;
			throw GMDSException("A vol. singularity point must only be associated to 3-cells");
		}
	}
	std::vector<gmds::Region> r;
	std::cout << "nb cells: " << m_cell_id_owner.size() << std::endl;

	for (unsigned int i = 0; i < m_cell_id_owner.size(); i++) {
		std::cout << "cell " << m_cell_id_owner[i] << std::endl;
	}
	for (unsigned int i = 0; i < m_cell_id_owner.size(); i++) {
		r.push_back(m_mesh->get<Region>(m_cell_id_owner[i]));
	}
	return r;
}
/*----------------------------------------------------------------------------*/
void
SingularityPoint::addMeshRegion(gmds::Region &ARegion)
{
	m_cell_id_owner.push_back(ARegion.id());
	m_dim_owner.push_back(3);
}
/*----------------------------------------------------------------------------*/
