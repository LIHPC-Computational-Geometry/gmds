/*----------------------------------------------------------------------------*/
/*
 * SingularityGraphBuilder2D.cpp
 *
 *  Created on: 13 juil. 2014
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
//#include <gmds/ig/IG.h>
#include "gmds/math/Constants.h"
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/io/GMSHWriter.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/MeditReader.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/math/AxisAngleRotation.h>
#include <gmds/math/Chart.h>
#include <gmds/math/Cross.h>
#include <gmds/math/Cross2D.h>
#include <gmds/math/Numerics.h>
#include <gmds/math/Quaternion.h>
#include <gmds/math/Ray.h>
#include <gmds/math/Segment.h>
#include <gmds/math/Triangle.h>
#include <gmds/singGraphBuild/SingularityGraphBuilder2D.h>
#include <gmds/singGraphBuild/Tools.h>
#include <gmds/utils/Variable.h>

#include <glpk.h>
#include <gmds/frame/LaplaceCross2D.h>
/*----------------------------------------------------------------------------*/
#include <chrono>
#include <list>
#include <sstream>
#include <unordered_map>
#include <vector>
typedef std::chrono::high_resolution_clock Clock;
/*----------------------------------------------------------------------------*/

namespace gmds {
namespace {

void
writeQuadMesh(Mesh *mesh, const std::string &AFileName)
{
	IGMeshIOService ioService(mesh);
	GMSHWriter gmshWriter(&ioService);
	gmshWriter.setCellOptions(N | F);
	gmshWriter.setDataOptions(0);
	gmshWriter.write(AFileName);
}

std::vector<double>
computeAngleAtSingularity(SingularityPoint *p)
{
	std::vector<double> angles;
	auto addAngle = [&angles](const SingularityPoint::Slot *slotA, const SingularityPoint::Slot *slotB) {
		if (slotA->isOnSurface || slotB->isOnSurface) {     // only angle with surface slot
			const double angle = slotA->line_direction.angle(slotB->line_direction);
			angles.push_back(angle);
		}
	};
	const auto slots = p->getSlots();
	for (int i = 0; i < slots.size() - 1; ++i) {
		addAngle(slots[i], slots[i + 1]);
	}
	addAngle(slots.back(), slots.front());
	return angles;
}

// Manipulate Singularity Point locations from a grapgh in order to optimize the ratio of linked edges
class LineRelaxator
{
 public:
	LineRelaxator(SingularityGraph *graph, const double meanEdgeLength);

	void run();
	double getWorstAngle();
	double getWorstEdgeRatio();

 private:
	void computeFixedGroupAndFixedLine();
	void computeMeanEdgeLengthByGroup();
	void applySpringForceOnPoints();
	void movePoints();
	void MoveIfNewLocationIsWorthTheAnglePenalty(SingularityPoint *singPoint, const gmds::math::Point newLocation);

	SingularityGraph *m_graph;

	double m_step;
	std::vector<double> m_lineSpringLenght;
	std::vector<double> m_minLengthByGroup;
	std::vector<double> m_maxLengthByGroup;
	std::vector<SingularityLine *> m_minLineByGroup;
	std::vector<SingularityLine *> m_maxLineByGroup;

	std::vector<size_t> m_numberOfLinePerGroup;
	std::vector<bool> m_groupIsFixed;
	std::unordered_map<size_t, bool> m_lineIsfixed;
	std::unordered_map<size_t, gmds::math::Vector3d> m_directions;

	std::unordered_map<size_t, std::vector<double>> m_angleAtSingularity;
};

class Timer
{
	using Clock = std::chrono::system_clock;

 public:
	Timer(const std::string &name)
	{
		start(name);
	}
	void stopAndRestart(const std::string &name)
	{
		stop();
		start(name);
	}
	void start(const std::string &name)
	{
		m_name = name;
		m_t0 = Clock::now();
	}
	void stop()
	{
		const auto tfinal = Clock::now();
		std::cout << m_name << " " << std::chrono::duration_cast<std::chrono::milliseconds>(tfinal - m_t0).count() << " milliseconds" << std::endl;
	}

 private:
	std::string m_name;
	std::chrono::system_clock::time_point m_t0;
};

}     // namespace

SingularityGraphBuilder2D::SingularityGraphBuilder2D(Mesh *AMesh, Variable<math::Cross2D> *AField, const bool ABuildGeomSing) :
  m_mesh(AMesh), m_field(AField), m_tool(AMesh, AField), m_output_directory_name(""), m_graph(AMesh)
{
	m_build_geometric_singularities = ABuildGeomSing;

	if (m_ATolerance < 0.01 || m_ATolerance > 0.1) throw GMDSException("SingularityGraphBuilder2D: Tolerance must be taken in [0.01,0.1]");

	double x_min, y_min, z_min, x_max, y_max, z_max;

	math::Point current_pnt = m_mesh->get<Node>(0).getPoint();
	x_min = current_pnt.X();
	x_max = current_pnt.X();
	y_min = current_pnt.Y();
	y_max = current_pnt.Y();
	z_min = current_pnt.Z();
	z_max = current_pnt.Z();

	for (auto n_id : m_mesh->nodes()) {
		math::Point current_pnt = m_mesh->get<Node>(n_id).getPoint();
		if (current_pnt.X() < x_min)
			x_min = current_pnt.X();
		else if (current_pnt.X() > x_max)
			x_max = current_pnt.X();

		if (current_pnt.Y() < y_min)
			y_min = current_pnt.Y();
		else if (current_pnt.Y() > y_max)
			y_max = current_pnt.Y();

		if (current_pnt.Z() < z_min)
			z_min = current_pnt.Z();
		else if (current_pnt.Z() > z_max)
			z_max = current_pnt.Z();
	}
	math::Point p_min(x_min, y_min, z_min);
	math::Point p_max(x_max, y_max, z_max);
	m_mesh_radius = p_min.distance(p_max);
	m_confusing_distance = m_ATolerance * m_mesh_radius;
}

void
SingularityGraphBuilder2D::initTestPrescribedSing()
{
	if (m_withGlobalComments) std::cout << "testPrescribedSing" << testPrescribedSing << std::endl;
	unsigned int singularity_type = 3;
	m_singularities_3.clear();
	m_singularities_5.clear();

	Face singularTri = m_mesh->get<Face>(4803);
	m_singularities_3.push_back(singularTri);
	m_mesh->mark(singularTri, m_mark_faces_with_sing_point);
	math::Point singleCenterTri;
	math::Cross2D singleCenterTriCross;
	constructOneCenterTriangleCross(singularTri, singleCenterTri, singleCenterTriCross);
	std::vector<math::Point> slot_points(singularity_type);
	std::vector<double> slot_param;
	std::vector<int> slot_cell_dim(singularity_type);
	std::vector<TCellID> slot_cell_id(singularity_type);
	// SINGULARITY POINT CREATION
	SurfaceSingularityPoint *singularity = m_graph.newSurfacePoint();
	singularity->setLocation(singleCenterTri);
	singularity->addMeshFace(singularTri);
	std::vector<math::Vector3d> vectors_i = singleCenterTriCross.componentVectors();
	vector<math::Vector3d> vectors_iNew(3);
	vectors_iNew[0] = vectors_i[0];
	vectors_iNew[1] = vectors_i[1];
	vectors_iNew[2] = vectors_i[3];

	math::Point intersectionPnt;
	double intersectionParam;
	for (unsigned int i = 0; i < singularity_type; i++) {
		math::Ray from_ray(singleCenterTri, vectors_iNew[i]);
		vector<Edge> currentEdges = singularTri.get<Edge>();
		for (unsigned int tt = 0; tt < 3; tt++) {
			vector<Node> currentNodes = currentEdges[tt].get<Node>();
			math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
			if (from_ray.SecondMetIntersect2D(oppSeg, intersectionPnt, intersectionParam, m_temp_epsilon)) {
				slot_points[i] = intersectionParam * currentNodes[0].getPoint() + (1 - intersectionParam) * currentNodes[1].getPoint();
				slot_cell_dim[i] = 1;
				slot_cell_id[i] = currentEdges[tt].id();
			}
		}
	}

	for (unsigned int i = 0; i < singularity_type; i++) {
		singularity->newSlot(slot_points[i], vectors_iNew[i], slot_cell_id[i], slot_cell_dim[i], true, 0);
	}

	visualizeCrossVectors();
}

void
SingularityGraphBuilder2D::visualizeCrossVectors()
{
	Mesh m(MeshModel(DIM3 | F | N | F2N));

	for (int i = 0; i < m_triangle_centers_cross.size(); i++) {

		const auto &componentVectors = m_triangle_centers_cross[i].componentVectors();
		const auto &faceCenter = m_triangle_centers[i];
		const Node nodeCenter = m.newNode(faceCenter);

		for (const auto &vect : componentVectors) {

			const auto p = faceCenter + vect * m_mean_edge_length * 0.7;
			Node node = m.newNode(p);
			m.newTriangle(nodeCenter, nodeCenter, node);
		}
	}
	IGMeshIOService meshIoServ(&m);
	VTKWriter writer(&meshIoServ);
	writer.setCellOptions(N | F);
	writer.setDataOptions(N | F);

	std::stringstream file_name;
	file_name << m_output_directory_name << "-CrossCenterTri.vtk";
	writer.write(file_name.str());
}

std::unique_ptr<Mesh>
SingularityGraphBuilder2D::getQuadMesh()
{
	const auto singPoints = m_graph.getPoints();
	const auto patchs = m_graph.getSurfacePatchs();
	auto mesh = std::make_unique<Mesh>(MeshModel(DIM3 | F | N | F2N));

	std::unordered_map<TCellID, TCellID> singID2NodeID;
	singID2NodeID.reserve(singPoints.size());
	for (const auto singPoint : singPoints) {
		const auto newNode = mesh.get()->newNode(singPoint->getLocation());
		singID2NodeID[singPoint->getNumber()] = newNode.id();
	}

	std::vector<SingularityPoint *> nodeBuffer(4, nullptr);
	for (const auto patch : patchs) {
		patch->getPoints(nodeBuffer);
		mesh.get()->newQuad(singID2NodeID[nodeBuffer[0]->getNumber()], singID2NodeID[nodeBuffer[1]->getNumber()],     //
		                    singID2NodeID[nodeBuffer[2]->getNumber()], singID2NodeID[nodeBuffer[3]->getNumber()]);
	}
	return mesh;
}

SingularityGraphBuilder2D &
SingularityGraphBuilder2D::setQuadMeshSmoothingEnable(bool enableSmoothing)
{
	m_enableQuadMeshSmoothing = enableSmoothing;
	return *this;
}

SingularityGraphBuilder2D &
SingularityGraphBuilder2D::setDebugFilesWritingEnable(bool enableDebugFilesWriting)
{
	m_enableDebugFilesWriting = enableDebugFilesWriting;
	return *this;
}

double
SingularityGraphBuilder2D::computeMeanEdgeLength()
{
	double mean_edge_length = 0.0;
	for (auto e_id : m_mesh->edges()) {
		mean_edge_length += m_mesh->get<Edge>(e_id).length();
	}
	return mean_edge_length / m_mesh->getNbEdges();
}

void
SingularityGraphBuilder2D::execute()
{
	std::cout << "========================================" << std::endl;
	std::cout << "Start singularity graph generation " << std::endl;

	m_graph = SingularityGraph(m_mesh);     // new fresh graph
	//==================================================================
	// Boolean marks initialization
	//==================================================================
	m_mark_faces_with_sing_point = m_mesh->newMark<Face>();
	m_mark_faces_with_sing_line = m_mesh->newMark<Face>();
	m_mark_forbiddenBdryEdge = m_mesh->newMark<Edge>();

	m_original_faces_number = m_mesh->getNbFaces();
	m_original_nodes_number = m_mesh->getNbNodes();

	//========================================================================
	// STEP 1 - Detection of singularity and slot creation
	//========================================================================

	/* We proceed as follows:
	- first the singular triangles are detected, having as input the frame/cross field for the mesh;
	- for each such singular triangle, the location of the singularity point is detected inside the triangle (for now
	the location is at the center of the triangle) and also the direction of the slot directions is computed */

	auto timer = Timer("Singular triangle detection and slot creation");

	detectSingularTriangles();

	for (const auto &currentFace : m_singularities_3)
		createSingPointAndSlots(currentFace);

	for (const auto &currentFace : m_singularities_5)
		createSingPointAndSlots(currentFace);

	timer.stopAndRestart("Cross field interpolation");
	constructCenterTriangleCrosses();

	timer.stopAndRestart("Mean edge length computation");
	m_mean_edge_length = computeMeanEdgeLength();
	m_temp_epsilon = 0.003 * m_mean_edge_length;

	//========================================================================
	// STEP 2 - Geometrical features are added in the singularity graph
	//========================================================================
	timer.stopAndRestart("Geometric singularities creation");

	/* for each interior boundary loop (if it exists), an artificial curve singularity point is created
	 * (at a random location, one one mesh vertex); this artificial singularity point must be removed
	 * after we have detected the singularity graph*/

	vector<CurveSingularityPoint *> artificialSingPointsCreated;
	addGeometryToSingularityGraph(artificialSingPointsCreated);

	if (m_build_geometric_singularities) {
		buildSlotsOfBoundarySingularities();
	}

	if (m_enableDebugFilesWriting) {
		writeSingularityPointsAndSlots();
		visualizeCrossVectors();
	}

	//========================================================================
	// STEP 3 - Singularity line building
	//========================================================================
	timer.stopAndRestart("Singularity lines creation");

	createSingularityLines();

	//========================================================================
	// STEP 4 - Detect singularity lines intersection and split them
	//========================================================================
	std::cout << "STEP 6 - Detect singularity lines intersection and split them" << std::endl;

	timer.stopAndRestart("Line intersections");
	detectLineIntersections();

	//========================================================================
	// STEP 5 - Build surface patchs
	//========================================================================

	timer.stopAndRestart("Quad Mesh building");

	deleteArtificialNodes(artificialSingPointsCreated);
	if (m_enableDebugFilesWriting) {
		writeOutputSingle("boundary_line_final");
	}

	m_graph.buildSurfacePatchs();

	if (m_enableQuadMeshSmoothing) {
		if (m_enableDebugFilesWriting) {
			writeQuadMesh(getQuadMesh().get(), m_output_directory_name + "gmshPatch.msh");
		}
		timer.stopAndRestart("Graph Optimization");
		LineRelaxator(&m_graph, m_mean_edge_length).run();
		timer.stop();
	}
	if (m_enableDebugFilesWriting) {
		writeOutputPatches("patchs");
		writeQuadMesh(getQuadMesh().get(), m_output_directory_name + "gmshPatch2.msh");
	}

	//========================================================================
	// cleaning
	//========================================================================
	m_mesh->unmarkAll<Face>(m_mark_faces_with_sing_point);
	m_mesh->unmarkAll<Face>(m_mark_faces_with_sing_line);
	m_mesh->unmarkAll<Edge>(m_mark_forbiddenBdryEdge);
	m_mesh->freeMark<Face>(m_mark_faces_with_sing_point);
	m_mesh->freeMark<Face>(m_mark_faces_with_sing_line);
	m_mesh->freeMark<Edge>(m_mark_forbiddenBdryEdge);

	std::cout << "--> Nb generated faces: " << m_graph.getNbSurfacePatches() << std::endl;
	std::cout << "========================================" << std::endl;
}

/*----------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::deleteArtificialNodes(vector<CurveSingularityPoint *> &artificialSingPointsCreated)
{
	// Delete the artifial nodes that have been created in addGeometryToSingularityGraph()
	for (unsigned int i = 0; i < artificialSingPointsCreated.size(); i++) {

		auto artificialSingPoint = artificialSingPointsCreated[i];
		auto slots = artificialSingPoint->getSlots();
		if (slots.size() > 2) continue;

		auto line1 = static_cast<CurveSingularityLine *>(slots[0]->line);
		auto line2 = static_cast<CurveSingularityLine *>(slots[1]->line);

		auto &points1 = line1->getDiscretizationPoints();
		auto &edges1 = line1->getMeshEdges();
		const auto &points2 = line2->getDiscretizationPoints();
		const auto &edges2 = line2->getMeshEdges();

		// line1 is kept
		if (points1.front() == artificialSingPoint->getLocation()) {
			auto &endSlots = line1->getSlots();
			std::reverse(endSlots.begin(), endSlots.end());
			std::reverse(points1.begin(), points1.end());
			std::reverse(edges1.begin(), edges1.end());
		}
		if (points1.back() != points2.front()) {
			points1.insert(points1.end(), points2.rbegin(), points2.rend());
			edges1.insert(edges1.end(), edges2.rbegin(), edges2.rend());
		}
		else {
			points1.insert(points1.end(), points2.begin(), points2.end());
			edges1.insert(edges1.end(), edges2.begin(), edges2.end());
		}

		// update lines sing points
		line1->removeSlot();
		const auto singPointsOfLine2 = line2->getEndPoints();
		auto singPointToUpdate = singPointsOfLine2[0] == artificialSingPoint ? singPointsOfLine2[1] : singPointsOfLine2[0];
		for (auto slotToUpdate : singPointToUpdate->getSlots()) {
			if (slotToUpdate->line && slotToUpdate->line->getNumber() == line2->getNumber()) {
				slotToUpdate->line = line1;
				line1->addSlot(slotToUpdate);
			}
		}

		m_graph.removeLine(line2);
		m_graph.removePoint(artificialSingPoint);
	}
}

/*---------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::detectLineIntersections()
{
	std::vector<math::Point> added_points;
	std::vector<Face> candidate_faces;
	if (m_withGlobalComments) cout << "detectLineIntersections " << endl;
	//========================================================================
	// We detect intersection faces before splitting any curves. After, it will be
	// too late to detect them.
	//=======================================================================
	vector<SurfaceSingularityLine *> surf_lines = m_graph.getSurfaceLines();
	for (auto f_id : m_mesh->faces()) {
		Face f = m_mesh->get<Face>(f_id);
		std::vector<SurfaceSingularityLine *> current_lines = getSingularityLinesIn(f);

		if (current_lines.size() > 1) {
			if (!m_mesh->isMarked(f, m_mark_faces_with_sing_point)) {
				candidate_faces.push_back(f);
			}
		}
	}     // for (; !it_faces.isDone(); it_faces.next())

	// writeTestMeshTriangles(candidate_faces, std::string("candidateFacesForLineSplit.vtk"));

	for (unsigned int i = 0; i < candidate_faces.size(); i++) {
		Face f = candidate_faces[i];
		if (m_withGlobalComments) cout << "f.id() " << f.id() << endl;
		std::vector<SurfaceSingularityLine *> current_lines = getSingularityLinesIn(f);

		if (current_lines.size() >= 2) {
			if (current_lines.size() == 2)
				createLineIntersection(current_lines[0], current_lines[1], f, added_points);
			else
				createLineIntersection(current_lines, f, added_points);
		}
	}
}

/*---------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::createLineIntersection(SurfaceSingularityLine *ALine1,
                                                  SurfaceSingularityLine *ALine2,
                                                  Face &AFace,
                                                  std::vector<math::Point> &AAddedPoints)
{
	// make sure that lines are not sharing an end point
	for (const auto &p1 : ALine1->getEndPoints())
		for (const auto &p2 : ALine2->getEndPoints())
			if (p1 == p2) return;

	// We look for the geometrical intersection point of the curve lines
	SurfaceSingularityLine *l0 = ALine1;
	SurfaceSingularityLine *l1 = ALine2;

	const math::Point face_center = AFace.center();
	const std::vector<Node> &face_nodes = AFace.get<Node>();
	const unsigned int nb_face_nodes = face_nodes.size();
	double face_radius = 0;
	for (unsigned int i = 0; i < nb_face_nodes; i++) {
		const math::Point &pi = face_nodes[i].getPoint();
		const math::Point &pj = face_nodes[(i + 1) % nb_face_nodes].getPoint();
		const double distance_ij = pi.distance2(pj);
		if (distance_ij > face_radius) face_radius = distance_ij;
	}

	math::Point p;
	if (l0->getIntersectionPoint(l1, face_center, face_radius, p)) {
		math::Triangle triangle(face_nodes[0].getPoint(), face_nodes[1].getPoint(), face_nodes[2].getPoint());

		if (triangle.isIn2ndMethod(p)) {
			for (const auto &singPoint : m_graph.getSurfacePoints()) {
				if (math::near(p.distance(singPoint->getLocation()), 0.0)) return;
			}
			for (const auto &point : AAddedPoints) {
				if (math::near(p.distance(point), 0.0)) return;
			}
			AAddedPoints.push_back(p);
			//==============================================================
			// Creation of the singularity point
			//==============================================================*
			SurfaceSingularityPoint *new_pnt = m_graph.newSurfacePoint();
			new_pnt->setLocation(p);
			new_pnt->addMeshFace(AFace);
			//==================================================================
			/* Creation of a new sing point and splitting of intersected sing. lines*/
			//==================================================================
			m_graph.splitSurfaceLine(new_pnt, l0);
			m_graph.splitSurfaceLine(new_pnt, l1);
		}
	}
}

/*---------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::createLineIntersection(std::vector<SurfaceSingularityLine *> &ALines, Face &AFace, std::vector<math::Point> &AAddedPoints)
{
	// We look for the geometrical intersection point of the curve lines
	std::vector<Node> face_nodes = AFace.get<Node>();
	double face_radius = 0;
	unsigned int nb_face_nodes = face_nodes.size();
	for (unsigned int i = 0; i < nb_face_nodes; i++) {
		math::Point pi = face_nodes[i].getPoint();
		math::Point pj = face_nodes[(i + 1) % nb_face_nodes].getPoint();
		double distance_ij = pi.distance2(pj);
		if (distance_ij > face_radius) face_radius = distance_ij;
	}

	const math::Point face_center = AFace.center();
	std::list<std::pair<SurfaceSingularityLine *, SurfaceSingularityPoint *>> linesToSplit;
	for (unsigned int i = 0; i < ALines.size() - 1; i++) {
		SurfaceSingularityLine *li = ALines[i];
		for (unsigned int j = i + 1; j < ALines.size(); j++) {
			SurfaceSingularityLine *lj = ALines[j];
			math::Point p;
			if (li->getIntersectionPoint(lj, face_center, face_radius, p)) {
				math::Triangle ATriangle(face_nodes[0].getPoint(), face_nodes[1].getPoint(), face_nodes[2].getPoint());
				if (ATriangle.isIn2ndMethod(p)) {
					bool already_added = false;
					vector<SurfaceSingularityPoint *> allSingPoints = m_graph.getSurfacePoints();
					for (unsigned int i = 0; i < allSingPoints.size(); i++) {
						if (math::near(p.distance(allSingPoints[i]->getLocation()), 0.0)) {
							already_added = true;
						}
					}
					for (unsigned int t = 0; !already_added && t < AAddedPoints.size(); t++) {
						if (math::near(p.distance(AAddedPoints[t]), 0.0)) already_added = true;
					}
					if (!already_added) {
						AAddedPoints.push_back(p);

						SurfaceSingularityPoint *new_pnt = m_graph.newSurfacePoint();
						new_pnt->setLocation(p);
						new_pnt->addMeshFace(AFace);

						linesToSplit.emplace_back(std::make_pair(li, new_pnt));
						linesToSplit.emplace_back(std::make_pair(lj, new_pnt));
					}
				}
			}
		}
	}
	while (!linesToSplit.empty()) {
		const auto pair = linesToSplit.back();
		linesToSplit.pop_back();
		auto line = pair.first;
		auto sing = pair.second;
		m_graph.splitSurfaceLine(sing, line);

		// the line we just cut might have had other intersections with other lines:
		const auto lineID = line->getNumber();
		const auto newLine = m_graph.getLines().back();
		for (auto &otherPair : linesToSplit) {
			auto &otherLine = otherPair.first;
			const auto intersectionPoint = otherPair.second->getLocation();
			if (otherLine->getNumber() == lineID) {
				// search on which of the two part the intersection is
				const std::vector<math::Point> &pnts = newLine->getDiscretizationPoints();
				for (unsigned int i = 0; i < pnts.size() - 1; i++) {
					math::Segment sij(pnts[i], pnts[i + 1]);
					if (sij.isIn2ndMethod(intersectionPoint)) {
						otherLine = dynamic_cast<SurfaceSingularityLine *>(newLine);
						break;
					}
				}
			}
		}
	}
}

/*---------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::detectSingularTriangles()
{
	m_singularities_3.clear();
	m_singularities_5.clear();

	for (auto f_id : m_mesh->faces()) {
		Face current = m_mesh->get<Face>(f_id);
		std::vector<TCellID> nodeIDs = current.getIDs<Node>();
		int ID1 = nodeIDs[0];
		int ID2 = nodeIDs[1];
		int ID3 = nodeIDs[2];

		// singularities detected along boundaries will almost always be unvalid
		if (m_mesh->isMarked<Node>(ID1, m_mark_nodes_on_curve)        //
		    || m_mesh->isMarked<Node>(ID2, m_mark_nodes_on_curve)     //
		    || m_mesh->isMarked<Node>(ID3, m_mark_nodes_on_curve)) {
			continue;
		}

		math::Cross2D c1 = (*m_field)[ID1];
		math::Cross2D c2 = (*m_field)[ID2];
		math::Cross2D c3 = (*m_field)[ID3];

		int index = math::Cross2D::index(c1, c2, c3);

		if (index == 1) {
			m_singularities_5.push_back(current);
			m_mesh->mark(current, m_mark_faces_with_sing_point);
		}
		else if (index == -1) {
			m_singularities_3.push_back(current);
			m_mesh->mark(current, m_mark_faces_with_sing_point);
		}
	}     // for (; !it_regions.isDone(); it_regions.next())
}

/*---------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::initMarks(const int AMarkNodePnt, const int AMarkNodeCrv, const int AMarkEdgeCrv, const int AMarkNodeForbiddenBdry)
{
	m_mark_nodes_on_point = AMarkNodePnt;
	m_mark_nodes_on_curve = AMarkNodeCrv;
	m_mark_edges_on_curve = AMarkEdgeCrv;
	m_mark_forbiddenBdryNode = AMarkNodeForbiddenBdry;
}

/*----------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::initConfusingBalls(SingularityPoint *APnt)
{
	math::Point sing_location = APnt->getLocation();
	bool found_face = false;
	//========================================================================
	// Faces
	//========================================================================
	// int nb_faces = 0;
	double m_confusing_distance_temp = 0.0;
	while (!found_face) {
		m_confusing_distance_temp = m_confusing_distance_temp + m_confusing_distance;
		for (auto f_id : m_mesh->faces()) {
			Face currentFace = m_mesh->get<Face>(f_id);
			math::Point center = currentFace.center();
			if (center.distance(sing_location) < m_confusing_distance_temp) {
				m_faces_to_singularity_on_surf[currentFace.id()] = APnt;
				found_face = true;
				// nb_faces++;
			}
		}
	}
	//========================================================================
	// Edges
	//========================================================================
	for (auto e_id : m_mesh->edges()) {
		Edge currentEdge = m_mesh->get<Edge>(e_id);
		math::Point center = currentEdge.center();
		if (center.distance(sing_location) < m_confusing_distance) {
			m_edges_to_singularity_on_surf[currentEdge.id()] = APnt;
		}
	}
	//========================================================================
	// Nodes
	//========================================================================
	for (auto n_id : m_mesh->nodes()) {
		Node currentNode = m_mesh->get<Node>(n_id);
		math::Point center = currentNode.getPoint();
		if (center.distance(sing_location) < m_confusing_distance) {
			m_nodes_to_singularity_on_surf[currentNode.id()] = APnt;
			vector<Face> adjacent_triangles = currentNode.get<Face>();
			for (unsigned int i = 0; i < adjacent_triangles.size(); i++) {
				math::Point centerFace = adjacent_triangles[i].center();
				// if(centerFace.distance(sing_location)< m_confusing_distance_temp) {
				m_faces_to_singularity_on_surf[adjacent_triangles[i].id()] = APnt;
				//}
			}
		}
	}
}

/*----------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::computeSingPointInfo(const Face &AFace, math::Point &APosSing)
{
	std::vector<Node> currentNodes = AFace.get<Node>();
	Node node1 = currentNodes[0];
	Node node2 = currentNodes[1];
	Node node3 = currentNodes[2];

	math::Vector3d v1 = (*m_field)[node1.id()].referenceVector();
	math::Vector3d v2 = (*m_field)[node2.id()].referenceVector();
	math::Vector3d v3 = (*m_field)[node3.id()].referenceVector();

	math::Point pointA = currentNodes[0].getPoint();
	math::Point pointB = currentNodes[1].getPoint();
	math::Point pointC = currentNodes[2].getPoint();

	// We solve a 2x2 Ax=b system with x = (alpha,beta)
	double alpha = 0, beta = 0, gamma = 0;
	math::Vector3d A = v1 - v3;
	math::Vector3d B = v2 - v3;
	math::Vector3d C = v3.opp();

	double dA = A[0] * B[1] - A[1] * C[0];
	if (dA != 0) {
		double Dx = C[0] * B[1] - C[1] * B[0];
		double Dy = A[0] * C[1] - A[1] * C[0];
		alpha = Dx / dA;
		beta = Dy / dA;
		gamma = 1 - alpha - beta;
	}
	else {
		throw GMDSException("Null Determinant in the computation of a singularity point location");
	}
	// WARNING: TODO the singularity point will always be placed at the traingle center...
	if (true) {     // dA<0){
		alpha = 0.333;
		beta = 0.333;
		gamma = 1 - alpha - beta;
	}
	//(alpha, beta, gamma) are the barycentric coordinates where the cross field vanishes
	// We can then compute the singularity point location
	double x = alpha * node1.X() + beta * node2.X() + gamma * node3.X();
	double y = alpha * node1.Y() + beta * node2.Y() + gamma * node3.Y();

	APosSing.setXYZ(x, y, 0);
}

/*----------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::createSingPointAndSlots(const Face &AFace)
{
	//=========================================================================
	// WE BEGIN BY COMPUTING THE LOCATION OF THE SINGULARITY POINT IN AFace
	//=========================================================================
	math::Point s;
	computeSingPointInfo(AFace, s);
	double xs = s.X();
	double ys = s.Y();
	//=========================================================================
	// NOW, WE COMPUTE THE INTERSECTIONS POINT ALONG EACH EDGE
	//=========================================================================
	// for each detected slot, we store its location, its direction, its out cell
	// id and the dimension of this out cell.
	std::vector<math::Point> slot_points;
	std::vector<double> slot_param;
	std::vector<math::Vector3d> slot_dirs;
	std::vector<int> slot_cell_dim;
	std::vector<TCellID> slot_cell_id;
	std::vector<Edge> edges = AFace.get<Edge>();
	for (int i = 0; i < 3; i++) {     // we walk along each edge of AFace
		Edge ei = edges[i];
		std::vector<Node> nodes_ei = ei.get<Node>();
		Node ni = nodes_ei[0];
		Node nj = nodes_ei[1];
		std::vector<math::Vector3d> vectors_i = (*m_field)[ni.id()].componentVectors();
		math::Cross2D cross_j = (*m_field)[nj.id()];

		math::Point pi = ni.getPoint();
		math::Point pj = nj.getPoint();

		double xi = pi.X();
		double yi = pi.Y();
		double xj = pj.X();
		double yj = pj.Y();

		double xji = xi - xj;
		double yji = yi - yj;

		double xsj = xj - xs;
		double ysj = yj - ys;

		for (int k = 0; k < 2; k++) {
			math::Vector3d cik = vectors_i[k];
			math::Vector3d cjk = cross_j.closestComponentVector(cik);

			double x_cik = cik.X();
			double y_cik = cik.Y();
			double x_cjk = cjk.X();
			double y_cjk = cjk.Y();

			double x_cijk = x_cik - x_cjk;
			double y_cijk = y_cik - y_cjk;

			double a = (xji * y_cijk) - (yji * x_cijk);
			double b = (y_cjk * xji) + (xsj * y_cijk) - (x_cjk * yji) - (ysj * x_cijk);
			double c = (y_cjk * xsj) - (x_cjk * ysj);

			std::vector<double> solutions;
			math::solve2ndDegreePolynomial(a, b, c, solutions);
			for (unsigned int i_sol = 0; i_sol < solutions.size(); i_sol++) {
				double alpha = solutions[i_sol];
				if (alpha > 1 || alpha < 0.0) continue;
				math::Point p = alpha * pi + (1 - alpha) * pj;

				slot_param.push_back(alpha);
				slot_points.push_back(p);

				slot_dirs.push_back(math::Vector3d(s, p));

				if (alpha == 0) {     // we go out from a node
					slot_cell_dim.push_back(0);
					slot_cell_id.push_back(ni.id());
				}
				else if (alpha == 0) {     // we go out from a node
					slot_cell_dim.push_back(0);
					slot_cell_id.push_back(nj.id());
				}
				else {     // general case, the edge
					slot_cell_dim.push_back(1);
					slot_cell_id.push_back(ei.id());
				}
			}
		}     // for(int k=0;k<2;k++)
	}        // for(int i=0;i<3;i++)

	// SINGULARITY POINT CREATION
	SurfaceSingularityPoint *singularity = m_graph.newSurfacePoint();
	singularity->setLocation(s);
	singularity->addMeshFace(AFace);
	for (unsigned int i = 0; i < slot_points.size(); i++) {
		singularity->newSlot(slot_points[i], slot_dirs[i], slot_cell_id[i], slot_cell_dim[i], true, 0);
	}
}

/*----------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::addGeometryToSingularityGraph(vector<CurveSingularityPoint *> &artificialSingPointsCreated)
{
	// Now we add singularity points and lines from geometric corners
	//============================================================================
	//	GEOMETRIC  SINGULARITY POINTS
	//============================================================================
	for (auto n_id : m_mesh->nodes()) {
		const Node &currentNode = m_mesh->get<Node>(n_id);
		if (m_mesh->isMarked(currentNode, m_mark_nodes_on_point)) {
			VertexSingularityPoint *sing_point = m_graph.newVertexPoint();
			sing_point->setXYZ(currentNode.X(), currentNode.Y(), currentNode.Z());
			sing_point->addMeshNode(currentNode);
		}
	}
	//============================================================================
	//	GEOMETRIC  SINGULARITY LINES
	//============================================================================
	// Now we have all the geometric corners viewed as singularity points
	// We create separatrices based on the geometric edges

	int mark_geom_edges = m_mesh->newMark<Node>();
	std::vector<CurveSingularityLine *> added_geom_lines;

	for (auto n_id : m_mesh->nodes()) {
		const Node &currentNode = m_mesh->get<Node>(n_id);

		if ((!m_mesh->isMarked(currentNode, m_mark_nodes_on_point)) && (m_mesh->isMarked(currentNode, m_mark_nodes_on_curve))
		    && (!m_mesh->isMarked(currentNode, mark_geom_edges))) {

			// new singularity line to create here
			std::vector<int> listOfNodesInSingLeft;
			std::vector<int> listOfNodesInSingRight;
			std::vector<Edge> edgesLeft;
			std::vector<Edge> edgesRight;

			m_mesh->mark(currentNode, mark_geom_edges);
			Node NodeLeft;
			Node NodeRight;
			bool gotFirstOne = false;
			for (const auto &edge_i : currentNode.get<Edge>()) {
				if (m_mesh->isMarked(edge_i, m_mark_edges_on_curve)) {
					if (!gotFirstOne) {
						gotFirstOne = true;
						NodeLeft = m_tool.getOpposedNodeOnEdge(currentNode, edge_i);
						edgesLeft.push_back(edge_i);
					}
					else {
						NodeRight = m_tool.getOpposedNodeOnEdge(currentNode, edge_i);
						edgesRight.push_back(edge_i);
						break;
					}
				}
			}
			// Here we have the two nodes left and right that continue the line
			//------------------------------------------------------------------------------------
			// first while-loop condition prevent from closing the curl
			Node Ncurr = currentNode;
			while ((NodeLeft.id() != NodeRight.id()) && (!m_mesh->isMarked(NodeLeft, m_mark_nodes_on_point))) {

				m_mesh->mark(NodeLeft, mark_geom_edges);
				listOfNodesInSingLeft.push_back(NodeLeft.id());

				Node NodeNext = NodeLeft;
				for (const auto &edge_i : NodeLeft.get<Edge>()) {
					if (m_mesh->isMarked(edge_i, m_mark_edges_on_curve)) {
						for (const auto &node_j : edge_i.get<Node>()) {
							if ((node_j.id() != NodeLeft.id()) && (node_j.id() != Ncurr.id())) {
								NodeNext = node_j;
								edgesLeft.push_back(edge_i);
							}
						}
					}
				}
				Ncurr = NodeLeft;
				NodeLeft = NodeNext;
			}
			m_mesh->mark(NodeLeft, mark_geom_edges);     // currently this could be marked as node on point
			listOfNodesInSingLeft.push_back(NodeLeft.id());

			if (NodeLeft.id() == NodeRight.id()) {
				listOfNodesInSingRight.push_back(NodeLeft.id());     // cycle
			}
			else {     // now we need to go on the right
				Node Ncurr = currentNode;
				while (!m_mesh->isMarked(NodeRight, m_mark_nodes_on_point)) {

					m_mesh->mark(NodeRight, mark_geom_edges);
					listOfNodesInSingRight.push_back(NodeRight.id());

					Node NodeNext = NodeRight;
					for (const auto &edge_i : NodeRight.get<Edge>()) {
						if (m_mesh->isMarked(edge_i, m_mark_edges_on_curve)) {
							for (const auto &node_j : edge_i.get<Node>()) {
								if ((node_j.id() != NodeRight.id()) && (node_j.id() != Ncurr.id())) {
									NodeNext = node_j;
									edgesRight.push_back(edge_i);
								}
							}
						}
					}
					Ncurr = NodeRight;
					NodeRight = NodeNext;
				}
				m_mesh->mark(NodeRight, mark_geom_edges);
				listOfNodesInSingRight.push_back(NodeRight.id());
			}
			//------------------------------------------------------------------------------------
			// create the singularity line
			//------------------------------------------------------------------------------------
			CurveSingularityLine *new_line = m_graph.newCurveLine();     //  [ reversedLeftNodes..., currentNode, rightNodes...]

			// add all edges
			std::vector<Edge> curve_edges;
			curve_edges.insert(curve_edges.end(), edgesLeft.rbegin(), edgesLeft.rend());
			curve_edges.insert(curve_edges.end(), edgesRight.begin(), edgesRight.end());
			new_line->setMeshEdges(curve_edges);

			// put all nodes into one vector
			std::vector<TCellID> curve_nodes;
			curve_nodes.insert(curve_nodes.end(), listOfNodesInSingLeft.rbegin(), listOfNodesInSingLeft.rend());
			curve_nodes.emplace_back(currentNode.id());
			curve_nodes.insert(curve_nodes.end(), listOfNodesInSingRight.begin(), listOfNodesInSingRight.end());

			// mark forbidden edges
			for (int i = 0; i < curve_nodes.size() - 1; ++i) {
				const TCellID id1 = curve_nodes[i];
				const TCellID id2 = curve_nodes[i + 1];
				if (m_mesh->isMarked<Node>(id1, m_mark_forbiddenBdryNode) && m_mesh->isMarked<Node>(id2, m_mark_forbiddenBdryNode)) {
					const TCellID edge = curve_edges[i].id();
					m_mesh->mark<Edge>(edge, m_mark_forbiddenBdryEdge);
				}
			}

			// add to node location to line discretization
			for (const TCellID id : curve_nodes) {
				Node node = m_mesh->get<Node>(id);
				new_line->addDiscretizationPoint(node.getPoint());
			}

			added_geom_lines.push_back(new_line);
		}
	}
	// Geometrical curves made of only one mesh edge are missing. We add them now.
	for (auto e_id : m_mesh->edges()) {
		const Edge &currentEdge = m_mesh->get<Edge>(e_id);
		// on ne traite que les aretes sur une courbe geometrique
		if (m_mesh->isMarked(currentEdge, m_mark_edges_on_curve)) {
			const std::vector<Node> &currentNodes = currentEdge.get<Node>();
			const Node &n1 = currentNodes[0];
			const Node &n2 = currentNodes[1];
			// on regarde si les deux noeuds correspondent ˆ des sommets geometriques
			if (m_mesh->isMarked(n1, m_mark_nodes_on_point) && m_mesh->isMarked(n2, m_mark_nodes_on_point)) {     // On cree donc la separatrice reliant les
				                                                                                                   // singularites associŽes ˆ n1 et n2

				CurveSingularityLine *new_line = m_graph.newCurveLine();
				std::vector<Edge> edges {currentEdge};
				new_line->setMeshEdges(edges);
				new_line->addDiscretizationPoint(n1.getPoint());
				new_line->addDiscretizationPoint(n2.getPoint());
				added_geom_lines.push_back(new_line);

				if (m_mesh->isMarked(n1, m_mark_forbiddenBdryNode) && m_mesh->isMarked(n2, m_mark_forbiddenBdryNode)) {
					m_mesh->mark<Edge>(e_id, m_mark_forbiddenBdryEdge);
				}
			}
		}
	}
	//============================================================================
	// Connection of singularity lines to singularity points
	//============================================================================
	std::vector<VertexSingularityPoint *> geom_points = m_graph.getVertexPoints();

	for (VertexSingularityPoint *pi : geom_points) {
		const Node &ni = pi->getMeshNode();
		for (CurveSingularityLine *lj : added_geom_lines) {

			std::vector<Edge> &lj_edges = lj->getMeshEdges();
			Node first_node, last_node;
			math::Vector3d first_vec, last_vec;
			if (lj_edges.size() == 1) {
				std::vector<Node> e_nodes = lj_edges[0].get<Node>();
				first_node = e_nodes[0];
				last_node = e_nodes[1];
				first_vec = math::Vector3d(first_node.getPoint(), last_node.getPoint());
				last_vec = math::Vector3d(last_node.getPoint(), first_node.getPoint());
			}
			else {
				Node nextNode = m_tool.getCommonNode(lj_edges[0], lj_edges[1]);
				first_node = m_tool.getOpposedNodeOnEdge(nextNode, lj_edges[0]);
				first_vec = math::Vector3d(first_node.getPoint(), nextNode.getPoint());

				Edge &last_but_one_edge = lj_edges[lj_edges.size() - 2];
				Edge &last_edge = lj_edges[lj_edges.size() - 1];

				Node penultiemNode = m_tool.getCommonNode(last_but_one_edge, last_edge);
				last_node = m_tool.getOpposedNodeOnEdge(penultiemNode, last_edge);
				last_vec = math::Vector3d(last_node.getPoint(), nextNode.getPoint());
			}
			// We compare the point node with the first and last node of the curve
			if (ni == first_node) {
				// pi and lj must be connected
				// WARNING SLOT DIMENSION IS DETERMINED BY THE LAST EDGE DIRECTION OF THE CURVE LINE STARTING FROM THIS
				// POINT;
				// PERHAPS THE ORIENTATION OF THE TOTAL CURVE LINE SHOULD BE TAKEN INTO ACC

				SingularityPoint::Slot *s = pi->newSlot(ni.getPoint(), first_vec, ni.id() /*starting cell id*/, 0 /*starting cell dim*/, true /*on surface*/, lj);
				s->isLaunched = true;
				s->isFreeze = true;
				lj->addSlot(s);
			}
			else if (ni == last_node) {
				// pi and lj must be connected
				// pi->newSlot(ni.getPoint(), last_vec,
				SingularityPoint::Slot *s = pi->newSlot(ni.getPoint(), last_vec, ni.id() /*starting cell id*/, 0 /*starting cell dim*/, true /*on surface*/, lj);
				s->isLaunched = true;
				s->isFreeze = true;
				lj->addSlot(s);
			}
		}
	}
	//============================================================================
	// Cycling boundary curves need to be connected to a singularity point.
	// Those will be created as "articial" and deleted at the end.
	//============================================================================
	std::vector<CurveSingularityLine *> geom_curves = m_graph.getCurveLines();

	for (CurveSingularityLine *ci : geom_curves) {
		if (ci->getSlots().empty()) {

			const std::vector<math::Point> &line_disc = ci->getDiscretizationPoints();
			const math::Point &p_begin = line_disc[0];
			const math::Point &p_end = line_disc[line_disc.size() - 1];
			const math::Vector3d slot_dir1 = math::Vector3d(p_begin, line_disc[1]);
			const math::Vector3d slot_dir2 = math::Vector3d(p_end, line_disc[line_disc.size() - 2]);

			CurveSingularityPoint *p = m_graph.newCurvePoint();
			artificialSingPointsCreated.push_back(p);
			p->setLocation(p_begin);

			SingularityPoint::Slot *s = p->newSlot(p->getLocation(),     // slot location
			                                       slot_dir1,            // slot direction
			                                       0,                    // No linked cell (id)
			                                       0,                    // No linked cell (dim)
			                                       true,                 // Always on surface (Maybe false in the future)
			                                       ci,                   // Connected line
			                                       slot_dir1.opp());     // Line direction is the slot direction
			s->isFreeze = true;
			ci->addSlot(s);

			s = p->newSlot(p->getLocation(), slot_dir2, 0, 0, true, ci, slot_dir2.opp());
			s->isFreeze = true;
			ci->addSlot(s);
		}
	}

	// slot must be in same order as discretization
	for (CurveSingularityLine *lj : added_geom_lines) {
		auto &slots = lj->getSlots();
		if (slots.front()->location != lj->getDiscretizationPoints().front()) std::reverse(slots.begin(), slots.end());
	}

	m_mesh->unmarkAll<Node>(mark_geom_edges);
	m_mesh->freeMark<Node>(mark_geom_edges);
}

/*----------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::buildSlotsOfBoundarySingularities()
{
	//================================================================================
	// Create the surface slots of vertex singularities (boundary singularity)
	//================================================================================)
	// Consider the boundaries of this potential rectangular trimesh, with a hole
	// and a pointy deformation at the bottom.
	//
	//  A##########################B
	//  #						 \__#     The nodes ABCDEFGHIJKL have all been marked as singularities.
	//	#   G######H           90°  #     They all have at least two boundary slots.
	//	#   #      #      			L	  (Here for example, the nodes ABCDEFGHIJK could have been marked
	//	#   # hole #     			#	   with gmds::BoundaryOperator class, and the node L manually)
	//	#   #	   #	    K		#
	//	#   I######J	   # #		#     To know how many slots they need, we need to compute the
	//  #                  # #		#	  solid angle made by the boundary at their location.
	//  C#################E   F####D
	//
	// - the nodes ABCDE, and F have a "solid angle" ~= 90° -> they are corners of the cross field,
	//	 no surface slot will be created for them
	// - GHIJ solid angles are ~= 270°, they need at least two surface slots to respect the cross field
	// - the node K is a spike on the boundary, its angle solid angle is high an it needs 3 slots
	// - the node L, surely marked by user, is on a flat boundary, and need only 1 surface slot
	//
	// However, we can't always be sure about the precise number of slots. Their validity will be check
	// right after they creation.

	for (VertexSingularityPoint *singularity : m_graph.getVertexPoints()) {

		Node currentNode = singularity->getMeshNode();
		math::Point currentPoint = currentNode.getPoint();

		// compute the solid angle around currentNode
		double angle_rad = 0;
		std::vector<math::Point> opposedPoint(2);
		for (const auto &cur_face : currentNode.get<Face>()) {

			int pointID = 0;
			for (const auto &ni : cur_face.get<Node>()) {
				if (ni.id() != currentNode.id()) {
					opposedPoint[pointID++] = ni.getPoint();
				}
			}
			math::Vector3d v1(currentPoint, opposedPoint[0]);
			math::Vector3d v2(currentPoint, opposedPoint[1]);
			angle_rad += v1.angle(v2);
		}
		double angle_deg = angle_rad * 180 / math::Constants::PI;

		// build slots. If the built slots appear to not respect the cross field,
		// buildSurfaceSlotsOfVertexSingularity() will call itself recursively, with 1 slot less to build.
		if (angle_deg > 275) {
			int nb_lines = 3;
			double single_line_angle = angle_rad / (nb_lines + 1);
			buildSlotsOfBoundarySingularity(singularity, single_line_angle, nb_lines);
		}
		else if (angle_deg > 180) {
			int nb_lines = 2;
			double single_line_angle = angle_rad / (nb_lines + 1);
			buildSlotsOfBoundarySingularity(singularity, single_line_angle, nb_lines);
		}
		else if (angle_deg > 125) {
			int nb_lines = 1;
			double single_line_angle = angle_rad / (nb_lines + 1);
			buildSlotsOfBoundarySingularity(singularity, single_line_angle, nb_lines);
		}
	}
}

/*----------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::buildSlotsOfBoundarySingularity(VertexSingularityPoint *AFrom, const double AAngle, const int ANbLines)
{
	const Node currentNode = AFrom->getMeshNode();

	// We get a first edge to start from
	const Edge firstBdryEdge = [&]() {
		for (const Edge &edge : currentNode.get<Edge>())
			if (m_mesh->isMarked(edge, m_mark_edges_on_curve)) return edge;
		throw std::runtime_error("unable to find boundary edge");
	}();

	const Node bdryNode = Tools::getOpposedNodeOnEdge(currentNode, firstBdryEdge);
	const Node surfNode = [&]() {
		const Face first_face = firstBdryEdge.get<Face>()[0];
		for (const Node &node : first_face.get<Node>())
			if (node.id() != currentNode.id() && node.id() != bdryNode.id()) return node;
		throw std::runtime_error("unable to find surface node");
	}();

	const auto firstBdryVector = math::Vector3d(currentNode.getPoint(), bdryNode.getPoint()).normalize();
	const auto surfaceVector = math::Vector3d(currentNode.getPoint(), surfNode.getPoint()).normalize();
	const auto axis = firstBdryVector.cross(surfaceVector).normalize();

	struct BaseSlot
	{
		math::Point location;
		math::Vector3d direction;
		TCellID cellId;
		int cellDim;     // a slot can be on a node or an edge
	};
	std::vector<BaseSlot> baseSlots;
	for (int i_line = 0; i_line < ANbLines && baseSlots.size() != ANbLines; i_line++) {

		const double lineAngle = AAngle * (i_line + 1);
		const math::AxisAngleRotation R(axis, lineAngle);
		const math::Vector3d lineDirection = (R * firstBdryVector).normalize();

		// Now for each face adjacent to currentNode, we look for an intersection with lineDirection
		for (const auto &currentFace : currentNode.get<Face>()) {

			const Edge oppositEdge = [&]() {
				for (const Edge &edge : currentFace.get<Edge>()) {
					const auto nodes = edge.get<Node>();
					if (currentNode != nodes[0] && currentNode != nodes[1]) return edge;
				}
				throw std::runtime_error("unable to find opposite edge");
			}();

			math::Vector3d slotDirection;     // = direction imposed by the local cross field

			// When really close to a node, its better to place the slot right on it rather than the edge
			// @ {
			const auto other_nodes = oppositEdge.get<Node>();
			if (!m_mesh->isMarked(other_nodes[0], m_mark_nodes_on_curve)) {

				const auto opp_node_loc1 = other_nodes[0].getPoint();
				const auto v_opp1 = math::Vector3d(currentNode.getPoint(), opp_node_loc1).normalize();
				if (v_opp1.dot(lineDirection) > 0.99) {     // 8°
					m_tool.computeOutVectorAtPoint(other_nodes[0], lineDirection, slotDirection);
					baseSlots.push_back({opp_node_loc1, slotDirection, other_nodes[0].id(), 0});
					break;
				}
			}
			if (!m_mesh->isMarked(other_nodes[1], m_mark_nodes_on_curve)) {

				const auto opp_node_loc2 = other_nodes[1].getPoint();
				const auto v_opp2 = math::Vector3d(currentNode.getPoint(), opp_node_loc2).normalize();
				if (v_opp2.dot(lineDirection) > 0.99) {     // 8°
					m_tool.computeOutVectorAtPoint(other_nodes[1], lineDirection, slotDirection);
					baseSlots.push_back({opp_node_loc2, slotDirection, other_nodes[1].id(), 0});
					break;
				}
			}
			// @ }

			// search for an intersection with the opposite edge
			double deviation = 0;
			math::Point slotPoint;
			if (m_tool.computeOutVectorFromRayAndEdge(oppositEdge, currentNode.getPoint(), lineDirection, slotPoint, slotDirection, deviation)) {
				baseSlots.push_back({slotPoint, slotDirection, oppositEdge.id(), 1});
				break;
			}
		}
	}

	if (baseSlots.size() != ANbLines) {
		throw GMDSException("unable to create the required number of slot," + std::to_string(baseSlots.size()) + " instead of " + std::to_string(ANbLines));
	}
	if (ANbLines == 3) {     // building 3 slots might have been a mistake if the cross field is grid like: (at least two slots with close directions)
		bool recompute = baseSlots[0].direction.angle(baseSlots[1].direction) < M_PI_4     //
		                 || baseSlots[1].direction.angle(baseSlots[2].direction) < M_PI_4;
		if (recompute) return buildSlotsOfBoundarySingularity(AFrom, AAngle * 4 / 3, 2);     // create 2 slots instead
	}
	else if (ANbLines == 2) {     // building 2 slots might also have been a mistake (close slots directions)
		bool recompute = baseSlots[0].direction.angle(baseSlots[1].direction) < M_PI_4;
		if (recompute) return buildSlotsOfBoundarySingularity(AFrom, AAngle * 3 / 2, 1);
	}
	else if (ANbLines == 1) {     // same for bulding only 1 slot: the singularity might be a corner of the cross field
		const math::AxisAngleRotation R(axis, 2 * AAngle);
		const math::Vector3d secondBdryVector = (R * firstBdryVector).normalize();
		bool doNotBuildSlot = baseSlots[0].direction.angle(firstBdryVector) < M_PI_4 || baseSlots[0].direction.angle(secondBdryVector) < M_PI_4;
		if (doNotBuildSlot) return;
	}
	for (const auto baseSlot : baseSlots)
		AFrom->newSlot(baseSlot.location, baseSlot.direction, baseSlot.cellId, baseSlot.cellDim, true);
}

/*----------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::createGeometricSingularityPoint(
   const math::Point &AInPnt, const math::Vector3d &AInVec, const int ACellDim, const TCellID ACellID, SingularityPoint *&APnt, SingularityPoint::Slot *&ASlot)
{
	if (ACellDim == 0) {

		Node n = m_mesh->get<Node>(ACellID);

		if (m_mesh->isMarked(n, m_mark_nodes_on_point)) {
			/*We arrive onto a geometric point !!!! It means that a singularity geometric point already exist*/
			for (VertexSingularityPoint *p : m_graph.getVertexPoints())
				if (p->getMeshNode().id() == ACellID) throw GMDSException("unimplemented scenario!");

			throw GMDSException("createGeometricSingularityPoint: Error, no geometric point to be connected to!");
		}
		else
			m_graph.splitCurveLine(AInPnt, AInVec, n, APnt, ASlot);
	}
	else {
		Edge e = m_mesh->get<Edge>(ACellID);
		m_graph.splitCurveLine(AInPnt, AInVec, e, APnt, ASlot);
	}
}

/*----------------------------------------------------------------------------*/

std::vector<SurfaceSingularityLine *>
SingularityGraphBuilder2D::getSingularityLinesIn(const Face &AFace)
{
	/*get the singularity lines that traverse triangle AFace*/
	std::vector<SurfaceSingularityLine *> found_lines;

	std::vector<SurfaceSingularityLine *> lines = m_graph.getSurfaceLines();

	TCellID face_id = AFace.id();
	for (unsigned int i = 0; i < lines.size(); i++) {
		SurfaceSingularityLine *line_i = lines[i];
		if (line_i->isTraversed(face_id)) {
			found_lines.push_back(line_i);
		}
	}

	return found_lines;
}

/*----------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::writeOutput(const std::string &AFileName)
{
	static int out = 0;
	std::string file_name = AFileName + "_" + std::to_string(out);

	writeOutputSingle(file_name);
	out++;
}

/*----------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::writeOutputSingle(const std::string &AFileName)
{
	std::string file_name = m_output_directory_name + "-" + AFileName + ".vtk";

	m_graph.createVTKOutputFile(file_name);
}

/*----------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::writeOutputPatches(const std::string &AFileName)
{
	std::string file_name = m_output_directory_name + "-" + AFileName + ".vtk";
	bool curvePatches = false;
	m_graph.createVTKOutputFile(file_name, curvePatches);
}

/*-----------------------------------------------------------------*/

void
SingularityGraphBuilder2D::writeSingularityPointsAndSlots()
{
	Mesh m(MeshModel(DIM3 | F | N | F2N));
	auto slotIDs = m.newVariable<int, GMDS_FACE>("slotIDs");

	int singID = 0;
	for (auto sing : m_graph.getPoints()) {
		int slotID = 0;
		for (auto slot : sing->getSlots()) {

			Node n0 = m.newNode(sing->getLocation());
			Node n1 = m.newNode(slot->location);
			Node n2 = m.newNode(slot->location + slot->direction.normalize() * m_mean_edge_length);

			auto f1 = m.newTriangle(n0, n1, n1);
			auto f2 = m.newTriangle(n1, n2, n2);

			int id = 5 * singID + slotID;
			(*slotIDs)[f1.id()] = id;
			(*slotIDs)[f2.id()] = id;
			++slotID;
		}
		++singID;
	}

	IGMeshIOService meshIoServ(&m);
	VTKWriter writer(&meshIoServ);
	writer.setCellOptions(N | F);
	writer.setDataOptions(N | F);

	writer.write(m_output_directory_name + "-points_and_slots.vtk");
}

/*----------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::writeConfusingBalls()
{
	// variable for debug purpose
	Variable<int> *ball_var = m_mesh->newVariable<int, GMDS_FACE>("sing_ball");
	for (auto f_id : m_mesh->faces()) {
		Face f = m_mesh->get<Face>(f_id);
		SingularityPoint *sing = m_faces_to_singularity_on_surf[f.id()];
		if (sing == 0)
			(*ball_var)[f.id()] = 0;
		else
			(*ball_var)[f.id()] = sing->index();
	}

	IGMeshIOService meshIoServ(m_mesh);
	VTKWriter writerB(&meshIoServ);
	writerB.setCellOptions(N | F);
	writerB.setDataOptions(N | F);
	std::stringstream file_name2;
	file_name2 << m_output_directory_name << "-confusing_balls.vtk";
	writerB.write(file_name2.str());
}

/*----------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::redefineOneConfusingBall(SingularityPoint *APnt,
                                                    vector<TCellID> &modifiedFaces,
                                                    double &previousDist,
                                                    const double increaseRadiusScale)
{
	if (m_withGlobalComments) cout << "redefineOneConfusingBall " << endl;
	// here it could be a problem if two sing are very close geometrically one to another

	unsigned int nbVerts = m_mesh->getNbNodes();

	vector<bool> visitedFaces(m_original_faces_number, false);
	vector<bool> visitedVerts(nbVerts, false);

	math::Point sing_location = APnt->getLocation();
	double newConfusingDistance = increaseRadiusScale * m_confusing_distance;

	std::vector<Face> neighbouringFaces;
	std::vector<Node> vertsToVisit;
	//  cout<<"APnt0>getNbMeshCells() "<<APnt->getNbMeshCells()<<endl;

	vector<Face> AtriVectFace = APnt->getMesh<Face>();

	if (AtriVectFace.size() != 0) {
		Face Atri = AtriVectFace[0];
		if (m_withGlobalComments) cout << "face " << Atri.id() << endl;
		std::vector<Node> nodes = Atri.get<Node>();
		for (unsigned int i = 0; i < nodes.size(); i++)
			vertsToVisit.push_back(nodes[i]);
	}
	else {     // geom point
		vector<Node> AtriVectNode = APnt->getMesh<Node>();

		if (AtriVectNode.size() != 0) {
			Node Atri = AtriVectNode[0];
			vertsToVisit.push_back(Atri);
			if (m_withGlobalComments) cout << "node " << Atri.id() << endl;
		}
		else {
			vector<Edge> AtriVectEdge = APnt->getMesh<Edge>();

			if (AtriVectEdge.size() != 0) {
				Edge Atri = AtriVectEdge[0];
				std::vector<Node> nodes = Atri.get<Node>();
				vertsToVisit.push_back(nodes[0]);
				vertsToVisit.push_back(nodes[1]);
				if (m_withGlobalComments) cout << "edge with nodes " << nodes[0].id() << " " << nodes[1] << endl;
			}
			// else curve sing point -> doesnt need confusing ball for now
		}
	}

	/*
	for(unsigned int i=0; i<originalConfusingBalls.size(); i++){
	   visitedFaces[originalConfusingBalls[i]] = true;
	   Face currentFace = m_mesh->get<Face> (originalConfusingBalls[i]);
	   std::vector<Node> nodes = currentFace.get<Node>();
	      for(unsigned int i=0; i<nodes.size(); i++){
	      visitedVerts[nodes[i].id()] = true;
	      }
	   }
	   */

	for (unsigned int i = 0; i < vertsToVisit.size(); i++) {
		visitedVerts[vertsToVisit[i].id()] = true;
	}

	while (!vertsToVisit.empty()) {
		Node currentNode = vertsToVisit[vertsToVisit.size() - 1];
		vertsToVisit.resize(vertsToVisit.size() - 1);
		currentNode.get<Face>(neighbouringFaces);
		for (unsigned int i = 0; i < neighbouringFaces.size(); i++) {
			if ((!visitedFaces[neighbouringFaces[i].id()])) {
				math::Point center = neighbouringFaces[i].center();
				double current_dist = center.distance(sing_location);
				if (current_dist < newConfusingDistance) {
					if (current_dist >= previousDist) {
						modifiedFaces.push_back(neighbouringFaces[i].id());
						m_faces_to_singularity_on_surf[neighbouringFaces[i].id()] = APnt;
					}
					vector<Node> currentFaceVerts;
					neighbouringFaces[i].get<Node>(currentFaceVerts);     // std::vector<TCellID> currentFaceVerts = neighbouringFaces[i].getIDs<Node>();
					for (unsigned int j = 0; j < currentFaceVerts.size(); j++) {
						if (!visitedVerts[currentFaceVerts[j].id()]) {
							visitedVerts[currentFaceVerts[j].id()] = true;
							vertsToVisit.push_back(currentFaceVerts[j]);
						}
					}
				}
				visitedFaces[neighbouringFaces[i].id()] = 1;
			}
		}
	}

	previousDist = newConfusingDistance;
}

/*----------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::removeSingularityLine(SingularityLine *ALine)
{
	/*
	if the line has one of its end points on the boundary, a geometric singularity
	has been added for it and it has also divided the boundary singularity line;
	therefore we have to remove the point created, merge back the 2 boundary singularity
	lines into one (and add it to the graph) and remove them*/
	vector<SurfaceSingularityLine *> surf_lines = m_graph.getSurfaceLines();
	// it should always be a SurfaceSingularityLine
	bool deletePoint = true;
	bool foundLine = false;
	for (unsigned int i = 0; i < surf_lines.size() && !foundLine; i++) {
		SingularityLine *current_line = surf_lines[i];
		if (current_line == ALine) {
			foundLine = true;
			vector<SingularityPoint *> singPointsOnSingLine = surf_lines[i]->getEndPoints();
			// cout<<"between ["<<singPointsOnSingLine[0]->getLocation().X()<<","<<singPointsOnSingLine[0]->getLocation().Y()<<"]"<<endl;
			// cout<<"and ["<<singPointsOnSingLine[1]->getLocation().X()<<","<<singPointsOnSingLine[1]->getLocation().Y()<<"]"<<endl;
			for (unsigned int t = 0; t < 2; t++) {
				SingularityPoint *current_sp = singPointsOnSingLine[t];

				if (current_sp->getType() == 1) { /* geom singularity point created on boundary (on edge or on vert)
					   if on vert, check its marked as NodeOnPoint (if so, we don't have to do anything)*/
					vector<Node> firstNodes = current_sp->getMesh<Node>();

					if (firstNodes.size() > 0) {
						// cout<<"the end point of the singularity line is a vertex"<<endl;
						if (m_mesh->isMarked(firstNodes[0], m_mark_nodes_on_point)) {
							deletePoint = false;
						}
					}

					if (deletePoint) {
						vector<SingularityPoint *> m_points = m_graph.getPoints();

						for (unsigned int j = 0; j < m_points.size(); j++) {
							if (m_points[j] == current_sp) {
								m_points.erase(m_points.begin() + j);
								break;
							}
						}
						vector<CurveSingularityLine *> curve_lines = m_graph.getCurveLines();
						bool foundCurveLine = false;
						unsigned int firstCurveLineIndex;
						// unsigned int current_spLineExtremity;
						bool foundBothLines = false;

						for (unsigned int j = 0; j < curve_lines.size() && !foundBothLines; j++) {
							vector<SingularityPoint *> singPointsOnCurveSingLine = curve_lines[j]->getEndPoints();
							// this should always be 2-sized (no intersections yet)
							for (unsigned int tt = 0; tt < 2; tt++) {
								if (singPointsOnCurveSingLine[tt] == current_sp) {

									if (foundCurveLine) {
										// this is the second curveline that we find (as well as the number)
										// keep the orientation of the first curve that we find
										// singular point is always the 2nd on the 1st line
										vector<math::Point> first_pnts = curve_lines[firstCurveLineIndex]->getDiscretizationPoints();
										vector<Edge> first_curve_edges = curve_lines[firstCurveLineIndex]->getMeshEdges();
										vector<math::Point> second_pnts = curve_lines[j]->getDiscretizationPoints();
										vector<Edge> second_curve_edges = curve_lines[j]->getMeshEdges();

										if (tt == 1) {
											std::reverse(second_pnts.begin(), second_pnts.end());
											std::reverse(second_curve_edges.begin(), second_curve_edges.end());
										}
										second_pnts.erase(second_pnts.begin());
										second_curve_edges.erase(second_curve_edges.begin());
										for (unsigned int ttt = 0; ttt < second_pnts.size(); ttt++)
											first_pnts.push_back(second_pnts[ttt]);

										for (unsigned int ttt = 0; ttt < second_curve_edges.size(); ttt++)
											first_curve_edges.push_back(second_curve_edges[ttt]);

										curve_lines[firstCurveLineIndex]->setMeshEdges(first_curve_edges);
										curve_lines[firstCurveLineIndex]->setDiscretizationPoints(first_pnts);

										curve_lines[firstCurveLineIndex]->removeSlot();
										curve_lines[firstCurveLineIndex]->addSingularityPoint(singPointsOnCurveSingLine[(tt + 1) % 2]);
										foundBothLines = true;
										m_graph.removeCurveLine(curve_lines[j]);
										break;
									}

									firstCurveLineIndex = j;

									foundCurveLine = true;
								}
							}
						}
					}
					break;
				}
			}
		}
	}

	if (!foundLine) {
		if (m_withGlobalComments) cout << "ALine is not a surface line" << endl;
	}
}

/*--------------------------------------------------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::constructCenterTriangleCrosses()
{
	m_triangle_centers.clear();
	m_triangle_centers.reserve(m_original_faces_number);
	m_triangle_centers_cross.clear();
	m_triangle_centers_cross.reserve(m_original_faces_number);

	for (int f_id = 0; f_id < m_original_faces_number; f_id++) {

		const auto &currentFace = m_mesh->get<Face>(f_id);
		m_triangle_centers.emplace_back(currentFace.center());

		const auto &nodesTri = currentFace.get<Node>();
		const auto &cross1 = (*m_field)[nodesTri[0].id()];
		const auto &cross2 = (*m_field)[nodesTri[1].id()];
		const auto &cross3 = (*m_field)[nodesTri[2].id()];

		const math::Vector3d vect1 = cross1.componentVectors()[0];
		const math::Vector3d vect2 = cross2.closestComponentVector(vect1);
		const math::Vector3d vect3 = cross3.closestComponentVector(vect1);

		math::Vector3d center_vect = (vect1 + vect2 + vect3) / 3;
		math::Vector3d center_vect_second = center_vect.getOneOrtho();

		m_triangle_centers_cross.emplace_back(math::Cross2D(center_vect, center_vect_second));
	}
}

/*--------------------------------------------------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::constructOneCenterTriangleCross(Face &ATriangle, math::Point &singleCenterTri, math::Cross2D &singleCenterTriCross)
{
	math::Point orig(0.0, 0.0, 0.0);

	std::vector<math::Cross2D> Tricrosses;
	vector<TCoord> AWeights(3, (double) 1 / 3);

	vector<Node> nodesTri = ATriangle.get<Node>();
	Tricrosses.push_back((*m_field)[nodesTri[0].id()]);
	Tricrosses.push_back((*m_field)[nodesTri[1].id()]);
	Tricrosses.push_back((*m_field)[nodesTri[2].id()]);
	vector<math::Vector3d> compV0 = (*m_field)[nodesTri[0].id()].componentVectors();
	vector<math::Vector3d> compV1 = (*m_field)[nodesTri[1].id()].componentVectors();
	vector<math::Vector3d> compV2 = (*m_field)[nodesTri[2].id()].componentVectors();

	singleCenterTri = ATriangle.center();

	TCoord ref_angle = Tricrosses[0].referenceAngle();

	math::Vector3d ref_vector = Tricrosses[0].referenceVector();

	TCoord pen_angle = 0.0;
	for (unsigned int i = 0; i < 3; i++) {
		const TCoord pen_i = Tricrosses[i].referenceAngle();
		pen_angle += pen_i;
	}
	pen_angle /= Tricrosses.size();

	ref_angle = math::modulo2PI(pen_angle);
	ref_vector = math::Cross2D(ref_angle).referenceVector();

	singleCenterTriCross = math::Cross2D(ref_angle);
}

/*---------------------------------------------------------------------------*/
void
SingularityGraphBuilder2D::growLine(SingularityPoint *AFromSingPnt,
                                    SingularityPoint::Slot *AFromSlot,
                                    SingularityPoint *&AToSingPnt,
                                    SingularityPoint::Slot *&AToSlot,
                                    math::Point &AFromPnt,
                                    math::Vector3d &AToDir,
                                    std::vector<math::Point> &APoints,
                                    std::vector<TCellID> &ATriangles,
                                    int &AToCellDim,
                                    TCellID &AToCellID,
                                    double &streamlineDeviation,
                                    double &accumulatedDistancePerSlotLine,
                                    double &currentSearchStep,
                                    bool &AEndOnBnd,
                                    bool &AToSlotIsFree)
{

	bool find_end = false;
	//========================================================================
	/* grow a line starting from the singularity point AFromSingPnt through the slot AFromSlot
	as long as the geometric distance travelled by this line (accumulatedDistancePerSlotLine)
	is smaller than currentSearchStep or until we meet a boundary vertex */
	//========================================================================
	SingularityPoint *found_pnt = 0;
	SingularityPoint::Slot *found_slot = 0;

	math::Point start_pnt = AFromPnt;     // starting point
	// APoints.push_back(start_pnt);
	math::Vector3d start_dir = AToDir;     // starting direction
	math::Vector3d prev_dir = AToDir;      // prev direction used in the
	TCellID start_cell_id = AToCellID;
	int start_cell_dim = AToCellDim;

	math::Point current_pnt = start_pnt;
	math::Vector3d current_vec = start_dir;

	if (start_cell_dim == 0)
		cout << "node: " << start_cell_id << endl;
	else if (start_cell_dim == 1) {
		vector<Node> currentNodes = (m_mesh->get<Edge>(start_cell_id)).get<Node>();
	}

	math::Point start_dirPnt(start_pnt.X() + start_dir.X(), start_pnt.Y() + start_dir.Y(), start_pnt.Z() + start_dir.Z());

	// We check some termination conditions on the boundary.
	if (start_cell_dim == 0) {
		Node currentNode = m_mesh->get<Node>(start_cell_id);
		if (m_mesh->isMarked(currentNode, m_mark_nodes_on_point) || m_mesh->isMarked(currentNode, m_mark_nodes_on_curve)) {
			AEndOnBnd = true;
		}
	}
	else {     // we have necessarry start_cell_dim=1
		Edge currentEdge = m_mesh->get<Edge>(start_cell_id);
		if (m_mesh->isMarked(currentEdge, m_mark_edges_on_curve)) {
			AEndOnBnd = true;
		}
	}

	while ((accumulatedDistancePerSlotLine < currentSearchStep) && (!AEndOnBnd) && (!find_end)) {
		TCellID next_cell_id = NullID;
		int next_cell_dim = -1;
		m_tool.findNextCell(start_pnt, start_dir, start_cell_dim, start_cell_id, next_cell_dim, next_cell_id);

		if (next_cell_dim == -1) {
			// The cell defined by (start_cell_dim, start_cell_id) is on the boundary.
			AEndOnBnd = true;
		}
		else if (next_cell_dim == 1) {
			cout << "next_cell_dim == 1" << endl;
			// we are going along an edge.
			// Our simple assumption is to follow this edge until reaching
			// one of its end points and to compute the next direction at
			// this point.
			Edge currentEdge = m_mesh->get<Edge>(next_cell_id);
			std::vector<TCellID> adj_faces = currentEdge.getIDs<Face>();
			ATriangles.insert(ATriangles.end(), adj_faces.begin(), adj_faces.end());

			std::vector<Node> currentNodes = currentEdge.get<Node>();
			math::Vector3d v0(start_pnt, currentNodes[0].getPoint());
			math::Vector3d v1(start_pnt, currentNodes[1].getPoint());
			Node next_node;
			if (math::near(v0.norm(), 0.0))
				next_node = currentNodes[1];
			else if (math::near(v1.norm(), 0.0))
				next_node = currentNodes[0];
			else if (v0.dot(start_dir) > v1.dot(start_dir))
				next_node = currentNodes[0];
			else
				next_node = currentNodes[1];

			math::Vector3d next_dir;
			math::Vector3d devVect(start_pnt, next_node.getPoint());
			devVect.normalize();
			m_tool.computeOutVectorAtPoint(next_node, start_dir, next_dir);
			start_dir.normalize();
			streamlineDeviation = streamlineDeviation + fabs(1.0 - start_dir.dot(devVect));

			// We assign the new value for the next step
			start_dir = next_dir;
			start_pnt = next_node.getPoint();
			start_cell_dim = 0;
			start_cell_id = next_node.id();
			accumulatedDistancePerSlotLine = accumulatedDistancePerSlotLine + start_pnt.distance(APoints[APoints.size() - 1]);
			APoints.push_back(start_pnt);
		}
		else {
			// general case, we are in a face
			Face currentFace = m_mesh->get<Face>(next_cell_id);
			ATriangles.push_back(currentFace.id());
			//==============================================================
			// CASE 1: ARE WE IN A FACE CONTAINING A SING. POINT?
			// VERY IMPROBABLE; GOT HERE BUT WE WEREN'T WITHIN THRESHOLD
			//==============================================================
			bool intersect_sing_point = false;

			SingularityPoint *next_sing_point = m_faces_to_singularity_on_surf[currentFace.id()];
			bool must_try_to_connect = false;
			if (next_sing_point != NULL &&             // face in a confusing ball ...
			    next_sing_point != AFromSingPnt) {     // ... of another singularity point
				must_try_to_connect = true;
			}
			else {
				if (next_sing_point != NULL &&             // face in the confusing ball ...
				    next_sing_point == AFromSingPnt) {     //... of the incoming singularity point
					cout << "face in a confusing ball of the incoming singularity point" << endl;
					// Warning: completly empiric, we just try to detect cyclic lines
					if (APoints.size() >= 100) must_try_to_connect = true;
				}
			}

			if (must_try_to_connect) {
				cout << "must_try_to_connect" << endl;
				math::Point start_dirPnt(start_pnt.X() + start_dir.X(), start_pnt.Y() + start_dir.Y(), start_pnt.Z() + start_dir.Z());

				// Now, we look for a compatible slot
				std::vector<SingularityPoint::Slot *> &cur_slots = next_sing_point->getSlots();
				//==============================================================
				// We look for a free slot and we connect the line to it
				//==============================================================
				bool found_compatible_slot = false;
				bool found_free_slot = false;
				found_pnt = next_sing_point;
				double slot_epsilon = 0.9;
				// we normalize the vector we arrive with
				math::Vector3d current_vec = start_dir;
				current_vec.normalize();
				while (!found_compatible_slot && slot_epsilon > 0.4) {
					double best_deviation = -2;
					double best_slot_id = 0;

					for (unsigned int i_slot = 0; i_slot < cur_slots.size(); i_slot++) {
						SingularityPoint::Slot *current_slot = cur_slots[i_slot];
						if (current_slot->isFreeze) continue;
						math::Vector3d slot_opp_dir = current_slot->direction.opp();
						slot_opp_dir.normalize();
						double slot_deviation = slot_opp_dir.dot(current_vec);

						if (slot_deviation > slot_epsilon && slot_deviation > best_deviation) {
							best_deviation = slot_deviation;
							best_slot_id = i_slot;
						}
					}
					if (best_deviation != -2 && !cur_slots.empty()) {
						SingularityPoint::Slot *best_slot = cur_slots[best_slot_id];
						math::Vector3d slot_opp_dir = best_slot->direction.opp();
						slot_opp_dir.normalize();
						if (best_slot->isLaunched) {
							// slot already assigned with a previous line (and direction)
							math::Vector3d prev_line_dir = best_slot->line_direction;
							prev_line_dir.normalize();
							double prev_deviation = slot_opp_dir.dot(prev_line_dir);

							if (best_deviation > prev_deviation) {     // the new alignment is better than the previous one
								found_free_slot = false;
								found_compatible_slot = true;
								found_slot = best_slot;
								intersect_sing_point = true;
							}
							else {     // We keep the previous association
								found_compatible_slot = false;
							}
						}
						else {     // WE HAVE A FREE SLOT
							// We keep the previous association
							found_free_slot = true;
							AToSlotIsFree = found_free_slot;
							found_compatible_slot = true;
							found_slot = best_slot;
							intersect_sing_point = true;
						}
						// HAVE WE FIND THE END OF THE LINE??
						if (found_compatible_slot) {
							find_end = true;
						}
					}
					slot_epsilon -= 0.1;
				}
			}

			//==============================================================
			// CASE 2: SO, WE JUST HAVE TO CROSS THE TRIANGLE (GENERAL CONF)
			//==============================================================
			// Does the current triangle has the same classif
			if (!intersect_sing_point) {
				math::Point out_pnt;
				math::Vector3d out_vec;
				TCellID out_cell_id;
				int out_cell_dim;

				m_tool.traverseTriangle(currentFace,          /* the face we work on*/
				                        start_pnt,            /* the geometric point we start from */
				                        start_dir,            /* the geometric direction to follow*/
				                        start_cell_dim,       /* the dimension of the cell start_pnt is located */
				                        start_cell_id,        /* the id of the cell start_pnt is located on*/
				                        out_pnt,              /* the geometric point where we go out */
				                        out_vec,              /* the geometric direction to follow after*/
				                        out_cell_dim,         /* the dimension of the out cell (0 or 1) */
				                        out_cell_id,          /* the id of the out cell*/
				                        streamlineDeviation); /*deviation of the streamline up to this point*/

				accumulatedDistancePerSlotLine = accumulatedDistancePerSlotLine + out_pnt.distance(APoints[APoints.size() - 1]);
				APoints.push_back(out_pnt);

				// we progress to the next point, next vector and so next face too
				prev_dir = start_dir;     // we store the prev direction for slot
				// reconnection with balls
				start_pnt = out_pnt;
				start_dir = out_vec;
				start_cell_dim = out_cell_dim;
				start_cell_id = out_cell_id;
			}
			if (intersect_sing_point) {
				find_end = true;
			}
			// post process, we just check if we haven't arrived onto a geometric boundary
			if (start_cell_dim == 0) {
				Node currentNode = m_mesh->get<Node>(start_cell_id);
				if (m_mesh->isMarked(currentNode, m_mark_nodes_on_point) || m_mesh->isMarked(currentNode, m_mark_nodes_on_curve)) {
					find_end = true;
					AEndOnBnd = true;
				}
			}
			else {     // we have necessarry start_cell_dim=1

				Edge currentEdge = m_mesh->get<Edge>(start_cell_id);
				if (m_mesh->isMarked(currentEdge, m_mark_edges_on_curve)) {
					find_end = true;
					AEndOnBnd = true;
				}
			}
		}
	}

	//==============================================================
	// Update of our parameters
	//==============================================================
	AToDir = start_dir;
	AToSingPnt = found_pnt;
	AToSlot = found_slot;
	AToCellDim = start_cell_dim;
	AToCellID = start_cell_id;
}
/*---------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::growLineRK4(const SingularityPoint::Slot *AFromSlot,
                                       SingularityPoint::Slot *&AToSlot,
                                       const math::Point &AFromPnt,
                                       math::Vector3d &AToDir,
                                       std::vector<math::Point> &APoints,
                                       std::vector<TCellID> &ATriangles,
                                       int &AToCellDim,
                                       TCellID &AToCellID,
                                       double &streamlineDeviation,
                                       const double &stepSize,
                                       bool &AEndOnBnd,
                                       bool &AToSlotIsFree,
                                       bool &must_try_to_connect)
{

	AEndOnBnd = false;
	//========================================================================
	/* grow a line starting from the singularity point AFromSingPnt through the slot AFromSlot
	 as long as the geometric distance travelled by this line
	 is smaller than stepSize or until we reach the boundary (or another slot)
	*/
	//========================================================================
	SingularityPoint *found_pnt = 0;
	SingularityPoint::Slot *found_slot = 0;

	math::Point start_pnt = AFromPnt;     // starting point
	// APoints.push_back(start_pnt);
	math::Vector3d start_dir = AToDir;     // starting direction
	math::Vector3d prev_dir = AToDir;      // prev direction used in the

	TCellID start_cell_id = AToCellID;
	int start_cell_dim = AToCellDim;

	vector<Face> AFaces;

	if (AToCellDim == 0) {
		AFaces = m_mesh->get<Node>(AToCellID).get<Face>();
		if (m_mesh->isMarked(m_mesh->get<Node>(AToCellID), m_mark_nodes_on_curve) || m_mesh->isMarked(m_mesh->get<Node>(AToCellID), m_mark_nodes_on_point)) {

			math::Point point_1 = start_pnt + start_dir * 2 * math::Constants::EPSILON;
			bool insideTri = false;
			for (unsigned int i = 0; i < AFaces.size(); i++) {
				std::vector<Node> f_nodes = AFaces[i].get<Node>();
				math::Triangle ATriangle(f_nodes[0].getPoint(), f_nodes[1].getPoint(), f_nodes[2].getPoint());

				if (ATriangle.isIn2ndMethod(point_1)) {
					insideTri = true;
					break;
				}
			}
			if (!insideTri) {
				AEndOnBnd = true;
				// APoints.push_back(start_pnt);
			}
		}
	}
	else {
		if (AToCellDim == 1) {
			AFaces = m_mesh->get<Edge>(AToCellID).get<Face>();
			if (m_mesh->isMarked(m_mesh->get<Edge>(AToCellID), m_mark_edges_on_curve)) {
				math::Point point_1 = start_pnt + start_dir * 2 * math::Constants::EPSILON;
				bool insideTri = false;
				for (unsigned int i = 0; i < AFaces.size(); i++) {
					std::vector<Node> f_nodes = AFaces[i].get<Node>();
					math::Triangle ATriangle(f_nodes[0].getPoint(), f_nodes[1].getPoint(), f_nodes[2].getPoint());

					if (ATriangle.isIn2ndMethod(point_1)) {
						insideTri = true;
						break;
					}
				}
				if (!insideTri) {
					AEndOnBnd = true;
					// APoints.push_back(start_pnt);
				}
			}
		}
		else
			AFaces.push_back(m_mesh->get<Face>(AToCellID));
	}

	if (!AEndOnBnd) {
		//==============================================================
		// CASE 1: ARE WE IN A FACE CONTAINING A SING. POINT?
		// VERY IMPROBABLE; GOT HERE BUT WE WEREN'T WITHIN THRESHOLD
		//==============================================================

		bool intersect_sing_point = false;
		must_try_to_connect = false;
		SingularityPoint *next_sing_point;

		for (unsigned int i = 0; i < AFaces.size(); i++) {
			next_sing_point = m_faces_to_singularity_on_surf[AFaces[i].id()];
			if (next_sing_point != NULL &&                      // face in a confusing ball ...
			    next_sing_point != AFromSlot->from_point) {     // ... of another singularity point
				must_try_to_connect = true;
				break;
			}
			else {
				if (next_sing_point != NULL &&                      // face in the confusing ball ...
				    next_sing_point == AFromSlot->from_point) {     //... of the incoming singularity point

					// Warning: completly empiric, we just try to detect cyclic lines
					if (APoints.size() >= 100) must_try_to_connect = true;
				}
			}
		}

		if (must_try_to_connect) {
			math::Point start_dirPnt(start_pnt.X() + start_dir.X(), start_pnt.Y() + start_dir.Y(), start_pnt.Z() + start_dir.Z());

			// Now, we look for a compatible slot
			vector<SingularityPoint::Slot *> &cur_slots = next_sing_point->getSlots();

			//==============================================================
			// We look for a free slot and we connect the line to it
			//==============================================================
			bool found_compatible_slot = false;
			bool found_free_slot = false;
			found_pnt = next_sing_point;
			double slot_epsilon = 0.9;
			// we normalize the vector we arrive with
			math::Vector3d current_vec = start_dir;
			current_vec.normalize();
			while (!found_compatible_slot && slot_epsilon > 0.4) {
				double best_deviation = -2;
				double best_slot_id = 0;
				for (unsigned int i_slot = 0; i_slot < cur_slots.size(); i_slot++) {
					// current_vec = start_dir;  						//current_vec.normalize();
					SingularityPoint::Slot *current_slot = cur_slots[i_slot];
					if (current_slot->isFreeze) continue;

					math::Vector3d resultingVec(start_dirPnt, current_slot->location);

					if (start_dir.angle(resultingVec) > math::Constants::PIDIV2) continue;

					math::Vector3d slot_opp_dir = current_slot->direction.opp();
					slot_opp_dir.normalize();
					// current_vec = (current_slot->location) - start_pnt;
					double slot_deviation = slot_opp_dir.dot(current_vec);

					if (slot_deviation > slot_epsilon && slot_deviation > best_deviation) {
						best_deviation = slot_deviation;
						best_slot_id = i_slot;
					}
				}

				if (best_deviation != -2 && !cur_slots.empty()) {
					SingularityPoint::Slot *best_slot = cur_slots[best_slot_id];
					math::Vector3d slot_opp_dir = best_slot->direction.opp();
					slot_opp_dir.normalize();
					if (best_slot->isLaunched) {
						// slot already assigned with a previous line (and direction)
						math::Vector3d prev_line_dir = best_slot->line_direction;
						prev_line_dir.normalize();
						double prev_deviation = slot_opp_dir.dot(prev_line_dir);
						if (best_deviation > prev_deviation) {
							// the new alignment is better than the previous one
							found_free_slot = false;
							found_compatible_slot = true;
							found_slot = best_slot;
							intersect_sing_point = true;
						}
						else {
							// We keep the previous association
							found_compatible_slot = false;
							// must_try_to_connect = false;
						}
					}
					else {     // WE HAVE A FREE SLOT
						// We keep the previous association
						found_free_slot = true;
						AToSlotIsFree = found_free_slot;
						found_compatible_slot = true;
						found_slot = best_slot;
						intersect_sing_point = true;
					}
					/*if(found_compatible_slot) {
					      find_end = true;
					}*/
				}
				slot_epsilon -= 0.1;
			}
		}
		//==============================================================
		// CASE 2: SO, WE JUST HAVE TO CROSS THE TRIANGLE (GENERAL CONF)
		//==============================================================

		if (!intersect_sing_point) {
			math::Point out_pnt;
			math::Vector3d out_vec;

			m_tool.RK4Computation(start_pnt, start_dir, out_pnt, out_vec, streamlineDeviation, stepSize, AEndOnBnd, ATriangles, start_cell_dim, start_cell_id);

			// accumulatedDistancePerSlotLine = accumulatedDistancePerSlotLine + out_pnt.distance(APoints[APoints.size()-1]);
			APoints.push_back(out_pnt);
			// if(!AEndOnBnd)
			//  start_cell_id = ATriangles[ATriangles.size()-1];
			// we progress to the next point, next vector and so next face too
			prev_dir = start_dir;     // we store the prev direction for slot
			start_pnt = out_pnt;
			start_dir = out_vec;
		}
		AToCellDim = start_cell_dim;
		AToCellID = start_cell_id;
	}
	//==============================================================
	// Update of out parameters
	//==============================================================
	// last followed direction
	AToDir = start_dir;
	// singularity point data if we found an end point
	if (found_slot != nullptr) {
		AToSlot->from_point = found_pnt;
		AToSlot = found_slot;
	}
}

/*---------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::writeTestMeshVerts(vector<Node> &ANodes, std::string &AFileName)
{
	Mesh m(MeshModel(DIM3 | F | N | F2N));
	math::Point APoint = ANodes[0].getPoint();
	for (unsigned int j = 1; j < ANodes.size(); j++) {
		math::Point prevPoint = APoint;
		APoint = ANodes[j].getPoint();
		Node n1 = m.newNode(prevPoint.X(), prevPoint.Y(), prevPoint.Z());
		Node n2 = m.newNode(APoint.X(), APoint.Y(), APoint.Z());
		Face f = m.newTriangle(n1, n1, n2);
	}
	IGMeshIOService ioService(&m);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N | F);
	vtkWriter.setDataOptions(N | F);
	vtkWriter.write(AFileName);
}

/*---------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::writeTestPoints(vector<math::Point> &APoints, std::string &AFileName)
{
	Mesh m(MeshModel(DIM3 | F | N | F2N));
	math::Point APoint = APoints[0];
	for (unsigned int j = 1; j < APoints.size(); j++) {
		math::Point prevPoint = APoint;
		APoint = APoints[j];
		Node n1 = m.newNode(prevPoint.X(), prevPoint.Y(), prevPoint.Z());
		Node n2 = m.newNode(APoint.X(), APoint.Y(), APoint.Z());
		Face f = m.newTriangle(n1, n1, n2);
	}
	IGMeshIOService ioService(&m);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N | F);
	vtkWriter.setDataOptions(N | F);
	vtkWriter.write(AFileName);
}
/*---------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::writeTestMeshTriangles(vector<Face> &ATriangles, std::string &AFileName)
{
	Mesh m(MeshModel(DIM3 | F | N | F2N));

	for (unsigned int j = 0; j < ATriangles.size(); j++) {
		vector<Node> currentNodes = ATriangles[j].get<Node>();
		math::Point current_point1 = currentNodes[0].getPoint();
		math::Point current_point2 = currentNodes[1].getPoint();
		math::Point current_point3 = currentNodes[2].getPoint();
		Node n1 = m.newNode(current_point1.X(), current_point1.Y(), current_point1.Z());
		Node n2 = m.newNode(current_point2.X(), current_point2.Y(), current_point2.Z());
		Node n3 = m.newNode(current_point3.X(), current_point3.Y(), current_point3.Z());
		Face f = m.newTriangle(n1, n2, n3);
	}
	IGMeshIOService ioService(&m);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N | F);
	vtkWriter.setDataOptions(N | F);
	vtkWriter.write(AFileName);
}

/*---------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::writeTestMeshTrianglesIds(vector<TCellID> &ATrianglesIds, std::string &AFileName)
{
	Mesh m(MeshModel(DIM3 | F | N | F2N));

	for (unsigned int j = 0; j < ATrianglesIds.size(); j++) {
		vector<Node> currentNodes = (m_mesh->get<Face>(ATrianglesIds[j])).get<Node>();
		math::Point current_point1 = currentNodes[0].getPoint();
		math::Point current_point2 = currentNodes[1].getPoint();
		math::Point current_point3 = currentNodes[2].getPoint();
		Node n1 = m.newNode(current_point1.X(), current_point1.Y(), current_point1.Z());
		Node n2 = m.newNode(current_point2.X(), current_point2.Y(), current_point2.Z());
		Node n3 = m.newNode(current_point3.X(), current_point3.Y(), current_point3.Z());
		Face f = m.newTriangle(n1, n2, n3);
	}
	IGMeshIOService ioService(&m);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N | F);
	vtkWriter.setDataOptions(N | F);
	vtkWriter.write(AFileName);
}

/*---------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::remeshTriangles(Mesh *newLocalMesh,
                                           vector<TCellID> &newLocalMesh_id_to_mesh_id_node,
                                           Variable<math::Cross2D> *local_cross_field_2D,
                                           vector<TCellID> &trianglesToRemesh)
{
	/* remesh the triangles in trianglesToRemesh; assumes that the original mesh is consistently oriented*/

	vector<TCellID> mesh_id_to_newLocalMesh_id_node(m_original_nodes_number);
	vector<bool> visitedVerts(m_original_nodes_number, false);

	math::Point current_point1, current_point2, current_point3, current_point_center;
	unsigned int NodeNumber = 0;
	for (unsigned int i = 0; i < trianglesToRemesh.size(); i++) {
		vector<Node> currentNodes = (m_mesh->get<Face>(trianglesToRemesh[i])).get<Node>();
		current_point1 = currentNodes[0].getPoint();
		current_point2 = currentNodes[1].getPoint();
		current_point3 = currentNodes[2].getPoint();

		Node n1;
		if (!visitedVerts[currentNodes[0].id()]) {
			n1 = newLocalMesh->newNode(current_point1.X(), current_point1.Y(), current_point1.Z());
			visitedVerts[currentNodes[0].id()] = true;
			mesh_id_to_newLocalMesh_id_node[currentNodes[0].id()] = NodeNumber;
			newLocalMesh_id_to_mesh_id_node.push_back(currentNodes[0].id());
			NodeNumber++;
		}
		else {
			n1 = newLocalMesh->get<Node>(mesh_id_to_newLocalMesh_id_node[currentNodes[0].id()]);
		}
		Node n2;
		if (!visitedVerts[currentNodes[1].id()]) {
			n2 = newLocalMesh->newNode(current_point2.X(), current_point2.Y(), current_point2.Z());
			visitedVerts[currentNodes[1].id()] = true;
			mesh_id_to_newLocalMesh_id_node[currentNodes[1].id()] = NodeNumber;
			newLocalMesh_id_to_mesh_id_node.push_back(currentNodes[1].id());
			NodeNumber++;
		}
		else {
			n2 = newLocalMesh->get<Node>(mesh_id_to_newLocalMesh_id_node[currentNodes[1].id()]);
		}
		Node n3;
		if (!visitedVerts[currentNodes[2].id()]) {
			n3 = newLocalMesh->newNode(current_point3.X(), current_point3.Y(), current_point3.Z());
			visitedVerts[currentNodes[2].id()] = true;
			mesh_id_to_newLocalMesh_id_node[currentNodes[2].id()] = NodeNumber;
			newLocalMesh_id_to_mesh_id_node.push_back(currentNodes[2].id());
			NodeNumber++;
		}
		else {
			n3 = newLocalMesh->get<Node>(mesh_id_to_newLocalMesh_id_node[currentNodes[2].id()]);
		}

		Node newn = newLocalMesh->newNode(m_triangle_centers[trianglesToRemesh[i]].X(), m_triangle_centers[trianglesToRemesh[i]].Y(),
		                                  m_triangle_centers[trianglesToRemesh[i]].Z());
		newLocalMesh_id_to_mesh_id_node.push_back(m_original_nodes_number + 1);
		NodeNumber++;
		Face f1 = newLocalMesh->newTriangle(n1, n2, newn);
		Face f2 = newLocalMesh->newTriangle(n2, n3, newn);
		Face f3 = newLocalMesh->newTriangle(n3, n1, newn);
	}

	IGMeshIOService ioService(newLocalMesh);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N | F);
	vtkWriter.setDataOptions(N | F);
	vtkWriter.write("localmesh.vtk");

	MeshDoctor doc(newLocalMesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	vector<Node> m_curve_nodes, m_surf_nodes;
	visitedVerts.clear();
	visitedVerts.resize(newLocalMesh->getNbNodes(), false);

	math::Vector3d OX(1, 0, 0);
	for (auto e_id : newLocalMesh->edges()) {
		Edge currentEdge = newLocalMesh->get<Edge>(e_id);
		vector<Node> currentNodes = (newLocalMesh->get<Edge>(e_id)).get<Node>();
		std::vector<Face> adj_faces = currentEdge.get<Face>();
		if (adj_faces.size() == 1) {
			if (!visitedVerts[currentNodes[0].id()]) {
				m_curve_nodes.push_back(currentNodes[0]);
				(*local_cross_field_2D)[currentNodes[0].id()] = (*m_field)[newLocalMesh_id_to_mesh_id_node[currentNodes[0].id()]];
				visitedVerts[currentNodes[0].id()] = true;
			}

			if (!visitedVerts[currentNodes[1].id()]) {
				m_curve_nodes.push_back(currentNodes[1]);
				(*local_cross_field_2D)[currentNodes[1].id()] = (*m_field)[newLocalMesh_id_to_mesh_id_node[currentNodes[1].id()]];
				visitedVerts[currentNodes[0].id()] = true;
			}
		}
	}

	for (auto n_id : newLocalMesh->nodes()) {
		if (!visitedVerts[n_id]) m_surf_nodes.push_back(newLocalMesh->get<Node>(n_id));
	}

	/*visualize crossVect for newLocalMesh
	Variable<math::Vector3d>* crossVect0 = newLocalMesh->newVariable<math::Vector3d,GMDS_NODE>("crossVect0");
	Variable<math::Vector3d>* crossVect1 = newLocalMesh->newVariable<math::Vector3d,GMDS_NODE>("crossVect1");
	Variable<math::Vector3d>* crossVect2 = newLocalMesh->newVariable<math::Vector3d,GMDS_NODE>("crossVect2");
	Variable<math::Vector3d>* crossVect3 = newLocalMesh->newVariable<math::Vector3d,GMDS_NODE>("crossVect3");
	for(auto f_id:newLocalMesh->faces()){
	   Face current = newLocalMesh->get<Face>(f_id);
	   std::vector<TCellID> nodeIDs = current.getIDs<Node>();
	   int ID1 = nodeIDs[0];
	   int ID2 = nodeIDs[1];
	   int ID3 = nodeIDs[2];

	   std::vector<math::Vector3d> crossV_ID1 = (*local_cross_field_2D)[ID1].componentVectors();
	   std::vector<math::Vector3d> crossV_ID2 = (*local_cross_field_2D)[ID2].componentVectors();
	   std::vector<math::Vector3d> crossV_ID3 = (*local_cross_field_2D)[ID3].componentVectors();
	   (*crossVect0)[ID1] = crossV_ID1[0];
	   (*crossVect0)[ID2] = crossV_ID2[0];
	   (*crossVect0)[ID3] = crossV_ID3[0];
	   (*crossVect1)[ID1] = crossV_ID1[1];
	   (*crossVect1)[ID2] = crossV_ID2[1];
	   (*crossVect1)[ID3] = crossV_ID3[1];
	   (*crossVect2)[ID1] = crossV_ID1[2];
	   (*crossVect2)[ID2] = crossV_ID2[2];
	   (*crossVect2)[ID3] = crossV_ID3[2];
	   (*crossVect3)[ID1] = crossV_ID1[3];
	   (*crossVect3)[ID2] = crossV_ID2[3];
	   (*crossVect3)[ID3] = crossV_ID3[3];
	}


	IGMeshIOService meshIoServref(newLocalMesh);
	VTKWriter writerref(&meshIoServref);
	writerref.setCellOptions(N|F);
	writerref.setDataOptions(N|F);

	std::stringstream file_nameref;
	file_nameref <<"before-newLocalMesh-crossVectors.vtk";
	   writerref.write(file_nameref.str());
	   */
	std::vector<Face> mesh_faces;
	mesh_faces.resize(newLocalMesh->getNbFaces());
	int f_index = 0;
	for (auto f_id : newLocalMesh->faces()) {
		mesh_faces[f_index++] = newLocalMesh->get<Face>(f_id);
	}

	LaplaceCross2D algo(newLocalMesh, local_cross_field_2D, m_curve_nodes, m_surf_nodes, mesh_faces);
	algo.execute();

	/*
	Variable<math::Vector3d>* crossVect0 = newLocalMesh->newVariable<math::Vector3d,GMDS_NODE>("crossVect0");
	Variable<math::Vector3d>* crossVect1 = newLocalMesh->newVariable<math::Vector3d,GMDS_NODE>("crossVect1");
	Variable<math::Vector3d>* crossVect2 = newLocalMesh->newVariable<math::Vector3d,GMDS_NODE>("crossVect2");
	Variable<math::Vector3d>* crossVect3 = newLocalMesh->newVariable<math::Vector3d,GMDS_NODE>("crossVect3");
	for(auto f_id:newLocalMesh->faces()){
	   Face current = newLocalMesh->get<Face>(f_id);
	   std::vector<TCellID> nodeIDs = current.getIDs<Node>();
	   int ID1 = nodeIDs[0];
	   int ID2 = nodeIDs[1];
	   int ID3 = nodeIDs[2];


	   std::vector<math::Vector3d> crossV_ID1 = (*local_cross_field_2D)[ID1].componentVectors();
	   std::vector<math::Vector3d> crossV_ID2 = (*local_cross_field_2D)[ID2].componentVectors();
	   std::vector<math::Vector3d> crossV_ID3 = (*local_cross_field_2D)[ID3].componentVectors();
	   (*crossVect0)[ID1] = crossV_ID1[0];
	   (*crossVect0)[ID2] = crossV_ID2[0];
	   (*crossVect0)[ID3] = crossV_ID3[0];
	   (*crossVect1)[ID1] = crossV_ID1[1];
	   (*crossVect1)[ID2] = crossV_ID2[1];
	   (*crossVect1)[ID3] = crossV_ID3[1];
	   (*crossVect2)[ID1] = crossV_ID1[2];
	   (*crossVect2)[ID2] = crossV_ID2[2];
	   (*crossVect2)[ID3] = crossV_ID3[2];
	   (*crossVect3)[ID1] = crossV_ID1[3];
	   (*crossVect3)[ID2] = crossV_ID2[3];
	   (*crossVect3)[ID3] = crossV_ID3[3];
	}


	IGMeshIOService meshIoServref(newLocalMesh);
	VTKWriter writerref(&meshIoServref);
	writerref.setCellOptions(N|F);
	writerref.setDataOptions(N|F);

	std::stringstream file_nameref;
	file_nameref <<"final-newLocalMesh-crossVectors.vtk";
	   writerref.write(file_nameref.str()); */
}

/*---------------------------------------------------------------------------*/

void
SingularityGraphBuilder2D::remeshTrianglesNewMesh(Mesh *newLocalMesh,
                                                  Variable<math::Cross2D> *newMesh_cross_field_2D,
                                                  vector<TCellID> &trianglesToRemesh,
                                                  vector<bool> &trianglesToRemeshBool,
                                                  vector<math::Point> &newTriangleCenters,
                                                  vector<math::Cross2D> &newTriangleCenterCrosses,
                                                  vector<math::Vector3d> &newBdryEdgeNormals,
                                                  vector<math::Vector3d> &newBdryNodeNormals,
                                                  vector<bool> &isCurveEdge,
                                                  vector<bool> &isCurveNode,
                                                  vector<TCellID> &newLocalMesh_id_to_mesh_id_node,
                                                  vector<TCellID> &mesh_id_to_newLocalMesh_id_node,
                                                  vector<TCellID> &newLocalMesh_id_to_mesh_id_face,
                                                  vector<TCellID> &mesh_id_to_newLocalMesh_id_face)
{
	/* function assumes that the original mesh is consistently oriented*/
	/* Method I. Introduce one single vertex in the center of each triangle to remesh ->
	 * will result 3 new triangles for each original one
	 * 	*easiest, least computationally expensive, but not as good
	 * Method II. Introduce 3 vertices for each triangle to remesh (in the center of each edge) ->
	 * result 4 new triangles for each original one
	 * - adjacent triangles remeshed as:
	 * a) introduce edge as [center edge (for remeshed edge), opossite original node]
	 * b) remesh each adjacent traingle as the trianglesToRemesh triangles
	 * => the entire mesh has to be remeshed this way...
	 * in this case, ideally the entire algorithm should be recomputed (starting with sing detection)
	 * simply put, would be easier to just remsh in GMSH...this is the method used
	 * MIIa) implemented here
	 */
	cout << "remeshTrianglesNewMesh" << endl;

	vector<bool> visitedVerts(m_original_nodes_number, false);

	math::Point current_point1, current_point2, current_point3, current_point_center;
	unsigned int NodeNumber = 0;
	unsigned int FaceNumber = 0;
	unsigned int EdgeNumber = 0;
	vector<TCellID> origEdge2NewPoint(m_mesh->getNbEdges());
	vector<bool> modifiedEdge(m_mesh->getNbEdges(), false);
	vector<Node> NodesToAdd(6);
	vector<TCoord> edgeWeights(2, 0.5);
	for (unsigned int i = 0; i < trianglesToRemesh.size(); i++) {
		Face currentFace = m_mesh->get<Face>(trianglesToRemesh[i]);
		vector<Node> currentNodes = currentFace.get<Node>();
		current_point1 = currentNodes[0].getPoint();
		current_point2 = currentNodes[1].getPoint();
		current_point3 = currentNodes[2].getPoint();

		if (!visitedVerts[currentNodes[0].id()]) {
			NodesToAdd[0] = newLocalMesh->newNode(current_point1.X(), current_point1.Y(), current_point1.Z());
			visitedVerts[currentNodes[0].id()] = true;
			mesh_id_to_newLocalMesh_id_node[currentNodes[0].id()] = NodeNumber;
			newLocalMesh_id_to_mesh_id_node.push_back(currentNodes[0].id());
			(*newMesh_cross_field_2D)[NodeNumber] = (*m_field)[currentNodes[0].id()];
			NodeNumber++;
		}
		else {
			NodesToAdd[0] = newLocalMesh->get<Node>(mesh_id_to_newLocalMesh_id_node[currentNodes[0].id()]);
		}

		if (!visitedVerts[currentNodes[1].id()]) {
			NodesToAdd[1] = newLocalMesh->newNode(current_point2.X(), current_point2.Y(), current_point2.Z());
			visitedVerts[currentNodes[1].id()] = true;
			mesh_id_to_newLocalMesh_id_node[currentNodes[1].id()] = NodeNumber;
			newLocalMesh_id_to_mesh_id_node.push_back(currentNodes[1].id());
			(*newMesh_cross_field_2D)[NodeNumber] = (*m_field)[currentNodes[1].id()];
			NodeNumber++;
		}
		else {
			NodesToAdd[1] = newLocalMesh->get<Node>(mesh_id_to_newLocalMesh_id_node[currentNodes[1].id()]);
		}

		if (!visitedVerts[currentNodes[2].id()]) {
			NodesToAdd[2] = newLocalMesh->newNode(current_point3.X(), current_point3.Y(), current_point3.Z());
			visitedVerts[currentNodes[2].id()] = true;
			mesh_id_to_newLocalMesh_id_node[currentNodes[2].id()] = NodeNumber;
			newLocalMesh_id_to_mesh_id_node.push_back(currentNodes[2].id());
			(*newMesh_cross_field_2D)[NodeNumber] = (*m_field)[currentNodes[2].id()];
			NodeNumber++;
		}
		else {
			NodesToAdd[2] = newLocalMesh->get<Node>(mesh_id_to_newLocalMesh_id_node[currentNodes[2].id()]);
		}

		vector<Edge> currentEdges = currentFace.get<Edge>();
		for (unsigned int k = 0; k < 3; k++) {
			vector<Node> edge_nodes = currentEdges[k].get<Node>();
			for (unsigned int j = 0; j < 3; j++) {     // get opposite node for each edge; put in this order
				if ((edge_nodes[0] != NodesToAdd[j]) && (edge_nodes[1] != NodesToAdd[j])) {
					if (!modifiedEdge[currentEdges[k].id()]) {
						math::Point middle_point = currentEdges[k].center();
						NodesToAdd[3 + fmod(j, 3)] = newLocalMesh->newNode(middle_point.X(), middle_point.Y(), middle_point.Z());

						vector<math::Vector3d> c_vectors0 = ((*m_field)[edge_nodes[0].id()]).componentVectors();
						vector<math::Vector3d> c_vectors1 = ((*m_field)[edge_nodes[1].id()]).componentVectors();

						math::Vector3d second_closest_vect = (*m_field)[edge_nodes[1].id()].closestComponentVector(c_vectors0[0]);
						math::Vector3d center_vect = c_vectors0[0] * edgeWeights[0] + second_closest_vect * edgeWeights[1];
						math::Vector3d center_vect_second = center_vect.getOneOrtho();

						(*newMesh_cross_field_2D)[NodeNumber] = math::Cross2D(center_vect, center_vect_second);
						origEdge2NewPoint[currentEdges[k].id()] = EdgeNumber;
						modifiedEdge[currentEdges[k].id()] = true;
						EdgeNumber++;
						NodeNumber++;
					}
					else {
						NodesToAdd[3 + fmod(j, 3)] = newLocalMesh->get<Node>(origEdge2NewPoint[currentEdges[k].id()]);
					}
				}
			}
		}
		Face f1 = newLocalMesh->newTriangle(NodesToAdd[3], NodesToAdd[4], NodesToAdd[5]);
		Face f2 = newLocalMesh->newTriangle(NodesToAdd[0], NodesToAdd[5], NodesToAdd[4]);
		Face f3 = newLocalMesh->newTriangle(NodesToAdd[1], NodesToAdd[3], NodesToAdd[5]);
		Face f4 = newLocalMesh->newTriangle(NodesToAdd[2], NodesToAdd[4], NodesToAdd[3]);

		newLocalMesh_id_to_mesh_id_face.push_back(currentFace.id());
		newLocalMesh_id_to_mesh_id_face.push_back(currentFace.id());
		newLocalMesh_id_to_mesh_id_face.push_back(currentFace.id());
		newLocalMesh_id_to_mesh_id_face.push_back(currentFace.id());
		FaceNumber = FaceNumber + 4;
		mesh_id_to_newLocalMesh_id_face[currentFace.id()] = FaceNumber;

		newTriangleCenters.push_back(f1.center());
		newTriangleCenters.push_back(f2.center());
		newTriangleCenters.push_back(f3.center());
		newTriangleCenters.push_back(f4.center());
	}

	// now treat neighbouring faces
	vector<bool> remeshedAdjFaces(m_original_faces_number, false);

	for (auto e_id : m_mesh->edges()) {
		Edge currentEdge = m_mesh->get<Edge>(e_id);

		if (modifiedEdge[e_id]) {
			vector<Face> adj_faces = currentEdge.get<Face>();
			for (unsigned int i = 0; i < adj_faces.size(); i++) {
				if (!trianglesToRemeshBool[adj_faces[i].id()]) {
					vector<Node> edge_nodes = currentEdge.get<Node>();
					Node midNode = newLocalMesh->get<Node>(origEdge2NewPoint[currentEdge.id()]);
					if (remeshedAdjFaces[adj_faces[i].id()]) {
						// TODO
					}
					else {
						vector<Node> face_nodes = adj_faces[i].get<Node>();

						for (unsigned int j = 0; j < 3; j++) {
							if ((face_nodes[j] != edge_nodes[0]) && (face_nodes[j] != edge_nodes[1])) {
								Node opposite_node = face_nodes[j];
								Node opposite_node_local;     // TODO = does it exist? or to add?
								// same for face_nodes!!!TODO
								Face f1 = newLocalMesh->newTriangle(opposite_node_local, face_nodes[fmod((j + 1), 3)], midNode);
								Face f2 = newLocalMesh->newTriangle(opposite_node_local, midNode, face_nodes[fmod((j + 2), 3)]);
								break;
							}
						}

						remeshedAdjFaces[adj_faces[i].id()] = true;
					}
				}
			}
		}
	}
	IGMeshIOService ioService(newLocalMesh);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N | F);
	vtkWriter.setDataOptions(N | F);
	vtkWriter.write("localmesh.vtk");

	MeshDoctor doc(newLocalMesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	isCurveEdge.resize(newLocalMesh->getNbEdges(), false);
	isCurveNode.resize(newLocalMesh->getNbNodes(), false);
	vector<Node> m_curve_nodes, m_surf_nodes;
	visitedVerts.clear();
	visitedVerts.resize(newLocalMesh->getNbNodes(), false);

	math::Vector3d OX(1, 0, 0);
	for (auto e_id : newLocalMesh->edges()) {
		Edge currentEdge = newLocalMesh->get<Edge>(e_id);
		vector<Node> currentNodes = (newLocalMesh->get<Edge>(e_id)).get<Node>();
		std::vector<Face> adj_faces = currentEdge.get<Face>();
		if (adj_faces.size() == 1) {
			isCurveEdge[e_id] = true;
			if (!visitedVerts[currentNodes[0].id()]) {
				isCurveNode[currentNodes[0].id()] = true;
				m_curve_nodes.push_back(currentNodes[0]);
				(*newMesh_cross_field_2D)[currentNodes[0].id()] = (*m_field)[newLocalMesh_id_to_mesh_id_node[currentNodes[0].id()]];
				visitedVerts[currentNodes[0].id()] = true;
			}

			if (!visitedVerts[currentNodes[1].id()]) {
				isCurveNode[currentNodes[1].id()] = true;
				m_curve_nodes.push_back(currentNodes[1]);
				(*newMesh_cross_field_2D)[currentNodes[1].id()] = (*m_field)[newLocalMesh_id_to_mesh_id_node[currentNodes[1].id()]];
				visitedVerts[currentNodes[0].id()] = true;
			}
		}
	}

	// unfortunately we recompute the frame algorithm for the entire new mesh

	newTriangleCenterCrosses.resize(newLocalMesh->getNbFaces());
	vector<TCoord> AWeights(3, (double) 1 / 3);
	for (auto f_id : newLocalMesh->faces()) {
		std::vector<math::Cross2D> Tricrosses;
		Face currentFace = newLocalMesh->get<Face>(f_id);
		vector<Node> nodesTri = currentFace.get<Node>();
		Tricrosses.push_back((*newMesh_cross_field_2D)[nodesTri[0].id()]);
		Tricrosses.push_back((*newMesh_cross_field_2D)[nodesTri[1].id()]);
		Tricrosses.push_back((*newMesh_cross_field_2D)[nodesTri[2].id()]);
		vector<math::Vector3d> c_vectors0 = Tricrosses[0].componentVectors();
		vector<math::Vector3d> c_vectors1 = Tricrosses[1].componentVectors();
		vector<math::Vector3d> c_vectors2 = Tricrosses[2].componentVectors();

		math::Vector3d second_closest_vect = Tricrosses[1].closestComponentVector(c_vectors0[0]);
		math::Vector3d third_closest_vect = Tricrosses[2].closestComponentVector(c_vectors0[0]);
		math::Vector3d center_vect = c_vectors0[0] * AWeights[0] + second_closest_vect * AWeights[1] + third_closest_vect * AWeights[2];
		math::Vector3d center_vect_second = center_vect.getOneOrtho();

		newTriangleCenterCrosses[f_id] = math::Cross2D(center_vect, center_vect_second);
	}

	newBdryEdgeNormals.resize(newLocalMesh->getNbEdges());
	math::Vector3d zeroVec(0.0, 0.0, 0.0);
	newBdryNodeNormals.resize(newLocalMesh->getNbNodes(), zeroVec);

	for (auto e_id : newLocalMesh->edges()) {
		Edge currentEdge = newLocalMesh->get<Edge>(e_id);
		if (isCurveEdge[currentEdge.id()]) {
			vector<Node> adjacent_nodes = currentEdge.get<Node>();
			math::Point p1 = adjacent_nodes[0].getPoint();
			math::Point p2 = adjacent_nodes[1].getPoint();
			math::Vector3d v1 = math::Vector3d(p1, p2);
			v1.normalize();
			vector<Face> adjacent_faces = currentEdge.get<Face>();
			newBdryEdgeNormals[e_id] = v1.cross(adjacent_faces[0].normal());
		}
	}

	for (auto n_id : newLocalMesh->nodes()) {
		Node currentNode = newLocalMesh->get<Node>(n_id);
		if (isCurveNode[currentNode.id()]) {
			vector<Edge> currentEdges = currentNode.get<Edge>();
			for (unsigned int i = 0; i < currentEdges.size(); i++) {
				if (isCurveEdge[currentEdges[i].id()]) {
					newBdryNodeNormals[n_id][0] = newBdryNodeNormals[n_id][0] + newBdryEdgeNormals[currentEdges[i].id()][0];
					newBdryNodeNormals[n_id][1] = newBdryNodeNormals[n_id][1] + newBdryEdgeNormals[currentEdges[i].id()][1];
					newBdryNodeNormals[n_id][2] = newBdryNodeNormals[n_id][2] + newBdryEdgeNormals[currentEdges[i].id()][2];
				}
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder2D::connectSingularityLines(math::Vector3d dir_slot_i,
                                                   math::Vector3d dir_slot_j,
                                                   math::Vector3d connection_dir,
                                                   math::Point start_pnt1,     // start_pnt = line_discretization[0]; //starting point
                                                   math::Point start_pnt2,
                                                   TCellID start_cell_id1,
                                                   TCellID start_cell_id2,
                                                   int start_cell_dim1,
                                                   int start_cell_dim2,
                                                   std::vector<math::Point> &APoints,
                                                   std::vector<TCellID> &ATriangles,
                                                   double &streamlineDeviation)
{
	cout << "connectSingularityLines " << endl;

	/* "transport" the last closest comp vect (to dir_slot_i) until singline2;
	 * if it corresponds with (opposite) dir_slot_j, we can connect.
	 * the measure could be the angle diff along the connection_dir */

	ATriangles.clear();
	APoints.clear();

	math::Point current_pnt = start_pnt1;
	math::Vector3d current_vec = connection_dir;

	math::Point start_pnt = start_pnt1;
	math::Vector3d start_dir = connection_dir;

	TCellID start_cell_id = start_cell_id1;
	int start_cell_dim = start_cell_dim1;

	math::Point start_dirPnt(start_pnt1.X() + connection_dir.X(), start_pnt1.Y() + connection_dir.Y(), start_pnt1.Z() + connection_dir.Z());

	std::vector<math::Vector3d> compVectors;
	math::Cross2D cross_cell;
	if (start_cell_dim1 == 0) {     // on node
		compVectors = (*m_field)[start_cell_id1].componentVectors();
		cross_cell = (*m_field)[start_cell_id1];
	}
	else {     // on edge
		vector<Node> currentNodes = (m_mesh->get<Edge>(start_cell_id1)).get<Node>();
		math::Vector3d AB(currentNodes[0].getPoint(), currentNodes[1].getPoint());
		math::Vector3d AC(currentNodes[0].getPoint(), start_pnt1);
		double alpha = AC.norm() / AB.norm();

		std::vector<math::Vector3d> compVectorsA, compVectorsB;
		math::Cross2D cross_cellA, cross_cellB;
		compVectorsA = (*m_field)[currentNodes[0].id()].componentVectors();
		compVectorsB = (*m_field)[currentNodes[1].id()].componentVectors();
		cross_cellA = (*m_field)[currentNodes[0].id()];
		cross_cellB = (*m_field)[currentNodes[1].id()];
		compVectors.resize(4);
		for (unsigned int i = 0; i < 4; i++) {
			compVectors[i] = alpha * compVectorsA[i] + (1.0 - alpha) * compVectorsB[i];
			compVectors[i].normalize();
		}
		cross_cell = math::Cross2D::mean(cross_cellA, alpha, cross_cellB, (1.0 - alpha));
	}

	math::Vector3d closestCompVect = compVectors[0];
	TCoord val = dir_slot_i.dot(closestCompVect);
	for (unsigned int i = 1; i < 4; i++) {
		math::Vector3d v_i = compVectors[i];
		TCoord val_i = connection_dir.dot(v_i);
		if (val_i > val) {
			val = val_i;
			closestCompVect = v_i;
		}
	}

	//========================================================================
	// Main loop to create the singularity line
	//========================================================================
	TCellID next_cell_id = NullID;
	while (start_cell_id2 != next_cell_id) {
		next_cell_id = NullID;
		int next_cell_dim = -1;
		m_tool.findNextCell(start_pnt, start_dir, start_cell_dim, start_cell_id, next_cell_dim, next_cell_id);

		if (next_cell_dim == -1) {
			throw GMDSException("connection tries to pass through concavity");
		}
		else {
			if (next_cell_dim == 1) {
				/* we are going along an edge. Our simple assumption is to follow this edge until
				 * reaching one of its end points and to compute the next direction at this point.*/
				Edge currentEdge = m_mesh->get<Edge>(next_cell_id);
				std::vector<TCellID> adj_faces = currentEdge.getIDs<Face>();
				ATriangles.insert(ATriangles.end(), adj_faces.begin(), adj_faces.end());

				std::vector<Node> currentNodes = currentEdge.get<Node>();
				math::Vector3d v0(start_pnt, currentNodes[0].getPoint());
				math::Vector3d v1(start_pnt, currentNodes[1].getPoint());
				Node next_node;
				if (math::near(v0.norm(), 0.0))
					next_node = currentNodes[1];
				else if (math::near(v1.norm(), 0.0))
					next_node = currentNodes[0];
				else if (v0.dot(start_dir) > v1.dot(start_dir))
					next_node = currentNodes[0];
				else
					next_node = currentNodes[1];

				math::Vector3d next_dir;

				start_pnt = next_node.getPoint();
				start_cell_dim = 0;
				start_cell_id = next_node.id();

				APoints.push_back(start_pnt);

				compVectors = (*m_field)[start_cell_id].componentVectors();
				cross_cell = (*m_field)[start_cell_id];
				math::Vector3d closestCompVect = compVectors[0];
				TCoord val = dir_slot_i.dot(closestCompVect);
				for (unsigned int i = 1; i < 4; i++) {
					math::Vector3d v_i = compVectors[i];
					TCoord val_i = connection_dir.dot(v_i);
					if (val_i > val) {
						val = val_i;
						closestCompVect = v_i;
					}
				}
			}
			else {     // general case, we are in a face
				Face currentFace = m_mesh->get<Face>(next_cell_id);
				ATriangles.push_back(currentFace.id());

				//==============================================================
				// CASE 2: SO, WE JUST HAVE TO CROSS THE TRIANGLE (GENERAL CONF)
				//==============================================================

				math::Point out_pnt, min_out_pnt;
				TCellID out_cell_id, min_out_cell_id;
				int out_cell_dim, min_out_cell_dim;

				vector<Node> currentNodes = currentFace.get<Node>();
				vector<Edge> currentEdges = currentFace.get<Edge>();

				double min_value = 7.7;
				bool found = false;

				for (unsigned int t = 0; t < 3; t++) {
					math::Point node_loc = currentNodes[t].getPoint();
					math::Vector3d trial_vector(start_pnt, node_loc);
					trial_vector.normalize();
					double temp_dot = trial_vector.dot(connection_dir);
					if (temp_dot > min_value) {
						min_value = temp_dot;
						min_out_cell_dim = 0;
						min_out_cell_id = currentNodes[t].id();
						min_out_pnt = node_loc;
					}

					if (math::near(temp_dot - 1, 0)) {
						out_cell_dim = 0;
						out_cell_id = currentNodes[t].id();
						out_pnt = node_loc;
						found = true;
						break;
					}

					math::Vector3d out_vec;
					double deviation = 0.0;

					bool intersectEdge = m_tool.computeOutVectorFromRayAndEdge(currentEdges[t], start_pnt, connection_dir, out_pnt, out_vec, deviation);
					if (intersectEdge) {
						out_cell_dim = 1;
						out_cell_id = currentEdges[t].id();
					}
				}

				if (!found) {
					out_cell_dim = min_out_cell_dim;
					out_cell_id = min_out_cell_id;
					out_pnt = min_out_pnt;
				}

				APoints.push_back(out_pnt);

				start_pnt = out_pnt;

				start_cell_dim = out_cell_dim;
				start_cell_id = out_cell_id;

				if (start_cell_dim == 0) {     // on node
					compVectors = (*m_field)[start_cell_id].componentVectors();
					cross_cell = (*m_field)[start_cell_id];
				}
				else {     // on edge
					vector<Node> currentNodes = (m_mesh->get<Edge>(start_cell_id)).get<Node>();
					math::Vector3d AB(currentNodes[0].getPoint(), currentNodes[1].getPoint());
					math::Vector3d AC(currentNodes[0].getPoint(), start_pnt);
					double alpha = AC.norm() / AB.norm();

					std::vector<math::Vector3d> compVectorsA, compVectorsB;
					math::Cross2D cross_cellA, cross_cellB;
					compVectorsA = (*m_field)[currentNodes[0].id()].componentVectors();
					compVectorsB = (*m_field)[currentNodes[1].id()].componentVectors();
					cross_cellA = (*m_field)[currentNodes[0].id()];
					cross_cellB = (*m_field)[currentNodes[1].id()];
					compVectors.resize(4);
					for (unsigned int i = 0; i < 4; i++) {
						compVectors[i] = alpha * compVectorsA[i] + (1.0 - alpha) * compVectorsB[i];
						compVectors[i].normalize();
					}
					cross_cell = math::Cross2D::mean(cross_cellA, alpha, cross_cellB, (1.0 - alpha));
				}
				math::Vector3d closestCompVect = compVectors[0];
				TCoord val = dir_slot_i.dot(closestCompVect);
				for (unsigned int i = 1; i < 4; i++) {
					math::Vector3d v_i = compVectors[i];
					TCoord val_i = connection_dir.dot(v_i);
					if (val_i > val) {
						val = val_i;
						closestCompVect = v_i;
					}
				}
			}
		}
	}
}

void
SingularityGraphBuilder2D::computeFace2FaceInfo()
{
	std::cout << "computeFace2FaceInfo" << std::endl;
	/*
	* WARNING if 3D, FIELD SHOULD BE DEFINED PER FACE! ; WE ALSO MUST TAKE INTO ACCOUNT
	* THE MISSMATCH BETWEEN 2 TRIANGLES (ACROSS THE EDGE) -
	// here we assume the mesh is manifold
	*/

	math::Vector3d zeroVec(0.0, 0.0, 0.0);
	if (m_face_normals.size() == 0) m_face_normals.resize(m_original_faces_number);

	// For every non-border edge
	for (auto e_id : m_mesh->edges()) {
		Edge currentEdge = m_mesh->get<Edge>(e_id);

		if (!m_mesh->isMarked(currentEdge, m_mark_edges_on_curve)) {

			vector<Face> adjacent_triangles = currentEdge.get<Face>();
			unsigned int fid0 = adjacent_triangles[0].id();
			unsigned int fid1 = adjacent_triangles[1].id();

			math::Vector3d N0 = m_face_normals[fid0];

			// find common edge on triangle 0 and 1
			int fid0_vc = -1;     // local id of common edge in f0 {0, 1 or 2}
			int fid1_vc = -1;     // local id of common edge in f1
			vector<Edge> currentEdges0 = adjacent_triangles[0].get<Edge>();
			vector<Edge> currentEdges1 = adjacent_triangles[1].get<Edge>();
			for (unsigned int i = 0; i < 3; i++) {
				if (currentEdges0[i].id() == e_id) fid0_vc = i;
				if (currentEdges1[i].id() == e_id) fid1_vc = i;
			}

			vector<Node> edge_nodes0 = currentEdges0[fid0_vc].get<Node>();
			math::Vector3d common_edge(edge_nodes0[0].getPoint(), edge_nodes0[1].getPoint());
			common_edge.normalize();

			// Common local Basis (CLB) Map the two triangles in a new space where the common edge is the x axis and the N0
			// the z axis

			Eigen::MatrixXd CLB(3, 3);
			math::Point commonOrig = edge_nodes0[0].getPoint();     // origin ; the first vert of common edge
			math::Vector3d tmp = common_edge.cross(m_face_normals[fid0]);

			for (unsigned int i = 0; i < 3; i++) {
				CLB(0, i) = common_edge[i];
				CLB(1, i) = tmp[i];
				CLB(2, i) = m_face_normals[fid0][i];
			}
			vector<Node> face_nodes0 = adjacent_triangles[0].get<Node>();
			Eigen::MatrixXd V0(3, 3);
			for (unsigned int i = 0; i < 3; i++) {
				V0(i, 0) = face_nodes0[i].X() - commonOrig.X();
				V0(i, 1) = face_nodes0[i].Y() - commonOrig.Y();
				V0(i, 2) = face_nodes0[i].Z() - commonOrig.Z();
			}

			V0 = (CLB * V0.transpose()).transpose();

			vector<Node> face_nodes1 = adjacent_triangles[1].get<Node>();
			Eigen::MatrixXd V1(3, 3);

			for (unsigned int i = 0; i < 3; i++) {
				V1(i, 0) = face_nodes1[i].X() - commonOrig.X();
				V1(i, 1) = face_nodes1[i].Y() - commonOrig.Y();
				V1(i, 2) = face_nodes1[i].Z() - commonOrig.Z();
			}
			V1 = (CLB * V1.transpose()).transpose();

			// compute rotation R such that R * N1 = N0
			// i.e. map both triangles to the same plane
			double alpha = -atan2(V1((fid1_vc + 2) % 3, 2), V1((fid1_vc + 2) % 3, 1));

			Eigen::MatrixXd R(3, 3);
			// Eigen::Matrix<typename DerivedV::Scalar, 3, 3> R;
			R << 1, 0, 0, 0, cos(alpha), -sin(alpha), 0, sin(alpha), cos(alpha);
			V1 = (R * V1.transpose()).transpose();

			// measure the angle between the reference frames
			// k_ij is the angle between the triangle on the left and the one on the right
			Eigen::MatrixXd ref0(1, 3);
			Eigen::MatrixXd ref1(1, 3);
			for (unsigned int i = 0; i < 3; i++) {
				ref1(0, i) = V1(i, 1) - V1(i, 0);
				ref0(0, i) = V0(0, 1) - V0(0, 0);
			}

			ref0.normalize();
			ref1.normalize();

			double ktemp = atan2(ref1(1), ref1(0)) - atan2(ref0(1), ref0(0));

			// just to be sure, rotate ref0 using angle ktemp...
			Eigen::MatrixXd R2(2, 2);
			R2 << cos(ktemp), -sin(ktemp), sin(ktemp), cos(ktemp);
		}
	}
}

/*--------------------------------------------------------------------------------------------------------------------*/

LineRelaxator::LineRelaxator(SingularityGraph *graph, const double meanEdgeLength)
{
	m_graph = graph;
	m_step = meanEdgeLength / 2;
	m_graph->computeLinkedEdges();
	computeFixedGroupAndFixedLine();
	for (auto line : m_graph->getLines()) {
		const auto slots = line->getSlots();
		const auto direction = math::Vector3d(slots[0]->from_point->getLocation(), slots[1]->from_point->getLocation());
		slots[0]->line_direction = direction;
		slots[1]->line_direction = direction.opp();
	}
	for (auto singPoint : m_graph->getPoints()) {
		if (singPoint->getSlots().size() > 3) {
			singPoint->sortSlots();
		}
		m_angleAtSingularity[singPoint->getNumber()] = computeAngleAtSingularity(singPoint);
	}
}

void
LineRelaxator::computeFixedGroupAndFixedLine()
{
	m_groupIsFixed = std::vector<bool>(m_graph->getNLinkedEdgeGroup(), false);
	m_lineIsfixed.clear();
	for (const auto line : m_graph->getLines()) {
		const auto points = line->getEndPoints();

		if (points[0]->getGeomType() == SingularityPoint::VERTEX && points[1]->getGeomType() == SingularityPoint::VERTEX) {
			m_lineIsfixed[line->getNumber()] = true;
			m_groupIsFixed[line->getLinkedEdgeID()] = true;
		}
		else if (!m_groupIsFixed[line->getLinkedEdgeID()]) {
			m_lineIsfixed[line->getNumber()] = false;
		}
	}

	m_numberOfLinePerGroup = std::vector<size_t>(m_graph->getNLinkedEdgeGroup(), 0);
	for (const auto line : m_graph->getLines())
		if (!m_groupIsFixed[line->getLinkedEdgeID()] || m_lineIsfixed[line->getNumber()]) m_numberOfLinePerGroup[line->getLinkedEdgeID()]++;
}

void
LineRelaxator::computeMeanEdgeLengthByGroup()
{
	m_lineSpringLenght = std::vector<double>(m_graph->getNLinkedEdgeGroup(), 0.0);
	for (const auto line : m_graph->getLines()) {
		const auto points = line->getEndPoints();
		const auto direction = math::Vector3d(points[0]->getLocation(), points[1]->getLocation());

		if (!m_groupIsFixed[line->getLinkedEdgeID()]) {
			const auto groupID = line->getLinkedEdgeID();
			m_lineSpringLenght[groupID] += direction.norm();
		}
		else if (m_lineIsfixed[line->getNumber()]) {
			const auto groupID = line->getLinkedEdgeID();
			m_lineSpringLenght[groupID] += direction.norm();
		}
	};
	for (int groupID = 0; groupID < m_graph->getNLinkedEdgeGroup(); groupID++) {
		m_lineSpringLenght[groupID] /= static_cast<double>(m_numberOfLinePerGroup[groupID]);
	};
}

void
LineRelaxator::applySpringForceOnPoints()
{
	// compute critical lines : those whose length is furthest from their group average length
	m_maxLengthByGroup = std::vector<double>(m_graph->getNLinkedEdgeGroup(), 0);
	m_minLengthByGroup = std::vector<double>(m_graph->getNLinkedEdgeGroup(), 1e12);
	m_maxLineByGroup = std::vector<SingularityLine *>(m_graph->getNLinkedEdgeGroup(), nullptr);
	m_minLineByGroup = std::vector<SingularityLine *>(m_graph->getNLinkedEdgeGroup(), nullptr);
	for (const auto line : m_graph->getLines()) {
		const auto points = line->getEndPoints();
		const double length = points[0]->getLocation().distance(points[1]->getLocation());
		const int groupID = line->getLinkedEdgeID();
		if (length > m_maxLengthByGroup[groupID]) {
			m_maxLengthByGroup[groupID] = length;
			m_maxLineByGroup[groupID] = line;
		}
		if (length < m_minLengthByGroup[groupID]) {
			m_minLengthByGroup[groupID] = length;
			m_minLineByGroup[groupID] = line;
		}
	}

	m_maxLineByGroup.insert(m_maxLineByGroup.end(), m_minLineByGroup.begin(), m_minLineByGroup.end());

	for (const auto point : m_graph->getPoints())
		m_directions[point->getNumber()] = math::Vector3d();

	// stack forces from those critical line onto their respective points
	for (const auto line : m_graph->getLines()) {     // m_graph->getLines()

		auto points = line->getEndPoints();
		auto direction = math::Vector3d(points[0]->getLocation(), points[1]->getLocation());

		const auto length = direction.norm();
		const auto L0 = m_lineSpringLenght[line->getLinkedEdgeID()];
		const double weight = 1.0 - length / L0;

		direction /= length;
		m_directions[points[1]->getNumber()] += weight * direction;     // shrink when weight is positive
		m_directions[points[0]->getNumber()] += weight * direction.opp();
	};
};

void
LineRelaxator::MoveIfNewLocationIsWorthTheAnglePenalty(SingularityPoint *singPoint, const math::Point newLocation)
{
	const auto oldLocation = singPoint->getLocation();
	singPoint->setLocation(newLocation);
	const auto slots = singPoint->getSlots();
	for (auto slot : slots) {
		const auto lineSlots = slot->line->getSlots();
		const auto direction = math::Vector3d(lineSlots[0]->from_point->getLocation(), lineSlots[1]->from_point->getLocation());
		// here we "direction" or "direction.opp()" doesn't matter, as only the absolute value of "90-angle" will be used.
		lineSlots[0]->line_direction = direction;
		lineSlots[1]->line_direction = direction;
	}

	const bool isBdry = singPoint->getGeomType() != SingularityPoint::SURFACE;
	// 4 slots : 20°, 3 slots : 45°; 5 slots : 40°
	const double maxGap = slots.size() || isBdry == 4 ? 0.34906 : (slots.size() == 3 ? 0.785398 : 0.69813);

	const auto newAngles = computeAngleAtSingularity(singPoint);
	const auto &oldAngles = m_angleAtSingularity[singPoint->getNumber()];
	for (int i = 0; i < newAngles.size(); ++i) {
		const double newGap = fabs(M_PI_2 - newAngles[i]);
		if (newGap > maxGap) {     // 20°
			const double oldGap = fabs(M_PI_2 - oldAngles[i]) + 1e-8;
			if (newGap > oldGap) {
				// the move is cancel
				singPoint->setLocation(oldLocation);
				return;
			}
		}
	}
	// the move is validated, the new angles are stored
	m_angleAtSingularity[singPoint->getNumber()] = newAngles;
}

void
constrainDirection(SingularityPoint *p, math::Vector3d &d)
{
	const double length = d.norm();
	double minLinelength = 1e12;
	for (auto *slot : p->getSlots()) {
		const auto lineSlots = slot->line->getSlots();
		const double lineLength = lineSlots.front()->location.distance(lineSlots.back()->location);
		if (minLinelength > lineLength) minLinelength = lineLength;
	}
	const double maxAuthorizedLength = 0.3 * minLinelength;
	if (length > maxAuthorizedLength) d *= (maxAuthorizedLength / length);
}

void
LineRelaxator::movePoints()
{
	for (const auto point : m_graph->getPoints()) {

		switch (point->getGeomType()) {
		case SingularityPoint::SURFACE: {
			auto direction = m_directions[point->getNumber()] * m_step;
			constrainDirection(point, direction);
			const auto newLocation = point->getLocation() + direction;
			MoveIfNewLocationIsWorthTheAnglePenalty(point, newLocation);
			break;
		}
		case SingularityPoint::CURVE: {     // TODO : curved (bent) line. (SingularityPoint::CURVE should be renamed ::BOUNDARY..)
			SingularityLine *ref_line = nullptr;
			for (const auto line : point->getLines()) {
				if (line->getType() == SingularityLine::CURVE) {
					ref_line = line;
					break;
				}
			}
			if (!ref_line) throw GMDSException("curve singularity without curved line");

			const auto &linePoints = ref_line->getEndPoints();
			const math::Vector3d lineDir(linePoints[0]->getLocation(), linePoints[1]->getLocation());
			auto projectedDirection = m_directions[point->getNumber()].dot(lineDir) * lineDir / lineDir.norm2();
			projectedDirection *= m_step;
			constrainDirection(point, projectedDirection);
			const auto newLocation = point->getLocation() + projectedDirection;
			MoveIfNewLocationIsWorthTheAnglePenalty(point, newLocation);
		}

		default: break;
		}
	};
}

double
LineRelaxator::getWorstEdgeRatio()
{
	auto maxLengthByGroup = std::vector<double>(m_graph->getNLinkedEdgeGroup(), 0);
	auto minLengthByGroup = std::vector<double>(m_graph->getNLinkedEdgeGroup(), 1e12);
	for (const auto line : m_graph->getLines()) {
		const auto points = line->getEndPoints();
		const double length = points[0]->getLocation().distance(points[1]->getLocation());
		const int groupID = line->getLinkedEdgeID();
		if (length > maxLengthByGroup[groupID]) maxLengthByGroup[groupID] = length;
		if (length < minLengthByGroup[groupID]) minLengthByGroup[groupID] = length;
	}
	double worstRatio = 10000;
	for (int i = 0; i < m_graph->getNLinkedEdgeGroup(); i++) {
		double maxLength = maxLengthByGroup[i];
		double minLength = minLengthByGroup[i];
		double ratio = minLength / maxLength;
		if (ratio < worstRatio) worstRatio = ratio;
	}
	return worstRatio;
}

double
LineRelaxator::getWorstAngle()
{
	double worstDiff90 = 0;
	double worstAngle = M_PI_2;
	math::Point worstPoint;
	for (const auto patch : m_graph->getSurfacePatchs()) {

		std::vector<SingularityPoint *> points;
		patch->getPoints(points);
		const auto v1 = math::Vector3d(points[0]->getLocation(), points[1]->getLocation());
		const auto v2 = math::Vector3d(points[1]->getLocation(), points[2]->getLocation());
		const auto v3 = math::Vector3d(points[2]->getLocation(), points[3]->getLocation());
		const auto v4 = math::Vector3d(points[3]->getLocation(), points[0]->getLocation());
		std::vector<double> angles;
		angles.push_back(v1.angle(v4.opp()));
		angles.push_back(v2.angle(v1.opp()));
		angles.push_back(v3.angle(v2.opp()));
		angles.push_back(v4.angle(v3.opp()));
		int count = 0;
		for (const auto angle : angles) {
			double diff90 = std::abs(M_PI_2 - angle);
			if (diff90 > worstDiff90) {
				worstDiff90 = diff90;
				worstAngle = angle;
				worstPoint = points[count]->getLocation();
			}
			count++;
		}
	}
	// std::cout << worstPoint << std::endl;
	return worstAngle * 180 / M_PI;
}

void
LineRelaxator::run()
{
	std::cout << "Relax mesh to optimize edge ratios: " << std::endl;
	std::cout << "	worst angle/ratio before relaxation : " << getWorstAngle() << ", " << getWorstEdgeRatio() << std::endl;

	computeMeanEdgeLengthByGroup();
	for (int i = 0; i < 99; i++) {
		if ((i - 1) % 20 == 0) {
			computeMeanEdgeLengthByGroup();
			m_step *= 0.75;
		}
		applySpringForceOnPoints();
		movePoints();
	}
	std::cout << "	worst angle/ratio after  relaxation : " << getWorstAngle() << ", " << getWorstEdgeRatio() << std::endl;
}

}     // namespace gmds
