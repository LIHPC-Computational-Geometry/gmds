/*----------------------------------------------------------------------------*/
/*
 * SingGraphBuilder2DShortestPath.cpp
 *
 *  Created on: 13 juil. 2014
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
#include "gmds/math/Constants.h"
#include <gmds/frame/LaplaceCross2D.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
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
/*----------------------------------------------------------------------------*/
#include <gmds/singGraphBuild/SingGraphBuilder2DShortestPath.h>
#include <gmds/singGraphBuild/SingularityGraphBuilder2D.h>
#include <gmds/singGraphBuild/Tools.h>
/*----------------------------------------------------------------------------*/
#include <glpk.h>
/*----------------------------------------------------------------------------*/
#include <chrono>
#include <queue>
#include <set>
#include <sstream>
#include <unordered_set>
#include <vector>

/*----------------------------------------------------------------------------*/
namespace gmds {

namespace {

constexpr double M_MAXDIST = 1000000.0;
constexpr double maxAngleWithoutPenalty = M_PI / 4.5;
constexpr double turnPenalizationTerm = 0.1 * M_MAXDIST;

struct DijkstraCellParam
{
	double orientedCost = M_MAXDIST;
	double absoluteCost = M_MAXDIST;
	unsigned int nVisitedCellsBefore = 1;
	double penalty = 0;
	double getFullCost() const
	{
		return absoluteCost / nVisitedCellsBefore + penalty;     // for oriented cost, add : 0.125 * ... + fabs(orientedCost)
	}
};

class Timer
{
	using Clock = std::chrono::system_clock;

 public:
	Timer(const std::string &name)
	{
		m_name = name;
		m_t0 = Clock::now();
	}
	~Timer()
	{
		const auto tfinal = Clock::now();
		std::cout << m_name << " " << std::chrono::duration_cast<std::chrono::milliseconds>(tfinal - m_t0).count() << " milliseconds" << std::endl;
	}

 private:
	std::string m_name;
	std::chrono::system_clock::time_point m_t0;
};

void
retraceShortestPath(const TCellID &lastFaceAdded, const vector<TCellID> &previous, const TCellID &EndCell, vector<TCellID> &OUTPath)
{
	TCellID tempSource = lastFaceAdded;
	OUTPath.clear();

	while (tempSource != EndCell) {
		OUTPath.emplace_back(tempSource);
		tempSource = previous[tempSource];
	}
}

void
writeCostPerSlot(const Mesh *const mesh, std::string fileName, const vector<DijkstraCellParam> &costByCell)
{
	Mesh m(MeshModel(DIM3 | F | N | F2N));

	for (const TCellID n_id : mesh->nodes()) {
		const auto &n = mesh->get<Node>(n_id);
		m.newNode(n.point());
	}

	Variable<double> *cost = m.newVariable<double, GMDS_FACE>("cost");
	for (auto f_id : mesh->faces()) {
		const Face &f = mesh->get<Face>(f_id);
		m.newFace(f.get<Node>());

		double c = costByCell[f_id].getFullCost();
		(*cost)[f_id] = std::min(c, 1.0);
	}

	IGMeshIOService meshIoServ(&m);
	VTKWriter writer(&meshIoServ);
	writer.setCellOptions(N | F);
	writer.setDataOptions(N | F);
	writer.write(fileName);
}

void
writeShortestPathMesh(vector<SurfaceSingularityLine *> surf_lines, std::string file_name_paths)
{
	Mesh meshSingPaths(MeshModel(DIM3 | F | N | F2N));
	for (int i = 0; i < surf_lines.size(); i++) {
		vector<math::Point> surf_line_discretization = surf_lines[i]->getDiscretizationPoints();
		Node mySing0 = meshSingPaths.newNode(surf_line_discretization[0].X(), surf_line_discretization[0].Y(), surf_line_discretization[0].Z());

		for (unsigned int j = 1; j < surf_line_discretization.size(); j++) {
			Node mySing = meshSingPaths.newNode(surf_line_discretization[j].X(), surf_line_discretization[j].Y(), surf_line_discretization[j].Z());
			meshSingPaths.newTriangle(mySing0, mySing, mySing);
			mySing0 = mySing;
		}
	}

	IGMeshIOService ioServiceSingPaths(&meshSingPaths);
	VTKWriter vtkWriterSingPaths(&ioServiceSingPaths);
	vtkWriterSingPaths.setCellOptions(N | F);
	vtkWriterSingPaths.setDataOptions(N | F);
	vtkWriterSingPaths.write(file_name_paths);
}

void
writeGraphPoint(std::vector<SingularityPoint *> &points, std::string file_name_paths)
{
	Mesh meshSingPaths(MeshModel(DIM3 | F | N | F2N));
	for (const auto p : points) {

		Node mySing0 = meshSingPaths.newNode(p->getLocation());
		meshSingPaths.newTriangle(mySing0, mySing0, mySing0);
	}

	IGMeshIOService ioServiceSingPaths(&meshSingPaths);
	VTKWriter vtkWriterSingPaths(&ioServiceSingPaths);
	vtkWriterSingPaths.setCellOptions(N | F);
	vtkWriterSingPaths.setDataOptions(N | F);
	vtkWriterSingPaths.write(file_name_paths);
}

bool
segmentsOverlap(const math::Segment &ASeg1, const math::Segment &ASeg2)
{
	const auto &p11 = ASeg1.getPoint(0);
	const auto &p12 = ASeg1.getPoint(1);
	const auto &p21 = ASeg2.getPoint(0);
	const auto &p22 = ASeg2.getPoint(1);
	if (p11 == p21 || p11 == p22 || p12 == p21 || p12 == p22) {
		return true;
	}
	return false;
}

}     // namespace

/*----------------------------------------------------------------------------------------------------*/
/*----------------------    initialization/util functions     ----------------------------------------*/
/*----------------------------------------------------------------------------------------------------*/

SingGraphBuilder2DShortestPath::SingGraphBuilder2DShortestPath(Mesh *AMesh, Variable<math::Cross2D> *AField, const bool ABuildGeomSing) :
  SingularityGraphBuilder2D::SingularityGraphBuilder2D(AMesh, AField, ABuildGeomSing)
{
}

void
SingGraphBuilder2DShortestPath::writeCrossCenterTri()
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

void
SingGraphBuilder2DShortestPath::initializeFieldsValue()
{
	m_singOrGeomFaces = std::vector<bool>(m_original_faces_number, false);
	if (m_build_geometric_singularities) {
		for (const auto current_point : m_graph.getVertexPoints()) {
			for (const Face &cur_face : current_point->getMeshNode().get<Face>()) {
				m_singOrGeomFaces[cur_face.id()] = true;
			}
		}
	}
	for (const Face &sing_it : m_singularities_3) {
		m_singOrGeomFaces[sing_it.id()] = true;
	}
	for (const Face &sing_it : m_singularities_5) {
		m_singOrGeomFaces[sing_it.id()] = true;
	}

	std::vector<SingularityPoint *> singularity_points = m_graph.getPoints();
	m_singPointNo = m_graph.getPoints().size();
	m_totalNumberOfSlots = m_singPointNo * 5;
	m_totalNumberOfVariables = m_totalNumberOfSlots + m_maxNBoundaryTarget;

	// build a flat vector containing all slots
	m_targets = vector<Slot *>(m_totalNumberOfSlots);
	for (unsigned int i = 0; i < m_singPointNo; i++) {
		SingularityPoint *pi = singularity_points[i];
		std::vector<Slot *> pi_slots = pi->getSlots();
		for (unsigned j = 0; j < pi_slots.size(); j++) {
			SourceID contSource = 5 * i + j;
			m_targets[contSource] = pi_slots[j];
			m_targets[contSource]->direction.normalize();     // make sure that all direction are normalized
		}
	}

	m_possibleVariablesSlotsLaunched = vector<bool>(m_totalNumberOfSlots, false);
	for (SourceID contSource = 0; contSource < m_totalNumberOfSlots; contSource++) {
		const Slot *slot = m_targets[contSource];
		if (slot && !slot->isLaunched) {
			m_possibleVariablesSlotsLaunched[contSource] = true;
		}
	}
	m_distances = vector<vector<double>>(m_totalNumberOfSlots, vector<double>(m_totalNumberOfVariables, M_MAXDIST));
	m_bdryPathEndParam.clear();
}

void
SingGraphBuilder2DShortestPath::computeFaceNeighboursInfo()
{
	auto timer = Timer("computeFaceNeighborsInfo");
	m_face2Face_neighbours_by_verts.clear();
	m_face2Face_neighbours_by_verts.resize(m_original_faces_number);
	m_is_bdry_face.clear();
	m_is_bdry_face.resize(m_original_faces_number, false);

	std::unordered_set<TCellID> visitedFaces;     // faster than a vector<bool> because there are few visits, and vector<bool>.clear takes time
	visitedFaces.reserve(25);

	for (int f_id = 0; f_id < m_original_faces_number; f_id++) {

		const Face currentFace = m_mesh->get<Face>(f_id);
		const vector<Node> currentNodes = currentFace.get<Node>();
		const vector<Edge> currentEdges = currentFace.get<Edge>();

		for (const Node &node : currentNodes) {     // assume only triangles
			if (m_mesh->isMarked(node, m_mark_nodes_on_curve) || m_mesh->isMarked(node, m_mark_nodes_on_point)) {
				m_is_bdry_face[f_id] = true;
				break;
			}
		}
		for (const Node &node : currentNodes) {
			for (const auto &incident_triangle : node.get<Face>()) {
				bool notVisited = visitedFaces.find(incident_triangle.id()) == visitedFaces.end();
				if (notVisited && (incident_triangle.id() != f_id)) {
					bool validNeigh = true;

					if (m_mesh->isMarked(node, m_mark_nodes_on_point) || m_mesh->isMarked(node, m_mark_nodes_on_curve)) {

						validNeigh = [&]() -> bool {
							const math::Segment from_seg(m_triangle_centers[f_id], m_triangle_centers[incident_triangle.id()]);
							math::Point intersectionPoint;
							double intersectionParam;
							for (const Edge &edge : currentEdges) {
								if (m_mesh->isMarked(edge, m_mark_edges_on_curve)) {
									if (from_seg.SecondMetIntersect2D(edge.segment(), intersectionPoint, intersectionParam, m_temp_epsilon)) {
										return false;
									}
								}
							}
							for (const Edge &edge : incident_triangle.get<Edge>()) {
								if (m_mesh->isMarked(edge, m_mark_edges_on_curve)) {
									if (from_seg.SecondMetIntersect2D(edge.segment(), intersectionPoint, intersectionParam, m_temp_epsilon)) {
										return false;
									}
								}
							}
							return true;
						}();
					}
					if (validNeigh) {
						m_face2Face_neighbours_by_verts[f_id].emplace_back(incident_triangle);
					}
					visitedFaces.insert(incident_triangle.id());
				}
			}
		}
		visitedFaces.clear();
	}
	// For every border edge
	if (m_bdry_edge_normals.size() == 0) m_bdry_edge_normals.resize(m_mesh->getNbEdges());
	for (const TCellID e_id : m_mesh->edges()) {

		const Edge currentEdge = m_mesh->get<Edge>(e_id);
		if (m_mesh->isMarked(currentEdge, m_mark_edges_on_curve)) {

			const vector<Node> adjacent_nodes = currentEdge.get<Node>();
			const math::Point p1 = adjacent_nodes[0].point();
			const math::Point p2 = adjacent_nodes[1].point();
			const math::Vector3d v1 = (p2- p1).normalize();
			const vector<Face> adjacent_faces = currentEdge.get<Face>();
			m_bdry_edge_normals[e_id] = v1.cross(adjacent_faces[0].normal());
		}
	}
}

bool
SingGraphBuilder2DShortestPath::targetIsBoundary(const TargetID &contTarget) const
{
	return contTarget >= m_totalNumberOfSlots;
}

void
SingGraphBuilder2DShortestPath::setRobustness(int level)
{
	if (level <= 0) {
		m_maxNBoundaryTarget = 2;
		m_maxNValidSurfaceTarget = 2;
	}
	else if (level == 1) {
		m_maxNBoundaryTarget = 3;
		m_maxNValidSurfaceTarget = 5;
	}
	else if (level == 2) {
		m_maxNBoundaryTarget = 5;
		m_maxNValidSurfaceTarget = 8;
	}
	else {
		m_maxNBoundaryTarget = 15;
		m_maxNValidSurfaceTarget = 100;
	}
}

void
SingGraphBuilder2DShortestPath::setRobustness(unsigned int nBoundaryTarget, unsigned int nSurfaceTarget, unsigned int glpkTimeLimit)
{
	m_maxNBoundaryTarget = nBoundaryTarget;
	m_maxNValidSurfaceTarget = nSurfaceTarget;
	m_glpkTimeLimit = glpkTimeLimit;
}

void
SingGraphBuilder2DShortestPath::setGLPKTimeLimit(int glpkTimeLimit)
{
	m_glpkTimeLimit = glpkTimeLimit;
}
/*----------------------------------------------------------------------------------------------------*/
/*----------------------    Graph building functions     ---------------------------------------------*/
/*----------------------------------------------------------------------------------------------------*/
void
SingGraphBuilder2DShortestPath::tryEarlySlotsConnection()
{
	// store faces close to target
	vector<set<TargetID>> neighboringTarget(m_original_faces_number);
	for (SourceID contSource = 0; contSource < m_totalNumberOfSlots; contSource++) {
		const TCellID slotFace = m_slotFaces[contSource];
		if (slotFace < m_original_faces_number) {
			neighboringTarget[slotFace].insert(contSource);

			// add all incident faces to the slot nodes
			std::vector<Node> slotNodes;
			if (m_targets[contSource]->starting_cell_dim == 1)
				slotNodes = m_mesh->get<Edge>(m_targets[contSource]->starting_cell_id).get<Node>();
			else
				slotNodes = {m_mesh->get<Node>(m_targets[contSource]->starting_cell_id)};
			for (const Node &slotNode : slotNodes) {
				for (const TCellID &neighborFaceID : slotNode.getIDs<Face>()) {
					neighboringTarget[neighborFaceID].insert(contSource);
				}
			}
		}
	}

	// try an early connection between close slots
	for (SourceID contSource = 0; contSource < m_totalNumberOfSlots; contSource++) {
		const Slot *target = m_targets[contSource];
		if (target && !target->isFreeze) {
			const TCellID slotFace = m_slotFaces[contSource];
			for (const TargetID &contTarget : neighboringTarget[slotFace]) {
				if (tryEarlySlot2SlotConnection(contSource, contTarget)) {
					break;
				}
			}
		}
	}

	// try an early connection between a slot and a close boundary
	for (SourceID contSource = 0; contSource < m_totalNumberOfSlots; contSource++) {
		const Slot *slot = m_targets[contSource];
		if (slot && !slot->isFreeze) {

			TCellID slotFace = m_slotFaces[contSource];
			if (m_is_bdry_face[slotFace]) {

				const math::Vector3d prevDir = m_targets[contSource]->direction;
				const math::Ray ray(m_targets[contSource]->location, m_targets[contSource]->location + prevDir);
				BoundaryPathEndParam bdryParam;

				if (findBoundary(ray, slotFace, bdryParam)) {

					m_targets[contSource]->isFreeze = true;
					m_bdryPathEndParam[contSource].emplace_back(bdryParam);
					m_fixedVariables.push_back(make_pair(contSource, contSource));
					m_distances[contSource][m_totalNumberOfSlots] = 0;
				}
			}
		}
	}
}

bool
SingGraphBuilder2DShortestPath::tryEarlySlot2SlotConnection(const SourceID contSource, const TargetID contTarget)
{
	Slot *sourceSlot = m_targets[contSource];
	Slot *targetSlot = m_targets[contTarget];
	if (targetSlot->isFreeze || sourceSlot->from_point == targetSlot->from_point) return false;

	const auto &sourcePoint = sourceSlot->from_point->getLocation();
	const auto &targetPoint = targetSlot->from_point->getLocation();

	const auto linkVector = targetPoint-sourcePoint;
	if ((linkVector.angle(sourceSlot->direction) < 0.523598) && (linkVector.angle(targetSlot->direction.opp()) < 0.523598)) {     // pi/6

		sourceSlot->isFreeze = true;
		targetSlot->isFreeze = true;
		sourceSlot->line_direction = linkVector;
		targetSlot->line_direction = linkVector.opp();
		m_distances[contSource][contTarget] = 0;
		m_fixedVariables.push_back(make_pair(contSource, contTarget));
		return true;
	}
	return false;
}

void
SingGraphBuilder2DShortestPath::computeStreamLineFirstPoints()
{
	m_slotFaces = vector<TCellID>(m_totalNumberOfSlots, m_original_faces_number);

	// First discretization from singularity and slot positions
	for (SourceID contSource = 0; contSource < m_totalNumberOfSlots; contSource++) {
		Slot *slot = m_targets[contSource];

		if (!slot || slot->isFreeze) continue;

		if (slot->starting_cell_dim == 0) {     // the slot is located right on a node (rare case)

			const Node &node = m_mesh->get<Node>(slot->starting_cell_id);
			TCellID bestFace = 0;
			double bestCos = -1;
			for (const TCellID &faceID : node.getIDs<Face>()) {
				const auto slotToFace = m_triangle_centers[faceID]-slot->location;
				const double cos = slotToFace.dot(slot->direction);
				if (cos > bestCos) {
					bestCos = cos;
					bestFace = faceID;
				}
			}
			m_slotFaces[contSource] = bestFace;
		}
		else {     // the slot is located right on an edge
			const Edge &currentEdge = m_mesh->get<Edge>(slot->starting_cell_id);
			const auto &adjFaceIDs = currentEdge.getIDs<Face>();
			const bool singPointIsInfirstFace = m_tool.isPntInTri(slot->from_point->getLocation(), m_mesh->get<Face>(adjFaceIDs[0]));
			if (adjFaceIDs.size() == 1) {
				if (singPointIsInfirstFace) {
					m_slotFaces[contSource] = adjFaceIDs[0];
				}
			}
			else {
				if (singPointIsInfirstFace)
					m_slotFaces[contSource] = adjFaceIDs[1];
				else
					m_slotFaces[contSource] = adjFaceIDs[0];
			}
		}
	}
}

bool
SingGraphBuilder2DShortestPath::findBoundary(const math::Ray &ray, const TCellID currentFace, BoundaryPathEndParam &bdryParam)
{
	math::Point endPoint;
	double intersectionParam;

	bool hasBdryEdge = false;
	bool reachedValidBdry = false;

	const Face &u = m_mesh->get<Face>(currentFace);
	for (const Edge &currentEdge : u.get<Edge>()) {
		if (m_mesh->isMarked(currentEdge, m_mark_edges_on_curve)) {

			hasBdryEdge = true;
			if (ray.SecondMetIntersect2D(currentEdge.segment(), endPoint, intersectionParam, m_temp_epsilon)) {

				double cosAngle = ray.getDirUnit().dot(m_bdry_edge_normals[currentEdge.id()]);
				if (cosAngle > 0.7071) {     // -> max angle = pi/4
					reachedValidBdry = true;
					bdryParam.finalFace = currentFace;
					bdryParam.cellId = currentEdge.id();
					bdryParam.dim = 1;
					bdryParam.point = endPoint;
					return true;
				}
			}
		}
	}
	if (!hasBdryEdge) {
		for (const Node &currentNode : u.get<Node>()) {
			if (m_mesh->isMarked(currentNode, m_mark_nodes_on_point) || m_mesh->isMarked(currentNode, m_mark_nodes_on_curve)) {

				for (const Edge &currentEdge : currentNode.get<Edge>()) {
					if (m_mesh->isMarked(currentEdge, m_mark_edges_on_curve)) {

						if (ray.SecondMetIntersect2D(currentEdge.segment(), endPoint, intersectionParam, m_temp_epsilon)) {

							double cosAngle = ray.getDirUnit().dot(m_bdry_edge_normals[currentEdge.id()]);
							if (cosAngle > 0.7071) {     // -> max angle = pi/4
								vector<Face> adj_faces1 = currentEdge.get<Face>();
								bdryParam.finalFace = adj_faces1[0].id();
								bdryParam.cellId = currentEdge.id();
								bdryParam.dim = 1;
								bdryParam.point = endPoint;
								return true;
							}
						}
					}
				}
			}
		}
	}
	return false;
}

void
SingGraphBuilder2DShortestPath::createSingularityLines()
{
	std::cout << "createSingularityLiesShortestPaths" << std::endl;

	// initialize variables
	computeFaceNeighboursInfo();
	initializeFieldsValue();
	writeCrossCenterTri();

	// begin the line discretization (first two points of each line)
	computeStreamLineFirstPoints();
	tryEarlySlotsConnection();

	// find shortest paths
	computeDijkstra();

	// check all illegal pair of line
	computeIllegalOverlappingPaths();

	// solve
	m_solutions = glpkSolve();

	// build singularityLines from shortest paths
	createSolvedSingularityLines();

	writeShortestPathMesh(m_graph.getSurfaceLines(), std::string(m_output_directory_name + "-ShortestPaths_result.vtk"));
}

void
SingGraphBuilder2DShortestPath::createSolvedSingularityLines()
{
	for (const auto &solution : m_solutions) {

		const SourceID contSource = solution.first;
		const TargetID contTarget = solution.second;

		const bool lineIsOneSegment = m_finalPaths[contSource][contTarget].empty();     // we ignore also slot points

		Slot *current_slot = m_targets[contSource];

		if (targetIsBoundary(contTarget)) {

			SurfaceSingularityLine *surf_line = m_graph.newSurfaceLine();

			current_slot->line = surf_line;
			current_slot->isLaunched = true;
			current_slot->isFreeze = true;
			current_slot->lineDeviation = m_distances[contSource][contTarget];

			SingularityPoint *geom_pnt = nullptr;
			Slot *incoming_slot = nullptr;
			const unsigned int bdryLocalID = contTarget - m_totalNumberOfSlots;
			const auto &boundaryParam = m_bdryPathEndParam[contSource][bdryLocalID];
			const auto &penultimatePoint =
			   lineIsOneSegment ? current_slot->from_point->getLocation() : m_triangle_centers[m_finalPaths[contSource][contTarget].back()];
			const auto lastLineDirection = (boundaryParam.point-penultimatePoint).normalize();

			createGeometricSingularityPoint(boundaryParam.point,      // the last point added
			                                lastLineDirection,        // the direction
			                                boundaryParam.dim,        // the dim. of the cell
			                                boundaryParam.cellId,     // the id of the cell
			                                geom_pnt,                 // the created point
			                                incoming_slot);           // and the slot

			incoming_slot->isLaunched = true;
			surf_line->addSlot(current_slot);
			surf_line->addSlot(incoming_slot);

			// add only two points in the discretization if the solution was fixed before GLPK
			surf_line->addDiscretizationPoint(current_slot->from_point->getLocation());
			if (!lineIsOneSegment) {
				surf_line->addDiscretizationPoint(current_slot->location);
				for (const TCellID &faceId : m_finalPaths[contSource][contTarget]) {
					surf_line->addDiscretizationPoint(m_triangle_centers[faceId]);
				}
			}
			surf_line->addDiscretizationPoint(boundaryParam.point);

			if (m_visualizeShortestPaths && surf_line->getDiscretizationPoints().size() >= 2) {
				std::string file_name = "ShortestPathsBdryEnd_" + std::to_string(contSource) + "-" + std::to_string(contTarget) + ".vtk";
				writeTestPoints(surf_line->getDiscretizationPoints(), file_name);
			}
		}
		else {     // not a boundary sing line

			SurfaceSingularityLine *surf_line = m_graph.newSurfaceLine();

			// first slot (line direction is not set here)
			current_slot->line = surf_line;
			current_slot->isLaunched = true;
			current_slot->isFreeze = true;
			current_slot->lineDeviation = m_distances[contSource][contTarget];

			// second slot
			Slot *to_slot = m_targets[contTarget];
			to_slot->line = surf_line;
			to_slot->isLaunched = true;
			to_slot->isFreeze = true;
			to_slot->lineDeviation = m_distances[contSource][contTarget];

			surf_line->addSlot(current_slot);
			surf_line->addSlot(to_slot);

			// add only two points in the discretization if the solution was fixed before GLPK
			surf_line->addDiscretizationPoint(current_slot->from_point->getLocation());
			if (!lineIsOneSegment) {
				surf_line->addDiscretizationPoint(current_slot->location);
				for (const TCellID &faceId : m_finalPaths[contSource][contTarget]) {
					surf_line->addDiscretizationPoint(m_triangle_centers[faceId]);
				}
				surf_line->addDiscretizationPoint(to_slot->location);
			}
			surf_line->addDiscretizationPoint(to_slot->from_point->getLocation());

			if (m_visualizeShortestPaths && surf_line->getDiscretizationPoints().size() >= 2) {
				std::string file_name = "ShortestPathsEnd_" + std::to_string(contSource) + "-" + std::to_string(contTarget) + ".vtk";
				writeTestPoints(surf_line->getDiscretizationPoints(), file_name);
			}
		}
	}
}

/*----------------------------------------------------------------------------------------------------*/
/*----------------------    Shortest path system functions     ---------------------------------------*/
/*----------------------------------------------------------------------------------------------------*/

void
SingGraphBuilder2DShortestPath::computeDijkstra()
{
	auto timer = Timer("graph construction");

	using TCellID2DMatrix = vector<vector<TCellID>>;
	m_finalPaths = vector<TCellID2DMatrix>(m_totalNumberOfSlots, TCellID2DMatrix(m_totalNumberOfVariables, vector<TCellID>()));

	vector<TCellID> targetFaceOpt(m_totalNumberOfSlots, m_original_faces_number);
	for (SourceID contSource = 0; contSource < m_totalNumberOfSlots; contSource++) {
		if (m_targets[contSource] && !m_targets[contSource]->isFreeze) targetFaceOpt[contSource] = m_slotFaces[contSource];
	}

	for (unsigned int i = 0; i < m_singPointNo; i++) {
		for (unsigned int j = 0; j < 5; j++) {

			unsigned int contSource = 5 * i + j;
			if (m_targets[contSource] && !m_targets[contSource]->isFreeze) {

				// prevent slot from connecting to a slot belonging to the same singularity
				vector<TCellID> tempTargetFaceOpt(5);
				for (int jk = 0; jk < 5; jk++) {
					tempTargetFaceOpt[jk] = targetFaceOpt[5 * i + jk];
					targetFaceOpt[5 * i + jk] = m_original_faces_number + 1;
				}

				getShortestPathBtwFacesOptimized(targetFaceOpt, contSource);

				for (int jk = 0; jk < 5; jk++) {
					targetFaceOpt[5 * i + jk] = tempTargetFaceOpt[jk];
				}
			}
		}
	}
}

void
SingGraphBuilder2DShortestPath::getShortestPathBtwFacesOptimized(const vector<TCellID> &targets, const SourceID &contSource)
{

	/*Description: get the shortest path between a singular triangle \param[in] source
	 * and all other singular triangles \param[in] targets by 'walking from face center to face center'*/

	const math::Point startPoint = m_targets[contSource]->location;
	const TCellID source = m_slotFaces[contSource];
	std::set<TargetID> validTargetReached;

	// compute face2target map : TODO  should be the same for all sources, when cycling singularity lines will be allowed
	vector<vector<TargetID>> neighbouringTarget(m_original_faces_number, vector<TargetID>(0));     // face2Target map
	for (TargetID i = 0; i < targets.size(); i++) {
		if (targets[i] < m_original_faces_number) {
			neighbouringTarget[targets[i]].push_back(i);
			for (const Face &neighbouringFace : m_face2Face_neighbours_by_verts[targets[i]])
				neighbouringTarget[neighbouringFace.id()].push_back(i);
		}
	}

	// previous direction
	const math::Vector3d sourceDirection = m_targets[contSource]->direction;

	// vector to keep track of each path
	vector<TCellID> previousFaceID(m_original_faces_number, m_original_faces_number);     // invalid cell ids for now
	vector<TCellID> retracedPath;
	vector<bool> visitedFaces(m_original_faces_number, false);
	vector<math::Vector3d> previousCrossDirection(m_original_faces_number, sourceDirection);

	// Dijkstra queue and cost
	std::set<std::pair<double, TCellID>> vertex_queue;
	vector<DijkstraCellParam> cellCosts(m_original_faces_number);

	// boundary lines set up
	std::priority_queue<std::pair<double, unsigned int>> boundaryQueue;
	std::vector<int> bdryFacesId;

	// search for nearby target, if found, only the singularities locations will be used,
	// similar as in tryEarlySlot2SlotConnection(), but less strict on angle
	for (const TargetID contTarget : neighbouringTarget[source]) {
		const auto linkVector = m_targets[contTarget]->from_point->getLocation()-m_targets[contSource]->from_point->getLocation();
		const auto targetOwnDirection = m_targets[contTarget]->direction.opp();
		const double srcAngle = linkVector.angle(sourceDirection);
		const double tgtAngle = linkVector.angle(targetOwnDirection);
		if (srcAngle < M_PI_4 && tgtAngle < M_PI_4) m_distances[contSource][contTarget] = srcAngle + tgtAngle;
	}

	// check if the slot face is actually a good candidate
	{
		const auto tri2tri = (m_triangle_centers[source]-startPoint).normalize();
		const double directionnalAngle = sourceDirection.angle(tri2tri);

		if (directionnalAngle < math::Constants::PIDIV2) {     // less strict here (pi/2), to manage turbulant cross field near singularities

			const math::Cross2D crossV = m_triangle_centers_cross[source];
			const math::Vector3d closestToPrevCross = crossV.closestComponentVector(sourceDirection);
			cellCosts[source].orientedCost = 0;
			cellCosts[source].absoluteCost = directionnalAngle + sourceDirection.angle(closestToPrevCross);
			cellCosts[source].penalty = directionnalAngle > math::Constants::PIDIV4 ? turnPenalizationTerm : 0.0;			
			previousCrossDirection[source] = closestToPrevCross;
			vertex_queue.insert(std::make_pair(cellCosts[source].getFullCost(), source));
		}
	}

	// Also initialize the source face neighbors
	for (const Face neighbour : m_face2Face_neighbours_by_verts[source]) {

		const TCellID v_id = neighbour.id();

		if (!m_mesh->isMarked<Face>(v_id, m_mark_faces_with_sing_point)) {

			const auto tri2tri = (m_triangle_centers[v_id]-startPoint).normalize();
			const double directionnalAngle = sourceDirection.angle(tri2tri);

			if (directionnalAngle < math::Constants::PIDIV4) {

				const math::Cross2D crossV = m_triangle_centers_cross[v_id];
				const math::Vector3d closestToPrevCross = crossV.closestComponentVector(sourceDirection);
				cellCosts[v_id].orientedCost = 0;
				cellCosts[v_id].absoluteCost = directionnalAngle + sourceDirection.angle(closestToPrevCross);
				cellCosts[v_id].penalty += m_singOrGeomFaces[v_id] ? turnPenalizationTerm : 0.0;
				// cellCosts[v_id].penalty = directionnalAngle > math::Constants::PIDIV4 ? turnPenalizationTerm : 0.0;
				previousCrossDirection[v_id] = closestToPrevCross;
				vertex_queue.insert(std::make_pair(cellCosts[v_id].getFullCost(), v_id));
			}
		}
	}

	while (!vertex_queue.empty()) {

		const TCellID u_id = vertex_queue.begin()->second;
		auto previousCostParam = cellCosts[u_id];

		vertex_queue.erase(vertex_queue.begin());
		visitedFaces[u_id] = true;

		const math::Vector3d prevCross_u = previousCrossDirection[u_id];

		// search for nearby target
		for (const TargetID contTarget : neighbouringTarget[u_id]) {

			const auto targetOwnDirection = m_targets[contTarget]->direction.opp();

			if (targetOwnDirection.angle(prevCross_u) > math::Constants::PIDIV4) continue;

			const auto tri2tri = (m_targets[contTarget]->location-m_triangle_centers[u_id]).normalize();
			const double directionnalAngle = prevCross_u.angle(tri2tri) + targetOwnDirection.angle(tri2tri);
			DijkstraCellParam targetCost = previousCostParam;     // ++nVisitedCellsBefore
			targetCost.nVisitedCellsBefore++;
			targetCost.absoluteCost += directionnalAngle + targetOwnDirection.angle(prevCross_u);
			if (directionnalAngle > maxAngleWithoutPenalty) targetCost.penalty += 3 * turnPenalizationTerm;

			const double newCost = targetCost.getFullCost();
			if (newCost < m_distances[contSource][contTarget]) {

				m_distances[contSource][contTarget] = targetCost.getFullCost();
				retraceShortestPath(u_id, previousFaceID, m_original_faces_number, retracedPath);
				m_finalPaths[contSource][contTarget].clear();
				if (!retracedPath.empty())
					m_finalPaths[contSource][contTarget].insert(m_finalPaths[contSource][contTarget].end(), retracedPath.rbegin(), retracedPath.rend());
				if (newCost < turnPenalizationTerm) {
					validTargetReached.insert(contTarget);
					if (validTargetReached.size() >= m_maxNValidSurfaceTarget) return;
				}
			}
			previousCostParam.penalty += turnPenalizationTerm / 100; // passing nearby a target incurs a small penalty
		}

		// search for boundary
		if (m_is_bdry_face[u_id]) {

			math::Ray ray(m_triangle_centers[u_id], prevCross_u);
			BoundaryPathEndParam bdryParam;
			bool reachedValidBdry = findBoundary(ray, u_id, bdryParam);

			if (reachedValidBdry) {
				const bool faceIsAlreadyEndOfAnotherLine = false;
				for (const auto face : bdryFacesId)
					if (bdryParam.finalFace == face) {
						reachedValidBdry = false;
						break;
					}
			}
			if (reachedValidBdry) {

				double newCost = previousCostParam.getFullCost();

				retraceShortestPath(u_id, previousFaceID, m_original_faces_number, retracedPath);
				for (const auto &faceId : retracedPath)
					if (m_singOrGeomFaces[faceId]) newCost += 2 * turnPenalizationTerm;

				bool addNewBoundary = false;
				unsigned int bdryLocalId;
				if (boundaryQueue.size() < m_maxNBoundaryTarget) {

					addNewBoundary = true;
					// update cost and boundary parameter
					bdryLocalId = boundaryQueue.size();
					m_bdryPathEndParam[contSource].emplace_back(bdryParam);
					bdryFacesId.push_back(bdryParam.finalFace);
				}
				else if (newCost < boundaryQueue.top().first) {

					addNewBoundary = true;
					bdryLocalId = boundaryQueue.top().second;
					m_bdryPathEndParam[contSource][bdryLocalId] = bdryParam;
					bdryFacesId[bdryLocalId] = bdryParam.finalFace;

					boundaryQueue.pop();
				}
				if (addNewBoundary) {
					// update cost
					m_distances[contSource][m_totalNumberOfSlots + bdryLocalId] = newCost;
					// update boudraryQueue
					boundaryQueue.emplace(newCost, bdryLocalId);
					// update path
					vector<TCellID> &path = m_finalPaths[contSource][m_totalNumberOfSlots + bdryLocalId];
					path.clear();
					if (retracedPath.size() > 0) path.insert(path.end(), retracedPath.rbegin(), retracedPath.rend());
				}
			}
		}

		// search among neighbours for next step
		for (const Face &neighbourFace : m_face2Face_neighbours_by_verts[u_id]) {

			const TCellID v_id = neighbourFace.id();
			const bool validNeighbour = m_mesh->isMarked<Face>(v_id, m_mark_faces_with_sing_point) ? false : true;

			if (!visitedFaces[v_id] && validNeighbour) {
				const auto tri2tri = (m_triangle_centers[v_id]-m_triangle_centers[u_id]).normalize();
				const double directionnalAngle = prevCross_u.orientedAngle(tri2tri);
				const double absDirectionnalAngle = fabs(directionnalAngle);

				DijkstraCellParam neighbourCellCost = previousCostParam;
				if (absDirectionnalAngle > M_PI_4 || m_singOrGeomFaces[v_id]) {
					continue;
				}
				else if (absDirectionnalAngle > maxAngleWithoutPenalty) {
					neighbourCellCost.penalty += turnPenalizationTerm;
				}

				neighbourCellCost.nVisitedCellsBefore++;
				neighbourCellCost.orientedCost += directionnalAngle;
				neighbourCellCost.absoluteCost += absDirectionnalAngle;

				const double newCost = neighbourCellCost.getFullCost();
				const double oldCost = cellCosts[v_id].getFullCost();
				if (newCost < oldCost) {

					if (cellCosts[v_id].orientedCost != M_MAXDIST) vertex_queue.erase(std::make_pair(oldCost, v_id));     // reprioritize
					cellCosts[v_id] = neighbourCellCost;
					previousFaceID[v_id] = u_id;
					const math::Cross2D &crossV = m_triangle_centers_cross[v_id];
					previousCrossDirection[v_id] = crossV.closestComponentVector(prevCross_u);
					vertex_queue.insert(std::make_pair(newCost, v_id));
				}
			}
		}
	}

	if (m_visualizeCost) writeCostPerSlot(m_mesh, m_output_directory_name + " -slotCost_" + std::to_string(contSource) + ".vtk", cellCosts);	
}

std::vector<std::pair<unsigned int, unsigned int>>
SingGraphBuilder2DShortestPath::glpkSolve()
{
	std::vector<std::pair<SourceID, TargetID>> solutions;
	auto timer = Timer("glpk");
	glp_prob *lp = glp_create_prob();

	int *ia, *ja;
	double *ar;

	glp_set_prob_name(lp, "SingGraphBuildOpt");
	glp_set_obj_dir(lp, GLP_MIN);

	// objective function //WARNING TODO if non launched slots -> their sum should be 0, not 2; or dont put condition
	unsigned int varNumber = m_totalNumberOfSlots * m_totalNumberOfVariables;
	glp_add_cols(lp, varNumber);

	vector<bool> fixedZeroVars(varNumber, false);

	for (unsigned int i = 0; i < m_singPointNo; i++) {
		for (unsigned int j = 0; j < 5; j++) {
			if (!m_possibleVariablesSlotsLaunched[5 * i + j]) {
				// entire lines fixed to 0
				for (unsigned int k = 0; k < m_totalNumberOfVariables; k++) {
					unsigned int tempVar = (5 * i + j) * m_totalNumberOfVariables + k;
					if (!fixedZeroVars[tempVar]) {

						fixedZeroVars[tempVar] = true;
						glp_set_col_bnds(lp, tempVar + 1, GLP_FX, 0.0, 0.0);
						glp_set_obj_coef(lp, tempVar + 1, 0);
					}
				}
				// columns fixed to 0
				for (unsigned int k = 0; k < m_totalNumberOfSlots; k++) {
					unsigned int tempVar = m_totalNumberOfVariables * k + (5 * i + j);
					if (!fixedZeroVars[tempVar]) {
						fixedZeroVars[tempVar] = true;
						glp_set_col_bnds(lp, tempVar + 1, GLP_FX, 0.0, 0.0);
						glp_set_obj_coef(lp, tempVar + 1, 0);
					}
				}
			}
		}
	}

	// CONDITION: slots corresponding to same singularity are not allowed to be connected in between themselves    ; however this condition should be removed
	for (unsigned int i = 0; i < m_totalNumberOfSlots; i++) {
		for (unsigned int j = 0; j < 5; j++) {
			unsigned int tempVar = m_totalNumberOfVariables * i + j;
			if (!fixedZeroVars[tempVar]) {
				fixedZeroVars[tempVar] = true;
				glp_set_col_bnds(lp, tempVar + 1, GLP_FX, 0.0, 0.0);
				glp_set_obj_coef(lp, tempVar + 1, 0);
			}
		}
	}

	// fix variable
	vector<pair<unsigned int, unsigned int>> fixedVariables(m_fixedVariables.size());
	for (unsigned int l1 = 0; l1 < m_fixedVariables.size(); l1++) {

		if (m_fixedVariables[l1].first == m_fixedVariables[l1].second) {
			fixedVariables[l1].first = m_fixedVariables[l1].first * m_totalNumberOfVariables + m_totalNumberOfSlots;
			fixedVariables[l1].second = fixedVariables[l1].first;
		}
		else {
			fixedVariables[l1].first = m_fixedVariables[l1].first * m_totalNumberOfVariables + m_fixedVariables[l1].second;
			fixedVariables[l1].second = m_fixedVariables[l1].second * m_totalNumberOfVariables + m_fixedVariables[l1].first;
		}
	}

	vector<double> w(m_totalNumberOfSlots * m_totalNumberOfVariables);
	for (unsigned int i = 0; i < m_distances.size(); i++) {
		for (unsigned int j = 0; j < m_distances[i].size(); j++) {
			w[i * m_totalNumberOfVariables + j] = m_distances[i][j];
		}
	}

	vector<bool> setFixedVars(fixedVariables.size(), false);
	for (unsigned int i = 0; i < m_totalNumberOfSlots * m_totalNumberOfVariables; i++) {
		bool setVal = false;

		for (unsigned int j = 0; j < fixedVariables.size(); j++) {
			if ((fixedVariables[j].first == i) || (fixedVariables[j].second == i)) {
				setVal = true;
				if (!setFixedVars[j]) {
					setFixedVars[j] = true;
					glp_set_col_bnds(lp, fixedVariables[j].first + 1, GLP_FX, 1.0, 1.0);
					glp_set_obj_coef(lp, fixedVariables[j].first + 1, 0.0001);
				}
			}
		}

		if (!setVal) {
			if (w[i] == M_MAXDIST) {
				glp_set_col_bnds(lp, i + 1, GLP_FX, 0.0, 0.0);
			}
			else {
				glp_set_col_kind(lp, i + 1, GLP_BV);     // GLP_IV - integer, but in our case we want binary
				glp_set_obj_coef(lp, i + 1, w[i]);
			}
		}
	}

	glp_add_rows(lp, m_totalNumberOfSlots + m_illegalOverlappingPaths.size());     //!!!!!!!!!!!

	// first equality constraints: all valid slots must have exactly 1 singularity line associated
	for (unsigned int i = 0; i < m_totalNumberOfSlots; i++) {
		if (m_possibleVariablesSlotsLaunched[i]) {
			glp_set_row_bnds(lp, i + 1, GLP_FX, 1, 1);
		}
		else
			glp_set_row_bnds(lp, i + 1, GLP_FX, 0, 0);
	}

	// second equality constraints: all IllegalCross pairs must have a maximum sum of 1 (at most 1 of them should be present in the final solution)
	for (unsigned int i = m_totalNumberOfSlots; i < m_totalNumberOfSlots + m_illegalOverlappingPaths.size(); i++) {
		glp_set_row_bnds(lp, i + 1, GLP_LO, 0.0, 0);
		glp_set_row_bnds(lp, i + 1, GLP_UP, 0.0, 1);
	}

	// dynamic allocation of the matrix where ia represents the row index, ja the col index and ar the entry
	ia = (int *) calloc(2 + (2 * varNumber + m_illegalOverlappingPaths.size() * 2), sizeof(int));
	ja = (int *) calloc(2 + (2 * varNumber + m_illegalOverlappingPaths.size() * 2), sizeof(int));
	ar = (double *) calloc(2 + (2 * varNumber + m_illegalOverlappingPaths.size() * 2), sizeof(double));

	unsigned int count = 1;

	for (unsigned int i = 0; i < m_totalNumberOfSlots; i++) {
		for (unsigned int j = 0; j < m_totalNumberOfVariables; j++) {
			ia[count] = i + 1;
			ja[count] = m_totalNumberOfVariables * i + j + 1;
			ar[count] = 1;
			count++;
			if ((j < m_totalNumberOfSlots) && (i != j)) {
				ia[count] = j + 1;
				ja[count] = m_totalNumberOfVariables * i + j + 1;
				ar[count] = 1;
				count++;
			}
		}
	}

	for (unsigned int i = 0; i < m_illegalOverlappingPaths.size(); i++) {
		ia[count] = m_totalNumberOfSlots + i + 1;
		ja[count] = m_illegalOverlappingPaths[i].first + 1;
		ar[count] = 1;
		count++;

		ia[count] = m_totalNumberOfSlots + i + 1;
		ja[count] = m_illegalOverlappingPaths[i].second + 1;
		ar[count] = 1;
		count++;
	}

	glp_load_matrix(lp, count - 1, ia, ja, ar);

	// the next line disables the default terminal report
	glp_term_out(GLP_ON);
	glp_iocp glpParams;
	glp_init_iocp(&glpParams);
	glpParams.tm_lim = m_glpkTimeLimit;
	glpParams.presolve = GLP_ON;

	// glp_write_lp(lp, NULL, "checkMe0.txt");

	switch (glp_intopt(lp, &glpParams)) {
	case GLP_ETMLIM:
	case 0: std::cout << "GLP OPT OK" << std::endl; break;
	default: throw GMDSException("SingularityGraphBuilder2D::SingGraphBuildOpt pb solving in GLPK.");
	}
	switch (glp_mip_status(lp)) {
	case GLP_UNDEF:
	case GLP_NOFEAS:
	case 0: throw GMDSException("SingularityGraphBuilder2D::SingGraphBuildOpt pb solving in GLPK.");
	default: std::cout << "MIP OK" << std::endl; break;
	}
	// glp_print_sol(lp, "SingGraphBuildOpt.txt");

	// for (unsigned int i = 0; i < varNumber; i++) {
	//	if (glp_mip_col_val(lp, i + 1) != 0) {
	//		std::cout << (int) (i / m_totalNumberOfVariables) << " " << fmod(i, m_totalNumberOfVariables) << std::endl;
	//	}
	//}

	// write solution
	solutions.reserve(m_totalNumberOfSlots);

	for (unsigned int i = 0; i < varNumber; i++) {
		if (glp_mip_col_val(lp, i + 1) != 0) {
			unsigned int contSource = (int) (i / m_totalNumberOfVariables);
			unsigned int contTarget = fmod(i, m_totalNumberOfVariables);
			solutions.emplace_back(std::make_pair(contSource, contTarget));
		}
	}
	free(ia);
	free(ja);
	free(ar);
	glp_delete_prob(lp);

	return solutions;
}

/*----------------------------------------------------------------------------------------------------*/
/*----------------------    illegal crossing functions     -------------------------------------------*/
/*----------------------------------------------------------------------------------------------------*/
void
SingGraphBuilder2DShortestPath::computeIllegalOverlappingPaths()
{
	auto illegalLineCrossingFinder = IllegalLineCrossingFinder(this);
	illegalLineCrossingFinder.registerAllLineSegmentsOnCells();
	m_illegalOverlappingPaths = illegalLineCrossingFinder.getIllegalOverlappingPaths();
	std::cout << "Number of illegal cross : " << m_illegalOverlappingPaths.size() << std::endl;
}

void
SingGraphBuilder2DShortestPath::IllegalLineCrossingFinder::registerAllLineSegmentsOnCells()
{
	auto timer = Timer("traversed faces");

	for (SourceID contSource = 0; contSource < m_graphBuilder->m_totalNumberOfSlots; contSource++) {

		const Slot *sourceSlot = m_graphBuilder->m_targets[contSource];
		if (!sourceSlot) continue;

		const TCellID sourceFace = m_graphBuilder->m_slotFaces[contSource];
		const auto startPoint = sourceSlot->location;

		for (TargetID contTarget = 0; contTarget < m_graphBuilder->m_totalNumberOfSlots; contTarget++) {

			if (m_graphBuilder->m_distances[contSource][contTarget] != M_MAXDIST) {

				const TCellID targetFace = m_graphBuilder->m_slotFaces[contTarget];
				registerOneLineSegment(startPoint, m_graphBuilder->m_targets[contTarget]->location, sourceFace, targetFace, contSource, contTarget);
			}
		}

		for (int bdryLocalId = 0; bdryLocalId < m_graphBuilder->m_bdryPathEndParam[contSource].size(); bdryLocalId++) {

			const TargetID contTarget = m_graphBuilder->m_totalNumberOfSlots + bdryLocalId;
			const auto &boundaryParam = m_graphBuilder->m_bdryPathEndParam[contSource][bdryLocalId];
			registerOneLineSegment(startPoint, boundaryParam.point, sourceFace, boundaryParam.finalFace, contSource, contTarget);
		}
	}
}

void
SingGraphBuilder2DShortestPath::IllegalLineCrossingFinder::registerOneLineSegment(const math::Point &startPnt,
                                                                                  const math::Point &endPnt,
                                                                                  const TCellID &startFace,
                                                                                  const TCellID &endFace,
                                                                                  const SourceID &contSource,
                                                                                  const TargetID &contTarget)
{
	const vector<TCellID> &facePath = m_graphBuilder->m_finalPaths[contSource][contTarget];

	if (facePath.empty()) {     // -> fixed variable
		const auto &p1 = m_graphBuilder->m_targets[contSource]->from_point->getLocation();
		math::Point p2;
		const int boundaryID = contTarget - m_graphBuilder->m_totalNumberOfSlots;
		if (boundaryID >= 0)
			p2 = m_graphBuilder->m_bdryPathEndParam[contSource][boundaryID].point;
		else
			p2 = m_graphBuilder->m_targets[contTarget]->from_point->getLocation();

		const auto segment = math::Segment(p1, p2);     // -> one and only segment of the path
		std::vector<TCellID> cellsToRegister;
		for (const Face &f : m_graphBuilder->m_face2Face_neighbours_by_verts[startFace])
			cellsToRegister.push_back(f.id());
		for (const Face &f : m_graphBuilder->m_face2Face_neighbours_by_verts[endFace])
			cellsToRegister.push_back(f.id());
		cellsToRegister.erase(unique(cellsToRegister.begin(), cellsToRegister.end()), cellsToRegister.end());
		for (const TCellID faceID : cellsToRegister)
			m_traversedFacesByLine[faceID].emplace_back(contSource, contTarget, segment);
		return;
	}

	// faces btw source slot and first face of the path
	registerLineSegmentBetweenTwoFaces(startFace, facePath[0], startPnt, m_graphBuilder->m_triangle_centers[facePath[0]], contSource, contTarget);

	// faces btw face path faces
	for (int i = 0; i < facePath.size() - 1; i++) {
		const TCellID f1 = facePath[i];
		const TCellID f2 = facePath[i + 1];
		const math::Point &p1 = m_graphBuilder->m_triangle_centers[f1];
		const math::Point &p2 = m_graphBuilder->m_triangle_centers[f2];
		registerLineSegmentBetweenTwoFaces(f1, f2, p1, p2, contSource, contTarget);
	}

	// faces btw last path face and target
	registerLineSegmentBetweenTwoFaces(facePath.back(), endFace, m_graphBuilder->m_triangle_centers[facePath.back()], endPnt, contSource, contTarget);

	// add surrounding faces for face on the path
	math::Segment firstSegment;
	if (facePath.size() == 1)
		firstSegment = math::Segment(startPnt, endPnt);
	else
		firstSegment = math::Segment(startPnt, m_graphBuilder->m_triangle_centers[facePath[1]]);

	if (facePath.front() != startFace) {
		firstSegment = math::Segment(startPnt, m_graphBuilder->m_triangle_centers[facePath[0]]);
		if (!m_graphBuilder->m_tool.isAdjacency(facePath.front(), startFace)) {
			const Node &n = Tools::getCommonNode(m_graphBuilder->m_mesh->get<Face>(facePath.front()), m_graphBuilder->m_mesh->get<Face>(startFace));
			for (const auto f : n.getIDs<Face>())     // only faces between f1 & f2 are needed in reality
				m_traversedFacesByLine[f].emplace_back(contSource, contTarget, firstSegment);
		}
		else {
			m_traversedFacesByLine[startFace].emplace_back(contSource, contTarget, firstSegment);
			m_traversedFacesByLine[facePath.front()].emplace_back(contSource, contTarget, firstSegment);
		}
	}
	else {
		m_traversedFacesByLine[startFace].emplace_back(contSource, contTarget, firstSegment);
	}

	for (int i = 1; i < facePath.size() - 1; i++) {
		const math::Point &p1 = m_graphBuilder->m_triangle_centers[facePath[i - 1]];
		const math::Point &p2 = m_graphBuilder->m_triangle_centers[facePath[i + 1]];
		m_traversedFacesByLine[facePath[i]].emplace_back(contSource, contTarget, math::Segment(p1, p2));
	}

	const math::Point &lastCenterFace = m_graphBuilder->m_triangle_centers[facePath.back()];
	const auto lastSegment = math::Segment(lastCenterFace, endPnt);

	if (facePath.back() != endFace) {
		m_traversedFacesByLine[endFace].emplace_back(contSource, contTarget, math::Segment(lastCenterFace, endPnt));

		if (!m_graphBuilder->m_tool.isAdjacency(facePath.back(), endFace)) {
			const Node &n = Tools::getCommonNode(m_graphBuilder->m_mesh->get<Face>(facePath.back()), m_graphBuilder->m_mesh->get<Face>(endFace));
			for (const auto f : n.getIDs<Face>())     // only faces between f1 & f2 are needed in reality
				m_traversedFacesByLine[f].emplace_back(contSource, contTarget, lastSegment);
		}
		else {
			m_traversedFacesByLine[endFace].emplace_back(contSource, contTarget, lastSegment);
			m_traversedFacesByLine[facePath.back()].emplace_back(contSource, contTarget, lastSegment);
		}
	}
	else {
		m_traversedFacesByLine[endFace].emplace_back(contSource, contTarget, lastSegment);
	}
}

void
SingGraphBuilder2DShortestPath::IllegalLineCrossingFinder::registerLineSegmentBetweenTwoFaces(
   const TCellID &f1Id, const TCellID &f2Id, const math::Point &p1, const math::Point &p2, const SourceID &contSource, const TargetID &contTarget)
{
	if (f1Id != f2Id || !m_graphBuilder->m_tool.isAdjacency(f1Id, f2Id)) {
		const Node &n = Tools::getCommonNode(m_graphBuilder->m_mesh->get<Face>(f1Id), m_graphBuilder->m_mesh->get<Face>(f2Id));
		const auto cellId = m_graphBuilder->m_original_faces_number + n.id();
		m_traversedFacesByLine[cellId].emplace_back(contSource, contTarget, math::Segment(p1, p2));
	}
}

std::vector<std::pair<unsigned int, unsigned int>>
SingGraphBuilder2DShortestPath::IllegalLineCrossingFinder::getIllegalOverlappingPaths()
{
	auto timer = Timer("illegal crossing");
	vector<pair<unsigned int, unsigned int>> illegalOverlappingPaths;

	// make pair of potential illegal intersection
	vector<std::tuple<TCellID, const LineSegmentOnCell *, const LineSegmentOnCell *>> candidateFaceByPair;
	for (const auto &candidate : m_traversedFacesByLine) {

		for (int i = 0; i < candidate.second.size() - 1; i++) {
			const LineSegmentOnCell *line_i = &candidate.second[i];
			const SourceID csrc_i = line_i->contSource;
			const TargetID ctgt_i = line_i->contTarget;

			for (int j = i + 1; j < candidate.second.size(); j++) {
				const LineSegmentOnCell *line_j = &candidate.second[j];
				const SourceID csrc_j = line_j->contSource;
				const TargetID ctgt_j = line_j->contTarget;

				if (csrc_i != csrc_j)     // only different src
				{
					if ((csrc_i != ctgt_j && ctgt_i != csrc_j && ctgt_i != ctgt_j)                                    // only different target
					    || (m_graphBuilder->targetIsBoundary(ctgt_i) && m_graphBuilder->targetIsBoundary(ctgt_j))     // or both targets on a boudary
					) {
						candidateFaceByPair.emplace_back(candidate.first, line_i, line_j);
					}
				}
			}
		}
	}
	// check illegal crossing among pair of intersection
	std::set<std::pair<unsigned int, unsigned int>> addedPair;
	for (const auto &lines : candidateFaceByPair) {

		const LineSegmentOnCell *line1 = std::get<1>(lines);
		const LineSegmentOnCell *line2 = std::get<2>(lines);
		const unsigned int line1ID = line1->contSource * m_graphBuilder->m_totalNumberOfVariables + line1->contTarget;
		const unsigned int line2ID = line2->contSource * m_graphBuilder->m_totalNumberOfVariables + line2->contTarget;
		const auto pairId = std::make_pair(line1ID, line2ID);

		if (addedPair.find(pairId) == addedPair.end()) {

			const TCellID cellID = std::get<0>(lines);
			const math::Segment &seg1 = line1->segment;
			const math::Segment &seg2 = line2->segment;
			math::Vector3d crossDir1, crossDir2;
			if (cellID > m_graphBuilder->m_original_faces_number) {
				crossDir1 = (*m_graphBuilder->m_field)[cellID - m_graphBuilder->m_original_faces_number].closestComponentVector(seg1.getDir());
				crossDir2 = (*m_graphBuilder->m_field)[cellID - m_graphBuilder->m_original_faces_number].closestComponentVector(seg2.getDir());
			}
			else {
				crossDir1 = m_graphBuilder->m_triangle_centers_cross[cellID].closestComponentVector(seg1.getDir());
				crossDir2 = m_graphBuilder->m_triangle_centers_cross[cellID].closestComponentVector(seg2.getDir());
			}

			if (math::near(std::fabs(crossDir1.dot(crossDir2)), 1.0)) {
				// check if paths are intersecting
				double intersectionParam;
				math::Point intersectionPnt;
				if (segmentsOverlap(seg1, seg2) || seg1.SecondMetIntersect2D(seg2, intersectionPnt, intersectionParam, m_graphBuilder->m_temp_epsilon)) {
					illegalOverlappingPaths.push_back(pairId);
					addedPair.insert(pairId);
				}
			}
		}
	}
	return illegalOverlappingPaths;
}

/*----------------------------------------------------------------------------------------------------*/
/*----------------------    line intersections functions     -----------------------------------------*/
/*----------------------------------------------------------------------------------------------------*/

void
SingGraphBuilder2DShortestPath::detectLineIntersections()
{
	auto lineIntersectionDetector = LineIntersectionDetector(this);
	lineIntersectionDetector.registerAllSolutionsOnCells();
	lineIntersectionDetector.detectAndCreateIntersections();
}

void
SingGraphBuilder2DShortestPath::LineIntersectionDetector::registerOneSingularityLine(
   SurfaceSingularityLine *singLine, const TCellID &startFace, const TCellID &endFace, const vector<TCellID> &facePath, const bool isBoundary)
{
	if (facePath.empty()) {     // -> fixed variable
		std::vector<TCellID> cellsToRegister {startFace};
		for (const Face &f : m_graphBuilder->m_face2Face_neighbours_by_verts[startFace])
			cellsToRegister.push_back(f.id());
		for (const Face &f : m_graphBuilder->m_face2Face_neighbours_by_verts[endFace])
			cellsToRegister.push_back(f.id());
		cellsToRegister.erase(unique(cellsToRegister.begin(), cellsToRegister.end()), cellsToRegister.end());
		for (const TCellID faceID : cellsToRegister)
			m_traversedCellsByLine[faceID].emplace_back(0, singLine);
		return;
	}

	// faces btw source slot and first face of the path
	if (facePath.front() != startFace) {
		registerLineOnNodesNearSlots(singLine, startFace, facePath.front(), 1, 0, false);
		registerLineOnFacesNearSlots(singLine, 1, 0, false);
	}
	// faces btw face path faces
	for (int i = 0; i < facePath.size() - 1; i++) {
		const TCellID &f1 = facePath[i];
		const TCellID &f2 = facePath[i + 1];
		registerLineBetweenTwoFaces(singLine, f1, f2, i + 2);
		m_traversedCellsByLine[f1].emplace_back(i + 1, singLine);
		if (m_singOrGeomFacesNeighbors[f1]) m_traversedCellsByLine[f1].emplace_back(i + 2, singLine);     // TODO explains this
	}
	const unsigned int slotDiscPointID = facePath.size() + 2;     // slot point index within the singLine discretization
	m_traversedCellsByLine[facePath.back()].emplace_back(slotDiscPointID - 1, singLine);

	// faces btw last path face and target
	if (facePath.back() != endFace) {
		registerLineOnNodesNearSlots(singLine, endFace, facePath.back(), slotDiscPointID, slotDiscPointID + 1, isBoundary);
		registerLineOnFacesNearSlots(singLine, slotDiscPointID, slotDiscPointID + 1, isBoundary);
	}
}

void
SingGraphBuilder2DShortestPath::LineIntersectionDetector::registerLineBetweenTwoFaces(SurfaceSingularityLine *singLine,
                                                                                      const TCellID &f1Id,
                                                                                      const TCellID &f2Id,
                                                                                      const unsigned int &prevDiscId)
{
	if (!m_graphBuilder->m_tool.isAdjacency(f1Id, f2Id)) {
		const Node &n = Tools::getCommonNode(m_graphBuilder->m_mesh->get<Face>(f1Id), m_graphBuilder->m_mesh->get<Face>(f2Id));
		const auto cellId = m_nFaces + n.id();
		m_traversedCellsByLine[cellId].emplace_back(prevDiscId, singLine);
	}
}

void
SingGraphBuilder2DShortestPath::LineIntersectionDetector::registerLineOnFacesNearSlots(SurfaceSingularityLine *singLine,
                                                                                       const unsigned int &slotDiscPointID,
                                                                                       const unsigned int &singDiscPointID,
                                                                                       const bool isBoundary)
{
	if (isBoundary) return;
	const int offset = singDiscPointID == 0 ? 0 : 1;     // when going from source to target
	const auto slot = singDiscPointID == 0 ? singLine->getSlots()[0] : singLine->getSlots()[1];
	if (slot->starting_cell_dim == 1) {     // edge
		for (const TCellID &faceID : m_graphBuilder->m_mesh->get<Edge>(slot->starting_cell_id).getIDs<Face>()) {
			m_traversedCellsByLine[faceID].emplace_back(singDiscPointID - offset, singLine);     // first or last segment
			m_traversedCellsByLine[faceID].emplace_back(slotDiscPointID - offset, singLine);     // 2nd or penultimate segment
		}
	}
	else {     // node (rare case)
		const Node node = m_graphBuilder->m_mesh->get<Node>(slot->starting_cell_id);
		for (const auto faceID : node.getIDs<Face>()) {
			m_traversedCellsByLine[faceID].emplace_back(singDiscPointID - offset, singLine);     // first or last segment
			m_traversedCellsByLine[faceID].emplace_back(slotDiscPointID - offset, singLine);     // 2nd or penultimate segment
		}
	}
}
void
SingGraphBuilder2DShortestPath::LineIntersectionDetector::registerLineOnNodesNearSlots(SurfaceSingularityLine *singLine,
                                                                                       const TCellID &slotFace,
                                                                                       const TCellID &finalPathFace,
                                                                                       const unsigned int &slotDiscPointID,
                                                                                       const unsigned int &singDiscPointID,
                                                                                       const bool isBoundary)
{
	if (isBoundary) {
		registerLineBetweenTwoFaces(singLine, finalPathFace, slotFace, slotDiscPointID - 1);
		return;
	}
	const int offset = singDiscPointID == 0 ? 0 : 1;     // when going from source to target
	if (!m_graphBuilder->m_tool.isAdjacency(slotFace, finalPathFace)) {
		const Node &n = Tools::getCommonNode(m_graphBuilder->m_mesh->get<Face>(slotFace), m_graphBuilder->m_mesh->get<Face>(finalPathFace));
		const auto &faces = n.getIDs<Face>();
		for (const auto f : faces)     // only faces between f1 & f2 are needed in reality
			m_traversedCellsByLine[f].emplace_back(slotDiscPointID - offset, singLine);
		m_traversedCellsByLine[n.id() + m_nFaces].emplace_back(slotDiscPointID - offset, singLine);     // normal registered node at 1
	}
	// register slot nodes at 0 and 1
	const auto slot = singDiscPointID == 0 ? singLine->getSlots()[0] : singLine->getSlots()[1];
	if (slot->starting_cell_dim == 1) {     // edge
		for (const auto &nodeID : m_graphBuilder->m_mesh->get<Edge>(slot->starting_cell_id).getIDs<Node>()) {
			m_traversedCellsByLine[nodeID + m_nFaces].emplace_back(singDiscPointID - offset, singLine);
			m_traversedCellsByLine[nodeID + m_nFaces].emplace_back(slotDiscPointID - offset, singLine);
		}
	}
	else {     // node (rare case)
		m_traversedCellsByLine[slot->starting_cell_id + m_nFaces].emplace_back(slotDiscPointID - offset, singLine);
		m_traversedCellsByLine[slot->starting_cell_id + m_nFaces].emplace_back(slotDiscPointID - offset, singLine);
	}
}

void
SingGraphBuilder2DShortestPath::LineIntersectionDetector::registerAllSolutionsOnCells()
{
	// targets neighbors
	m_singOrGeomFacesNeighbors = std::vector<bool>(m_nFaces, false);
	for (const TCellID &faceId : m_graphBuilder->m_slotFaces) {
		if (faceId < m_nFaces) {
			m_singOrGeomFacesNeighbors[faceId] = true;
			for (const Face &neighborFace : m_graphBuilder->m_face2Face_neighbours_by_verts[faceId]) {
				m_singOrGeomFacesNeighbors[neighborFace.id()] = true;
			}
		}
	}

	for (const auto &solution : m_graphBuilder->m_solutions) {

		const SourceID contSource = solution.first;
		const TargetID contTarget = solution.second;

		auto *surf_line = dynamic_cast<SurfaceSingularityLine *>(m_graphBuilder->m_targets[contSource]->line);
		const TCellID slotFace = m_graphBuilder->m_slotFaces[contSource];
		const auto faceIDsPath = m_graphBuilder->m_finalPaths[contSource][contTarget];

		if (m_graphBuilder->targetIsBoundary(contTarget)) {
			const unsigned int bdryLocalID = contTarget - m_graphBuilder->m_totalNumberOfSlots;
			const auto &boundaryParam = m_graphBuilder->m_bdryPathEndParam[contSource][bdryLocalID];
			registerOneSingularityLine(surf_line, slotFace, boundaryParam.finalFace, faceIDsPath, true);
		}
		else {
			registerOneSingularityLine(surf_line, slotFace, m_graphBuilder->m_slotFaces[contTarget], faceIDsPath, false);
		}
	}
}

void
SingGraphBuilder2DShortestPath::LineIntersectionDetector::detectAndCreateIntersections()
{
	// make pair of line that could intersect each other near a Cell
	vector<std::tuple<TCellID, TraversedCell, TraversedCell>> candidateLinesByPair;
	for (const auto &candidate : m_traversedCellsByLine) {
		const auto cellId = candidate.first;
		const auto lines = candidate.second;
		for (int i = 0; i < lines.size() - 1; i++) {
			const TraversedCell &line_i = lines[i];
			for (int j = i + 1; j < lines.size(); j++) {
				if (line_i.singLine != lines[j].singLine) candidateLinesByPair.emplace_back(cellId, line_i, lines[j]);
			}
		}
	}
	struct intersectionStruct
	{
		SurfaceSingularityLine *singLine;
		unsigned int prevPointID;
		double param;
		SurfaceSingularityPoint *singPoint;
		intersectionStruct(SurfaceSingularityLine *singLine, unsigned int prevPointID, double param, SurfaceSingularityPoint *singPoint) :
		  singLine(singLine), prevPointID(prevPointID), param(param), singPoint(singPoint)
		{
		}
		bool operator>(const intersectionStruct &other) const
		{
			const auto id = singLine->getNumber();
			const auto other_id = other.singLine->getNumber();
			return std::tie(id, prevPointID, param) > std::tie(other_id, other.prevPointID, other.param);
		}
	};

	// assess those potential intersections and sort them
	std::vector<intersectionStruct> sorted_intersections;
	// tricky issue : an intersection can be found near more than one cell, but it needs to be added only once
	std::vector<math::Point> added_points;
	for (const auto &sing : m_graphBuilder->m_graph.getPoints())
		added_points.emplace_back(sing->getLocation());

	for (const auto &lines : candidateLinesByPair) {

		const TraversedCell &line1 = std::get<1>(lines);
		const TraversedCell &line2 = std::get<2>(lines);
		const TCellID cellID = std::get<0>(lines);

		// 2 variables useless for now
		const auto singLineID1 = line1.singLine->getNumber();
		const auto singLineID2 = line2.singLine->getNumber();

		// traversed cell : a node or a face
		const auto &path1 = line1.singLine->getDiscretizationPoints();
		const auto &path2 = line2.singLine->getDiscretizationPoints();
		const auto seg1 = math::Segment(path1[line1.prevPnt], path1[line1.prevPnt + 1]);
		const auto seg2 = math::Segment(path2[line2.prevPnt], path2[line2.prevPnt + 1]);
		double intersectionParam1, intersectionParam2;
		math::Point intersectionPnt;
		double eps = 1e-6;
		if (seg1.intersect3D(seg2, intersectionPnt, intersectionParam2, intersectionParam1, eps)) {
			auto alreadyAdded = [&]() -> bool {
				if (cellID > m_nFaces && segmentsOverlap(seg1, seg2)) return true;     // *
				//* ->not a subtle solution, but it allows to not take into account the first segment of lines belonging to the same singularity
				for (const auto &point : added_points)
					if (math::near(intersectionPnt.distance2(point), 0.0)) return true;
				return false;
			};
			if (alreadyAdded()) continue;

			SurfaceSingularityPoint *new_sing_pnt = m_graphBuilder->m_graph.newSurfacePoint();
			new_sing_pnt->setLocation(intersectionPnt);
			added_points.push_back(intersectionPnt);
			sorted_intersections.emplace_back(line1.singLine, line1.prevPnt, intersectionParam1, new_sing_pnt);
			sorted_intersections.emplace_back(line2.singLine, line2.prevPnt, intersectionParam2, new_sing_pnt);
		}
	}

	// sort intersection points on lines and split them
	std::sort(sorted_intersections.begin(), sorted_intersections.end(), std::greater<>());
	for (const auto &intersection : sorted_intersections) {
		m_graphBuilder->m_graph.splitSurfaceLineSimple(intersection.singPoint, intersection.singLine, intersection.prevPointID);
	}

	writeGraphPoint(m_graphBuilder->m_graph.getPoints(), m_graphBuilder->m_output_directory_name + std::string("gaphPointsAfterLineIntersection.vtk"));
}

}     // namespace gmds
