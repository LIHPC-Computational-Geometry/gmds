/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux (2015)
 *
 * franck.ledoux@cea.fr
 *
 * The FRAME software is a computer program whose purpose is to provide a set
 * of algorithms to build 2D and 3D meshes using frame field concept. The main
 * focus of these algorithms is quadrilateral and hexahedral meshing.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*----------------------------------------------------------------------------*/
/*
 * SingGraphBuilder2DSimultStartRK4.cpp
 *
 *  Created on: 13 juil. 2014
 *      Author: bibi
 */
/*----------------------------------------------------------------------------*/
#include "gmds/math/Constants.h"
#include <gmds/frame/LaplaceCross2D.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator.h>
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
/*----------------------------------------------------------------------------*/
#include <gmds/singGraphBuild/SingGraphBuilder2DSimultStartRK4.h>
#include <gmds/singGraphBuild/SingularityGraphBuilder2D.h>
#include <gmds/singGraphBuild/Tools.h>
/*----------------------------------------------------------------------------*/
#include <chrono>
#include <set>
#include <sstream>
#include <vector>
typedef std::chrono::high_resolution_clock Clock;
/*----------------------------------------------------------------------------*/
using namespace gmds;

SingGraphBuilder2DSimultStartRK4::SingGraphBuilder2DSimultStartRK4(Mesh *AMesh, Variable<math::Cross2D> *AField, const bool ABuildGeomSing) :
  SingularityGraphBuilder2D::SingularityGraphBuilder2D(AMesh, AField, ABuildGeomSing)
{
}

void
SingGraphBuilder2DSimultStartRK4::createSingularityLines()
{
	for (SingularityPoint *pi : m_graph.getPoints()) {
		// WARNING if sing points too close it might "overide the singular triangles"
		initConfusingBalls(pi);
	}
	writeConfusingBalls();
	//========================================================================
	/*We start simultaneously from all singularity points (through all slots) to compute streamlines in a
	 * direction following as close as possible the cross field. At each iteration all such streamlines
	 * advance until they have reached a certain geometric length step (searchStep) - which
	 * is increased with each iteration. If such a streamline reaches the boundary, the respective line is
	 * logges and the slot is being frozen.
	 * If two streamlines (parting from two different singularity points) advance until they are
	 * within a certain distance from each other (thresholdStreamLineDist); in this case they are connected only
	 * if the angle between their directions is not so close to pi/2.
	 */
	//========================================================================

	std::cout << "createSingularityLinesSimultaneousStartRK4" << std::endl;
	/* when deciding whether to connect 2 singularity lines we have two options:
	1. connectByField = false; connect depending on the geometrical distance between the extremities of the singularity
	lines and depending on the angle made by those
	2. connectByField = true; connect depending on the frame field; consider the straight line between the extremities of
	the 2 singularity lines - this line will intersect edges (and vertices rarely); calculate the closest component
	vector at the intersection point with regard to the straight line's direction; if the closest component vector
	coincides (end to end) - we can connect; (by coincides -  consider direction 0-2 and direction 1-3)

	Computing sinsgularity lines: 2 Methods:
	I. Heun's method -> per triangle
	II. RK4 method; 2 approaches:
	A. compute the original method (with the prescribed stepSize)
	   - advantage: fixed stepSize - can be used directly as stopping condition for simultaneous
	   - advantage: could be cheaper computationally (depending on stepSize)
	   - *advantage perhaps: in case of a coarse mesh, perhaps choosing a small stepSize and could improve
	   following the field and making the entire algorithm more robust
	   - advantage perhaps: if the mesh is non-uniform, it will perform better
	   - disadv: extra work for finding traverseTriangle
	   - disadv: if stepSize is not very small -> it doesnt follow the field as close as with B. approach
	B. after RK4 orig we have an outputPoint, which can be outside the current triangle;
	   find intersection of (inputPoint, outputPoint) with the current triangle edges and continue as in Heun
	   - advantage follows the field better
	   - more expensive - proportional to number of traverseTriangle
	*/
	bool visualizeSingBallSim = true;
	bool connectByField = false;

	double stepSize;                    // search step for streamlines
	double thresholdStreamLineDist;     // if 2 lines are within this threshold from one another => connect them
	double weightAngle = 0.75;
	double weightDistance = 0.25;

	stepSize = m_mean_edge_length;
	// WARNING stepSize should be sufficiently small
	// WARNING thresholdStreamLineDist should be sufficiently big
	thresholdStreamLineDist = 100 * m_mean_edge_length * m_ATolerance;

	vector<SingularityPoint::Slot *> searchSlots;
	vector<SingularityPoint *> singularity_points = m_graph.getPoints();

	for (unsigned int i = 0; i < singularity_points.size(); i++) {
		SingularityPoint *pi = singularity_points[i];
		if (pi->getType() == 0) {
			std::vector<SingularityPoint::Slot *> pi_slots = pi->getSlots();
			for (unsigned int j = 0; j < pi_slots.size(); j++) {
				if (!pi_slots[j]->isLaunched) searchSlots.push_back(pi_slots[j]);
			}
		}
	}

	unsigned int NoSimLines = searchSlots.size();
	vector<vector<math::Point>> line_discretization(NoSimLines, vector<math::Point>(0)), copy_line_discretization(NoSimLines, vector<math::Point>(0));
	vector<vector<TCellID>> line_triangles(NoSimLines, vector<TCellID>(0)), copy_line_triangles(NoSimLines, vector<TCellID>(0));
	vector<double> accumulatedDistancePerSlotLine(NoSimLines, 0.0);
	bool find_end_bdry = false;
	bool end_on_free_slot = false;
	vector<TCellID> to_cell_id(NoSimLines);
	vector<int> to_cell_dim(NoSimLines);
	vector<math::Vector3d> to_dir(NoSimLines);
	vector<double> streamlineDeviation(NoSimLines, 0.0);

	unsigned int noFrozenSLots = 0;

	for (unsigned int i = 0; i < NoSimLines; i++) {
		line_discretization[i].push_back(searchSlots[i]->from_point->getLocation());
		line_discretization[i].push_back(searchSlots[i]->location);
		to_cell_id[i] = searchSlots[i]->starting_cell_id;
		to_cell_dim[i] = searchSlots[i]->starting_cell_dim;
		to_dir[i] = searchSlots[i]->direction;     // starting direction
	}

	// max contor stop condition should take into account the m_mean_edge_length relative to the mesh...
	int contor = 0;

	math::Vector3d surfNormal(1.0, 1.0, 1.0);
	vector<pair<pair<unsigned int, unsigned int>, pair<double, double>>>
	   possibleConnectingLines; /*vector storing ((slotA, slotB) , (distance between them, angle of intersection))*/
	vector<pair<unsigned int, unsigned int>> forbiddenConnectingLines;
	double distBtwSlots;
	/*
	rather than adding the streamlines as soon as two meet (or meet the boundary), we will add the pair
	to the possibleConnectingLines and decide later;
	possibleConnectingLines contains also the (geometrical) distance and the angle in between two;
	   for bry streamlines add -1 as second vertex
	   for bry streamlines - if no other streamline has been met in the following two iterations
	   (with good angles and distances = to define what "good" is), construct BdryLine
	   for surface streamlines - if no other streamline has been met in the following two iterations
	   (with better angles and distances), connect; however, stop growLine() for the slot at the
	   iteration where it reaches the vicinity of another slot line*/

	vector<bool> stopIncrease(NoSimLines, false);

	Variable<int> *ball_var_sim = m_mesh->newVariable<int, GMDS_FACE>("sing_ball_sim");

	for (auto f_id : m_mesh->faces()) {
		(*ball_var_sim)[f_id] = 0;
	}

	while (noFrozenSLots < NoSimLines) {
		contor++;
		possibleConnectingLines.clear();
		stopIncrease.clear();
		stopIncrease.resize(NoSimLines, false);
		for (unsigned int i = 0; i < NoSimLines; i++) {
			SingularityPoint *to_sing_pnt = 0;
			SingularityPoint::Slot *to_slot = 0;
			SingularityPoint::Slot *current_slot = searchSlots[i];
			if ((!current_slot->isFreeze) && (!stopIncrease[i])) {
				find_end_bdry = false;
				end_on_free_slot = false;
				bool find_end = false;

				if (contor == 1) {
					// here it can also go out through a node<!!! add triangles
					Edge currentEdge = m_mesh->get<Edge>(to_cell_id[i]);
					std::vector<TCellID> adj_faces = currentEdge.getIDs<Face>();
					bool AOnEdge0, AOnEdge1, AOnEdge2;
					vector<double> lambdas(3);

					if (m_tool.isPntInTri(searchSlots[i]->from_point->getLocation(), m_mesh->get<Face>(adj_faces[0]), AOnEdge0, AOnEdge1, AOnEdge2, lambdas[0],
					                      lambdas[1])) {
						// to_cell_id[i] = adj_faces[1];
						line_triangles[i].push_back(adj_faces[0]);
						line_triangles[i].push_back(adj_faces[1]);
					}
					else {
						// to_cell_id[i] = adj_faces[0];
						line_triangles[i].push_back(adj_faces[1]);
						line_triangles[i].push_back(adj_faces[0]);
					}
				}

				growLineRK4(current_slot, to_slot, line_discretization[i][line_discretization[i].size() - 1], to_dir[i], line_discretization[i], line_triangles[i],
				            to_cell_dim[i], to_cell_id[i], streamlineDeviation[i], stepSize, find_end_bdry, end_on_free_slot, find_end);

				if (to_slot != 0) {
					std::cout << "current_slot->location " << current_slot->location << "; to_slot->location " << to_slot->location << std::endl;
				}
				// first check if find_end_bdry or end_on_free_slot(unlikely)
				if (find_end_bdry) {
					//======================================================================
					// CASE 1 - We finish on the boundary. A geometric point must be created
					// The cell defined by (start_cell_dim[i], start_cell_id[i]) is on the boundary.
					// We have to create a geometric singularity point.
					// however we could do this at the end of the entire iteration(for the current searchStep)
					// line creation
					if (m_withGlobalComments)
						std::cout << "!!!find_end_bdry for current_slot->from_point " << current_slot->from_point->getLocation().X() << " "
						          << current_slot->from_point->getLocation().Y() << " and with dir " << current_slot->direction[0] << " "
						          << current_slot->direction[1] << std::endl;
					SurfaceSingularityLine *surf_line = m_graph.newSurfaceLine();
					int sepNumberTmp = m_graph.getNbLines();
					surf_line->setNumber(sepNumberTmp);
					SingularityPoint *from_sing_pnt = current_slot->from_point;
					surf_line->addSlot(current_slot);
					surf_line->addDiscretizationPoint(from_sing_pnt->getLocation());
					current_slot->line = surf_line;

					math::Vector3d firstDir = math::Vector3d(from_sing_pnt->getLocation(), line_discretization[i][0]);

					current_slot->line_direction = firstDir;
					current_slot->isLaunched = true;
					if (to_slot != 0) {
						streamlineDeviation[i] = (streamlineDeviation[i]
						                          + fabs(1
						                                 - to_slot->direction.dot(math::Vector3d(line_discretization[i][line_discretization[i].size() - 1],
						                                                                         line_discretization[i][line_discretization[i].size() - 2]))));
						streamlineDeviation[i] = streamlineDeviation[i] / (line_discretization[i].size() + 1);
						current_slot->lineDeviation = streamlineDeviation[i];
					}
					// Insertion of line points
					for (unsigned int j = 0; j < line_discretization[i].size(); j++) {
						surf_line->addDiscretizationPoint(line_discretization[i][j]);
					}

					for (unsigned int j = 0; j < line_triangles[i].size(); j++) {
						surf_line->addTraversedFace(line_triangles[i][j]);
					}

					SingularityPoint *geom_pnt;
					SingularityPoint::Slot *incoming_slot;
					createGeometricSingularityPoint(line_discretization[i][line_discretization[i].size() - 1],     // the last point added
					                                to_dir[i],                                                     // the direction we come from
					                                to_cell_dim[i],                                                // the dim. of the cell
					                                to_cell_id[i],                                                 // the id of the cell
					                                geom_pnt,                                                      // the created point
					                                incoming_slot);                                                // and the slot

					current_slot->lineDeviation = streamlineDeviation[i] / line_discretization[i].size();

					surf_line->addSlot(incoming_slot);
					noFrozenSLots++;
					current_slot->isFreeze = true;
					current_slot->isLaunched = true;

					for (unsigned int j = 0; j < possibleConnectingLines.size(); j++) {
						if ((i == possibleConnectingLines[j].first.first) || ((i == possibleConnectingLines[j].first.second))) {
							possibleConnectingLines.erase(possibleConnectingLines.begin() + j);
						}
					}

					gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
					for (unsigned int j = 0; j < line_discretization[i].size() - 1; j++) {
						gmds::Node n1 = m.newNode(line_discretization[i][j].X(), line_discretization[i][j].Y(), line_discretization[i][j].Z());
						gmds::Node n2 = m.newNode(line_discretization[i][j + 1].X(), line_discretization[i][j + 1].Y(), line_discretization[i][j + 1].Z());
						gmds::Face f = m.newTriangle(n1, n1, n2);
					}
					gmds::IGMeshIOService ioService(&m);
					gmds::VTKWriter vtkWriter(&ioService);
					vtkWriter.setCellOptions(gmds::N | gmds::F);
					vtkWriter.setDataOptions(gmds::N | gmds::F);
					std::string file_name = "SimultaneousStartRK4_SingToBdry" + to_string(i) + ".vtk";
					vtkWriter.write(file_name);
					std::cout << "has written SingToBdry" << file_name << std::endl;
				}
				else {     // not find_end_bdry
					if (end_on_free_slot) {
						if (m_withGlobalComments) {
							std::cout << "end_on_free_slot; to_slot " << to_slot << std::endl;
							std::cout << "current_slot->location " << current_slot->location << std::endl;
							std::cout << "line disret " << line_discretization[i][line_discretization[i].size() - 2];
							std::cout << " " << line_discretization[i][line_discretization[i].size() - 1] << std::endl;
						}
						SurfaceSingularityLine *surf_line = m_graph.newSurfaceLine();
						int sepNumberTmp = m_graph.getNbLines();
						surf_line->setNumber(sepNumberTmp);

						SingularityPoint *from_sing_pnt = current_slot->from_point;
						surf_line->addSlot(current_slot);
						surf_line->addDiscretizationPoint(from_sing_pnt->getLocation());

						current_slot->line = surf_line;
						math::Vector3d firstDir = math::Vector3d(from_sing_pnt->getLocation(), line_discretization[i][0]);

						current_slot->line_direction = firstDir;
						current_slot->isLaunched = true;

						current_slot->lineDeviation = streamlineDeviation[i];
						for (unsigned int j = 0; j < line_discretization[i].size(); j++) {
							surf_line->addDiscretizationPoint(line_discretization[i][j]);
						}

						for (unsigned int j = 0; j < line_triangles[i].size(); j++) {
							surf_line->addTraversedFace(line_triangles[i][j]);
						}

						current_slot->isLaunched = true;
						to_slot->isLaunched = true;
						current_slot->isFreeze = true;
						to_slot->isFreeze = true;
						to_slot->line = surf_line;
						to_slot->line_direction = to_dir[i];
						to_slot->lineDeviation = streamlineDeviation[i];

						surf_line->addSlot(to_slot);
						surf_line->addDiscretizationPoint(to_sing_pnt->getLocation());
						to_slot->line_direction = to_dir[i];
						noFrozenSLots = noFrozenSLots + 2;

						gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
						for (unsigned int j = 0; j < line_discretization[i].size() - 1; j++) {
							gmds::Node n1 = m.newNode(line_discretization[i][j].X(), line_discretization[i][j].Y(), line_discretization[i][j].Z());
							gmds::Node n2 = m.newNode(line_discretization[i][j + 1].X(), line_discretization[i][j + 1].Y(), line_discretization[i][j + 1].Z());
							gmds::Face f = m.newTriangle(n1, n1, n2);
						}
						gmds::IGMeshIOService ioService(&m);
						gmds::VTKWriter vtkWriter(&ioService);
						vtkWriter.setCellOptions(gmds::N | gmds::F);
						vtkWriter.setDataOptions(gmds::N | gmds::F);
						std::string file_name = "SimultaneousStartRK4_SingToSingBall" + to_string(i) + "-x.vtk";
						vtkWriter.write(file_name);
						std::cout << "has written SingToSing " << file_name << std::endl;
					}
					else {
						std::cout << "see if we got in the vecinity of another line" << std::endl;
						//======================================================================
						// CASE 2 - We encounter another line within the thresholdStreamLineDist
						// for now, we will simply connect the last points added to the 2 lines (which are within
						// thresholdStreamLineDist); connect only if deviation for the connecting line is sufficiently small;
						// otherwise,let it arrive at boundary
						for (unsigned int j = 0; j < NoSimLines; j++) {
							if ((searchSlots[i]->from_point != searchSlots[j]->from_point) && (line_discretization[j].size() > 2) && (!searchSlots[j]->isFreeze)) {
								distBtwSlots =
								   line_discretization[i][line_discretization[i].size() - 1].distance(line_discretization[j][line_discretization[j].size() - 1]);
								if (distBtwSlots <= thresholdStreamLineDist) {
									bool forbiddenPair = false;
									for (unsigned int ti0 = 0; ti0 < forbiddenConnectingLines.size(); ti0++) {
										if ((forbiddenConnectingLines[ti0].first == i) && (forbiddenConnectingLines[ti0].second == j)) {
											forbiddenPair = true;
											break;
										}
									}
									if (!forbiddenPair) {
										math::Point last_added_pnt = line_discretization[i][line_discretization[i].size() - 1];
										math::Vector3d start_dir = math::Vector3d(last_added_pnt, line_discretization[j][line_discretization[j].size() - 1]);
										if (!connectByField) {
											double tempDev = 0.0;
											vector<gmds::Face> candidate_faces;
											vector<bool> visitedFaces(m_original_faces_number, false), addedFaces(m_original_faces_number, false);
											if (to_cell_dim[i] == 0) {     // node
												gmds::Node currentNode = m_mesh->get<gmds::Node>(to_cell_id[i]);
												candidate_faces = currentNode.get<gmds::Face>();
												for (unsigned int t = 0; t < candidate_faces.size(); t++) {
													visitedFaces[candidate_faces[t].id()] = true;
												}
											}
											else {
												if (to_cell_dim[i] == 1) {     // edge
													gmds::Edge currentEdge = m_mesh->get<gmds::Edge>(to_cell_id[i]);
													candidate_faces = currentEdge.get<gmds::Face>();
													for (unsigned int t = 0; t < candidate_faces.size(); t++) {
														visitedFaces[candidate_faces[t].id()] = true;
													}
												}
												else {     // to_cell_dim[i]==2 - face
													vector<gmds::Node> currentNodes = m_mesh->get<gmds::Face>(to_cell_id[i]).get<Node>();
													for (unsigned int tt = 0; tt < 3; tt++) {
														vector<gmds::Face> temp_faces = currentNodes[tt].get<gmds::Face>();
														for (unsigned int tt2 = 0; tt2 < temp_faces.size(); tt2++) {
															if (!visitedFaces[temp_faces[tt2].id()]) {
																visitedFaces[temp_faces[tt2].id()] = true;
																candidate_faces.push_back(temp_faces[tt2]);
															}
														}
													}
												}
											}

											for (unsigned int t = 0; t < line_triangles[i].size(); t++)
												addedFaces[line_triangles[i][t]] = true;
											for (unsigned int t = 0; t < line_triangles[j].size(); t++)
												addedFaces[line_triangles[j][t]] = true;

											math::Segment seg1(last_added_pnt, line_discretization[j][line_discretization[j].size() - 1]);
											math::Ray from_ray(last_added_pnt, line_discretization[j][line_discretization[j].size() - 1]);
											math::Vector3d connLineDir = from_ray.getDirUnit();
											math::Vector3d closest2Current, closest0, closest1;
											unsigned int t = 0;
											bool stopAdding = false;
											vector<bool> visitedEdges(m_mesh->getNbEdges(), false);
											unsigned int noIntEdges = 0;

											copy_line_triangles[i] = line_triangles[i];

											while (t < candidate_faces.size()) {
												Face currentFace = candidate_faces[t];
												vector<gmds::Edge> currentEdges = currentFace.get<gmds::Edge>();
												for (unsigned int tt = 0; tt < 3; tt++) {
													vector<gmds::Node> currentNodes = currentEdges[tt].get<gmds::Node>();
													math::Segment oppSeg(currentNodes[0].point(),
                                                                         currentNodes[1].point());
													math::Point intersectionPnt;
													double intersectionParam;
													if (seg1.SecondMetIntersect2D(oppSeg, intersectionPnt, intersectionParam, m_temp_epsilon)) {
														if (from_ray.SecondMetIntersect2D(oppSeg, intersectionPnt, intersectionParam, m_temp_epsilon)) {
															if (!visitedEdges[currentEdges[tt].id()]) {
																math::Cross2D cross_0 = (*m_field)[currentNodes[0].id()];
																math::Cross2D cross_1 = (*m_field)[currentNodes[1].id()];
																closest0 = cross_0.closestComponentVector(connLineDir);
																closest1 = cross_1.closestComponentVector(connLineDir);
																closest2Current = intersectionParam * closest0 + (1 - intersectionParam) * closest1;
																tempDev = tempDev + fabs(1.0 - connLineDir.dot(closest2Current));
																visitedEdges[currentEdges[tt].id()] = true;
																noIntEdges++;
															}
															if (currentFace.id() == line_triangles[j][line_triangles[j].size() - 1]) stopAdding = true;
															if (!addedFaces[currentFace.id()]) {
																copy_line_triangles[i].push_back(currentFace.id());
																addedFaces[currentFace.id()] = true;
															}
															if (!stopAdding) {
																vector<Face> adj_faces0 = currentNodes[0].get<Face>();
																vector<Face> adj_faces1 = currentNodes[1].get<Face>();
																adj_faces0.insert(adj_faces0.end(), adj_faces1.begin(), adj_faces1.end());
																for (unsigned int ttt = 0; ttt < adj_faces0.size(); ttt++) {
																	if (!visitedFaces[adj_faces0[ttt].id()]) {
																		visitedFaces[adj_faces0[ttt].id()] = true;
																		candidate_faces.push_back(adj_faces0[ttt]);
																	}
																}
															}
														}
													}
												}
												t++;
											}

											std::cout << "slots " << i << " and " << j << " to be connected with deviation " << fmod(tempDev, 1.5707963267948966)
											          << std::endl;
											math::Vector3d test1(line_discretization[i][line_discretization[i].size() - 2],
											                     line_discretization[i][line_discretization[i].size() - 1]);
											math::Vector3d test2(line_discretization[j][line_discretization[j].size() - 2],
											                     line_discretization[j][line_discretization[j].size() - 1]);
											test1.normalize();
											test2.normalize();
											// if close to perpendicular=> shouldn't connect
											tempDev = test1.angle(test2);     // tempDevin [0,PI]]
											// we must also take into account the distance; 2 streamlines could be parallel but far
											// from one another(although inside the proximityRadius)
											math::Vector3d tempCross = test1.cross(test2);
											// we connect only if the dot product with the normal (in 2d case - (1,1,1) is higher
											// that 0) - in order to account for orientation
											if (tempCross.dot(surfNormal) >= 0) {
												tempDev = tempDev / (double) (math::Constants::PI);

												if ((tempDev >= 0.75) && (tempDev <= 1.25)) {
													possibleConnectingLines.push_back(make_pair(make_pair(i, j), make_pair(distBtwSlots, std::fabs(1.0 - tempDev))));
													// ideally the angle btw the lines should be 0; therefore tempDev = 1
													stopIncrease[j] = true;
												}
											}
										}
										else {
											// connectByField hasnt been implemented
										}
									}
								}
							}
						}
					}
				}

				gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
				for (unsigned int ti = 0; ti < NoSimLines; ti++) {
					for (unsigned int ti1 = 0; ti1 < line_discretization[ti].size() - 1; ti1++) {
						gmds::Node n1 = m.newNode(line_discretization[ti][ti1].X(), line_discretization[ti][ti1].Y(), line_discretization[ti][ti1].Z());
						gmds::Node n2 = m.newNode(line_discretization[ti][ti1 + 1].X(), line_discretization[ti][ti1 + 1].Y(), line_discretization[ti][ti1 + 1].Z());
						gmds::Face f = m.newTriangle(n1, n1, n2);
					}
				}

				gmds::IGMeshIOService ioService(&m);
				gmds::VTKWriter vtkWriter(&ioService);
				vtkWriter.setCellOptions(gmds::N | gmds::F);
				vtkWriter.setDataOptions(gmds::N | gmds::F);
				std::string test_file_name = m_output_directory_name + "SimultaneousStartRK4_-lastStep" + ".vtk";
				vtkWriter.write(test_file_name);
			}
		}

		// if one slot is within radius with 2 other slots, choose the one with the lowest (geometric) deviation

		if (possibleConnectingLines.size() > 0) {
			unsigned int t1 = 0;
			unsigned int t2;
			double connectionTerm1, connectionTerm2;

			while (t1 < possibleConnectingLines.size() - 1) {
				t2 = t1 + 1;
				while (t2 < possibleConnectingLines.size()) {
					if (((possibleConnectingLines[t1].first.first == possibleConnectingLines[t2].first.first)
					     && (possibleConnectingLines[t1].first.second == possibleConnectingLines[t2].first.second))
					    || ((possibleConnectingLines[t1].first.first == possibleConnectingLines[t2].first.second)
					        && (possibleConnectingLines[t1].first.second == possibleConnectingLines[t2].first.first))) {

						std::cout << "connection between " << possibleConnectingLines[t1].first.first << " and " << possibleConnectingLines[t1].first.second
						          << std::endl;
						connectionTerm1 =
						   weightDistance
						   * (possibleConnectingLines[t1].second.first / (possibleConnectingLines[t1].second.first + possibleConnectingLines[t2].second.first));
						connectionTerm1 = connectionTerm1
						                  + weightAngle
						                       * (possibleConnectingLines[t1].second.second
						                          / (possibleConnectingLines[t1].second.second + possibleConnectingLines[t2].second.second));
						connectionTerm2 =
						   weightDistance
						   * (possibleConnectingLines[t2].second.first / (possibleConnectingLines[t1].second.first + possibleConnectingLines[t2].second.first));
						connectionTerm2 = connectionTerm2
						                  + weightAngle
						                       * (possibleConnectingLines[t2].second.second
						                          / (possibleConnectingLines[t1].second.second + possibleConnectingLines[t2].second.second));

						if (connectionTerm1 > connectionTerm2) {     // remove t1
							possibleConnectingLines.erase(possibleConnectingLines.begin() + t1);
							t1--;
							break;
						}
						else {
							possibleConnectingLines.erase(possibleConnectingLines.begin() + t2);
							t2--;
						}
					}
					t2++;
				}
				t1++;
			}
		}

		for (unsigned int t2 = 0; t2 < possibleConnectingLines.size(); t2++) {
			unsigned int i = possibleConnectingLines[t2].first.first;
			unsigned int j = possibleConnectingLines[t2].first.second;
			if (visualizeSingBallSim) {
				gmds::math::Point centerPoint = line_discretization[j].back() + line_discretization[i].back();
				centerPoint.X() = centerPoint.X() / 2;
				centerPoint.Y() = centerPoint.Y() / 2;
				for (unsigned int t5 = 0; t5 < m_original_faces_number - 1; t5++) {
					gmds::Face testFace1 = m_mesh->get<Face>(t5);
					if (m_triangle_centers[t5].distance(centerPoint) < thresholdStreamLineDist) {
						for (unsigned int t6 = t5 + 1; t6 < m_original_faces_number; t6++) {
							gmds::Face testFace2 = m_mesh->get<Face>(t6);
							if (m_triangle_centers[t6].distance(centerPoint) < thresholdStreamLineDist) {
								if (m_triangle_centers[t5].distance(m_triangle_centers[t6]) < thresholdStreamLineDist) {
									(*ball_var_sim)[t5] = 1;
									(*ball_var_sim)[t6] = 1;
								}
							}
						}
					}
				}
			}

			// the 2 lines are within radius, they must be connected
			for (int t = line_discretization[j].size() - 1; t >= 0; t--) {
				line_discretization[i].push_back(line_discretization[j][t]);
			}
			if (m_withGlobalComments) {
				std::cout << "creates line between " << i << " and " << j << " with angle 1+ " << possibleConnectingLines[t2].second.second << std::endl;
				std::cout << "and distance " << possibleConnectingLines[t2].second.first << std::endl;
				std::cout << "these should be diff(if small triangle count, not necessarilly) " << line_triangles[i].size() << "<" << copy_line_triangles[i].size()
				          << std::endl;
			}
			line_triangles[i] = copy_line_triangles[i];
			for (int t = line_triangles[j].size() - 1; t >= 0; t--) {
				line_triangles[i].push_back(line_triangles[j][t]);
			}

			SurfaceSingularityLine *surf_line = m_graph.newSurfaceLine();
			int sepNumberTmp = m_graph.getNbLines();
			surf_line->setNumber(sepNumberTmp);
			// connect line to initial singularity point
			SingularityPoint *from_sing_pnt = searchSlots[i]->from_point;
			SingularityPoint::Slot *the_other_slot = searchSlots[j];
			SingularityPoint *towards_sing_pnt = the_other_slot->from_point;
			surf_line->addSlot(searchSlots[i]);
			surf_line->addSlot(the_other_slot);
			// surf_line->addDiscretizationPoint(from_sing_pnt->getLocation());
			math::Vector3d firstDir = math::Vector3d(line_discretization[i][0], line_discretization[i][1]);
			searchSlots[i]->line_direction = firstDir;

			firstDir = math::Vector3d(line_discretization[j][0], line_discretization[j][1]);
			searchSlots[j]->line_direction = firstDir;

			streamlineDeviation[i] = streamlineDeviation[i] + streamlineDeviation[j];
			streamlineDeviation[i] = streamlineDeviation[i] / (line_discretization[i].size() + 1);
			streamlineDeviation[j] = streamlineDeviation[i];
			// WARNING also compute streamlineDeviation locally, inside thresholdStreamLineDist
			searchSlots[i]->lineDeviation = streamlineDeviation[i];
			searchSlots[j]->lineDeviation = streamlineDeviation[i];

			for (unsigned int t = 0; t < line_discretization[i].size(); t++) {
				surf_line->addDiscretizationPoint(line_discretization[i][t]);
			}

			for (unsigned int t = 0; t < line_triangles[i].size(); t++) {
				surf_line->addTraversedFace(line_triangles[i][t]);
			}
			searchSlots[i]->line = surf_line;
			searchSlots[j]->line = surf_line;     //  the inverse of this; although in the original algorithm it's also surf_line;

			searchSlots[i]->isLaunched = true;
			searchSlots[i]->isFreeze = true;
			searchSlots[j]->isLaunched = true;
			searchSlots[j]->isFreeze = true;

			/*gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
			 *for (unsigned int j = 0; j < line_discretization[i].size()-1; j++){
			      gmds::Node n1 = m.newNode(line_discretization[i][j].X(), line_discretization[i][j].Y(),
			line_discretization[i][j].Z()); gmds::Node n2 = m.newNode(line_discretization[i][j + 1].X(),
			line_discretization[i][j + 1].Y(), line_discretization[i][j + 1].Z()); gmds::Face f  =  m.newTriangle(n1, n1,
			n2);
			}
			gmds::IGMeshIOService ioService(&m);
			gmds::VTKWriter vtkWriter(&ioService);
			vtkWriter.setCellOptions(gmds::N|gmds::F);
			vtkWriter.setDataOptions(gmds::N|gmds::F);
			std::string file_name = "SingToSing"+to_string(i)+"-"+to_string(j)+".vtk";
			vtkWriter.write(file_name);
			std::cout<<"has written SingToSing"<<file_name<<std::endl;*/

			noFrozenSLots = noFrozenSLots + 2;
		}

		/*only for visualization purposes
		   for(unsigned int j0 = 0; j0 < line_discretization.size(); j0++){
		   gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
		   std::cout<<"line_discretization["<<j0<<"].size() "<<line_discretization[j0].size()<<std::endl;
		   for (unsigned int j = 0; j < line_discretization[j0].size()-1; j++){
		      gmds::Node n1 = m.newNode(line_discretization[j0][j].X(), line_discretization[j0][j].Y(),
		line_discretization[j0][j].Z()); gmds::Node n2 = m.newNode(line_discretization[j0][j + 1].X(),
		line_discretization[j0][j + 1].Y(), line_discretization[j0][j + 1].Z()); gmds::Face f  =  m.newTriangle(n1, n1,
		n2);
		   }

		   gmds::IGMeshIOService ioService(&m);
		   gmds::VTKWriter vtkWriter(&ioService);
		   vtkWriter.setCellOptions(gmds::N|gmds::F);
		   vtkWriter.setDataOptions(gmds::N|gmds::F);
		   std::string file_name = "SimultaneuosRK4_lastStep"+to_string(j0)+".vtk";
		   vtkWriter.write(file_name);
		}    */
		writeOutputSingle("boundary_line");

		if (contor == m_mesh->getNbEdges()) {
			std::cout << "m_mesh->getNbEdges() " << m_mesh->getNbEdges() << std::endl;
			throw GMDSException("contor==total number of edges");
		}
	}

	/*for the sp-sp streamlines, remesh the area around them (all neighbouring faces of the
	 *traversed triangles), recompute cross field in order to have a final refined streamline*/
	/*vector<bool> visitedFaces(m_original_faces_number, false);

	vector<SurfaceSingularityLine*> surface_lines = m_graph.getSurfaceLines();

	for(unsigned int i=0; i<surface_lines.size(); i++){
	   vector<SingularityPoint*> endPoints = surface_lines[i]->getEndPoints();
	   if((endPoints[0]->getType()==0)&&(endPoints[1]->getType()==0)){
	      //streamline in between two singularity points
	      vector<gmds::TCellID> trianglesToRemesh = surface_lines[i]->getTraversedFaces();
	      unsigned int originalTraversedTrianglesNo = trianglesToRemesh.size();
	      for(unsigned int j=0; j<originalTraversedTrianglesNo; j++){
	         visitedFaces[trianglesToRemesh[j]] = true;
	      }
	      for(unsigned int j=0; j<originalTraversedTrianglesNo; j++){
	         vector<gmds::Node> currentNodes = m_mesh->get<Face>(trianglesToRemesh[j]).get<Node>();
	         for(unsigned int t=0; t< 3; t++){
	            vector<Face> adj_faces = currentNodes[t].get<Face>();
	            for(unsigned int t1=0; t1<adj_faces.size(); t1++){
	               if(!visitedFaces[adj_faces[t1].id()]){
	                  trianglesToRemesh.push_back(adj_faces[t1].id());
	                  visitedFaces[adj_faces[t1].id()] = true;
	               }
	            }
	         }
	      }
	      gmds::Mesh newLocalMesh(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
	      gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

	      vector<TCellID> newLocalMesh_id_to_mesh_id_node;
	      gmds::Variable<gmds::math::Cross2D>* local_cross_field_2D =
	newLocalMesh.newVariable<math::Cross2D,GMDS_NODE>("cross_X");
	      //at the begining all seem to be intialized as OX axis

	      remeshTriangles(&newLocalMesh,
	                     newLocalMesh_id_to_mesh_id_node
	                     local_cross_field_2D,
	                     trianglesToRemesh);

	      //visualize newly computed cross field for the last sp-sp streamline
	      Variable<math::Vector3d>* crossVect0 = newLocalMesh.newVariable<math::Vector3d,GMDS_NODE>("crossVect0");
	      Variable<math::Vector3d>* crossVect1 = newLocalMesh.newVariable<math::Vector3d,GMDS_NODE>("crossVect1");
	      Variable<math::Vector3d>* crossVect2 = newLocalMesh.newVariable<math::Vector3d,GMDS_NODE>("crossVect2");
	      Variable<math::Vector3d>* crossVect3 = newLocalMesh.newVariable<math::Vector3d,GMDS_NODE>("crossVect3");
	      for(auto f_id:newLocalMesh.faces()){
	         Face current = newLocalMesh.get<Face>(f_id);
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

	      gmds::IGMeshIOService meshIoServref(&newLocalMesh);
	      gmds::VTKWriter writerref(&meshIoServref);
	      writerref.setCellOptions(gmds::N|gmds::F);
	      writerref.setDataOptions(gmds::N|gmds::F);
	         std::stringstream file_nameref;
	      file_nameref <<"final-newLocalMesh-crossVectors.vtk";
	         writerref.write(file_nameref.str());
	      visitedFaces.clear();
	         visitedFaces.resize(m_original_faces_number, false);
	   }
	   }
	*/

	/*WARNING Although improbable, for the lines that have arrived within the confusing ball of a
	 * different singularity, we haven't done anything; at the end, inspect if we have such cases and
	 * if so, treat them*/

	if (visualizeSingBallSim) {

		gmds::IGMeshIOService meshIoServ(m_mesh);
		gmds::VTKWriter writerB(&meshIoServ);
		writerB.setCellOptions(gmds::N | gmds::F);
		writerB.setDataOptions(gmds::N | gmds::F);
		std::stringstream file_name2;
		file_name2 << m_output_directory_name << "-confusing_balls_sim.vtk";
		writerB.write(file_name2.str());
		std::cout << "wrote confusing_balls_sim.vtk" << std::endl;
	}
}
