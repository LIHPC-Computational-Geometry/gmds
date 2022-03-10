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
 * SingGraphBuilder2DOriginal.cpp
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
#include <gmds/singGraphBuild/SingGraphBuilder2DOriginal.h>
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

SingGraphBuilder2DOriginal::SingGraphBuilder2DOriginal(Mesh *AMesh, Variable<math::Cross2D> *AField, const bool ABuildGeomSing) :
  SingularityGraphBuilder2D::SingularityGraphBuilder2D(AMesh, AField, ABuildGeomSing)
{
}

void
SingGraphBuilder2DOriginal::createSingularityLines()
{
	for (SingularityPoint *pi : m_graph.getPoints()) {
		// WARNING if sing points too close it might "overide the singular triangles"
		initConfusingBalls(pi);
	}
	writeConfusingBalls();

	// Now we need to create separatrices on the faces
	// we will start by launching them from the slot of the singularities
	// for each one we have the pos and dir to launch the sep into, as well as the first triangle

	//========================================================================
	// At the beginning, only the singularity points of the cross field are
	// known. Some other points may be created, when we intersect
	// boundary, singularity lines, etc.
	//========================================================================
	std::vector<SingularityPoint *> singularity_points = m_graph.getPoints();

	for (unsigned int i = 0; i < singularity_points.size(); i++) {
		SingularityPoint *pi = singularity_points[i];
		std::vector<SingularityPoint::Slot *> pi_slots = pi->getSlots();

		for (unsigned int j = 0; j < pi_slots.size(); j++) {
			if (!pi_slots[j]->isLaunched) m_free_slots.push_back(pi_slots[j]);
		}
	}
	//========================================================================
	// Creation of singularity lines from cross field singular points
	//========================================================================

	unsigned int cont = 0;
	while (!m_free_slots.empty()) {
		SingularityPoint::Slot *current_slot = m_free_slots.front();
		m_free_slots.pop_front();
		if (m_withGlobalComments) {
			std::cout << "current_slot->location " << current_slot->location << std::endl;
			if (current_slot->line) std::cout << "begining current_slot->line has line " << current_slot->line->getNumber() << std::endl;
		}
		if (current_slot->isLaunched != true) {
			current_slot->isLaunched = true;
			if (m_withGlobalComments) std::cout << "computeSingularityLine" << std::endl;
			SingularityPoint::Slot *removed_slot = 0;
			computeSingularityLine(current_slot->from_point, current_slot, cont, removed_slot);

			if (removed_slot != 0) {
				m_free_slots.push_back(removed_slot);
				removed_slot->isLaunched = false;
				if (m_withGlobalComments) {
					std::cout << "connection to slot has been removed; removed_slot->location " << removed_slot->location << std::endl;
					if (removed_slot->line)
						std::cout << "connection to slot has been removed removed_slot->line has line " << removed_slot->line->getNumber() << std::endl;
				}
			}
			cont++;
			writeOutput("boundary_line");
		}
		if (m_withGlobalComments) {
			if (current_slot->line) std::cout << "end current_slot->line has line " << current_slot->line->getNumber() << std::endl;

			std::cout << "each step m_free_slots detectLineIntersections sing points " << std::endl;
			vector<SingularityPoint *> pi = m_graph.getPoints();
			for (unsigned int i = 0; i < pi.size(); i++)
				std::cout << pi[i]->getLocation().X() << " " << pi[i]->getLocation().Y() << std::endl;
		}
	}
}

void
SingGraphBuilder2DOriginal::computeSingularityLine(SingularityPoint *AFromPoint,
                                                   SingularityPoint::Slot *AFromSlot,
                                                   unsigned int &cont,
                                                   SingularityPoint::Slot *&ARemovedSlot)
{
	//========================================================================
	// Data initialization for line building
	//========================================================================
	bool end_on_bnd = false;
	bool end_on_free_slot = false;
	bool must_create_pnt = false;
	SingularityPoint *to_sing_pnt = 0;
	SingularityPoint::Slot *to_slot = 0;
	TCellID to_cell_id;
	int to_cell_dim;
	math::Point to_pnt;
	math::Vector3d to_dir;
	std::vector<math::Point> line_discretization;
	std::vector<TCellID> line_triangles;
	double streamlineDeviation = 0.0;
	//========================================================================
	// Stream line computation
	//========================================================================
	computeStreamLine(AFromPoint, AFromSlot, to_sing_pnt, to_slot, to_pnt, to_dir, line_discretization, line_triangles, to_cell_dim, to_cell_id,
	                  streamlineDeviation, end_on_bnd, end_on_free_slot, must_create_pnt);
	/*gmds::Mesh meshSing(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));

	for(unsigned int t=0; t<line_triangles.size(); t++) {
	   vector<gmds::Node> currentTriNodes = m_mesh->get<Face>(line_triangles[t]).get<Node>();
	   std::cout<<"currentTriNodes,isze() "<<currentTriNodes.size()<<std::endl;
	   gmds::Node mySing1 = meshSing.newNode(currentTriNodes[0].getPoint());
	   gmds::Node mySing2 = meshSing.newNode(currentTriNodes[1].getPoint());
	   gmds::Node mySing3 = meshSing.newNode(currentTriNodes[2].point());
	   meshSing.newTriangle(mySing1, mySing2, mySing3);
	}
	gmds::IGMeshIOService ioServiceSing(&meshSing);
	gmds::VTKWriter vtkWriterSing(&ioServiceSing);
	vtkWriterSing.setCellOptions(gmds::N|gmds::F);
	vtkWriterSing.setDataOptions(gmds::N|gmds::F);
	std::string file_name = "HolesInSquare-allStreamLines_"+std::to_string(cont)+".vtk";
	stkWriterSing.write(file_name);  */

	//========================================================================
	// Stream line analysis and creation
	//========================================================================
	if (m_withGlobalComments) {
		std::cout << "line_discretization[0] " << line_discretization[0].X() << " " << line_discretization[0].Y() << std::endl;
		std::cout << "line_discretization[1] " << line_discretization[1].X() << " " << line_discretization[1].Y() << std::endl;
	}
	// line creation
	SurfaceSingularityLine *surf_line = m_graph.newSurfaceLine();
	int sepNumberTmp = m_graph.getNbLines();
	surf_line->setNumber(sepNumberTmp);

	// connect line to initial singularity point
	SingularityPoint *from_sing_pnt = AFromSlot->from_point;
	surf_line->addSlot(AFromSlot);
	surf_line->addDiscretizationPoint(from_sing_pnt->getLocation());
	if (m_withGlobalComments) std::cout << "surf_line also adds disc pt " << from_sing_pnt->getLocation().X() << from_sing_pnt->getLocation().Y() << std::endl;

	AFromSlot->line = surf_line;
	// AFromSlot->line_direction =  AFromSlot->direction;
	math::Vector3d firstDir = math::Vector3d(from_sing_pnt->getLocation(), line_discretization[0]);
	if (m_withGlobalComments) std::cout << "firstDir " << firstDir.X() << " " << firstDir.Y() << std::endl;
	AFromSlot->line_direction = firstDir;
	AFromSlot->isLaunched = true;
	if (to_slot != 0) {
		streamlineDeviation = (streamlineDeviation
		                       + fabs(1
		                              - to_slot->direction.dot(math::Vector3d(line_discretization[line_discretization.size() - 1],
		                                                                      line_discretization[line_discretization.size() - 2]))));
		streamlineDeviation = streamlineDeviation / (line_discretization.size() + 1);
		AFromSlot->lineDeviation = streamlineDeviation;
	}
	// Insertion of line points
	for (unsigned int i = 0; i < line_discretization.size(); i++) {
		surf_line->addDiscretizationPoint(line_discretization[i]);
	}

	for (unsigned int i = 0; i < line_triangles.size(); i++) {
		surf_line->addTraversedFace(line_triangles[i]);
	}

	//========================================================================
	// Termination of the streamline
	//========================================================================
	if (end_on_bnd) {
		if (m_withGlobalComments) std::cout << "!!!!!!!!!!!!!!!!!end_on_bnd " << std::endl;
		//======================================================================
		// CASE 1 - We finish on the boundary. A geometric point must be created
		// The cell defined by (start_cell_dim, start_cell_id) is on the boundary.
		// We have to create a geometric singularity point so.
		// Here it always assumes we reach a non nodeOnPoint vertex; otherwise just error

		SingularityPoint *geom_pnt;
		SingularityPoint::Slot *incoming_slot;
		createGeometricSingularityPoint(to_pnt,             // the point we are
		                                to_dir,             // the direction we come from
		                                to_cell_dim,        // the dim. of the cell
		                                to_cell_id,         // the id of the cell
		                                geom_pnt,           // the created point
		                                incoming_slot);     // and the slot

		AFromSlot->lineDeviation = streamlineDeviation / line_discretization.size();
		if (m_withGlobalComments) std::cout << "AFromSlot->lineDeviation " << AFromSlot->lineDeviation << std::endl;

		surf_line->addSlot(incoming_slot);
		// surf_line->addDiscretizationPoint(geom_pnt->getLocation());/* This has already been added, but check what
		// happens if end on node on point*/
		if (m_withGlobalComments)
			std::cout << "from " << AFromPoint->getMesh<Face>()[0] << " towards bdry point " << to_pnt.X() << " " << to_pnt.Y() << std::endl;
		//    geom_pnt->connectLine(surf_line, to_dir);
	}     // if(to_cell_dim==0)
	else {
		if (m_withGlobalComments) std::cout << "!!!!!!!!!!!!!!!!! not end_on_bnd " << std::endl;

		if (end_on_free_slot) {
			//======================================================================
			// CASE 2 - We finish on a field singularity where a slot is free
			// I don't need to backtrack, i already have the line
			// TODO connect it here; the same situation for an occupied slot for which the new line is better aligned
			if (m_withGlobalComments) std::cout << "end_on_free_slot, to_slot->isLaunched " << to_slot->isLaunched << std::endl;
		}
		// else{
		//============================================
		if (!end_on_free_slot) {
			//======================================================================
			// CASE 3 - We finish on a field singularity where a slot is not free
			// The slot must be freed before connection.
			if (m_withGlobalComments) std::cout << "!end_on_free_slot" << std::endl;
			SingularityPoint *singPointToRemove = 0;
			if (to_slot != 0) {
				SingularityLine *removedLine = to_slot->line;
				SingularityPoint *otherSingPnt;

				vector<SingularityPoint *> removedLineSingPoints = removedLine->getEndPoints();
				if (removedLineSingPoints[0] == to_sing_pnt) {
					otherSingPnt = removedLineSingPoints[1];
				}
				else
					otherSingPnt = removedLineSingPoints[0];

				if (otherSingPnt->getType() == 0) {
					std::vector<SingularityPoint::Slot *> otherSingPnt_slots = otherSingPnt->getSlots();

					for (unsigned int t = 0; t < otherSingPnt_slots.size(); t++) {
						if (otherSingPnt_slots[t]->line->getNumber() == to_slot->line->getNumber()) ARemovedSlot = otherSingPnt_slots[t];
					}
					if (m_withGlobalComments) {
						std::cout << "to_slot!=0; to_slot->line " << to_slot->line->getNumber() << " will be removed " << std::endl;
						std::cout << "this singline is between [" << to_slot->line->getEndPoints()[0]->getLocation().X() << ","
						          << to_slot->line->getEndPoints()[0]->getLocation().Y() << "]" << std::endl;
						std::cout << "and [" << to_slot->line->getEndPoints()[1]->getLocation().X() << "," << to_slot->line->getEndPoints()[1]->getLocation().Y()
						          << "]" << std::endl;
						std::cout << "to_slot->starting_cell_id " << to_slot->starting_cell_id << std::endl;
					}
				}
				else { /* the other sing point is a geometric one;
					 if it has been created as streamline-bdry intersection => it should be removed */
					vector<SingularityPoint::Slot *> curvePointSlots = otherSingPnt->getSlots();
					for (unsigned int t = 0; t < curvePointSlots.size(); t++) {
						if (!curvePointSlots[t]->isFreeze) {
							singPointToRemove = otherSingPnt;
							break;
						}
					}
				}
			}
			// vector<SingularityLine*> linesAdded = m_graph.getCurveLines();
			removeSingularityLine(to_slot->line);
			m_graph.removeLine(to_slot->line);

			if (singPointToRemove != 0) {
				m_graph.removePoint(singPointToRemove);
			}
		}

		to_slot->isLaunched = true;
		to_slot->line = surf_line;
		to_slot->line_direction = to_dir;
		to_slot->lineDeviation = streamlineDeviation;

		/* We need to compute the discretization points, which are in the
		triangles of the confusing ball.*/
		surf_line->addSlot(to_slot);
		surf_line->addDiscretizationPoint(to_sing_pnt->getLocation());

		/* vector<math::Point> Newline_discretization = surf_line->getDiscretizationPoints();
		gmds::Mesh meshSing(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
		std::cout<<"first elem Newline_discretization "<<std::endl;
		for(unsigned int j=0; j<Newline_discretization.size(); j++) {
		   gmds::Node mySing = meshSing.newNode(Newline_discretization[j].X(), Newline_discretization[j].Y(),
		Newline_discretization[j].Z()); meshSing.newTriangle(mySing, mySing, mySing);
		}
		gmds::IGMeshIOService ioServiceSing(&meshSing);
		gmds::VTKWriter vtkWriterSing(&ioServiceSing);
		vtkWriterSing.setCellOptions(gmds::N|gmds::F);
		vtkWriterSing.setDataOptions(gmds::N|gmds::F);
		std::string file_name =
		"HolesInSquare-discretizationPointsBeforeBacktrack_"+std::to_string(m_graph.getNbLines())+".vtk";
		vtkWriterSing.write(file_name);  */

		streamlineDeviation = 0.0;
		bool foundPath = false;
		if (m_withGlobalComments)
			std::cout << "backtrackSingularityLine from " << to_sing_pnt->getMesh<Face>()[0] << ", backtrackSingularityLine towards "
			          << AFromPoint->getMesh<Face>()[0] << std::endl;
		backtrackSingularityLine(surf_line,       // the line we modify
		                         to_sing_pnt,     // the point we start from
		                         to_slot,         // the slot we start from
		                         AFromPoint,      // the point we go to
		                         AFromSlot,       // the slot we go to
		                         to_dir,          // the direction of the streamline for to_slot
		                         foundPath);      // boolean value indicating if we have found a path

		if (foundPath) {
			to_slot->line_direction = to_dir;
			streamlineDeviation = 0.0;

			gmds::math::Point current_point, previous_point;
			vector<gmds::math::Point> finalLineDiscretization = surf_line->getDiscretizationPoints();
			vector<TCellID> finalLineTraversedTriangles = surf_line->getTraversedFaces();
			current_point = finalLineDiscretization[0];
			unsigned int contTri = 0;
			bool foundTri, temp;
			double lambda0, lambda1;
			math::Vector3d previous_dir, next_dir;
			math::Vector3d closest0, closest1, closest2, closest2Current;

			// finalLineDiscretization from AFromSlot towards to_slot
			for (unsigned int i = 1; i < finalLineDiscretization.size() - 1; i++) {
				previous_point = current_point;
				current_point = finalLineDiscretization[i];
				previous_dir = math::Vector3d(previous_point, current_point);
				next_dir = math::Vector3d(current_point, finalLineDiscretization[i + 1]);
				if (i == 1) {
					streamlineDeviation = streamlineDeviation + fabs(1.0 - previous_dir.dot(AFromSlot->direction));
				}
				if (i == finalLineDiscretization.size() - 2) {
					streamlineDeviation = streamlineDeviation + fabs(1.0 - next_dir.dot(to_slot->direction.opp()));
				}
				foundTri = false;
				while (contTri < finalLineTraversedTriangles.size() && !foundTri) {
					if (m_tool.isPntInTri(current_point, m_mesh->get<Face>(finalLineTraversedTriangles[contTri]), temp, temp, temp, lambda0, lambda1)) {
						vector<gmds::Node> current_verts = m_mesh->get<Face>(finalLineTraversedTriangles[contTri]).get<Node>();
						math::Cross2D cross_0 = (*m_field)[current_verts[0].id()];
						math::Cross2D cross_1 = (*m_field)[current_verts[1].id()];
						math::Cross2D cross_2 = (*m_field)[current_verts[2].id()];
						closest0 = cross_0.closestComponentVector(previous_dir);
						closest1 = cross_1.closestComponentVector(previous_dir);
						closest2 = cross_2.closestComponentVector(previous_dir);
						closest2Current = lambda0 * closest0 + lambda1 * closest1 + (1 - lambda0 - lambda1) * closest2;
						streamlineDeviation = streamlineDeviation + fabs(1.0 - previous_dir.dot(closest2Current));

						closest0 = cross_0.closestComponentVector(next_dir);
						closest1 = cross_1.closestComponentVector(next_dir);
						closest2 = cross_2.closestComponentVector(next_dir);
						closest2Current = lambda0 * closest0 + lambda1 * closest1 + (1 - lambda0 - lambda1) * closest2;
						streamlineDeviation = streamlineDeviation + fabs(1.0 - next_dir.dot(closest2Current));
					}
					contTri++;
				}
			}
			streamlineDeviation = streamlineDeviation / finalLineDiscretization.size();

			to_slot->lineDeviation = streamlineDeviation;
			AFromSlot->lineDeviation = streamlineDeviation;
		}
	}
}

void
SingGraphBuilder2DOriginal::computeStreamLine(SingularityPoint *AFromPnt,
                                              SingularityPoint::Slot *AFromSlot,
                                              SingularityPoint *&AToSingPnt,
                                              SingularityPoint::Slot *&AToSlot,
                                              math::Point &AToPnt,
                                              math::Vector3d &AToDir,
                                              std::vector<math::Point> &APoints,
                                              std::vector<TCellID> &ATriangles,
                                              int &AToCellDim,
                                              TCellID &AToCellID,
                                              double &streamlineDeviation,
                                              bool &AEndOnBnd,
                                              bool &AToSlotIsFree,
                                              bool &APntToCreate)
{
	if (m_withGlobalComments) cout << "begining computeStreamLine from point " << AFromPnt->getLocation().X() << " , " << AFromPnt->getLocation().Y() << endl;

	ATriangles.clear();
	APoints.clear();

	math::Point start_pnt = AFromSlot->location;     // starting point
	APoints.push_back(start_pnt);
	math::Vector3d start_dir = AFromSlot->direction;     // starting direction
	math::Vector3d prev_dir = AFromSlot->direction;      // prev direction used in the
	/*termination of extrapolation process (when we get into a
	confusing ball) */

	TCellID start_cell_id = AFromSlot->starting_cell_id;
	int start_cell_dim = AFromSlot->starting_cell_dim;

	math::Point current_pnt = start_pnt;
	math::Vector3d current_vec = start_dir;

	if (start_cell_dim == 0)
		cout << "node: " << start_cell_id << endl;
	else if (start_cell_dim == 1) {
		vector<Node> currentNodes = (m_mesh->get<Edge>(start_cell_id)).get<Node>();
		if (m_withGlobalComments) cout << "edge: " << start_cell_id << " between " << currentNodes[0].id() << " " << currentNodes[1].id() << endl;
	}

	math::Point start_dirPnt(start_pnt.X() + start_dir.X(), start_pnt.Y() + start_dir.Y(), start_pnt.Z() + start_dir.Z());

	//========================================================================
	// Singularity point we will be connecting to. It can be another
	// singularity point of the cross field or a geometric singularity point
	// created on the fly.
	//========================================================================
	SingularityPoint *found_pnt = 0;
	SingularityPoint::Slot *found_slot = 0;

	bool find_end = false;
	/* indicates that we reach a boundary point or line */
	bool end_on_boundary = false;
	/* indicates that we reach an existing singularity point*/
	// bool end_on_field_singularity = false;

	// We check some termination conditions on the boundary.
	if (start_cell_dim == 0) {
		Node currentNode = m_mesh->get<Node>(start_cell_id);
		if (m_mesh->isMarked(currentNode, m_mark_nodes_on_point) || m_mesh->isMarked(currentNode, m_mark_nodes_on_curve)) {
			find_end = true;
			end_on_boundary = true;
		}
	}
	else {     // we have necessarry start_cell_dim=1
		Edge currentEdge = m_mesh->get<Edge>(start_cell_id);
		if (m_mesh->isMarked(currentEdge, m_mark_edges_on_curve)) {
			find_end = true;
			end_on_boundary = true;
		}
	}
	//========================================================================
	// Main loop to create the singularity line
	//========================================================================
	// int nbwalk=0;
	while (!find_end) {
		// nbwalk++;
		TCellID next_cell_id = NullID;
		int next_cell_dim = -1;
		m_tool.findNextCell(start_pnt, start_dir, start_cell_dim, start_cell_id, next_cell_dim, next_cell_id);

		if (next_cell_dim == -1) {
			if (m_withGlobalComments) cout << "next_cell_dim == -1" << endl;
			// The cell defined by (start_cell_dim, start_cell_id) is on the boundary.
			find_end = true;
			end_on_boundary = true;
		}
		else if (next_cell_dim == 1) {
			if (m_withGlobalComments) cout << "next_cell_dim == 1" << endl;
			/* we are going along an edge.
			Our simple assumption is to follow this edge until reaching
			one of its end points and to compute the next direction at
			this point.*/
			// WARNING, we do not check the confusing ball!!!!!
			Edge currentEdge = m_mesh->get<Edge>(next_cell_id);
			std::vector<TCellID> adj_faces = currentEdge.getIDs<Face>();
			ATriangles.insert(ATriangles.end(), adj_faces.begin(), adj_faces.end());

			std::vector<Node> currentNodes = currentEdge.get<Node>();
			math::Vector3d v0(start_pnt, currentNodes[0].point());
			math::Vector3d v1(start_pnt, currentNodes[1].point());
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
			math::Vector3d devVect(start_pnt, next_node.point());
			devVect.normalize();
			m_tool.computeOutVectorAtPoint(next_node, start_dir, next_dir);
			start_dir.normalize();
			streamlineDeviation = streamlineDeviation + fabs(1.0 - start_dir.dot(devVect));

			// We assign the new value for the next step
			start_dir = next_dir;
			start_pnt = next_node.point();
			start_cell_dim = 0;
			start_cell_id = next_node.id();
			find_end = false;

			APoints.push_back(start_pnt);
		}
		else {     // general case, we are in a face
			Face currentFace = m_mesh->get<Face>(next_cell_id);
			// cout<<"next_cell_id "<<next_cell_id<<endl;
			ATriangles.push_back(currentFace.id());
			//==============================================================
			// CASE 1: DO WE ARE IN A FACE CONTAINING A SING. POINT?
			//==============================================================
			bool intersect_sing_point = false;

			SingularityPoint *next_sing_point = m_faces_to_singularity_on_surf[currentFace.id()];
			bool must_try_to_connect = false;
			if (next_sing_point != NULL &&         // face in a confusing ball ...
			    next_sing_point != AFromPnt) {     // ... of another singularity point
				if (m_withGlobalComments) cout << "face in a confusing ball of another singularity point" << endl;
				must_try_to_connect = true;
			}
			else if (next_sing_point != NULL &&         // face in the confusing ball ...
			         next_sing_point == AFromPnt) {     //... of the incoming singularity point
				if (m_withGlobalComments) cout << "face in a confusing ball of the incoming singularity point" << endl;
				// Warning: completly empiric, we just try to detect cyclic lines
				if (APoints.size() >= 100) must_try_to_connect = true;
			}

			if (must_try_to_connect) {
				if (m_withGlobalComments) cout << "must_try_to_connect" << endl;
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
						}     // if (slot_deviation < slot_epsilon) {
					}        // for (unsigned int i_slot = 0; !found_free_slot && i_slot < ....

					if (best_deviation != -2 && !cur_slots.empty()) {
						if (m_withGlobalComments) cout << "best_deviation!=-2 &&  !cur_slots.empty()" << endl;
						SingularityPoint::Slot *best_slot = cur_slots[best_slot_id];
						math::Vector3d slot_opp_dir = best_slot->direction.opp();
						if (m_withGlobalComments) {
							cout << " best_slot->direction " << best_slot->direction.X() << " " << best_slot->direction.Y() << endl;
							cout << "current_vec " << current_vec.X() << " " << current_vec.Y() << endl;
						}
						slot_opp_dir.normalize();
						if (best_slot->isLaunched) {
							// slot already assigned with a previous line (and direction)
							if (m_withGlobalComments) cout << "slot already assigned with a previous line (and direction)" << endl;

							math::Vector3d prev_line_dir = best_slot->line_direction;
							if (m_withGlobalComments)
								cout << " best_slot->line_direction " << best_slot->line_direction.X() << " " << best_slot->line_direction.Y() << endl;
							prev_line_dir.normalize();
							double prev_deviation = slot_opp_dir.dot(prev_line_dir);
							if (m_withGlobalComments) cout << "prev_deviation " << prev_deviation << endl;
							if (best_deviation > prev_deviation) {
								// the new alignment is better than the previous one
								if (m_withGlobalComments)
									cout << "best_deviation " << best_deviation << " > "
									     << "prev_deviation " << prev_deviation << ", the new alignment is better than the previous one" << endl;

								found_free_slot = false;
								found_compatible_slot = true;
								found_slot = best_slot;
								intersect_sing_point = true;
							}
							else {
								// We keep the previous association
								if (m_withGlobalComments) cout << "We keep the previous association" << endl;

								found_compatible_slot = false;
							}
						}          // if ( current_slot>isLaunched) {
						else {     // WE HAVE A FREE SLOT
							// We keep the previous association
							if (m_withGlobalComments) cout << "WE HAVE A FREE SLOT" << endl;

							found_free_slot = true;
							AToSlotIsFree = found_free_slot;
							found_compatible_slot = true;
							found_slot = best_slot;
							intersect_sing_point = true;
						}
						// HAVE WE FOUND THE END OF THE LINE??
						if (found_compatible_slot) {
							if (m_withGlobalComments) cout << "last    found_compatible_slot" << endl;
							// COMPATIBLE AND FREE SLOT
							find_end = true;
						}
					}
					slot_epsilon -= 0.1;
				}     // while (!found_compatible_slot && slot_epsilon < -0.4)
			}        //	if(must_try_to_connect) {

			//==============================================================
			// CASE 2: SO, WE JUST HAVE TO CROSS THE TRIANGLE (GENERAL CONF)
			//==============================================================
			// Does the current triangle has the same classif
			if (!intersect_sing_point) {
				math::Point out_pnt;
				math::Vector3d out_vec;
				TCellID out_cell_id;
				int out_cell_dim;
				if (m_withGlobalComments) {
					cout << "start_pnt " << start_pnt.X() << " " << start_pnt.Y() << endl;
					cout << "currentFace " << currentFace.id() << endl;
				}
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
				if (m_withGlobalComments) cout << "after traverseTriangle out_pnt= " << out_pnt.X() << " " << out_pnt.Y() << endl;
				APoints.push_back(out_pnt);

				// we progress to the next point, next vector and so next face too
				prev_dir = start_dir;     // we store the prev direction for slot
				// reconnection with balls
				start_pnt = out_pnt;
				start_dir = out_vec;
				start_cell_dim = out_cell_dim;
				start_cell_id = out_cell_id;
			}     // if (!intersect_line && !intersect_sing_point) {

			if (intersect_sing_point) {
				if (m_withGlobalComments) cout << "intersect_sing_point" << endl;
				find_end = true;
			}

			// post process, we just check whether we have arrived onto a geometric boundary
			if (start_cell_dim == 0) {
				if (m_withGlobalComments) cout << "post process start_cell_dim==0" << endl;
				Node currentNode = m_mesh->get<Node>(start_cell_id);
				if (m_mesh->isMarked(currentNode, m_mark_nodes_on_point) || m_mesh->isMarked(currentNode, m_mark_nodes_on_curve)) {
					find_end = true;
					end_on_boundary = true;
				}
			}
			else {     // we have necessarry start_cell_dim=1
				if (m_withGlobalComments) cout << "post process we have necessarry start_cell_dim=1" << endl;
				Edge currentEdge = m_mesh->get<Edge>(start_cell_id);
				if (m_mesh->isMarked(currentEdge, m_mark_edges_on_curve)) {
					find_end = true;
					end_on_boundary = true;
				}
			}
		}     // else { //general case, we are in a face
	}        // while(!find_end)

	//==============================================================
	// Update of out parameters
	//==============================================================
	// last followed direction
	AToDir = start_dir;
	AToPnt = start_pnt;
	AEndOnBnd = end_on_boundary;
	// the end point must be created if it has not been found
	APntToCreate = (found_pnt == 0);
	// singularity point data if we found an end point
	AToSingPnt = found_pnt;
	AToSlot = found_slot;
	AToCellDim = start_cell_dim;
	AToCellID = start_cell_id;

	/*
	if(found_slot!=0){
	   gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));

	   for (unsigned int j = 0; j < APoints.size()-1; j++){
	       gmds::Node n1 = m.newNode(APoints[j].X(), APoints[j].Y(), APoints[j].Z());
	       gmds::Node n2 = m.newNode(APoints[j + 1].X(), APoints[j + 1].Y(), APoints[j + 1].Z());
	       gmds::Face f  =  m.newTriangle(n1, n1, n2);
	   }
	   gmds::IGMeshIOService ioService(&m);
	   gmds::VTKWriter vtkWriter(&ioService);
	   vtkWriter.setCellOptions(gmds::N|gmds::F);
	   vtkWriter.setDataOptions(gmds::N|gmds::F);
	   std::string file_name = "testNow"+to_string(found_slot->direction[0])+".vtk";
	   vtkWriter.write(file_name);
	} */
}

void
SingGraphBuilder2DOriginal::backtrackSingularityLine(SurfaceSingularityLine *ALine,
                                                     SingularityPoint *AFromPnt,
                                                     SingularityPoint::Slot *AFromSlot,
                                                     SingularityPoint *AToPnt,
                                                     SingularityPoint::Slot *AToSlot,
                                                     gmds::math::Vector3d &to_dir,
                                                     bool &foundBackTrackPath)
{
	//==============================================================
	// We keep in mind the length of the line
	//==============================================================
	if (m_withGlobalComments) cout << "backtrackSingularityLine" << endl;
	// double line_length = ALine->length();

	SingularityPoint *arrival_sing = 0;
	SingularityPoint::Slot *arrival_slot;
	math::Point arrival_pnt;
	math::Vector3d arrival_dir;
	int arrival_cell_dim;
	TCellID arrival_cell_id;
	std::vector<math::Point> line_pnts;
	std::vector<TCellID> new_traversed_triangles;
	bool arrival_on_bnd;
	bool arrival_on_free_slot;
	bool arrival_pnt_to_create;

	vector<gmds::TCellID> modifiedFaces;
	double previousRad = m_confusing_distance;
	double redefRadius = 1.5;
	int cont = 0;     // here take care not to bump into a new confusing ball
	double streamlineDeviation = 0.0;

	while ((arrival_sing != AToPnt) && (cont < 10)) {
		computeStreamLine(AFromPnt, AFromSlot, arrival_sing, arrival_slot, arrival_pnt, arrival_dir, line_pnts, new_traversed_triangles, arrival_cell_dim,
		                  arrival_cell_id, streamlineDeviation, arrival_on_bnd, arrival_on_free_slot, arrival_pnt_to_create);

		if (arrival_sing != AToPnt) {
			cont++;
			redefineOneConfusingBall(AToPnt, modifiedFaces, previousRad, redefRadius);
			redefRadius = redefRadius + 0.5;
		}
	}

	if (arrival_sing != AToPnt) {
		if (m_withGlobalComments) {
			cout << "we have started from triangle " << AFromPnt->getMesh<Face>()[0].id() << " and we want to reach " << AToPnt->getMesh<Face>()[0].id() << endl;
			cout << "backtrack didn't reach destination point" << endl;
		}
		// throw GMDSException("Backtraking issue: We don't reach the departure point");
	}
	else {
		foundBackTrackPath = true;
		if (m_withGlobalComments) cout << "backtrack we have reached destination point" << endl;

		if (modifiedFaces.size() != 0) {
			/*std::string name = "ball_varredefReached"+std::to_string(cont);
			Variable<int>* modif_sing_ball = m_mesh->newVariable<int,GMDS_FACE>(name);
			for(unsigned int t=0; t<modifiedFaces.size(); t++){
			   Face f = m_mesh->get<Face>(modifiedFaces[t]);

			   SingularityPoint* sing = m_faces_to_singularity_on_surf[f.id()];
			   if (sing==0)
			      (*modif_sing_ball)[f.id()] = 0;
			   else
			      (*modif_sing_ball)[f.id()] = sing->index();
			}
			gmds::IGMeshIOService meshIoServ(m_mesh);
			gmds::VTKWriter writerB(&meshIoServ);
			writerB.setCellOptions(gmds::N|gmds::F);
			writerB.setDataOptions(gmds::N|gmds::F);
			std::stringstream file_name2;
			file_name2 <<m_output_directory_name<<"-confusing_balls_Reached.vtk";
			writerB.write(file_name2.str());  */
		}
	}

	if (cont > 0) {
		for (unsigned int i = 0; i < modifiedFaces.size(); i++) {
			m_faces_to_singularity_on_surf[modifiedFaces[i]] = 0;
		}
	}

	if (foundBackTrackPath) {
		std::vector<gmds::math::Point> old_pnts = ALine->getDiscretizationPoints();
		/*gmds::Mesh meshSing(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
		gmds::Node firstNo = meshSing.newNode(line_pnts[0].X(), line_pnts[0].Y(), line_pnts[0].Z());
		for(unsigned int itsp=1; itsp<line_pnts.size(); itsp++){
		    gmds::Node secondNo = meshSing.newNode(line_pnts[itsp].X(), line_pnts[itsp].Y(), line_pnts[itsp].Z());
  meshSing.newTriangle(firstNo, secondNo, secondNo);
 firstNo = secondNo;
		}

		gmds::IGMeshIOService ioServiceSing(&meshSing);
		gmds::VTKWriter vtkWriterSing(&ioServiceSing);
		vtkWriterSing.setCellOptions(gmds::N|gmds::F);
		vtkWriterSing.setDataOptions(gmds::N|gmds::F);

		std::string file_name = "HolesInSquareNewPoints-_"+std::to_string(cont)+"_"+std::to_string(AToSlot->line->getNumber())+".vtk";

		vtkWriterSing.write(file_name);
		//throw GMDSException("old_pnts"); */

		// old_pnts and line_pnts are traversed in opposite directions. Moreover,
		// they don't terminate at the same location. One goes from the singularity
		// point 1 to the boundary of the confusing ball 2 and the other does the
		// opposite.
		std::vector<gmds::math::Point> new_pnts;

		new_pnts.push_back(old_pnts[0]);
		math::Point last_point = old_pnts[old_pnts.size() - 1];

		for (unsigned int i_old = 1; i_old < old_pnts.size() - 1; i_old++) {
			math::Point current_pnt = old_pnts[i_old];
			math::Segment seg(line_pnts[0], line_pnts[1]);
			math::Point proj_pnt = seg.project(current_pnt);
			double proj_dist = current_pnt.distance(proj_pnt);

			for (unsigned int i_new = 2; i_new < line_pnts.size(); i_new++) {
				math::Segment seg_i(line_pnts[i_new - 1], line_pnts[i_new]);
				math::Point proj_pnt_i = seg_i.project(current_pnt);
				double proj_dist_i = current_pnt.distance(proj_pnt_i);
				if (proj_dist_i < proj_dist) {
					proj_dist = proj_dist_i;
					proj_pnt = proj_pnt_i;
				}
			}
			// We have our new point
			new_pnts.push_back(0.5 * current_pnt + 0.5 * proj_pnt);
		}

		new_pnts.push_back(last_point);
		ALine->setDiscretizationPoints(new_pnts);

		//==============================================================
		// Now we recompute the intersected faces.
		//==============================================================

		// We keep in mind the faces traversed the first time and during the
		// backtracking process. Indeed, we don't know which faces will be really
		// traversed at the end.
		std::vector<TCellID> all_traversed_faces = ALine->getTraversedFaces();
		all_traversed_faces.insert(all_traversed_faces.end(), new_traversed_triangles.begin(), new_traversed_triangles.end());
		std::set<TCellID> candidates;
		for (unsigned int i = 0; i < all_traversed_faces.size(); i++) {
			Face f = m_mesh->get<Face>(all_traversed_faces[i]);
			std::vector<Node> f_nodes = f.get<Node>();
			for (unsigned int j = 0; j < f_nodes.size(); j++) {
				Node nj = f_nodes[j];
				std::vector<TCellID> nj_faces = nj.getIDs<Face>();
				for (unsigned int k = 0; k < nj_faces.size(); k++) {
					candidates.insert(nj_faces[k]);
				}
			}
		}

		std::set<TCellID> set_of_traversed_faces;

		for (std::set<TCellID>::iterator it = candidates.begin(); it != candidates.end(); it++) {
			Face f = m_mesh->get<Face>(*it);
			std::vector<Node> f_nodes = f.get<Node>();
			math::Triangle t(f_nodes[0].point(), f_nodes[1].point(), f_nodes[2].point());
			bool found_pnt = false;
			for (unsigned int j = 0; j < new_pnts.size() && !found_pnt; j++) {
				math::Point pj = new_pnts[j];
				if (t.isIn(pj)) {
					found_pnt = true;
					set_of_traversed_faces.insert(f.id());
				}
				if (j != 0) {
					math::Point pk = new_pnts[j - 1];
					if (t.intersect(math::Segment(pj, pk))) {
						found_pnt = true;
						set_of_traversed_faces.insert(f.id());
					}
				}
			}
		}

		std::vector<TCellID> final_traversed_faces;
		final_traversed_faces.insert(final_traversed_faces.end(), set_of_traversed_faces.begin(), set_of_traversed_faces.end());
		ALine->setTraversedFaces(final_traversed_faces);
		to_dir = arrival_dir;
	}
}
