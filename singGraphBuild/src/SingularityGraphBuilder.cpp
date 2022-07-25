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
 * In this respect, the user's attention is drawn to the risks apassoiressociated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulatetNe,  and  that  also
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
 * SingularityGraphBuilder.cpp
 *
 *  Created on: 13 juil. 2014
 *      Author: bibi
 */
/*----------------------------------------------------------------------------*/
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/MeditReader.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/math/Chart.h>
#include <gmds/math/Cross.h>
#include <gmds/math/Quaternion.h>
/*----------------------------------------------------------------------------*/
#include <gmds/singGraphBuild/SingularityGraphBuilder.h>
/*----------------------------------------------------------------------------*/
#include <set>
#include <sstream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::execute()
{

	//==================================================================
	// BOOLEAN MARKS ARE INITIALIZED
	//==================================================================
	initMarks();

	//==================================================================
	// MARK BOUNDARY CELLS
	//==================================================================
	std::cout << "=============================================================" << std::endl;
	std::cout << "Mark boundary " << std::endl;
	BoundaryOperator boundaryOp(m_mesh);
	if (!boundaryOp.isValid()) {
		std::cout << "Invalid model for boundary operations" << std::endl;
		throw GMDSException("Invalid model for boundary operations");
	}

	boundaryOp.markCellOnGeometry(m_markFacesOnSurf, m_markEdgesOnSurf, m_markNodesOnSurf, m_markEdgesOnCurv, m_markNodesOnCurv, m_markNodesOnVert,
	                              m_markAloneNodes);

	Variable<int> *geom_var = m_mesh->newVariable<int, GMDS_NODE>("geometry");

	for (auto n_id : m_mesh->nodes()) {
		Node n = m_mesh->get<Node>(n_id);
		if (m_mesh->isMarked(n, m_markNodesOnVert))
			(*geom_var)[n.id()] = 3;
		else if (m_mesh->isMarked(n, m_markNodesOnCurv))
			(*geom_var)[n.id()] = 2;
		else if (m_mesh->isMarked(n, m_markNodesOnSurf))
			(*geom_var)[n.id()] = 1;
		else
			(*geom_var)[n.id()] = 0;
	}
	// VTKWriter<IGMesh> writer(*m_mesh);
	gmds::IGMeshIOService ioService(m_mesh);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N | gmds::F | gmds::R);
	vtkWriter.setDataOptions(gmds::N | gmds::F | gmds::R);
	vtkWriter.write("geom_mesh");
	// writer.write("geom_mesh", DIM3 | F | N | R);
	std::cout << "\t DONE" << std::endl;

	try {
		m_var_distance = m_mesh->getVariable<double, GMDS_NODE>("distance");
	}
	catch (GMDSException &e) {
		std::cout << "ERROR: NO DISTANCE!!! " << e.what() << std::endl;
		exit(0);
	}
	//==================================================================
	// CREATION OF QUATERNIONS (OR WE GET THEM)
	//==================================================================
	std::cout << "=============================================================" << std::endl;
	std::cout << "Initialization of the quaternion field" << std::endl;
	extractQuaternions();
	std::cout << "\t DONE" << std::endl;

	//==================================================================
	// COLOR FACES TO GATHER THEM ACCORDING TO THE GEOM. CLASSIFICATION
	//==================================================================
	std::cout << "=============================================================" << std::endl;
	std::cout << "color faces" << std::endl;
	colorFaces(m_markFacesOnSurf, m_markEdgesOnCurv);
	std::cout << "\t DONE" << std::endl;

	//========================================================================
	// STEP 1 - Detection of singular tetrahedra and storage
	//========================================================================
	std::cout << "=============================================================" << std::endl;
	std::cout << "Detection of singular tetrahedra" << std::endl;
	detectSingularTetrahedra();
	std::cout << "\t DONE" << std::endl;

	//========================================================================
	// STEP 1 - Singularity points inside the volume
	//========================================================================
	std::cout << "=============================================================" << std::endl;
	std::cout << "Creation of volume singularity points" << std::endl;
	createVolumeSingularityPoints();
	std::cout << "\t DONE" << std::endl;

	//========================================================================
	// STEP 1 - Singularity points on the surface
	//========================================================================
	std::cout << "=============================================================" << std::endl;
	std::cout << "Creation of Surface singularity points" << std::endl;
	initSingularityPointOnSurf();
	std::cout << "\t DONE" << std::endl;

	Variable<SurfaceSingularityPoint *> *ball_var = m_mesh->newVariable<SurfaceSingularityPoint *, GMDS_FACE>("sing_ball");

	for (auto f_id : m_mesh->faces()) {
		Face f = m_mesh->get<Face>(f_id);

		SurfaceSingularityPoint *sing = m_face_to_singularity_on_surf[f.id()];

		(*ball_var)[f.id()] = sing;
	}
	// VTKWriter<IGMesh> writerB(*m_mesh);
	// gmds::IGMeshIOService ioService(&m);
	gmds::VTKWriter writerB(&ioService);
	writerB.setCellOptions(gmds::N | gmds::F);
	writerB.setDataOptions(gmds::N | gmds::F);
	writerB.write("confusing_balls");
	// writerB.write("confusing_balls", DIM3 | F | N);
	std::cout << "\t DONE" << std::endl;
	//========================================================================
	// STEP 2 - Singularity lines inside the volume
	//========================================================================
	// at this point we have made all cluster of sing. points
	// now we need to compute the lines of sing into separatrices
	std::cout << "=============================================================" << std::endl;
	std::cout << "Creation of volume singularity lines" << std::endl;
	createVolumeSingularityLines();
	std::cout << "\t DONE" << std::endl;

	writeOutputSingle("inner_skeleton");
	std::cout << "INNER SKELETON BUILD" << std::endl;

	//========================================================================
	// STEP 3 - Geometrical features are added in the singularity graph
	//========================================================================
	std::cout << "=============================================================" << std::endl;
	std::cout << "Addition of geometrical features" << std::endl;
	addGeometryToSingularitGraph();
	std::cout << "\t DONE" << std::endl;
	writeOutputSingle("with_geom");

	//========================================================================
	// STEP 4 - Boundary lines from the frame field are built
	//========================================================================
	std::cout << "=============================================================" << std::endl;
	std::cout << "Creation of surface frame lines" << std::endl;
	createBoundarySingularityGraphFromField();
	std::cout << "\t DONE" << std::endl;
	writeOutputSingle("with_surface_field");

	//========================================================================
	// STEP 5 - Other boundary lines to be built
	//========================================================================

	//========================================================================
	// STEP 6 - New inner lines
	//========================================================================

	//
	//	std::cout<<"on est sorti de la phase de creation des sep au bord issues de singularites du champs"<<std::endl;
	//	std::cout<<"on entre dans la phase de creation des autres separatrices au bord"<<std::endl;
	//
	//	// on recupere les separatrices qui ne terminent pas sur deux singularites, on les range dans un vecteur
	//	// de separatrices a traiter
	//	std::vector<Separatrix3D> unfinished_sep;
	//
	//	int m_markEndEdgeOfSepatrices = mesh.newMark();
	//	for (unsigned int i = 0;i < separatrices.size(); i++){
	//		Separatrix3D sepI = separatrices[i];
	//		if (sepI.isToBeAssociatedWithAnotherSep){
	//			unfinished_sep.push_back(sepI);
	//		}
	//		else{
	//			//on marque les aretes terminant les separatrices
	//			Node n0 = mesh.get<Node>(sepI.NodeEndID[0]);
	//			Node n1 = mesh.get<Node>(sepI.NodeEndID[1]);
	//			Edge endEdge =0;
	//			endEdge= getEdge(n0,n1);
	//			if(endEdge!=0)
	//				mesh.m_mark(endEdge,m_markEndEdgeOfSepatrices);
	//		}
	//	}
	//
	//	while(!unfinished_sep.empty()){
	//		std::cout<<"on a "<<unfinished_sep.size()<<" separatrices incompletes"<<std::endl;
	//		//on recupere la derniere separatrice
	//		Separatrix3D current_sep = unfinished_sep[unfinished_sep.size()-1];
	//		unfinished_sep.pop_back();
	//		std::cout<<"Nb singularities = "<<current_sep.SingNumber.size()<<std::endl;
	//		if(current_sep.SingNumber.size()>1){
	//			continue;
	//		}
	//		for(unsigned int k=0;k<current_sep.SingNumber.size();k++)
	//			std::cout<<"  --> "<<current_sep.SingNumber[k]<<std::endl;
	//
	//		//on recupere les noeuds de l'arete de maillage contenant la fin de la separatrice
	//
	//		// on regarde si on part d'une arete ou d'un sommet
	//		//on teste si une separatrice match avec
	//		int matchingIndex;
	//		if(getMatchingBoundaryFreeSep(current_sep,unfinished_sep,
	//				matchingIndex,m_markBordFace,m_markBordEdge,m_markBordNode,
	//				m_markEndEdgeOfSepatrices,
	//				singularities)==true){
	//			std::cout<<"===>> La separatrice "<<matchingIndex<<" matche"<<std::endl;
	//
	//			// On recupere les deux separatrices qui nous interesse
	//			Separatrix3D smatch;
	//			for(unsigned int i=0;i<separatrices.size();i++){
	//				Separatrix3D si = separatrices[i];
	//				if(si.getSepNumber()==matchingIndex)
	//					smatch = si;
	//			}
	//
	//			std::cout<<"Point de depart ("<<current_sep.Sens1X[current_sep.Sens1X.size()-1]<<", "
	//										<<current_sep.Sens1Y[current_sep.Sens1Y.size()-1]<<", "
	//										<<current_sep.Sens1Z[current_sep.Sens1Z.size()-1]<<") "<<std::endl;
	//
	//			std::cout<<"Point final ("<<smatch.Sens1X[smatch.Sens1X.size()-1]<<", "
	//										<<smatch.Sens1Y[smatch.Sens1Y.size()-1]<<", "
	//										<<smatch.Sens1Z[smatch.Sens1Z.size()-1]<<") "<<std::endl;
	//
	//			std::cout<<" --> realisation du matching "<<std::endl;
	//			matchCurves(current_sep,smatch,
	//					singularities, separatrices,
	//					m_markBordFace,m_markBordEdge,m_markBordNode,m_markbordSing,
	//					m_markEndEdgeOfSepatrices);
	//
	//			//on retire smatch des separatrices a traiter aussi
	//			int rem_index =-1;
	//			for(unsigned int i=0;i<unfinished_sep.size();i++){
	//				if(unfinished_sep[i].getSepNumber()==smatch.getSepNumber())
	//					rem_index = i;
	//			}
	//			if(rem_index == unfinished_sep.size()-1)
	//				unfinished_sep.pop_back();
	//			else{
	//				//on deplace le dernier elt
	//				unfinished_sep[rem_index] = unfinished_sep[unfinished_sep.size()-1];
	//				unfinished_sep.pop_back();
	//			}
	//			std::cout<<" --> realisation du matching DONE"<<std::endl;
	//		}
	//		else{
	//			std::cout<<"===>> aucune separatrice ne correspond"<<std::endl;
	//			Separatrix3D newSep = spreadCurve(current_sep,
	//					singularities, separatrices,
	//					m_markBordFace,m_markBordEdge,m_markBordNode,m_markbordSing,m_markEndEdgeOfSepatrices);
	//			//on ajoute la nouvelle separatrice au debut pour la traiter en dernier
	//			unfinished_sep.push_back(unfinished_sep[unfinished_sep.size()-1]);
	//			Separatrix3D prev0 = unfinished_sep[0];
	//			for(unsigned int is=unfinished_sep.size()-1;is>0;is--){
	//				unfinished_sep[is]=unfinished_sep[is-1];
	//			}
	//			unfinished_sep[0]=newSep;
	//		}
	//		std::cout<<">>>>>>>>>>>>>>>>>>>>>>>>> MESH SKELETON"<<std::endl;
	//		buildMeshSkeleton(separatrices);
	//		writeSkeleton();
	//
	//		//creation d'une singularite
	////		std::cout<<"on a une pos de: "<<separatrices[i].Sens1X[separatrices[i].Sens1X.size() - 1]<<";
	///"<<separatrices[i].Sens1Y[separatrices[i].Sens1X.size() - 1]<<"; "<<separatrices[i].Sens1Z[separatrices[i].Sens1X.size() - 1]<<std::endl;
	//
	//	}

	//==================================================================
	// BOOLEAN MARKS ARE FINALLY CLEANED
	//==================================================================
	std::cout << "=============================================================" << std::endl;
	std::cout << "Boolean marks cleaning" << std::endl;
	cleanMarks();
	std::cout << "\t DONE" << std::endl;

	std::cout << "=============================================================" << std::endl;
	writeOutput("output_FRAMEFIELD");
	std::cout << "Skeleton done!" << std::endl;
	std::cout << "=============================================================" << std::endl;
}

/*---------------------------------------------------------------------------*/
void
SingularityGraphBuilder::detectSingularTetrahedra()
{

	m_2SingTetIDsAlongSurf.clear();
	m_3SingTetIDs.clear();
	m_region_singularity_type.clear();

	for (auto r_id : m_mesh->regions()) {
		Region current = m_mesh->get<Region>(r_id);
		std::vector<TCellID> nodeIDs = current.getIDs<Node>();
		int ID1 = nodeIDs[0];
		int ID2 = nodeIDs[1];
		int ID3 = nodeIDs[2];
		int ID4 = nodeIDs[3];

		math::Quaternion q1 = (*m_var_quatern)[ID1];
		math::Quaternion q2 = (*m_var_quatern)[ID2];
		math::Quaternion q3 = (*m_var_quatern)[ID3];
		math::Quaternion q4 = (*m_var_quatern)[ID4];
		int singularity_type = math::Quaternion::testSingularity(q1, q2, q3, q4);

		m_region_singularity_type[current.id()] = singularity_type;

		if (singularity_type == 3) {
			m_3SingTetIDs.push_back(current.id());
		}
		else if (singularity_type == 2) {
			std::vector<Face> current_faces = current.get<Face>();
			bool along_surf = false;
			for (unsigned int i = 0; i < current_faces.size(); i++) {
				Face current_face = current_faces[i];
				if (m_mesh->isMarked(current_face, m_markFacesOnSurf)) along_surf = true;
			}

			if (along_surf) m_2SingTetIDsAlongSurf.push_back(current.id());
		}     // else if (singularity_type == 2)

	}     // for (; !it_regions.isDone(); it_regions.next())

	std::cout << "Nb 2-singular tet (along surf): " << m_2SingTetIDsAlongSurf.size() << std::endl;
	std::cout << "Nb 3-singular tet: " << m_3SingTetIDs.size() << std::endl;
}
/*---------------------------------------------------------------------------*/
void
SingularityGraphBuilder::initMarks()
{

	// set of marks for knowing which cells are already processed in the
	// sing. graph building algorithm.
	m_markClusterSingDone = m_mesh->newMark<Region>();
	m_markBordVolSingForFaces = m_mesh->newMark<Face>();
	m_markBordVolSingForEdges = m_mesh->newMark<Edge>();
	m_markBordVolSingForNodes = m_mesh->newMark<Node>();

	m_markBordSurfSingForFaces = m_mesh->newMark<Face>();
	m_markBordSurfSingForEdges = m_mesh->newMark<Edge>();
	m_markBordSurfSingForNodes = m_mesh->newMark<Node>();
	// set of marks for the geometric classification
	m_markFacesOnSurf = m_mesh->newMark<Face>();
	m_markEdgesOnSurf = m_mesh->newMark<Edge>();
	m_markNodesOnSurf = m_mesh->newMark<Node>();

	m_markEdgesOnCurv = m_mesh->newMark<Edge>();
	m_markNodesOnCurv = m_mesh->newMark<Node>();

	m_markNodesOnVert = m_mesh->newMark<Node>();

	m_markAloneNodes = m_mesh->newMark<Node>();
}
/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::cleanMarks()
{
	m_mesh->unmarkAll<Node>(m_markNodesOnVert);
	m_mesh->unmarkAll<Node>(m_markNodesOnCurv);
	m_mesh->unmarkAll<Node>(m_markNodesOnSurf);
	m_mesh->unmarkAll<Node>(m_markAloneNodes);
	m_mesh->unmarkAll<Edge>(m_markEdgesOnSurf);
	m_mesh->unmarkAll<Edge>(m_markEdgesOnCurv);
	m_mesh->unmarkAll<Face>(m_markFacesOnSurf);
	m_mesh->unmarkAll<Region>(m_markClusterSingDone);

	m_mesh->unmarkAll<Face>(m_markBordVolSingForFaces);
	m_mesh->unmarkAll<Edge>(m_markBordVolSingForEdges);
	m_mesh->unmarkAll<Node>(m_markBordVolSingForNodes);

	m_mesh->unmarkAll<Face>(m_markBordSurfSingForFaces);
	m_mesh->unmarkAll<Edge>(m_markBordSurfSingForEdges);
	m_mesh->unmarkAll<Node>(m_markBordSurfSingForNodes);

	m_mesh->freeMark<Node>(m_markNodesOnVert);
	m_mesh->freeMark<Node>(m_markNodesOnCurv);
	m_mesh->freeMark<Node>(m_markNodesOnSurf);
	m_mesh->freeMark<Node>(m_markAloneNodes);
	m_mesh->freeMark<Edge>(m_markEdgesOnSurf);
	m_mesh->freeMark<Edge>(m_markEdgesOnCurv);

	m_mesh->freeMark<Face>(m_markFacesOnSurf);

	m_mesh->freeMark<Region>(m_markClusterSingDone);

	m_mesh->freeMark<Face>(m_markBordVolSingForFaces);
	m_mesh->freeMark<Edge>(m_markBordVolSingForEdges);
	m_mesh->freeMark<Node>(m_markBordVolSingForNodes);

	m_mesh->freeMark<Face>(m_markBordSurfSingForFaces);
	m_mesh->freeMark<Edge>(m_markBordSurfSingForEdges);
	m_mesh->freeMark<Node>(m_markBordSurfSingForNodes);
}
///*----------------------------------------------------------------------------*/
// math::Vector3d<TCoord,3> SingularityGraphBuilder::
// findStableDirection(math::Chart& AT1,math::Chart& AT2,math::Chart& AT3)
//{
//
//}
///*----------------------------------------------------------------------------*/
// void SingularityGraphBuilder::
// alignAlongDirection(math::Chart& AT1,math::Vector3d<TCoord,3>& ADir)
//{
//
//}
/*----------------------------------------------------------------------------*/
SingularityPoint *
SingularityGraphBuilder::createSurfaceSingularityPointFromVolumeInfo(Face &FCurr, const bool ACreatedFromVolumeLine)
{

	SurfaceSingularityPoint *newSingPoint = m_graph.newSurfacePoint();
	m_mesh->mark(FCurr, m_markBordSurfSingForFaces);

	// the face FCurr contains the singularity but one of its incident surface
	// face can belong to the singularity cluster too.
	// We have to check that at first

	std::vector<Edge> current_edges = FCurr.get<Edge>();
	for (unsigned int i_edges = 0; i_edges < current_edges.size(); i_edges++)
		m_mesh->mark(current_edges[i_edges], m_markBordSurfSingForEdges);

	std::vector<Node> current_nodes = FCurr.get<Node>();
	for (unsigned int i_nodes = 0; i_nodes < current_nodes.size(); i_nodes++)
		m_mesh->mark(current_nodes[i_nodes], m_markBordSurfSingForNodes);

	// we need to create the directions of separatrices, the number of separatrices, etc...
	// We start by computing the normal to the triangle and reprojecting triads to this normal
	math::Vector3d normal = getInputNormal(FCurr, FCurr.get<Region>()[0]);

	// Now we reproject triads
	TCellID ID1 = FCurr.get<Node>()[0].id();
	TCellID ID2 = FCurr.get<Node>()[1].id();
	TCellID ID3 = FCurr.get<Node>()[2].id();

	//==================================================================
	//==================================================================
	//==================================================================
	math::Chart tq[3];
	tq[0] = math::Chart((*m_var_quatern)[ID1]);
	tq[1] = math::Chart((*m_var_quatern)[ID2]);
	tq[2] = math::Chart((*m_var_quatern)[ID3]);

	std::vector<math::Vector3d> vecRepLoc;
	// for each chart k (0-->2), and each of its vector (i:0-->2)
	for (unsigned int k = 0; k < 3; k++) {
		for (unsigned int i = 0; i < 3; i++) {

			double dot_ki = tq[k].get(i).dot(normal);

			if (dot_ki < 0.9 && dot_ki > (-0.9)) {
				math::Vector3d v_ki = tq[k].get(i);
				vecRepLoc.push_back(v_ki);
			}
		}
	}
	//==================================================================
	// We have here a vec orthogonal for each
	// we need to launch a function taking the vector of rep, the normal, the ID of nodes,
	// and get back the pos of the sing, and vectors for the sep slots
	math::Point posSing;
	std::vector<math::Point> posSep;
	std::vector<math::Vector3d> dirSep;
	std::vector<int> sepInFaceID;

	// std::cout << "got in through first one" << std::endl;
	computeBoundSingInfo(m_markFacesOnSurf, FCurr, vecRepLoc, normal, posSing, posSep, dirSep, sepInFaceID);

	//==================================================================
	// We udpate the singularity point data
	//==================================================================
	newSingPoint->setLocation(posSing);
	newSingPoint->addMeshFace(FCurr);
	//==================================================================
	// The slot for the line inside the volume is initialized first
	if (!ACreatedFromVolumeLine) {
		// We have to define the direction vector to follow inside the
		// volume.

		// TODO CHANGE      newSingPoint->newSlot(posSing, normal.opp(), FCurr, false);
	}
	//==================================================================
	// Slots for surface lines are initialized
	for (unsigned int i = 0; i < posSep.size(); i++) {
		Face f = m_mesh->get<Face>(sepInFaceID[i]);
		// TODO CHANGEnewSingPoint->newSlot(posSep[i], dirSep[i], f, true);
	}

	//==================================================================
	// now we computed the confusion ball around this point.
	// For that, we use the distance field
	initConfusionBallOfASurfaceSingularityPoint(newSingPoint);
	// Now we traverse adjacent faces until getting faces which are in the good distance

	std::cout << "New surface singularity point starting from: " << posSing << "!!" << std::endl;

	return newSingPoint;
}

/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::initConfusionBallOfASurfaceSingularityPoint(SurfaceSingularityPoint *ASingPoint)
{
	Face f = ASingPoint->getMeshFace();
	std::vector<Node> current_nodes = f.get<Node>();

	math::Point current_point = current_nodes[0].point();

	double max_d = (*m_var_distance)[current_nodes[0].id()];

	for (unsigned int i = 1; i < current_nodes.size(); i++) {
		double dist_i = (*m_var_distance)[current_nodes[i].id()];
		if (dist_i > max_d) max_d = dist_i;
	}

	math::Point ref_point = f.center();

	std::vector<Face> done_faces;
	std::vector<Face> to_do_faces;
	int mark_done = m_mesh->newMark<Face>();
	done_faces.push_back(f);
	to_do_faces.push_back(f);
	m_mesh->mark(f, mark_done);

	while (!to_do_faces.empty()) {
		Face current = to_do_faces.back();
		to_do_faces.pop_back();
		// We get adjacent boundary faces
		std::vector<Face> adj_faces = getAdjacentFacesByNodes(current, m_markFacesOnSurf);

		for (unsigned int i = 0; i < adj_faces.size(); i++) {
			Face adj_face_i = adj_faces[i];
			math::Point center_i = adj_face_i.center();
			if (!m_mesh->isMarked(adj_face_i, mark_done) &&     // not already added
			    center_i.distance(ref_point) < max_d)           // not too far
			{
				to_do_faces.push_back(adj_face_i);
				m_mesh->mark(adj_face_i, mark_done);
				done_faces.push_back(adj_face_i);
			}

		}     // for (unsigned int i = 0; i < adj_faces.size(); i++)

	}     // while (!to_do_faces.empty())

	//====================================================================================
	// Now, we can build the confusing information
	//====================================================================================
	for (unsigned int i = 0; i < done_faces.size(); i++) {
		Face current = done_faces[i];
		m_face_to_singularity_on_surf[current.id()] = ASingPoint;

		m_mesh->unmark(current, mark_done);
	}
	//====================================================================================
	// Faces were just unmmarked and so, boolean mark must be made free
	//====================================================================================

	m_mesh->freeMark<Face>(mark_done);
}

/*----------------------------------------------------------------------------*/
std::vector<Face>
SingularityGraphBuilder::getAdjacentFacesByNodes(Face &AFace, const int AMark)
{
	int init_color = (*m_var_color)[AFace.id()];
	std::set<Face> adj_faces;
	std::vector<Node> adj_nodes = AFace.get<Node>();

	for (unsigned int j = 0; j < adj_nodes.size(); j++) {
		Node n_j = adj_nodes[j];
		std::vector<Face> faces_nj = n_j.get<Face>();
		for (unsigned int j = 0; j < faces_nj.size(); j++) {
			Face f_j = faces_nj[j];
			int color_j = (*m_var_color)[f_j.id()];
			if (f_j.id() != AFace.id() &&           // fj is not AFace
			    m_mesh->isMarked(f_j, AMark) &&     // fj is on the surface
			    color_j == init_color)              // fj and AFace are on the same surface
			{
				adj_faces.insert(f_j);
			}
		}
	}
	std::vector<Face> adj;
	adj.insert(adj.end(), adj_faces.begin(), adj_faces.end());
	return adj;
}
/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::createOneSingularityLineFrom(SingularityPoint *ASingPoint, SingularityPoint::Slot *ASlot, const int AMarkFaceUsedForSep)
{
	std::cout << "===== START OF AN INNER SEP CREATION =====" << std::endl;

	// we have a new sep to create here starting with Ftmp
	// first we look if we have to create a new sing, on the boundary,
	// or to retrieve the sing inside the mesh volume

	SingularityPoint *starting_singularity_pnt = ASingPoint;
	Face starting_face;     // TODO CHANGE = ASlot->starting_triangle;

	std::vector<Region> adj_regions = starting_face.get<Region>();
	Region firstReg;
	for (unsigned int i_reg = 0; i_reg < adj_regions.size(); i_reg++) {
		if (!m_mesh->isMarked(adj_regions[i_reg], m_markClusterSingDone)) firstReg = adj_regions[i_reg];
	}

	SingularityLine *new_line = m_graph.newVolumeLine();
	new_line->setNumber(m_graph.getNbLines());
	new_line->addSingularityPoint(starting_singularity_pnt);
	new_line->addDiscretizationPoint(starting_singularity_pnt->getLocation());
	ASlot->isLaunched = true;
	ASlot->line = new_line;

	Face FCurr = starting_face;
	// we get the face center
	math::Point currentPnt = ASlot->location;

	// tetrahedral element in which we currently propagate/build the sing. line
	Region currentTet = firstReg;

	//=====================================================================
	// at this point we have a sep, a point to insert, and a region to
	// continue the sep with we will now iteratively continue until reaching
	// either a boundary or a singularity point
	//=====================================================================
	bool completeSingLine = false;
	std::cout << "starting ReachedEnd" << std::endl;

	while (!completeSingLine) {
		new_line->addDiscretizationPoint(currentPnt);
		m_mesh->mark(FCurr, AMarkFaceUsedForSep);
		std::cout << "We go from face " << FCurr.id() << " into tetrahedron " << currentTet.id() << std::endl;
		// Now look for the opposite face
		int nbOfNewFacesOkay = 0;

		std::vector<Face> current_tet_faces = currentTet.get<Face>();
		int StopForLoop = 0;
		for (unsigned int i = 0; i < 4; i++) {
			Face current_tet_face = currentTet.get<Face>()[i];
			if (!StopForLoop && current_tet_face.id() != FCurr.id()) {

				if (isSingularFace(current_tet_face)) {
					StopForLoop = 1;
					// We have the new face
					nbOfNewFacesOkay++;
					FCurr = current_tet_face;
					std::vector<Node> current_face_nodes = FCurr.get<Node>();
					std::cout << "We are in face " << FCurr.id() << " and tet" << currentTet.id() << std::endl;

					// we get the center of the face as next discretization point
					currentPnt = FCurr.center();

					// Now, we distinguish the case of a boundary face and an inner face
					if (FCurr.get<Region>().size() == 1) {     // BOUNDARY
						std::cout << "== Case 1 - Geom. boundary face ==" << std::endl;
						// We reached the domain boundary, so we will have to stop after this point creation
						completeSingLine = true;
						std::vector<SurfaceSingularityPoint *> surf_sing = m_graph.getSurfacePoints();

						for (unsigned int i_sing = 0; i_sing < surf_sing.size(); i_sing++) {
							SurfaceSingularityPoint *surf_sing_i = surf_sing[i_sing];
							Face sing_face = surf_sing_i->getMeshFace();
							if (sing_face.id() == FCurr.id()) {

								// IMPORTANT Fill the singularity slot !!!!
								bool canAdd = surf_sing_i->addLineFromVolume(new_line);
								if (canAdd) {
									// to finish the new_line on this sing. point, we must have a compatible free slot
									new_line->addSingularityPoint(surf_sing_i);
									new_line->addDiscretizationPoint(surf_sing_i->getLocation());
									// the line is finished
									completeSingLine = true;
								}
							}
						}

						m_mesh->mark(FCurr, AMarkFaceUsedForSep);
					}     // if (FCurr.get<Region>().size() == 1)
					else if (m_mesh->isMarked(current_tet_face, m_markBordVolSingForFaces))
					//((m_mesh->isMarked(current_tet_face->get<Region>()[0], m_markClusterSingDone)) &&
					//(!m_mesh->isMarked(current_tet_face->get<Region>()[1], m_markClusterSingDone)))
					//||
					//((!m_mesh->isMarked(current_tet_face->get<Region>()[0], m_markClusterSingDone)) &&
					//(m_mesh->isMarked(current_tet_face->get<Region>()[1], m_markClusterSingDone))))
					{
						// The current_tet_face is incident to one region in a cluster and one out
						// We reached a sing. cluster, so we will have to stop after this point creation
						std::cout << "== Case 2 - Cluster boundary face ==" << std::endl;

						new_line->addDiscretizationPoint(currentPnt);
						m_mesh->mark(FCurr, AMarkFaceUsedForSep);
						std::cout << "Nb adj regions: " << FCurr.get<Region>().size() << std::endl;
						// need to find the corresponding cluster
						int IDOfTetToTest;
						if (FCurr.get<Region>()[0].id() == currentTet.id())
							currentTet = FCurr.get<Region>()[1];
						else
							currentTet = FCurr.get<Region>()[0];
						IDOfTetToTest = currentTet.id();
						int point_index = 0;
						std::vector<VolumeSingularityPoint *> singularities = m_graph.getVolumePoints();
						std::cout << "Nb singularities " << singularities.size() << std::endl;
						for (unsigned int i = 0; i < singularities.size(); i++) {
							// for each singularity point, we check its mesh cells
							std::cout << "-- i: " << i << std::endl;
							std::cout << "-- sing: " << singularities[i] << std::endl;
							std::vector<Region> current_cells = singularities[i]->getMeshRegions();

							std::cout << "-- nb reg: " << current_cells.size() << std::endl;
							for (unsigned int j = 0; j < current_cells.size(); j++) {
								if (current_cells[j].type() == GMDS_TETRA && current_cells[j].id() == IDOfTetToTest) {
									// on a trouve le bon cluster
									point_index = i;
								}
							}
						}
						std::cout << "point index     : " << point_index << std::endl;
						std::cout << "nb singularities: " << singularities.size() << std::endl;

						// new_line->addDiscretizationPoint(singularities[point_index]->getLocation());

						// IMPORTANT Fill the singularity slot !!!!
						bool canAdd = singularities[point_index]->addLine(new_line, FCurr);
						if (canAdd) {
							// to finish the new_line on this sing. point, we must have a compatible free slot
							new_line->addSingularityPoint(singularities[point_index]);
							new_line->addDiscretizationPoint(singularities[point_index]->getLocation());
							// the line is finished
							completeSingLine = true;
						}
					}
					else {
						std::cout << "== Case 3 - Inner normal face ==" << std::endl;
						// continue business as usual
						if (FCurr.get<Region>()[0].id() == currentTet.id()) {
							currentTet = FCurr.get<Region>()[1];
						}
						else {
							currentTet = FCurr.get<Region>()[0];
						}
					}
				}
			}     // if (!StopForLoop && current_tet_face.id() != FCurr.id())
		}

		std::cout << "|good faces| = " << nbOfNewFacesOkay << std::endl;
	}     // while (!completeSingLine)
	std::cout << "A sing. line has been built" << std::endl;
	// writeOutput("singularity_line");
}

/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::computeBoundSingInfo(const int AFaceMark,
                                              Face &AFace,
                                              std::vector<math::Vector3d> &AVecRep,
                                              math::Vector3d &ADirNormal,
                                              math::Point &APosSing,
                                              std::vector<math::Point> &posSep,
                                              std::vector<math::Vector3d> &dirSep,
                                              std::vector<int> &sepInFaceID)
{
	std::cout << "== Entree dans computeboundSingInfo ==" << std::endl;

	// le vecteur orthogonal a vecRepresentationLoc[0] dans le plan orthogonal a dirMoy
	math::Vector3d orthogBase = ADirNormal.cross(AVecRep[0]);
	std::vector<Node> current_nodes = AFace.get<Node>();
	Node node1 = current_nodes[0];
	Node node2 = current_nodes[1];
	Node node3 = current_nodes[2];

	std::cout << "Node 1 " << node1.id() << ": (" << node1.X() << ", " << node1.Y() << ", " << node1.Z() << ")" << std::endl;
	std::cout << "Node 2 " << node2.id() << ": (" << node2.X() << ", " << node2.Y() << ", " << node2.Z() << ")" << std::endl;
	std::cout << "Node 3 " << node3.id() << ": (" << node3.X() << ", " << node3.Y() << ", " << node3.Z() << ")" << std::endl;

	std::cout << "Normal Vec: " << ADirNormal << std::endl;
	std::cout << "First  Vec: " << AVecRep[0] << std::endl;
	std::cout << "Second Vec: " << orthogBase << std::endl;
	math::Vector3d vecRepresentationLoc[3];
	for (unsigned int i = 0; i < 3; i++) {
		vecRepresentationLoc[i] = AVecRep[i];
	}

	for (unsigned int i = 1; i < 3; i++) {
		// on s'occupe du point i
		double prodScalTmp = AVecRep[0].dot(AVecRep[i]);
		double angleTmp = acos(prodScalTmp);
		double prodScalTmpBis = orthogBase.dot(AVecRep[i]);
		if (prodScalTmpBis < 0.0) {
			angleTmp = 4.0 * acos(0.0) - angleTmp;
		}
		angleTmp = 4.0 * angleTmp;
		while (angleTmp >= 4 * acos(0.0)) {
			angleTmp -= 4 * acos(0.0);
		}
		while (angleTmp < 0.0) {
			angleTmp += 4 * acos(0.0);
		}
		std::cout << "Pour le point " << i << " on a un angle final de " << angleTmp << std::endl;
		// ici on doit calculer la matrice de rotation
		double cTmp = cos(angleTmp);
		double sTmp = sin(angleTmp);
		std::cout << "Pour le point " << i << " on a un cos de " << cTmp << " et sin de " << sTmp << std::endl;
		double matrixRot[3][3];
		matrixRot[0][0] = ADirNormal.X() * ADirNormal.X() + (1 - ADirNormal.X() * ADirNormal.X()) * cTmp;
		matrixRot[0][1] = ADirNormal.X() * ADirNormal.Y() * (1 - cTmp) - ADirNormal.Z() * sTmp;
		matrixRot[0][2] = ADirNormal.X() * ADirNormal.Z() * (1 - cTmp) + ADirNormal.Y() * sTmp;
		matrixRot[1][0] = ADirNormal.X() * ADirNormal.Y() * (1 - cTmp) + ADirNormal.Z() * sTmp;
		matrixRot[1][1] = ADirNormal.Y() * ADirNormal.Y() + (1 - ADirNormal.Y() * ADirNormal.Y()) * cTmp;
		matrixRot[1][2] = ADirNormal.Z() * ADirNormal.Y() * (1 - cTmp) - ADirNormal.X() * sTmp;
		matrixRot[2][0] = ADirNormal.X() * ADirNormal.Z() * (1 - cTmp) - ADirNormal.Y() * sTmp;
		matrixRot[2][1] = ADirNormal.Z() * ADirNormal.Y() * (1 - cTmp) + ADirNormal.X() * sTmp;
		matrixRot[2][2] = ADirNormal.Z() * ADirNormal.Z() + (1 - ADirNormal.Z() * ADirNormal.Z()) * cTmp;

		// multiplication du vecteur avec la matrice de rotation
		double vecResTmp[3];
		for (unsigned int j = 0; j < 3; j++) {
			double sumTmp = 0.0;
			for (unsigned int k = 0; k < 3; k++) {
				// sumTmp += matrixRotTmp[j][k] * vecRepresentationLoc[i][k];
				sumTmp += matrixRot[j][k] * AVecRep[0][k];
			}
			vecResTmp[j] = sumTmp;
		}

		vecRepresentationLoc[i].setX(vecResTmp[0]);
		vecRepresentationLoc[i].setY(vecResTmp[1]);
		vecRepresentationLoc[i].setZ(vecResTmp[2]);
	}

	std::cout << "On a vecRepresentationLoc[0]: " << vecRepresentationLoc[0] << std::endl;
	std::cout << "On a vecRepresentationLoc[1]: " << vecRepresentationLoc[1] << std::endl;
	std::cout << "On a vecRepresentationLoc[2]: " << vecRepresentationLoc[2] << std::endl;
	double prodScalTests[3];

	for (unsigned int i = 0; i < 3; i++)
		prodScalTests[i] = vecRepresentationLoc[i].dot(ADirNormal);

	std::cout << "Dot products: " << prodScalTests[0] << "; " << prodScalTests[1] << "; " << prodScalTests[2] << std::endl;

	for (unsigned int i = 0; i < 3; i++)
		prodScalTests[i] = AVecRep[i].dot(ADirNormal);

	std::cout << "Init. Dot products: " << prodScalTests[0] << "; " << prodScalTests[1] << "; " << prodScalTests[2] << std::endl;

	// At this point we have the 3 representation vector.
	// Now we need to compute the position of the sing
	// usage des moindres carres
	double coeff[3];
	double matrixA[3][2];
	double matrixb[3];
	double matrixAT[2][3];
	double matrixATA[2][2];
	double matrixATb[2];

	TCoord vecRepresentationTab[3][3];
	for (unsigned int i = 0; i < 3; i++) {
		vecRepresentationTab[i][0] = vecRepresentationLoc[i].X();
		vecRepresentationTab[i][1] = vecRepresentationLoc[i].Y();
		vecRepresentationTab[i][2] = vecRepresentationLoc[i].Z();
	}

	for (unsigned int i = 0; i < 3; i++) {
		for (unsigned int j = 0; j < 2; j++) {
			matrixA[i][j] = vecRepresentationTab[j][i] - vecRepresentationTab[2][i];
			matrixAT[j][i] = matrixA[i][j];
		}
		matrixb[i] = -vecRepresentationTab[2][i];
	}
	for (unsigned int i = 0; i < 2; i++) {
		for (unsigned int j = 0; j < 2; j++) {
			matrixATA[i][j] = 0.0;
			for (unsigned int k = 0; k < 3; k++) {
				matrixATA[i][j] += matrixAT[i][k] * matrixA[k][j];
			}
		}
		matrixATb[i] = 0.0;
		for (unsigned int j = 0; j < 3; j++) {
			matrixATb[i] += matrixAT[i][j] * matrixb[j];
		}
	}
	double det = matrixATA[0][0] * matrixATA[1][1] - matrixATA[0][1] * matrixATA[1][0];
	if (det == 0.0) {
		std::cout << "gros pb det nul" << std::endl;
		throw GMDSException("Determinant null!!!");
	}
	else {
		coeff[0] = matrixATA[1][1] * matrixATb[0] / det - matrixATA[0][1] * matrixATb[1] / det;
		coeff[1] = matrixATA[0][0] * matrixATb[1] / det - matrixATA[1][0] * matrixATb[0] / det;
		coeff[2] = 1.0 - coeff[0] - coeff[1];
	}
	// on a recup dans coeff les coefficient barycentriques de la singularite
	double posSingX = coeff[0] * node1.X() + coeff[1] * node2.X() + coeff[2] * node3.X();
	double posSingY = coeff[0] * node1.Y() + coeff[1] * node2.Y() + coeff[2] * node3.Y();
	double posSingZ = coeff[0] * node1.Z() + coeff[1] * node2.Z() + coeff[2] * node3.Z();

	APosSing.setXYZ(posSingX, posSingY, posSingZ);

	// maintenant il faut recuperer les positions des separatrices
	// on va appeler ComputeSingInfoOnEdge sur chacune des trois edges
	std::cout << "We are in face " << AFace.id() << std::endl;

	//============================================================
	//			FIRST EDGE
	//============================================================
	std::cout << "Info on the first edge" << std::endl;
	computeSingInfoAlongEdgeOfTriangle(AFaceMark, AFace, node1, node2,     // two nodes of AFace
	                                   AVecRep[0], AVecRep[1],             // vectors in previous nodes
	                                   ADirNormal,                         // normal to AFace
	                                   APosSing,                           // field singularity in AFace
	                                   posSep, dirSep, sepInFaceID);

	//============================================================
	//			SECOND EDGE
	//============================================================
	std::cout << "Info on the second edge" << std::endl;
	computeSingInfoAlongEdgeOfTriangle(AFaceMark, AFace, node2, node3,     // two nodes of AFace
	                                   AVecRep[1], AVecRep[2],             // vectors in previous nodes
	                                   ADirNormal,                         // normal to AFace
	                                   APosSing,                           // field singularity in AFace
	                                   posSep, dirSep, sepInFaceID);

	//============================================================
	//			THIRD EDGE
	//============================================================
	std::cout << "Info on the third edge" << std::endl;

	computeSingInfoAlongEdgeOfTriangle(AFaceMark, AFace, node3, node1,     // two nodes of AFace
	                                   AVecRep[2], AVecRep[0],             // vectors in previous nodes
	                                   ADirNormal,                         // normal to AFace
	                                   APosSing,                           // field singularity in AFace
	                                   posSep, dirSep, sepInFaceID);

	// trois aretes faites
	std::cout << "sortie de computeboundSing" << std::endl << std::endl;
}
/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::computeSingInfoAlongEdgeOfTriangle(const int AFaceMark,
                                                            Face &AFace,
                                                            Node &ANode1,
                                                            Node &ANode2,
                                                            math::Vector3d &AV1,
                                                            math::Vector3d &AV2,
                                                            math::Vector3d &ADirNormal,
                                                            const math::Point &AFaceSingularityPoint,
                                                            std::vector<math::Point> &posSep,
                                                            std::vector<math::Vector3d> &dirSep,
                                                            std::vector<int> &sepInFaceID)
{
	double cross1Tmp[2][3];
	double cross2Tmp[2][3];
	double pos1Tmp[3];
	double pos2Tmp[3];
	double posSingTmp[3];
	std::vector<math::Point> posSepTmp;
	std::vector<math::Vector3d> dirSepTmp;
	std::vector<int> sepInFaceIDTmp;

	math::Vector3d orthogForCross;

	orthogForCross = ADirNormal.cross(AV1);
	cross1Tmp[0][0] = AV1.X();
	cross1Tmp[0][1] = AV1.Y();
	cross1Tmp[0][2] = AV1.Z();
	cross1Tmp[1][0] = orthogForCross.X();
	cross1Tmp[1][1] = orthogForCross.Y();
	cross1Tmp[1][2] = orthogForCross.Z();

	orthogForCross = ADirNormal.cross(AV2);
	cross2Tmp[0][0] = AV2.X();
	cross2Tmp[0][1] = AV2.Y();
	cross2Tmp[0][2] = AV2.Z();
	cross2Tmp[1][0] = orthogForCross.X();
	cross2Tmp[1][1] = orthogForCross.Y();
	cross2Tmp[1][2] = orthogForCross.Z();

	pos1Tmp[0] = ANode1.X();
	pos1Tmp[1] = ANode1.Y();
	pos1Tmp[2] = ANode1.Z();
	pos2Tmp[0] = ANode2.X();
	pos2Tmp[1] = ANode2.Y();
	pos2Tmp[2] = ANode2.Z();

	posSingTmp[0] = AFaceSingularityPoint.X();
	posSingTmp[1] = AFaceSingularityPoint.Y();
	posSingTmp[2] = AFaceSingularityPoint.Z();

	posSepTmp.clear();
	dirSepTmp.clear();
	sepInFaceIDTmp.clear();

	computeSingInfoOnEdge(AFaceMark, ANode1, ANode2, AFace, cross1Tmp, cross2Tmp, ADirNormal, pos1Tmp, pos2Tmp, posSingTmp, posSepTmp, dirSepTmp,
	                      sepInFaceIDTmp);

	std::cout << "For an edge, we have size of posSepTmp: " << posSepTmp.size() << "; and dirSepTmp: " << dirSepTmp.size() << "; "
	          << ", and sepInFaceIDTmp: " << sepInFaceIDTmp.size() << std::endl;

	for (unsigned int i = 0; i < posSepTmp.size(); i++) {
		std::cout << "For an edge, we insert at location " << posSepTmp[i] << " the beginning of a sep. in direction " << dirSepTmp[i]
		          << " starting from triangle " << sepInFaceIDTmp[i] << std::endl;
		posSep.push_back(posSepTmp[i]);
		dirSep.push_back(dirSepTmp[i]);
		sepInFaceID.push_back(sepInFaceIDTmp[i]);
	}
}
/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::computeSingInfoOnEdge(const int AFaceMark,
                                               Node &ANode1,
                                               Node &ANode2,
                                               Face &AFace,
                                               double cross1[2][3],
                                               double cross2[2][3],
                                               math::Vector3d &ANormalToFace,
                                               double pos1[3],
                                               double pos2[3],
                                               double posSing[3],
                                               std::vector<math::Point> &posSep,
                                               std::vector<math::Vector3d> &dirSep,
                                               std::vector<int> &sepInFaceID)
{
	// Pour commencer on calcule la matrice de rotation M amenant la normale en 0,0,1
	// ensuite on applique cette matrice aux deux croix et aux points
	// puis on resoud les equations comme en 2D
	// et on applique linverse de M aux resultats
	// CALCUL de la matrice de rotation
	math::Vector3d orthogBase({ANormalToFace.Y(), -ANormalToFace.X(), 0});     // le vecteur orthogonal a normaleTriang et 0,0,1

	double matrixRotation[3][3];
	double normeOrthog = orthogBase.norm();

	if (normeOrthog == 0.0) {
		for (unsigned int i = 0; i < 3; i++) {
			for (unsigned int j = 0; j < 3; j++) {
				if (i == j) {
					matrixRotation[i][j] = 1.0;
				}
				else {
					matrixRotation[i][j] = 0.0;
				}
			}
		}
	}
	else {
		orthogBase.normalize();

		// matrix non identite
		double prodScalTmp = ANormalToFace.Z();
		double angleTmp = acos(prodScalTmp);
		while (angleTmp >= 4 * acos(0.0)) {
			angleTmp -= 4 * acos(0.0);
		}
		while (angleTmp < 0.0) {
			angleTmp += 4 * acos(0.0);
		}
		double cTmp = cos(angleTmp);
		double sTmp = sin(angleTmp);
		matrixRotation[0][0] = orthogBase.X() * orthogBase.X() + (1 - orthogBase.X() * orthogBase.X()) * cTmp;
		matrixRotation[0][1] = orthogBase.X() * orthogBase.Y() * (1 - cTmp) - orthogBase.Z() * sTmp;
		matrixRotation[0][2] = orthogBase.X() * orthogBase.Z() * (1 - cTmp) + orthogBase.Y() * sTmp;
		matrixRotation[1][0] = orthogBase.X() * orthogBase.Y() * (1 - cTmp) + orthogBase.Z() * sTmp;
		matrixRotation[1][1] = orthogBase.Y() * orthogBase.Y() + (1 - orthogBase.Y() * orthogBase.Y()) * cTmp;
		matrixRotation[1][2] = orthogBase.Z() * orthogBase.Y() * (1 - cTmp) - orthogBase.X() * sTmp;
		matrixRotation[2][0] = orthogBase.X() * orthogBase.Z() * (1 - cTmp) - orthogBase.Y() * sTmp;
		matrixRotation[2][1] = orthogBase.Z() * orthogBase.Y() * (1 - cTmp) + orthogBase.X() * sTmp;
		matrixRotation[2][2] = orthogBase.Z() * orthogBase.Z() + (1 - orthogBase.Z() * orthogBase.Z()) * cTmp;
	}
	double cross1Rotated[2][3];
	double cross2Rotated[2][3];
	double pos1Rotated[3];
	double pos2Rotated[3];
	double posSingRotated[3];
	std::vector<double> posSepXRotated;
	std::vector<double> posSepYRotated;
	std::vector<double> dirSepXRotated;
	std::vector<double> dirSepYRotated;
	for (unsigned int i = 0; i < 3; i++) {
		cross1Rotated[0][i] = 0.0;
		cross1Rotated[1][i] = 0.0;
		cross2Rotated[0][i] = 0.0;
		cross2Rotated[1][i] = 0.0;
		pos1Rotated[i] = 0.0;
		pos2Rotated[i] = 0.0;
		posSingRotated[i] = 0.0;
		for (unsigned int j = 0; j < 3; j++) {
			cross1Rotated[0][i] += matrixRotation[i][j] * cross1[0][j];
			cross1Rotated[1][i] += matrixRotation[i][j] * cross1[1][j];
			cross2Rotated[0][i] += matrixRotation[i][j] * cross2[0][j];
			cross2Rotated[1][i] += matrixRotation[i][j] * cross2[1][j];
			pos1Rotated[i] += matrixRotation[i][j] * pos1[j];
			pos2Rotated[i] += matrixRotation[i][j] * pos2[j];
			posSingRotated[i] += matrixRotation[i][j] * posSing[j];
		}
	}

	double xs = posSing[0];
	double ys = posSing[1];
	double x1 = pos1[0];
	double y1 = pos1[1];
	double z1 = pos1[2];
	double x2 = pos2[0];
	double y2 = pos2[1];
	double z2 = pos2[2];

	double x2MINx1 = x2 - x1;
	double y2MINy1 = y2 - y1;
	double z2MINz1 = pos2[2] - pos1[2];
	double x1MINxs = x1 - xs;
	double y1MINys = y1 - ys;
	double z1MINzs = pos1[2] - posSing[2];
	double x2MINx1Rotated = 0.0;
	double y2MINy1Rotated = 0.0;
	double x1MINxsRotated = 0.0;
	double y1MINysRotated = 0.0;

	x2MINx1Rotated += matrixRotation[0][0] * x2MINx1;
	x2MINx1Rotated += matrixRotation[0][1] * y2MINy1;
	x2MINx1Rotated += matrixRotation[0][2] * z2MINz1;
	y2MINy1Rotated += matrixRotation[1][0] * x2MINx1;
	y2MINy1Rotated += matrixRotation[1][1] * y2MINy1;
	y2MINy1Rotated += matrixRotation[1][2] * z2MINz1;
	x1MINxsRotated += matrixRotation[0][0] * x1MINxs;
	x1MINxsRotated += matrixRotation[0][1] * y1MINys;
	x1MINxsRotated += matrixRotation[0][2] * z1MINzs;
	y1MINysRotated += matrixRotation[1][0] * x1MINxs;
	y1MINysRotated += matrixRotation[1][1] * y1MINys;
	y1MINysRotated += matrixRotation[1][2] * z1MINzs;

	//=============================================================
	// first the first vector of the 4 RoSy
	//=============================================================
	double v1 = cross1Rotated[0][0];
	double w1 = cross1Rotated[0][1];
	double v2 = cross2Rotated[0][0];
	double w2 = cross2Rotated[0][1];
	double v2test = cross2Rotated[1][0];
	double w2test = cross2Rotated[1][1];
	if ((v1 * v2 + w1 * w2) < (v1 * v2test + w1 * w2test)) {
		v2 = v2test;
		w2 = w2test;
	}
	v2test = -cross2Rotated[0][0];
	w2test = -cross2Rotated[0][1];
	if ((v1 * v2 + w1 * w2) < (v1 * v2test + w1 * w2test)) {
		v2 = v2test;
		w2 = w2test;
	}
	v2test = -cross2Rotated[1][0];
	w2test = -cross2Rotated[1][1];
	if ((v1 * v2 + w1 * w2) < (v1 * v2test + w1 * w2test)) {
		v2 = v2test;
		w2 = w2test;
	}

	double a = (w2 - w1) * x2MINx1Rotated - (v2 - v1) * y2MINy1Rotated;
	double b = (w2 - w1) * x1MINxsRotated + x2MINx1Rotated * w1 - ((v2 - v1) * y1MINysRotated + y2MINy1Rotated * v1);
	double c = w1 * x1MINxsRotated - v1 * y1MINysRotated;

	solveSecondOrderEq(a, b, c, AFaceMark, AFace, ANode1, ANode2, x1, y1, z1, x2, y2, v1, w1, v2, w2, x2MINx1, y2MINy1, z2MINz1, x1MINxsRotated, x2MINx1Rotated,
	                   y1MINysRotated, y2MINy1Rotated, posSep, dirSepXRotated, dirSepYRotated, sepInFaceID);

	//=============================================================
	// second vector
	//=============================================================

	v1 = cross1Rotated[1][0];
	w1 = cross1Rotated[1][1];
	v2 = cross2Rotated[0][0];
	w2 = cross2Rotated[0][1];
	v2test = cross2Rotated[1][0];
	w2test = cross2Rotated[1][1];
	if ((v1 * v2 + w1 * w2) < (v1 * v2test + w1 * w2test)) {
		v2 = v2test;
		w2 = w2test;
	}
	v2test = -cross2Rotated[0][0];
	w2test = -cross2Rotated[0][1];
	if ((v1 * v2 + w1 * w2) < (v1 * v2test + w1 * w2test)) {
		v2 = v2test;
		w2 = w2test;
	}
	v2test = -cross2Rotated[1][0];
	w2test = -cross2Rotated[1][1];
	if ((v1 * v2 + w1 * w2) < (v1 * v2test + w1 * w2test)) {
		v2 = v2test;
		w2 = w2test;
	}

	a = (w2 - w1) * x2MINx1Rotated - (v2 - v1) * y2MINy1Rotated;
	b = (w2 - w1) * x1MINxsRotated + x2MINx1Rotated * w1;
	b = b - ((v2 - v1) * y1MINysRotated + y2MINy1Rotated * v1);
	c = w1 * x1MINxsRotated - v1 * y1MINysRotated;
	solveSecondOrderEq(a, b, c, AFaceMark, AFace, ANode1, ANode2, x1, y1, z1, x2, y2, v1, w1, v2, w2, x2MINx1, y2MINy1, z2MINz1, x1MINxsRotated, x2MINx1Rotated,
	                   y1MINysRotated, y2MINy1Rotated, posSep, dirSepXRotated, dirSepYRotated, sepInFaceID);

	// premiere arete deuxieme vecteur fait
	// premiere arete finie
	// maintenant on a besoin de multiplier par
	// linverse de la matrice de rotation
	// qui est sa transposee
	for (unsigned int i = 0; i < dirSepXRotated.size(); i++) {
		double dirSepTmp[3];
		dirSepTmp[0] = dirSepXRotated[i];
		dirSepTmp[1] = dirSepYRotated[i];
		dirSepTmp[2] = 0.0;

		double dirSepRotTmp[3];
		for (unsigned int j = 0; j < 3; j++) {

			dirSepRotTmp[j] = 0.0;
			for (unsigned int k = 0; k < 3; k++) {

				dirSepRotTmp[j] += matrixRotation[k][j] * dirSepTmp[k];
			}
		}
		math::Vector3d v({dirSepRotTmp[0], dirSepRotTmp[1], dirSepRotTmp[2]});
		dirSep.push_back(v);
	}
}

/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::solveSecondOrderEq(const double a,
                                            const double b,
                                            const double c,
                                            const int &AFaceMark,
                                            Face &AFace,
                                            Node &ANode1,
                                            Node &ANode2,
                                            const double x1,
                                            const double y1,
                                            const double z1,
                                            const double x2,
                                            const double y2,
                                            const double v1,
                                            const double w1,
                                            const double v2,
                                            const double w2,
                                            const double x2MINx1,
                                            const double y2MINy1,
                                            const double z2MINz1,
                                            const double x1MINxsRotated,
                                            const double x2MINx1Rotated,
                                            const double y1MINysRotated,
                                            const double y2MINy1Rotated,
                                            std::vector<math::Point> &posSep,
                                            std::vector<double> &dirSepXRotated,
                                            std::vector<double> &dirSepYRotated,
                                            std::vector<int> &sepInFaceID)
{

	double soleq = 0.0;

	if (a == 0.0) {
		if (b != 0.0) {
			soleq = -c / b;
			if ((soleq >= 0.0) && (soleq <= 1.0)) {
				buildInfoAlongEdge(AFaceMark, AFace, ANode1, ANode2, soleq, x1, y1, z1, x2, y2, v1, w1, v2, w2, x2MINx1, y2MINy1, z2MINz1, x1MINxsRotated,
				                   x2MINx1Rotated, y1MINysRotated, y2MINy1Rotated, posSep, dirSepXRotated, dirSepYRotated, sepInFaceID);
			}
		}
		else {
			std::cout << "Erreur sur equation pas de solution" << std::endl;
			//	throw GMDSException("Erreur sur equation pas de solution");
		}
	}
	else     // a!=0
	{
		// traditional case, second-order equation
		double det = b * b - 4 * a * c;
		if (det < 0.0) {
			std::cout << "erreur sur eq pas de solution car det neg" << std::endl;
			//	throw GMDSException("erreur sur eq pas de solution car det neg");
		}
		else     // det >=0
		{
			soleq = (-b - sqrt(det)) / (2 * a);
			if ((soleq >= 0.0) && (soleq <= 1.0)) {
				buildInfoAlongEdge(AFaceMark, AFace, ANode1, ANode2, soleq, x1, y1, z1, x2, y2, v1, w1, v2, w2, x2MINx1, y2MINy1, z2MINz1, x1MINxsRotated,
				                   x2MINx1Rotated, y1MINysRotated, y2MINy1Rotated, posSep, dirSepXRotated, dirSepYRotated, sepInFaceID);
			}
			if (det != 0.0)     // a second solution can be computed
			{
				soleq = (-b + sqrt(det)) / (2 * a);
				if ((soleq >= 0.0) && (soleq <= 1.0)) {
					buildInfoAlongEdge(AFaceMark, AFace, ANode1, ANode2, soleq, x1, y1, z1, x2, y2, v1, w1, v2, w2, x2MINx1, y2MINy1, z2MINz1, x1MINxsRotated,
					                   x2MINx1Rotated, y1MINysRotated, y2MINy1Rotated, posSep, dirSepXRotated, dirSepYRotated, sepInFaceID);
				}
			}
		}
	}
}
/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::buildInfoAlongEdge(const int &AFaceMark,
                                            Face &AFace,
                                            Node &ANode1,
                                            Node &ANode2,
                                            const double soleq,
                                            const double x1,
                                            const double y1,
                                            const double z1,
                                            const double x2,
                                            const double y2,
                                            const double v1,
                                            const double w1,
                                            const double v2,
                                            const double w2,
                                            const double x2MINx1,
                                            const double y2MINy1,
                                            const double z2MINz1,
                                            const double x1MINxsRotated,
                                            const double x2MINx1Rotated,
                                            const double y1MINysRotated,
                                            const double y2MINy1Rotated,
                                            std::vector<math::Point> &posSep,
                                            std::vector<double> &dirSepXRotated,
                                            std::vector<double> &dirSepYRotated,
                                            std::vector<int> &sepInFaceID)
{
	double xFound = x1 + soleq * x2MINx1;
	double yFound = y1 + soleq * y2MINy1;
	double zFound = z1 + soleq * z2MINz1;
	double vFound = v1 + soleq * (v2 - v1);
	double wFound = w1 + soleq * (w2 - w1);

	if (((x1MINxsRotated + soleq * x2MINx1Rotated) * vFound + (y1MINysRotated + soleq * y2MINy1Rotated) * wFound) < 0) {
		// vecteur dans le mauvais sens, vers singularite
		vFound = -vFound;
		wFound = -wFound;
	}
	double norme = sqrt(vFound * vFound + wFound * wFound);
	if (norme != 0.0) {
		vFound = vFound / norme;
		wFound = wFound / norme;
	}
	else {
		vFound = 1.0;
		wFound = 0.0;
	}
	Node Node1Tmp = ANode1;
	Node Node2Tmp = ANode2;
	for (int i = 0; i < Node1Tmp.nbFaces(); i++) {
		for (int j = 0; j < Node2Tmp.nbFaces(); j++) {
			if ((Node1Tmp.get<Face>()[i].id() == Node2Tmp.get<Face>()[j].id()) && (Node1Tmp.get<Face>()[i].id() != AFace.id())
			    && (m_mesh->isMarked(Node1Tmp.get<Face>()[i], AFaceMark))) {
				sepInFaceID.push_back(Node1Tmp.get<Face>()[i].id());
			}
		}
	}
	posSep.push_back(math::Point(xFound, yFound, zFound));
	dirSepXRotated.push_back(vFound);
	dirSepYRotated.push_back(wFound);
}
/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::createVolumeSingularityLines()
{

	// at this point we have made all cluster of sing into Singularity3D
	// now we need to compute the lines of sing into separatrices
	int markFaceUsedForSep = m_mesh->newMark<Face>();

	/* We go through all the volume singularity points and we create a sing line
	 * for every available free slot.
	 */
	//	std::vector<VolumeSingularityPoint*> singularity_points = m_graph.olumePoints();
	std::vector<SingularityPoint *> singularity_points = m_graph.getPoints();
	for (unsigned int i = 0; i < singularity_points.size(); i++) {
		std::cout << "POINT " << i << std::endl;
		SingularityPoint *pnt_i = singularity_points[i];

		// We have a volume singularity point
		std::vector<SingularityPoint::Slot *> slots_i = pnt_i->getSlots();
		std::cout << "==> We create sing line from vol. sing. pnt " << pnt_i->getLocation() << ", having " << slots_i.size() << std::endl;
		for (unsigned int j = 0; j < slots_i.size(); j++) {
			std::cout << " --> slot " << j << std::endl;
			SingularityPoint::Slot *s_j = slots_i[j];
			if (s_j->isOnSurface)
				std::cout << "on surface" << std::endl;
			else {
				std::cout << "  <got> " << j << " --> sj adress " << s_j << std::endl;
				if (!s_j->isLaunched) {
					std::cout << " not launched" << std::endl;
					createOneSingularityLineFrom(pnt_i, s_j, markFaceUsedForSep);
				}
			}
		}
		writeOutput("lines");
	}

	//=========================================================================
	// SMOOTHING OF SEPARATRICES
	//=========================================================================
	std::vector<SingularityLine *> sing_lines = m_graph.getLines();
	std::cout << "Nb sing lines = " << sing_lines.size() << std::endl;
	for (unsigned int i = 0; i < sing_lines.size(); i++) {
		SingularityLine *line_i = sing_lines[i];
		std::cout << i << ") " << line_i << std::endl;
		line_i->smooth();
	}

	m_mesh->unmarkAll<Face>(markFaceUsedForSep);
	m_mesh->freeMark<Face>(markFaceUsedForSep);
}
/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::initSingularityPointOnSurf()
{
	//=========================================================================
	// First, we initialize the m_face_to_singularity_on_surf value to NULL
	// for every boundary face
	//=========================================================================

	for (auto f_id : m_mesh->faces()) {
		// while (!it_face.isDone()){
		Face f = m_mesh->get<Face>(f_id);
		if (m_mesh->isMarked(f, m_markFacesOnSurf)) m_face_to_singularity_on_surf[f.id()] = NULL;
	}
	//=========================================================================
	// For each boudary region, we create a singularity point with
	// slots oriented inside the volume and along the surface.
	//=========================================================================
	std::list<TCellID>::iterator it = m_2SingTetIDsAlongSurf.begin();
	for (; it != m_2SingTetIDsAlongSurf.end(); it++) {
		Region current = m_mesh->get<Region>(*it);
		std::vector<Face> current_faces = current.get<Face>();
		std::vector<Face> singular_bnd_faces;
		// We have to detect which face is traversed by the sing line.
		// Indeed it can be some faces of current that are not boundary faces.
		// As a consequence, we have no surface sing. point to introduce
		// in this case.
		for (unsigned int j = 0; j < current_faces.size(); j++) {
			Face face_j = current_faces[j];
			if (isSingularFace(face_j) && m_mesh->isMarked(face_j, m_markFacesOnSurf)) {
				singular_bnd_faces.push_back(face_j);
			}
		}

		if (singular_bnd_faces.size() == 1) {
			// We have the face the singularity comes through
			createSurfaceSingularityPointFromVolumeInfo(singular_bnd_faces[0], false);
		}
		else if (singular_bnd_faces.size() > 1) {
			// Not yet handle
			throw GMDSException("Problem to define a frame singulary point on the surface: too many mesh faces are candidates!");
		}
		// Otherwise no, traversed boundary face, nothing to do.
	}
}
/*----------------------------------------------------------------------------*/
bool
SingularityGraphBuilder::isSingularFace(Face &AF)
{
	std::vector<TCellID> node_ids = AF.getIDs<Node>();
	int ID1 = node_ids[0];
	int ID2 = node_ids[1];
	int ID3 = node_ids[2];
	int ID4 = node_ids[2];
	math::Quaternion q1 = (*m_var_quatern)[ID1];
	math::Quaternion q2 = (*m_var_quatern)[ID2];
	math::Quaternion q3 = (*m_var_quatern)[ID3];
	math::Quaternion q4 = (*m_var_quatern)[ID4];
	int singularity_type = math::Quaternion::testSingularity(q1, q2, q3, q4);

	return (singularity_type == 2);
}
/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::createVolumeSingularityPoints()
{

	int nbOfCluster = 0;

	//=========================================================================
	// INTERN SKELETON CREATION
	//=========================================================================
	// FIRST LOOP ON REGIONS TO GET ALL THE 3-SING. TETS
	//=========================================================================
	std::cout << "Start detection" << std::endl;
	for (auto r_id : m_mesh->regions()) {
		Region Rtmp = m_mesh->get<Region>(r_id);

		int singTypeTmp = m_region_singularity_type[Rtmp.id()];
		if (singTypeTmp == 3 && (!m_mesh->isMarked(Rtmp, m_markClusterSingDone))) {
			nbOfCluster++;
			// we are currently in a region that is part of a cluster and has
			// not been used yet, we are going to use it now

			m_mesh->mark(Rtmp, m_markClusterSingDone);
			std::vector<Region> vecRegions;

			vecRegions.push_back(Rtmp);
			int stopAlready = 0;
			while (!stopAlready) {
				stopAlready = 1;
				for (unsigned int i = 0; i < vecRegions.size(); i++) {
					for (int j = 0; j < vecRegions[i].nbFaces(); j++) {
						Region ROpposite;
						if (vecRegions[i].get<Face>()[j].get<Region>()[0].id() == vecRegions[i].id()) {
							ROpposite = vecRegions[i].get<Face>()[j].get<Region>()[1];
						}
						else {
							ROpposite = vecRegions[i].get<Face>()[j].get<Region>()[0];
						}

						int singTypeTmpOpp = m_region_singularity_type[ROpposite.id()];
						if (singTypeTmpOpp == 3) {
							if (!m_mesh->isMarked(ROpposite, m_markClusterSingDone)) {
								vecRegions.push_back(ROpposite);
								m_mesh->mark(ROpposite, m_markClusterSingDone);
								stopAlready = 0;
							}
						}
					}
				}
			}     // while (!stopAlready)
			//==================================================================
			// We create the volume singularity point
			//==================================================================
			createOneVolumeSingularityPoint(vecRegions);
		}     // if (singTypeTmp == 3 && (!m_mesh->isMarked(Rtmp, m_markClusterSingDone)))
	}

	//=========================================================================
	std::cout << "We have " << nbOfCluster << " clusters of 3-sing. tetrahedra, ";
	std::cout << "and so " << m_graph.getNbVolumePoints() << " sing. graph points" << std::endl;
}

/*----------------------------------------------------------------------------*/

void
SingularityGraphBuilder::createOneVolumeSingularityPoint(std::vector<gmds::Region> &ACluster)

// m_markClusterSingDone, m_markBordVolSingForFaces, m_markBordVolSingForEdges
// const int AMarkCluster,
// const int AMarkBordSingFace, const int AMarkBordSingEdge)
{
	VolumeSingularityPoint *singularity = m_graph.newVolumePoint();
	//==================================================================
	// First, we compute the center of mass of the cluster
	//==================================================================
	int nbTetInClust = ACluster.size();
	std::cout << "In the current cluster we have " << nbTetInClust << " tet" << std::endl;
	double posXTmp = 0.0;
	double posYTmp = 0.0;
	double posZTmp = 0.0;
	for (int i = 0; i < nbTetInClust; i++) {
		for (unsigned int j = 0; j < 4; j++) {
			posXTmp += ACluster[i].get<Node>()[j].X();
			posYTmp += ACluster[i].get<Node>()[j].Y();
			posZTmp += ACluster[i].get<Node>()[j].Z();
		}
	}
	posXTmp = posXTmp / (4.0 * ((double) nbTetInClust));
	posYTmp = posYTmp / (4.0 * ((double) nbTetInClust));
	posZTmp = posZTmp / (4.0 * ((double) nbTetInClust));
	//==================================================================
	// location
	math::Point sing_pnt(posXTmp, posYTmp, posZTmp);
	singularity->setLocation(sing_pnt);
	//==================================================================
	// mesh association
	for (unsigned int i = 0; i < ACluster.size(); i++)
		singularity->addMeshRegion(ACluster[i]);

	//==================================================================
	// slots definition

	// We get the outer faces of the cluster
	std::vector<Face> cluster_bnd_faces;

	for (int i = 0; i < nbTetInClust; i++) {
		std::vector<Face> local_faces = ACluster[i].get<Face>();
		for (unsigned int j = 0; j < local_faces.size(); j++) {
			Face current_face = local_faces[j];
			std::vector<Region> local_regions = current_face.get<Region>();
			if (local_regions.size() == 1)
				cluster_bnd_faces.push_back(current_face);
			else {
				Region r0 = local_regions[0];
				Region r1 = local_regions[1];
				if ((m_mesh->isMarked(r0, m_markClusterSingDone) && !m_mesh->isMarked(r1, m_markClusterSingDone))
				    || (!m_mesh->isMarked(r0, m_markClusterSingDone) && m_mesh->isMarked(r1, m_markClusterSingDone)))
					cluster_bnd_faces.push_back(current_face);
			}
		}
	}
	for (unsigned int i = 0; i < cluster_bnd_faces.size(); i++) {
		Face bnd_face = cluster_bnd_faces[i];

		m_mesh->mark(bnd_face, m_markBordVolSingForFaces);
		std::vector<Edge> current_edges = bnd_face.get<Edge>();
		for (unsigned int i_edges = 0; i_edges < current_edges.size(); i_edges++)
			m_mesh->mark(current_edges[i_edges], m_markBordVolSingForEdges);

		std::vector<Node> current_nodes = bnd_face.get<Node>();
		for (unsigned int i_nodes = 0; i_nodes < current_nodes.size(); i_nodes++)
			m_mesh->mark(current_nodes[i_nodes], m_markBordVolSingForNodes);

		if (isSingularFace(bnd_face)) {
			// A Singularty line goes through bnd_face
			math::Point face_center = bnd_face.center();
			math::Vector3d line_vec= face_center-sing_pnt;
			// TODO CHANGE singularity->newSlot(face_center, line_vec, bnd_face, false);
		}
	}

	std::cout << "Volume sing. point created at " << singularity->getLocation() << " in a cluster of " << singularity->getNbMeshCells() << " with "
	          << singularity->getSlots().size() << " slots " << std::endl;
}
/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::InvertMatrix(double matrixInit[3][3],
                                      double &matrixInv00,
                                      double &matrixInv01,
                                      double &matrixInv02,
                                      double &matrixInv10,
                                      double &matrixInv11,
                                      double &matrixInv12,
                                      double &matrixInv20,
                                      double &matrixInv21,
                                      double &matrixInv22)
{
	double a = matrixInit[0][0];
	double b = matrixInit[0][1];
	double c = matrixInit[0][2];
	double d = matrixInit[1][0];
	double e = matrixInit[1][1];
	double f = matrixInit[1][2];
	double g = matrixInit[2][0];
	double h = matrixInit[2][1];
	double i = matrixInit[2][2];
	double detM = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
	matrixInv00 = (e * i - f * h) / detM;
	matrixInv01 = (c * h - b * i) / detM;
	matrixInv02 = (b * f - c * e) / detM;
	matrixInv10 = (f * g - d * i) / detM;
	matrixInv11 = (a * i - c * g) / detM;
	matrixInv12 = (c * d - a * f) / detM;
	matrixInv20 = (d * h - e * g) / detM;
	matrixInv21 = (b * g - a * h) / detM;
	matrixInv22 = (a * e - b * d) / detM;
}
/*----------------------------------------------------------------------------*/
bool
SingularityGraphBuilder::intersecDirectionWithSegment(double P1[3],
                                                      double P2[3],
                                                      double PStart[3],
                                                      double dirStart[3],
                                                      double &coeff,
                                                      double &coeffBis,
                                                      double &PEnd0,
                                                      double &PEnd1,
                                                      double &PEnd2,
                                                      double epsilon)
{
	// PRINCIPE: on verifie si meme direction, dans ce cas pas intersection
	// sinon on prend prod vectoriel des deux directions
	// donne une direction normale
	// puis on construit matrice trois par trois du pb
	// std::cout << std::endl;
	// std::cout << "Get into intersecDirectionWithSegment" << std::endl;
	// std::cout << "with P1: " << P1[0] << "; " << P1[1] << "; " << P1[2] << std::endl;
	// std::cout << "     P2: " << P2[0] << "; " << P2[1] << "; " << P2[2] << std::endl;
	// std::cout << "     PStart: " << PStart[0] << "; " << PStart[1] << "; " << PStart[2] << std::endl;
	// std::cout << "     dirStart: " << dirStart[0] << "; " << dirStart[1] << "; " << dirStart[2] << std::endl;
	double dirPoints[3];
	dirPoints[0] = P2[0] - P1[0];
	dirPoints[1] = P2[1] - P1[1];
	dirPoints[2] = P2[2] - P1[2];
	// std::cout << "on a dirPoints: " << dirPoints[0] << "; " << dirPoints[1] << "; " << dirPoints[2] << std::endl;
	double normeTmp = sqrt(dirPoints[0] * dirPoints[0] + dirPoints[1] * dirPoints[1] + dirPoints[2] * dirPoints[2]);
	if (normeTmp > 0.0000001) {
		//		for (unsigned int i = 0;i < 3;i++){
		//			dirPoints[i] = dirPoints[i] / normeTmp;
		//		}
		double dirOrthog[3];
		dirOrthog[0] = dirStart[1] * dirPoints[2] - dirStart[2] * dirPoints[1];
		dirOrthog[1] = dirStart[2] * dirPoints[0] - dirStart[0] * dirPoints[2];
		dirOrthog[2] = dirStart[0] * dirPoints[1] - dirStart[1] * dirPoints[0];
		//	std::cout<<"on a dirOrthog: "<<dirOrthog[0]<<"; "<<dirOrthog[1]<<"; "<<dirOrthog[2]<<std::endl;
		normeTmp = sqrt(dirOrthog[0] * dirOrthog[0] + dirOrthog[1] * dirOrthog[1] + dirOrthog[2] * dirOrthog[2]);
		for (unsigned int i = 0; i < 3; i++) {
			dirOrthog[i] = dirOrthog[i] / normeTmp;
		}
		// on a les trois vecteurs qui forment une base, on va creer la matrice du pb et linverser
		double matrixInit[3][3];
		for (unsigned int i = 0; i < 3; i++) {
			matrixInit[i][0] = dirStart[i];
			matrixInit[i][1] = dirPoints[i];
			matrixInit[i][2] = dirOrthog[i];
		}
		double matrixInv[3][3];
		double matrixInv00;
		double matrixInv01;
		double matrixInv02;
		double matrixInv10;
		double matrixInv11;
		double matrixInv12;
		double matrixInv20;
		double matrixInv21;
		double matrixInv22;
		InvertMatrix(matrixInit, matrixInv00, matrixInv01, matrixInv02, matrixInv10, matrixInv11, matrixInv12, matrixInv20, matrixInv21, matrixInv22);
		matrixInv[0][0] = matrixInv00;
		matrixInv[0][1] = matrixInv01;
		matrixInv[0][2] = matrixInv02;
		matrixInv[1][0] = matrixInv10;
		matrixInv[1][1] = matrixInv11;
		matrixInv[1][2] = matrixInv12;
		matrixInv[2][0] = matrixInv20;
		matrixInv[2][1] = matrixInv21;
		matrixInv[2][2] = matrixInv22;
		for (unsigned int i = 0; i < 3; i++) {
			for (unsigned int j = 0; j < 3; j++) {
				//				std::cout<<"pour i: "<<i;
				//				std::cout<<" et pour j: "<<j;
				//				std::cout<<std::endl<<"on a matrixInit["<<i<<"]["<<j<<"] = "<<matrixInit[i][j]<<std::endl;
				//				std::cout<<"on a matrixInv["<<i<<"]["<<j<<"] = "<<matrixInv[i][j]<<std::endl;
				double sumTmp = 0.0;
				double sumTmpBis = 0.0;
				for (unsigned int k = 0; k < 3; k++) {
					sumTmp += matrixInit[i][k] * matrixInv[k][j];
					sumTmpBis += matrixInv[i][k] * matrixInit[k][j];
				}
				//			std::cout<<"on a une somme de "<<sumTmp<<" et somme bis de "<<sumTmpBis<<std::endl;;
			}
		}
		double memberB[3];
		for (unsigned int i = 0; i < 3; i++) {
			memberB[i] = P1[i] - PStart[i];
			//		std::cout<<"on a memberB de "<<i<<" = "<<memberB[i]<<std::endl;
		}
		double solut[3];
		for (unsigned int i = 0; i < 3; i++) {
			solut[i] = 0.0;
			for (unsigned int j = 0; j < 3; j++) {
				solut[i] += matrixInv[i][j] * memberB[j];
			}
		}
		// std::cout << "on a des solutions: " << solut[0] << "; " << solut[1] << "; " << solut[2] << std::endl;
		coeff = -solut[1];
		coeffBis = solut[2] * solut[2];
		PEnd0 = P1[0] - solut[1] * dirPoints[0];
		PEnd1 = P1[1] - solut[1] * dirPoints[1];
		PEnd2 = P1[2] - solut[1] * dirPoints[2];
		//		for (unsigned int i = 0;i < 3;i++){
		//			PEnd[i] = P1[i] - solut[1] * dirPoints[i];
		//		}
		if (solut[1] < 0.0) {
			if (solut[1] > (-1.0)) {
				if (solut[0] > 0.000001) {
					return true;
				}
				else {
					return false;
				}
			}
			else {
				return false;
			}
		}
		else {
			return false;
		}
	}
	else {
		return false;
	}

	bool result = false;

	return result;
}

/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::addGeometryToSingularitGraph()
{
	std::cout << "== Geometry Incorporation ==" << std::endl;
	// Now we add singularty points and lines for the corner and edges of the geometry

	std::cout << "NB inner singularities = " << m_graph.getNbPoints() << std::endl;

	for (auto n_id : m_mesh->nodes()) {
		Node current_node = m_mesh->get<Node>(n_id);
		if (m_mesh->isMarked(current_node, m_markNodesOnVert)) {
			VertexSingularityPoint *sing_point = m_graph.newVertexPoint();
			sing_point->setXYZ(current_node.X(), current_node.Y(), current_node.Z());
			sing_point->addMeshNode(current_node);

			// Volume singularity lines must be added for non-convex surface points
			std::vector<Edge> current_edges = current_node.get<Edge>();
			std::vector<Face> current_faces = current_node.get<Face>();
			for (unsigned int i = 0; i < current_edges.size(); i++) {
				Edge currentE = current_edges[i];

				// only edges classified on geometrical curves are taken into account
				if (!m_mesh->isMarked(currentE, m_markEdgesOnCurv)) continue;

				std::vector<Node> edge_nodes = currentE.get<Node>();
				Node other_node;
				if (edge_nodes[0] == current_node)
					other_node = edge_nodes[1];
				else
					other_node = edge_nodes[0];

				math::Vector3d vec= other_node.point()-current_node.point();
				vec.normalize();

				// Maintenant on va regarder pour toutes les faces adjacentes
				// a current_node si une d'elles est intersectee par vec.
				bool found = false;     // indique si on a trouve une intersection
				for (unsigned int j = 0; j < current_faces.size() && !found; j++) {
					Face currentF = current_faces[j];
					if (!m_mesh->isMarked(currentF, m_markFacesOnSurf)) continue;     // on ne traite que les faces au bord

					// on commence par regarder si on est dans une face succeptible d'etre
					// intersectee par la courbe a prolonger
					math::Vector3d fNormale = currentF.normal();
					// dans le cas qui suit on est quasi orthogonal a la face
					double dotProduct = vec.dot(fNormale);

					if (fabs(dotProduct) > 0.6) continue;

					// on projete le vecteur vec sur la face
					math::Vector3d proj;
					if (dotProduct > 0)
						proj = vec - dotProduct * fNormale;
					else
						proj = vec + dotProduct * fNormale;

					// on recupere les noeuds de l'arete opposee dans la face pour
					// tester si on intersecte
					std::vector<Node> current_nodes = currentF.get<Node>();
					Node oppositeNodes[2];
					if (current_nodes[0] == current_node) {
						oppositeNodes[0] = current_nodes[1];
						oppositeNodes[1] = current_nodes[2];
					}
					else if (current_nodes[1] == current_node) {
						oppositeNodes[0] = current_nodes[0];
						oppositeNodes[1] = current_nodes[2];
					}
					else {
						oppositeNodes[0] = current_nodes[0];
						oppositeNodes[1] = current_nodes[1];
					}

					double P1[3], P2[3], PStart[3], dirStart[3];
					P1[0] = oppositeNodes[0].X();
					P1[1] = oppositeNodes[0].Y();
					P1[2] = oppositeNodes[0].Z();

					P2[0] = oppositeNodes[1].X();
					P2[1] = oppositeNodes[1].Y();
					P2[2] = oppositeNodes[1].Z();

					PStart[0] = current_node.X();
					PStart[1] = current_node.Y();
					PStart[2] = current_node.Z();

					dirStart[0] = vec.X();
					dirStart[1] = vec.Y();
					dirStart[2] = vec.Z();
					double coeff;
					double coeffBis;
					double PEnd0;
					double PEnd1;
					double PEnd2;
					double epsilon = 0.0000000000001;

					if (intersecDirectionWithSegment(P1, P2, PStart, dirStart, coeff, coeffBis, PEnd0, PEnd1, PEnd2, epsilon)) {
						// on ajoute alors une separatrice au depart de cette singularite
						std::cout << " --> A geometrical singularity line should be created from node " << current_node.id() << std::endl;
						// TODO GESTION DES COURBES DE SINGULARITES GEOMETRIQUES
						// found = true;
						// SingularityLine* new_line = m_graph.newVolumeLine();
						///*						int sepNumber = separatrices.size();
						//						NewSep.setSepNumber(sepNumber);
						//						NewSep.SingNumber.push_back(singularities.size());
						//						NewSep.SlotSingNumber.push_back(1);
						//						NewSep.addPointInCurrentSens(SNode.getPosX(),SNode.getPosY(),SNode.getPosZ());
						//						separatrices.push_back(NewSep);
						//						std::cout<<"La separatrice "<<sepNumber<<" est cree"<<std::endl;
						//						std::cout<<"a partir du vecteur ("<<dirStart[0]<<", "<<dirStart[1]<<", "<<dirStart[2]<<")"<<std::endl;
						//						std::cout<<" et de la face "<<currentF.id()<<std::endl;
						//						*/					//mise a jour de la singularite
						// sing_point->addLine(new_line,(dirStart[0], dirStart[1], dirStart[2]);
						// SNode.setSepTriangle(currentF);

						// SNode.setSepXYZ(SNode.getPosX(), SNode.getPosY(), SNode.getPosZ());
						////		std::cout<<"NB SEP issues du noeud"<<Ntmp.id()<<" = "<<SNode.getNbSepToLaunch()<<std::endl;
						////				SNode.setSepNumber(SNode.getNbSepToLaunch(),NewSep.getSepNumber());
						////			std::cout<<"=> NB SEP issues du noeud"<<Ntmp.id()<<" = "<<SNode.getNbSepToLaunch()<<std::endl;
					}
				}
			}

			// FIN test ajout des slots pour les lancement de sep pour concavite
		}
	}
	std::cout << "NB singularities (in+geom): " << m_graph.getNbPoints() << std::endl;

	//===================================================================================================
	//		GEOMETRIC  SINGULARITY LINES
	//===================================================================================================
	// Now we have all the corner of the geom as singularity
	// We will create separatrices based on the geometric edges
	int markGeomEdgeUsedAsSing = m_mesh->newMark<Node>();
	// it_node = m_mesh->nodes_begin();
	std::vector<SingularityLine *> added_geom_lines;
	for (auto n_id : m_mesh->nodes()) {
		Node current_node = m_mesh->get<Node>(n_id);
		if ((!m_mesh->isMarked(current_node, m_markNodesOnVert)) && (m_mesh->isMarked(current_node, m_markNodesOnCurv))
		    && (!m_mesh->isMarked(current_node, markGeomEdgeUsedAsSing))) {
			// new singularity line to create here
			std::vector<int> listOfNodesInSingLeft;
			std::vector<int> listOfNodesInSingRight;
			m_mesh->mark(current_node, markGeomEdgeUsedAsSing);
			Node NodeLeft;
			Node NodeRight;
			int gotFirstOne = 0;
			for (unsigned int i = 0; i < current_node.get<Edge>().size(); i++) {
				if (m_mesh->isMarked(current_node.get<Edge>()[i], m_markEdgesOnCurv)) {
					for (unsigned int j = 0; j < 2; j++) {
						if (current_node.get<Edge>()[i].get<Node>()[j].id() != current_node.id()) {
							if (!gotFirstOne) {
								gotFirstOne = 1;
								NodeLeft = current_node.get<Edge>()[i].get<Node>()[j];
							}
							else {
								NodeRight = current_node.get<Edge>()[i].get<Node>()[j];
							}
						}
					}
				}
			}
			Node Ncurr = current_node;
			// int numberOfTheFirstSing = 0;
			// int numberOfTheSecondSing = 0;
			// Here we have the two nodes left and right that continue the line
			//		std::cout << "entree dans premiere boucle while" << std::endl;
			//------------------------------------------------------------------------------------
			// la premiere condition de la boucle while empeche la fermeture des cercles
			while ((NodeLeft.id() != NodeRight.id()) && (!m_mesh->isMarked(NodeLeft, m_markNodesOnVert))) {
				//	std::cout << "first boucle: on est dans le node " << NodeLeft.id() << std::endl;
				m_mesh->mark(NodeLeft, markGeomEdgeUsedAsSing);
				listOfNodesInSingLeft.push_back(NodeLeft.id());
				//	if (m_mesh->isMarked(NodeLeft, m_markNodesOnVert)){
				//		std::cout << "this is a corner" << std::endl;
				//	}
				/*	if (m_mesh->isMarked(NodeLeft, m_markNodesOnCurv)){
			  std::cout << "this is an edge" << std::endl;
			  }
			  std::cout << "on a un nombre de edges voisines de "
			  << NodeLeft.get<Edge>().size() << std::endl;*/
				Node NodeNext = NodeLeft;
				for (unsigned int i = 0; i < NodeLeft.get<Edge>().size(); i++) {
					if (m_mesh->isMarked(NodeLeft.get<Edge>()[i], m_markEdgesOnCurv)) {
						// std::cout << "start of edge" << std::endl;
						for (unsigned int j = 0; j < NodeLeft.get<Edge>()[i].get<Node>().size(); j++) {
							//		std::cout << "on a nodeLeft qui est " << NodeLeft.id() << " et on test le node " << NodeLeft.get<Edge>()[i].get<Node>()[j].id() <<
							//std::endl;
							if (NodeLeft.get<Edge>()[i].get<Node>()[j].id() != NodeLeft.id()) {
								if (NodeLeft.get<Edge>()[i].get<Node>()[j].id() != Ncurr.id()) {
									NodeNext = NodeLeft.get<Edge>()[i].get<Node>()[j];
								}
							}
						}
						//		std::cout << "out of edge" << std::endl;
					}
				}
				Ncurr = NodeLeft;
				NodeLeft = NodeNext;
			}
			//	std::cout << "sortie de premiere boucle while" << std::endl;
			//------------------------------------------------------------------------------------
			m_mesh->mark(NodeLeft, markGeomEdgeUsedAsSing);
			std::cout << "A" << std::endl;
			listOfNodesInSingLeft.push_back(NodeLeft.id());
			std::cout << "A" << std::endl;
			if (NodeLeft.id() != NodeRight.id()) {
				// we need to get on the right too
				Node Ncurr = current_node;
				//	std::cout << "entree dans seconde boucle while" << std::endl;
				while (!m_mesh->isMarked(NodeRight, m_markNodesOnVert)) {
					//	std::cout << "seconde boucle: on est dans le node " << NodeRight.id() << std::endl;
					m_mesh->mark(NodeRight, markGeomEdgeUsedAsSing);
					listOfNodesInSingRight.push_back(NodeRight.id());
					// if (m_mesh->isMarked(NodeRight, m_markNodesOnVert))
					//{
					//	std::cout << "this is a corner" << std::endl;
					//}
					Node NodeNext = NodeRight;
					for (unsigned int i = 0; i < NodeRight.get<Edge>().size(); i++) {
						if (m_mesh->isMarked(NodeRight.get<Edge>()[i], m_markEdgesOnCurv)) {
							for (unsigned int j = 0; j < 2; j++) {
								// std::cout << "on test le node " << NodeRight.get<Edge>()[i].get<Node>()[j].id() << std::endl;
								if (NodeRight.get<Edge>()[i].get<Node>()[j].id() != NodeRight.id()) {
									if (NodeRight.get<Edge>()[i].get<Node>()[j].id() != Ncurr.id()) {
										NodeNext = NodeRight.get<Edge>()[i].get<Node>()[j];
									}
								}
							}
						}
					}
					Ncurr = NodeRight;
					NodeRight = NodeNext;
				}
				//	std::cout << "sortie de seconde boucle while" << std::endl;
				m_mesh->mark(NodeRight, markGeomEdgeUsedAsSing);
				listOfNodesInSingRight.push_back(NodeRight.id());
			}
			else {     // cycle
				listOfNodesInSingRight.push_back(NodeLeft.id());
			}
			// Here we have the list of points to insert in the separatrix
			CurveSingularityLine *new_line = m_graph.newCurveLine();
			// we give the line id
			new_line->setNumber(m_graph.getNbLines());

			std::vector<Node> curve_nodes;
			for (unsigned int i = 0; i < listOfNodesInSingLeft.size(); i++) {
				Node NAtThisPoint = m_mesh->get<Node>(listOfNodesInSingLeft[listOfNodesInSingLeft.size() - 1 - i]);
				new_line->addDiscretizationPoint(NAtThisPoint.point());
				curve_nodes.push_back(NAtThisPoint);
			}

			new_line->addDiscretizationPoint(current_node.point());
			curve_nodes.push_back(current_node);

			for (unsigned int i = 0; i < listOfNodesInSingRight.size(); i++) {
				Node NAtThisPoint = m_mesh->get<Node>(listOfNodesInSingRight[i]);
				new_line->addDiscretizationPoint(NAtThisPoint.point());
				curve_nodes.push_back(NAtThisPoint);
			}

			// now we add the edges from the nodes in the ordered way
			std::vector<Edge> curve_edges;
			for (unsigned int i = 0; i < curve_nodes.size() - 1; i++) {
				Node current = curve_nodes[i];
				Node next = curve_nodes[i + 1];
				std::vector<Edge> current_edges = current.get<Edge>();

				bool found_edge = false;
				for (unsigned int j = 0; j < current_edges.size() && !found_edge; j++) {
					Edge ej = current_edges[j];
					if (m_mesh->isMarked(ej, m_markEdgesOnCurv)) {
						std::vector<Node> ej_nodes = ej.get<Node>();
						if (ej_nodes[0].id() == next.id() || ej_nodes[1].id() == next.id()) {
							curve_edges.push_back(ej);
							found_edge = true;
						}
					}
				}
			}
			// TODO Attention aux courbes cycliques, a priori non traitees
			new_line->setMeshEdges(curve_edges);
			added_geom_lines.push_back(new_line);
		}
	}

	// Geometrical curves made of only one mesh edges are missing. We add them now.
	// Mesh::edge_iterator itEdgeGeo = m_mesh->edges_begin();

	for (auto e_id : m_mesh->edges()) {
		Edge currentEdge = m_mesh->get<Edge>(e_id);
		// on ne traite que les aretes sur une courbe geometrique
		if (m_mesh->isMarked(currentEdge, m_markEdgesOnCurv)) {
			std::vector<Node> currentNodes = currentEdge.get<Node>();
			Node n1 = currentNodes[0];
			Node n2 = currentNodes[1];
			// on regarde si les deux noeuds correspondent ˆ des sommets geometriques
			if (m_mesh->isMarked(n1, m_markNodesOnVert)
			    && m_mesh->isMarked(n2, m_markNodesOnVert)) {     // On cree donc la separatrice reliant les singularites associŽes ˆ n1 et n2
				std::cout << "singularity line having only one mesh edge [" << n1.id() << ", " << n2.id() << "]" << std::endl;
				SingularityLine *new_line = m_graph.newCurveLine();
				new_line->setNumber(m_graph.getNbLines());

				new_line->addDiscretizationPoint(n1.point());
				new_line->addDiscretizationPoint(n2.point());
				added_geom_lines.push_back(new_line);
			}
		}
	}
	//===================================================================================================
	//	FINALLY, GEOM SINGULARITY POINTS AND GEOM SINGULARITY LINES MUST BE CONNECTED
	//===================================================================================================
	std::vector<VertexSingularityPoint *> geom_points = m_graph.getVertexPoints();
	for (unsigned int i = 0; i < geom_points.size(); i++) {
		VertexSingularityPoint *pi = geom_points[i];
		gmds::math::Point pi_point = pi->getLocation();
		for (unsigned int j = 0; j < added_geom_lines.size(); j++) {
			SingularityLine *lj = added_geom_lines[j];
			std::vector<SingularityPoint *> lj_points = lj->getEndPoints();
			if (lj_points.size() == 2) continue;     // the line is connected to its both end points

			std::vector<gmds::math::Point> &lj_discretization = lj->getDiscretizationPoints();

			gmds::math::Point lj_end_loc[2];
			gmds::math::Vector3d lj_end_dir[2];
			lj_end_loc[0] = lj_discretization[0];
			lj_end_loc[1] = lj_discretization[lj_discretization.size() - 1];
			lj_end_dir[0] = lj_discretization[1]-lj_end_loc[0];
			lj_end_dir[1] = lj_discretization[lj_discretization.size() - 2]-lj_end_loc[1];
			for (int k = 0; k < 2; k++) {
				gmds::math::Point pk = lj_end_loc[k];
				gmds::math::Vector3d v= pi_point-pk;
				if (v.norm2() < 1e-5) {
					Face f;
					// pi and lj must be connected
					// TODO CHANGE   pi->newSlot(lj_end_loc[k], lj_end_dir[k], f, lj);
					lj->addSingularityPoint(pi);
				}
			}
		}
	}
	m_mesh->unmarkAll<Node>(markGeomEdgeUsedAsSing);
	m_mesh->freeMark<Node>(markGeomEdgeUsedAsSing);
}
/*----------------------------------------------------------------------------*/
bool
SingularityGraphBuilder::isOnBoundary(const math::Point &APnt, Face &AFace)
{
	std::vector<Edge> e = AFace.get<Edge>();
	for (unsigned int i = 0; i < e.size(); i++) {
		if (isIn(APnt, e[i])) return true;
	}
	return false;
}
/*----------------------------------------------------------------------------*/
bool
SingularityGraphBuilder::isIn(const math::Point &APnt, Edge &AEdge)
{
	std::vector<Node> n = AEdge.get<Node>();
	std::cout << APnt << " in ?" << std::endl;
	std::cout << "\t " << n[0].point() << std::endl;
	std::cout << "\t " << n[1].point() << std::endl;

	math::Point p = APnt;
	math::Point p0 = n[0].point();
	math::Point p1 = n[1].point();
	// first, we look if the p is colinear with p0,p1
	math::Vector3d v01=p1-p0;
	math::Vector3d v0=p-p0;
	math::Vector3d v1=p-p1;

	math::Vector3d v01n = v01;
	v01n.normalize();
	math::Vector3d v0n = v0;
	v0n.normalize();
	std::cout << "CHECK VAL = " << fabs(v0n.dot(v01n)) << std::endl;
	if (fabs(v0n.dot(v01n) - 1) > 1e-5) {
		std::cout << "\t --> NOT ALIGNED" << std::endl;
		return false;
	}

	TCoord norme0 = v0.norm2();
	TCoord norme1 = v1.norm2();
	TCoord norm01 = v01.norm2();

	if (norme0 <= norm01 && norme1 <= norm01)
		std::cout << "\t --> YES" << std::endl;
	else
		std::cout << "\t --> NO" << std::endl;
	return (norme0 <= norm01 && norme1 <= norm01);
}
/*----------------------------------------------------------------------------*/
bool
SingularityGraphBuilder::isInside(const math::Point &APnt, const math::Vector3d &AVec, Face &AFace)
{
	std::vector<Edge> e = AFace.get<Edge>();
	std::vector<Node> nf = AFace.get<Node>();
	Edge start_edge;
	for (unsigned int i = 0; i < e.size(); i++) {
		if (isIn(APnt, e[i])) start_edge = e[i];
	}
	if (start_edge.id() == NullID) throw GMDSException("isInside(): Error, the point is not on the triangle boundary");

	// APnt is on start_edge
	// We take the point of AFace that is not incident to AEdge
	std::vector<Node> ne = start_edge.get<Node>();
	Node opp_node;
	for (unsigned int i = 0; i < ne.size(); i++) {
		if (nf[i] != ne[0] && nf[i] != ne[1]) opp_node = nf[i];
	}
	if (opp_node.id() == NullID) throw GMDSException("isInside(): Error, no opposite node");

	math::Point opp_pnt = opp_node.point();
	math::Vector3d v_test=opp_pnt-APnt;

	return (AVec.dot(v_test) >= 0.0);
}
/*----------------------------------------------------------------------------*/
bool
SingularityGraphBuilder::tryToconnectToExistingCurveSingularity(int m_markBordNode,     // boundary nodes are m_marked with m_markBordNode
                                                                int m_markBordEdge,     // boundary edges are m_marked with m_markBordEdge
                                                                SurfaceSingularityLine *ALine,
                                                                Edge &AEdge,
                                                                gmds::math::Point &APnt,
                                                                gmds::math::Vector3d &AVec)
{
	std::cout << "In tryToconnectToExistingCurveSingularity" << std::endl;
	std::vector<Node> edge_nodes = AEdge.get<Node>();
	math::Point p0 = edge_nodes[0].point();
	math::Point p1 = edge_nodes[1].point();

	std::vector<CurveSingularityPoint *> candidates = m_graph.getCurvePoints();

	std::vector<CurveSingularityPoint *> matching;
	for (unsigned int i = 0; i < candidates.size(); i++) {
		CurveSingularityPoint *current_pnt = candidates[i];
		if (APnt.distance(current_pnt->getLocation())) matching.push_back(current_pnt);
	}

	// We haven't found any candidate
	if (matching.empty()) return false;

	if (matching.size() > 1) throw GMDSException("Unexpected configuration in SingularityGraphBuilder::tryToconnectToExistingCurveSingularity");

	std::vector<Edge> adj_edges;

	// We have one candidate
	CurveSingularityPoint *sing_point = matching[0];

	// we get the line
	// TODO Check this addition without any vector
	sing_point->addLine(ALine);
	ALine->addSingularityPoint(sing_point);
	ALine->addDiscretizationPoint(sing_point->getLocation());
	return true;
}
/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::computeSurfaceSingularityLine(SingularityPoint *AFromPoint, SingularityPoint::Slot *AFromSlot)
{
	std::cout << "==== Computation of a surface singularity line ====" << std::endl;
	math::Point startPnt = AFromSlot->location;         // starting point
	math::Vector3d startDir = AFromSlot->direction;     // starting direction
	Face startFace;                                     // TODO CHANGE = AFromSlot->starting_triangle;			//starting triangle

	if (startFace.id() != NullID)
		std::cout << "Starting face = " << startFace.id() << std::endl;
	else
		std::cout << "Starting face = NULL" << std::endl;
	//========================================================================
	// Initialization of the line we create
	//========================================================================
	SurfaceSingularityLine *surf_line = m_graph.newSurfaceLine();
	int sepNumberTmp = m_graph.getNbLines();
	surf_line->setNumber(sepNumberTmp);

	Face current_face = startFace;
	math::Point current_pnt = startPnt;
	math::Vector3d current_vec = startDir;

	std::cout << "START PNT: " << startPnt << std::endl;
	std::cout << "START DIR: " << startDir << std::endl;
	math::Point startDirPnt(startPnt.X() + startDir.X(), startPnt.Y() + startDir.Y(), startPnt.Z() + startDir.Z());

	std::cout << "START DIR PNT: " << startDirPnt << std::endl;
	math::Vector3d prev_vec;
	Face prev_face;
	SingularityPoint *from_sing_pnt = AFromSlot->from_point;
	//========================================================================
	// connect line to singularity point
	//========================================================================
	surf_line->addSingularityPoint(from_sing_pnt);
	surf_line->addDiscretizationPoint(from_sing_pnt->getLocation());
	//========================================================================
	// connect point (i.e. its slot) to the new line
	AFromSlot->line = surf_line;
	AFromSlot->isLaunched = true;

	if (from_sing_pnt->getGeomType() == SingularityPoint::SURFACE) {
		SurfaceSingularityPoint *from_surf_sing_pnt = dynamic_cast<SurfaceSingularityPoint *>(from_sing_pnt);
		prev_face = from_surf_sing_pnt->getMeshFace();
		std::cout << "we have a surface sing point" << std::endl;
	}

	// first we add the singularity point of the slot ,and we define prev_vec
	if (prev_face.id() == NullID) {
		std::cout << "we have a prev face" << std::endl;
		prev_vec = current_pnt-AFromSlot->from_point->getLocation();
		std::cout << " PREV VEC: " << prev_vec << std::endl;
		math::Point prevPnt(startPnt.X() + prev_vec.X(), startPnt.Y() + prev_vec.Y(), startPnt.Z() + prev_vec.Z());

		std::cout << "PREV PNT: " << prevPnt << std::endl;
	}
	std::vector<Node> prev_face_nodes, cur_face_nodes;
	prev_face_nodes = prev_face.get<Node>();
	cur_face_nodes = current_face.get<Node>();
	std::cout << "====== CURRENT FACE:";
	for (unsigned in = 0; in < cur_face_nodes.size(); in++) {
		std::cout << " " << cur_face_nodes[in].id();
	}
	std::cout << std::endl;
	std::cout << "====== PREV   FACE:";
	for (unsigned in = 0; in < prev_face_nodes.size(); in++) {
		std::cout << " " << prev_face_nodes[in].id();
	}
	std::cout << std::endl;

	bool find_end = false;
	do {

		// We arrive in triangle current_face from fromPnt following vector fromVec
		// fromPnt must belong to the boundary of current_face and
		// fromDir must be oriented inside current_face
		// if (!isOnBoundary(current_pnt, current_face))
		//	throw GMDSException("ERROR: fromPnt is not on the boundary of current_face");
		// if (!isInside(current_pnt, current_vec, current_face))
		//	throw GMDSException("ERROR: (fromPnt,fromVec) is not inside current_face");

		//==============================================================
		// CASE 1: DOES THE FACE BE CLASSIFIED ONTO ANOTHER GEOM SURF?
		//==============================================================
		bool reach_geom_boundary = false;
		if (prev_face.id() != NullID) {
			std::cout << "\n --> we check if we are on a geometric curve" << std::endl;
			// At the beginning, from_face can be null, meaning we start
			// from a geometrical singularity point

			// We check the classification
			// if different, we stop and build a geom_sing_point
			if ((*m_var_color)[prev_face.id()] != (*m_var_color)[current_face.id()]) {
				std::cout << "\t YES WE ARE" << std::endl;

				reach_geom_boundary = true;
				// we get the edge common to from_face and current_face
				Edge common_edge;
				bool found_common = false;
				std::vector<Edge> prev_edges = prev_face.get<Edge>();
				std::vector<Edge> current_edges = current_face.get<Edge>();

				for (unsigned int from_i = 0; from_i < prev_edges.size() && !found_common; from_i++) {
					for (unsigned int cur_i = 0; cur_i < current_edges.size() && !found_common; cur_i++) {
						if (current_edges[cur_i] == prev_edges[from_i]) {
							found_common = true;
							common_edge = current_edges[cur_i];
						}
					}
				}
				std::cout << "common edge = " << common_edge.get<Node>()[0].id() << ", " << common_edge.get<Node>()[1].id() << std::endl;

				if (!found_common) throw GMDSException("ERROR, no common edge between two adjacent faces");

				// we compute slots info
				// we create a geom sing.
				// such a singularity has 4 connected singularity lines
				// - 2 for the geometric curves that is split in two parts
				// - 1 per adjacent surface
				// But the two first do not open a slot (they are fixed)

				// math::Vector3d opp_prev_vec(-prev_vec.X(), -prev_vec.Y(), -prev_vec.Z());

				// if (!tryToconnectToExistingCurveSingularity(m_markBordNode, m_markBordEdge, surf_line,
				//	common_edge, current_pnt, opp_prev_vec))
				//{
				//	m_graph.splitCurveLine(surf_line, current_pnt, common_edge,
				//		prev_face, current_face,
				//		opp_prev_vec, current_vec);

				//}
				// In both cases, we reach the end of the line
				find_end = true;

			}     // if ((*m_var_color)[prev_face.id()] != (*m_var_color)[current_face.id()])
		}        // if (prev_face.id() != NullID)

		//==============================================================
		// CASE 2: DO WE ENTER IN A FACE CONTAINING ANOTHER SING. LINE?
		//==============================================================
		bool intersect_line = false;
		if (!reach_geom_boundary) {
			// TODO Le cas ou on intersecte une autre sing. line n'est pas
			// encore gere
			////the current_face can contain another singularity line
			// std::vector<SurfaceSingularityLine*> intersected_lines =
			//	getSurfaceSingularityLineIn(current_face);
			// if (!intersected_lines.empty())
			//{
			//	//We compute if we intersect a line of intersected_lines or not
			//	bool do_intersect = false;
			//	//YES
			//	if (do_intersect)
			//	{
			//		intersect_line = true;
			//		//computation of the singularity point and split of the
			//		// intersected sing. line
			//	}
			//}
		}

		//==============================================================
		// CASE 3: DO WE ENTER IN A FACE CONTAINING A SING. POINT?
		//==============================================================
		bool intersect_sing_point = false;
		if (!reach_geom_boundary && !intersect_line) {
			SurfaceSingularityPoint *next_sing_point = m_face_to_singularity_on_surf[current_face.id()];
			bool must_connect = false;
			if (next_sing_point != NULL &&         // face in a confusing ball ...
			    next_sing_point != AFromPoint)     // ... of another singularity point
			{
				must_connect = true;
			}
			else if (next_sing_point != NULL &&         // face in the confusing ball ...
			         next_sing_point == AFromPoint)     //... of the incoming singularity point
			{
				// Warning: completly empiric
				if (surf_line->getDiscretizationPoints().size() >= 10) must_connect = true;
			}

			if (must_connect) {
				std::cout << "\n --> triangle " << current_face.id() << " with sing. point" << std::endl;
				intersect_sing_point = true;
				// We get the surface singularity located in current_face
				std::vector<SurfaceSingularityPoint *> surf_points = m_graph.getSurfacePoints();
				SurfaceSingularityPoint *sing = 0;
				std::cout << "nb. surf. points = " << surf_points.size() << std::endl;
				for (unsigned int i_pnt = 0; sing == 0 && i_pnt < surf_points.size(); i_pnt++) {
					SurfaceSingularityPoint *cur_surf_pnt = surf_points[i_pnt];
					Face cur_sing_face = cur_surf_pnt->getMeshFace();
					std::cout << i_pnt << " sing. point -> " << cur_sing_face.id() << std::endl;
					if (cur_sing_face == current_face) sing = cur_surf_pnt;
				}

				if (sing == 0) {
					std::cout << "ERROR, a bordSing-m_marked face doesn't contain a sing point" << std::endl;

					throw GMDSException("ERROR, a bordSing-m_marked face doesn't contain a sing point");
				}

				// Now, we look for an available slot
				std::vector<SingularityPoint::Slot *> &cur_slots = sing->getSlots();

				bool found_free_slot = false;

				double slot_connection_tol = -0.9;
				while (!found_free_slot) {
					for (unsigned int i_slot = 0; !found_free_slot && i_slot < cur_slots.size(); i_slot++) {
						SingularityPoint::Slot *current_slot = cur_slots[i_slot];
						if (current_slot->isLaunched == false) {
							// It's a free slot, we test it
							math::Vector3d slot_dir = current_slot->direction;
							slot_dir.normalize();
							current_vec.normalize();
							std::cout << "SLOT DIR: " << slot_dir << std::endl;
							std::cout << "CURR VEC: " << current_vec << std::endl;
							if (slot_dir.dot(current_vec) < slot_connection_tol) {
								found_free_slot = true;
								current_slot->isLaunched = true;
								current_slot->line = surf_line;
								surf_line->addSingularityPoint(sing);
								surf_line->addDiscretizationPoint(sing->getLocation());
							}
						}
					}
					slot_connection_tol = slot_connection_tol + 0.1;
				}
				find_end = true;
				// SECOND
				////We compute if we intersect a line of intersected_lines or not
				// bool do_intersect = false;
				////YES NOT YET HANDLED
				// if (do_intersect)
				//{
				//	intersect_sing_point = true;
				//	//computation of the singularity point and split of the
				//	// intersected sing. line
				//}
			}
		}     // if (!reach_geom_boundary && !intersect_line)
		//==============================================================
		// CASE 4: SO, WE JUST HAVE TO CROSS THE TRIANGLE (GENERAL CONF)
		//==============================================================
		// Does the current triangle has the same classif
		if (!reach_geom_boundary && !intersect_line && !intersect_sing_point) {
			std::cout << "\n --> we are in a smooth field triangle" << std::endl;

			math::Point toPnt;
			math::Vector3d toVec;
			Face toFace;

			computeOutTriangleOnTheSurface(current_face, current_pnt, current_vec, toFace, toPnt, toVec);
			// we keep the point toPnt to define the line

			surf_line->addDiscretizationPoint(toPnt);
			// we progress to the next point, next vector and so next face too
			prev_face = current_face;
			prev_vec = current_vec;
			current_pnt = toPnt;
			current_vec = toVec;
			current_face = toFace;
		}

		if (reach_geom_boundary || intersect_line || intersect_sing_point) find_end = true;

	} while (!find_end);

	std::cout << "Getting out from ComputeSepLine" << std::endl;
}
/*----------------------------------------------------------------------------*/
bool
SingularityGraphBuilder::computeIntersection(
   gmds::math::Point &APnt, gmds::math::Vector3d &AVec, gmds::Edge AEdge, gmds::math::Point &AOUTPnt, double &AOUTParam, double &AOUTParamRay)
{
	std::cout << "We compute an intersection between an edge and a ray" << std::endl;

	double epsilon = 1e-4;
	std::vector<Node> edge_nodes = AEdge.get<Node>();

	math::Point q0 = APnt;
	math::Point p0 = edge_nodes[0].point();
	math::Point p1 = edge_nodes[1].point();

	math::Vector3d q0q1(AVec);
	math::Vector3d q0p0=p0-q0;
	math::Vector3d q0p1=p1-q0;
	math::Vector3d p0p1=p1-p0;

	// by default (APnt,AVec) and AEdge should not be colinear,
	// and APnt, p1 and p2 are assumed as being coplanar
	math::Vector3d col1(q0q1);
	math::Vector3d col2(p0p1);
	col1.normalize();
	col2.normalize();
	double colinear = fabs(col1.dot(col2));
	if (fabs(1 - colinear) < epsilon) return false;


	// We have an intersection point to compute
	// coeff of the line L0(Q0, V0) equation sytem
	double vx = q0q1.X();
	double vy = q0q1.Y();
	double vz = q0q1.Z();

	// coeff of the line L1(P0, U0) equation sytem
	double p0x = p0.X();
	double p0y = p0.Y();
	double p0z = p0.Z();
	double ux = p0p1.X();
	double uy = p0p1.Y();
	double uz = p0p1.Z();

	double wx = q0p0.X();
	double wy = q0p0.Y();
	double wz = q0p0.Z();

	double u1 = 0;
	double u2 = 0;
	double v1 = 0;
	double v2 = 0;
	double w1 = 0;
	double w2 = 0;
	if (fabs(ux * vy - uy * vx) > epsilon)     // Work in XY
	{
		u1 = ux;
		u2 = uy;
		v1 = vx;
		v2 = vy;
		w1 = wx;
		w2 = wy;
	}
	else if (fabs(ux * vz - uz * vx) > epsilon)     // work in XZ
	{
		u1 = ux;
		u2 = uz;
		v1 = vx;
		v2 = vz;
		w1 = wx;
		w2 = wz;
	}
	else if (fabs(uz * vy - uy * vz) > epsilon)     // work in YZ
	{
		u1 = uy;
		u2 = uz;
		v1 = vy;
		v2 = vz;
		w1 = wy;
		w2 = wz;
	}
	else {
		std::cout << "No projection plane found to compute a (ray,segment) intersection" << std::endl;
		throw GMDSException("ERROR: No projection plane found to compute a(ray, segment) intersection");
	}

	// parameter of the intersection pnt for L0
	double t0 = (u1 * w2 - u2 * w1) / (u1 * v2 - u2 * v1);
	// parameter of the intersection pnt for L1
	double t1 = (v2 * w1 - v1 * w2) / (v1 * u2 - v2 * u1);
	AOUTParam = t1;
	AOUTParamRay = t0;
	if (t0 < 0) {
		std::cout << "OUT NO SOLUTION due to backward option" << std::endl;
		return false;
	}
	if (t1 < 0 || t1 > 1) {
		std::cout << "OUT NO SOLUTION wit paramater" << t1 << std::endl;
		return false;
	}

	// POINT
	AOUTPnt.setXYZ(p0x + t1 * ux, p0y + t1 * uy, p0z + t1 * uz);
	return true;
}
/*----------------------------------------------------------------------------*/
int
SingularityGraphBuilder::computeOutDataOnTheSurface(const gmds::math::Vector3d &ANormal,
                                                    gmds::math::Point &AINPnt,
                                                    gmds::math::Vector3d &AINVec,
                                                    gmds::Node &AINNode1,     // val=1
                                                    gmds::Node &AINNode2,     // val=2
                                                    gmds::Node &AINNode3,     // val=3
                                                    gmds::Edge &AINEdge1,     // val=4
                                                    gmds::Edge &AINEdge2,     // val=5
                                                    gmds::math::Point &AOUTPnt,
                                                    gmds::math::Vector3d &AOUTVec)
{
	int val = -1;
	std::cout << "===================================" << std::endl;
	std::cout << "compute OUT QUAT from" << std::endl;
	std::cout << "\t AIN Pnt: " << AINPnt << std::endl;
	std::cout << "\t AIN Vec: " << AINVec << std::endl;
	math::Point pout = AINPnt + AINVec;
	std::cout << "\t \t out pnt: " << pout << std::endl;
	std::cout << "\t Node candidate : " << AINNode1.id() << std::endl;
	std::cout << "\t Edge candidate1: " << AINEdge1.get<Node>()[0].id() << ",  " << AINEdge1.get<Node>()[1].id() << std::endl;
	std::cout << "\t Edge candidate2: " << AINEdge2.get<Node>()[0].id() << ",  " << AINEdge2.get<Node>()[1].id() << std::endl;

	std::cout << "Node " << AINNode1.id() << " - " << AINNode1.point() << std::endl;
	std::cout << "Node " << AINNode2.id() << " - " << AINNode2.point() << std::endl;
	std::cout << "Node " << AINNode3.id() << " - " << AINNode3.point() << std::endl;
	math::Vector3d v_in = AINVec;
	v_in.normalize();
	math::Point node_loc = AINNode1.point();
	math::Vector3d v1=node_loc-AINPnt;
	v1.normalize();
	math::Quaternion q;
	if (fabs(v1.dot(v_in) - 1) < 1e-5) {
		std::cout << "INTERSECT NODE --> out in a node" << std::endl;
		// throw GMDSException("ERROR: OUT INTO A NODE");
		// intersected point = node
		AOUTPnt = node_loc;
		val = 1;
		computeVecInSurfacePoint(AINNode1, ANormal, AINVec, AOUTVec);
	}
	else {
		// We have to intersect (AINPnt,AINVec) with each edge

		gmds::math::Point pnt;
		Edge out_edge;
		double param = 0.0;
		double param1 = 0.0;
		double param2 = 0.0;
		double paramFrom1 = 0.0;
		double paramFrom2 = 0.0;

		bool intersectEdge1 = false;
		intersectEdge1 = computeIntersection(AINPnt, AINVec, AINEdge1, pnt, param1, paramFrom1);

		bool intersectEdge2 = false;

		if (!intersectEdge1) intersectEdge2 = computeIntersection(AINPnt, AINVec, AINEdge2, pnt, param1, paramFrom2);

		if (intersectEdge1) {
			std::cout << "intersect the first edge" << std::endl;
			val = 4;
			out_edge = AINEdge1;
			param = param1;
		}
		else if (intersectEdge2) {
			std::cout << "intersect the seconde edge" << std::endl;
			val = 5;
			out_edge = AINEdge2;
			param = param2;
		}
		else     // We do not intersect an edge but a node
		{
			gmds::TCellID nodeID = NullID;
			// no direct intersection
			if (paramFrom1 < 0 && paramFrom2 > 0) {
				std::vector<Node> e_nodes = AINEdge2.get<Node>();
				// We use nodes of curve 2
				if (fabs(param2) < 0.2) {
					// we are close to the first edge point
					math::Point node_loc = e_nodes[0].point();
					AOUTPnt = node_loc;
					nodeID = e_nodes[0].id();
				}
				else {
					// we are close to the first edge point
					math::Point node_loc = e_nodes[1].point();
					AOUTPnt = node_loc;
					nodeID = e_nodes[1].id();
				}
			}
			else if (paramFrom2 < 0 && paramFrom1 > 0) {
				std::vector<Node> e_nodes = AINEdge1.get<Node>();
				// We use nodes of curve 1
				if (fabs(param1) < 0.2) {
					// we are close to the first edge point
					math::Point node_loc = e_nodes[0].point();
					AOUTPnt = node_loc;
					nodeID = e_nodes[0].id();
				}
				else {
					// we are close to the first edge point
					math::Point node_loc = e_nodes[0].point();
					AOUTPnt = node_loc;
					nodeID = e_nodes[1].id();
				}
			}
			else {
				std::cout << "ERROR: No out point in SingularityGraphBuilder::computeOutDirection" << std::endl;
				throw GMDSException("ERROR: No out point in SingularityGraphBuilder::computeOutDirection");
			}

			Node outNode;
			if (AINNode1.id() == nodeID) {
				outNode = AINNode1;
				val = 1;
			}
			else if (AINNode2.id() == nodeID) {
				outNode = AINNode2;
				val = 2;
			}
			else if (AINNode3.id() == nodeID) {
				outNode = AINNode3;
				val = 3;
			}
			else {
				std::cout << "ERROR: a bad id is used for a node" << std::endl;
				throw GMDSException("ERROR: No out point in SingularityGraphBuilder::computeOutDirection");
			}
			std::cout << "Param (" << paramFrom1 << ", " << paramFrom2 << ") -> " << val << std::endl;
			std::cout << "INTERSECT NODE ALONG EDGE --> out in a node" << std::endl;
			computeVecInSurfacePoint(outNode, ANormal, AINVec, AOUTVec);
		}     // else //We do not intersect an edge but a node

		if (val == 4 || val == 5) {

			AOUTPnt = pnt;
			Node n0 = out_edge.get<Node>()[0];
			Node n1 = out_edge.get<Node>()[1];
			std::cout << "face normal: " << ANormal << std::endl;
			// We get the cross in each extrem point
			math::Cross c0((*m_var_quatern)[n0.id()], ANormal);
			math::Cross c1((*m_var_quatern)[n1.id()], ANormal);
			std::cout << "cross in " << n0.id() << ": " << c0 << std::endl;
			std::cout << "cross in " << n1.id() << ": " << c1 << std::endl;

			math::Vector3d c0_closestC1X = c0.closestVector(c1.X());

			math::Vector3d v_out1 = (1 - param) * c0_closestC1X + param * c1.X();

			math::Vector3d v_normal = ANormal;
			v_normal.normalize();
			v_out1.normalize();

			math::Vector3d v_out2 = v_out1.cross(v_normal);
			math::Vector3d v_out3 = v_out2.cross(v_normal);
			math::Vector3d v_out4 = v_out3.cross(v_normal);

			double val1 = AINVec.dot(v_out1);
			double val2 = AINVec.dot(v_out2);
			double val3 = AINVec.dot(v_out3);
			double val4 = AINVec.dot(v_out4);
			if (val1 >= val2 && val1 >= val3 && val1 >= val4) {
				AOUTVec = v_out1;
			}
			else if (val2 >= val1 && val2 >= val3 && val2 >= val4) {
				AOUTVec = v_out2;
			}
			else if (val3 >= val1 && val3 >= val1 && val3 >= val4) {
				AOUTVec = v_out3;
			}
			else {
				AOUTVec = v_out4;
			}
		}     // if (val == 4 || val == 5)
	}

	return val;
}
/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::computeVecInSurfacePoint(const Node &ANode, const math::Vector3d &ANormal, math::Vector3d &AFromVec, math::Vector3d &AToVec)
{

	math::Vector3d v_normal = ANormal;
	v_normal.normalize();
	math::Cross cross_in_node((*m_var_quatern)[ANode.id()], ANormal);

	math::Vector3d v1 = cross_in_node.X();
	v1.normalize();
	math::Vector3d v2 = v1.cross(v_normal);
	math::Vector3d v3 = v2.cross(v_normal);
	math::Vector3d v4 = v3.cross(v_normal);

	double val1 = AFromVec.dot(v1);
	double val2 = AFromVec.dot(v2);
	double val3 = AFromVec.dot(v3);
	double val4 = AFromVec.dot(v4);
	if (val1 >= val2 && val1 >= val3 && val1 >= val4) {
		AToVec = v1;
	}
	else if (val2 >= val1 && val2 >= val3 && val2 >= val4) {
		AToVec = v2;
	}
	else if (val3 >= val1 && val3 >= val1 && val3 >= val4) {
		AToVec = v3;
	}
	else {
		AToVec = v4;
	}
}
/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::computeOutTriangleOnTheSurface(gmds::Face &ACurrentFace,
                                                        gmds::math::Point &AINPnt,
                                                        gmds::math::Vector3d &AINVec,
                                                        gmds::Face &AOUTFace,
                                                        gmds::math::Point &AOUTPnt,
                                                        gmds::math::Vector3d &AOUTVec)
{
	//=====================================================================
	// STEP 1: We find the edge of ACurrentFace containing AINPnt
	//=====================================================================
	std::vector<Edge> current_edges = ACurrentFace.get<Edge>();
	std::vector<Node> current_nodes = ACurrentFace.get<Node>();

	std::cout << "CURRENT FACE " << ACurrentFace.id() << ": " << std::endl;
	for (unsigned int i = 0; i < current_nodes.size(); i++)
		std::cout << " " << current_nodes[i].id();
	std::cout << std::endl;

	Edge edge_in;
	for (unsigned int i = 0; /*edge_in.id() == NullID &&*/ i < current_edges.size(); i++) {
		if (isIn(AINPnt, current_edges[i])) {
			std::cout << "in edge " << i << std::endl;
			edge_in = current_edges[i];
		}
	}

	if (edge_in.id() == NullID) {
		std::cout << "ERROR in SingularityGraphBuilder::computeOutTriangle - No edge contains AINPnt" << std::endl;
		std::cout << "pnt: " << AINPnt << std::endl;
		for (unsigned int i = 0; i < current_edges.size(); i++) {
			std::cout << "Edge " << i << std::endl;
			Edge ei = current_edges[i];
			std::vector<Node> nei = ei.get<Node>();
			for (unsigned int j = 0; j < nei.size(); j++) {
				Node nj = nei[j];
				std::cout << "    -> " << nj.point() << std::endl;
			}
		}
		throw GMDSException("ERROR in SingularityGraphBuilder::computeOutTriangle - No edge contains AINPnt");
	}

	Edge other_edges[2];
	int nb_edges = 0;
	for (unsigned int i = 0; i < current_edges.size(); i++) {
		if (current_edges[i] != edge_in) other_edges[nb_edges++] = current_edges[i];
	}

	Node other_node;
	std::vector<Node> nodes_in = edge_in.get<Node>();
	for (unsigned int i = 0; i < current_nodes.size(); i++) {
		Node current_node = current_nodes[i];
		if (current_node != nodes_in[0] && current_node != nodes_in[1]) other_node = current_node;
	}

	std::cout << "EDGE IN  : " << nodes_in[0].id() << ", " << nodes_in[1].id() << std::endl;
	std::vector<Node> nodes_out0 = other_edges[0].get<Node>();
	std::vector<Node> nodes_out1 = other_edges[1].get<Node>();
	// std::cout << "EDGE OUT0: " << nodes_out0[0].id() << ", " << nodes_out0[1].id() << std::endl;
	// std::cout << "EDGE OUT1: " << nodes_out1[0].id() << ", " << nodes_out1[1].id() << std::endl;
	std::cout << "OTHER TRIANGLE NODE: " << other_node.id() << std::endl;
	//=====================================================================
	// STEP 2: We compute a first out pnt and vector
	//=====================================================================
	std::cout << "========================= Heun's 1 ====================" << std::endl;
	math::Point out_pnt;
	math::Vector3d out_vec;
	AINVec.normalize();
	computeOutDataOnTheSurface(ACurrentFace.normal(), AINPnt, AINVec, other_node, edge_in.get<Node>()[0], edge_in.get<Node>()[1], other_edges[0], other_edges[1],
	                           out_pnt, out_vec);

	std::cout << "First computed out point: " << out_pnt << std::endl;
	std::cout << "First computed out vect: " << out_vec << std::endl;
	//=====================================================================
	// STEP 3: We recompute out pnt and vector using Heun's method
	//=====================================================================
	std::cout << "========================= Heun's 2 ====================" << std::endl;
	std::cout << "Recomputation via Heun's method" << std::endl;
	// math::Vector3d out_vec = cross(out_q, ACurrentFace.normal()).closestVector(AINVec);

	gmds::math::Vector3d v_Heun = 0.5 * AINVec + 0.5 * out_vec;
	v_Heun.normalize();

	int out_id = computeOutDataOnTheSurface(ACurrentFace.normal(), AINPnt, v_Heun, other_node, edge_in.get<Node>()[0], edge_in.get<Node>()[1], other_edges[0],
	                                        other_edges[1], out_pnt, out_vec);
	// We recompute with 1/2 vi + 1/2 out_vec
	// out_q = computeOutmath::quaternion(ACurrentFace.normal(),
	//	AINPnt, v_Heun,
	//	other_node, other_edges[0], other_edges[1],
	//	out_pnt, out_element);

	std::cout << "Second computed out point: " << out_pnt << std::endl;

	// we leave the current face along the first or the second edge
	AOUTPnt = out_pnt;

	if (out_id == 4 || out_id == 5) {

		Edge out_edge;
		if (out_id == 4)
			out_edge = other_edges[0];
		else if (out_id == 5)
			out_edge = other_edges[1];

		std::cout << "OUT EDGE = " << out_edge.get<Node>()[0].id() << ", " << out_edge.get<Node>()[1].id() << std::endl;
		std::vector<Face> candidate_out_faces = out_edge.get<Face>();
		std::wcout << "Nb adj faces = " << candidate_out_faces.size() << std::endl;
		for (unsigned int i = 0; i < candidate_out_faces.size(); i++) {
			Face current_face = candidate_out_faces[i];
			if (current_face != ACurrentFace && m_mesh->isMarked(current_face, m_markFacesOnSurf)) AOUTFace = current_face;
		}
		if (AOUTFace.id() == NullID) std::cout << "No out face" << std::endl;

		math::Vector3d normal = AOUTFace.normal();
		// cross c(out_q, AOUTFace.normal());
		// std::cout << "CROSS INTERPOLATED: " << c << std::endl;
	}
	else {
		Node out_node;

		if (out_id == 1) {
			out_node = other_node;
		}
		else if (out_id == 2) {
			out_node = edge_in.get<Node>()[0];
		}
		else if (out_id == 3) {
			out_node = edge_in.get<Node>()[1];
		}
		TCellID out_node_id = out_node.id();

		std::vector<Face> candidate_out_faces = out_node.get<Face>();
		std::vector<Face> candidate_out_bnd_faces;
		std::wcout << "Nb adj faces = " << candidate_out_faces.size() << std::endl;
		for (unsigned int i = 0; i < candidate_out_faces.size(); i++) {
			Face current_face = candidate_out_faces[i];
			if (current_face != ACurrentFace && m_mesh->isMarked(current_face, m_markFacesOnSurf)) candidate_out_bnd_faces.push_back(current_face);
		}

		Node other_nodes[2];
		for (unsigned int i = 0; i < candidate_out_bnd_faces.size(); i++) {
			Face current_face = candidate_out_bnd_faces[i];
			std::vector<Node> current_nodes = current_face.get<Node>();

			int other_id = 0;
			for (unsigned int j = 0; j < current_nodes.size(); j++) {
				Node nj = current_nodes[j];
				if (nj.id() != out_node_id) {
					other_nodes[other_id++] = nj;
				}
			}
		}
		math::Vector3d v0=other_nodes[0].point()-out_node.point();
		math::Vector3d v1=other_nodes[1].point()-out_node.point();

		math::Vector3d normal = AOUTFace.normal();
		// cross c(out_q, AOUTFace.normal());
		// std::cout << "CROSS INTERPOLATED: " << c << std::endl;
	}

	// AOUTPnt and AOUTFace are computed, we have now to project out_q onto AOUTFace
	AOUTVec = out_vec;     // c.closestVector(AINVec);//v_Heun); //TODO VERIFY
}
/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::computeShortSepLine(const math::Point &startPnt,
                                             const math::Vector3d &startDir,
                                             Face &startFace,     // next face
                                             const math::Point &singPnt,
                                             Face &singFace,
                                             std::vector<math::Point> &points,
                                             int sepNumber)
{
	std::cout << "entree dans ComputeShortSepLine" << std::endl;

	std::cout << "On va de sing  a start, raf" << std::endl;
	points.push_back(singPnt);
	points.push_back(startPnt);

	//	separatrices[sepNumber].isToBeAssociatedWithAnotherSep = 1;
	std::vector<TCellID> common = m_mesh->getCommonNodes(startFace, singFace);
	if (common.size() != 2) {
		std::cout << "ERREUR -> " << common.size() << std::endl;
		exit(1);
	}
	//	separatrices[sepNumber].NodeEndID[0] = common[0].id();
	//	separatrices[sepNumber].NodeEndID[1] = common[1].id();
	std::cout << "sortie de ComputeShortSepLine" << std::endl;
}
/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::createBoundarySingularityGraphFromField()
{
	// Now we need to create separatrices on the faces
	// we will start by launching them from the slot of the singularities
	// for each one we have the pos and dir to launch the sep into, as well as the first triangle
	std::cout << "== Boundary singularity graph extraction from frame field ==" << std::endl;

	std::vector<SingularityPoint *> singularity_points = m_graph.getPoints();
	std::cout << "Nb sing points = " << singularity_points.size() << std::endl;
	int nb_done = 0;
	for (unsigned int i = 0; i < singularity_points.size(); i++) {

		SingularityPoint *pi = singularity_points[i];

		if (pi->getGeomType() == SingularityPoint::VERTEX)
			std::cout << i << " -> VERTEX" << std::endl;
		else if (pi->getGeomType() == SingularityPoint::CURVE)
			std::cout << i << " -> CURVE" << std::endl;
		else if (pi->getGeomType() == SingularityPoint::SURFACE)
			std::cout << i << " -> SURFACE" << std::endl;
		else if (pi->getGeomType() == SingularityPoint::VOLUME)
			std::cout << i << " -> VOLUME" << std::endl;
		else
			std::cout << i << " -> UNDEF UNDEF UNDEF UNDEF" << std::endl;

		// Volume singularities have not to be handled
		if (pi->getGeomType() == SingularityPoint::VOLUME) continue;

		// Curve and vertex singularities will be handled in a second step
		if (pi->getGeomType() == SingularityPoint::VERTEX || pi->getGeomType() == SingularityPoint::CURVE) continue;

		if (pi->getNbMeshCells() != 1) continue;
		std::vector<SingularityPoint::Slot *> pi_slots = pi->getSlots();

		for (unsigned int j = 0; j < pi_slots.size(); j++) {
			SingularityPoint::Slot *slot_j = pi_slots[j];
			if (!slot_j->isLaunched && slot_j->isOnSurface) {
				slot_j->isLaunched = true;

				computeSurfaceSingularityLine(pi, slot_j);
			}
			writeOutput("boundary_line");
			nb_done++;
			// if (nb_done == 2)
			//	return;
		}
	}
	std::cout << "== DONE ==" << std::endl;
}

/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::extractQuaternions()
{

	std::cout << "Start math::quaternion definition" << std::endl;

	Variable<math::Vector3d> *varX = m_mesh->getVariable<math::Vector3d, GMDS_NODE>("quaternion_X");
	Variable<math::Vector3d> *varY = m_mesh->getVariable<math::Vector3d, GMDS_NODE>("quaternion_Y");
	Variable<math::Vector3d> *varZ = m_mesh->getVariable<math::Vector3d, GMDS_NODE>("quaternion_Z");

	m_var_quatern = m_mesh->newVariable<math::Quaternion, GMDS_NODE>("quaternion_field");

	for (auto n_id : m_mesh->nodes()) {
		Node n = m_mesh->get<Node>(n_id);
		math::Vector3d vx = (*varX)[n.id()];
		math::Vector3d vy = (*varY)[n.id()];
		math::Vector3d vz = (*varZ)[n.id()];
		math::Quaternion q(math::Chart(vx, vy, vz));
		(*m_var_quatern)[n.id()] = q;
	}
	std::cout << "\t -> quaternions created!" << std::endl;
}
/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::colorFaces(const int m_markF, const int m_markE)
{

	try {
		m_var_color = m_mesh->getVariable<int, GMDS_FACE>("color");
	}
	catch (GMDSException &e) {
		std::cout << e.what() << std::endl;
		m_var_color = m_mesh->newVariable<int, GMDS_FACE>("color");
	}
	int color = 1;     // par defaut 0 pour tout le monde je pense
	int m_markDone = m_mesh->newMark<Face>();

	for (auto f_id : m_mesh->faces()) {
		Face f = m_mesh->get<Face>(f_id);
		// on ne considere que les faces au bord
		if (!m_mesh->isMarked(f, m_markF)) continue;

		// qui n'ont pas encore ete traitees
		if (m_mesh->isMarked(f, m_markDone)) continue;

		// nouvelle surface, yep
		color++;     // donc nouvelle couleur
		m_mesh->mark(f, m_markDone);
		(*m_var_color)[f.id()] = color;

		// on se propage et on marque
		std::vector<Face> next;
		next.push_back(f);

		while (!next.empty()) {
			Face current = next[next.size() - 1];
			next.pop_back();
			// recuperation des faces voisines non traitees et appartenant a la surface
			std::vector<Edge> current_edges = current.get<Edge>();
			for (unsigned int ie = 0; ie < current_edges.size(); ie++) {
				Edge ei = current_edges[ie];
				if (!m_mesh->isMarked(ei, m_markE))     // si ce n'est pas une arete au bord
				{
					std::vector<Face> f_edges = ei.get<Face>();
					for (unsigned int ifa = 0; ifa < f_edges.size(); ifa++) {
						Face fi = f_edges[ifa];
						if (m_mesh->isMarked(fi, m_markF) && !m_mesh->isMarked(fi, m_markDone)) {
							m_mesh->mark(fi, m_markDone);
							(*m_var_color)[fi.id()] = color;
							next.push_back(fi);
						}
					}
				}
			}
		}
	}
	m_mesh->unmarkAll<Face>(m_markDone);
	m_mesh->freeMark<Face>(m_markDone);
}
/*----------------------------------------------------------------------------*/
math::Vector3d
SingularityGraphBuilder::getOutputNormal(Face &AFace, Region &ARegion)
{
	std::vector<Node> region_nodes = ARegion.get<Node>();
	std::vector<Node> face_nodes = AFace.get<Node>();

	if (region_nodes.size() != 4) throw GMDSException("SingularityGraphBuilder::getOutputNormal can only be used on tetrahedral regions");
	if (face_nodes.size() != 3) throw GMDSException("SingularityGraphBuilder::getOutputNormal can only be used on triangular faces");

	// we go through all the nodes of ARegion to find the one that do not belong
	// to AFAce
	for (unsigned int i = 0; i < region_nodes.size(); i++) {
		Node n = region_nodes[i];
		if (n != face_nodes[0] && n != face_nodes[1] && n != face_nodes[2]) {
			// n is the node opposite to the face AFace
			Node n0 = face_nodes[0];
			Node n1 = face_nodes[1];
			Node n2 = face_nodes[2];
			math::Vector3d normal_to_face = AFace.normal();
			math::Vector3d in_vector=n.point()-n0.point();
			if (normal_to_face.dot(in_vector) > 0.0) {
				return math::Vector3d({-normal_to_face[0], -normal_to_face[1], -normal_to_face[2]});
				/*return math::Vector3d(-normal_to_face.get(0),
				                              -normal_to_face.get(1),
				                              -normal_to_face.get(2)); */
			}
			else {
				return normal_to_face;
			}
		}
	}
	throw GMDSException("SingularityGraphBuilder::getOutputNormal unexpected behaviour");
}
/*----------------------------------------------------------------------------*/
math::Vector3d
SingularityGraphBuilder::getOutputNormalOfABoundaryFace(Face &AFace)
{
	std::vector<Region> adj_regions = AFace.get<Region>();
	if (adj_regions.size() != 1) throw GMDSException("A boundary face must be adjacent to 1 region!!!");

	return getOutputNormal(AFace, adj_regions[0]);
}
/*----------------------------------------------------------------------------*/
math::Vector3d
SingularityGraphBuilder::getInputNormal(Face &AFace, Region &ARegion)
{
	math::Vector3d outVec = getOutputNormal(AFace, ARegion);
	return outVec.opp();
}
/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::writeOutput(const std::string &AFileName)
{
	static int out = 0;
	std::stringstream file_name;
	file_name << m_output_directory_name << "/" << AFileName << "_" << out;
	writeOutputSingle(file_name.str());
	out++;
}
/*----------------------------------------------------------------------------*/
void
SingularityGraphBuilder::writeOutputSingle(const std::string &AFileName)
{
	m_graph.createVTKOutputFile(AFileName);
	// std::stringstream file_name;
	// file_name << AFileName<<"_FROM_MESH";
	// VTKWriter<IGMesh> writer(*m_mesh);
	gmds::IGMeshIOService ioService(m_mesh);
	gmds::VTKWriter writer(&ioService);
	writer.setCellOptions(gmds::N | gmds::F);
	writer.setDataOptions(gmds::N | gmds::F);
	writer.write("AFileName");
	// writer.write(file_name.str(), DIM3 | R | F | N);
}
/*----------------------------------------------------------------------------*/
math::Quaternion
SingularityGraphBuilder::computeLinearQuaternion(gmds::math::Point &APnt, gmds::Edge &AEdge)
{
	std::vector<gmds::Node> nodes = AEdge.get<Node>();

	Node n1 = nodes[0];
	Node n2 = nodes[1];

	int ID1 = n1.id();
	int ID2 = n2.id();

	math::Point p1 = n1.point();
	math::Point p2 = n2.point();

	// ration pour l'interpolation lineaire
	double edge_length = p1.distance(p2);
	double ratio_length = p1.distance(APnt);
	double param = ratio_length / edge_length;

	std::cout << "Linear Interpolation between quaternions of nodes " << ID1 << "and " << ID2 << ", with param " << param << std::endl;

	math::Quaternion q1 = (*m_var_quatern)[ID1];
	math::Quaternion q2 = (*m_var_quatern)[ID2];
	return math::Quaternion::mean(q1, 1 - param, q2, param);
}
/*----------------------------------------------------------------------------*/
math::Quaternion
SingularityGraphBuilder::buildLinearQuaternion(gmds::Edge &AEdge, const double AParam)
{
	std::vector<gmds::Node> nodes = AEdge.get<Node>();

	Node n1 = nodes[0];
	Node n2 = nodes[1];

	int ID1 = n1.id();
	int ID2 = n2.id();

	math::Quaternion q1 = (*m_var_quatern)[ID1];
	math::Quaternion q2 = (*m_var_quatern)[ID2];

	std::cout << "Linear Interpolation between quaternions of nodes " << ID1 << "and " << ID2 << std::endl;

	return math::Quaternion::mean(q1, 1 - AParam, q2, AParam);
}
/*----------------------------------------------------------------------------*/
