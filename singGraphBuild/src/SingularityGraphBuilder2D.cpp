/*----------------------------------------------------------------------------*/
/*
 * SingularityGraphBuilder2D.cpp
 *
 *  Created on: 13 juil. 2014
 *      Author: bibi
 */
/*----------------------------------------------------------------------------*/
//#include <gmds/ig/IG.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/io/MeditReader.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/math/Cross.h>
#include <gmds/math/Cross2D.h>
#include <gmds/math/AxisAngleRotation.h>
#include <gmds/math/Quaternion.h>
#include <gmds/math/Chart.h>
#include <gmds/math/Ray.h>
#include <gmds/math/Triangle.h>
#include <gmds/math/Segment.h>
#include <gmds/math/Numerics.h>
#include "gmds/math/Constants.h"
#include <gmds/singGraphBuild/SingularityGraphBuilder2D.h>
#include <gmds/singGraphBuild/Tools.h>

#include <gmds/frame/LaplaceCross2D.h>
#include <glpk.h>
/*----------------------------------------------------------------------------*/
#include <sstream>
#include <set>
#include<vector>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace std;
/*----------------------------------------------------------------------------*/
SingularityGraphBuilder2D::SingularityGraphBuilder2D(Mesh*                    AMesh,
                                                     Variable<math::Cross2D>* AField,
                                                     const bool               ABuildGeomSing)
                                                     :m_mesh(AMesh), m_field(AField),
                                                      m_tool(AMesh,AField),
                                                      m_output_directory_name(""),
                                                      m_graph(AMesh)
{
	m_build_geometric_singularities= ABuildGeomSing;

	if(ATolerance<0.01 || ATolerance>0.1)
		throw GMDSException("SingularityGraphBuilder2D: Tolerance must be taken in [0.01,0.1]");
  
	double x_min, y_min, z_min, x_max, y_max, z_max;
    
	math::Point current_pnt = m_mesh->get<Node>(0).getPoint();
	x_min = current_pnt.X();
	x_max = current_pnt.X();
	y_min = current_pnt.Y();
	y_max = current_pnt.Y();
	z_min = current_pnt.Z();
	z_max = current_pnt.Z();
  
	for(auto n_id:m_mesh->nodes()){
		math::Point current_pnt = m_mesh->get<Node>(n_id).getPoint();
		if(current_pnt.X()<x_min)
			x_min = current_pnt.X();
		else  if(current_pnt.X()>x_max)
			x_max = current_pnt.X();

	if(current_pnt.Y()<y_min)
		y_min = current_pnt.Y();
	else  if(current_pnt.Y()>y_max)
		y_max = current_pnt.Y();

	if(current_pnt.Z()<z_min)
		z_min = current_pnt.Z();
	else  if(current_pnt.Z()>z_max)
		z_max = current_pnt.Z();        
	}
	math::Point p_min(x_min, y_min, z_min);
	math::Point p_max(x_max, y_max, z_max);
	m_mesh_radius = p_min.distance(p_max);
	m_confusing_distance = ATolerance*m_mesh_radius;
}

/*----------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::execute(const Strategy AStrategy,  unsigned int& number_of_control_points)
{
	auto t0 = Clock::now();
	std::cout <<"========================================"<< std::endl;
	std::cout<< "Start singularity graph generation "<<std::endl;
	//==================================================================
	// Boolean marks initialization
	//==================================================================
	m_mark_faces_with_sing_point = m_mesh->newMark<Face>();
	m_mark_faces_with_sing_line  = m_mesh->newMark<Face>();

	original_faces_number = m_mesh->getNbFaces();
	original_nodes_number = m_mesh->getNbNodes();
	//==================================================================
	// GEOMETRY VARIABLE FOR DEBUG
	//==================================================================
	m_index = m_mesh->newVariable<int,GMDS_FACE>("index");

/*   Variable<int>* geom_var = m_mesh->newVariable<int,GMDS_NODE>("geometry");


    	for(auto n_id:m_mesh->nodes()){
        Node n = m_mesh->get<Node>(n_id);
        if (m_mesh->isMarked(n, m_mark_nodes_on_point))
            (*geom_var)[n.id()] = 2;
        else if (m_mesh->isMarked(n,m_mark_nodes_on_curve))
            (*geom_var)[n.id()] = 1;
        else
            (*geom_var)[n.id()] = 0;
    }

    gmds::IGMeshIOService meshIoServ(m_mesh);
    gmds::VTKWriter writer(&meshIoServ);
    writer.setCellOptions(gmds::N|gmds::F);
    writer.setDataOptions(gmds::N|gmds::F);
    std::stringstream file_name;
    file_name <<m_output_directory_name<<"-geom_mesh.vtk";
    writer.write(file_name.str());
*/
/*
    	//visualize refVect
	Variable<math::Vector3d>* RefVect = m_mesh->newVariable<math::Vector3d,GMDS_NODE>("refVect");
	for(auto f_id:m_mesh->faces()){
		Face current = m_mesh->get<Face>(f_id);
	  	std::vector<TCellID> nodeIDs = current.getIDs<Node>();
	  	int ID1 = nodeIDs[0];
	  	int ID2 = nodeIDs[1];
	  	int ID3 = nodeIDs[2];

	  	math::Cross2D c1 = (*m_field)[ID1];
	  	math::Cross2D c2 = (*m_field)[ID2];
	  	math::Cross2D c3 = (*m_field)[ID3];

	  	(*RefVect)[ID1] = c1.referenceVector();
	  	(*RefVect)[ID2] = c2.referenceVector();
	  	(*RefVect)[ID3] = c3.referenceVector();
	}

	gmds::IGMeshIOService meshIoServref1(m_mesh);
	gmds::VTKWriter writerref1(&meshIoServref1);
	writerref1.setCellOptions(gmds::N|gmds::F);
	writerref1.setDataOptions(gmds::N|gmds::F);
	std::stringstream file_nameref1;
    	file_nameref1 <<m_output_directory_name<<"-ReferenceVectors.vtk";
	writerref1.write(file_nameref1.str());
   */
    //========================================================================
    // STEP 1 - Detection of singular tetrahedra and storage
    //========================================================================

	unsigned int cont = 0;
	if(AStrategy!=testPrescribedSing){
		/*if the chosen strategy is different than testPrescribedSing, we proceed as follows:
		- first the singular triangles are detected, having as input the frame/cross field for the mesh;
		- for each such singular triangle, the location of the singularity point is detected inside the triangle (for now the location is at the center of the triangle) and also the direction of the slot directions is computed */
		if(withGlobalComments)
			cout<<"AStrategy!=testPrescribedSing "<<endl;
		detectSingularTriangles();

		//visualize crossVect
		Variable<math::Vector3d>* crossVect0 = m_mesh->newVariable<math::Vector3d,GMDS_NODE>("crossVect0");
		Variable<math::Vector3d>* crossVect1 = m_mesh->newVariable<math::Vector3d,GMDS_NODE>("crossVect1");
		Variable<math::Vector3d>* crossVect2 = m_mesh->newVariable<math::Vector3d,GMDS_NODE>("crossVect2");
		Variable<math::Vector3d>* crossVect3 = m_mesh->newVariable<math::Vector3d,GMDS_NODE>("crossVect3");
		for(auto f_id:m_mesh->faces()){
			Face current = m_mesh->get<Face>(f_id);
			std::vector<TCellID> nodeIDs = current.getIDs<Node>();
			int ID1 = nodeIDs[0];
			int ID2 = nodeIDs[1];
			int ID3 = nodeIDs[2];

			std::vector<math::Vector3d> crossV_ID1 = (*m_field)[ID1].componentVectors();
			std::vector<math::Vector3d> crossV_ID2 = (*m_field)[ID2].componentVectors();
			std::vector<math::Vector3d> crossV_ID3 = (*m_field)[ID3].componentVectors();
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

		gmds::IGMeshIOService meshIoServref(m_mesh);
		gmds::VTKWriter writerref(&meshIoServref);
		writerref.setCellOptions(gmds::N|gmds::F);
		writerref.setDataOptions(gmds::N|gmds::F);

		std::stringstream file_nameref;
		file_nameref <<m_output_directory_name<<"-crossVectors.vtk";
		writerref.write(file_nameref.str());

		//========================================================================
		// STEP 2 - Creation of singularity points and slots
		//========================================================================
		// For 3-valent singular points ...
		std::list<gmds::Face>::iterator sing_it;
		for(auto sing_it:m_singularities_3){
			Face currentFace = sing_it;
			createSingPointAndSlots(currentFace, cont);
			cont++;
		}
		// ... and 5-valent singular points
		for(auto sing_it:m_singularities_5){
			Face currentFace = sing_it;
			createSingPointAndSlots(currentFace, cont);
			cont++;
		}

	}
	else{ //just testPrescribedSing
		/*this strategy assumes that here we will impose the singularity triangles by providing their ID: Face singularTri = m_mesh->get<Face>(triangle id);
		this has been implemented with the main scope of testing (ex for a cane-shaped mesh, where the singular triangle is either not detected or detected in a different location that needed for testing;
		note: the detection is highly dependent on the computation of the frame field;
		this strategy could also be used if the user would have available a visualization platform and could simply click on a triangle -> it's ID could be parsed directly)*/
		if(withGlobalComments)
			cout<<"testPrescribedSing"<<testPrescribedSing<<endl;
		unsigned int singularity_type = 3;
		m_singularities_3.clear();
		m_singularities_5.clear();

		Face singularTri = m_mesh->get<Face>(4803);
		m_singularities_3.push_back(singularTri);
		m_mesh->mark(singularTri, m_mark_faces_with_sing_point);
		gmds::math::Point singleCenterTri;
		gmds::math::Cross2D singleCenterTriCross;
		constructOneCenterTriangleCross(singularTri, singleCenterTri, singleCenterTriCross);
		std::vector<math::Point>  slot_points(singularity_type);
		std::vector<double>       slot_param;
		std::vector<int>          slot_cell_dim(singularity_type);
		std::vector<TCellID>      slot_cell_id(singularity_type);
		// SINGULARITY POINT CREATION
		SurfaceSingularityPoint*  singularity = m_graph.newSurfacePoint();
		singularity->setLocation(singleCenterTri);
		singularity->addMeshFace(singularTri);
		std::vector<math::Vector3d> vectors_i = singleCenterTriCross.componentVectors();
		vector<math::Vector3d> vectors_iNew(3);
		vectors_iNew[0] = vectors_i[0];
		vectors_iNew[1] = vectors_i[1];
		vectors_iNew[2] = vectors_i[3];

		math::Point intersectionPnt;
		double intersectionParam;
		for(unsigned int i=0; i<singularity_type; i++){
			math::Ray from_ray(singleCenterTri, vectors_iNew[i]);
			vector<gmds::Edge> currentEdges = singularTri.get<gmds::Edge>();
			for(unsigned int tt=0; tt<3; tt++){
				vector<gmds::Node> currentNodes = currentEdges[tt].get<gmds::Node>();
				math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
				if(from_ray.SecondMetIntersect2D(oppSeg, intersectionPnt, intersectionParam, temp_epsilon)){
					slot_points[i] = intersectionParam * currentNodes[0].getPoint() + (1 - intersectionParam) * currentNodes[1].getPoint();
					slot_cell_dim[i] = 1;
					slot_cell_id[i] = currentEdges[tt].id();
				}
			}
		}

		for(unsigned int i=0;i<singularity_type; i++){
			singularity->newSlot(slot_points[i],
			                     vectors_iNew[i],
			                     slot_cell_id[i],
			                     slot_cell_dim[i],
			                     true,
			                     0);
		}

		Variable<math::Vector3d>* crossVect0 = m_mesh->newVariable<math::Vector3d,GMDS_NODE>("crossVect0");
		Variable<math::Vector3d>* crossVect1 = m_mesh->newVariable<math::Vector3d,GMDS_NODE>("crossVect1");
		Variable<math::Vector3d>* crossVect2 = m_mesh->newVariable<math::Vector3d,GMDS_NODE>("crossVect2");
		Variable<math::Vector3d>* crossVect3 = m_mesh->newVariable<math::Vector3d,GMDS_NODE>("crossVect3");
		for(auto f_id:m_mesh->faces()){
			Face current = m_mesh->get<Face>(f_id);
			std::vector<TCellID> nodeIDs = current.getIDs<Node>();
			int ID1 = nodeIDs[0];
			int ID2 = nodeIDs[1];
			int ID3 = nodeIDs[2];

			std::vector<math::Vector3d> crossV_ID1 = (*m_field)[ID1].componentVectors();
			std::vector<math::Vector3d> crossV_ID2 = (*m_field)[ID2].componentVectors();
			std::vector<math::Vector3d> crossV_ID3 = (*m_field)[ID3].componentVectors();
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

		gmds::IGMeshIOService meshIoServref(m_mesh);
		gmds::VTKWriter writerref(&meshIoServref);
		writerref.setCellOptions(gmds::N|gmds::F);
		writerref.setDataOptions(gmds::N|gmds::F);

		std::stringstream file_nameref;
		file_nameref <<m_output_directory_name<<"-crossVectors.vtk";
		writerref.write(file_nameref.str());
	}

	//========================================================================
	// STEP 3 - Geometrical features are added in the singularity graph
	//========================================================================
	vector<CurveSingularityPoint*> artificialSingPointsCreated;
	/*for each interior boundary loop (if it exists), an artificial curve singularity point is created
	*(at a random location, one one mesh vertex); this artificial singularity point must be removed
	* after we have detected the singularity graph*/

	addGeometryToSingularityGraph(artificialSingPointsCreated);
	if(withGlobalComments){
		cout<<"addGeometryToSingularityGraph"<<endl;
		std::vector<SingularityPoint*> singularity_points1 = m_graph.getPoints();
		cout<<" singularity_points.size() "<<singularity_points1.size()<<endl;
		for (unsigned int i = 0; i < singularity_points1.size(); i++) {
			SingularityPoint* pi = singularity_points1[i];
			std::vector<SingularityPoint::Slot*> pi_slots = pi->getSlots();
			cout<<" sing point "<<i<<" has a total number of slots= "<< pi_slots.size()<<endl;
			for(unsigned int j=0; j<pi_slots.size(); j++){
				cout<<"slot j has direction ";
				cout<<pi_slots[j]->from_point->getLocation()<<" -> "<<(pi_slots[j]->from_point->getLocation()).X()+(pi_slots[j]->direction)[0]<<" , "<<(pi_slots[j]->from_point->getLocation()).Y()+(pi_slots[j]->direction)[1]<<endl;
			}
		}
	}
	vector<bool> singOrGeomFaces(original_faces_number, false);
	cout<<"m_build_geometric_singularities "<<m_build_geometric_singularities<<endl;

	if (m_build_geometric_singularities){
		/*if this flag has been assigned the value "true", geometric singularities and slots will be placed at "concave" corners */
		buildGeometricSlots(singOrGeomFaces);
	}
	//writeOutputSingle("with_geom");

	/*vector<CurveSingularityLine* >   curve_Sing_lines = m_graph.getCurveLines();

	for(unsigned int i=0; i<curve_Sing_lines.size(); i++){
		vector<gmds::math::Point> discrPts = curve_Sing_lines[i]->getDiscretizationPoints();
		if(discrPts.size()>0){
		  	gmds::Mesh meshSing(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
			gmds::Node firstPt = meshSing.newNode(discrPts[0].X(), discrPts[0].Y(), discrPts[0].Z());
			for(unsigned int j=1; j<discrPts.size(); j++){
	 			 gmds::Node nextPt = meshSing.newNode(discrPts[j].X(), discrPts[j].Y(), discrPts[j].Z());
	  			meshSing.newTriangle(firstPt, nextPt, nextPt);
				firstPt = nextPt;
			}

			gmds::IGMeshIOService ioServiceSing(&meshSing);
			gmds::VTKWriter vtkWriterSing(&ioServiceSing);
			vtkWriterSing.setCellOptions(gmds::N|gmds::F);
			vtkWriterSing.setDataOptions(gmds::N|gmds::F);

			std::stringstream vtk_fileSing;
			vtk_fileSing <<m_output_directory_name<<"-curve_Sing_lines" <<i<<".vtk";
			vtkWriterSing.write(vtk_fileSing.str());
		}
	}

	vector<VertexSingularityPoint *> vertex_Sing_points = m_graph.getVertexPoints();
	cout<<"vertex_Sing_points --- Type --- geomType --- number of slots"<<endl;
	for(unsigned int i =0; i<vertex_Sing_points.size(); i++)
		cout<<vertex_Sing_points[i]->getMeshNode().id()<<" "<<vertex_Sing_points[i]->getType()<<" "
		<<vertex_Sing_points[i]->getGeomType()<<" "<<vertex_Sing_points[i]->getSlots().size()<<endl;
	cout<<endl;

	vector<CurveSingularityPoint *> curve_Sing_points = m_graph.getCurvePoints();
	cout<<"curve_Sing_points (location) --- Type --- geomType --- number of slots"<<endl;
    	for(unsigned int i =0; i<curve_Sing_points.size(); i++)
 		cout<<curve_Sing_points[i]->getLocation()<<" "<<curve_Sing_points[i]->getType()<<" "
	 		<<curve_Sing_points[i]->getGeomType()<<" "<<curve_Sing_points[i]->getSlots().size()<<endl;
    	cout<<endl;

    	vector<SurfaceSingularityPoint *> surface_Sing_points = m_graph.getSurfacePoints();
     cout<<"surface_Sing_points (face) --- Type --- geomType --- number of slots"<<endl;
    	for(unsigned int i =0; i<surface_Sing_points.size(); i++)
 		cout<<surface_Sing_points[i]->getMeshFace().id()<<" "<<surface_Sing_points[i]->getType()<<" "
	 		<<surface_Sing_points[i]->getGeomType()<<" "<<surface_Sing_points[i]->getSlots().size()<<endl;
    	cout<<endl; */
    	/* faces consistently oriented
    	math::Vector3d prevNormal(0,0,-1);
    	cout<<"face normals"<<endl;
    	for(auto f_id:m_mesh->faces()){
     	gmds::Face currentFace = m_mesh->get<Face>(f_id);
      	vector<gmds::Node> currentNodes = currentFace.get<Node>();
      	math::Vector3d AB(currentNodes[0].getPoint(),currentNodes[1].getPoint());
      	math::Vector3d AC(currentNodes[0].getPoint(),currentNodes[2].getPoint());
       	math::Vector3d myNormal = AB.cross(AC);
       	myNormal.normalize();
       	if(prevNormal!=myNormal)
	  		cout<<"face "<<f_id<<" has normal "<<myNormal[0]<<","<<myNormal[1]<<","<<myNormal[2]<<endl;
       	prevNormal = myNormal;
    	}*/

	writeSingularityPointsAndSlots();

	vector<vector<gmds::TCellID>> newNode2NewNodeNeighbours;
	vector<vector<double>> face2FaceTransport;
	vector<vector<unsigned int>> face2FaceDeviation;

	constructCenterTriangleCrosses();
	computeFace2FaceInfo(newNode2NewNodeNeighbours, face2FaceTransport, face2FaceDeviation);

	//========================================================================
	// STEP 4 - Initialization of confusing balls around singularity points
	//========================================================================

	std::vector<SingularityPoint*> singularity_points = m_graph.getPoints();
	for (unsigned int i = 0; i < singularity_points.size(); i++) {
		SingularityPoint* pi = singularity_points[i];
		//WARNING if sing points too close it might "overide the singular triangles"
		initConfusingBalls(pi);
	}

	//variable for debug purpose
	Variable<int>* ball_var = m_mesh->newVariable<int,GMDS_FACE>("sing_ball");
	for(auto f_id:m_mesh->faces()){
		Face f = m_mesh->get<Face>(f_id);
		SingularityPoint* sing = m_faces_to_singularity_on_surf[f.id()];
		if (sing==0)
			(*ball_var)[f.id()] = 0;
		else
			(*ball_var)[f.id()] = sing->index();
	}

	gmds::IGMeshIOService meshIoServ(m_mesh);
	gmds::VTKWriter writerB(&meshIoServ);
	writerB.setCellOptions(gmds::N|gmds::F);
	writerB.setDataOptions(gmds::N|gmds::F);
	std::stringstream file_name2;
	file_name2 <<m_output_directory_name<<"-confusing_balls.vtk";
	writerB.write(file_name2.str());
	/*
	* vector<gmds::TCellID> modifiedFaces;
	//get the orginal confusing balls; these are for all singularities, not just for the current one
	vector<gmds::TCellID> originalConfusingBalls;
	for(auto f_id:m_mesh->faces()){
		gmds::Face currentFace = m_mesh->get<Face>(f_id);
		if(m_faces_to_singularity_on_surf[f_id]!=0){
			originalConfusingBalls.push_back(f_id);
		}
	}

	for (unsigned int i = 0; i < singularity_points.size(); i++) {
		SingularityPoint* pi = singularity_points[i];
		redefineOneConfusingBall(pi, modifiedFaces, m_confusing_distance);
	}
	//variable for debug purpose
	Variable<int>* ball_varnew = m_mesh->newVariable<int,GMDS_FACE>("sing_ballNew");
	for(auto f_id:m_mesh->faces()){
		Face f = m_mesh->get<Face>(f_id);
		SingularityPoint* sing = m_faces_to_singularity_on_surf[f.id()];
		if (sing==0)
			(*ball_varnew)[f.id()] = 0;
		else
			(*ball_varnew)[f.id()] = sing->index();
	}

	gmds::VTKWriter writerB2(&meshIoServ);
	writerB2.setCellOptions(gmds::N|gmds::F);
	writerB2.setDataOptions(gmds::N|gmds::F);
	std::stringstream file_name3;
	file_name3 <<m_output_directory_name<<"-redefinedconfusing_balls.vtk";
	writerB2.write(file_name3.str());
	*/
	//========================================================================
	// STEP 5 - Singularity line building
	//========================================================================

	auto t1 = Clock::now();
	cout << "setup before choosing strategy "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
              << " milliseconds" << std::endl;


	mean_edge_length = 0.0;
	for(auto e_id:m_mesh->edges()){
		Edge e = m_mesh->get<Edge>(e_id);
		mean_edge_length = mean_edge_length + e.length();
	}

	mean_edge_length = mean_edge_length/m_mesh->getNbEdges();
	temp_epsilon = 0.03*mean_edge_length;
	if(withGlobalComments){
		cout<<"mean_edge_length "<<mean_edge_length<<endl;
		cout<<"temp_epsilon= "<<temp_epsilon<<endl;
		cout<<"math::Constants::EPSILON "<<math::Constants::EPSILON<<endl;
	}

	if(AStrategy==original)
		createSingularityLines();
	else{
		if(AStrategy==simultaneousStartHeun)
			createSingularityLinesSimultaneousStart();
		else{
			if(AStrategy==simultaneousStartRK4)
				createSingularityLinesSimultaneousStartRK4();
			else{
				if((AStrategy==shortestPaths)||(testPrescribedSing))
					createSingularityLinesShortestPaths(newNode2NewNodeNeighbours, face2FaceTransport, face2FaceDeviation, singOrGeomFaces);
				else
					throw GMDSException("Wrong strategy value for singularity graph computation");
			}
		}
	}
	//========================================================================
	// STEP 6 - Detect singularity lines intersection and split them
	//========================================================================
	cout<<"STEP 6 - Detect singularity lines intersection and split them"<<endl;
	// if(AStrategy==simultaneousStartRK4)
	//writeOutputSingle("boundary_line_semifinalRK4");
  	// else
	// writeOutputSingle("boundary_line_semifinal");

	detectLineIntersections(AStrategy);
  	//  if(AStrategy==simultaneousStartRK4)
	//	writeOutput("boundary_line_intersectionsRK4");
	//else
	// writeOutput("boundary_line_intersections");
	//========================================================================
	// STEP 7 - Build surface patchs
	//========================================================================

	//Delete the artifial nodes that have been created at the begining
	vector<math::Point> firstLineFoundDiscrPts, secondLineFoundDiscrPts;
	vector<gmds::Edge> firstLineMeshEdges, secondLineMeshEdges;
	cont = 0;
	unsigned int firstLineFound;
	bool artificialPointBeginingLine;
	SingularityPoint *theOtherSingPointFirstLine;
	SingularityPoint *theOtherSingPointSecondLine;

	for(unsigned int i=0; i<artificialSingPointsCreated.size(); i++){
		if(withGlobalComments)
			cout<<(artificialSingPointsCreated[i]->getLocation()).X()<<" "<<(artificialSingPointsCreated[i]->getLocation()).Y()<<endl;

		vector<CurveSingularityLine*> curveLines = m_graph.getCurveLines();
		for(unsigned int j=0; j<curveLines.size(); j++){
			vector<SingularityPoint* > singPts = curveLines[j]->getEndPoints();
			for(unsigned int t=0; t<singPts.size(); t++){
				if(singPts[t]==artificialSingPointsCreated[i]){
					if(cont==0){
						theOtherSingPointFirstLine = singPts[(t+1)%2];
						firstLineFound = j;
						firstLineFoundDiscrPts = curveLines[j]->getDiscretizationPoints();
						firstLineMeshEdges = curveLines[j]->getMeshEdges();
						if((curveLines[j]->getDiscretizationPoints()[0].X()==artificialSingPointsCreated[i]->getLocation().X()) &&
						(curveLines[j]->getDiscretizationPoints()[0].Y()==artificialSingPointsCreated[i]->getLocation().Y()))
							artificialPointBeginingLine = true;
						else
							artificialPointBeginingLine = false;
					}
					else{
						theOtherSingPointSecondLine = singPts[(t+1)%2];
						secondLineFoundDiscrPts = curveLines[j]->getDiscretizationPoints();

						secondLineMeshEdges = curveLines[j]->getMeshEdges();
						if(artificialPointBeginingLine){
							if((curveLines[j]->getDiscretizationPoints()[0].X()==artificialSingPointsCreated[i]->getLocation().X())&&
							(curveLines[j]->getDiscretizationPoints()[0].Y()==artificialSingPointsCreated[i]->getLocation().Y())){
								reverse(secondLineFoundDiscrPts.begin(), secondLineFoundDiscrPts.end());
								reverse(secondLineMeshEdges.begin(), secondLineMeshEdges.end());
							}
							for(unsigned int t1=1; t1<firstLineFoundDiscrPts.size(); t1++){
								secondLineFoundDiscrPts.push_back(firstLineFoundDiscrPts[t1]);
							}
							for(unsigned int t1=1; t1<firstLineMeshEdges.size(); t1++){
								secondLineMeshEdges.push_back(firstLineMeshEdges[t1]);
							}
							curveLines[firstLineFound]->setMeshEdges(secondLineMeshEdges);
							curveLines[firstLineFound]->setDiscretizationPoints(secondLineFoundDiscrPts);
							curveLines[firstLineFound]->removeAllSingularityPoints();
							curveLines[firstLineFound]->addSingularityPoint(theOtherSingPointSecondLine);
							curveLines[firstLineFound]->addSingularityPoint(theOtherSingPointFirstLine);
						}
						else{
							if((curveLines[j]->getDiscretizationPoints()[0].X()!=artificialSingPointsCreated[i]->getLocation().X())||
							(curveLines[j]->getDiscretizationPoints()[0].Y()!=artificialSingPointsCreated[i]->getLocation().Y())){
								reverse(secondLineFoundDiscrPts.begin(), secondLineFoundDiscrPts.end());
								reverse(secondLineMeshEdges.begin(), secondLineMeshEdges.end());
							}
							for(unsigned int t1=1; t1<secondLineFoundDiscrPts.size(); t1++){
								curveLines[firstLineFound]->addDiscretizationPoint(secondLineFoundDiscrPts[t1]);
							}
							for(unsigned int t1=1; t1<secondLineMeshEdges.size(); t1++){
								firstLineMeshEdges.push_back(secondLineMeshEdges[t1]);
							}

							curveLines[firstLineFound]->setMeshEdges(firstLineMeshEdges);
							curveLines[firstLineFound]->removeSingularityPoint();
							curveLines[firstLineFound]->addSingularityPoint(theOtherSingPointSecondLine);
						}
						vector<SingularityPoint::Slot*> pi_slots = theOtherSingPointSecondLine->getSlots();
						for(unsigned int t2=0; t2< pi_slots.size(); t2++){
							SingularityLine* currentsingline = pi_slots[t2]->line;
							if(currentsingline){
								if(currentsingline->getType()==0){// curve line
									if(pi_slots[t2]->line->getNumber()==curveLines[j]->getNumber()){
										pi_slots[t2]->line = curveLines[firstLineFound];
									}
								}
							}
						}

						m_graph.removeLine(curveLines[j]);
						m_graph.removePoint(artificialSingPointsCreated[i]);
					}
					cont++;
				}
			}
		}
	}

	if(AStrategy==simultaneousStartRK4)
		writeOutputSingle("boundary_line_finalRK4");
	else
		writeOutputSingle("boundary_line_final");

	//m_graph.buildSurfacePatchs();

	m_graph.buildCurveSurfacePatchs(number_of_control_points);

	if(AStrategy==simultaneousStartRK4)
		writeOutputPatches("patchsRK4");
	else
		writeOutputPatches("patchs");
	std::vector<SingularityPatch* >  patchs = m_graph.getSurfacePatchs();

	//========================================================================
	// Boolean marks cleaning
	//========================================================================
	m_mesh->unmarkAll<Face>(m_mark_faces_with_sing_point);
	m_mesh->unmarkAll<Face>(m_mark_faces_with_sing_line);
	m_mesh->freeMark<Face>(m_mark_faces_with_sing_point);
	m_mesh->freeMark<Face>(m_mark_faces_with_sing_line);

	std::cout<<"\t --> Nb generated faces: "<<patchs.size()<<std::endl;
	std::cout <<"========================================"<< std::endl;

}

/*---------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::detectLineIntersections(const Strategy AStrategy)
{
	std::vector<math::Point> added_points;
	std::vector<Face> candidate_faces;
	if(withGlobalComments)
		cout<<"detectLineIntersections "<<endl;
	//========================================================================
	// We detect intersection faces before splitting any curves. After, it will be
	// too late to detect them.
	//=======================================================================
	vector<SurfaceSingularityLine* > surf_lines = m_graph.getSurfaceLines();
	for(auto f_id:m_mesh->faces()){
		Face f = m_mesh->get<Face>(f_id);
		std::vector<SurfaceSingularityLine*> current_lines = getSingularityLinesIn(f);

		if (current_lines.size()>1){
			if(!m_mesh->isMarked(f, m_mark_faces_with_sing_point)){
				candidate_faces.push_back(f);
			}
		}
	} //for (; !it_faces.isDone(); it_faces.next())

	//std::string file_name = "bdyLine";
	// writeOutput(file_name);

	for(unsigned int i=0;i<candidate_faces.size();i++){
		Face f = candidate_faces[i];
		if(withGlobalComments)
			cout<<"f.id() "<<f.id()<<endl;
		std::vector<SurfaceSingularityLine*> current_lines = getSingularityLinesIn(f);

		if(current_lines.size()>=2){
			if(current_lines.size()==2)
				createLineIntersection(AStrategy, current_lines[0],current_lines[1],f, added_points);
			else
				createLineIntersection(AStrategy, current_lines,f, added_points);
		}
	}

}

/*---------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::createLineIntersection(const Strategy                  AStrategy,
                                                       SurfaceSingularityLine*         ALine1,
                                                       SurfaceSingularityLine*         ALine2,
                                                       Face&                           AFace,
                                                       std::vector<gmds::math::Point>& AAddedPoints)
{
	// We look for the geometrical intersection point of the curve lines
	SurfaceSingularityLine *l0 = ALine1;
	SurfaceSingularityLine *l1 = ALine2;
	math::Point p;

	math::Point face_center = AFace.center();
	std::vector<Node> face_nodes = AFace.get<Node>();
	double face_radius = 0;
	unsigned int nb_face_nodes = face_nodes.size();
	for(unsigned int i=0; i < nb_face_nodes; i++) {
		math::Point pi = face_nodes[i].getPoint();
		math::Point pj = face_nodes[(i+1)%nb_face_nodes].getPoint();
		double distance_ij = pi.distance(pj);
		if(distance_ij>face_radius)
			face_radius = distance_ij;
	}

	if(AStrategy==shortestPaths){
		face_radius = 2 * face_radius;
	}

	if(l0->getIntersectionPoint(l1,face_center,face_radius,p)) {
		gmds::math::Triangle ATriangle(face_nodes[0].getPoint(), face_nodes[1].getPoint(), face_nodes[2].getPoint());
		//cout<<"ATriangle.isIn(p) "<<ATriangle.isIn(p)<<endl;
		if(ATriangle.isIn2ndMethod(p)){
			bool already_added = false;
			vector<SurfaceSingularityPoint*> allSingPoints = m_graph.getSurfacePoints();
			for(unsigned int i=0; i<allSingPoints.size();i++){
				//if(math::near(p.distance(allSingPoints[i]->getLocation()),0.0)){
				if(p.distance(allSingPoints[i]->getLocation())<mean_edge_length/2.0){
					already_added = true;
				}
			}

			for(unsigned int i=0;!already_added && i<AAddedPoints.size();i++){
				if(math::near(p.distance(AAddedPoints[i]),0.0))
					already_added = true;
			}

			if(!already_added){
				AAddedPoints.push_back(p);
				//==============================================================
				// Creation of the singularity point
				//==============================================================*
				SurfaceSingularityPoint* new_pnt = m_graph.newSurfacePoint();
				new_pnt->setLocation(p);
				new_pnt->addMeshFace(AFace);
				//==================================================================
				/* Creation of a new sing point and splitting of intersected sing. lines*/
				//==================================================================
				if(withGlobalComments)
					cout<<"line intersection should be at "<<p.X()<<" "<<p.Y()<<endl;
				m_graph.splitSurfaceLine(new_pnt,l0);
				m_graph.splitSurfaceLine(new_pnt,l1);
			}
		}
	}

}

/*---------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::createLineIntersection(const Strategy                         AStrategy,
                                                       std::vector<SurfaceSingularityLine *>& ALines,
                                                       Face&                                  AFace,
                                                       std::vector<gmds::math::Point>&        AAddedPoints)
{
	// We look for the geometrical intersection point of the curve lines
	math::Point p;
	math::Point face_center = AFace.center();
	std::vector<Node> face_nodes = AFace.get<Node>();
	double face_radius =0;
	unsigned int nb_face_nodes = face_nodes.size();
	for(unsigned int i=0; i < nb_face_nodes; i++) {
		math::Point pi = face_nodes[i].getPoint();
		math::Point pj = face_nodes[(i+1)%nb_face_nodes].getPoint();
		double distance_ij = pi.distance(pj);
		if(distance_ij>face_radius)
			face_radius = distance_ij;
		}

		if(AStrategy==shortestPaths){
			face_radius = 2 * face_radius;
		}

		// WARNING THIS DOUBLE LOOP DOES NOT ENSURE TO GET ALL THE POSSIBLE
		// INTERSECTIONS!!!! THIS SHOULD BE IMPROVED.
		for(unsigned int i=0;i<ALines.size()-1; i++){
			SurfaceSingularityLine *li = ALines[i];
			for(unsigned int j=i+1; j<ALines.size(); j++){
			SurfaceSingularityLine *lj = ALines[j];
			if(li->getIntersectionPoint(lj,face_center,face_radius,p)) {
				if(withGlobalComments)
					cout<<"intersection point p "<<p<<endl;
				gmds::math::Triangle ATriangle(face_nodes[0].getPoint(), face_nodes[1].getPoint(), face_nodes[2].getPoint());
				if(ATriangle.isIn2ndMethod(p)){
					bool already_added = false;
					vector<SurfaceSingularityPoint*> allSingPoints = m_graph.getSurfacePoints();
					for(unsigned int i=0; i<allSingPoints.size();i++){
						//if(math::near(p.distance(allSingPoints[i]->getLocation()),0.0)){
						if(p.distance(allSingPoints[i]->getLocation())<mean_edge_length/2.0){
							already_added = true;
						}
					}

					for(unsigned int t=0;!already_added && t<AAddedPoints.size();t++){
						if(math::near(p.distance(AAddedPoints[t]),0.0))
							already_added = true;
					}

					if(!already_added){
						AAddedPoints.push_back(p);
						//==============================================================
						// Creation of the singularity point
						//==============================================================*
						SurfaceSingularityPoint* new_pnt = m_graph.newSurfacePoint();
						new_pnt->setLocation(p);
						new_pnt->addMeshFace(AFace);
						//==================================================================
						// Creation of a new sing point
						// And splitting of intersected sing. lines
						//==================================================================
						m_graph.splitSurfaceLine(new_pnt,li);
						m_graph.splitSurfaceLine(new_pnt,lj);
						if(withGlobalComments)
							cout<<"line intersection should be at "<<p.X()<<" "<<p.Y()<<endl;
					}
				}
			}
		}
	}

}

/*---------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::detectSingularTriangles()
{
	m_singularities_3.clear();
	m_singularities_5.clear();

	for(auto f_id:m_mesh->faces()){
		Face current = m_mesh->get<Face>(f_id);
		std::vector<TCellID> nodeIDs = current.getIDs<Node>();
		int ID1 = nodeIDs[0];
		int ID2 = nodeIDs[1];
		int ID3 = nodeIDs[2];

		math::Cross2D c1 = (*m_field)[ID1];
		math::Cross2D c2 = (*m_field)[ID2];
		math::Cross2D c3 = (*m_field)[ID3];

		int index = math::Cross2D::index(c1,c2,c3);

		(*m_index)[current.id()] = index;

		if (index ==1) {
			m_singularities_5.push_back(current);
			m_mesh->mark(current, m_mark_faces_with_sing_point);
		}
		else if (index == -1) {
			m_singularities_3.push_back(current);
			m_mesh->mark(current, m_mark_faces_with_sing_point);
		}
	} //for (; !it_regions.isDone(); it_regions.next())

}

/*---------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::initMarks(const int AMarkNodePnt,
                                          const int AMarkNodeCrv,
                                          const int AMarkEdgeCrv)
{
	m_mark_nodes_on_point = AMarkNodePnt;
	m_mark_nodes_on_curve = AMarkNodeCrv;
	m_mark_edges_on_curve = AMarkEdgeCrv;
}

/*----------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::initConfusingBalls(SingularityPoint* APnt)
{
	math::Point sing_location = APnt->getLocation();
	bool found_face = false;
	//========================================================================
	// Faces
	//========================================================================
	// int nb_faces = 0;
	double m_confusing_distance_temp = 0.0;
	while(!found_face){
		m_confusing_distance_temp = m_confusing_distance_temp + m_confusing_distance;
		for(auto f_id:m_mesh->faces()){
			Face currentFace = m_mesh->get<Face>(f_id);
			math::Point center = currentFace.center();
			if(center.distance(sing_location)< m_confusing_distance_temp) {
				m_faces_to_singularity_on_surf[currentFace.id()] = APnt;
				found_face = true;
				//nb_faces++;
			}
		}
	}
	//========================================================================
	// Edges
	//========================================================================
	for(auto e_id:m_mesh->edges()){
		Edge currentEdge = m_mesh->get<Edge>(e_id);
		math::Point center = currentEdge.center();
		if(center.distance(sing_location)< m_confusing_distance) {
			m_edges_to_singularity_on_surf[currentEdge.id()] = APnt;
		}
	}
	//========================================================================
	// Nodes
	//========================================================================
	for(auto n_id:m_mesh->nodes()){
		Node currentNode = m_mesh->get<Node>(n_id);
		math::Point center = currentNode.getPoint();
		if(center.distance(sing_location)< m_confusing_distance) {
			m_nodes_to_singularity_on_surf[currentNode.id()] = APnt;
			vector<gmds::Face> adjacent_triangles = currentNode.get<gmds::Face>();
			for(unsigned int i=0; i<adjacent_triangles.size(); i++){
				math::Point centerFace = adjacent_triangles[i].center();
				//if(centerFace.distance(sing_location)< m_confusing_distance_temp) {
				m_faces_to_singularity_on_surf[adjacent_triangles[i].id()] = APnt;
				//}
			}
		}
	}

}

/*----------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::computeSingPointInfo(Face&        AFace,
                                                     math::Point& APosSing)
{
	std::vector<Node> currentNodes = AFace.get<Node>();
	Node node1 = currentNodes[0];
	Node node2 = currentNodes[1];
	Node node3 = currentNodes[2];

	math::Vector3d v1 =(*m_field)[node1.id()].referenceVector();
	math::Vector3d v2 =(*m_field)[node2.id()].referenceVector();
	math::Vector3d v3 =(*m_field)[node3.id()].referenceVector();

	math::Point pointA = currentNodes[0].getPoint();
	math::Point pointB = currentNodes[1].getPoint();
	math::Point pointC = currentNodes[2].getPoint();

	//We solve a 2x2 Ax=b system with x = (alpha,beta)
	double alpha=0, beta=0, gamma=0;
	math::Vector3d A = v1-v3;
	math::Vector3d B = v2-v3;
	math::Vector3d C = v3.opp();

	double dA = A[0]*B[1]-A[1]*C[0];
	if(dA!=0){
		double Dx = C[0]*B[1]-C[1]*B[0];
		double Dy = A[0]*C[1]-A[1]*C[0];
		alpha = Dx/dA;
		beta  = Dy/dA;
		gamma = 1-alpha-beta;
	}
	else{
		throw GMDSException("Null Determinant in the computation of a singularity point location");
	}
	// WARNING: TODO the singularity point will always be placed at the traingle center...
	if(true){//dA<0){
		alpha = 0.333;
		beta  = 0.333;
		gamma = 1-alpha-beta;
	}
	//(alpha, beta, gamma) are the barycentric coordinates where the cross field vanishes
	//We can then compute the singularity point location
	double x = alpha * node1.X() + beta * node2.X() + gamma * node3.X();
	double y = alpha * node1.Y() + beta * node2.Y() + gamma * node3.Y();

	APosSing.setXYZ(x,y,0);

}

/*----------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::createSingPointAndSlots(Face& AFace, unsigned int& cont)
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
	//for each detected slot, we store its location, its direction, its out cell
	// id and the dimension of this out cell.
	std::vector<math::Point>  slot_points;
	std::vector<double>       slot_param;
	std::vector<math::Vector3d> slot_dirs;
	std::vector<int>          slot_cell_dim;
	std::vector<TCellID>      slot_cell_id;
	std::vector<Edge> edges = AFace.get<Edge>();
	for(int i=0; i<3; i++){ //we walk along each edge of AFace
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

		for(int k=0;k<2;k++){
			math::Vector3d cik =  vectors_i[k];
			math::Vector3d cjk =  cross_j.closestComponentVector(cik);

			double x_cik = cik.X();
			double y_cik = cik.Y();
			double x_cjk = cjk.X();
			double y_cjk = cjk.Y();

			double x_cijk = x_cik - x_cjk;
			double y_cijk = y_cik - y_cjk;

			double a = (xji* y_cijk) - (yji   * x_cijk);
			double b = (y_cjk * xji) + (xsj   * y_cijk) - (x_cjk * yji) - (ysj * x_cijk);
			double c = (y_cjk * xsj) - (x_cjk * ysj);

			std::vector<double> solutions;
			math::solve2ndDegreePolynomial(a,b,c,solutions);
			for(unsigned int i_sol=0; i_sol<solutions.size();i_sol++){
				double alpha = solutions[i_sol];
				if(alpha>1 || alpha<0.0)
					continue;
				math::Point p = alpha*pi + (1-alpha)*pj;

				slot_param.push_back(alpha);
				slot_points.push_back(p);

				slot_dirs.push_back(math::Vector3d(s,p));

				if(alpha==0){ //we go out from a node
					slot_cell_dim.push_back(0);
					slot_cell_id.push_back(ni.id());
				}
				else if(alpha==0){ //we go out from a node
					slot_cell_dim.push_back(0);
					slot_cell_id.push_back(nj.id());
				}
				else{ // general case, the edge
					slot_cell_dim.push_back(1);
					slot_cell_id.push_back(ei.id());
				}
			}
		}//for(int k=0;k<2;k++)
	}// for(int i=0;i<3;i++)

	// SINGULARITY POINT CREATION
	SurfaceSingularityPoint*  singularity = m_graph.newSurfacePoint();
	singularity->setLocation(s);
	singularity->addMeshFace(AFace);
	for(unsigned int i=0;i<slot_points.size(); i++){
		singularity->newSlot(slot_points[i],
		                     slot_dirs[i],
		                     slot_cell_id[i],
		                     slot_cell_dim[i],
		                     true,
		                     0);
	}

}

/*----------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::addGeometryToSingularityGraph(vector<CurveSingularityPoint*> &artificialSingPointsCreated) {
	//Now we add singularity points and lines from geometric corners
	//============================================================================
	//	GEOMETRIC  SINGULARITY POINTS
	//============================================================================
	for(auto n_id:m_mesh->nodes()){
		Node currentNode = m_mesh->get<Node>(n_id);
		if (m_mesh->isMarked(currentNode, m_mark_nodes_on_point)) {
			VertexSingularityPoint *sing_point = m_graph.newVertexPoint();
			sing_point->setXYZ(currentNode.X(), currentNode.Y(), currentNode.Z());
			sing_point->addMeshNode(currentNode);
		}
	}
	//============================================================================
	//	GEOMETRIC  SINGULARITY LINES
	//============================================================================
	//Now we have all the geometric corners viewed as singularity points
	//We create separatrices based on the geometric edges

	int mark_geom_edges = m_mesh->newMark<Node>();
	std::vector<CurveSingularityLine *> added_geom_lines;

	for(auto n_id:m_mesh->nodes()){
		Node currentNode = m_mesh->get<Node>(n_id);

		if ((!m_mesh->isMarked(currentNode, m_mark_nodes_on_point)) &&
		(m_mesh->isMarked(currentNode, m_mark_nodes_on_curve)) &&
		(!m_mesh->isMarked(currentNode, mark_geom_edges))) {
			//new singularity line to create here
			std::vector<int> listOfNodesInSingLeft;
			std::vector<int> listOfNodesInSingRight;
			m_mesh->mark(currentNode, mark_geom_edges);
			Node NodeLeft;
			Node NodeRight;
			int gotFirstOne = 0;
			for (unsigned int i = 0; i < currentNode.get<Edge>().size(); i++) {
				if (m_mesh->isMarked(currentNode.get<Edge>()[i], m_mark_edges_on_curve)) {
					for (unsigned int j = 0; j < 2; j++) {
						if (currentNode.get<Edge>()[i].get<Node>()[j].id() != currentNode.id()) {
							if (!gotFirstOne) {
								gotFirstOne = 1;
								NodeLeft = currentNode.get<Edge>()[i].get<Node>()[j];
							}
							else {
								NodeRight = currentNode.get<Edge>()[i].get<Node>()[j];
							}
						}
					}
				}
			}
			Node Ncurr = currentNode;
			//     int numberOfTheFirstSing = 0;
			//     int numberOfTheSecondSing = 0;
			//Here we have the two nodes left and right that continue the line
			//------------------------------------------------------------------------------------
			//la premiere condition de la boucle while empeche la fermeture des cercles
			while ((NodeLeft.id() != NodeRight.id()) &&
			(!m_mesh->isMarked(NodeLeft, m_mark_nodes_on_point))) {
				m_mesh->mark(NodeLeft, mark_geom_edges);
				listOfNodesInSingLeft.push_back(NodeLeft.id());

				Node NodeNext = NodeLeft;
				for (unsigned int i = 0; i < NodeLeft.get<Edge>().size(); i++) {
					if (m_mesh->isMarked(NodeLeft.get<Edge>()[i], m_mark_edges_on_curve)) {
						for (unsigned int j = 0; j < NodeLeft.get<Edge>()[i].get<Node>().size(); j++) {
							if (NodeLeft.get<Edge>()[i].get<Node>()[j].id() != NodeLeft.id()) {
								if (NodeLeft.get<Edge>()[i].get<Node>()[j].id() != Ncurr.id()) {
									NodeNext = NodeLeft.get<Edge>()[i].get<Node>()[j];
								}
							}
						}
					}
				}
				Ncurr = NodeLeft;
				NodeLeft = NodeNext;
			}
			//------------------------------------------------------------------------------------
			m_mesh->mark(NodeLeft, mark_geom_edges);// currently this could be marked as node on point
			listOfNodesInSingLeft.push_back(NodeLeft.id());
			if (NodeLeft.id() != NodeRight.id()) { //we need to get on the right too
				Node Ncurr = currentNode;
				while (!m_mesh->isMarked(NodeRight, m_mark_nodes_on_point)) {
					m_mesh->mark(NodeRight, mark_geom_edges);
					listOfNodesInSingRight.push_back(NodeRight.id());

					Node NodeNext = NodeRight;
					for (unsigned int i = 0; i < NodeRight.get<Edge>().size(); i++) {
						if (m_mesh->isMarked(NodeRight.get<Edge>()[i], m_mark_edges_on_curve)) {
							for (unsigned int j = 0; j < 2; j++) {
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
				m_mesh->mark(NodeRight, mark_geom_edges);
				listOfNodesInSingRight.push_back(NodeRight.id());
			}
			else { //cycle
				listOfNodesInSingRight.push_back(NodeLeft.id());
			}
			CurveSingularityLine *new_line = m_graph.newCurveLine();
			new_line->setNumber(m_graph.getNbLines()); //redundant now

			std::vector<Node> curve_nodes;
			for (unsigned int i = 0; i < listOfNodesInSingLeft.size(); i++) {
				Node NAtThisPoint = m_mesh->get<Node>(listOfNodesInSingLeft[listOfNodesInSingLeft.size() - 1 - i]);

				new_line->addDiscretizationPoint(NAtThisPoint.getPoint());
				curve_nodes.push_back(NAtThisPoint);
			}

			new_line->addDiscretizationPoint(currentNode.getPoint());
			curve_nodes.push_back(currentNode);

			for (unsigned int i = 0; i < listOfNodesInSingRight.size(); i++) {
				Node NAtThisPoint = m_mesh->get<Node>(listOfNodesInSingRight[i]);

				new_line->addDiscretizationPoint(NAtThisPoint.getPoint());
				curve_nodes.push_back(NAtThisPoint);
			}
			std::vector<Edge> curve_edges;
			for (unsigned int i = 0; i < curve_nodes.size() - 1; i++) {
				Node current = curve_nodes[i];
				Node next = curve_nodes[i + 1];
				std::vector<Edge> currentEdges = current.get<Edge>();

				bool found_edge = false;
				for (unsigned int j = 0; j < currentEdges.size() && !found_edge; j++) {
					Edge ej = currentEdges[j];
					if (m_mesh->isMarked(ej, m_mark_edges_on_curve)) {
						std::vector<Node> ej_nodes = ej.get<Node>();
						if (ej_nodes[0].id() == next.id() || ej_nodes[1].id() == next.id()) {
							curve_edges.push_back(ej);
							found_edge = true;
						}
					}
				}
			}

			//TODO Attention aux courbes cycliques, a priori non traitees
			new_line->setMeshEdges(curve_edges);
			added_geom_lines.push_back(new_line);
		}
	}

	// Geometrical curves made of only one mesh edge are missing. We add them now.

	for(auto e_id:m_mesh->edges()){
		Edge currentEdge = m_mesh->get<Edge>(e_id);
		// on ne traite que les aretes sur une courbe geometrique
		if (m_mesh->isMarked(currentEdge, m_mark_edges_on_curve)) {
			std::vector<Node> currentNodes = currentEdge.get<Node>();
			Node n1 = currentNodes[0];
			Node n2 = currentNodes[1];
			// on regarde si les deux noeuds correspondent ˆ des sommets geometriques
			if (m_mesh->isMarked(n1, m_mark_nodes_on_point) &&
			m_mesh->isMarked(n2, m_mark_nodes_on_point)) { // On cree donc la separatrice reliant les singularites associŽes ˆ n1 et n2

				CurveSingularityLine *new_line = m_graph.newCurveLine();
				new_line->setNumber(m_graph.getNbLines());
				std::vector<Edge> edges;
				edges.push_back(currentEdge);
				new_line->setMeshEdges(edges);
				//new_line->setNumber(1);
				new_line->addDiscretizationPoint(n1.getPoint());
				new_line->addDiscretizationPoint(n2.getPoint());
				added_geom_lines.push_back(new_line);
			}
		}
	}
	//============================================================================
	// GEOM SINGULARITY POINTS AND GEOM SINGULARITY LINES MUST BE CONNECTED
	//============================================================================
	std::vector<VertexSingularityPoint *> geom_points = m_graph.getVertexPoints();
	for (unsigned int i = 0; i < geom_points.size(); i++) {
		VertexSingularityPoint *pi = geom_points[i];
		Node ni = pi->getMeshNode();
		for (unsigned int j = 0; j < added_geom_lines.size(); j++) {
			CurveSingularityLine *lj = added_geom_lines[j];
			std::vector<Edge>& lj_edges= lj->getMeshEdges();
			Node first_node, last_node;
			math::Vector3d first_vec, last_vec;
			if(lj_edges.size()==1){
				std::vector<Node> e_nodes = lj_edges[0].get<Node>();
				first_node = e_nodes[0];
				last_node = e_nodes[1];
				first_vec = math::Vector3d(first_node.getPoint(), last_node.getPoint());
				last_vec  = math::Vector3d(last_node.getPoint() , first_node.getPoint());
			}
			else {
				Edge first_edge = lj_edges[0];
				Edge second_edge = lj_edges[1];
				std::vector<Node> ref_nodes = first_edge.get<Node>();
				std::vector<Node> next_nodes = second_edge.get<Node>();
				Node forbidden_node;
				Node next_node;
				for(auto ni:ref_nodes){
					for (auto nj:next_nodes){
						if(ni==nj){
							forbidden_node = ni;
						}
					}
				}
				for(auto ni:ref_nodes){
					if(ni!=forbidden_node){
						first_node= ni;
					}
					else{
						next_node=ni;
					}
				}

				first_vec = math::Vector3d(first_node.getPoint(), next_node.getPoint());
				Edge last_but_one_edge = lj_edges[lj_edges.size()-2];
				Edge last_edge = lj_edges[lj_edges.size()-1];

				next_nodes = last_but_one_edge.get<Node>();
				ref_nodes = last_edge.get<Node>();
				for(auto ni:ref_nodes){
					for (auto nj:next_nodes){
						if(ni==nj){
							forbidden_node = ni;
						}
					}
				}
				for(auto ni:ref_nodes){
					if(ni!=forbidden_node){
						last_node= ni;
					}
					else{
						next_node=ni;
					}
				}

				last_vec = math::Vector3d(last_node.getPoint(), next_node.getPoint());

			}
			//We compare the point node with the first and last node of the curve
			if(ni==first_node){
				//pi and lj must be connected
				//WARNING SLOT DIMENSION IS DETERMINED BY THE LAST EDGE DIRECTION OF THE CURVE LINE STARTING FROM THIS POINT;
				// PERHAPS THE ORIENTATION OF THE TOTAL CURVE LINE SHOULD BE TAKEN INTO ACC
				//WARNING HERE WAS LAST_VEC INSTEAD OF FIRST_VEC
				SingularityPoint::Slot *s = pi->newSlot(ni.getPoint(),
				                                        first_vec,
				                                        ni.id() /*starting cell id*/,
				                                        0 /*starting cell dim*/,
				                                        true /*on surface*/,
				                                        lj);
				s->isLaunched = true;
				s->isFreeze = true;

				vector <gmds::Edge>myEdges=lj->getMeshEdges();

				for(unsigned int tt=0;tt<myEdges.size();tt++){
					vector<Node> e_nodes = myEdges[tt].get<Node>();
				}

				lj->addSingularityPoint(pi);
			}
			else if(ni==last_node){
				//pi and lj must be connected
				//pi->newSlot(ni.getPoint(), last_vec,
				SingularityPoint::Slot *s = pi->newSlot(ni.getPoint(),
				                                        last_vec,
				                                        ni.id() /*starting cell id*/,
				                                        0 /*starting cell dim*/,
				                                        true /*on surface*/,
				                                        lj);
				s->isLaunched = true;
				s->isFreeze = true;

				vector <gmds::Edge>myEdges=lj->getMeshEdges();

				for(unsigned int tt=0;tt<myEdges.size();tt++){
					vector<Node> e_nodes = myEdges[tt].get<Node>();
				}

				lj->addSingularityPoint(pi);
			}
		}
	}
	cout<<"CYCLE GEOMETRIC CURVE WILL NOT BE CONNECTED TO GEOM POINT. WE CREATE SUCH"<<endl;
	//============================================================================
	// CYCLE GEOMETRIC CURVE WILL NOT BE CONNECTED TO GEOM POINT. WE CREATE SUCH
	// POINT IN AN ARTIFICAL WAY
	//============================================================================

	std::vector<CurveSingularityLine *> geom_curves = m_graph.getCurveLines();
	for (unsigned int i = 0; i < geom_curves.size(); i++) {
		CurveSingularityLine *ci = geom_curves[i];
		if (ci->getEndPoints().empty()) {
			gmds::math::Point loc = ci->getDiscretizationPoints()[0];
			CurveSingularityPoint *p = m_graph.newCurvePoint();
			artificialSingPointsCreated.push_back(p);
			p->setLocation(loc);
			std::vector<math::Point> line_disc = ci->getDiscretizationPoints();
			math::Point p_begin = line_disc[0];
			math::Point p_end = line_disc[line_disc.size() - 1];

			math::Vector3d slot_dir1 = math::Vector3d(p_begin, line_disc[1]);
			math::Vector3d slot_dir2 = math::Vector3d(p_end, line_disc[line_disc.size() - 2]);

			SingularityPoint::Slot *s = p->newSlot(p->getLocation(), // slot location
			                                       slot_dir1, // slot direction
			                                       0,    // No linked cell (id)
			                                       0,    // No linked cell (dim)
			                                       true, // Always on surface (Maybe false in the future)
			                                       ci,// Connected line
			                                       slot_dir1.opp());// Line direction is the slot direction
			s->isFreeze = true;
			s = p->newSlot(p->getLocation(), // slot location
			               slot_dir2, // slot direction
			               0,    // No linked cell (id)
			               0,    // No linked cell (dim)
			               true, // Always on surface (Maybe false in the future)
			               ci,// Connected line
			               slot_dir2.opp());// Line direction is the slot direction
			s->isFreeze = true;
			ci->addSingularityPoint(p);
			ci->addSingularityPoint(p);
		}
	}
	m_mesh->unmarkAll<Node>(mark_geom_edges);
	m_mesh->freeMark<Node>(mark_geom_edges);

}

/*----------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::buildGeometricSlots(vector<bool>& singOrGeomFaces)
{
	//================================================================================
	//WE TRAVERSE GEOM SING POINTS (ASSIGNED TO GEOMETRIC CORNERS) TO DEFINE THEIR
	// SLOTS AT NON-CONVEX AREAS
	//================================================================================
	std::vector<VertexSingularityPoint *> vertex_points = m_graph.getVertexPoints();
	for (unsigned int i = 0; i < vertex_points.size(); i++) {
		VertexSingularityPoint *current_point = vertex_points[i];
		Node currentNode = current_point->getMeshNode();
		math::Point in_pnt = currentNode.getPoint();
		//=========================================================
		// First, we compute the solid angle around currentNode
		std::vector<Face> currentFaces = currentNode.get<Face>();
		std::vector<Edge> currentEdges = currentNode.get<Edge>();
		if(withGlobalComments)
			cout<<"sing_point "<<i<<endl;
		double angle_rad=0;
		for (auto cur_face:currentFaces) {
			singOrGeomFaces[cur_face.id()] = true;
			std::vector<Node> face_nodes = cur_face.get<Node>();
			Node other_nodes[2];
			int i_on = 0;
			for (auto ni:face_nodes) {
				if (ni != currentNode) {
					other_nodes[i_on++] = ni;
				}
			}

			math::Vector3d v1(in_pnt, other_nodes[0].getPoint());
			math::Vector3d v2(in_pnt, other_nodes[1].getPoint());
			angle_rad += v1.angle(v2);
			if(withGlobalComments){
				cout<<"in_pnt "<<in_pnt<<endl;
				cout<<"other_nodes[0].getPoint() "<<other_nodes[0].getPoint()<<endl;
				cout<<"other_nodes[1].getPoint() "<<other_nodes[1].getPoint()<<endl;
			}
		}
		double angle_deg = angle_rad * 180 / math::Constants::PI;
		int nb_lines = 0;
		double single_line_angle = 0;
		if(withGlobalComments)
			cout<<"angle_deg "<<angle_deg<<endl;

		if (angle_deg > 285) {
			//THIS CASE IS NOT WELL IMPLEMENTED
			//Reverse configuration
			//3 lines to throw
			nb_lines = 3;
			single_line_angle = angle_rad / (nb_lines+1);
			createLineFrom(current_point,single_line_angle,3);
		}
		else if (angle_deg > 255) {
			//2 lines to throw
			nb_lines = 2;
			single_line_angle = angle_rad / (nb_lines+1);
			createLineFrom(current_point,single_line_angle,2);
		}
	}
}

/*----------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::createLineFrom(VertexSingularityPoint* AFrom,
                                               const double            AAngle,
                                               const int               ANbLines)
{
	if(withGlobalComments)
		cout<<"createLineFrom, ANbLines "<<ANbLines<<endl;
	Node currentNode = AFrom->getMeshNode();
	std::vector<Face> currentFaces = currentNode.get<Face>();
	std::vector<Edge> currentEdges = currentNode.get<Edge>();
	//We get a first edge to start from
	Edge first_edge;

	for(unsigned int e_cont=0; e_cont<currentEdges.size(); e_cont++){
		gmds::Edge ei = currentEdges[e_cont];
		if(m_mesh->isMarked(ei,m_mark_edges_on_curve)){
			first_edge = ei;
			break;
		}
	}
	for(auto i_line=0; i_line<ANbLines; i_line++){
		double line_angle = AAngle*(i_line+1);
		bool found_out = false;

		Face first_face = first_edge.get<Face>()[0];
		std::vector<Node> first_face_nodes = first_face.get<Node>();
		Node from_n, bnd_n, to_n;
		vector<unsigned int> contNode(2,0);//from_n, bnd_n, to_n

		for(unsigned int ni=0; ni<3; ni++){
			if (first_face_nodes[ni]==currentNode){
				from_n = first_face_nodes[ni];
				contNode[0] = ni;
			}
			else{
				if((m_mesh->isMarked(first_face_nodes[ni],m_mark_nodes_on_curve))||(m_mesh->isMarked(first_face_nodes[ni],m_mark_nodes_on_point))){
					/* for curve reduce to one mesh edge*/
					bnd_n = first_face_nodes[ni];
					contNode[1] = ni;
				}
			}
		}

		for(auto ni:first_face_nodes){
			if (ni!=from_n && ni!=bnd_n) {
				to_n = ni;
			}
		}

		math::Vector3d v_bnd(from_n.getPoint(), bnd_n.getPoint());
		math::Vector3d v_to (from_n.getPoint(), to_n.getPoint() );

		v_bnd.normalize();
		v_to.normalize();
		math::Vector3d axis = v_bnd.cross(v_to);//surf normal
		axis.normalize();

		math::AxisAngleRotation R(axis,line_angle);
		math::Vector3d line_dir = R*v_bnd;
		math::Vector3d vec(line_dir.X(), line_dir.Y(), line_dir.Z());
		vec.normalize();
		//cout<<"vec "<<vec<<endl;
		math::Point op = currentNode.getPoint()+vec;
		//cout<<"op "<<op<<endl;
		// Now for each face adjacent to currentNode, we look for one intersected
		// by vec

		for (unsigned int k = 0; k < currentFaces.size() && !found_out; k++){
			Face fk = currentFaces[k];
			Node node_from = currentNode;
			//=====================================================================
			// We look for the opposite edge
			//=====================================================================
			Edge opposite_edge;
			std::vector<Edge> fk_edges = fk.get<Edge>();

			for (unsigned int i = 0; i < fk_edges.size(); i++) {
				Edge edgei = fk_edges[i];
				std::vector<Node> ei_nodes = edgei.get<Node>();
				if (node_from != ei_nodes[0] && node_from != ei_nodes[1]){
					opposite_edge = edgei;
					//cout<<"node_from "<<node_from<<endl;
					//cout<<"opposite_edge "<<ei_nodes[0]<<" "<<ei_nodes[1]<<endl;
				}
			}

			//=====================================================================
			// We compute an out pnt and vector
			//=====================================================================
			math::Point  out_pnt;
			math::Vector3d out_vec;
			int out_dim =0;
			TCellID out_id = NullID;
			math::Vector3d in_vec = vec;
			std::vector<Node> other_nodes = opposite_edge.get<Node>();
			//===========================================================
			// Go through the first opposite node if it is not on a curve ?
			//===========================================================
		 	if(!m_mesh->isMarked(other_nodes[0],m_mark_nodes_on_curve)) {

				math::Point opp_node_loc1 = other_nodes[0].getPoint();
				math::Vector3d v_opp1(currentNode.getPoint(), opp_node_loc1);
				v_opp1.normalize();
				if (math::near(v_opp1.dot(in_vec) - 1,0)) {
					found_out = true;
					out_pnt = opp_node_loc1;
					m_tool.computeOutVectorAtPoint(other_nodes[0], in_vec, out_vec);
					out_dim = 0;
					out_id  = other_nodes[0].id();
				}
			}
			if(!found_out && !m_mesh->isMarked(other_nodes[1],m_mark_nodes_on_curve)) {

				math::Point opp_node_loc2 = other_nodes[1].getPoint();
				math::Vector3d v_opp2(currentNode.getPoint(), opp_node_loc2);
				v_opp2.normalize();

				if (math::near(v_opp2.dot(in_vec) - 1,0)) {
					found_out = true;
					out_pnt = opp_node_loc2;
					m_tool.computeOutVectorAtPoint(other_nodes[1], in_vec, out_vec);
					out_dim = 0;
					out_id  = other_nodes[1].id();
				}
			}

			//================================================
			// Go through the opposite edge
			// And not through one of its end points due to
			// previous tests.
			//================================================
			if(!found_out) {
				double deviation = 0;
				found_out = m_tool.computeOutVectorFromRayAndEdge(opposite_edge,
				                                                  currentNode.getPoint(),
				                                                  in_vec,
				                                                  out_pnt,
				                                                  out_vec,
				                                                  deviation);
				if(found_out){
					out_dim = 1;
					out_id  = opposite_edge.id();
				}
			}

			if(found_out){
				//cout<<"fk.id() "<<fk.id()<<endl;
				// problem when inside circular hole
				//cout<<"slot dir "<<currentNode.getPoint()<<" -> "<<out_pnt<<endl;
				//cout<<"slot dir "<<currentNode.getPoint()<<" -> "<<currentNode.getPoint().X()+out_vec[0]<<" , "<<currentNode.getPoint().Y()+out_vec[1]<<endl;

				for(unsigned int e_contTemp=0; e_contTemp<currentEdges.size(); e_contTemp++){
					gmds::Edge eiTemp = currentEdges[e_contTemp];
					if(m_mesh->isMarked(eiTemp,m_mark_edges_on_curve)){
						std::vector<Node> currentNodes_mine = eiTemp.get<Node>();

						Node from_nTemp, bnd_nTemp;
						for(unsigned int ni=0; ni<2; ni++){
							if (currentNodes_mine[ni]==currentNode){
								from_nTemp = currentNodes_mine[ni];
							}
							else{
							bnd_nTemp = currentNodes_mine[ni];
							}
						}

						//cout<<"v_bndTemp "<<from_nTemp.getPoint()<<" -> "<<bnd_nTemp.getPoint()<<endl;
						math::Vector3d v_bndTemp(from_nTemp.getPoint(), bnd_nTemp.getPoint());
						//cout<<"out_vec.angle(v_bndTemp) "<<out_vec.angle(v_bndTemp)<<endl;
					}
				}


				AFrom->newSlot(out_pnt, out_vec, out_id, out_dim, false);
				break;// Circle_with_concavities_coarse
			}
		}
	}

}

/*----------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::createSingularityLines()
{
	//Now we need to create separatrices on the faces
	//we will start by launching them from the slot of the singularities
	//for each one we have the pos and dir to launch the sep into, as well as the first triangle

	//========================================================================
	// At the beginning, only the singularity points of the cross field are
	// known. Some other points may be created, when we intersect
	// boundary, singularity lines, etc.
	//========================================================================
	std::vector<SingularityPoint*> singularity_points = m_graph.getPoints();

	for (unsigned int i = 0; i < singularity_points.size(); i++) {
		SingularityPoint* pi = singularity_points[i];
		std::vector<SingularityPoint::Slot*> pi_slots = pi->getSlots();

		for (unsigned int j = 0; j < pi_slots.size(); j++) {
			if (!pi_slots[j]->isLaunched)
				m_free_slots.push_back(pi_slots[j]);
		}
	}
	//========================================================================
	// Creation of singularity lines from cross field singular points
	//========================================================================

	unsigned int cont = 0;
	while(!m_free_slots.empty()){
		SingularityPoint::Slot* current_slot = m_free_slots.front();
		m_free_slots.pop_front();
		if(withGlobalComments){
			cout<<"current_slot->location "<<current_slot->location<<endl;
			if(current_slot->line)
				cout<<"begining current_slot->line has line "<<current_slot->line->getNumber()<<endl;
		}
		if(current_slot->isLaunched!=true) {
			current_slot->isLaunched = true;
			if(withGlobalComments)
				cout<<"computeSingularityLine"<<endl;
			SingularityPoint::Slot* removed_slot = 0;
			computeSingularityLine(current_slot->from_point, current_slot, cont, removed_slot);

			if(removed_slot!=0){
				m_free_slots.push_back(removed_slot);
				removed_slot->isLaunched = false;
				if(withGlobalComments){
					cout<<"connection to slot has been removed; removed_slot->location "<<removed_slot->location<<endl;
					if(removed_slot->line)
						cout<<"connection to slot has been removed removed_slot->line has line "<<removed_slot->line->getNumber()<<endl;
				}
			}
			cont++;
			writeOutput("boundary_line");
		}
		if(withGlobalComments){
			if(current_slot->line)
				cout<<"end current_slot->line has line "<<current_slot->line->getNumber()<<endl;

			cout<<"each step m_free_slots detectLineIntersections sing points "<<endl;
			vector<SingularityPoint*> pi = m_graph.getPoints();
			for(unsigned int i=0;i<pi.size(); i++)
				cout<<pi[i]->getLocation().X()<<" "<<pi[i]->getLocation().Y()<<endl;
		}
	}

}

/*----------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::computeSingularityLine(SingularityPoint*        AFromPoint,
                                                       SingularityPoint::Slot*  AFromSlot,
                                                       unsigned int&            cont,
                                                       SingularityPoint::Slot*& ARemovedSlot)
{
	//========================================================================
	// Data initialization for line building
	//========================================================================
	bool end_on_bnd = false;
	bool end_on_free_slot = false;
	bool must_create_pnt = false;
	SingularityPoint       *to_sing_pnt  = 0;
	SingularityPoint::Slot *to_slot = 0;
	TCellID to_cell_id;
	int to_cell_dim;
	math::Point  to_pnt;
	math::Vector3d to_dir;
	std::vector<math::Point> line_discretization;
	std::vector<TCellID>     line_triangles;
	double                   streamlineDeviation = 0.0;
	//========================================================================
	// Stream line computation
	//========================================================================
	computeStreamLine(AFromPoint,
	                  AFromSlot,
	                  to_sing_pnt,
	                  to_slot,
	                  to_pnt,
	                  to_dir,
	                  line_discretization,
	                  line_triangles,
	                  to_cell_dim,
	                  to_cell_id,
	                  streamlineDeviation,
	                  end_on_bnd,
	                  end_on_free_slot,
	                  must_create_pnt);
	/*gmds::Mesh meshSing(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));

	for(unsigned int t=0; t<line_triangles.size(); t++) {
		vector<gmds::Node> currentTriNodes = m_mesh->get<Face>(line_triangles[t]).get<Node>();
		cout<<"currentTriNodes,isze() "<<currentTriNodes.size()<<endl;
		gmds::Node mySing1 = meshSing.newNode(currentTriNodes[0].getPoint());
		gmds::Node mySing2 = meshSing.newNode(currentTriNodes[1].getPoint());
		gmds::Node mySing3 = meshSing.newNode(currentTriNodes[2].getPoint());
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
	if(withGlobalComments){
		cout<<"line_discretization[0] "<<line_discretization[0].X()<<" "<<line_discretization[0].Y()<<endl;
		cout<<"line_discretization[1] "<<line_discretization[1].X()<<" "<<line_discretization[1].Y()<<endl;
	}
	//line creation
	SurfaceSingularityLine* surf_line = m_graph.newSurfaceLine();
	int sepNumberTmp = m_graph.getNbLines();
	surf_line->setNumber(sepNumberTmp);

	//connect line to initial singularity point
	SingularityPoint* from_sing_pnt = AFromSlot->from_point;
	surf_line->addSingularityPoint(from_sing_pnt);
	surf_line->addDiscretizationPoint(from_sing_pnt->getLocation());
	if(withGlobalComments)
		cout<<"surf_line also adds disc pt "<<from_sing_pnt->getLocation().X()<<from_sing_pnt->getLocation().Y()<<endl;

	AFromSlot->line = surf_line;
	//AFromSlot->line_direction =  AFromSlot->direction;
	math::Vector3d firstDir = math::Vector3d(from_sing_pnt->getLocation(), line_discretization[0]);
	if(withGlobalComments)
		cout<<"firstDir "<<firstDir.X()<<" "<<firstDir.Y()<<endl;
	AFromSlot->line_direction = firstDir;
	AFromSlot->isLaunched = true;
	if(to_slot!=0){
		streamlineDeviation = (streamlineDeviation + fabs(1 - to_slot->direction.dot(math::Vector3d(line_discretization[line_discretization.size()-1],line_discretization[line_discretization.size()-2]))));
		streamlineDeviation = streamlineDeviation/(line_discretization.size()+1);
		AFromSlot->lineDeviation = streamlineDeviation;
	}
	//Insertion of line points
	for(unsigned int i=0; i<line_discretization.size(); i++) {
		surf_line->addDiscretizationPoint(line_discretization[i]);
	}

	for(unsigned int i=0; i<line_triangles.size(); i++) {
		surf_line->addTraversedFace(line_triangles[i]);
	}

	//========================================================================
	// Termination of the streamline
	//========================================================================
	if(end_on_bnd) {
		if(withGlobalComments)
			cout<<"!!!!!!!!!!!!!!!!!end_on_bnd "<<endl;
		//======================================================================
		// CASE 1 - We finish on the boundary. A geometric point must be created
		// The cell defined by (start_cell_dim, start_cell_id) is on the boundary.
		// We have to create a geometric singularity point so.
		// Here it always assumes we reach a non nodeOnPoint vertex; otherwise just error

		SingularityPoint* geom_pnt;
		SingularityPoint::Slot* incoming_slot;
		createGeometricSingularityPoint(to_pnt,      // the point we are
		                                to_dir,      // the direction we come from
		                                to_cell_dim, // the dim. of the cell
		                                to_cell_id,  // the id of the cell
		                                geom_pnt,       // the created point
		                                incoming_slot); // and the slot

		AFromSlot->lineDeviation = streamlineDeviation/line_discretization.size();
		if(withGlobalComments)
			cout<<"AFromSlot->lineDeviation "<<AFromSlot->lineDeviation<<endl;

		surf_line->addSingularityPoint(geom_pnt);
		// surf_line->addDiscretizationPoint(geom_pnt->getLocation());/* This has already been added, but check what happens if end on node on point*/
		if(withGlobalComments)
			cout<<"from "<<AFromPoint->getMesh<Face>()[0]<<" towards bdry point "<<to_pnt.X()<<" "<<to_pnt.Y()<<endl;
		//    geom_pnt->connectLine(surf_line, to_dir);
	} // if(to_cell_dim==0)
	else {
		if(withGlobalComments)
			cout<<"!!!!!!!!!!!!!!!!! not end_on_bnd "<<endl;

		if(end_on_free_slot){
			//======================================================================
			// CASE 2 - We finish on a field singularity where a slot is free
			// I don't need to backtrack, i already have the line
			// TODO connect it here; the same situation for an occupied slot for which the new line is better aligned
			if(withGlobalComments)
				cout<<"end_on_free_slot, to_slot->isLaunched "<<to_slot->isLaunched<<endl;
		}
		// else{
		//============================================
		if(!end_on_free_slot){
			//======================================================================
			// CASE 3 - We finish on a field singularity where a slot is not free
			//The slot must be freed before connection.
			if(withGlobalComments)
				cout<<"!end_on_free_slot"<<endl;
			SingularityPoint* singPointToRemove = 0;
			if(to_slot!=0){
				SingularityLine* removedLine = to_slot->line;
				SingularityPoint* otherSingPnt;

				vector<SingularityPoint*> removedLineSingPoints = removedLine->getEndPoints();
				if(removedLineSingPoints[0]==to_sing_pnt){
					otherSingPnt = removedLineSingPoints[1];
				}
				else
					otherSingPnt = removedLineSingPoints[0];

				if(otherSingPnt->getType()==0){
					std::vector<SingularityPoint::Slot*> otherSingPnt_slots = otherSingPnt->getSlots();

					for(unsigned int t=0; t<otherSingPnt_slots.size(); t++){
						if(otherSingPnt_slots[t]->line->getNumber() == to_slot->line->getNumber())
							ARemovedSlot = otherSingPnt_slots[t];
					}
					if(withGlobalComments){
						cout<<"to_slot!=0; to_slot->line "<<to_slot->line->getNumber()<<" will be removed "<<endl;
						cout<<"this singline is between ["<<to_slot->line->getEndPoints()[0]->getLocation().X()<<","<<to_slot->line->getEndPoints()[0]->getLocation().Y()<<"]"<<endl;
						cout<<"and ["<<to_slot->line->getEndPoints()[1]->getLocation().X()<<","<<to_slot->line->getEndPoints()[1]->getLocation().Y()<<"]"<<endl;
						cout<<"to_slot->starting_cell_id "<<to_slot->starting_cell_id<<endl;
					}
				}
				else{ /* the other sing point is a geometric one;
					if it has been created as streamline-bdry intersection => it should be removed */
					vector<SingularityPoint::Slot*> curvePointSlots = otherSingPnt->getSlots();
					for(unsigned int t=0; t<curvePointSlots.size(); t++){
						if(!curvePointSlots[t]->isFreeze){
							singPointToRemove = otherSingPnt;
							break;
						}
					}
				}
			}
			//vector<SingularityLine*> linesAdded = m_graph.getCurveLines();
			removeSingularityLine(to_slot->line);
			m_graph.removeLine(to_slot->line);

			if(singPointToRemove!=0){
				m_graph.removePoint(singPointToRemove);
			}
		}

		to_slot->isLaunched = true;
		to_slot->line = surf_line;
		to_slot->line_direction = to_dir;
		to_slot->lineDeviation = streamlineDeviation;

		/* We need to compute the discretization points, which are in the
		triangles of the confusing ball.*/
		surf_line->addSingularityPoint(to_sing_pnt);
		surf_line->addDiscretizationPoint(to_sing_pnt->getLocation());

		/* vector<math::Point> Newline_discretization = surf_line->getDiscretizationPoints();
		gmds::Mesh meshSing(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
		cout<<"first elem Newline_discretization "<<endl;
	 	for(unsigned int j=0; j<Newline_discretization.size(); j++) {
			gmds::Node mySing = meshSing.newNode(Newline_discretization[j].X(), Newline_discretization[j].Y(), Newline_discretization[j].Z());
		    	meshSing.newTriangle(mySing, mySing, mySing);
		}
		gmds::IGMeshIOService ioServiceSing(&meshSing);
		gmds::VTKWriter vtkWriterSing(&ioServiceSing);
		vtkWriterSing.setCellOptions(gmds::N|gmds::F);
		vtkWriterSing.setDataOptions(gmds::N|gmds::F);
		std::string file_name = "HolesInSquare-discretizationPointsBeforeBacktrack_"+std::to_string(m_graph.getNbLines())+".vtk";
		vtkWriterSing.write(file_name);  */

		streamlineDeviation = 0.0;
		bool foundPath=false;
		if(withGlobalComments)
			cout<<"backtrackSingularityLine from "<<to_sing_pnt->getMesh<Face>()[0]<<", backtrackSingularityLine towards "<<AFromPoint->getMesh<Face>()[0]<<endl;
		backtrackSingularityLine(surf_line, // the line we modify
		                         to_sing_pnt, // the point we start from
		                         to_slot,// the slot we start from
		                         AFromPoint,// the point we go to
		                         AFromSlot,// the slot we go to
		                         to_dir, // the direction of the streamline for to_slot
		                         foundPath);  // boolean value indicating if we have found a path

		if(foundPath){
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
			for(unsigned int i=1; i<finalLineDiscretization.size()-1; i++){
				previous_point = current_point;
				current_point = finalLineDiscretization[i];
				previous_dir = math::Vector3d(previous_point, current_point);
				next_dir = math::Vector3d(current_point, finalLineDiscretization[i+1]);
				if(i==1){
					streamlineDeviation = streamlineDeviation + fabs(1.0 - previous_dir.dot(AFromSlot->direction));
				}
				if(i==finalLineDiscretization.size()-2){
					streamlineDeviation = streamlineDeviation + fabs(1.0 - next_dir.dot(to_slot->direction.opp()));
				}
				foundTri = false;
				while(contTri<finalLineTraversedTriangles.size() && !foundTri){
					if(m_tool.isPntInTri(current_point, m_mesh->get<Face>(finalLineTraversedTriangles[contTri]), temp, temp, temp, lambda0, lambda1)){
						vector<gmds::Node> current_verts = m_mesh->get<Face>(finalLineTraversedTriangles[contTri]).get<Node>();
						math::Cross2D cross_0 = (*m_field)[current_verts[0].id()];
						math::Cross2D cross_1  = (*m_field)[current_verts[1].id()];
						math::Cross2D cross_2  = (*m_field)[current_verts[2].id()];
						closest0 =  cross_0.closestComponentVector(previous_dir);
						closest1 =  cross_1.closestComponentVector(previous_dir);
						closest2 =  cross_2.closestComponentVector(previous_dir);
						closest2Current = lambda0 * closest0 + lambda1 * closest1 + (1 - lambda0 - lambda1) * closest2;
						streamlineDeviation = streamlineDeviation + fabs(1.0 - previous_dir.dot(closest2Current));

						closest0 =  cross_0.closestComponentVector(next_dir);
						closest1 =  cross_1.closestComponentVector(next_dir);
						closest2 =  cross_2.closestComponentVector(next_dir);
						closest2Current = lambda0 * closest0 + lambda1 * closest1 + (1 - lambda0 - lambda1) * closest2;
						streamlineDeviation = streamlineDeviation + fabs(1.0 - next_dir.dot(closest2Current));
					}
					contTri++;
				}
			}
			streamlineDeviation = streamlineDeviation/finalLineDiscretization.size();

			to_slot->lineDeviation = streamlineDeviation;
			AFromSlot->lineDeviation = streamlineDeviation;
		}
	}
}

/*----------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::computeStreamLine(SingularityPoint*          AFromPnt,
                                                  SingularityPoint::Slot*    AFromSlot,
                                                  SingularityPoint*&         AToSingPnt,
                                                  SingularityPoint::Slot*&   AToSlot,
                                                  math::Point&               AToPnt,
                                                  math::Vector3d&            AToDir,
                                                  std::vector<math::Point>&  APoints,
                                                  std::vector<TCellID>&      ATriangles,
                                                  int&                       AToCellDim,
                                                  TCellID&                   AToCellID,
                                                  double&                    streamlineDeviation,
                                                  bool&                      AEndOnBnd,
                                                  bool&                      AToSlotIsFree,
                                                  bool&                      APntToCreate)
{
	if(withGlobalComments)
		cout<<"begining computeStreamLine from point "<<AFromPnt->getLocation().X()<<" , "<<AFromPnt->getLocation().Y()<<endl;

	ATriangles.clear();
	APoints.clear();

	math::Point  start_pnt = AFromSlot->location ; //starting point
	APoints.push_back(start_pnt);
	math::Vector3d start_dir = AFromSlot->direction; //starting direction
	math::Vector3d prev_dir  = AFromSlot->direction; //prev direction used in the
	/*termination of extrapolation process (when we get into a
	confusing ball) */

	TCellID start_cell_id  = AFromSlot->starting_cell_id ;
	int start_cell_dim = AFromSlot->starting_cell_dim;

	math::Point  current_pnt = start_pnt;
	math::Vector3d current_vec = start_dir;

	if(start_cell_dim==0)
		cout<<"node: "<<start_cell_id<<endl;
	else
		if(start_cell_dim==1){
			vector<Node> currentNodes = (m_mesh->get<Edge>(start_cell_id)).get<Node>();
			if(withGlobalComments)
				cout<<"edge: "<<start_cell_id<<" between "<<currentNodes[0].id()<<" "<<currentNodes[1].id()<<endl;
		}

	math::Point start_dirPnt(start_pnt.X() + start_dir.X(),
	                         start_pnt.Y() + start_dir.Y(),
	                         start_pnt.Z() + start_dir.Z());

	//========================================================================
	// Singularity point we will be connecting to. It can be another
	// singularity point of the cross field or a geometric singularity point
	// created on the fly.
	//========================================================================
	SingularityPoint*       found_pnt = 0;
	SingularityPoint::Slot* found_slot =0;

	bool find_end = false;
	/* indicates that we reach a boundary point or line */
	bool end_on_boundary = false;
	/* indicates that we reach an existing singularity point*/
	//bool end_on_field_singularity = false;

	// We check some termination conditions on the boundary.
	if (start_cell_dim==0){
		Node currentNode = m_mesh->get<Node>(start_cell_id);
		if(m_mesh->isMarked(currentNode, m_mark_nodes_on_point) ||
		m_mesh->isMarked(currentNode, m_mark_nodes_on_curve)){
			find_end        = true;
			end_on_boundary = true;
		}
	}
	else { //we have necessarry start_cell_dim=1
		Edge currentEdge = m_mesh->get<Edge>(start_cell_id);
		if(m_mesh->isMarked(currentEdge, m_mark_edges_on_curve)){
			find_end        = true;
			end_on_boundary = true;
		}
	}
	//========================================================================
	// Main loop to create the singularity line
	//========================================================================
	// int nbwalk=0;
	while(!find_end) {
		//nbwalk++;
		TCellID next_cell_id  = NullID;
		int next_cell_dim = -1;
		m_tool.findNextCell(start_pnt,
		                    start_dir,
		                    start_cell_dim,
		                    start_cell_id,
		                    next_cell_dim,
		                    next_cell_id);

		if (next_cell_dim == -1){
			if(withGlobalComments)
				cout<<"next_cell_dim == -1"<<endl;
			// The cell defined by (start_cell_dim, start_cell_id) is on the boundary.
			find_end = true;
			end_on_boundary = true;
		}
		else if (next_cell_dim == 1){
			if(withGlobalComments)
				cout<<"next_cell_dim == 1"<<endl;
			/* we are going along an edge.
			Our simple assumption is to follow this edge until reaching
			one of its end points and to compute the next direction at
			this point.*/
			// WARNING, we do not check the confusing ball!!!!!
			Edge currentEdge = m_mesh->get<Edge>(next_cell_id);
			std::vector<TCellID> adj_faces = currentEdge.getIDs<Face>();
			ATriangles.insert(ATriangles.end(),adj_faces.begin(),adj_faces.end());

			std::vector<Node> currentNodes = currentEdge.get<Node>();
			math::Vector3d v0(start_pnt, currentNodes[0].getPoint());
			math::Vector3d v1(start_pnt, currentNodes[1].getPoint());
			Node next_node;
			if(math::near(v0.norm(),0.0))
				next_node = currentNodes[1];
			else if(math::near(v1.norm(),0.0))
				next_node = currentNodes[0];
			else if(v0.dot(start_dir) > v1.dot(start_dir))
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
			find_end = false;

			APoints.push_back(start_pnt);
		}
		else { //general case, we are in a face
			Face currentFace = m_mesh->get<Face>(next_cell_id);
			//cout<<"next_cell_id "<<next_cell_id<<endl;
			ATriangles.push_back(currentFace.id());
			//==============================================================
			// CASE 1: DO WE ARE IN A FACE CONTAINING A SING. POINT?
			//==============================================================
			bool intersect_sing_point = false;

			SingularityPoint* next_sing_point = m_faces_to_singularity_on_surf[currentFace.id()];
			bool must_try_to_connect = false;
			if (next_sing_point != NULL && //face in a confusing ball ...
			next_sing_point != AFromPnt) {// ... of another singularity point
				if(withGlobalComments)
					cout<<"face in a confusing ball of another singularity point"<<endl;
				must_try_to_connect = true;
			}
			else if (next_sing_point != NULL && //face in the confusing ball ...
			next_sing_point == AFromPnt) {//... of the incoming singularity point
				if(withGlobalComments)
					cout<<"face in a confusing ball of the incoming singularity point"<<endl;
				//Warning: completly empiric, we just try to detect cyclic lines
				if (APoints.size() >= 100)
					must_try_to_connect = true;
			}

			if(must_try_to_connect) {
				if(withGlobalComments)
					cout<<"must_try_to_connect"<<endl;
				math::Point start_dirPnt(start_pnt.X() + start_dir.X(),
				                         start_pnt.Y() + start_dir.Y(),
				                         start_pnt.Z() + start_dir.Z());

				//Now, we look for a compatible slot
				std::vector<SingularityPoint::Slot*>& cur_slots = next_sing_point->getSlots();
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
						SingularityPoint::Slot* current_slot = cur_slots[i_slot];
						if(current_slot->isFreeze)
							continue;
						math::Vector3d slot_opp_dir = current_slot->direction.opp();
						slot_opp_dir.normalize();
						double slot_deviation = slot_opp_dir.dot(current_vec);
						if (slot_deviation > slot_epsilon && slot_deviation > best_deviation) {
							best_deviation = slot_deviation;
							best_slot_id = i_slot;
						} //if (slot_deviation < slot_epsilon) {
					} //for (unsigned int i_slot = 0; !found_free_slot && i_slot < ....

					if(best_deviation!=-2 &&  !cur_slots.empty()){
						if(withGlobalComments)
							cout<<"best_deviation!=-2 &&  !cur_slots.empty()"<<endl;
						SingularityPoint::Slot* best_slot = cur_slots[best_slot_id];
						math::Vector3d            slot_opp_dir = best_slot->direction.opp();
						if(withGlobalComments){
							cout<<" best_slot->direction "<<best_slot->direction.X()<<" "<<best_slot->direction.Y()<<endl;
							cout<<"current_vec "<<current_vec.X()<<" "<<current_vec.Y()<<endl;
						}
						slot_opp_dir.normalize();
						if ( best_slot->isLaunched) {
							//slot already assigned with a previous line (and direction)
							if(withGlobalComments)
								cout<<"slot already assigned with a previous line (and direction)"<<endl;

							math::Vector3d prev_line_dir = best_slot->line_direction;
							if(withGlobalComments)
								cout<<" best_slot->line_direction "<<best_slot->line_direction.X()<<" "<<best_slot->line_direction.Y()<<endl;
							prev_line_dir.normalize();
							double prev_deviation =slot_opp_dir.dot(prev_line_dir);
							if(withGlobalComments)
								cout<<"prev_deviation "<<prev_deviation<<endl;
							if(best_deviation>prev_deviation) {
								//the new alignment is better than the previous one
								if(withGlobalComments)
									cout<<"best_deviation "<<best_deviation<<" > "<<"prev_deviation "<<prev_deviation<<", the new alignment is better than the previous one"<<endl;

								found_free_slot = false;
								found_compatible_slot = true;
								found_slot = best_slot;
								intersect_sing_point = true;
							}
							else{
							 // We keep the previous association
								if(withGlobalComments)
									cout<<"We keep the previous association"<<endl;

								found_compatible_slot = false;
							}
						} //if ( current_slot>isLaunched) {
						else{  //WE HAVE A FREE SLOT
						 // We keep the previous association
							if(withGlobalComments)
								cout<<"WE HAVE A FREE SLOT"<<endl;

							found_free_slot = true;
							AToSlotIsFree = found_free_slot;
							found_compatible_slot = true;
							found_slot = best_slot;
							intersect_sing_point = true;
						}
						// HAVE WE FOUND THE END OF THE LINE??
						if(found_compatible_slot) {
							if(withGlobalComments)
								cout<<"last    found_compatible_slot"<<endl;
							// COMPATIBLE AND FREE SLOT
							find_end = true;
						}
					}
					slot_epsilon -=0.1;
				}//while (!found_compatible_slot && slot_epsilon < -0.4)
			}//	if(must_try_to_connect) {

			//==============================================================
			// CASE 2: SO, WE JUST HAVE TO CROSS THE TRIANGLE (GENERAL CONF)
			//==============================================================
			//Does the current triangle has the same classif
			if ( !intersect_sing_point) {
				math::Point  out_pnt;
				math::Vector3d out_vec;
				TCellID out_cell_id;
				int out_cell_dim;
				if(withGlobalComments){
					cout<<"start_pnt "<<start_pnt.X()<<" "<<start_pnt.Y()<<endl;
					cout<<"currentFace "<<currentFace.id()<<endl;
				}
				m_tool.traverseTriangle(currentFace,   /* the face we work on*/
				                        start_pnt,      /* the geometric point we start from */
				                        start_dir,      /* the geometric direction to follow*/
				                        start_cell_dim, /* the dimension of the cell start_pnt is located */
				                        start_cell_id,  /* the id of the cell start_pnt is located on*/
				                        out_pnt,        /* the geometric point where we go out */
				                        out_vec,        /* the geometric direction to follow after*/
				                        out_cell_dim,   /* the dimension of the out cell (0 or 1) */
				                        out_cell_id,    /* the id of the out cell*/
				                        streamlineDeviation);   /*deviation of the streamline up to this point*/
				if(withGlobalComments)
					cout<<"after traverseTriangle out_pnt= "<<out_pnt.X()<<" "<<out_pnt.Y()<<endl;
				APoints.push_back(out_pnt);

				// we progress to the next point, next vector and so next face too
				prev_dir = start_dir; //we store the prev direction for slot
				//reconnection with balls
				start_pnt = out_pnt;
				start_dir = out_vec;
				start_cell_dim = out_cell_dim;
				start_cell_id = out_cell_id;
			} //if (!intersect_line && !intersect_sing_point) {

			if (intersect_sing_point){
				if(withGlobalComments)
					cout<<"intersect_sing_point"<<endl;
				find_end = true;
			}

			//post process, we just check whether we have arrived onto a geometric boundary
			if (start_cell_dim==0){
				if(withGlobalComments)
					cout<<"post process start_cell_dim==0"<<endl;
				Node currentNode = m_mesh->get<Node>(start_cell_id);
				if(m_mesh->isMarked(currentNode, m_mark_nodes_on_point) ||
				m_mesh->isMarked(currentNode, m_mark_nodes_on_curve)){
					find_end        = true;
					end_on_boundary = true;
				}
			}
			else { //we have necessarry start_cell_dim=1
				if(withGlobalComments)
					cout<<"post process we have necessarry start_cell_dim=1"<<endl;
				Edge currentEdge = m_mesh->get<Edge>(start_cell_id);
				if(m_mesh->isMarked(currentEdge, m_mark_edges_on_curve)){
					find_end        = true;
					end_on_boundary = true;
				}
			}
		} // else { //general case, we are in a face
	} //while(!find_end)

	//==============================================================
	// Update of out parameters
	//==============================================================
	//last followed direction
	AToDir = start_dir;
	AToPnt = start_pnt;
	AEndOnBnd = end_on_boundary;
	//the end point must be created if it has not been found
	APntToCreate = (found_pnt==0);
	//singularity point data if we found an end point
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

/*----------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::createGeometricSingularityPoint(const math::Point&       AInPnt,
                                                                const math::Vector3d&    AInVec,
                                                                const int                ACellDim,
                                                                const TCellID            ACellID,
                                                                SingularityPoint*&       APnt,
                                                                SingularityPoint::Slot*& ASlot)
{
	if(withGlobalComments){
		cout<<"createGeometricSingularityPoint"<<endl;
		cout<<"ACellDim "<<ACellDim<<endl;
	}
	if(ACellDim==0){
		if(withGlobalComments)
			cout<<"ACellDim==0"<<endl;
		Node n = m_mesh->get<Node>(ACellID);

		if(m_mesh->isMarked(n,m_mark_nodes_on_point)){
			if(withGlobalComments)
				cout<<"We arrive onto a geometric point"<<endl;
			/*We arrive onto a geometric point !!!! It means that a singularity geometric point already exist for it
			We look for it*/
			std::vector<VertexSingularityPoint* >
			geom_sing_points = m_graph.getVertexPoints();

			bool found_pnt = false;
			for(unsigned int i=0; i<geom_sing_points.size() && !found_pnt; i++){
				VertexSingularityPoint* current_pnt = geom_sing_points[i];
				Node currentNode = current_pnt->getMeshNode();
				if(currentNode.id()==ACellID){
					found_pnt = true;
				}
			}
			if(!found_pnt){
				throw GMDSException("createGeometricSingularityPoint: Error, no geometric point to be connected to!");
			}
		} // if(m_mesh->isMarked(n,m_mark_nodes_on_point))
		else{     	// We are on a geometrical curve
			m_graph.splitCurveLine(AInPnt, AInVec, n, APnt, ASlot);
			cout<<"has split curve"<<endl;
		}
 	} // if(ACellDim!=0)
	else {
		Edge e = m_mesh->get<Edge>(ACellID);
		vector<gmds::Node> myTempNodes = e.get<gmds::Node>();
		if(withGlobalComments)
			cout<<"ACellDim!=0, edge "<<e.id()<<" in between nodes"<<myTempNodes[0]<<" - "<<myTempNodes[1]<<endl;
		m_graph.splitCurveLine(AInPnt, AInVec, e, APnt, ASlot);
	}
}

/*----------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::backtrackSingularityLine(SurfaceSingularityLine*    ALine,
                                                         SingularityPoint*          AFromPnt,
                                                         SingularityPoint::Slot*    AFromSlot,
                                                         SingularityPoint*          AToPnt,
                                                         SingularityPoint::Slot*    AToSlot,
                                                         gmds::math::Vector3d&      to_dir,
                                                         bool&                      foundBackTrackPath)
{
	//==============================================================
	// We keep in mind the length of the line
	//==============================================================
	if(withGlobalComments)
		cout<<"backtrackSingularityLine"<<endl;
		//double line_length = ALine->length();

	SingularityPoint*        arrival_sing;
	SingularityPoint::Slot*  arrival_slot;
	math::Point              arrival_pnt;
	math::Vector3d           arrival_dir;
	int                      arrival_cell_dim;
	TCellID                  arrival_cell_id;
	std::vector<math::Point> line_pnts;
	std::vector<TCellID>     new_traversed_triangles;
	bool arrival_on_bnd;
	bool arrival_on_free_slot;
	bool arrival_pnt_to_create;

	vector<gmds::TCellID> modifiedFaces;
	double previousRad = m_confusing_distance;
	double redefRadius=1.5;
	int cont = 0; // here take care not to bump into a new confusing ball
	double streamlineDeviation = 0.0;

	while((arrival_sing!=AToPnt)&&(cont<10)){
		computeStreamLine(AFromPnt,
		                  AFromSlot,
		                  arrival_sing,
		                  arrival_slot,
		                  arrival_pnt,
		                  arrival_dir,
		                  line_pnts,
		                  new_traversed_triangles,
		                  arrival_cell_dim,
		                  arrival_cell_id,
		                  streamlineDeviation,
		                  arrival_on_bnd,
		                  arrival_on_free_slot,
		                  arrival_pnt_to_create);

		if(arrival_sing!=AToPnt){
			cont++;
			redefineOneConfusingBall(AToPnt, modifiedFaces, previousRad, redefRadius);
			redefRadius=redefRadius+0.5;
		}
	}

	if(arrival_sing!=AToPnt){
		if(withGlobalComments){
			cout<<"we have started from triangle "<<AFromPnt->getMesh<Face>()[0].id()<<" and we want to reach "<<AToPnt->getMesh<Face>()[0].id()<<endl;
			cout<<"backtrack didn't reach destination point"<<endl;
		}
            //throw GMDSException("Backtraking issue: We don't reach the departure point");
	}
	else{
		foundBackTrackPath = true;
		if(withGlobalComments)
			cout<<"backtrack we have reached destination point"<<endl;

		if(modifiedFaces.size()!=0){
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

	if(cont>0){
		for(unsigned int i=0; i<modifiedFaces.size();i++){
			m_faces_to_singularity_on_surf[modifiedFaces[i]] = 0;
		}
	}

	if(foundBackTrackPath){
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
		math::Point last_point = old_pnts[old_pnts.size()-1];

		for(unsigned int i_old = 1; i_old<old_pnts.size()-1;i_old++){
			math::Point current_pnt = old_pnts[i_old];
			math::Segment seg(line_pnts[0],line_pnts[1]);
			math::Point proj_pnt = seg.project(current_pnt);
			double proj_dist = current_pnt.distance(proj_pnt);

			for(unsigned int i_new = 2; i_new<line_pnts.size();i_new++){
				math::Segment seg_i(line_pnts[i_new-1],line_pnts[i_new]);
				math::Point proj_pnt_i = seg_i.project(current_pnt);
				double proj_dist_i = current_pnt.distance(proj_pnt_i);
				if(proj_dist_i<proj_dist){
					proj_dist = proj_dist_i;
					proj_pnt = proj_pnt_i;
				}
			}
			// We have our new point
			new_pnts.push_back(0.5*current_pnt+0.5*proj_pnt);
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
		for(unsigned int i=0; i<all_traversed_faces.size();i++) {
			Face f=m_mesh->get<Face>(all_traversed_faces[i]);
			std::vector<Node> f_nodes = f.get<Node>();
			for(unsigned int j=0;j<f_nodes.size();j++){
				Node nj = f_nodes[j];
				std::vector<TCellID> nj_faces = nj.getIDs<Face>();
				for(unsigned int k=0;k<nj_faces.size();k++){
					candidates.insert(nj_faces[k]);
				}
			}
		}

		std::set<TCellID> set_of_traversed_faces;

		for(std::set<TCellID>::iterator it = candidates.begin(); it!=candidates.end(); it++) {
			Face f = m_mesh->get<Face>(*it);
			std::vector<Node> f_nodes = f.get<Node>();
			math::Triangle t(f_nodes[0].getPoint(), f_nodes[1].getPoint(), f_nodes[2].getPoint());
			bool found_pnt = false;
			for(unsigned int j=0;j<new_pnts.size() && !found_pnt; j++){
				math::Point pj = new_pnts[j];
				if(t.isIn(pj)){
					found_pnt = true;
					set_of_traversed_faces.insert(f.id());
				}
				if(j!=0){
					math::Point pk = new_pnts[j-1];
					if(t.intersect(math::Segment(pj,pk))){
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

/*----------------------------------------------------------------------------*/

std::vector<SurfaceSingularityLine*> SingularityGraphBuilder2D::getSingularityLinesIn(const Face& AFace)
{
 /*get the singularity lines that traverse triangle AFace*/
	std::vector<SurfaceSingularityLine* > found_lines;

	std::vector<SurfaceSingularityLine* > lines = m_graph.getSurfaceLines();

	TCellID face_id = AFace.id();
	for(unsigned int i=0; i<lines.size(); i++){
		SurfaceSingularityLine* line_i = lines[i];
		if(line_i->isTraversed(face_id)){
			found_lines.push_back(line_i);
		}
	}

	return found_lines;
}

/*----------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::writeOutput(const std::string& AFileName)
{
	static int out = 0;
	std::string file_name = AFileName+"_"+std::to_string(out);

	writeOutputSingle(file_name);
	out++;
}

/*----------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::writeOutputSingle(const std::string& AFileName)
{
	std::string file_name = m_output_directory_name +"-"+AFileName+".vtk";

	m_graph.createVTKOutputFile(file_name);

}

/*----------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::writeOutputPatches(const std::string& AFileName)
{
	std::string file_name = AFileName+".vtk"; //= m_output_directory_name +"-"+AFileName+".vtk";
	bool curvePatches = true;
	m_graph.createVTKOutputFile(file_name, curvePatches);

}

/*-----------------------------------------------------------------*/

void SingularityGraphBuilder2D::writeSingularityPointsAndSlots()
{
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
	std::vector<SingularityPoint*> points = m_graph.getPoints();
	std::vector<SingularityPoint*>::iterator it = points.begin();
	for(;it!=points.end();it++){
		SingularityPoint* sing = *it;
		math::Point loc = sing->getLocation();
		gmds::Node center = m.newNode(loc.X(), loc.Y(), loc.Z());
		std::vector< SingularityPoint::Slot*>& slots = sing->getSlots();

		for (unsigned int i = 0; i <slots.size(); i++){

			SingularityPoint::Slot* si = slots[i];

			math::Point si_loc = si->location;
			math::Vector3d si_dir = si->direction;
			si_dir.normalize();
			double kappa = 10;
			gmds::Node si_departure = m.newNode(si_loc.X(), si_loc.Y(), si_loc.Z());
			gmds::Node si_end = m.newNode(si_loc.X() + si_dir.X()/kappa,
			                              si_loc.Y() + si_dir.Y()/kappa,
			                              si_loc.Z() + si_dir.Z()/kappa );

			m.newTriangle(center, si_departure, si_departure);
			m.newTriangle(si_end, si_departure, si_departure);
		}
	}

	gmds::IGMeshIOService meshIoServ(&m);
	gmds::VTKWriter writer(&meshIoServ);
	writer.setCellOptions(gmds::N|gmds::F);
	writer.setDataOptions(gmds::N|gmds::F);

	std::stringstream file_name;
	file_name <<m_output_directory_name<<"-points_and_slots.vtk";
	writer.write(file_name.str());

}

/*----------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::redefineOneConfusingBall(SingularityPoint*      APnt,
                                                         vector<gmds::TCellID>& modifiedFaces,
                                                         double                &previousDist,
                                                         const double           increaseRadiusScale)
{
	if(withGlobalComments)
		cout<<"redefineOneConfusingBall "<<endl;
	// here it could be a problem if two sing are very close geometrically one to another

	unsigned int nbVerts = m_mesh->getNbNodes();

	vector<bool> visitedFaces(original_faces_number, false);
	vector<bool> visitedVerts(nbVerts, false);

	math::Point sing_location = APnt->getLocation();
	double newConfusingDistance = increaseRadiusScale*m_confusing_distance;

	std::vector<gmds::Face> neighbouringFaces;
	std::vector<gmds::Node> vertsToVisit;
	//  cout<<"APnt0>getNbMeshCells() "<<APnt->getNbMeshCells()<<endl;

	vector<gmds::Face>  AtriVectFace = APnt->getMesh<Face>();

	if(AtriVectFace.size()!=0){
		Face Atri = AtriVectFace[0];
		if(withGlobalComments)
			cout<<"face "<<Atri.id()<<endl;
		std::vector<Node> nodes = Atri.get<Node>();
		for(unsigned int i=0; i<nodes.size(); i++)
			vertsToVisit.push_back(nodes[i]);
	}
	else{//geom point
		vector<gmds::Node> AtriVectNode = APnt->getMesh<Node>();

		if(AtriVectNode.size()!=0){
			gmds::Node Atri = AtriVectNode[0];
			vertsToVisit.push_back(Atri);
			if(withGlobalComments)
				cout<<"node "<<Atri.id()<<endl;
		}
		else{
			vector<gmds::Edge> AtriVectEdge = APnt->getMesh<Edge>();

			if(AtriVectEdge.size()!=0){
				gmds::Edge Atri = AtriVectEdge[0];
				std::vector<Node> nodes = Atri.get<Node>();
				vertsToVisit.push_back(nodes[0]);
				vertsToVisit.push_back(nodes[1]);
				if(withGlobalComments)
					cout<<"edge with nodes "<<nodes[0].id()<<" "<<nodes[1]<<endl;
			}
			//else curve sing point -> doesnt need confusing ball for now
		}
	}

  	/*
	for(unsigned int i=0; i<originalConfusingBalls.size(); i++){
     	visitedFaces[originalConfusingBalls[i]] = true;
     	gmds::Face currentFace = m_mesh->get<Face> (originalConfusingBalls[i]);
     	std::vector<Node> nodes = currentFace.get<Node>();
       	for(unsigned int i=0; i<nodes.size(); i++){
	 		visitedVerts[nodes[i].id()] = true;
       	}
    	}
   	*/

	for(unsigned int i=0; i<vertsToVisit.size(); i++){
		visitedVerts[vertsToVisit[i].id()] = true;
	}

	while(!vertsToVisit.empty()){
		gmds::Node currentNode = vertsToVisit[vertsToVisit.size()-1];
		vertsToVisit.resize(vertsToVisit.size()-1);
		currentNode.get<Face>(neighbouringFaces);
		for(unsigned int i=0; i<neighbouringFaces.size(); i++){
			if((!visitedFaces[neighbouringFaces[i].id()])){
				math::Point center = neighbouringFaces[i].center();
				double current_dist = center.distance(sing_location);
				if(current_dist < newConfusingDistance) {
					if(current_dist >= previousDist){
						modifiedFaces.push_back(neighbouringFaces[i].id());
						m_faces_to_singularity_on_surf[neighbouringFaces[i].id()] = APnt;
					}
					vector<gmds::Node> currentFaceVerts;
					neighbouringFaces[i].get<Node>(currentFaceVerts); //std::vector<TCellID> currentFaceVerts = neighbouringFaces[i].getIDs<Node>();
					for(unsigned int j=0; j<currentFaceVerts.size(); j++){
						if(!visitedVerts[currentFaceVerts[j].id()]){
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

void SingularityGraphBuilder2D::getVertNeigh(){

	unsigned int vertNo = m_mesh->getNbNodes();
	vertNeighFaces.resize(vertNo, vector<gmds::TCellID>(0));

	for(auto f_id:m_mesh->faces()){
		gmds::Face currentFace = m_mesh->get<Face>(f_id);
		vector<Node> nodesTri = currentFace.get<Node>();
		for(unsigned int i=0; i < nodesTri.size(); i++){
			vertNeighFaces[nodesTri[i].id()].push_back(f_id);
		}
	}
   	// the n2n is simpler to go through edges
  	/*
  	try{
    		if(m_mesh->getModel().has(N2N))
    			cout<<"m_mesh->getModel().has(N2N)"<<endl;
  	}
  	catch(GMDSException& e){
    		cout << e.what() << std::endl;
  	}
  	try{
    		if(m_mesh->getModel().has(N2F))
    			cout<<"m_mesh->getModel().has(N2F)"<<endl;
    		for(auto n_id:m_mesh->nodes()){
      		gmds::Node n = m_mesh->get<Node>(n_id);
    		}
  	}
  	catch(GMDSException& e){
    		cout << e.what() << std::endl;
    		cout<<"no m_mesh->getModel().has(N2F)"<<endl;
  	}*/

}

/*----------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::removeSingularityLine(SingularityLine* ALine)
{
	/*
	if the line has one of its end points on the boundary, a geometric singularity
	has been added for it and it has also divided the boundary singularity line;
	therefore we have to remove the point created, merge back the 2 boundary singularity
	lines into one (and add it to the graph) and remove them*/
	vector<SurfaceSingularityLine* > surf_lines = m_graph.getSurfaceLines();
	// it should always be a SurfaceSingularityLine
	bool deletePoint = true;
	bool foundLine = false;
	for (unsigned int i = 0; i < surf_lines.size() && !foundLine; i++){
		SingularityLine* current_line = surf_lines[i];
		if(current_line==ALine){
			foundLine = true;
			vector<SingularityPoint* > singPointsOnSingLine = surf_lines[i]->getEndPoints();
			//cout<<"between ["<<singPointsOnSingLine[0]->getLocation().X()<<","<<singPointsOnSingLine[0]->getLocation().Y()<<"]"<<endl;
			//cout<<"and ["<<singPointsOnSingLine[1]->getLocation().X()<<","<<singPointsOnSingLine[1]->getLocation().Y()<<"]"<<endl;
			for(unsigned int t=0; t<2; t++){
				SingularityPoint* current_sp = singPointsOnSingLine[t];

				if(current_sp->getType() == 1){/* geom singularity point created on boundary (on edge or on vert)
					if on vert, check its marked as NodeOnPoint (if so, we don't have to do anything)*/
					vector<gmds::Node> firstNodes = current_sp->getMesh<gmds::Node>();

					if(firstNodes.size()>0){
						//cout<<"the end point of the singularity line is a vertex"<<endl;
						if(m_mesh->isMarked(firstNodes[0], m_mark_nodes_on_point)){
							deletePoint = false;
						}
					}

					if(deletePoint){
						vector<SingularityPoint*> m_points = m_graph.getPoints();

						for(unsigned int j=0; j<m_points.size(); j++){
							if(m_points[j]==current_sp){
								m_points.erase(m_points.begin()+j);
								break;
							}
						}
						vector<CurveSingularityLine* > curve_lines = m_graph.getCurveLines();
						bool foundCurveLine = false;
						unsigned int firstCurveLineIndex;
						//unsigned int current_spLineExtremity;
						bool foundBothLines = false;

						for(unsigned int j=0; j<curve_lines.size() && !foundBothLines; j++){
							vector<SingularityPoint* > singPointsOnCurveSingLine = curve_lines[j]->getEndPoints();
							//this should always be 2-sized (no intersections yet)
							for(unsigned int tt=0; tt<2; tt++){
								if(singPointsOnCurveSingLine[tt]==current_sp){

									if(foundCurveLine){
										//this is the second curveline that we find (as well as the number)
										// keep the orientation of the first curve that we find
										// singular point is always the 2nd on the 1st line
										vector<gmds::math::Point> first_pnts = curve_lines[firstCurveLineIndex]->getDiscretizationPoints();
										vector<gmds::Edge> first_curve_edges = curve_lines[firstCurveLineIndex]->getMeshEdges();
										vector<gmds::math::Point> second_pnts = curve_lines[j]->getDiscretizationPoints();
										vector<gmds::Edge> second_curve_edges = curve_lines[j]->getMeshEdges();

										if(tt==1){
											std::reverse(second_pnts.begin(),second_pnts.end());
											std::reverse(second_curve_edges.begin(),second_curve_edges.end());
										}
										second_pnts.erase(second_pnts.begin());
										second_curve_edges.erase(second_curve_edges.begin());
										for(unsigned int ttt=0; ttt<second_pnts.size(); ttt++)
											first_pnts.push_back(second_pnts[ttt]);

										for(unsigned int ttt=0; ttt<second_curve_edges.size(); ttt++)
											first_curve_edges.push_back(second_curve_edges[ttt]);

										curve_lines[firstCurveLineIndex]->setMeshEdges(first_curve_edges);
										curve_lines[firstCurveLineIndex]->setDiscretizationPoints(first_pnts);

										curve_lines[firstCurveLineIndex]->removeSingularityPoint();
										curve_lines[firstCurveLineIndex]->addSingularityPoint(singPointsOnCurveSingLine[(tt+1)%2]);
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

	if(!foundLine){
		if(withGlobalComments)
			cout<<"ALine is not a surface line"<<endl;
	}
}

/*--------------------------------------------------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::constructCenterTriangleCrosses()
{
	if(withGlobalComments)
		cout<<"constructCenterTriangleCrosses"<<endl;

	gmds::math::Point orig(0.0, 0.0, 0.0);
	if(triangle_centers.size()==0)
		triangle_centers.resize(original_faces_number, orig);

	if(triangle_centers_cross.size()==0)
		triangle_centers_cross.resize(original_faces_number);

	vector<TCoord> AWeights(3, (double)1/3);
	for(auto f_id:m_mesh->faces()){
		std::vector<math::Cross2D> Tricrosses;
		gmds::Face currentFace = m_mesh->get<Face>(f_id);
		vector<Node> nodesTri = currentFace.get<Node>();
		Tricrosses.push_back((*m_field)[nodesTri[0].id()]);
		Tricrosses.push_back((*m_field)[nodesTri[1].id()]);
		Tricrosses.push_back((*m_field)[nodesTri[2].id()]);
		vector<math::Vector3d> c_vectors0 = Tricrosses[0].componentVectors();
		vector<math::Vector3d> c_vectors1 = Tricrosses[1].componentVectors();
		vector<math::Vector3d> c_vectors2 = Tricrosses[2].componentVectors();
		triangle_centers[f_id] = currentFace.center();
		//mean is more like a median, not actually mean
		// triangle_centers_cross[f_id] = gmds::math::Cross2D::mean(Tricrosses, AWeights);
		gmds::math::Vector3d first_vect = c_vectors0[0];
		gmds::math::Vector3d second_closest_vect =  Tricrosses[1].closestComponentVector(first_vect);
		gmds::math::Vector3d third_closest_vect =  Tricrosses[2].closestComponentVector(first_vect);
		gmds::math::Vector3d center_vect = first_vect * AWeights[0] + second_closest_vect * AWeights[1] + third_closest_vect * AWeights[2];
		gmds::math::Vector3d center_vect_second = center_vect.getOneOrtho();

		triangle_centers_cross[f_id] = gmds::math::Cross2D(center_vect, center_vect_second);

	}

	gmds::Mesh centerTriMesh2(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));

	Variable<math::Vector3d>* crossVectC0 = centerTriMesh2.newVariable<math::Vector3d,GMDS_NODE>("crossVectC0");
	Variable<math::Vector3d>* crossVectC1 = centerTriMesh2.newVariable<math::Vector3d,GMDS_NODE>("crossVectC1");
	Variable<math::Vector3d>* crossVectC2 = centerTriMesh2.newVariable<math::Vector3d,GMDS_NODE>("crossVectC2");
	Variable<math::Vector3d>* crossVectC3 = centerTriMesh2.newVariable<math::Vector3d,GMDS_NODE>("crossVectC3");
	for(auto f_id:m_mesh->faces()){
		Face current = m_mesh->get<Face>(f_id);
		gmds::Node n = centerTriMesh2.newNode(triangle_centers[f_id].X(), triangle_centers[f_id].Y(), triangle_centers[f_id].Z());
		gmds::Face f  = centerTriMesh2.newTriangle(n, n, n);
		std::vector<math::Vector3d> crossV_centerTri = triangle_centers_cross[f_id].componentVectors();
		(*crossVectC0)[n.id()] = crossV_centerTri[0];
		(*crossVectC1)[n.id()] = crossV_centerTri[1];
		(*crossVectC2)[n.id()] = crossV_centerTri[2];
		(*crossVectC3)[n.id()] = crossV_centerTri[3];
	}

	gmds::IGMeshIOService meshIoServCenterTri2(&centerTriMesh2);
	gmds::VTKWriter writerCenterTri2(&meshIoServCenterTri2);
	writerCenterTri2.setCellOptions(gmds::N|gmds::F);
	writerCenterTri2.setDataOptions(gmds::N|gmds::F);

	std::stringstream file_nameCenterTri2;
	file_nameCenterTri2 <<m_output_directory_name<<"-CentreTri-crossVectors.vtk";
	writerCenterTri2.write(file_nameCenterTri2.str());

}

/*--------------------------------------------------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::constructOneCenterTriangleCross(gmds::Face&          ATriangle,
                                                                gmds::math::Point&   singleCenterTri,
                                                                gmds::math::Cross2D& singleCenterTriCross)
{
	gmds::math::Point orig(0.0, 0.0, 0.0);

	std::vector<math::Cross2D> Tricrosses;
	vector<TCoord> AWeights(3, (double)1/3);

	vector<Node> nodesTri = ATriangle.get<Node>();
	Tricrosses.push_back((*m_field)[nodesTri[0].id()]);
	Tricrosses.push_back((*m_field)[nodesTri[1].id()]);
	Tricrosses.push_back((*m_field)[nodesTri[2].id()]);
	vector<math::Vector3d> compV0 = (*m_field)[nodesTri[0].id()].componentVectors();
	vector<math::Vector3d> compV1 = (*m_field)[nodesTri[1].id()].componentVectors();
	vector<math::Vector3d> compV2 = (*m_field)[nodesTri[2].id()].componentVectors();

	singleCenterTri = ATriangle.center();

	TCoord ref_angle =  Tricrosses[0].referenceAngle();

	gmds::math::Vector3d ref_vector =  Tricrosses[0].referenceVector();

	TCoord pen_angle=0.0;
	for (unsigned int i = 0; i <3; i++){
		const TCoord pen_i = Tricrosses[i].referenceAngle();
		pen_angle += pen_i;
	}
	pen_angle /=Tricrosses.size();

	ref_angle  = gmds::math::modulo2PI(pen_angle);
	ref_vector = gmds::math::Cross2D(ref_angle).referenceVector();

	singleCenterTriCross = gmds::math::Cross2D(ref_angle);

}

/*---------------------------------------------------------------------------*/

int SingularityGraphBuilder2D::getShortestPathBtwNewNodes(unsigned int&                                             numberOfNewNodes,
                                                          gmds::TCellID&                                            source,
                                                          vector<unsigned int>&                                     targets,
                                                          vector<pair<double, unsigned int>>&                       min_distance,
                                                          vector<int>&                                              previous,
                                                          int&                                                      found,
                                                          vector<vector<gmds::TCellID>>&                            newNode2NewNodeNeighbours,
                                                          vector<vector<double>>&                                   face2FaceTransport,
                                                          vector<vector<unsigned int>>&                             face2FaceDeviation,
                                                          pair<gmds::math::Vector3d, gmds::math::Vector3d>&         prevDirCross,
                                                          vector<pair<gmds::math::Vector3d, gmds::math::Vector3d>>& targetDirCross,
                                                          vector<bool>&                                             forbiddenFaces,
                                                          bool&                                                     targetBdry,
                                                          gmds::math::Point&                                        startPoint,
                                                          gmds::math::Point&                                        endPoint,
                                                          double&                                                   maxDist)
{
	/*function detected the shortest paths between a source triangle and a vector of target triangles; walk through triangle centers and also through nodes */
	if(withGlobalComments)
		cout<<"getShortestPathBtwNewNodes between source "<<source<<" and "<<targets[0]<<endl;
	min_distance.clear();

	min_distance.resize(numberOfNewNodes, make_pair(maxDist, 1));
	min_distance[source] = make_pair(0.0,1);
	previous.clear();
	previous.resize(numberOfNewNodes, -1);

	double turnPenalizationTerm = 0.5 * maxDist;

	std::set<std::pair<double,unsigned int> > vertex_queue;
	math::Vector3d closestToPrevDir, closestToPrevCross, prevDir;

	vector<pair<gmds::math::Vector3d, gmds::math::Vector3d>> tempPrevDirCross(numberOfNewNodes, prevDirCross);

	for(unsigned int i=0; i<numberOfNewNodes; i++){
		vertex_queue.insert(std::make_pair(min_distance[i].first,0));
	}
	vertex_queue.erase(std::make_pair(min_distance[source].first,0));
	vertex_queue.insert(std::make_pair(0.0, source));
	cout<<"source "<<source<<endl;
	std::vector<unsigned int>::iterator it;
	gmds::Face v_face, u_face;
	gmds::Node v_node, u_node;
	gmds::TCellID u_id, v_id;

	vector<bool> visitedNodesAndFaces(numberOfNewNodes, false);

	u_id = vertex_queue.begin()->second;
	if(u_id<original_nodes_number)
		u_node = m_mesh->get<gmds::Node>(vertex_queue.begin()->second);
	else
		u_face = m_mesh->get<gmds::Face>(vertex_queue.begin()->second - original_nodes_number);

	double total_dist = vertex_queue.begin()->first;
	unsigned int visitedFacesNo = min_distance[vertex_queue.begin()->second].second;

	vertex_queue.erase(vertex_queue.begin());
	visitedNodesAndFaces[u_id] = true;

	it = std::find (targets.begin(), targets.end(), u_id);
	if(it!=targets.end()){

		found = std::distance( targets.begin(), it );
		prevDir =  tempPrevDirCross[u_id].first;
		math::Vector3d prevDirCross =  tempPrevDirCross[u_id].second;
		if(!targetBdry){
			if((targetDirCross[it - targets.begin()].second).angle(prevDirCross)>=gmds::math::Constants::PIDIV4)
				min_distance[u_id].first =  total_dist + maxDist + 1.0;
			else{
				min_distance[u_id].first =  total_dist + 10*(targetDirCross[it - targets.begin()].second).angle(prevDirCross);
				return u_id;
			}
		}
		else{
			double intersectionParam;
			math::Ray from_ray(startPoint, startPoint + prevDirCross);

			vector<gmds::Edge> currentEdges = u_face.get<gmds::Edge>();
			bool hasBdryEdge = false;
			for(unsigned int i=0; i<currentEdges.size(); i++){
				if (m_mesh->isMarked(currentEdges[i], m_mark_edges_on_curve)){
					hasBdryEdge = true;
					vector<gmds::Node> currentNodes =  currentEdges[i].get<Node>();
					math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
					if(from_ray.SecondMetIntersect2D(oppSeg, endPoint, intersectionParam, temp_epsilon)){
					  	double tempAngleVal = 1.0 - fabs(prevDirCross.dot(bdry_edge_normals[currentEdges[i].id()]));

						if(tempAngleVal<gmds::math::Constants::PIDIV4){
							min_distance[u_id].first =  total_dist + tempAngleVal;
							return u_id;
						}
					}
				}
			}
			if(!hasBdryEdge){
				// VERTEX-> angle with vert normal
				vector<gmds::Node> currentNodes = u_face.get<gmds::Node>();
				for(unsigned int i=0; i<currentNodes.size(); i++){
					if(m_mesh->isMarked(currentNodes[i], m_mark_nodes_on_point) ||
					m_mesh->isMarked(currentNodes[i], m_mark_nodes_on_curve)){
 						vector<gmds::Edge> currentEdges = currentNodes[i].get<gmds::Edge>();
						for(unsigned int j=0; j<currentEdges.size(); j++){
							if(m_mesh->isMarked(currentEdges[j], m_mark_edges_on_curve)){
								vector<gmds::Node> edge_nodes = currentEdges[j].get<gmds::Node>();
								math::Segment oppSeg(edge_nodes[0].getPoint(), edge_nodes[1].getPoint());
								if(from_ray.SecondMetIntersect2D(oppSeg, endPoint, intersectionParam, temp_epsilon)){
									double tempAngleVal = 1.0 - fabs(prevDirCross.dot(bdry_edge_normals[currentEdges[j].id()]));
									if(tempAngleVal<gmds::math::Constants::PIDIV4){
										min_distance[u_id].first =  total_dist + tempAngleVal;
										return u_id;
									}
								}
							}
						}
					}
				}
			}
		}
	}

	bool toCheck = false;
	if(is_bdry_face[u_id])
		toCheck = true;

	for(unsigned int i=0; i<newNode2NewNodeNeighbours[u_id].size(); i++){
		v_id = newNode2NewNodeNeighbours[u_id][i];
		if(v_id<original_nodes_number)
			v_node = m_mesh->get<Node>(v_id);
		else
			v_face = m_mesh->get<Face>(v_id - original_nodes_number);

		bool validNeighbour = true;
		if(toCheck==true){
			if(is_bdry_face[v_id]){
				math::Point intersectionPnt;
				double intersectionParam;
				math::Ray from_ray(triangle_centers[u_id], triangle_centers[v_id]);
				vector<gmds::Edge> currentEdges = u_face.get<gmds::Edge>();
				for(unsigned int j=0; j<3; j++){
					if (m_mesh->isMarked(currentEdges[j], m_mark_edges_on_curve)){
						vector<gmds::Node> currentNodes =  currentEdges[j].get<Node>();
						math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
						if(from_ray.SecondMetIntersect2D(oppSeg, intersectionPnt, intersectionParam, temp_epsilon)){
							if((intersectionParam>0.000000001)||(intersectionParam<0.999999999)){
								validNeighbour = false;
							}
						}
					}
				}
			}
		}

		if((!visitedNodesAndFaces[v_id])&&(!forbiddenFaces[v_id])&&(validNeighbour)){

			math::Vector3d tri2tri(startPoint, triangle_centers[v_id]);
			tri2tri.normalize();

			math::Cross2D crossV = triangle_centers_cross[v_id];
			prevDir =  tempPrevDirCross[u_id].first;
			closestToPrevDir =  crossV.closestComponentVector(prevDir);
			closestToPrevCross =  crossV.closestComponentVector(tempPrevDirCross[u_id].second);

			double distance_through_u;

			if(tri2tri.angle(tempPrevDirCross[u_id].second)>=gmds::math::Constants::PIDIV4)
				distance_through_u =  total_dist + maxDist + closestToPrevCross.angle(tri2tri);
			else
				distance_through_u =  total_dist + //(tempPrevDirCross[u_id].second).angle(tri2tri)
				                      + closestToPrevCross.angle(tri2tri)
				                      + (floor)((tempPrevDirCross[u_id].second.angle(closestToPrevCross))/gmds::math::Constants::PIDIV2)*10.0;

			if (distance_through_u/(visitedFacesNo+1) < min_distance[v_id].first/min_distance[v_id].second) {

				vertex_queue.erase(std::make_pair(min_distance[v_id].first,v_id));
				min_distance[v_id].first = distance_through_u;
				min_distance[v_id].second = visitedFacesNo + 1;
				previous[v_id] = u_id;
				vertex_queue.insert(std::make_pair(min_distance[v_id].first,v_id));
				tempPrevDirCross[v_id].first = tri2tri;
				tempPrevDirCross[v_id].second = closestToPrevCross;
			}
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////

	while (!vertex_queue.empty()){

		u_id = vertex_queue.begin()->second;
		if(u_id<original_nodes_number)
			u_node = m_mesh->get<gmds::Node>(vertex_queue.begin()->second);
		else
			u_face = m_mesh->get<gmds::Face>(vertex_queue.begin()->second - original_nodes_number);

		double total_dist = vertex_queue.begin()->first;
		unsigned int visitedFacesNo = min_distance[vertex_queue.begin()->second].second;

		vertex_queue.erase(vertex_queue.begin());
		visitedNodesAndFaces[u_id] = true;

		it = std::find (targets.begin(), targets.end(), u_id);
		if(it!=targets.end()){

			found = std::distance( targets.begin(), it );
			prevDir =  tempPrevDirCross[u_id].first;
			math::Vector3d prevDirCross =  tempPrevDirCross[u_id].second;
			if(!targetBdry){
				if((targetDirCross[it - targets.begin()].second).angle(prevDirCross)>=gmds::math::Constants::PIDIV4)
					min_distance[u_id].first =  total_dist + turnPenalizationTerm + 1.0;
				else{
					min_distance[u_id].first =  total_dist + 10*(targetDirCross[it - targets.begin()].second).angle(prevDirCross);
					return u_id;
				}
			}
			else{
				double intersectionParam;
				math::Ray from_ray(triangle_centers[u_id], prevDirCross);

				vector<gmds::Edge> currentEdges = u_face.get<gmds::Edge>();
				bool hasBdryEdge = false;
				for(unsigned int i=0; i<currentEdges.size(); i++){
					if (m_mesh->isMarked(currentEdges[i], m_mark_edges_on_curve)){
						hasBdryEdge = true;
						vector<gmds::Node> currentNodes =  currentEdges[i].get<Node>();
						math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
						if(from_ray.SecondMetIntersect2D(oppSeg, endPoint, intersectionParam, temp_epsilon)){
							double tempAngleVal = 1.0 - fabs(prevDirCross.dot(bdry_edge_normals[currentEdges[i].id()]));

							if(tempAngleVal<gmds::math::Constants::PIDIV4){
								min_distance[u_id].first =  total_dist + tempAngleVal;
								return u_id;
							}
						}
					}
				}
				if(!hasBdryEdge){
					// VERTEX-> angle with vert normal
					vector<gmds::Node> currentNodes = u_face.get<gmds::Node>();
					for(unsigned int i=0; i<currentNodes.size(); i++){
						if(m_mesh->isMarked(currentNodes[i], m_mark_nodes_on_point) ||
						m_mesh->isMarked(currentNodes[i], m_mark_nodes_on_curve)){
 							vector<gmds::Edge> currentEdges = currentNodes[i].get<gmds::Edge>();
							for(unsigned int j=0; j<currentEdges.size(); j++){
								if(m_mesh->isMarked(currentEdges[j], m_mark_edges_on_curve)){
									vector<gmds::Node> edge_nodes = currentEdges[j].get<gmds::Node>();
									math::Segment oppSeg(edge_nodes[0].getPoint(), edge_nodes[1].getPoint());
									if(from_ray.SecondMetIntersect2D(oppSeg, endPoint, intersectionParam, temp_epsilon)){
										double tempAngleVal = 1.0 - fabs(prevDirCross.dot(bdry_edge_normals[currentEdges[j].id()]));

										if(tempAngleVal<gmds::math::Constants::PIDIV4){
											min_distance[u_id].first =  total_dist + tempAngleVal;
											return u_id;
										}
									}
								}
							}
						}
					}
				}
			}
		}

		toCheck = false;
		if(is_bdry_face[u_id])
		  	toCheck = true;


		for(unsigned int i=0; i<newNode2NewNodeNeighbours[u_id].size(); i++){
			v_id = newNode2NewNodeNeighbours[u_id][i];
			if(v_id<original_nodes_number)
				v_node = m_mesh->get<Node>(v_id);
			else
				v_face = m_mesh->get<Face>(v_id - original_nodes_number);

			bool validNeighbour = true;
			if(toCheck==true){
				if(is_bdry_face[v_id]){
					math::Point intersectionPnt;
					double intersectionParam;
					math::Ray from_ray(triangle_centers[u_id], triangle_centers[v_id]);
					vector<gmds::Edge> currentEdges = u_face.get<gmds::Edge>();
					for(unsigned int j=0; j<3; j++){
						if (m_mesh->isMarked(currentEdges[j], m_mark_edges_on_curve)){
							vector<gmds::Node> currentNodes =  currentEdges[j].get<Node>();
							math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
							if(from_ray.SecondMetIntersect2D(oppSeg, intersectionPnt, intersectionParam, temp_epsilon)){
								if((intersectionParam>0.000000001)||(intersectionParam<0.999999999)){
									validNeighbour = false;
								}
							}
						}
					}
				}
			}
			if((!visitedNodesAndFaces[v_id])&&(!forbiddenFaces[v_id])&&(validNeighbour)){

				math::Vector3d tri2tri(triangle_centers[u_id], triangle_centers[v_id]);
				tri2tri.normalize();
				math::Cross2D crossV = triangle_centers_cross[v_id];
				prevDir =  tempPrevDirCross[u_id].first;
				closestToPrevDir =  crossV.closestComponentVector(prevDir);
				closestToPrevCross =  crossV.closestComponentVector(tempPrevDirCross[u_id].second);
				double distance_through_u;

				if(tri2tri.angle(tempPrevDirCross[u_id].second)>=gmds::math::Constants::PIDIV4){
					distance_through_u =  total_dist + turnPenalizationTerm + closestToPrevCross.angle(tri2tri);
				}
 				else{
					distance_through_u = total_dist + //(tempPrevDirCross[u_id].second).angle(tri2tri)
					                     + closestToPrevCross.angle(tri2tri)
					                     + (floor)((tempPrevDirCross[u_id].second.angle(closestToPrevCross))/gmds::math::Constants::PIDIV2)*10.0;

				}

				//WARNING division par visitedFacesNo not particularly ok; especially for very irregular meshes
				if (distance_through_u/(visitedFacesNo+1) < min_distance[v_id].first/min_distance[v_id].second) {
				//if (distance_through_u < min_distance[v_id].first) {
					vertex_queue.erase(std::make_pair(min_distance[v_id].first,v_id));

					min_distance[v_id].first = distance_through_u;
					min_distance[v_id].second = visitedFacesNo + 1;
					previous[v_id] = u_id;
					vertex_queue.insert(std::make_pair(min_distance[v_id].first,v_id));

					tempPrevDirCross[v_id].first = tri2tri;
					tempPrevDirCross[v_id].second = closestToPrevCross;

				}
			}
		}
	}

	//cout<<"we should never get here"<<endl;
	return -1;
}

/*---------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::createSingularityLinesSimultaneousStart()
{
	//========================================================================
	/*We start simultaneously from all the singularity points (through all slots) to compute streamlines in a
	* direction following as close as possible the cross field. At each iteration all such streamlines
	* advance until they have reached a certain geometric length step (searchStep) - which
	* is increased with each iteration. If such a streamline reaches the boundary, the respective line is
	* logged and the slot is being frozen.
	* If two streamlines (departing from two different singularity points) advance until they are
	* within a certain distance from each other (thresholdStreamLineDist), they will be connected only
	* if their directions (at the tip of the streamlines) are close to parallel
	*/
	//========================================================================

	cout<<"createSingularityLinesSimultaneousStart"<<endl;
	/* when deciding whether to connect 2 singularity lines we have two options:
	1. connectByField = false; connect depending on the geometrical distance between the extremities of the singularity lines and depending on the angle made by those
	2. connectByField = true; connect depending on the frame field; consider the straight line between the extremities of the 2 singularity lines - this line will intersect edges (and vertices rarely); calculate the closest component vector at the intersection point with regard to the straight line's direction; if the closest component vector coincides (end to end) - we can connect; (by coincides -  consider direction 0-2 and direction 1-3)

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
	bool connectByField = false;

	double searchStep; // search step for streamlines
	double thresholdStreamLineDist; // if 2 lines are within this threshold from one another => connect them
	double weightAngle = 0.75;
	double weightDistance = 0.25;

	searchStep = mean_edge_length;
	//searchStep = 3*searchStep; // WARNING searchStep should be sufficiently small
	//WARNING thresholdStreamLineDist should be sufficiently big
	thresholdStreamLineDist = 3*searchStep;

	vector<SingularityPoint::Slot*> searchSlots;
	vector<SingularityPoint*> singularity_points = m_graph.getPoints();

	for (unsigned int i = 0; i < singularity_points.size(); i++) {
		SingularityPoint* pi = singularity_points[i];
		if(pi->getType()==0){
			std::vector<SingularityPoint::Slot*> pi_slots = pi->getSlots();
			for (unsigned int j = 0; j < pi_slots.size(); j++) {
				if (!pi_slots[j]->isLaunched)
					searchSlots.push_back(pi_slots[j]);
			}
		}
	}

	unsigned int NoSimLines = searchSlots.size();
	vector<vector<math::Point>> line_discretization(NoSimLines, vector<math::Point>(0)), copy_line_discretization(NoSimLines, vector<math::Point>(0));
	vector<vector<TCellID>>     line_triangles(NoSimLines, vector<TCellID>(0)), copy_line_triangles(NoSimLines, vector<TCellID>(0));
	vector<double> accumulatedDistancePerSlotLine(NoSimLines, 0.0);
	bool 					find_end_bdry = false;
	bool 					end_on_free_slot = false;
	vector<TCellID> 		to_cell_id(NoSimLines);
	vector<int> 			to_cell_dim(NoSimLines);
	vector<math::Vector3d> 	to_dir(NoSimLines);
	vector<double>            streamlineDeviation(NoSimLines, 0.0);

	double currentAccDist = 0.0;
	unsigned int noFrozenSLots = 0;
	double currentSearchStep;

	for(unsigned int i=0; i<NoSimLines; i++){
		line_discretization[i].push_back(searchSlots[i]->from_point->getLocation());
		line_discretization[i].push_back(searchSlots[i]->location);
		to_cell_id[i]  = searchSlots[i]->starting_cell_id ;
		to_cell_dim[i] = searchSlots[i]->starting_cell_dim;
		to_dir[i] = searchSlots[i]->direction; //starting direction
	}

	currentSearchStep = searchStep;

	int contor = 0;

	math::Vector3d surfNormal(1.0, 1.0, 1.0);
	// vector storing ((slotA, slotB) , (distance between them, angle of intersection))
	vector<pair<pair<unsigned int, unsigned int>, pair<double, double>>> possibleConnectingLines;
	vector<pair<unsigned int, unsigned int>> forbiddenConnectingLines;
	double distBtwSlots;
	/*
	rather than adding the streamlines as soon as two of them meet (or meet the boundary),
	i will add them to the possibleConnectingLines and decide later;
	possibleConnectingLines contains also the (geometrical) distance and the angle in between the two;
	for bry verts add -1 as second vertex
	for bry lines - if no other streamline has been met in the following two iterations
	(with good angles and distances = to define what "good" is), construct BdryLine
	for surface line - if no other streamline has been met in the following two iterations
	(with better angles and distances), connect; however, stop growLine() for the slot at the
	iteration where it reaches the vicinity of another slot line*/

	vector<bool> stopIncrease(NoSimLines, false);
	unsigned int testCont=0;

	while(noFrozenSLots<NoSimLines){
		cout<<"noFrozenSLots= "<<noFrozenSLots<<endl;
		contor++;
		possibleConnectingLines.clear();
		stopIncrease.clear();
		stopIncrease.resize(NoSimLines, false);
		for(unsigned int i=0; i<NoSimLines; i++){

			SingularityPoint       *to_sing_pnt  = 0;
			SingularityPoint::Slot *to_slot = 0;
			SingularityPoint::Slot* current_slot = searchSlots[i];
			if((!current_slot->isFreeze)&&(!stopIncrease[i])){
				find_end_bdry = false;
				end_on_free_slot = false;
				while((!find_end_bdry)&&(!end_on_free_slot)&&(accumulatedDistancePerSlotLine[i]<currentSearchStep)){

					/* grow line until its length equals currentSearchStep or until find_end_bdry or until end_on_free_slot*/
					growLine(current_slot->from_point,
								current_slot,
								to_sing_pnt,
								to_slot,
								line_discretization[i][line_discretization[i].size()-1],
								to_dir[i],
								line_discretization[i],
								line_triangles[i],
								to_cell_dim[i],
								to_cell_id[i],
								streamlineDeviation[i],
								accumulatedDistancePerSlotLine[i],
								currentSearchStep,
								find_end_bdry,
								end_on_free_slot);

					if(to_slot!=0){
						cout<<"current_slot->location "<<current_slot->location<<"; to_slot->location "<<to_slot->location<<endl;
					}
					if(find_end_bdry){
						//======================================================================
						// CASE 1 - We finish on the boundary. A geometric point must be created
						// The cell defined by (start_cell_dim[i], start_cell_id[i]) is on the boundary.
						// We have to create a geometric singularity point.
						// however we could do this at the end of the entire iteration(for the current searchStep)
						//line creation

						SurfaceSingularityLine* surf_line = m_graph.newSurfaceLine();
						int sepNumberTmp = m_graph.getNbLines();
						surf_line->setNumber(sepNumberTmp);
						//connect line to initial singularity point
						SingularityPoint* from_sing_pnt = current_slot->from_point;
						surf_line->addSingularityPoint(from_sing_pnt);
						surf_line->addDiscretizationPoint(from_sing_pnt->getLocation());
						current_slot->line = surf_line;

						math::Vector3d firstDir = math::Vector3d(from_sing_pnt->getLocation(), line_discretization[i][0]);

						current_slot->line_direction = firstDir;
						current_slot->isLaunched = true;
						if(to_slot!=0){
							streamlineDeviation[i] = (streamlineDeviation[i] + fabs(1 - to_slot->direction.dot(math::Vector3d(line_discretization[i][line_discretization[i].size()-1],line_discretization[i][line_discretization[i].size()-2]))));
							streamlineDeviation[i] = streamlineDeviation[i]/(line_discretization.size()+1);
							current_slot->lineDeviation = streamlineDeviation[i];
						}
						//Insertion of line points
						for(unsigned int j=0; j<line_discretization[i].size(); j++) {
							surf_line->addDiscretizationPoint(line_discretization[i][j]);
						}

						for(unsigned int j=0; j<line_triangles[i].size(); j++) {
							surf_line->addTraversedFace(line_triangles[i][j]);
						}

						SingularityPoint* geom_pnt;
						SingularityPoint::Slot* incoming_slot;
						createGeometricSingularityPoint(line_discretization[i].back(),      // the last point added
						                                to_dir[i],      // the direction we come from
						                                to_cell_dim[i], // the dim. of the cell
						                                to_cell_id[i],  // the id of the cell
						                                geom_pnt,       // the created point
						                                incoming_slot); // and the slot

						current_slot->lineDeviation = streamlineDeviation[i];

						surf_line->addSingularityPoint(geom_pnt);
						noFrozenSLots++;
						current_slot->isFreeze = true;
						current_slot->isLaunched = true;

						for(unsigned int j=0; j<possibleConnectingLines.size(); j++){
							if((i==possibleConnectingLines[j].first.first)||((i==possibleConnectingLines[j].first.second))){
								possibleConnectingLines.erase(possibleConnectingLines.begin()+j);
							}
						}

						gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));

						for (unsigned int j = 0; j < line_discretization[i].size()-1; j++){
							gmds::Node n1 = m.newNode(line_discretization[i][j].X(), line_discretization[i][j].Y(), line_discretization[i][j].Z());
							gmds::Node n2 = m.newNode(line_discretization[i][j + 1].X(), line_discretization[i][j + 1].Y(), line_discretization[i][j + 1].Z());
							gmds::Face f  =  m.newTriangle(n1, n1, n2);
						}
						gmds::IGMeshIOService ioService(&m);
						gmds::VTKWriter vtkWriter(&ioService);
						vtkWriter.setCellOptions(gmds::N|gmds::F);
						vtkWriter.setDataOptions(gmds::N|gmds::F);
						std::string file_name = "SimultaneousRK4_SingToBdry"+to_string(i)+".vtk";
						vtkWriter.write(file_name);
					}
					else{// not find_end_bdry
						if(end_on_free_slot){

							SurfaceSingularityLine* surf_line = m_graph.newSurfaceLine();
							int sepNumberTmp = m_graph.getNbLines();
							surf_line->setNumber(sepNumberTmp);

							SingularityPoint* from_sing_pnt = current_slot->from_point;
							surf_line->addSingularityPoint(from_sing_pnt);
							surf_line->addDiscretizationPoint(from_sing_pnt->getLocation());

							current_slot->line = surf_line;

							math::Vector3d firstDir = math::Vector3d(from_sing_pnt->getLocation(), line_discretization[i][0]);

							current_slot->line_direction = firstDir;
							current_slot->isLaunched = true;

							current_slot->lineDeviation = streamlineDeviation[i];

							for(unsigned int j=0; j<line_discretization[i].size(); j++) {
								surf_line->addDiscretizationPoint(line_discretization[i][j]);
							}

							for(unsigned int j=0; j<line_triangles.size(); j++) {
								surf_line->addTraversedFace(line_triangles[i][j]);
							}

							to_slot->isLaunched = true;
							current_slot->isFreeze = true;
							to_slot->isFreeze = true;

							to_slot->line = surf_line;
							to_slot->line_direction = to_dir[i];

							to_slot->lineDeviation = streamlineDeviation[i];

							surf_line->addSingularityPoint(to_sing_pnt);
							surf_line->addDiscretizationPoint(to_sing_pnt->getLocation());

							noFrozenSLots = noFrozenSLots + 2;

							gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
							for (unsigned int j = 0; j < line_discretization[i].size()-1; j++){
								gmds::Node n1 = m.newNode(line_discretization[i][j].X(), line_discretization[i][j].Y(), line_discretization[i][j].Z());
								gmds::Node n2 = m.newNode(line_discretization[i][j + 1].X(), line_discretization[i][j + 1].Y(), line_discretization[i][j + 1].Z());
								gmds::Face f  =  m.newTriangle(n1, n1, n2);
							}
							gmds::IGMeshIOService ioService(&m);
							gmds::VTKWriter vtkWriter(&ioService);
							vtkWriter.setCellOptions(gmds::N|gmds::F);
							vtkWriter.setDataOptions(gmds::N|gmds::F);
							std::string file_name = "SimultaneousRK4_SingToSingBall"+to_string(i)+"-x.vtk";
							vtkWriter.write(file_name);
							cout<<"has written SingToSing "<<file_name<<endl;
						}
						else{
							cout<<"got in the vecinity of another line"<<endl;
							//======================================================================
							//CASE 2 - We encounter another line within the thresholdStreamLineDist
							// for now, we will simply connect the last points added to the 2 lines (which are within thresholdStreamLineDist); this can and should be improved
							// connect only if deviation for the connecting line is sufficiently small; otherwise,let it to arrive to boundary
							for(unsigned int j=0; j<NoSimLines; j++){
								if((searchSlots[i]->from_point!=searchSlots[j]->from_point)&&(line_discretization[j].size()>2)&&(!searchSlots[j]->isFreeze)){
									distBtwSlots = line_discretization[i][line_discretization[i].size()-1].distance(line_discretization[j][line_discretization[j].size()-1]);
									if(distBtwSlots<=thresholdStreamLineDist){
										bool forbiddenPair = false;
										for(unsigned int ti0=0; ti0<forbiddenConnectingLines.size(); ti0++){
											if((forbiddenConnectingLines[ti0].first==i)&&(forbiddenConnectingLines[ti0].second==j)){
												forbiddenPair = true;
												break;
											}
										}
										if(!forbiddenPair){
											math::Point last_added_pnt = line_discretization[i][line_discretization[i].size()-1];
											math::Vector3d start_dir = math::Vector3d(last_added_pnt, line_discretization[j][line_discretization[j].size()-1]);
											if(!connectByField){
												double tempDev = 0.0;

												vector<gmds::Face> candidate_faces;
												vector<bool> visitedFaces(original_faces_number,false), addedFaces(original_faces_number,false);
												if(to_cell_dim[i]==0){
													gmds::Node currentNode = m_mesh->get<gmds::Node>(to_cell_id[i]);
													candidate_faces = currentNode.get<gmds::Face>();
												}
												else{
													if(to_cell_dim[i]==1){
														gmds::Edge currentEdge = m_mesh->get<gmds::Edge>(to_cell_id[i]);
														candidate_faces = currentEdge.get<gmds::Face>();
													}
													else{
														gmds::Face currentFace = m_mesh->get<gmds::Face>(to_cell_id[i]);
														candidate_faces = currentFace.get<gmds::Face>();
													}
												}

												for(unsigned int t=0; t<candidate_faces.size(); t++){
													visitedFaces[candidate_faces[t].id()] = true;
												}

												for(unsigned int t=0; t<line_triangles[i].size(); t++)
													addedFaces[line_triangles[i][t]] = true;
												for(unsigned int t=0; t<line_triangles[j].size(); t++)
													addedFaces[line_triangles[j][t]] = true;

												math::Segment seg1(last_added_pnt, line_discretization[j][line_discretization[j].size()-1]);
												math::Ray from_ray(last_added_pnt,line_discretization[j][line_discretization[j].size()-1]);
												math::Vector3d connLineDir =	from_ray.getDirUnit();
												math::Vector3d closest2Current, closest0, closest1;
												unsigned int t = 0;
												bool stopAdding = false;
												vector<bool> visitedEdges(m_mesh->getNbEdges(),false);
												unsigned int noIntEdges = 0;

												copy_line_triangles[i] = line_triangles[i];

												while(t<candidate_faces.size()){
													Face currentFace = candidate_faces[t];
													vector<gmds::Edge> currentEdges = currentFace.get<gmds::Edge>();
													for(unsigned int tt=0; tt<3; tt++){
														vector<gmds::Node> currentNodes = currentEdges[tt].get<gmds::Node>();
														math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
														math::Point intersectionPnt;
														double intersectionParam;
														if(seg1.SecondMetIntersect2D(oppSeg, intersectionPnt, intersectionParam, temp_epsilon)){
															if(from_ray.SecondMetIntersect2D(oppSeg, intersectionPnt, intersectionParam, temp_epsilon)){
																if(!visitedEdges[currentEdges[tt].id()]){
																	math::Cross2D cross_0 = (*m_field)[currentNodes[0].id()];
																	math::Cross2D cross_1  = (*m_field)[currentNodes[1].id()];
																	closest0 =  cross_0.closestComponentVector(connLineDir);
																	closest1 =  cross_1.closestComponentVector(connLineDir);
																	closest2Current = intersectionParam * closest0 + (1 - intersectionParam) * closest1;
																	tempDev = tempDev + fabs(1.0 - connLineDir.dot(closest2Current));
																	visitedEdges[currentEdges[tt].id()] = true;
																	noIntEdges++;
																}

																if(currentFace.id() == line_triangles[j][line_triangles[j].size()-1])
																	stopAdding = true;
																if(!addedFaces[currentFace.id()]){
																	copy_line_triangles[i].push_back(currentFace.id());
																	addedFaces[currentFace.id()] = true;
																}
																if(!stopAdding){
																	vector<Face> adj_faces0 = currentNodes[0].get<Face>();
																	vector<Face> adj_faces1 = currentNodes[1].get<Face>();
																	adj_faces0.insert(adj_faces0.end(), adj_faces1.begin(), adj_faces1.end());
																	for(unsigned int ttt=0; ttt<adj_faces0.size(); ttt++){
																		if(!visitedFaces[adj_faces0[ttt].id()]){
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

												math::Vector3d test1(line_discretization[i][line_discretization[i].size()-2], line_discretization[i][line_discretization[i].size()-1]);
												math::Vector3d test2(line_discretization[j][line_discretization[j].size()-2], line_discretization[j][line_discretization[j].size()-1]);
												test1.normalize();
												test2.normalize();
												// if close to perpendicular=> shouldn't connect
												tempDev = test1.angle(test2); // tempDevin [0,PI]]

												math::Vector3d   tempCross = test1.cross(test2);
												// we connect only if the dot product with the normal (in 2d case - (1,1,1) is higher that 0) - in order to account for orientation
												if(tempCross.dot(surfNormal)>=0){
													//we want (test1.angle(test2))/(double)(math::Constants::PI)) to be in between [0.75 , 1.25]; if so -> we can connect 2 streamlines, otherwise not
													tempDev = tempDev/(double)(math::Constants::PI);
													if((tempDev>=0.75)&&(tempDev<=1.25)){
														possibleConnectingLines.push_back(make_pair(make_pair(i, j),make_pair(distBtwSlots, std::fabs(1.0 - tempDev))));
														//ideally the angle btw the lines should be 0; therefore tempDev = 1
														stopIncrease[j] = true;
													}
												}
											}
											else{

											}
										}
									}
								}
							}
						}
					}

					gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
					for(unsigned int ti=0;ti<NoSimLines;ti++){
						for(unsigned int ti1=0;ti1<line_discretization[ti].size()-1;ti1++){
							gmds::Node n1 = m.newNode(line_discretization[ti][ti1].X(), line_discretization[ti][ti1].Y(), line_discretization[ti][ti1].Z());
							gmds::Node n2 = m.newNode(line_discretization[ti][ti1+ 1].X(), line_discretization[ti][ti1 + 1].Y(), line_discretization[ti][ti1 + 1].Z());
							gmds::Face f  =  m.newTriangle(n1, n1, n2);
						}
					}

					gmds::IGMeshIOService ioService(&m);
					gmds::VTKWriter vtkWriter(&ioService);
					vtkWriter.setCellOptions(gmds::N|gmds::F);
					vtkWriter.setDataOptions(gmds::N|gmds::F);
					std::string test_file_name = m_output_directory_name +"SimultaneousRK4-eachStep"+to_string(testCont)+".vtk";
					vtkWriter.write(test_file_name);
					testCont++;
				}
			}
		}

		// if one slot is within radius with 2 other slots, choose the one with the lowest deviation
		if(possibleConnectingLines.size()>0){
			unsigned int t1=0;
			unsigned int t2;
			double connectionTerm1, connectionTerm2;

			while(t1<possibleConnectingLines.size()-1){
				t2 = t1+1;
				while(t2<possibleConnectingLines.size()){
					if(((possibleConnectingLines[t1].first.first==possibleConnectingLines[t2].first.first)&&
					(possibleConnectingLines[t1].first.second==possibleConnectingLines[t2].first.second))||
					((possibleConnectingLines[t1].first.first==possibleConnectingLines[t2].first.second)&&
					(possibleConnectingLines[t1].first.second==possibleConnectingLines[t2].first.first))){

						connectionTerm1 = weightDistance * (possibleConnectingLines[t1].second.first/(possibleConnectingLines[t1].second.first + possibleConnectingLines[t2].second.first));
						connectionTerm1 = connectionTerm1 + weightAngle * (possibleConnectingLines[t1].second.second/(possibleConnectingLines[t1].second.second + possibleConnectingLines[t2].second.second));

						connectionTerm2 = weightDistance * (possibleConnectingLines[t2].second.first/(possibleConnectingLines[t1].second.first + possibleConnectingLines[t2].second.first));
						connectionTerm2 = connectionTerm2 + weightAngle * (possibleConnectingLines[t2].second.second/(possibleConnectingLines[t1].second.second + possibleConnectingLines[t2].second.second));
						if(connectionTerm1>connectionTerm2){  // remove t1
							possibleConnectingLines.erase(possibleConnectingLines.begin()+t1);
							t1--;
							break;
						}
						else{
							possibleConnectingLines.erase(possibleConnectingLines.begin()+t2);
							t2--;
						}
					}
					t2++;
				}
				t1++;
			}
		}

		for(unsigned int t2=0; t2<possibleConnectingLines.size(); t2++){
			unsigned int i = possibleConnectingLines[t2].first.first;
			unsigned int j = possibleConnectingLines[t2].first.second;
			// the 2 lines are within radius, they must be connected
			for(int t=line_discretization[j].size()-1; t>=0; t--){
				line_discretization[i].push_back(line_discretization[j][t]);
			}

			line_triangles[i] = copy_line_triangles[i];
			for(int t=line_triangles[j].size()-1; t>=0; t--){
				line_triangles[i].push_back(line_triangles[j][t]);
			}

			//line creation
			SurfaceSingularityLine* surf_line = m_graph.newSurfaceLine();
			int sepNumberTmp = m_graph.getNbLines();
			surf_line->setNumber(sepNumberTmp);

			//connect line to initial singularity point
			SingularityPoint* from_sing_pnt = searchSlots[i]->from_point;
			SingularityPoint::Slot* the_other_slot = searchSlots[j];
			SingularityPoint* towards_sing_pnt = the_other_slot->from_point;
			surf_line->addSingularityPoint(from_sing_pnt);
			surf_line->addSingularityPoint(towards_sing_pnt);
			//surf_line->addDiscretizationPoint(from_sing_pnt->getLocation());

			math::Vector3d firstDir = math::Vector3d(line_discretization[i][0], line_discretization[i][1]);
			searchSlots[i]->line_direction = firstDir;

			firstDir = math::Vector3d(line_discretization[j][0], line_discretization[j][1]);
			searchSlots[j]->line_direction = firstDir;

			streamlineDeviation[i] = streamlineDeviation[i] + streamlineDeviation[j];
			streamlineDeviation[i] = streamlineDeviation[i]/(line_discretization[i].size()+1);
			streamlineDeviation[j] = streamlineDeviation[i] ;
			//WARNING also compute streamlineDeviation locally, inside thresholdStreamLineDist
			searchSlots[i]->lineDeviation = streamlineDeviation[i];
			searchSlots[j]->lineDeviation = streamlineDeviation[i];

			for(unsigned int t=0; t<line_discretization[i].size(); t++) {
				surf_line->addDiscretizationPoint(line_discretization[i][t]);
			}

			for(unsigned int t=0; t<line_triangles[i].size(); t++) {
				surf_line->addTraversedFace(line_triangles[i][t]);
			}

			searchSlots[i]->line = surf_line;
			searchSlots[j]->line = surf_line;
			searchSlots[i]->isLaunched = true;
			searchSlots[i]->isFreeze = true;
			searchSlots[j]->isLaunched = true;
			searchSlots[j]->isFreeze = true;
			gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));

			for (unsigned int j = 0; j < line_discretization[i].size()-1; j++){
				gmds::Node n1 = m.newNode(line_discretization[i][j].X(), line_discretization[i][j].Y(), line_discretization[i][j].Z());
				gmds::Node n2 = m.newNode(line_discretization[i][j + 1].X(), line_discretization[i][j + 1].Y(), line_discretization[i][j + 1].Z());
				gmds::Face f  =  m.newTriangle(n1, n1, n2);
			}
			gmds::IGMeshIOService ioService(&m);
			gmds::VTKWriter vtkWriter(&ioService);
			vtkWriter.setCellOptions(gmds::N|gmds::F);
			vtkWriter.setDataOptions(gmds::N|gmds::F);
			std::string file_name = "SingToSing"+to_string(i)+"-"+to_string(j)+".vtk";
			vtkWriter.write(file_name);
			cout<<"has written SingToSing"<<file_name<<endl;

			noFrozenSLots = noFrozenSLots + 2;

		}

		for(unsigned int j0 = 0; j0 < line_discretization.size(); j0++){
		  	if(!searchSlots[j0]->isFreeze)
				cout<<"!searchSlots["<<j0<<"]->isFreeze"<<endl;
		}
		for(unsigned int j0 = 0; j0 < line_discretization.size(); j0++){
			gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
			cout<<"line_discretization["<<j0<<"].size() "<<line_discretization[j0].size()<<endl;
			for (unsigned int j = 0; j < line_discretization[j0].size()-1; j++){
				gmds::Node n1 = m.newNode(line_discretization[j0][j].X(), line_discretization[j0][j].Y(), line_discretization[j0][j].Z());
				gmds::Node n2 = m.newNode(line_discretization[j0][j + 1].X(), line_discretization[j0][j + 1].Y(), line_discretization[j0][j + 1].Z());
				gmds::Face f  =  m.newTriangle(n1, n1, n2);
			}

			gmds::IGMeshIOService ioService(&m);
			gmds::VTKWriter vtkWriter(&ioService);
			vtkWriter.setCellOptions(gmds::N|gmds::F);
			vtkWriter.setDataOptions(gmds::N|gmds::F);
			std::string file_name = "SimultaneousRK4_lastStep"+to_string(j0)+".vtk";
			vtkWriter.write(file_name);
		}

		writeOutputSingle("boundary_line");
		cout<<"wrote writeOutputSingle(\"boundary_line\"); "<<endl;
		currentSearchStep = currentSearchStep + searchStep;
		if(contor==m_mesh->getNbEdges()){
			cout<<"m_mesh->getNbEdges() "<<m_mesh->getNbEdges()<<endl;
			throw GMDSException("contor==total number of edges");
		}
	}

	/* remeshing
    	constructCenterTriangleCrosses();

   	vector<bool> visitedFaces(original_faces_number, false);

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
		  	gmds::Variable<gmds::math::Cross2D>* local_cross_field_2D = newLocalMesh.newVariable<math::Cross2D,GMDS_NODE>("cross_X");
			//at the begining all seem to be intialized as OX axis

            	remeshTriangles(&newLocalMesh, newLocalMesh_id_to_mesh_id_node,
					   local_cross_field_2D, trianglesToRemesh);

		  	//TODO also check exactly how the singularity slots have that direction (angle between last added 'vector' and the opoosite of slot direction)

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

        		//aiciaicithrow GMDSException("stop here");
		  	visitedFaces.clear();
        		visitedFaces.resize(original_faces_number, false);
        	}
   	}
*/

}

/*---------------------------------------------------------------------------*/
void SingularityGraphBuilder2D::growLine(SingularityPoint*               AFromSingPnt,
                                         SingularityPoint::Slot*         AFromSlot,
                                         SingularityPoint*&              AToSingPnt,
                                         SingularityPoint::Slot*&        AToSlot,
                                         gmds::math::Point&              AFromPnt,
                                         gmds::math::Vector3d&           AToDir,
                                         std::vector<gmds::math::Point>& APoints,
                                         std::vector<gmds::TCellID>&     ATriangles,
                                         int&                            AToCellDim,
                                         gmds::TCellID&                  AToCellID,
                                         double&                         streamlineDeviation,
                                         double&                         accumulatedDistancePerSlotLine,
                                         double&                         currentSearchStep,
                                         bool&                           AEndOnBnd,
                                         bool&                           AToSlotIsFree)
{

	bool find_end = false;
	//========================================================================
	/* grow a line starting from the singularity point AFromSingPnt through the slot AFromSlot
	as long as the geometric distance travelled by this line (accumulatedDistancePerSlotLine)
	is smaller than currentSearchStep or until we meet a boundary vertex */
	//========================================================================
	SingularityPoint*       found_pnt = 0;
	SingularityPoint::Slot* found_slot =0;

	math::Point  start_pnt = AFromPnt ; //starting point
	//APoints.push_back(start_pnt);
	math::Vector3d start_dir = AToDir; //starting direction
	math::Vector3d prev_dir  = AToDir; //prev direction used in the
	TCellID start_cell_id  = AToCellID ;
	int     start_cell_dim = AToCellDim;

	math::Point  current_pnt = start_pnt;
	math::Vector3d current_vec = start_dir;

	if(start_cell_dim==0)
		cout<<"node: "<<start_cell_id<<endl;
	else if(start_cell_dim==1){
		vector<Node> currentNodes = (m_mesh->get<Edge>(start_cell_id)).get<Node>();
	}

	math::Point start_dirPnt(start_pnt.X() + start_dir.X(),
	                         start_pnt.Y() + start_dir.Y(),
	                         start_pnt.Z() + start_dir.Z());

	// We check some termination conditions on the boundary.
	if (start_cell_dim==0){
		Node currentNode = m_mesh->get<Node>(start_cell_id);
		if(m_mesh->isMarked(currentNode, m_mark_nodes_on_point) ||
		m_mesh->isMarked(currentNode, m_mark_nodes_on_curve)){
			AEndOnBnd       = true;
		}
	}
	else { //we have necessarry start_cell_dim=1
		Edge currentEdge = m_mesh->get<Edge>(start_cell_id);
		if(m_mesh->isMarked(currentEdge, m_mark_edges_on_curve)){
			AEndOnBnd       = true;
		}
	}

	while((accumulatedDistancePerSlotLine<currentSearchStep)&&(!AEndOnBnd)&&(!find_end)){
		TCellID next_cell_id  = NullID;
		int     next_cell_dim = -1;
		m_tool.findNextCell(start_pnt, start_dir,
		start_cell_dim, start_cell_id,
		next_cell_dim, next_cell_id);

		if (next_cell_dim == -1){
			// The cell defined by (start_cell_dim, start_cell_id) is on the boundary.
			AEndOnBnd = true;
		}
		else if (next_cell_dim == 1){
			cout<<"next_cell_dim == 1"<<endl;
			// we are going along an edge.
			// Our simple assumption is to follow this edge until reaching
			// one of its end points and to compute the next direction at
			// this point.
			Edge currentEdge = m_mesh->get<Edge>(next_cell_id);
			std::vector<TCellID> adj_faces = currentEdge.getIDs<Face>();
			ATriangles.insert(ATriangles.end(),adj_faces.begin(),adj_faces.end());

			std::vector<Node> currentNodes = currentEdge.get<Node>();
			math::Vector3d v0(start_pnt, currentNodes[0].getPoint());
			math::Vector3d v1(start_pnt, currentNodes[1].getPoint());
			Node next_node;
			if(math::near(v0.norm(),0.0))
				next_node = currentNodes[1];
			else if(math::near(v1.norm(),0.0))
				next_node = currentNodes[0];
			else if(v0.dot(start_dir) > v1.dot(start_dir))
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
			accumulatedDistancePerSlotLine = accumulatedDistancePerSlotLine + start_pnt.distance(APoints[APoints.size()-1]);
			APoints.push_back(start_pnt);
		}
		else {
			//general case, we are in a face
			Face currentFace = m_mesh->get<Face>(next_cell_id);
			ATriangles.push_back(currentFace.id());
			//==============================================================
			// CASE 1: ARE WE IN A FACE CONTAINING A SING. POINT?
			//VERY IMPROBABLE; GOT HERE BUT WE WEREN'T WITHIN THRESHOLD
			//==============================================================
			bool intersect_sing_point = false;

			SingularityPoint* next_sing_point = m_faces_to_singularity_on_surf[currentFace.id()];
			bool must_try_to_connect = false;
			if (next_sing_point != NULL && //face in a confusing ball ...
				next_sing_point != AFromSingPnt) {// ... of another singularity point
				must_try_to_connect = true;
			}
			else{
				if (next_sing_point != NULL && //face in the confusing ball ...
				next_sing_point == AFromSingPnt) {//... of the incoming singularity point
					cout<<"face in a confusing ball of the incoming singularity point"<<endl;
					//Warning: completly empiric, we just try to detect cyclic lines
					if (APoints.size() >= 100)
						must_try_to_connect = true;
				}
			}

			if(must_try_to_connect) {
				cout<<"must_try_to_connect"<<endl;
				math::Point start_dirPnt(start_pnt.X() + start_dir.X(),
				                         start_pnt.Y() + start_dir.Y(),
				                         start_pnt.Z() + start_dir.Z());

				//Now, we look for a compatible slot
				std::vector<SingularityPoint::Slot*>& cur_slots = next_sing_point->getSlots();
				//==============================================================
				// We look for a free slot and we connect the line to it
				//==============================================================
				bool found_compatible_slot = false;
				bool found_free_slot = false;
				found_pnt = next_sing_point;
				double slot_epsilon = 0.9;
				//we normalize the vector we arrive with
				math::Vector3d current_vec = start_dir;
				current_vec.normalize();
				while (!found_compatible_slot && slot_epsilon > 0.4) {
					double best_deviation = -2;
					double best_slot_id = 0;

					for (unsigned int i_slot = 0; i_slot < cur_slots.size(); i_slot++) {
						SingularityPoint::Slot* current_slot = cur_slots[i_slot];
						if(current_slot->isFreeze)
							continue;
						math::Vector3d slot_opp_dir = current_slot->direction.opp();
						slot_opp_dir.normalize();
						double slot_deviation =slot_opp_dir.dot(current_vec);

						if (slot_deviation > slot_epsilon && slot_deviation > best_deviation) {
							best_deviation = slot_deviation;
							best_slot_id=i_slot;
						}
					}
					if(best_deviation!=-2 &&  !cur_slots.empty()){
						SingularityPoint::Slot* best_slot = cur_slots[best_slot_id];
						math::Vector3d            slot_opp_dir = best_slot->direction.opp();
						slot_opp_dir.normalize();
						if ( best_slot->isLaunched) {
						 //slot already assigned with a previous line (and direction)
							math::Vector3d prev_line_dir = best_slot->line_direction;
							prev_line_dir.normalize();
							double prev_deviation =slot_opp_dir.dot(prev_line_dir);

							if(best_deviation>prev_deviation) { //the new alignment is better than the previous one
								found_free_slot = false;
								found_compatible_slot = true;
								found_slot = best_slot;
								intersect_sing_point = true;
							}
							else {// We keep the previous association
								found_compatible_slot = false;
							}
						}
						else {  //WE HAVE A FREE SLOT
							// We keep the previous association
							found_free_slot = true;
							AToSlotIsFree = found_free_slot;
							found_compatible_slot = true;
							found_slot = best_slot;
							intersect_sing_point = true;
						}
						// HAVE WE FIND THE END OF THE LINE??
						if(found_compatible_slot) {
							find_end = true;
						}
					}
					slot_epsilon -=0.1;
				}
			}

			//==============================================================
			// CASE 2: SO, WE JUST HAVE TO CROSS THE TRIANGLE (GENERAL CONF)
			//==============================================================
			//Does the current triangle has the same classif
			if ( !intersect_sing_point) {
				math::Point  out_pnt;
				math::Vector3d out_vec;
				TCellID out_cell_id;
				int out_cell_dim;

				m_tool.traverseTriangle(currentFace,   /* the face we work on*/
				                        start_pnt,      /* the geometric point we start from */
				                        start_dir,      /* the geometric direction to follow*/
				                        start_cell_dim, /* the dimension of the cell start_pnt is located */
				                        start_cell_id,  /* the id of the cell start_pnt is located on*/
				                        out_pnt,        /* the geometric point where we go out */
				                        out_vec,        /* the geometric direction to follow after*/
				                        out_cell_dim,   /* the dimension of the out cell (0 or 1) */
				                        out_cell_id,    /* the id of the out cell*/
				                        streamlineDeviation);   /*deviation of the streamline up to this point*/

				accumulatedDistancePerSlotLine = accumulatedDistancePerSlotLine + out_pnt.distance(APoints[APoints.size()-1]);
				APoints.push_back(out_pnt);

				// we progress to the next point, next vector and so next face too
				prev_dir = start_dir; //we store the prev direction for slot
				//reconnection with balls
				start_pnt = out_pnt;
				start_dir = out_vec;
				start_cell_dim = out_cell_dim;
				start_cell_id = out_cell_id;
			}
			if (intersect_sing_point){
				find_end = true;
			}
			//post process, we just check if we haven't arrived onto a geometric boundary
			if (start_cell_dim==0){
				Node currentNode = m_mesh->get<Node>(start_cell_id);
				if(m_mesh->isMarked(currentNode, m_mark_nodes_on_point) ||
				m_mesh->isMarked(currentNode, m_mark_nodes_on_curve)){
					find_end = true;
					AEndOnBnd = true;
				}
			}
			else { //we have necessarry start_cell_dim=1

				Edge currentEdge = m_mesh->get<Edge>(start_cell_id);
				if(m_mesh->isMarked(currentEdge, m_mark_edges_on_curve)){
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

void SingularityGraphBuilder2D::growLineRK4(SingularityPoint*               AFromSingPnt,
                                            SingularityPoint::Slot*         AFromSlot,
                                            SingularityPoint*&              AToSingPnt,
                                            SingularityPoint::Slot*&        AToSlot,
                                            gmds::math::Point&              AFromPnt,
                                            gmds::math::Vector3d&           AToDir,
                                            std::vector<gmds::math::Point>& APoints,
                                            std::vector<gmds::TCellID>&     ATriangles,
                                            int&                            AToCellDim,
                                            gmds::TCellID&                  AToCellID,
                                            double&                         streamlineDeviation,
                                            double&                         stepSize,
                                            bool&                           AEndOnBnd,
                                            bool&                           AToSlotIsFree,
                                            bool&                           must_try_to_connect)
{

	AEndOnBnd = false;
	//========================================================================
	/* grow a line starting from the singularity point AFromSingPnt through the slot AFromSlot
	 as long as the geometric distance travelled by this line
	 is smaller than stepSize or until we reach the boundary (or another slot)
	*/
	//========================================================================
	SingularityPoint*       found_pnt = 0;
	SingularityPoint::Slot* found_slot =0;

	math::Point  start_pnt = AFromPnt ; //starting point
	//APoints.push_back(start_pnt);
	math::Vector3d start_dir = AToDir; //starting direction
	math::Vector3d prev_dir  = AToDir; //prev direction used in the

	TCellID start_cell_id  = AToCellID ;
	int     start_cell_dim = AToCellDim;

	vector<gmds::Face> AFaces;

	if(AToCellDim==0){
		AFaces = m_mesh->get<gmds::Node>(AToCellID).get<gmds::Face>();
		if(m_mesh->isMarked(m_mesh->get<gmds::Node>(AToCellID),m_mark_nodes_on_curve) ||
		m_mesh->isMarked(m_mesh->get<gmds::Node>(AToCellID),m_mark_nodes_on_point)){

			gmds::math::Point point_1 = start_pnt + start_dir*2*math::Constants::EPSILON;
			bool insideTri = false;
			for(unsigned int i=0; i<AFaces.size(); i++){
				std::vector<Node> f_nodes = AFaces[i].get<Node>();
				gmds::math::Triangle ATriangle(f_nodes[0].getPoint(), f_nodes[1].getPoint(), f_nodes[2].getPoint());

				if(ATriangle.isIn2ndMethod(point_1)){
					insideTri = true;
					break;
				}
			}
			if(!insideTri){
				AEndOnBnd = true;
			//APoints.push_back(start_pnt);
			}
		}
	}
	else{
		if(AToCellDim==1){
			AFaces = m_mesh->get<gmds::Edge>(AToCellID).get<gmds::Face>();
			if(m_mesh->isMarked(m_mesh->get<gmds::Edge>(AToCellID),m_mark_edges_on_curve)){
				gmds::math::Point point_1 = start_pnt + start_dir*2*math::Constants::EPSILON;
				bool insideTri = false;
				for(unsigned int i=0; i<AFaces.size(); i++){
					std::vector<Node> f_nodes = AFaces[i].get<Node>();
					gmds::math::Triangle ATriangle(f_nodes[0].getPoint(), f_nodes[1].getPoint(), f_nodes[2].getPoint());

					if(ATriangle.isIn2ndMethod(point_1)){
						insideTri = true;
						break;
					}
				}
				if(!insideTri){
					AEndOnBnd = true;
					//APoints.push_back(start_pnt);
				}
			}
		}
		else
			AFaces.push_back(m_mesh->get<gmds::Face>(AToCellID));
	}

	if(!AEndOnBnd){
		//==============================================================
		// CASE 1: ARE WE IN A FACE CONTAINING A SING. POINT?
		//VERY IMPROBABLE; GOT HERE BUT WE WEREN'T WITHIN THRESHOLD
		//==============================================================

		bool intersect_sing_point = false;
		must_try_to_connect = false;
		SingularityPoint* next_sing_point;

		for(unsigned int i=0; i<AFaces.size(); i++){
			next_sing_point = m_faces_to_singularity_on_surf[AFaces[i].id()];
			if(next_sing_point != NULL && //face in a confusing ball ...
			next_sing_point != AFromSingPnt) {// ... of another singularity point
				must_try_to_connect = true;
				break;
			}
			else{
				if(next_sing_point != NULL && //face in the confusing ball ...
				next_sing_point == AFromSingPnt) {//... of the incoming singularity point

					//Warning: completly empiric, we just try to detect cyclic lines
					if (APoints.size() >= 100)
						must_try_to_connect = true;
				}
			}
		}

		if(must_try_to_connect) {
			math::Point start_dirPnt(start_pnt.X() + start_dir.X(),
			                         start_pnt.Y() + start_dir.Y(),
			                         start_pnt.Z() + start_dir.Z());

			//Now, we look for a compatible slot
			vector<SingularityPoint::Slot*>& cur_slots = next_sing_point->getSlots();

			//==============================================================
			// We look for a free slot and we connect the line to it
			//==============================================================
			bool found_compatible_slot = false;
			bool found_free_slot = false;
			found_pnt = next_sing_point;
			double slot_epsilon = 0.9;
			//we normalize the vector we arrive with
			math::Vector3d current_vec = start_dir;
			current_vec.normalize();
			while (!found_compatible_slot && slot_epsilon > 0.4) {
				double best_deviation = -2;
				double best_slot_id = 0;
				for(unsigned int i_slot = 0; i_slot < cur_slots.size(); i_slot++) {
					//current_vec = start_dir;  						//current_vec.normalize();
					SingularityPoint::Slot* current_slot = cur_slots[i_slot];
					if(current_slot->isFreeze)
						continue;

					math::Vector3d resultingVec(start_dirPnt, current_slot->location);

					if(start_dir.angle(resultingVec)>gmds::math::Constants::PIDIV2)
						continue;

					math::Vector3d slot_opp_dir = current_slot->direction.opp();
					slot_opp_dir.normalize();
					//current_vec = (current_slot->location) - start_pnt;
					double slot_deviation = slot_opp_dir.dot(current_vec);

					if (slot_deviation > slot_epsilon && slot_deviation > best_deviation) {
						best_deviation = slot_deviation;
						best_slot_id=i_slot;
					}
				}

				if(best_deviation!=-2 &&  !cur_slots.empty()){
					SingularityPoint::Slot* best_slot = cur_slots[best_slot_id];
					math::Vector3d            slot_opp_dir = best_slot->direction.opp();
					slot_opp_dir.normalize();
					if ( best_slot->isLaunched) {
						//slot already assigned with a previous line (and direction)
						math::Vector3d prev_line_dir = best_slot->line_direction;
						prev_line_dir.normalize();
						double prev_deviation =slot_opp_dir.dot(prev_line_dir);
						if(best_deviation>prev_deviation) {
							//the new alignment is better than the previous one
							found_free_slot = false;
							found_compatible_slot = true;
							found_slot = best_slot;
							intersect_sing_point = true;
						}
						else{
							// We keep the previous association
							found_compatible_slot = false;
							//must_try_to_connect = false;
						}
					}
					else{  //WE HAVE A FREE SLOT
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
				slot_epsilon -=0.1;
			}
		}
		//==============================================================
		// CASE 2: SO, WE JUST HAVE TO CROSS THE TRIANGLE (GENERAL CONF)
		//==============================================================

		if(!intersect_sing_point) {
			math::Point  out_pnt;
			math::Vector3d out_vec;

			m_tool.RK4Computation(start_pnt,
			                      start_dir,
			                      out_pnt,
			                      out_vec,
			                      streamlineDeviation,
			                      stepSize,
			                      AEndOnBnd,
			                      ATriangles,
			                      start_cell_dim,
			                      start_cell_id);

			// accumulatedDistancePerSlotLine = accumulatedDistancePerSlotLine + out_pnt.distance(APoints[APoints.size()-1]);
			APoints.push_back(out_pnt);
			//if(!AEndOnBnd)
			//  start_cell_id = ATriangles[ATriangles.size()-1];
			// we progress to the next point, next vector and so next face too
			prev_dir = start_dir; //we store the prev direction for slot
			start_pnt = out_pnt;
			start_dir = out_vec;
		}
		AToCellDim = start_cell_dim;
		AToCellID = start_cell_id;
	}
	//==============================================================
	// Update of out parameters
	//==============================================================
	//last followed direction
	AToDir = start_dir;
	//singularity point data if we found an end point
	AToSingPnt = found_pnt;
	AToSlot = found_slot;

}

/*---------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::writeTestMeshVerts(vector<gmds::Node>& ANodes,
                                                   std::string&        AFileName)
{
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
	math::Point APoint = ANodes[0].getPoint();
	for (unsigned int j = 1; j < ANodes.size(); j++){
		math::Point prevPoint = APoint;
		APoint = ANodes[j].getPoint();
		gmds::Node n1 = m.newNode(prevPoint.X(), prevPoint.Y(), prevPoint.Z());
		gmds::Node n2 = m.newNode(APoint.X(), APoint.Y(), APoint.Z());
		gmds::Face f  =  m.newTriangle(n1, n1, n2);
	}
	gmds::IGMeshIOService ioService(&m);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write(AFileName);
}

/*---------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::writeTestPoints(vector<gmds::math::Point>& APoints,
                                                std::string&               AFileName)
{
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
	math::Point APoint = APoints[0];
	for (unsigned int j = 1; j < APoints.size(); j++){
		math::Point prevPoint = APoint;
		APoint = APoints[j];
		gmds::Node n1 = m.newNode(prevPoint.X(), prevPoint.Y(), prevPoint.Z());
		gmds::Node n2 = m.newNode(APoint.X(), APoint.Y(), APoint.Z());
		gmds::Face f  =  m.newTriangle(n1, n1, n2);
	}
	gmds::IGMeshIOService ioService(&m);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write(AFileName);
}
/*---------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::writeTestMeshTriangles(vector<gmds::Face>& ATriangles,
                                                       std::string&        AFileName)
{
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));

	for (unsigned int j = 0; j < ATriangles.size(); j++){
		vector<gmds::Node> currentNodes = ATriangles[j].get<Node>();
		math::Point current_point1 = currentNodes[0].getPoint();
		math::Point current_point2 = currentNodes[1].getPoint();
		math::Point current_point3 = currentNodes[2].getPoint();
		gmds::Node n1 = m.newNode(current_point1.X(), current_point1.Y(), current_point1.Z());
		gmds::Node n2 = m.newNode(current_point2.X(), current_point2.Y(), current_point2.Z());
		gmds::Node n3 = m.newNode(current_point3.X(), current_point3.Y(), current_point3.Z());
		gmds::Face f  =  m.newTriangle(n1, n2, n3);
	}
	gmds::IGMeshIOService ioService(&m);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write(AFileName);
}

/*---------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::writeTestMeshTrianglesIds(vector<gmds::TCellID>& ATrianglesIds,
                                                          std::string&           AFileName)
{
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));

	for (unsigned int j = 0; j < ATrianglesIds.size(); j++){
		vector<gmds::Node> currentNodes = (m_mesh->get<gmds::Face>(ATrianglesIds[j])).get<Node>();
		math::Point current_point1 = currentNodes[0].getPoint();
		math::Point current_point2 = currentNodes[1].getPoint();
		math::Point current_point3 = currentNodes[2].getPoint();
		gmds::Node n1 = m.newNode(current_point1.X(), current_point1.Y(), current_point1.Z());
		gmds::Node n2 = m.newNode(current_point2.X(), current_point2.Y(), current_point2.Z());
		gmds::Node n3 = m.newNode(current_point3.X(), current_point3.Y(), current_point3.Z());
		gmds::Face f  =  m.newTriangle(n1, n2, n3);
	}
	gmds::IGMeshIOService ioService(&m);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write(AFileName);
}

/*---------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::remeshTriangles(gmds::Mesh*                          newLocalMesh,
                                                vector<gmds::TCellID>&               newLocalMesh_id_to_mesh_id_node,
                                                gmds::Variable<gmds::math::Cross2D>* local_cross_field_2D,
                                                vector<gmds::TCellID>&               trianglesToRemesh)
{
	/* remesh the triangles in trianglesToRemesh; assumes that the original mesh is consistently oriented*/

	vector<gmds::TCellID> mesh_id_to_newLocalMesh_id_node(original_nodes_number);
	vector<bool> visitedVerts(original_nodes_number, false);

	math::Point current_point1, current_point2, current_point3, current_point_center;
	unsigned int NodeNumber = 0;
	for (unsigned int i = 0; i < trianglesToRemesh.size(); i++){
		vector<gmds::Node> currentNodes = (m_mesh->get<gmds::Face>(trianglesToRemesh[i])).get<Node>();
		current_point1 = currentNodes[0].getPoint();
		current_point2 = currentNodes[1].getPoint();
		current_point3 = currentNodes[2].getPoint();

		gmds::Node n1;
		if(!visitedVerts[currentNodes[0].id()]){
			n1 = newLocalMesh->newNode(current_point1.X(), current_point1.Y(), current_point1.Z());
			visitedVerts[currentNodes[0].id()] = true;
			mesh_id_to_newLocalMesh_id_node[currentNodes[0].id()] = NodeNumber;
			newLocalMesh_id_to_mesh_id_node.push_back(currentNodes[0].id());
			NodeNumber++;
		}
		else{
			n1 = newLocalMesh->get<gmds::Node>(mesh_id_to_newLocalMesh_id_node[currentNodes[0].id()]);
		}
		gmds::Node n2;
		if(!visitedVerts[currentNodes[1].id()]){
			n2 = newLocalMesh->newNode(current_point2.X(), current_point2.Y(), current_point2.Z());
			visitedVerts[currentNodes[1].id()] = true;
			mesh_id_to_newLocalMesh_id_node[currentNodes[1].id()] = NodeNumber;
			newLocalMesh_id_to_mesh_id_node.push_back(currentNodes[1].id());
			NodeNumber++;
		}
		else{
			n2 = newLocalMesh->get<gmds::Node>(mesh_id_to_newLocalMesh_id_node[currentNodes[1].id()]);
		}
		gmds::Node n3;
		if(!visitedVerts[currentNodes[2].id()]){
			n3 = newLocalMesh->newNode(current_point3.X(), current_point3.Y(), current_point3.Z());
			visitedVerts[currentNodes[2].id()] = true;
			mesh_id_to_newLocalMesh_id_node[currentNodes[2].id()] = NodeNumber;
			newLocalMesh_id_to_mesh_id_node.push_back(currentNodes[2].id());
			NodeNumber++;
		}
		else{
			n3 = newLocalMesh->get<gmds::Node>(mesh_id_to_newLocalMesh_id_node[currentNodes[2].id()]);
		}

		gmds::Node newn = newLocalMesh->newNode(triangle_centers[trianglesToRemesh[i]].X(), triangle_centers[trianglesToRemesh[i]].Y(), triangle_centers[trianglesToRemesh[i]].Z());
		newLocalMesh_id_to_mesh_id_node.push_back(original_nodes_number+1);
		NodeNumber++;
		gmds::Face  f1 = newLocalMesh->newTriangle(n1, n2, newn);
		gmds::Face  f2 = newLocalMesh->newTriangle(n2, n3, newn);
		gmds::Face  f3 = newLocalMesh->newTriangle(n3, n1, newn);
	}

	gmds::IGMeshIOService ioService(newLocalMesh);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("localmesh.vtk");

	gmds::MeshDoctor doc(newLocalMesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	vector<gmds::Node> m_curve_nodes, m_surf_nodes;
	visitedVerts.clear();
	visitedVerts.resize(newLocalMesh->getNbNodes(), false);

	math::Vector3d OX(1,0,0);
	for(auto e_id:newLocalMesh->edges()){
		gmds::Edge currentEdge = newLocalMesh->get<gmds::Edge>(e_id);
		vector<gmds::Node> currentNodes = (newLocalMesh->get<gmds::Edge>(e_id)).get<gmds::Node>();
		std::vector<gmds::Face> adj_faces = currentEdge.get<gmds::Face>();
		if(adj_faces.size()==1){
			if(!visitedVerts[currentNodes[0].id()]){
				m_curve_nodes.push_back(currentNodes[0]);
				(*local_cross_field_2D)[currentNodes[0].id()] = (*m_field)[newLocalMesh_id_to_mesh_id_node[currentNodes[0].id()]];
				visitedVerts[currentNodes[0].id()] = true;
			}

			if(!visitedVerts[currentNodes[1].id()]){
				m_curve_nodes.push_back(currentNodes[1]);
				(*local_cross_field_2D)[currentNodes[1].id()] = (*m_field)[newLocalMesh_id_to_mesh_id_node[currentNodes[1].id()]];
				visitedVerts[currentNodes[0].id()] = true;
			}
		}
	}

	for(auto n_id:newLocalMesh->nodes()){
		if(!visitedVerts[n_id])
			m_surf_nodes.push_back(newLocalMesh->get<gmds::Node>(n_id));
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


	gmds::IGMeshIOService meshIoServref(newLocalMesh);
	gmds::VTKWriter writerref(&meshIoServref);
	writerref.setCellOptions(gmds::N|gmds::F);
	writerref.setDataOptions(gmds::N|gmds::F);

	std::stringstream file_nameref;
	file_nameref <<"before-newLocalMesh-crossVectors.vtk";
    	writerref.write(file_nameref.str());
		*/
	std::vector<Face> mesh_faces;
	mesh_faces.resize(newLocalMesh->getNbFaces());
	int f_index=0;
	for(auto f_id:newLocalMesh->faces()){
		mesh_faces[f_index++]= newLocalMesh->get<Face>(f_id);
	}

	LaplaceCross2D algo(newLocalMesh,
	                    local_cross_field_2D,
	                    m_curve_nodes,
	                    m_surf_nodes,
	                    mesh_faces);
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


	gmds::IGMeshIOService meshIoServref(newLocalMesh);
	gmds::VTKWriter writerref(&meshIoServref);
	writerref.setCellOptions(gmds::N|gmds::F);
	writerref.setDataOptions(gmds::N|gmds::F);

	std::stringstream file_nameref;
	file_nameref <<"final-newLocalMesh-crossVectors.vtk";
    	writerref.write(file_nameref.str()); */
}

/*---------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::remeshTrianglesNewMesh(gmds::Mesh*                          newLocalMesh,
                                                       gmds::Variable<gmds::math::Cross2D>* newMesh_cross_field_2D,
                                                       vector<gmds::TCellID>&               trianglesToRemesh,
                                                       vector<bool>&                        trianglesToRemeshBool,
                                                       vector<gmds::math::Point>&           newTriangleCenters,
                                                       vector<gmds::math::Cross2D>&         newTriangleCenterCrosses,
                                                       vector<gmds::math::Vector3d>&        newBdryEdgeNormals,
                                                       vector<gmds::math::Vector3d>&        newBdryNodeNormals,
                                                       vector<bool>&                        isCurveEdge,
                                                       vector<bool>&                        isCurveNode,
                                                       vector<gmds::TCellID>&               newLocalMesh_id_to_mesh_id_node,
                                                       vector<gmds::TCellID>&               mesh_id_to_newLocalMesh_id_node,
                                                       vector<gmds::TCellID>&               newLocalMesh_id_to_mesh_id_face,
                                                       vector<gmds::TCellID>&               mesh_id_to_newLocalMesh_id_face)
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
	cout<<"remeshTrianglesNewMesh"<<endl;

	vector<bool> visitedVerts(original_nodes_number, false);

	math::Point current_point1, current_point2, current_point3, current_point_center;
	unsigned int NodeNumber = 0;
	unsigned int FaceNumber = 0;
	unsigned int EdgeNumber = 0;
	vector<gmds::TCellID> origEdge2NewPoint(m_mesh->getNbEdges());
	vector<bool> modifiedEdge(m_mesh->getNbEdges(), false);
	vector<gmds::Node> NodesToAdd(6);
	vector<TCoord> edgeWeights(2, 0.5);
	for(unsigned int i=0; i<trianglesToRemesh.size(); i++){
		gmds::Face currentFace = m_mesh->get<gmds::Face>(trianglesToRemesh[i]);
		vector<gmds::Node> currentNodes = currentFace.get<Node>();
		current_point1 = currentNodes[0].getPoint();
		current_point2 = currentNodes[1].getPoint();
		current_point3 = currentNodes[2].getPoint();

		if(!visitedVerts[currentNodes[0].id()]){
			NodesToAdd[0] = newLocalMesh->newNode(current_point1.X(), current_point1.Y(), current_point1.Z());
			visitedVerts[currentNodes[0].id()] = true;
			mesh_id_to_newLocalMesh_id_node[currentNodes[0].id()] = NodeNumber;
			newLocalMesh_id_to_mesh_id_node.push_back(currentNodes[0].id());
			(*newMesh_cross_field_2D)[NodeNumber] = (*m_field)[currentNodes[0].id()];
			NodeNumber++;
		}
		else{
			NodesToAdd[0] = newLocalMesh->get<gmds::Node>(mesh_id_to_newLocalMesh_id_node[currentNodes[0].id()]);
		}

		if(!visitedVerts[currentNodes[1].id()]){
			NodesToAdd[1] = newLocalMesh->newNode(current_point2.X(), current_point2.Y(), current_point2.Z());
			visitedVerts[currentNodes[1].id()] = true;
			mesh_id_to_newLocalMesh_id_node[currentNodes[1].id()] = NodeNumber;
			newLocalMesh_id_to_mesh_id_node.push_back(currentNodes[1].id());
			(*newMesh_cross_field_2D)[NodeNumber] = (*m_field)[currentNodes[1].id()];
			NodeNumber++;
		}
		else{
			NodesToAdd[1] = newLocalMesh->get<gmds::Node>(mesh_id_to_newLocalMesh_id_node[currentNodes[1].id()]);
		}

		if(!visitedVerts[currentNodes[2].id()]){
			NodesToAdd[2] = newLocalMesh->newNode(current_point3.X(), current_point3.Y(), current_point3.Z());
			visitedVerts[currentNodes[2].id()] = true;
			mesh_id_to_newLocalMesh_id_node[currentNodes[2].id()] = NodeNumber;
			newLocalMesh_id_to_mesh_id_node.push_back(currentNodes[2].id());
			(*newMesh_cross_field_2D)[NodeNumber] = (*m_field)[currentNodes[2].id()];
			NodeNumber++;
		}
		else{
			NodesToAdd[2] = newLocalMesh->get<gmds::Node>(mesh_id_to_newLocalMesh_id_node[currentNodes[2].id()]);
		}

		vector<gmds::Edge> currentEdges = currentFace.get<gmds::Edge>();
		for(unsigned int k=0; k<3; k++){
			vector<gmds::Node> edge_nodes = currentEdges[k].get<gmds::Node>();
			for(unsigned int j=0; j<3; j++){// get opposite node for each edge; put in this order
				if((edge_nodes[0]!=NodesToAdd[j])&&(edge_nodes[1]!=NodesToAdd[j])){
					if(!modifiedEdge[currentEdges[k].id()]){
						gmds::math::Point middle_point = currentEdges[k].center();
 						NodesToAdd[3+fmod(j,3)] = newLocalMesh->newNode(middle_point.X(), middle_point.Y(), middle_point.Z());

						vector<math::Vector3d> c_vectors0 = ((*m_field)[edge_nodes[0].id()]).componentVectors();
						vector<math::Vector3d> c_vectors1 = ((*m_field)[edge_nodes[1].id()]).componentVectors();

						gmds::math::Vector3d second_closest_vect =  (*m_field)[edge_nodes[1].id()].closestComponentVector(c_vectors0[0]);
						gmds::math::Vector3d center_vect = c_vectors0[0] * edgeWeights[0] + second_closest_vect * edgeWeights[1];
						gmds::math::Vector3d center_vect_second = center_vect.getOneOrtho();

						(*newMesh_cross_field_2D)[NodeNumber] = gmds::math::Cross2D(center_vect, center_vect_second);
						origEdge2NewPoint[currentEdges[k].id()] = EdgeNumber;
						modifiedEdge[currentEdges[k].id()] = true;
						EdgeNumber++;
						NodeNumber++;
					}
					else{
						NodesToAdd[3+fmod(j,3)] = newLocalMesh->get<gmds::Node>(origEdge2NewPoint[currentEdges[k].id()]);
					}

				}
			}
		}
		gmds::Face  f1 = newLocalMesh->newTriangle(NodesToAdd[3], NodesToAdd[4], NodesToAdd[5]);
		gmds::Face  f2 = newLocalMesh->newTriangle(NodesToAdd[0], NodesToAdd[5], NodesToAdd[4]);
		gmds::Face  f3 = newLocalMesh->newTriangle(NodesToAdd[1], NodesToAdd[3], NodesToAdd[5]);
		gmds::Face  f4 = newLocalMesh->newTriangle(NodesToAdd[2], NodesToAdd[4], NodesToAdd[3]);

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
	vector<bool> remeshedAdjFaces(original_faces_number, false);

	for(auto e_id:m_mesh->edges()){
		gmds::Edge currentEdge = m_mesh->get<gmds::Edge>(e_id);

		if(modifiedEdge[e_id]){
			vector<gmds::Face> adj_faces = currentEdge.get<gmds::Face>();
			for(unsigned int i=0; i<adj_faces.size(); i++){
				if(!trianglesToRemeshBool[adj_faces[i].id()]){
					vector<gmds::Node> edge_nodes = currentEdge.get<gmds::Node>();
					gmds::Node midNode = newLocalMesh->get<gmds::Node>(origEdge2NewPoint[currentEdge.id()]);
					if(remeshedAdjFaces[adj_faces[i].id()]){
					  // TODO

					}
					else{
						vector<gmds::Node> face_nodes = adj_faces[i].get<gmds::Node>();

						for(unsigned int j=0; j<3; j++){
							if((face_nodes[j]!=edge_nodes[0])&&(face_nodes[j]!=edge_nodes[1])){
								gmds::Node opposite_node = face_nodes[j];
								gmds::Node opposite_node_local;//TODO = does it exist? or to add?
								//same for face_nodes!!!TODO
								gmds::Face  f1 = newLocalMesh->newTriangle(opposite_node_local, face_nodes[fmod((j+1),3)], midNode);
								gmds::Face  f2 = newLocalMesh->newTriangle(opposite_node_local, midNode, face_nodes[fmod((j+2),3)]);
								break;
							}
						}

						remeshedAdjFaces[adj_faces[i].id()] = true;
					}
				}
			}
		}
	}
	gmds::IGMeshIOService ioService(newLocalMesh);
 	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("localmesh.vtk");

	gmds::MeshDoctor doc(newLocalMesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	isCurveEdge.resize(newLocalMesh->getNbEdges(), false);
	isCurveNode.resize(newLocalMesh->getNbNodes(), false);
	vector<gmds::Node> m_curve_nodes, m_surf_nodes;
	visitedVerts.clear();
	visitedVerts.resize(newLocalMesh->getNbNodes(), false);

	math::Vector3d OX(1,0,0);
	for(auto e_id:newLocalMesh->edges()){
		gmds::Edge currentEdge = newLocalMesh->get<gmds::Edge>(e_id);
		vector<gmds::Node> currentNodes = (newLocalMesh->get<gmds::Edge>(e_id)).get<gmds::Node>();
		std::vector<gmds::Face> adj_faces = currentEdge.get<gmds::Face>();
		if(adj_faces.size()==1){
			isCurveEdge[e_id] = true;
			if(!visitedVerts[currentNodes[0].id()]){
				isCurveNode[currentNodes[0].id()] = true;
				m_curve_nodes.push_back(currentNodes[0]);
				(*newMesh_cross_field_2D)[currentNodes[0].id()] = (*m_field)[newLocalMesh_id_to_mesh_id_node[currentNodes[0].id()]];
				visitedVerts[currentNodes[0].id()] = true;
			}

			if(!visitedVerts[currentNodes[1].id()]){
				isCurveNode[currentNodes[1].id()] = true;
				m_curve_nodes.push_back(currentNodes[1]);
				(*newMesh_cross_field_2D)[currentNodes[1].id()] = (*m_field)[newLocalMesh_id_to_mesh_id_node[currentNodes[1].id()]];
				visitedVerts[currentNodes[0].id()] = true;
			}
		}
	}

	//unfortunately we recompute the frame algorithm for the entire new mesh

	newTriangleCenterCrosses.resize(newLocalMesh->getNbFaces());
	vector<TCoord> AWeights(3, (double)1/3);
	for(auto f_id:newLocalMesh->faces()){
		std::vector<math::Cross2D> Tricrosses;
		gmds::Face currentFace = newLocalMesh->get<Face>(f_id);
		vector<Node> nodesTri = currentFace.get<Node>();
		Tricrosses.push_back((*newMesh_cross_field_2D)[nodesTri[0].id()]);
		Tricrosses.push_back((*newMesh_cross_field_2D)[nodesTri[1].id()]);
		Tricrosses.push_back((*newMesh_cross_field_2D)[nodesTri[2].id()]);
		vector<math::Vector3d> c_vectors0 = Tricrosses[0].componentVectors();
		vector<math::Vector3d> c_vectors1 = Tricrosses[1].componentVectors();
		vector<math::Vector3d> c_vectors2 = Tricrosses[2].componentVectors();

		gmds::math::Vector3d second_closest_vect =  Tricrosses[1].closestComponentVector(c_vectors0[0]);
		gmds::math::Vector3d third_closest_vect =  Tricrosses[2].closestComponentVector(c_vectors0[0]);
		gmds::math::Vector3d center_vect = c_vectors0[0] * AWeights[0] + second_closest_vect * AWeights[1] + third_closest_vect * AWeights[2];
		gmds::math::Vector3d center_vect_second = center_vect.getOneOrtho();

		newTriangleCenterCrosses[f_id] = gmds::math::Cross2D(center_vect, center_vect_second);
	}

	newBdryEdgeNormals.resize(newLocalMesh->getNbEdges());
	gmds::math::Vector3d zeroVec(0.0, 0.0, 0.0);
	newBdryNodeNormals.resize(newLocalMesh->getNbNodes(), zeroVec);

	for (auto e_id:newLocalMesh->edges()){
		gmds::Edge currentEdge = newLocalMesh->get<gmds::Edge>(e_id);
		if (isCurveEdge[currentEdge.id()]) {
			vector<gmds::Node> adjacent_nodes = currentEdge.get<gmds::Node>();
			math::Point p1 = adjacent_nodes[0].getPoint();
			math::Point p2 = adjacent_nodes[1].getPoint();
			math::Vector3d v1 = math::Vector3d(p1, p2);
			v1.normalize();
			vector<gmds::Face> adjacent_faces = currentEdge.get<gmds::Face>();
			newBdryEdgeNormals[e_id] = v1.cross(adjacent_faces[0].normal());
		}
	}

	for(auto n_id:newLocalMesh->nodes()){
		gmds::Node currentNode = newLocalMesh->get<gmds::Node>(n_id);
		if(isCurveNode[currentNode.id()]){
			vector<gmds::Edge> currentEdges = currentNode.get<gmds::Edge>();
			for(unsigned int i=0; i<currentEdges.size(); i++){
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
void SingularityGraphBuilder2D::connectSingularityLines(math::Vector3d             dir_slot_i,
                                                        math::Vector3d             dir_slot_j,
                                                        math::Vector3d             connection_dir,
                                                        math::Point                start_pnt1,            // start_pnt = line_discretization[0]; //starting point
                                                        math::Point                start_pnt2,
                                                        TCellID                    start_cell_id1,
                                                        TCellID                    start_cell_id2,
                                                        int                        start_cell_dim1,
                                                        int                        start_cell_dim2,
                                                        std::vector<math::Point>&  APoints,
                                                        std::vector<TCellID>&      ATriangles,
                                                        double&                    streamlineDeviation)
{
	cout<<"connectSingularityLines "<<endl;

	/* "transport" the last closest comp vect (to dir_slot_i) until singline2;
	* if it corresponds with (opposite) dir_slot_j, we can connect.
	* the measure could be the angle diff along the connection_dir */

	ATriangles.clear();
	APoints.clear();

	math::Point  current_pnt = start_pnt1;
	math::Vector3d current_vec = connection_dir;

	math::Point  start_pnt = start_pnt1;
	math::Vector3d start_dir = connection_dir;

	TCellID             start_cell_id = start_cell_id1;
	int                 start_cell_dim = start_cell_dim1;

	math::Point start_dirPnt(start_pnt1.X() + connection_dir.X(),
	                         start_pnt1.Y() + connection_dir.Y(),
	                         start_pnt1.Z() + connection_dir.Z());

	std::vector<math::Vector3d> compVectors;
	math::Cross2D cross_cell;
	if(start_cell_dim1==0){// on node
		compVectors = (*m_field)[start_cell_id1].componentVectors();
		cross_cell = (*m_field)[start_cell_id1];
	}
	else{//on edge
		vector<Node> currentNodes = (m_mesh->get<Edge>(start_cell_id1)).get<Node>();
		math::Vector3d AB(currentNodes[0].getPoint(), currentNodes[1].getPoint());
		math::Vector3d AC(currentNodes[0].getPoint(), start_pnt1);
		double alpha = AC.norm()/AB.norm();

		std::vector<math::Vector3d> compVectorsA, compVectorsB;
		math::Cross2D cross_cellA, cross_cellB;
		compVectorsA = (*m_field)[currentNodes[0].id()].componentVectors();
		compVectorsB = (*m_field)[currentNodes[1].id()].componentVectors();
		cross_cellA = (*m_field)[currentNodes[0].id()];
		cross_cellB = (*m_field)[currentNodes[1].id()];
		compVectors.resize(4);
		for(unsigned int i=0;i<4;i++){
			compVectors[i] = alpha * compVectorsA[i] + (1.0-alpha) * compVectorsB[i];
			compVectors[i].normalize();
		}
		cross_cell = math::Cross2D::mean(cross_cellA, alpha, cross_cellB, (1.0-alpha));
	}

	math::Vector3d closestCompVect = compVectors[0];
	TCoord val = dir_slot_i.dot(closestCompVect);
	for(unsigned int i=1;i<4;i++){
		math::Vector3d v_i = compVectors[i];
		TCoord val_i=connection_dir.dot(v_i);
		if(val_i> val){
			val = val_i;
			closestCompVect = v_i;
		}
	}

	//========================================================================
	// Main loop to create the singularity line
	//========================================================================
	TCellID next_cell_id  = NullID;
	while(start_cell_id2!=next_cell_id) {
		next_cell_id  = NullID;
		int     next_cell_dim = -1;
		m_tool.findNextCell(start_pnt, start_dir,
		                    start_cell_dim, start_cell_id,
		                    next_cell_dim, next_cell_id);

		if(next_cell_dim == -1){
			throw GMDSException("connection tries to pass through concavity");
		}
		else{
			if (next_cell_dim == 1){
				/* we are going along an edge. Our simple assumption is to follow this edge until
				* reaching one of its end points and to compute the next direction at this point.*/
				Edge currentEdge = m_mesh->get<Edge>(next_cell_id);
				std::vector<TCellID> adj_faces = currentEdge.getIDs<Face>();
				ATriangles.insert(ATriangles.end(),adj_faces.begin(),adj_faces.end());

				std::vector<Node> currentNodes = currentEdge.get<Node>();
				math::Vector3d v0(start_pnt, currentNodes[0].getPoint());
				math::Vector3d v1(start_pnt, currentNodes[1].getPoint());
				Node next_node;
				if(math::near(v0.norm(),0.0))
					next_node = currentNodes[1];
				else if(math::near(v1.norm(),0.0))
					next_node = currentNodes[0];
				else if(v0.dot(start_dir) > v1.dot(start_dir))
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
				for(unsigned int i=1;i<4;i++){
					math::Vector3d v_i = compVectors[i];
					TCoord val_i=connection_dir.dot(v_i);
					if(val_i> val){
						val = val_i;
						closestCompVect = v_i;
					}
				}
			}
			else { //general case, we are in a face
				Face currentFace = m_mesh->get<Face>(next_cell_id);
				ATriangles.push_back(currentFace.id());

				//==============================================================
				// CASE 2: SO, WE JUST HAVE TO CROSS THE TRIANGLE (GENERAL CONF)
				//==============================================================

				math::Point  out_pnt, min_out_pnt;
				TCellID out_cell_id, min_out_cell_id;
				int out_cell_dim, min_out_cell_dim;

				vector<gmds::Node> currentNodes = currentFace.get<gmds::Node>();
				vector<gmds::Edge> currentEdges = currentFace.get<gmds::Edge>();

				double min_value = 7.7;
				bool found = false;

				for(unsigned int t=0; t<3; t++){
					math::Point node_loc = currentNodes[t].getPoint();
					math::Vector3d trial_vector(start_pnt, node_loc);
					trial_vector.normalize();
					double temp_dot = trial_vector.dot(connection_dir);
					if(temp_dot>min_value){
						min_value = temp_dot;
						min_out_cell_dim = 0;
						min_out_cell_id = currentNodes[t].id();
						min_out_pnt = node_loc;
					}

					if (math::near(temp_dot - 1,0)) {
						out_cell_dim = 0;
						out_cell_id = currentNodes[t].id();
						out_pnt = node_loc;
						found = true;
						break;
					}

					math::Vector3d out_vec;
					double deviation = 0.0;

					bool intersectEdge = m_tool.computeOutVectorFromRayAndEdge(currentEdges[t],
					                                                           start_pnt,
					                                                           connection_dir,
					                                                           out_pnt,
					                                                           out_vec,
					                                                           deviation);
					if(intersectEdge){
						out_cell_dim = 1;
						out_cell_id = currentEdges[t].id();
					}
				}

				if(!found){
					out_cell_dim = min_out_cell_dim;
					out_cell_id = min_out_cell_id;
 					out_pnt = min_out_pnt;
				 }

				APoints.push_back(out_pnt);

				start_pnt = out_pnt;

				start_cell_dim = out_cell_dim;
				start_cell_id = out_cell_id;

				if(start_cell_dim==0){// on node
					compVectors = (*m_field)[start_cell_id].componentVectors();
					cross_cell = (*m_field)[start_cell_id];
	 			}
	 			else{//on edge
					vector<Node> currentNodes = (m_mesh->get<Edge>(start_cell_id)).get<Node>();
					math::Vector3d AB(currentNodes[0].getPoint(), currentNodes[1].getPoint());
					math::Vector3d AC(currentNodes[0].getPoint(), start_pnt);
					double alpha = AC.norm()/AB.norm();

					std::vector<math::Vector3d> compVectorsA, compVectorsB;
					math::Cross2D cross_cellA, cross_cellB;
					compVectorsA = (*m_field)[currentNodes[0].id()].componentVectors();
					compVectorsB = (*m_field)[currentNodes[1].id()].componentVectors();
					cross_cellA = (*m_field)[currentNodes[0].id()];
					cross_cellB = (*m_field)[currentNodes[1].id()];
					compVectors.resize(4);
					for(unsigned int i=0;i<4;i++){
						compVectors[i] = alpha * compVectorsA[i] + (1.0-alpha) * compVectorsB[i];
						compVectors[i].normalize();
					}
					cross_cell = math::Cross2D::mean(cross_cellA, alpha, cross_cellB, (1.0-alpha));
	 			}
				math::Vector3d closestCompVect = compVectors[0];
				TCoord val = dir_slot_i.dot(closestCompVect);
				for(unsigned int i=1;i<4;i++){
					math::Vector3d v_i = compVectors[i];
					TCoord val_i = connection_dir.dot(v_i);
					if(val_i> val){
						val = val_i;
						closestCompVect = v_i;
					}
				}
			}
		}
	}
}

/*---------------------------------------------------------------------------*/
void SingularityGraphBuilder2D::createSingularityLinesSimultaneousStartRK4()
{
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

	cout<<"createSingularityLinesSimultaneousStartRK4"<<endl;
	/* when deciding whether to connect 2 singularity lines we have two options:
	1. connectByField = false; connect depending on the geometrical distance between the extremities of the singularity lines and depending on the angle made by those
	2. connectByField = true; connect depending on the frame field; consider the straight line between the extremities of the 2 singularity lines - this line will intersect edges (and vertices rarely); calculate the closest component vector at the intersection point with regard to the straight line's direction; if the closest component vector coincides (end to end) - we can connect; (by coincides -  consider direction 0-2 and direction 1-3)

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

	double stepSize; // search step for streamlines
	double thresholdStreamLineDist; // if 2 lines are within this threshold from one another => connect them
	double weightAngle = 0.75;
	double weightDistance = 0.25;

	stepSize = mean_edge_length;
	// WARNING stepSize should be sufficiently small
	//WARNING thresholdStreamLineDist should be sufficiently big
	thresholdStreamLineDist = 100*mean_edge_length*ATolerance;

	vector<SingularityPoint::Slot*> searchSlots;
	vector<SingularityPoint*> singularity_points = m_graph.getPoints();

	for (unsigned int i = 0; i < singularity_points.size(); i++) {
		SingularityPoint* pi = singularity_points[i];
		if(pi->getType()==0){
			std::vector<SingularityPoint::Slot*> pi_slots = pi->getSlots();
			for (unsigned int j = 0; j < pi_slots.size(); j++) {
				if (!pi_slots[j]->isLaunched)
					searchSlots.push_back(pi_slots[j]);
			}
		}
	}

	unsigned int NoSimLines = searchSlots.size();
	vector<vector<math::Point>> line_discretization(NoSimLines, vector<math::Point>(0)), copy_line_discretization(NoSimLines, vector<math::Point>(0));
	vector<vector<TCellID>>     line_triangles(NoSimLines, vector<TCellID>(0)), copy_line_triangles(NoSimLines, vector<TCellID>(0));
	vector<double> accumulatedDistancePerSlotLine(NoSimLines, 0.0);
	bool 					find_end_bdry = false;
	bool 					end_on_free_slot = false;
	vector<TCellID> 		to_cell_id(NoSimLines);
	vector<int> 			to_cell_dim(NoSimLines);
	vector<math::Vector3d> 	to_dir(NoSimLines);
	vector<double>            streamlineDeviation(NoSimLines, 0.0);

	unsigned int noFrozenSLots = 0;

	for(unsigned int i=0; i<NoSimLines; i++){
		line_discretization[i].push_back(searchSlots[i]->from_point->getLocation());
		line_discretization[i].push_back(searchSlots[i]->location);
		to_cell_id[i]  = searchSlots[i]->starting_cell_id ;
		to_cell_dim[i] = searchSlots[i]->starting_cell_dim;
		to_dir[i] = searchSlots[i]->direction; //starting direction
	}

	// max contor stop condition should take into account the mean_edge_length relative to the mesh...
	int contor = 0;

	math::Vector3d surfNormal(1.0, 1.0, 1.0);
	vector<pair<pair<unsigned int, unsigned int>, pair<double, double>>> possibleConnectingLines; /*vector storing ((slotA, slotB) , (distance between them, angle of intersection))*/
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

	Variable<int>* ball_var_sim = m_mesh->newVariable<int,GMDS_FACE>("sing_ball_sim");

	for(auto f_id:m_mesh->faces()){
		(*ball_var_sim)[f_id] = 0;
	}

	while(noFrozenSLots<NoSimLines){
		contor++;
		possibleConnectingLines.clear();
		stopIncrease.clear();
		stopIncrease.resize(NoSimLines, false);
		for(unsigned int i=0; i<NoSimLines; i++){
			SingularityPoint       *to_sing_pnt  = 0;
			SingularityPoint::Slot *to_slot = 0;
			SingularityPoint::Slot* current_slot = searchSlots[i];
			if((!current_slot->isFreeze)&&(!stopIncrease[i])){
				find_end_bdry = false;
				end_on_free_slot = false;
				bool find_end = false;

				if(contor==1){
					// here it can also go out through a node<!!! add triangles
					Edge currentEdge = m_mesh->get<Edge>(to_cell_id[i]);
					std::vector<TCellID> adj_faces = currentEdge.getIDs<Face>();
					bool AOnEdge0, AOnEdge1, AOnEdge2;
					vector<double> lambdas(3);

					if(m_tool.isPntInTri(searchSlots[i]->from_point->getLocation(), m_mesh->get<Face>(adj_faces[0]), AOnEdge0, AOnEdge1, AOnEdge2, lambdas[0], lambdas[1])){
						//to_cell_id[i] = adj_faces[1];
						line_triangles[i].push_back(adj_faces[0]);
						line_triangles[i].push_back(adj_faces[1]);
					}
					else{
						// to_cell_id[i] = adj_faces[0];
						line_triangles[i].push_back(adj_faces[1]);
						line_triangles[i].push_back(adj_faces[0]);
					}
				}

				growLineRK4(current_slot->from_point,
				            current_slot,
				            to_sing_pnt,
				            to_slot,
				            line_discretization[i][line_discretization[i].size()-1],
				            to_dir[i],
				            line_discretization[i],
				            line_triangles[i],
				            to_cell_dim[i],
				            to_cell_id[i],
				            streamlineDeviation[i],
				            stepSize,
				            find_end_bdry,
				            end_on_free_slot,
				            find_end);

				if(to_slot!=0){
					cout<<"current_slot->location "<<current_slot->location<<"; to_slot->location "<<to_slot->location<<endl;
				}
				// first check if find_end_bdry or end_on_free_slot(unlikely)
				if(find_end_bdry){
					//======================================================================
					// CASE 1 - We finish on the boundary. A geometric point must be created
					// The cell defined by (start_cell_dim[i], start_cell_id[i]) is on the boundary.
					// We have to create a geometric singularity point.
					// however we could do this at the end of the entire iteration(for the current searchStep)
					//line creation
					if(withGlobalComments)
						cout<<"!!!find_end_bdry for current_slot->from_point "<<current_slot->from_point->getLocation().X()<<" "<<current_slot->from_point->getLocation().Y()<<
						" and with dir "<<current_slot->direction[0]<<" "<<current_slot->direction[1]<<endl;
					SurfaceSingularityLine* surf_line = m_graph.newSurfaceLine();
					int sepNumberTmp = m_graph.getNbLines();
					surf_line->setNumber(sepNumberTmp);
					SingularityPoint* from_sing_pnt = current_slot->from_point;
					surf_line->addSingularityPoint(from_sing_pnt);
					surf_line->addDiscretizationPoint(from_sing_pnt->getLocation());
					current_slot->line = surf_line;

					math::Vector3d firstDir = math::Vector3d(from_sing_pnt->getLocation(), line_discretization[i][0]);

					current_slot->line_direction = firstDir;
					current_slot->isLaunched = true;
					if(to_slot!=0){
						streamlineDeviation[i] = (streamlineDeviation[i] + fabs(1 - to_slot->direction.dot(math::Vector3d(line_discretization[i][line_discretization[i].size()-1],line_discretization[i][line_discretization[i].size()-2]))));
						streamlineDeviation[i] = streamlineDeviation[i]/(line_discretization[i].size()+1);
						current_slot->lineDeviation = streamlineDeviation[i];
					}
					//Insertion of line points
					for(unsigned int j=0; j<line_discretization[i].size(); j++) {
						surf_line->addDiscretizationPoint(line_discretization[i][j]);
					}

					for(unsigned int j=0; j<line_triangles[i].size(); j++) {
						surf_line->addTraversedFace(line_triangles[i][j]);
					}

					SingularityPoint* geom_pnt;
					SingularityPoint::Slot* incoming_slot;
					createGeometricSingularityPoint(line_discretization[i][line_discretization[i].size()-1],      // the last point added
					                                to_dir[i],      // the direction we come from
					                                to_cell_dim[i], // the dim. of the cell
					                                to_cell_id[i],  // the id of the cell
					                                geom_pnt,       // the created point
					                                incoming_slot); // and the slot

					current_slot->lineDeviation = streamlineDeviation[i]/line_discretization[i].size();

					surf_line->addSingularityPoint(geom_pnt);
					noFrozenSLots++;
					current_slot->isFreeze = true;
					current_slot->isLaunched = true;

					for(unsigned int j=0; j<possibleConnectingLines.size(); j++){
						if((i==possibleConnectingLines[j].first.first)||((i==possibleConnectingLines[j].first.second))){
							possibleConnectingLines.erase(possibleConnectingLines.begin()+j);
						}
					}

					gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
					for (unsigned int j = 0; j < line_discretization[i].size()-1; j++){
						gmds::Node n1 = m.newNode(line_discretization[i][j].X(), line_discretization[i][j].Y(), line_discretization[i][j].Z());
						gmds::Node n2 = m.newNode(line_discretization[i][j + 1].X(), line_discretization[i][j + 1].Y(), line_discretization[i][j + 1].Z());
						gmds::Face f  =  m.newTriangle(n1, n1, n2);
					}
					gmds::IGMeshIOService ioService(&m);
					gmds::VTKWriter vtkWriter(&ioService);
					vtkWriter.setCellOptions(gmds::N|gmds::F);
					vtkWriter.setDataOptions(gmds::N|gmds::F);
					std::string file_name = "SimultaneousStartRK4_SingToBdry"+to_string(i)+".vtk";
					vtkWriter.write(file_name);
					cout<<"has written SingToBdry"<<file_name<<endl;
				}
				else{// not find_end_bdry
		  			if(end_on_free_slot){
						if(withGlobalComments){
							cout<<"end_on_free_slot; to_slot "<<to_slot<<endl;
							cout<<"current_slot->location "<<current_slot->location<<endl;
							cout<<"line disret "<<line_discretization[i][line_discretization[i].size()-2];
							cout<<" "<<line_discretization[i][line_discretization[i].size()-1]<<endl;
						}
						SurfaceSingularityLine* surf_line = m_graph.newSurfaceLine();
						int sepNumberTmp = m_graph.getNbLines();
						surf_line->setNumber(sepNumberTmp);

						SingularityPoint* from_sing_pnt = current_slot->from_point;
						surf_line->addSingularityPoint(from_sing_pnt);
						surf_line->addDiscretizationPoint(from_sing_pnt->getLocation());

						current_slot->line = surf_line;
						math::Vector3d firstDir = math::Vector3d(from_sing_pnt->getLocation(), line_discretization[i][0]);

						current_slot->line_direction = firstDir;
						current_slot->isLaunched = true;

						current_slot->lineDeviation = streamlineDeviation[i];
						for(unsigned int j=0; j<line_discretization[i].size(); j++) {
							surf_line->addDiscretizationPoint(line_discretization[i][j]);
						}

						for(unsigned int j=0; j<line_triangles[i].size(); j++) {
							surf_line->addTraversedFace(line_triangles[i][j]);
						}

						current_slot->isLaunched = true;
						to_slot->isLaunched = true;
						current_slot->isFreeze = true;
						to_slot->isFreeze = true;
						to_slot->line = surf_line;
						to_slot->line_direction = to_dir[i];
						to_slot->lineDeviation = streamlineDeviation[i];

						surf_line->addSingularityPoint(to_sing_pnt);
						surf_line->addDiscretizationPoint(to_sing_pnt->getLocation());
						to_slot->line_direction = to_dir[i];
						noFrozenSLots = noFrozenSLots + 2;

						gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
						for(unsigned int j = 0; j < line_discretization[i].size()-1; j++){
							gmds::Node n1 = m.newNode(line_discretization[i][j].X(), line_discretization[i][j].Y(), line_discretization[i][j].Z());
							gmds::Node n2 = m.newNode(line_discretization[i][j + 1].X(), line_discretization[i][j + 1].Y(), line_discretization[i][j + 1].Z());
							gmds::Face f  =  m.newTriangle(n1, n1, n2);
						}
						gmds::IGMeshIOService ioService(&m);
						gmds::VTKWriter vtkWriter(&ioService);
						vtkWriter.setCellOptions(gmds::N|gmds::F);
						vtkWriter.setDataOptions(gmds::N|gmds::F);
						std::string file_name = "SimultaneousStartRK4_SingToSingBall"+to_string(i)+"-x.vtk";
						vtkWriter.write(file_name);
						cout<<"has written SingToSing "<<file_name<<endl;
					}
					else{
						cout<<"see if we got in the vecinity of another line"<<endl;
						//======================================================================
						//CASE 2 - We encounter another line within the thresholdStreamLineDist
						// for now, we will simply connect the last points added to the 2 lines (which are within thresholdStreamLineDist);
						// connect only if deviation for the connecting line is sufficiently small; otherwise,let it arrive at boundary
						for(unsigned int j=0; j<NoSimLines; j++){
							if((searchSlots[i]->from_point!=searchSlots[j]->from_point)&&(line_discretization[j].size()>2)&&(!searchSlots[j]->isFreeze)){
							distBtwSlots = line_discretization[i][line_discretization[i].size()-1].distance(line_discretization[j][line_discretization[j].size()-1]);
								if(distBtwSlots<=thresholdStreamLineDist){
									bool forbiddenPair = false;
									for(unsigned int ti0=0; ti0<forbiddenConnectingLines.size(); ti0++){
										if((forbiddenConnectingLines[ti0].first==i)&&(forbiddenConnectingLines[ti0].second==j)){
											forbiddenPair = true;
											break;
										}
									}
									if(!forbiddenPair){
										math::Point last_added_pnt = line_discretization[i][line_discretization[i].size()-1];
										math::Vector3d start_dir = math::Vector3d(last_added_pnt, line_discretization[j][line_discretization[j].size()-1]);
										if(!connectByField){
											double tempDev = 0.0;
											vector<gmds::Face> candidate_faces;
											vector<bool> visitedFaces(original_faces_number,false), addedFaces(original_faces_number,false);
											if(to_cell_dim[i]==0){//node
												gmds::Node currentNode = m_mesh->get<gmds::Node>(to_cell_id[i]);
												candidate_faces = currentNode.get<gmds::Face>();
												for(unsigned int t=0; t<candidate_faces.size(); t++){
													visitedFaces[candidate_faces[t].id()] = true;
												}
											}
											else{
												if(to_cell_dim[i]==1){//edge
													gmds::Edge currentEdge = m_mesh->get<gmds::Edge>(to_cell_id[i]);
													candidate_faces = currentEdge.get<gmds::Face>();
													for(unsigned int t=0; t<candidate_faces.size(); t++){
														visitedFaces[candidate_faces[t].id()] = true;
													}
												}
												else{//to_cell_dim[i]==2 - face
													vector<gmds::Node> currentNodes = m_mesh->get<gmds::Face>(to_cell_id[i]).get<Node>();
													for(unsigned int tt =0; tt<3; tt++){
														vector<gmds::Face> temp_faces = currentNodes[tt].get<gmds::Face>();
														for(unsigned int tt2 =0; tt2<temp_faces.size(); tt2++){
															if(!visitedFaces[temp_faces[tt2].id()]){
																visitedFaces[temp_faces[tt2].id()] = true;
																candidate_faces.push_back(temp_faces[tt2]);
															}
														}
													}
												}
											}

											for(unsigned int t=0; t<line_triangles[i].size(); t++)
												addedFaces[line_triangles[i][t]] = true;
											for(unsigned int t=0; t<line_triangles[j].size(); t++)
												addedFaces[line_triangles[j][t]] = true;

											math::Segment seg1(last_added_pnt, line_discretization[j][line_discretization[j].size()-1]);
											math::Ray from_ray(last_added_pnt,line_discretization[j][line_discretization[j].size()-1]);
											math::Vector3d connLineDir =	from_ray.getDirUnit();
											math::Vector3d closest2Current, closest0, closest1;
											unsigned int t = 0;
											bool stopAdding = false;
											vector<bool> visitedEdges(m_mesh->getNbEdges(),false);
											unsigned int noIntEdges = 0;

											copy_line_triangles[i] = line_triangles[i];

											while(t<candidate_faces.size()){
												Face currentFace = candidate_faces[t];
												vector<gmds::Edge> currentEdges = currentFace.get<gmds::Edge>();
												for(unsigned int tt=0; tt<3; tt++){
													vector<gmds::Node> currentNodes = currentEdges[tt].get<gmds::Node>();
													math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
													math::Point intersectionPnt;
													double intersectionParam;
													if(seg1.SecondMetIntersect2D(oppSeg, intersectionPnt, intersectionParam, temp_epsilon)){
														if(from_ray.SecondMetIntersect2D(oppSeg, intersectionPnt, intersectionParam, temp_epsilon)){
															if(!visitedEdges[currentEdges[tt].id()]){
																math::Cross2D cross_0 = (*m_field)[currentNodes[0].id()];
																math::Cross2D cross_1  = (*m_field)[currentNodes[1].id()];
																closest0 =  cross_0.closestComponentVector(connLineDir);
																closest1 =  cross_1.closestComponentVector(connLineDir);
																closest2Current = intersectionParam * closest0 + (1 - intersectionParam) * closest1;
																tempDev = tempDev + fabs(1.0 - connLineDir.dot(closest2Current));
																visitedEdges[currentEdges[tt].id()] = true;
																noIntEdges++;
															}
															if(currentFace.id() == line_triangles[j][line_triangles[j].size()-1])
																stopAdding = true;
															if(!addedFaces[currentFace.id()]){
																copy_line_triangles[i].push_back(currentFace.id());
																addedFaces[currentFace.id()] = true;
															}
															if(!stopAdding){
																vector<Face> adj_faces0 = currentNodes[0].get<Face>();
																vector<Face> adj_faces1 = currentNodes[1].get<Face>();
																adj_faces0.insert(adj_faces0.end(), adj_faces1.begin(), adj_faces1.end());
																for(unsigned int ttt=0; ttt<adj_faces0.size(); ttt++){
																	if(!visitedFaces[adj_faces0[ttt].id()]){
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

											cout<<"slots "<<i<<" and "<<j<<" to be connected with deviation "<<fmod(tempDev, 1.5707963267948966)<<endl;
											math::Vector3d test1(line_discretization[i][line_discretization[i].size()-2], line_discretization[i][line_discretization[i].size()-1]);
											math::Vector3d test2(line_discretization[j][line_discretization[j].size()-2], line_discretization[j][line_discretization[j].size()-1]);
											test1.normalize();
											test2.normalize();
											// if close to perpendicular=> shouldn't connect
											tempDev = test1.angle(test2); // tempDevin [0,PI]]
											//we must also take into account the distance; 2 streamlines could be parallel but far from one another(although inside the proximityRadius)
											math::Vector3d   tempCross = test1.cross(test2);
											// we connect only if the dot product with the normal (in 2d case - (1,1,1) is higher that 0) - in order to account for orientation
											if(tempCross.dot(surfNormal)>=0){
												tempDev = tempDev/(double)(math::Constants::PI);

												if((tempDev>=0.75)&&(tempDev<=1.25)){
													possibleConnectingLines.push_back(make_pair(make_pair(i, j),make_pair(distBtwSlots, std::fabs(1.0 - tempDev))));
													//ideally the angle btw the lines should be 0; therefore tempDev = 1
													stopIncrease[j] = true;
												}
											}
										}
										else{
										  //connectByField hasnt been implemented
										}
									}
								}
							}
						}
					}
				}

				gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
				for(unsigned int ti=0;ti<NoSimLines;ti++){
					for(unsigned int ti1=0;ti1<line_discretization[ti].size()-1;ti1++){
						gmds::Node n1 = m.newNode(line_discretization[ti][ti1].X(), line_discretization[ti][ti1].Y(), line_discretization[ti][ti1].Z());
						gmds::Node n2 = m.newNode(line_discretization[ti][ti1+ 1].X(), line_discretization[ti][ti1 + 1].Y(), line_discretization[ti][ti1 + 1].Z());
						gmds::Face f  =  m.newTriangle(n1, n1, n2);
					}
				}

				gmds::IGMeshIOService ioService(&m);
				gmds::VTKWriter vtkWriter(&ioService);
				vtkWriter.setCellOptions(gmds::N|gmds::F);
				vtkWriter.setDataOptions(gmds::N|gmds::F);
				std::string test_file_name = m_output_directory_name +"SimultaneousStartRK4_-lastStep"+".vtk";
				vtkWriter.write(test_file_name);
			}
		}

		// if one slot is within radius with 2 other slots, choose the one with the lowest (geometric) deviation

		if(possibleConnectingLines.size()>0){
			unsigned int t1=0;
			unsigned int t2;
			double connectionTerm1, connectionTerm2;

			while(t1<possibleConnectingLines.size()-1){
				t2 = t1+1;
				while(t2<possibleConnectingLines.size()){
					if(((possibleConnectingLines[t1].first.first==possibleConnectingLines[t2].first.first)&&
					(possibleConnectingLines[t1].first.second==possibleConnectingLines[t2].first.second))||
					((possibleConnectingLines[t1].first.first==possibleConnectingLines[t2].first.second)&&
					(possibleConnectingLines[t1].first.second==possibleConnectingLines[t2].first.first))){

						cout<<"connection between "<<possibleConnectingLines[t1].first.first<<" and "
						<<  possibleConnectingLines[t1].first.second<<endl;
						connectionTerm1 = weightDistance * (possibleConnectingLines[t1].second.first/(possibleConnectingLines[t1].second.first + possibleConnectingLines[t2].second.first));
						connectionTerm1 = connectionTerm1 + weightAngle * (possibleConnectingLines[t1].second.second/(possibleConnectingLines[t1].second.second + possibleConnectingLines[t2].second.second));
						connectionTerm2 = weightDistance * (possibleConnectingLines[t2].second.first/(possibleConnectingLines[t1].second.first + possibleConnectingLines[t2].second.first));
						connectionTerm2 = connectionTerm2 + weightAngle * (possibleConnectingLines[t2].second.second/(possibleConnectingLines[t1].second.second + possibleConnectingLines[t2].second.second));

						if(connectionTerm1>connectionTerm2){  // remove t1
							possibleConnectingLines.erase(possibleConnectingLines.begin()+t1);
							t1--;
							break;
						}
						else{
							possibleConnectingLines.erase(possibleConnectingLines.begin()+t2);
							t2--;
						}
					}
					t2++;
				}
				t1++;
			}
		}

		for(unsigned int t2=0; t2<possibleConnectingLines.size(); t2++){
			unsigned int i = possibleConnectingLines[t2].first.first;
			unsigned int j = possibleConnectingLines[t2].first.second;
			if(visualizeSingBallSim){
				gmds::math::Point centerPoint = line_discretization[j].back() + line_discretization[i].back();
				centerPoint.X() = centerPoint.X()/2;
				centerPoint.Y() = centerPoint.Y()/2;
				for(unsigned int t5=0; t5<original_faces_number-1; t5++){
					gmds::Face testFace1 = m_mesh->get<Face>(t5);
					if(triangle_centers[t5].distance(centerPoint)<thresholdStreamLineDist){
						for(unsigned int t6=t5+1; t6<original_faces_number; t6++){
							gmds::Face testFace2 = m_mesh->get<Face>(t6);
							if(triangle_centers[t6].distance(centerPoint)<thresholdStreamLineDist){
								if(triangle_centers[t5].distance(triangle_centers[t6])<thresholdStreamLineDist){
									(*ball_var_sim)[t5] = 1;
									(*ball_var_sim)[t6] = 1;
								}
							}
						}
					}
				}
			}

			// the 2 lines are within radius, they must be connected
			for(int t=line_discretization[j].size()-1; t>=0; t--){
				line_discretization[i].push_back(line_discretization[j][t]);
			}
			if(withGlobalComments){
				cout<<"creates line between "<<i<<" and "<<j<<" with angle 1+ "<<possibleConnectingLines[t2].second.second<<endl;
				cout<<"and distance "<<possibleConnectingLines[t2].second.first<<endl;
				cout<<"these should be diff(if small triangle count, not necessarilly) "<<line_triangles[i].size()<<"<"<<copy_line_triangles[i].size()<<endl;
			}
			line_triangles[i] = copy_line_triangles[i];
			for(int t=line_triangles[j].size()-1; t>=0; t--){
				line_triangles[i].push_back(line_triangles[j][t]);
			}

			SurfaceSingularityLine* surf_line = m_graph.newSurfaceLine();
			int sepNumberTmp = m_graph.getNbLines();
			surf_line->setNumber(sepNumberTmp);
			//connect line to initial singularity point
			SingularityPoint* from_sing_pnt = searchSlots[i]->from_point;
			SingularityPoint::Slot* the_other_slot = searchSlots[j];
			SingularityPoint* towards_sing_pnt = the_other_slot->from_point;
			surf_line->addSingularityPoint(from_sing_pnt);
			surf_line->addSingularityPoint(towards_sing_pnt);
			//surf_line->addDiscretizationPoint(from_sing_pnt->getLocation());
			math::Vector3d firstDir = math::Vector3d(line_discretization[i][0], line_discretization[i][1]);
			searchSlots[i]->line_direction = firstDir;

			firstDir = math::Vector3d(line_discretization[j][0], line_discretization[j][1]);
			searchSlots[j]->line_direction = firstDir;

			streamlineDeviation[i] = streamlineDeviation[i] + streamlineDeviation[j];
			streamlineDeviation[i] = streamlineDeviation[i]/(line_discretization[i].size()+1);
			streamlineDeviation[j] = streamlineDeviation[i] ;
			//WARNING also compute streamlineDeviation locally, inside thresholdStreamLineDist
			searchSlots[i]->lineDeviation = streamlineDeviation[i];
			searchSlots[j]->lineDeviation = streamlineDeviation[i];

			for(unsigned int t=0; t<line_discretization[i].size(); t++) {
				surf_line->addDiscretizationPoint(line_discretization[i][t]);
			}

			for(unsigned int t=0; t<line_triangles[i].size(); t++) {
				surf_line->addTraversedFace(line_triangles[i][t]);
			}
			searchSlots[i]->line = surf_line;
			searchSlots[j]->line = surf_line;//  the inverse of this; although in the original algorithm it's also surf_line;

			searchSlots[i]->isLaunched = true;
			searchSlots[i]->isFreeze = true;
			searchSlots[j]->isLaunched = true;
			searchSlots[j]->isFreeze = true;

			/*gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
			 *for (unsigned int j = 0; j < line_discretization[i].size()-1; j++){
		    		gmds::Node n1 = m.newNode(line_discretization[i][j].X(), line_discretization[i][j].Y(), line_discretization[i][j].Z());
		    		gmds::Node n2 = m.newNode(line_discretization[i][j + 1].X(), line_discretization[i][j + 1].Y(), line_discretization[i][j + 1].Z());
		    		gmds::Face f  =  m.newTriangle(n1, n1, n2);
			}
			gmds::IGMeshIOService ioService(&m);
			gmds::VTKWriter vtkWriter(&ioService);
			vtkWriter.setCellOptions(gmds::N|gmds::F);
			vtkWriter.setDataOptions(gmds::N|gmds::F);
			std::string file_name = "SingToSing"+to_string(i)+"-"+to_string(j)+".vtk";
			vtkWriter.write(file_name);
			cout<<"has written SingToSing"<<file_name<<endl;*/

			noFrozenSLots = noFrozenSLots + 2;
		}

		/*only for visualization purposes
    		for(unsigned int j0 = 0; j0 < line_discretization.size(); j0++){
			gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
			cout<<"line_discretization["<<j0<<"].size() "<<line_discretization[j0].size()<<endl;
			for (unsigned int j = 0; j < line_discretization[j0].size()-1; j++){
				gmds::Node n1 = m.newNode(line_discretization[j0][j].X(), line_discretization[j0][j].Y(), line_discretization[j0][j].Z());
		    		gmds::Node n2 = m.newNode(line_discretization[j0][j + 1].X(), line_discretization[j0][j + 1].Y(), line_discretization[j0][j + 1].Z());
		    		gmds::Face f  =  m.newTriangle(n1, n1, n2);
			}

			gmds::IGMeshIOService ioService(&m);
			gmds::VTKWriter vtkWriter(&ioService);
			vtkWriter.setCellOptions(gmds::N|gmds::F);
			vtkWriter.setDataOptions(gmds::N|gmds::F);
			std::string file_name = "SimultaneuosRK4_lastStep"+to_string(j0)+".vtk";
			vtkWriter.write(file_name);
		}    */
		writeOutputSingle("boundary_line");

		if(contor==m_mesh->getNbEdges()){
			cout<<"m_mesh->getNbEdges() "<<m_mesh->getNbEdges()<<endl;
			throw GMDSException("contor==total number of edges");
		}
	}

	/*for the sp-sp streamlines, remesh the area around them (all neighbouring faces of the
	 *traversed triangles), recompute cross field in order to have a final refined streamline*/
	/*vector<bool> visitedFaces(original_faces_number, false);

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
			gmds::Variable<gmds::math::Cross2D>* local_cross_field_2D = newLocalMesh.newVariable<math::Cross2D,GMDS_NODE>("cross_X");
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
        		visitedFaces.resize(original_faces_number, false);
     	}
   	}
	*/

	/*WARNING Although improbable, for the lines that have arrived within the confusing ball of a
	* different singularity, we haven't done anything; at the end, inspect if we have such cases and
	* if so, treat them*/

	if(visualizeSingBallSim){

		gmds::IGMeshIOService meshIoServ(m_mesh);
		gmds::VTKWriter writerB(&meshIoServ);
		writerB.setCellOptions(gmds::N|gmds::F);
		writerB.setDataOptions(gmds::N|gmds::F);
		std::stringstream file_name2;
		file_name2 <<m_output_directory_name<<"-confusing_balls_sim.vtk";
		writerB.write(file_name2.str());
		cout<<"wrote confusing_balls_sim.vtk"<<endl;
	}
}

/*--------------------------------------------------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::retraceShortestPath(gmds::TCellID&        foundTarget,
                                                    vector<int>&          previous,
                                                    vector<unsigned int>& path)
{
	// between foundTarget and the source (from previously called getShortestPathBtwFaces)
	gmds::TCellID tempSource = foundTarget;
	path.clear();

	while(tempSource!=original_faces_number){//WARNING check; before while(tempSource!=-1){
		path.push_back(tempSource);
		tempSource = previous[tempSource];
	}
}

/*--------------------------------------------------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::retraceShortestPath(int&                  foundTarget,
                                                    vector<int>&          previous,
                                                    vector<unsigned int>& path)
{
	// between foundTarget and the source (from previously called getShortestPathBtwFaces)
	int tempSource = foundTarget;
	path.clear();

	while(tempSource!=-1){
		path.push_back(tempSource);
		tempSource = previous[tempSource];
	}
}

/*--------------------------------------------------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::retraceShortestPath(gmds::TCellID&        foundTarget,
                                                    vector<int>&          previous,
                                                    vector<unsigned int>& path,
                                                    unsigned int&         tempFacesNo)
{
	//between foundTarget and the source (from previously called getShortestPathBtwFaces)
	gmds::TCellID tempSource = foundTarget;
	path.clear();

	while(tempSource<original_faces_number){
		path.push_back(tempSource);
		tempSource = previous[tempSource];
	}
}

/*--------------------------------------------------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::computeFace2FaceInfo(vector<vector<gmds::TCellID>>& newNode2NewNodeNeighbours,
                                                     vector<vector<double>>&        face2FaceTransport,
                                                     vector<vector<unsigned int>>&  face2FaceDeviation)
{
	cout<<"computeFace2FaceInfo"<<endl;
	/*
	* WARNING if 3D, FIELD SHOULD BE DEFINED PER FACE! ; WE ALSO MUST TAKE INTO ACCOUNT
	* THE MISSMATCH BETWEEN 2 TRIANGLES (ACROSS THE EDGE) - face2FaceDeviation
	* first the local basis is constructed for each triangle (local_basis)
	local_basis[triangle T] has 2 components: 2 orthogonal vectors in the triangle plane
	that are orthonormal to the normal of the triangle
	// here we assume the mesh is manifold
	face2FaceTransport - the angle between the two reference frames across and edge (shared by two adjacent triangles)
	if two vectors across this edge differ by this angle they are considered to be parallel.
	*/
	unsigned int noEdges = m_mesh->getNbEdges();
	gmds::math::Vector3d zeroVec(0.0, 0.0, 0.0);
	local_basis.clear();
	if(face_normals.size()==0)
		face_normals.resize(original_faces_number);
	if(face2Face_neighbours.size()==0)
		face2Face_neighbours.resize(original_faces_number);
	if(face2Face_neighbours_by_verts.size()==0)
		face2Face_neighbours_by_verts.resize(original_faces_number);
	if(face2FaceTransport.size()==0)
	 	face2FaceTransport.resize(original_faces_number, vector<double> (3));
	if(face2FaceDeviation.size()==0)
	 	face2FaceDeviation.resize(original_faces_number, vector<unsigned int> (3, -1));
	if(is_bdry_face.size()==0)
		is_bdry_face.resize(original_faces_number,false);
	if(bdry_edge_normals.size()==0)
		bdry_edge_normals.resize(noEdges);
	if(bdry_node_normals.size()==0)
		bdry_node_normals.resize(original_nodes_number, zeroVec);
	if(newNode2NewNodeNeighbours.size()==0)
		newNode2NewNodeNeighbours.resize(original_nodes_number + original_faces_number);

	vector<bool> isBorderEdge(noEdges, false);
	vector<bool> visitedFaces(original_faces_number, false);


	for(auto f_id:m_mesh->faces()){
		gmds::Face currentFace = m_mesh->get<Face>(f_id);
		vector<gmds::Node> currentNodes = currentFace.get<Node>();
		gmds::math::Vector3d firstEdge(currentNodes[0].getPoint(), currentNodes[1].getPoint());
		firstEdge.normalize();
		gmds::math::Vector3d lastEdgeInv(currentNodes[0].getPoint(), currentNodes[2].getPoint());
		face_normals[f_id] = currentFace.normal();
		gmds::math::Vector3d BasisY = (face_normals[f_id].cross(firstEdge));
		BasisY.normalize();
		local_basis.push_back(make_pair(firstEdge, BasisY));

		vector<gmds::Edge> currentEdges = currentFace.get<gmds::Edge>();

		visitedFaces.clear();
		visitedFaces.resize(original_faces_number, false);
		for(unsigned int i=0; i<3; i++){ //assume only triangles
			vector<gmds::Face> adjacent_triangles = currentEdges[i].get<gmds::Face>();
			if (adjacent_triangles.size()==2){
				if(adjacent_triangles[0].id()!=f_id)
					face2Face_neighbours[f_id].push_back(adjacent_triangles[0]);
				else
					face2Face_neighbours[f_id].push_back(adjacent_triangles[1]);
			}
			else{
				isBorderEdge[currentEdges[i].id()] = true;
			}

			if(m_mesh->isMarked(currentNodes[i],m_mark_nodes_on_curve) ||
			m_mesh->isMarked(currentNodes[i],m_mark_nodes_on_point)){
				is_bdry_face[f_id] = true;
			}
		}

		//TODO i have to compute this cheaper - disjoint intervals

		for(unsigned int i=0; i<3; i++){ //assume only triangles
			vector<gmds::Face> adjacent_triangles = currentNodes[i].get<gmds::Face>();
			newNode2NewNodeNeighbours[original_nodes_number + f_id].push_back(currentNodes[i].id());
			for(unsigned int j=0; j<adjacent_triangles.size(); j++){
				if((!visitedFaces[adjacent_triangles[j].id()])&&(adjacent_triangles[j].id()!=f_id)){
					bool validNeigh = true;

					if(m_mesh->isMarked(currentNodes[i], m_mark_nodes_on_point) ||
					m_mesh->isMarked(currentNodes[i], m_mark_nodes_on_curve)){

						math::Segment from_seg(triangle_centers[f_id], triangle_centers[adjacent_triangles[j].id()]);
						for(unsigned int k=0; k<3; k++){
							if (m_mesh->isMarked(currentEdges[k], m_mark_edges_on_curve)){
								vector<gmds::Node> temp_nodes =  currentEdges[k].get<gmds::Node>();
								math::Segment oppSeg(temp_nodes[0].getPoint(), temp_nodes[1].getPoint());
								gmds::math::Point intersectionPoint;
								double intersectionParam;
								if(from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon)){
									validNeigh = false;
								}
							}
						}
						if(validNeigh){
							vector<gmds::Edge> adj_edges = adjacent_triangles[j].get<gmds::Edge>();
							for(unsigned int k=0; k<3; k++){
								if (m_mesh->isMarked(adj_edges[k], m_mark_edges_on_curve)){
									vector<gmds::Node> temp_nodes =  adj_edges[k].get<gmds::Node>();
									math::Segment oppSeg(temp_nodes[0].getPoint(), temp_nodes[1].getPoint());
									gmds::math::Point intersectionPoint;
									double intersectionParam;
									if(from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon)){
										validNeigh = false;
									}
								}
							}
						}
					}
					if(validNeigh){
						face2Face_neighbours_by_verts[f_id].push_back(adjacent_triangles[j]);
						newNode2NewNodeNeighbours[original_nodes_number + f_id].push_back(original_nodes_number + adjacent_triangles[j].id());
					}
					visitedFaces[adjacent_triangles[j].id()] = true;
				}
			}
		}
		visitedFaces.clear();
		visitedFaces.resize(original_faces_number, false);

		//face2FaceTransport[f_id].resize(face2Face_neighbours[f_id].size(), 0.0);
		//face2FaceDeviation[f_id].resize(face2Face_neighbours[f_id].size());
	}

	//vector<double> K(noEdges, 0.0);
	// For every non-border edge

	for (auto e_id:m_mesh->edges()){
		gmds::Edge currentEdge = m_mesh->get<gmds::Edge>(e_id);
		if (!m_mesh->isMarked(currentEdge, m_mark_edges_on_curve)) {

			vector<gmds::Face> adjacent_triangles = currentEdge.get<gmds::Face>();
			unsigned int fid0 = adjacent_triangles[0].id();
			unsigned int fid1 = adjacent_triangles[1].id();

			math::Vector3d N0 = face_normals[fid0];

			// find common edge on triangle 0 and 1
			int fid0_vc = -1;//local id of common edge in f0 {0, 1 or 2}
			int fid1_vc = -1;//local id of common edge in f1
			vector<gmds::Edge> currentEdges0 = adjacent_triangles[0].get<gmds::Edge>();
			for(unsigned int i=0; i<3; i++){
				if(currentEdges0[i].id()==e_id)
					fid0_vc = i;
			}

			vector<gmds::Edge> currentEdges1 = adjacent_triangles[1].get<gmds::Edge>();
			for(unsigned int i=0; i<3; i++){
				if(currentEdges1[i].id()==e_id)
					fid1_vc = i;
			}

			vector<gmds::Node> edge_nodes0 = currentEdges0[fid0_vc].get<gmds::Node>();
			gmds::math::Vector3d common_edge(edge_nodes0[0].getPoint(), edge_nodes0[1].getPoint());
			common_edge.normalize();

			//Common local Basis (CLB) Map the two triangles in a new space where the common edge is the x axis and the N0 the z axis

			Eigen::MatrixXd CLB(3, 3);
			gmds::math::Point commonOrig = edge_nodes0[0].getPoint();// origin ; the first vert of common edge
			gmds::math::Vector3d tmp = common_edge.cross(face_normals[fid0]);

			for(unsigned int i=0; i<3; i++){
				CLB(0,i) = common_edge[i];
				CLB(1,i) = tmp[i];
				CLB(2,i) = face_normals[fid0][i];
			}
			vector<gmds::Node> face_nodes0 = adjacent_triangles[0].get<gmds::Node>();
			Eigen::MatrixXd V0(3, 3);
			for(unsigned int i=0; i<3; i++){
				V0(i,0) = face_nodes0[i].X() - commonOrig.X();
				V0(i,1) = face_nodes0[i].Y() - commonOrig.Y();
				V0(i,2) = face_nodes0[i].Z() - commonOrig.Z();
			}

			V0 = (CLB*V0.transpose()).transpose();

			vector<gmds::Node> face_nodes1 = adjacent_triangles[1].get<gmds::Node>();
			Eigen::MatrixXd V1(3, 3);

			for(unsigned int i=0; i<3; i++){
				V1(i,0) = face_nodes1[i].X() - commonOrig.X();
				V1(i,1) = face_nodes1[i].Y() - commonOrig.Y();
				V1(i,2) = face_nodes1[i].Z() - commonOrig.Z();
			}
			V1 = (CLB*V1.transpose()).transpose();

			// compute rotation R such that R * N1 = N0
			// i.e. map both triangles to the same plane
			double alpha = -atan2(V1((fid1_vc+2)%3,2),V1((fid1_vc+2)%3,1));

			Eigen::MatrixXd R(3, 3);
			//Eigen::Matrix<typename DerivedV::Scalar, 3, 3> R;
			R << 1,          0,            0,
				0, cos(alpha), -sin(alpha) ,
				0, sin(alpha),  cos(alpha);
			V1 = (R*V1.transpose()).transpose();

			// measure the angle between the reference frames
			// k_ij is the angle between the triangle on the left and the one on the right
			Eigen::MatrixXd ref0(1,3);
			Eigen::MatrixXd ref1(1,3);
	 		for(unsigned int i=0; i<3; i++){
				ref1(0,i) = V1(i,1) - V1(i,0);
				ref0(0,i) = V0(0,1) - V0(0,0);
			}

			ref0.normalize();
			ref1.normalize();

			double ktemp = atan2(ref1(1),ref1(0)) - atan2(ref0(1),ref0(0));

			// just to be sure, rotate ref0 using angle ktemp...
			Eigen::MatrixXd R2(2,2);
			R2 << cos(ktemp), -sin(ktemp), sin(ktemp), cos(ktemp);

			for(unsigned int i=0; i<face2Face_neighbours[fid0].size(); i++){
				if(face2Face_neighbours[fid0][i].id()==fid1){
					face2FaceTransport[fid0][i] = ktemp;
					//cout<<"face2FaceTransport["<<fid0<<"]["<<fid1<<"] ="<<ktemp<<endl;
					break;
				}
			}
			for(unsigned int i=0; i<face2Face_neighbours[fid1].size(); i++){
				if(face2Face_neighbours[fid1][i].id()==fid0){
					face2FaceTransport[fid1][i] = ktemp;
					break;
				}
			}
			//K[e_id] = ktemp;
		}
		else{
			vector<gmds::Node> adjacent_nodes = currentEdge.get<gmds::Node>();
			math::Point p1 = adjacent_nodes[0].getPoint();
			math::Point p2 = adjacent_nodes[1].getPoint();
			math::Vector3d v1 = math::Vector3d(p1, p2);
			v1.normalize();
			vector<gmds::Face> adjacent_faces = currentEdge.get<gmds::Face>();
			bdry_edge_normals[e_id] = v1.cross(adjacent_faces[0].normal());
		}
	}

	//edges
	for(auto n_id:m_mesh->nodes()){
		gmds::Node currentNode = m_mesh->get<gmds::Node>(n_id);
		if(m_mesh->isMarked(currentNode,m_mark_nodes_on_curve) ||
		m_mesh->isMarked(currentNode,m_mark_nodes_on_point)){
			vector<gmds::Edge> currentEdges = currentNode.get<gmds::Edge>();
			for(unsigned int i=0; i<currentEdges.size(); i++){
				if(m_mesh->isMarked(currentEdges[i],m_mark_edges_on_curve)){
					bdry_node_normals[n_id][0] = bdry_node_normals[n_id][0] + bdry_edge_normals[currentEdges[i].id()][0];
					bdry_node_normals[n_id][1] = bdry_node_normals[n_id][1] + bdry_edge_normals[currentEdges[i].id()][1];
					bdry_node_normals[n_id][2] = bdry_node_normals[n_id][2] + bdry_edge_normals[currentEdges[i].id()][2];
				}
			}
		}
	}
	//nodes
/*


  	for(auto f_id:m_mesh->faces()){
		Face current = m_mesh->get<Face>(f_id);
      	for (int j=0;j<3;j++){
		     if(j>face2Face_neighbours[f_id].size())
				face2FaceDeviation[f_id][j] = 0;
        		else{

		 		std::vector<math::Vector3d> compVectors;
	  			compVectors = (*m_faceField)[f_id].componentVectors();

			  	gmds::math::Vector3d dir0 = compVectors[0];// 1st comp vect of CF in f_id

			  	compVectors.clear();
			  	compVectors = (*m_faceField)[ face2Face_neighbours[f_id][j].id()].componentVectors();
     			gmds::math::Vector3d dir1 = compVectors[0]; //1st comp vect of CF in f1
      			gmds::math::Vector3d n0 = face_normals[f_id];
      			gmds::math::Vector3d n1 = face_normals[face2Face_neighbours[f_id][j].id()];

				///find the axis of rotation from v0 to v1
				 math::Matrix<3,3,double> rotMatrix;
				double innerProdV0V1 = (n0.normalize()).dot(n1.normalize());
  				math::Vector3d axis = n0.cross(n1);
  				axis.normalize();

  				///assemble the rotation matrix
  				double u=axis[0];
  				double v=axis[1];
  				double w=axis[2];
  				double phi=std::acos(innerProdV0V1);
  				double phicos = std::cos(phi);
  				double phisin = std::sin(phi);

  				rotMatrix.set(0,0,phicos + u*u*(1.0 - phicos));
  				rotMatrix.set(1,0, w * phisin + v*u*(1.0 - phicos));
  				rotMatrix.set(2,0, -v * phisin + w*u*(1.0 - phicos));
  				rotMatrix.set(0,1, -w * phisin + u*v*(1.0 - phicos));
  				rotMatrix.set(1,1, phicos + v*v*(1.0 - phicos));
  				rotMatrix.set(2,1, u * phisin + w*v*(1.0 - phicos));
  				rotMatrix.set(0,2, v * phisin + u*w*(1.0 - phicos));
  				rotMatrix.set(1,2, -u * phisin + v*w*(1.0 - phicos));
				rotMatrix.set(2,2, phicos + w*w*(1.0-phicos));

       			gmds::math::Vector3d dir1Rot = rotMatrix*dir1;
      			dir1Rot.normalize();

       			compVectors.clear();
	  			compVectors = (*m_faceField)[f_id].componentVectors();

       			double angle_diff = atan2(dir1Rot.dot(compVectors[1]),dir1Rot.dot(compVectors[0]));

      			double step=gmds::math::Constants::PIDIV2;
      			int missMatch=(int)std::floor((angle_diff/step)+0.5);
      			face2FaceDeviation[f_id][j] = 0;
      			if (missMatch>=0)
					 face2FaceDeviation[f_id][j] = missMatch%4;
      			else
					 face2FaceDeviation[f_id][j] = (-(3*missMatch))%4;
	   		}
      	}
	}
*/

}

/*---------------------------------------------------------------------------------------------*/

void SingularityGraphBuilder2D::createSingularityLinesShortestPaths(vector<vector<gmds::TCellID>>& newNode2NewNodeNeighbours,
                                                                    vector<vector<double>>&        face2FaceTransport,
                                                                    vector<vector<unsigned int>>&  face2FaceDeviation,
                                                                    vector<bool>&                  singOrGeomFaces)
{
	cout<<"createSingularityLinesShortestPaths"<<endl;

	bool remeshConflictingTriangles = false;

	bool visualizeAllSingLines = false;

	bool visualizeTraversedTriangles = false;

	auto tSetupInitial = Clock::now();
	std::vector<SingularityPoint*> singularity_points = m_graph.getPoints();
	unsigned int singPointNo = singularity_points.size();
	cout<<"singPointNo "<<singPointNo<<endl;
	unsigned int totalNumberOfSlots = singPointNo*5;
	//WARNING consider maximum number of variables for singularity  as equal to 5
	gmds::math::Vector3d zeroVec(0.0, 0.0, 0.0);
	gmds::math::Point zeroPoint(0.0, 0.0, 0.0);

	vector<unsigned int> originalSlotsNoPerSing(singPointNo, 0);
	vector<vector<TCellID>>  possibleTargetTriangles(singPointNo, vector<TCellID> (5, original_faces_number+2));
	vector<vector<gmds::math::Point>>  last_discr_point(singPointNo, vector<gmds::math::Point> (5, original_faces_number+2));
	vector<bool> possibleVariablesSlotsLaunched(singPointNo*5, false);
	vector<vector<bool>> nonFrozenSlots(singPointNo, vector<bool> (5,false));
	unsigned int contSource = singPointNo;

	for (unsigned int i = 0; i < singPointNo; i++) {
		SingularityPoint* pi = singularity_points[i];
		std::vector<SingularityPoint::Slot*> pi_slots = pi->getSlots();
		originalSlotsNoPerSing[i] = pi_slots.size();
		bool validSing = false;
		for (unsigned int j = 0; j < pi_slots.size(); j++){
			//slotToVarBdryCorrespondence[i][j] = tempCont;
			if (!pi_slots[j]->isLaunched){
				possibleVariablesSlotsLaunched[5*i+j] = true;
				validSing = true;
			}
//			else
//				cout<<"not valid sing i= "<<i<<" slot j= "<<j<<endl;
			if (!pi_slots[j]->isFreeze)
				nonFrozenSlots[i][j] = true;
		}
	}

	vector<vector<vector<math::Point>>> line_discretization(singPointNo, vector<vector<math::Point>>(5, vector<math::Point>() ));
	vector<vector<vector<TCellID>>>     line_triangles(singPointNo, vector<vector<TCellID>>(5, vector<TCellID>()));
	vector<vector<bool>>     addedAlready(singPointNo, vector<bool>(5, false));
	vector<vector<double>> streamlineDeviation(singPointNo, vector<double>(5, 0.0));

	vector<vector<TCellID>> 	to_cell_id(singPointNo, vector<TCellID>(5));
	vector<vector<int>> 	to_cell_dim(singPointNo, vector<int>(5));
	vector<vector<math::Vector3d>> 	to_dir(singPointNo, vector<math::Vector3d>(5));
	vector<vector<TCellID>>  targetTriangles(singPointNo, vector<TCellID>(5));

	/*//visualize bdry_edge_normals

	gmds::Mesh meshSing(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));

	for(auto e_id:m_mesh->edges()){
		gmds::Edge currentEdge = m_mesh->get<Edge>(e_id);
	  	if (m_mesh->isMarked(currentEdge, m_mark_edges_on_curve)) {
	  		vector<gmds::Node> currentNodes = currentEdge.get<Node>();
	  		gmds::math::Point myPoint0 = (currentNodes[0].getPoint() + currentNodes[1].getPoint())*0.5;
	  		gmds::Node mySing0 = meshSing.newNode(myPoint0.X(), myPoint0.Y(), myPoint0.Z());
	  		gmds::math::Point myPoint1 = myPoint0 + bdry_edge_normals[e_id];
	  		gmds::Node mySing1 = meshSing.newNode(myPoint1.X(), myPoint1.Y(), myPoint1.Z());
			meshSing.newTriangle(mySing0, mySing1, mySing1);
	  	}
   }

   gmds::IGMeshIOService ioServiceSing(&meshSing);
   gmds::VTKWriter vtkWriterSing(&ioServiceSing);
   vtkWriterSing.setCellOptions(gmds::N|gmds::F);
   vtkWriterSing.setDataOptions(gmds::N|gmds::F);
   std::string file_name = m_output_directory_name+"-bdry_edge_normals.vtk";
   vtkWriterSing.write(file_name);*/

	vector<TCellID> bdryFaces;

	for(unsigned int tt=0; tt< original_faces_number; tt++){
		if(is_bdry_face[tt]){
			bdryFaces.push_back(tt);
		}
	}
	/*std::string file_name_bdry_faces = "bdryFaces.vtk";
	writeTestMeshTrianglesIds(bdryFaces,file_name_bdry_faces);*/

	for(unsigned int i=0; i<singPointNo; i++){
		SingularityPoint* pi = singularity_points[i];
		std::vector<SingularityPoint::Slot*> pi_slots = pi->getSlots();

		for(unsigned j = 0; j<originalSlotsNoPerSing[i]; j++){
			line_discretization[i][j].push_back(pi_slots[j]->from_point->getLocation());
			line_discretization[i][j].push_back(pi_slots[j]->location);
			to_cell_id[i][j]  = pi_slots[j]->starting_cell_id ;
			to_cell_dim[i][j] = pi_slots[j]->starting_cell_dim;
			to_dir[i][j] = pi_slots[j]->direction;
		}
	}

	cout<<"------------------------------------------------------------------------"<<endl;

	double searchStep = mean_edge_length;

	bool AOnEdge0, AOnEdge1, AOnEdge2;
	vector<double> lambdas(3);

	vector<pair<unsigned int, unsigned int>> fixedVariables;
	vector<pair<unsigned int, unsigned int>> fixedVariablesOpt;
	unsigned int maxNoSing = 5;

	vector<pair<vector<unsigned int>, vector<unsigned int >>> isTraversedFaceCompVector(original_faces_number);
	vector<pair<vector<unsigned int>, vector<unsigned int >>> isTraversedNodeCompVector(original_nodes_number);
	vector<pair<vector<unsigned int>, vector<unsigned int >>> isTraversedFaceCompVector_SegmentPathCont(original_faces_number);
	vector<pair<vector<unsigned int>, vector<unsigned int >>> isTraversedFaceCompVector_SegmentPathCode(original_faces_number);

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	/*for each slot pi_slots[j] of singularity point singularity_points[i], advance along the slots directions (using Runge-Kutta 4) until the landing faces are different.*/
	for (unsigned int i = 0; i < singPointNo; i++) {
		SingularityPoint* pi = singularity_points[i];

		std::vector<SingularityPoint::Slot*> pi_slots = pi->getSlots();

		bool needToLaunch = false;
		for(unsigned int j=0; j<pi_slots.size(); j++){
			if(!pi_slots[j]->isLaunched)
				needToLaunch = true;
		}
		if(needToLaunch){
			gmds::math::Point intersectionPoint;
			double intersectionParam;
			bool diffFaces = false;
			bool firstPass = true;
			while(!diffFaces){
				for (unsigned int j = 0; j < pi_slots.size(); j++){
					SingularityPoint::Slot* current_slot = pi_slots[j];
					if((current_slot->isFreeze!=true)&&(!addedAlready[i][j])) {
						SingularityPoint       *to_sing_pnt  = 0;
						SingularityPoint::Slot *to_slot = 0;

						bool find_end_bdry = false;
						bool end_on_free_slot = false;
						bool find_end = false;

						if(firstPass){
							Edge currentEdge = m_mesh->get<Edge>(to_cell_id[i][j]);
							std::vector<TCellID> adj_faces = currentEdge.getIDs<Face>();
							if(adj_faces.size()==1){
								if(m_tool.isPntInTri(current_slot->from_point->getLocation(), m_mesh->get<Face>(adj_faces[0]), AOnEdge0, AOnEdge1, AOnEdge2, lambdas[0], lambdas[1]))
									line_triangles[i][j].push_back(adj_faces[0]);
							}
							else{
								if(m_tool.isPntInTri(current_slot->from_point->getLocation(), m_mesh->get<Face>(adj_faces[0]), AOnEdge0, AOnEdge1, AOnEdge2, lambdas[0], lambdas[1])){
									line_triangles[i][j].push_back(adj_faces[0]);
									line_triangles[i][j].push_back(adj_faces[1]);
								}
								else{
									line_triangles[i][j].push_back(adj_faces[1]);
									line_triangles[i][j].push_back(adj_faces[0]);
								}
							}
						}

						growLineRK4(current_slot->from_point,
						            current_slot,
						            to_sing_pnt,
						            to_slot,
						            line_discretization[i][j].back(),
						            to_dir[i][j],
						            line_discretization[i][j],
						            line_triangles[i][j],
						            to_cell_dim[i][j],
						            to_cell_id[i][j],
						            streamlineDeviation[i][j],
						            searchStep,
						            find_end_bdry,
						            end_on_free_slot,
						            find_end);

						possibleTargetTriangles[i][j] = line_triangles[i][j].back();
						last_discr_point[i][j] = line_discretization[i][j].back();

						if(find_end_bdry){
							//======================================================================
							// CASE 1 - We finish on the boundary. A geometric point must be created
							// The cell defined by (start_cell_dim[j], start_cell_id[j]) is on the boundary.
							// We have to create a geometric singularity point.
							// however we could do this at the end of the entire iteration(for the current searchStep)
							//line creation
							cout<<"!!!find_end_bdry for current_slot->from_point "<<current_slot->from_point->getLocation().X()<<" "<<current_slot->from_point->getLocation().Y()<<
							" and with dir "<<current_slot->direction[0]<<" "<<current_slot->direction[1]<<endl;

							if(m_build_geometric_singularities){
								/*WARNING TODO: we must check if the boundary intersection point is very
								 *close to a geometric singularity and if the streamline is aligned with one
								 * of the geometric singularity's slots. If so, that won't be a boundary
								 * singularity line, but a singularity line between two slots of two singularities. */

								for(auto sing_it:m_singularities_3){
									if(sing_it.id()==line_triangles[i][j].back()){

										cout<<"sing_it.id()==line_triangles[i][j][k].back()= "<<sing_it.id()<<endl;
									}
								}
							}
							SurfaceSingularityLine* surf_line = m_graph.newSurfaceLine();
							int sepNumberTmp = m_graph.getNbLines();
							surf_line->setNumber(sepNumberTmp);
							//connect line to initial singularity point
							SingularityPoint* from_sing_pnt = current_slot->from_point;
							surf_line->addSingularityPoint(from_sing_pnt);
							surf_line->addDiscretizationPoint(from_sing_pnt->getLocation());
							current_slot->line = surf_line;
							math::Vector3d firstDir = math::Vector3d(from_sing_pnt->getLocation(), line_discretization[i][j][1]);

							current_slot->line_direction = firstDir;
							current_slot->isLaunched = true;
							if(to_slot!=0){
								streamlineDeviation[i][j] = (streamlineDeviation[i][j] + fabs(1 - to_slot->direction.dot(math::Vector3d(line_discretization[i][j].back(),line_discretization[i][j][line_discretization[i][j].size()-2]))));
								streamlineDeviation[i][j] = streamlineDeviation[i][j]/(line_discretization[i][j].size()+1);
								current_slot->lineDeviation = streamlineDeviation[i][j];
							}
							//Insertion of line points
							for(unsigned int k=0; k<line_discretization[i][j].size(); k++) {
								surf_line->addDiscretizationPoint(line_discretization[i][j][k]);
							}

							for(unsigned int k=0; k<line_triangles[i][j].size(); k++) {
								surf_line->addTraversedFace(line_triangles[i][j][k]);
							}

							SingularityPoint* geom_pnt;
							SingularityPoint::Slot* incoming_slot;

							createGeometricSingularityPoint(line_discretization[i][j].back(),      // the last point added
							                                to_dir[i][j],      // the direction we come from
							                                to_cell_dim[i][j], // the dim. of the cell
							                                to_cell_id[i][j],  // the id of the cell
							                                geom_pnt,       // the created point
							                                incoming_slot); // and the slot
							current_slot->lineDeviation = streamlineDeviation[i][j];
							surf_line->addDiscretizationPoint(geom_pnt->getLocation());
							surf_line->addSingularityPoint(geom_pnt);
							current_slot->isFreeze = true;
							current_slot->isLaunched = true;
							possibleTargetTriangles[i][j] = original_faces_number + 1;
							addedAlready[i][j] = true;
							contSource = 5*i + j;
							for(unsigned int k0=1; k0<line_discretization[i][j].size(); k0++){
								math::Segment from_seg(line_discretization[i][j][k0-1], line_discretization[i][j][k0]);
								math::Vector3d tempDir(line_discretization[i][j][k0-1], line_discretization[i][j][k0]);
								for(unsigned int k=1; k<line_triangles[i][j].size(); k++) {
									vector<gmds::Edge> temp_edges = m_mesh->get<gmds::Face>(line_triangles[i][j][k]).get<gmds::Edge>();
									for(unsigned int j1=0; j1<temp_edges.size(); j1++){
										vector<gmds::Node> currentNodes =  temp_edges[j1].get<gmds::Node>();
										math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
										if(from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon)){
											math::Cross2D crossV = triangle_centers_cross[line_triangles[i][j][k]];
											unsigned int contMatrix = contSource*(totalNumberOfSlots + 1) + totalNumberOfSlots;
											unsigned int compVId = crossV.closestComponentVectorAsIndex(tempDir)%2;
											if(compVId==0){
												isTraversedFaceCompVector[line_triangles[i][j][k]].first.push_back(contMatrix);
												isTraversedFaceCompVector_SegmentPathCode[line_triangles[i][j][k]].first.push_back(0);
												isTraversedFaceCompVector_SegmentPathCont[line_triangles[i][j][k]].first.push_back(k0-1);
											}
											else{
												isTraversedFaceCompVector[line_triangles[i][j][k]].second.push_back(contMatrix);
												isTraversedFaceCompVector_SegmentPathCode[line_triangles[i][j][k]].second.push_back(0);
												isTraversedFaceCompVector_SegmentPathCont[line_triangles[i][j][k]].second.push_back(k0-1);
											}
										}
									}
								}
							}

							cout<<"addedAlready["<<i<<"]["<<j<<"]"<<endl;
							vector<math::Point> Newline_discretization = surf_line->getDiscretizationPoints();
							gmds::Mesh meshSing(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
							gmds::Node mySing0 = meshSing.newNode(Newline_discretization[0].X(), Newline_discretization[0].Y(), Newline_discretization[0].Z());
							for(unsigned int j5=1; j5<Newline_discretization.size(); j5++) {
								gmds::Node mySing = meshSing.newNode(Newline_discretization[j5].X(), Newline_discretization[j5].Y(), Newline_discretization[j5].Z());
								meshSing.newTriangle(mySing0, mySing, mySing);
								mySing0 = mySing;
							}
							gmds::IGMeshIOService ioServiceSing(&meshSing);
							gmds::VTKWriter vtkWriterSing(&ioServiceSing);
							vtkWriterSing.setCellOptions(gmds::N|gmds::F);
							vtkWriterSing.setDataOptions(gmds::N|gmds::F);
							std::string file_name = m_output_directory_name+"-111discretizationPointsFoundBdry_"+std::to_string(5*i+j)+".vtk";
							vtkWriterSing.write(file_name);
							//file_name = "bdyLine";
							//writeOutput(file_name);
							cout<<"!!!has written "<<file_name<<endl;
							fixedVariables.push_back(make_pair(i*5+j, i*5+j));
							if(visualizeTraversedTriangles){
								std::string file_name_trav_tri = "trav_tri_"+std::to_string(5*i+j)+".vtk";

								writeTestMeshTrianglesIds(line_triangles[i][j], file_name_trav_tri);
							}
						}
						else{// not find_end_bdry
							if(end_on_free_slot){
								if(withGlobalComments){
									cout<<"!!!end_on_free_slot; current_slot->location "<<current_slot->location<<endl;
									cout<<"line disret "<<line_discretization[i][j][line_discretization[i][j].size()-2];
									cout<<" "<<line_discretization[i][j][line_discretization[i][j].size()-1]<<endl;
								}
								SurfaceSingularityLine* surf_line = m_graph.newSurfaceLine();
								int sepNumberTmp = m_graph.getNbLines();
								surf_line->setNumber(sepNumberTmp);

								SingularityPoint* from_sing_pnt = current_slot->from_point;
								surf_line->addSingularityPoint(from_sing_pnt);
								surf_line->addDiscretizationPoint(from_sing_pnt->getLocation());

								current_slot->line = surf_line;

								math::Vector3d firstDir = math::Vector3d(from_sing_pnt->getLocation(), line_discretization[i][j][1]);
								current_slot->line_direction = firstDir;
								current_slot->isLaunched = true;
								current_slot->lineDeviation = streamlineDeviation[i][j];

								for(unsigned int k=0; k<line_discretization[i][j].size();k++) {
									surf_line->addDiscretizationPoint(line_discretization[i][j][k]);
								}

								for(unsigned int k=0; k<line_triangles[i][j].size(); k++) {
									surf_line->addTraversedFace(line_triangles[i][j][k]);
								}

								unsigned int temp_i, temp_j;

								SingularityPoint* pi2 = to_slot->from_point;
								std::vector<SingularityPoint::Slot*> pi2_slots = pi2->getSlots();
								for (unsigned int t1 = 0; t1 < singPointNo; t1++) {
									if(singularity_points[t1]==pi2){
										temp_i = t1;
										for (unsigned int t2 = 0; t2 < originalSlotsNoPerSing[t1]; t2++){
											if(pi2_slots[t2]==to_slot)
												temp_j = t2;
										}
									}
								}

								to_slot->isLaunched = true;
								current_slot->isFreeze = true;
								to_slot->isFreeze = true;

								to_slot->line = surf_line;
								to_slot->line_direction = to_dir[i][j];

								to_slot->lineDeviation = streamlineDeviation[i][j];

								surf_line->addSingularityPoint(to_sing_pnt);
								surf_line->addDiscretizationPoint(to_sing_pnt->getLocation());

								possibleTargetTriangles[i][j] = original_faces_number+1;//no+1
								possibleTargetTriangles[temp_i][temp_j] = original_faces_number+1;
								addedAlready[i][j] = true;
								addedAlready[temp_i][temp_j] = true;

								for(unsigned int k0=1; k0<line_discretization[i][j].size(); k0++){
									math::Segment from_seg(line_discretization[i][j][k0-1], line_discretization[i][j][k0]);
									math::Vector3d tempDir(line_discretization[i][j][k0-1], line_discretization[i][j][k0]);
									for(unsigned int k=1; k<line_triangles[i][j].size(); k++) {
										vector<gmds::Edge> temp_edges = m_mesh->get<gmds::Face>(line_triangles[i][j][k]).get<gmds::Edge>();
										for(unsigned int j1=0; j1<temp_edges.size(); j1++){
											vector<gmds::Node> currentNodes =  temp_edges[j1].get<gmds::Node>();
											math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
											if(from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon)){
												math::Cross2D crossV = triangle_centers_cross[line_triangles[i][j][k]];
												unsigned int contMatrix = contSource*(totalNumberOfSlots + 1) + 5*temp_i+temp_j;
												unsigned int compVId = crossV.closestComponentVectorAsIndex(tempDir)%2;
												if(compVId==0){
													isTraversedFaceCompVector[line_triangles[i][j][k]].first.push_back(contMatrix);
													isTraversedFaceCompVector_SegmentPathCode[line_triangles[i][j][k]].first.push_back(0);
													isTraversedFaceCompVector_SegmentPathCont[line_triangles[i][j][k]].first.push_back(k0-1);
												}
												else{
													isTraversedFaceCompVector[line_triangles[i][j][k]].second.push_back(contMatrix);
													isTraversedFaceCompVector_SegmentPathCode[line_triangles[i][j][k]].second.push_back(0);
													isTraversedFaceCompVector_SegmentPathCont[line_triangles[i][j][k]].second.push_back(k0-1);
												}
											}
										}
								 	}
								}

								cout<<"addedAlready["<<i<<"]["<<j<<"]"<<endl;
								cout<<"addedAlready["<<temp_i<<"]["<<temp_j<<"]"<<endl;
								vector<math::Point> Newline_discretization = surf_line->getDiscretizationPoints();
								gmds::Mesh meshSing(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
								gmds::Node mySing0 = meshSing.newNode(Newline_discretization[0].X(), Newline_discretization[0].Y(), Newline_discretization[0].Z());
			 					for(unsigned int j5=1; j5<Newline_discretization.size(); j5++) {
									gmds::Node mySing = meshSing.newNode(Newline_discretization[j5].X(), Newline_discretization[j5].Y(), Newline_discretization[j5].Z());
									meshSing.newTriangle(mySing0, mySing, mySing);
									mySing0 = mySing;
								}
								gmds::IGMeshIOService ioServiceSing(&meshSing);
								gmds::VTKWriter vtkWriterSing(&ioServiceSing);
								vtkWriterSing.setCellOptions(gmds::N|gmds::F);
								vtkWriterSing.setDataOptions(gmds::N|gmds::F);
								std::string file_name = m_output_directory_name+"-discretizationPointsFreeSlot_"+std::to_string(i)+"_"+std::to_string(j)+"--"+std::to_string(temp_i)+"_"+std::to_string(temp_j)+".vtk";
								vtkWriterSing.write(file_name);
								cout<<"!!!has written "<<file_name<<endl;
								fixedVariables.push_back(make_pair(i*5+j, temp_i*5 + temp_j));

								//file_name = "bdyLine";
								//writeOutput(file_name);

								if(visualizeTraversedTriangles){
									std::string file_name_trav_tri = "trav_tri_"+std::to_string(5*i+j)+"--"+std::to_string(5*temp_i+temp_j)+".vtk";
									writeTestMeshTrianglesIds(line_triangles[i][j], file_name_trav_tri);
								}

							}
							else{/*we reached a slot that has already been connected; highly improbable*/
								if(withGlobalComments){
									cout<<"not end_on_free_slot; possibleTargetTriangles[i][j] "<<possibleTargetTriangles[i][j]<<endl;
									cout<<"last_discr_point[i][j] "<<last_discr_point[i][j]<<endl;
								}
								//connect if within threshold distance or check neigh triangles
								math::Vector3d sourceDirTemp  = to_dir[i][j];
								sourceDirTemp.normalize();
								math::Vector3d sourceCrossTemp = triangle_centers_cross[line_triangles[i][j].back()].closestComponentVector(sourceDirTemp);

								for(unsigned int i3=0; (i3!=i) && i3<singPointNo && (!addedAlready[i][j]); i3++){
									SingularityPoint* pi3 = singularity_points[i3];
									std::vector<SingularityPoint::Slot*> pi3_slots = pi3->getSlots();
									for(unsigned int j3=0; j3<pi3_slots.size() && (!addedAlready[i][j]); j3++){
										if((!addedAlready[i3][j3])&&(!line_triangles[i3][j3].empty())){
											for(unsigned int k3=0; k3<face2Face_neighbours_by_verts[line_triangles[i][j].back()].size(); k3++){
												if(line_triangles[i3][j3].back() == face2Face_neighbours_by_verts[line_triangles[i][j].back()][k3].id()){

													math::Vector3d targetDirTemp  = to_dir[i3][j3];
													targetDirTemp.normalize();
													math::Vector3d targetCrossTemp = triangle_centers_cross[line_triangles[i3][j3].back()].closestComponentVector(targetDirTemp);
													targetCrossTemp = -targetCrossTemp;
													math::Vector3d tempConnection(line_discretization[i][j].back(), line_discretization[i3][j3].back());
													//cout<<"tempConnection.angle(sourceCrossTemp) "<<tempConnection.angle(sourceCrossTemp)<<endl;
													//cout<<"tempConnection.angle(targetCrossTemp) "<<tempConnection.angle(targetCrossTemp)<<endl;
													if((targetCrossTemp.angle(sourceCrossTemp)<gmds::math::Constants::PIDIV4)&&(tempConnection.angle(sourceCrossTemp)<gmds::math::Constants::PIDIV4)){

														SurfaceSingularityLine* surf_line = m_graph.newSurfaceLine();
														int sepNumberTmp = m_graph.getNbLines();
														surf_line->setNumber(sepNumberTmp);
														SingularityPoint* from_sing_pnt = current_slot->from_point;
														surf_line->addSingularityPoint(from_sing_pnt);
														surf_line->addDiscretizationPoint(from_sing_pnt->getLocation());

														current_slot->line = surf_line;
														math::Vector3d firstDir = math::Vector3d(from_sing_pnt->getLocation(), line_discretization[i][j][1]);

														current_slot->line_direction = firstDir;
														current_slot->lineDeviation = streamlineDeviation[i][j];

														fixedVariables.push_back(make_pair(i*5+j, i3*5 + j3));

														math::Segment from_seg_temp(line_discretization[i][j].back(), line_discretization[i3][j3].back());

														gmds::TCellID prevFace;
														gmds::TCellID nextFaceTemp = line_triangles[i3][j3].back();
														gmds::TCellID currentFaceTemp = line_triangles[i][j].back();

														if(line_triangles[i][j].size()>1)
															prevFace = line_triangles[i][j][line_triangles[i][j].size()-2];
															gmds::math::Point intersectionPoint;

															double intersectionParam;
															bool moved2nextFace = false;
															while(currentFaceTemp!=nextFaceTemp){
																TCellID tempFace;
																vector<gmds::Edge> currentEdges = (m_mesh->get<gmds::Face>(currentFaceTemp)).get<gmds::Edge>();
																for(unsigned int j7=0; j7<3; j7++){
																vector<gmds::Node> currentNodes =  currentEdges[j7].get<gmds::Node>();
																math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());

																if(from_seg_temp.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon)){

																	if((intersectionParam<=temp_epsilon)||(intersectionParam>=1.0-temp_epsilon)){
																		//PASSES through a node
																		gmds::Node passedNode;
																	if((intersectionParam<=temp_epsilon))
																			passedNode = currentNodes[0];
																		else
																			passedNode = currentNodes[1];

																		vector<gmds::Face> nodeAdjFaces = passedNode.get<gmds::Face>();
																		for(unsigned int j9=0; j9<nodeAdjFaces.size(); j9++){
																			vector<gmds::Edge> tempEdges = nodeAdjFaces[j9].get<gmds::Edge>();
																			for(unsigned int j4=0; j4<3;j4++){
																				vector<gmds::Node> tempEdgeNodes = tempEdges[j4].get<gmds::Node>();
																				if((tempEdgeNodes[0]!=passedNode)&&(tempEdgeNodes[1]!=passedNode)){
																					math::Segment tempOppSeg(tempEdgeNodes[0].getPoint(), tempEdgeNodes[1].getPoint());

																					if(from_seg_temp.SecondMetIntersect2D(tempOppSeg, intersectionPoint, intersectionParam, temp_epsilon)){

																						prevFace = currentFaceTemp;
																						currentFaceTemp = nodeAdjFaces[j9].id();
																						moved2nextFace = true;
																						break;
																					}
																				}
																			}
																		}
																	}
																	else{
																		vector<gmds::Face> adj_faces = currentEdges[j7].get<gmds::Face>();

																		if(adj_faces[0].id()!=currentFaceTemp){
																			tempFace = adj_faces[0].id();
																			if(adj_faces[0].id()!=prevFace){
																				prevFace = currentFaceTemp;
																				currentFaceTemp = adj_faces[0].id();
																				line_triangles[i][j].push_back(currentFaceTemp);
																				moved2nextFace = true;
																				break;
																			}
																		}
																		if(adj_faces.size()>1){
																			if(adj_faces[1].id()!=currentFaceTemp){
																				tempFace = adj_faces[1].id();
																				if(adj_faces[1].id()!=prevFace){
																					prevFace = currentFaceTemp;
																					currentFaceTemp = adj_faces[1].id();
																					line_triangles[i][j].push_back(currentFaceTemp);
																					moved2nextFace = true;
																					break;
																				}
																			}
																		}
																	}
																}
															}
															if(!moved2nextFace){
																prevFace = currentFaceTemp;
																currentFaceTemp = tempFace;
																line_triangles[i][j].push_back(currentFaceTemp);
															}
														}

														reverse(line_triangles[i3][j3].begin(),line_triangles[i3][j3].end());
														line_triangles[i][j].insert(line_triangles[i][j].end(),line_triangles[i3][j3].begin(), line_triangles[i3][j3].end());
														reverse(line_triangles[i3][j3].begin(),line_triangles[i3][j3].end());
														reverse(line_discretization[i3][j3].begin(), line_discretization[i3][j3].end());
														line_discretization[i][j].insert(line_discretization[i][j].end(),line_discretization[i3][j3].begin(), line_discretization[i3][j3].end());
														reverse(line_discretization[i3][j3].begin(), line_discretization[i3][j3].end());

														for(unsigned int k=0; k<line_discretization[i][j].size();k++) {
															surf_line->addDiscretizationPoint(line_discretization[i][j][k]);
														}

														for(unsigned int k=0; k<line_triangles[i][j].size(); k++) {
															surf_line->addTraversedFace(line_triangles[i][j][k]);
														}
														current_slot->isLaunched = true;
														pi3_slots[j3]->isLaunched = true;
														current_slot->isFreeze = true;
														pi3_slots[j3]->isFreeze = true;
														pi3_slots[j3]->line = surf_line;
														pi3_slots[j3]->line_direction = to_dir[i3][j3];
														pi3_slots[j3]->lineDeviation = streamlineDeviation[i][j];
														surf_line->addSingularityPoint(pi3);
														surf_line->addDiscretizationPoint(pi3->getLocation());
														possibleTargetTriangles[i][j] = original_faces_number+1;
														possibleTargetTriangles[i3][j3] = original_faces_number+1;
														addedAlready[i][j] = true;
														addedAlready[i3][j3] = true;
														/* //unsigned int contMatrix = 5*i+j*(totalNumberOfSlots + 1) + 5*i3+j3;
														unsigned int compVId = crossV.closestComponentVectorAsIndex(tempDir)%2;
														if(compVId==0)
															isTraversedFaceCompVector[currentFace_id].first.push_back(contMatrix);
														else
															isTraversedFaceCompVector[currentFace_id].first.push_back(contMatrix);
														*/
														cout<<"addedAlready["<<i<<"]["<<j<<"]"<<endl;
														cout<<"addedAlready["<<i3<<"]["<<j3<<"]"<<endl;
														vector<math::Point> Newline_discretization = surf_line->getDiscretizationPoints();
														gmds::Mesh meshSing(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
														gmds::Node mySing0 = meshSing.newNode(Newline_discretization[0].X(), Newline_discretization[0].Y(), Newline_discretization[0].Z());
														for(unsigned int j5=1; j5<Newline_discretization.size(); j5++) {
															gmds::Node mySing = meshSing.newNode(Newline_discretization[j5].X(), Newline_discretization[j5].Y(), Newline_discretization[j5].Z());
															meshSing.newTriangle(mySing0, mySing, mySing);
															mySing0 = mySing;
														}
														gmds::IGMeshIOService ioServiceSing(&meshSing);
														gmds::VTKWriter vtkWriterSing(&ioServiceSing);
														vtkWriterSing.setCellOptions(gmds::N|gmds::F);
														vtkWriterSing.setDataOptions(gmds::N|gmds::F);
														std::string file_name = m_output_directory_name+"-discretizationPointsNoFreeSlotThresoldDistSlot_"+std::to_string(5*i+j)+" and "+std::to_string(5*i3+j3)+".vtk";
														vtkWriterSing.write(file_name);
														cout<<"!!!has written "<<file_name<<endl;
														//file_name = "bdyLine";
     													//writeOutput(file_name);
														break;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}

				bool sameTris = false;

				for (unsigned int j1 = 0; j1 < originalSlotsNoPerSing[i]-1; j1++){
					if((nonFrozenSlots[i][j1]) && (!addedAlready[i][j1])){
						for (unsigned int j2 = j1+1; j2 < originalSlotsNoPerSing[i]; j2++){
							if((!addedAlready[i][j2])&&(nonFrozenSlots[i][j2])){
								if(line_triangles[i][j1].back()==line_triangles[i][j2].back()){
									sameTris = true;
									break;
								}
								else{
									vector<gmds::Node> tempNodes1 = (m_mesh->get<gmds::Face>(line_triangles[i][j1].back())).get<gmds::Node>();
									vector<gmds::Node> tempNodes2 = (m_mesh->get<gmds::Face>(line_triangles[i][j2].back())).get<gmds::Node>();
									for(unsigned int j3 = 0; j3<3; j3++){
										for(unsigned int j4 = 0; j4<3; j4++){
											if(tempNodes1[j3]==tempNodes2[j4]){
												sameTris = true;
												break;
											}
										}
									}
								}
							}
						}
					}
				}

				if(!sameTris)
					diffFaces = true;
				firstPass = false;

			}
		}
	}

//	for(unsigned int i=0; i<fixedVariables.size(); i++){
//		cout<<"fixedVariables["<<i<<"] "<<fixedVariables[i].first<<" (sing "<<(int)(fixedVariables[i].first/5)<<") slot (, "<<fmod(fixedVariables[i].first,5)<<")"<<endl;
//		cout<<"towards fixedVariables["<<i<<"] "<<fixedVariables[i].second<<" (sing "<<(int)(fixedVariables[i].second/5)<<") slot (, "<<fmod(fixedVariables[i].second,5)<<")"<<endl;
//	}

	vector<unsigned int> faceNo2Cont(original_faces_number+2);
	singPointNo = possibleTargetTriangles.size();
	for (unsigned int i = 0; i < singPointNo; i++) {
		for (unsigned int j = 0; j < originalSlotsNoPerSing[i]; j++){	//bef slotsNoPerSing
//			cout<<"possibleTargetTriangles["<<i<<"]["<<j<<"]= "<<possibleTargetTriangles[i][j]<<endl;
			if(possibleTargetTriangles[i][j]<=original_faces_number+1){
				faceNo2Cont[possibleTargetTriangles[i][j]] = 5*i+j;
			}
		}
	}

	unsigned int totalNumberOfVariables = totalNumberOfSlots + 1;

	double maxDist = 1000000.0;
	vector<vector<double>> distances(totalNumberOfSlots, vector<double>(totalNumberOfVariables, maxDist));

	vector<vector<vector<gmds::TCellID>>> traversedSPTris(totalNumberOfSlots, vector<vector<gmds::TCellID>>(totalNumberOfVariables, vector<gmds::TCellID>()));
	if(withGlobalComments){
		cout<<"original_faces_number "<<original_faces_number<<endl;
		cout<<"traversedSPTris.size() "<<traversedSPTris.size()<<endl;
		cout<<"totalNumberOfSlots "<<totalNumberOfSlots<<endl;
		cout<<"distances.size() "<<distances.size()<<endl;
	}

	vector<vector<vector<gmds::TCellID>>> finalPaths(totalNumberOfSlots, vector<vector<gmds::TCellID>>(totalNumberOfVariables, vector<gmds::TCellID>()));
	vector<unsigned int> paths2(totalNumberOfSlots);

	vector<vector<pair<vector<gmds::math::Point>, vector<gmds::math::Point>>>> pointPaths(totalNumberOfSlots, vector<pair<vector<gmds::math::Point>, vector<gmds::math::Point>>>(totalNumberOfVariables));
	vector<gmds::math::Point> endPoint(totalNumberOfSlots);
	vector<vector<vector<gmds::math::Point>>> finalCenterPoints(totalNumberOfSlots, vector<vector<gmds::math::Point>>(totalNumberOfSlots, vector<gmds::math::Point>()));

	contSource = 0;
	vector<bool> forbiddenFaces(original_faces_number, false); //- all slot faces(+singularity faces) except source and target
	vector<bool> forbiddenFacesJustSing(original_faces_number, false);

	for (unsigned int i = 0; i < singPointNo; i++) {
		for (unsigned int j = 0; j < originalSlotsNoPerSing[i]; j++){	 //bef slotsNoPerSing
			if(possibleTargetTriangles[i][j]<original_faces_number+1){
				forbiddenFaces[possibleTargetTriangles[i][j]] = true;

				for(unsigned int t=0; t<line_triangles[i][j].size();t++){
					forbiddenFaces[line_triangles[i][j][t]] = true;
				}
				cout<<endl;
			}
		}
	}

	for(auto sing_it:m_singularities_3){
		forbiddenFaces[sing_it.id()] = true;
		forbiddenFacesJustSing[sing_it.id()] = true;
		singOrGeomFaces[sing_it.id()] = true;
	}
	for(auto sing_it:m_singularities_5){
		forbiddenFaces[sing_it.id()] = true;
		forbiddenFacesJustSing[sing_it.id()] = true;
		singOrGeomFaces[sing_it.id()] = true;
	}

	bool walkThroughNodes = false;

	bool targetBdry = false;

	vector<unsigned int> contSourceToSing;

	/////////////////////////////////////////////////////////////////////////////////////////////
	/*vector<CurveSingularityLine* >   curve_Sing_lines = m_graph.getCurveLines();

	for(unsigned int i=0; i<curve_Sing_lines.size(); i++){
		vector<gmds::math::Point> discrPts = curve_Sing_lines[i]->getDiscretizationPoints();
		if(discrPts.size()>0){
		  	gmds::Mesh meshSing(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
			gmds::Node firstPt = meshSing.newNode(discrPts[0].X(), discrPts[0].Y(), discrPts[0].Z());
			for(unsigned int j=1; j<discrPts.size(); j++){
	 			 gmds::Node nextPt = meshSing.newNode(discrPts[j].X(), discrPts[j].Y(), discrPts[j].Z());
	  			meshSing.newTriangle(firstPt, nextPt, nextPt);
				firstPt = nextPt;
			}

			gmds::IGMeshIOService ioServiceSing(&meshSing);
			gmds::VTKWriter vtkWriterSing(&ioServiceSing);
			vtkWriterSing.setCellOptions(gmds::N|gmds::F);
			vtkWriterSing.setDataOptions(gmds::N|gmds::F);

			std::stringstream vtk_fileSing;
			vtk_fileSing <<m_output_directory_name<<"-curve_Sing_lines" <<i<<".vtk";
			vtkWriterSing.write(vtk_fileSing.str());
		}
	}
	vector<SurfaceSingularityLine* >   surf_Sing_lines = m_graph.getSurfaceLines();

	for(unsigned int i=0; i<surf_Sing_lines.size(); i++){

		vector<gmds::math::Point> discrPts = surf_Sing_lines[i]->getDiscretizationPoints();
		if(discrPts.size()>0){
		  	gmds::Mesh meshSing1(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
			gmds::Node firstPt = meshSing1.newNode(discrPts[0].X(), discrPts[0].Y(), discrPts[0].Z());
			for(unsigned int j=1; j<discrPts.size(); j++){
	 			 gmds::Node nextPt = meshSing1.newNode(discrPts[j].X(), discrPts[j].Y(), discrPts[j].Z());
	  			meshSing1.newTriangle(firstPt, nextPt, nextPt);
				firstPt = nextPt;
			}

			gmds::IGMeshIOService ioServiceSing(&meshSing1);
			gmds::VTKWriter vtkWriterSing(&ioServiceSing);
			vtkWriterSing.setCellOptions(gmds::N|gmds::F);
			vtkWriterSing.setDataOptions(gmds::N|gmds::F);

			std::stringstream vtk_fileSing;
			vtk_fileSing <<m_output_directory_name<<"-surf_Sing_lines" <<i<<".vtk";
			vtkWriterSing.write(vtk_fileSing.str());
		}
	}
	*/

	vector<pair<unsigned int, unsigned int>> IllegalCross;
	vector<pair<unsigned int, unsigned int>> contSourceToSingularity(totalNumberOfSlots*totalNumberOfVariables);

//	auto tSetupFinal = Clock::now();
//	cout << "setup "
//	<< std::chrono::duration_cast<std::chrono::milliseconds>(tSetupFinal - tSetupInitial).count()
//	<< " milliseconds" << std::endl;
	/////////////////////////////////////////////////////////////////////////////////////////////
	auto t0 = Clock::now();
	double timerTravTri = 0.0;
	fixedVariablesOpt.resize(fixedVariables.size());

	double intersectionParam;

	gmds::TCellID lastVisTri = original_faces_number;
	vector<vector<int>>  to_cell_dim_sp(to_cell_dim);
	vector<vector<TCellID>> 	to_cell_id_sp(to_cell_id);

	vector<TCellID> targetFaceOpt(5*singPointNo, original_faces_number);
	vector<gmds::math::Point> targetPoints(5*singPointNo);
	vector<pair<gmds::math::Vector3d, gmds::math::Vector3d>> alltargetDirCross(5*singPointNo);
	for (unsigned int i = 0; i < singPointNo; i++) {
		for (unsigned int j = 0; j < originalSlotsNoPerSing[i]; j++){
			targetFaceOpt[5*i+j] = possibleTargetTriangles[i][j];
			targetPoints[5*i+j] = line_discretization[i][j].back();

			if(possibleTargetTriangles[i][j]<original_faces_number+1){
				gmds::math::Vector3d prevDir = to_dir[i][j];
				prevDir.normalize();
				gmds::math::Vector3d prevCross = triangle_centers_cross[targetFaceOpt[5*i+j]].closestComponentVector(prevDir);
				alltargetDirCross[5*i+j].first = -prevDir;
				alltargetDirCross[5*i+j].second = -prevCross;
			}
		}
	}
	for (unsigned int i = 0; i < singPointNo; i++) {
		for (unsigned int j = 0; j < originalSlotsNoPerSing[i]; j++){

			if(possibleTargetTriangles[i][j]<original_faces_number+1){
				gmds::TCellID source_id = possibleTargetTriangles[i][j];
				gmds::Face source = m_mesh->get<Face>(source_id);

				vector<unsigned int> tempTargetFaces;
				vector<gmds::math::Vector3d> tempTargetDir;
				gmds::math::Vector3d prevDir = to_dir[i][j];
				prevDir.normalize();
				vector<pair<double, unsigned int>> min_distance;
				vector<int> previous, previousBdry;
				int found;
				vector<vector<double>> face2FaceMatch;
				pair<gmds::math::Vector3d, gmds::math::Vector3d> prevDirCross;
				gmds::math::Vector3d prevCross;

				prevDirCross.first = prevDir;
				prevDirCross.second = prevDir;
				vector<pair<gmds::math::Vector3d, gmds::math::Vector3d>> targetDirCross(1);

				contSource = 5*i+j;

				bool detectedBdry = false;
				vector<math::Point> tempBdry;
				if(is_bdry_face[source_id]){
					math::Point intersectionPnt;
					double intersectionParam;
					math::Segment from_ray(line_discretization[i][j].back(), line_discretization[i][j].back() + prevDirCross.second);

					vector<gmds::Edge> currentEdges = source.get<gmds::Edge>();
					bool hasBdryEdge = false;
					for(unsigned int t7=0; t7<currentEdges.size(); t7++){
						if (m_mesh->isMarked(currentEdges[t7], m_mark_edges_on_curve)){
							hasBdryEdge = true;
							to_cell_dim_sp[i][j] = 1;
							to_cell_id_sp[i][j] = currentEdges[t7].id();
							vector<gmds::Node> currentNodes =  currentEdges[t7].get<Node>();
							math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());

							if(from_ray.SecondMetIntersect2D(oppSeg, intersectionPnt, intersectionParam, temp_epsilon)){
								double tempAngleVal = 1.0 - fabs((prevDirCross.second).dot(bdry_edge_normals[currentEdges[t7].id()]));
								if(tempAngleVal<gmds::math::Constants::PIDIV4){
									detectedBdry = true;
									distances[contSource][totalNumberOfSlots] = tempAngleVal;
									tempBdry.push_back(intersectionPnt);
									break;
								}
							}
						}
					}
					if(!hasBdryEdge){
						// VERTEX-> angle with vert normal
						vector<gmds::Node> currentNodes = source.get<gmds::Node>();
						for(unsigned int t7=0; t7<currentNodes.size(); t7++){
							if(m_mesh->isMarked(currentNodes[t7], m_mark_nodes_on_point) ||
							m_mesh->isMarked(currentNodes[t7], m_mark_nodes_on_curve)){
								to_cell_dim_sp[i][j] = 0;
								to_cell_id_sp[i][j] = currentNodes[t7].id();
 								vector<gmds::Edge> currentEdges = currentNodes[t7].get<gmds::Edge>();
								for(unsigned int t9=0; t9<currentEdges.size(); t9++){
									if(m_mesh->isMarked(currentEdges[t9], m_mark_edges_on_curve)){
										vector<gmds::Node> edge_nodes = currentEdges[t9].get<gmds::Node>();
										math::Segment oppSeg(edge_nodes[0].getPoint(), edge_nodes[1].getPoint());
										if(from_ray.SecondMetIntersect2D(oppSeg, intersectionPnt, intersectionParam, temp_epsilon)){
											double tempAngleVal = 1.0 - fabs((prevDirCross.second).dot(bdry_edge_normals[currentEdges[t9].id()]));
											if(tempAngleVal<gmds::math::Constants::PIDIV4){
												vector<gmds::Face> adj_faces_temp = currentEdges[t9].get<gmds::Face>();
												line_triangles[i][j].push_back(adj_faces_temp[0].id());
												detectedBdry = true;
												distances[contSource][totalNumberOfSlots] = tempAngleVal;
												//pointPaths[contSource][totalNumberOfSlots].first.push_back(intersectionPnt);
												tempBdry.push_back(intersectionPnt);
												break;
											}
										}
									}
								}
							}
						}
					}
				}

				if(detectedBdry){
					if(withGlobalComments)
						cout<<"detectedBdry before graph; contSource "<<contSource<<endl;
					pointPaths[contSource][totalNumberOfSlots].first.insert(pointPaths[contSource][totalNumberOfSlots].first.end(),line_discretization[i][j].begin(),line_discretization[i][j].end());
					pointPaths[contSource][totalNumberOfSlots].first.insert(pointPaths[contSource][totalNumberOfSlots].first.end(),tempBdry.begin(),tempBdry.end());
					if(visualizeTraversedTriangles){
						std::string file_name_trav_tri = "trav_tri_BDRY_"+std::to_string(5*i+j)+".vtk";
						writeTestMeshTrianglesIds(line_triangles[i][j], file_name_trav_tri);
					}

					SurfaceSingularityLine* surf_line = m_graph.newSurfaceLine();
					int sepNumberTmp = m_graph.getNbLines();
					surf_line->setNumber(sepNumberTmp);

					SingularityPoint* from_sing_pnt = singularity_points[i];
					std::vector<SingularityPoint::Slot*> pi_slots = from_sing_pnt->getSlots();

					SingularityPoint::Slot* current_slot = pi_slots[j];
					surf_line->addSingularityPoint(from_sing_pnt);

					current_slot->line = surf_line;
					math::Vector3d firstDir = math::Vector3d(from_sing_pnt->getLocation(), line_discretization[i][j][1]);
					current_slot->line_direction = firstDir;
					current_slot->isLaunched = true;

					for(unsigned int k3=0; k3<line_triangles[i][j].size(); k3++) {
						surf_line->addTraversedFace(line_triangles[i][j][k3]);
					}

					vector<math::Point> Newline_discretization = pointPaths[contSource][totalNumberOfSlots].first;
					gmds::Mesh meshSing(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
					gmds::Node mySing0 = meshSing.newNode(Newline_discretization[0].X(), Newline_discretization[0].Y(), Newline_discretization[0].Z());
					for(unsigned int j5=1; j5<Newline_discretization.size(); j5++) {
						gmds::Node mySing = meshSing.newNode(Newline_discretization[j5].X(), Newline_discretization[j5].Y(), Newline_discretization[j5].Z());
						meshSing.newTriangle(mySing0, mySing, mySing);
						mySing0 = mySing;
					}
					gmds::IGMeshIOService ioServiceSing(&meshSing);
					gmds::VTKWriter vtkWriterSing(&ioServiceSing);
					vtkWriterSing.setCellOptions(gmds::N|gmds::F);
					vtkWriterSing.setDataOptions(gmds::N|gmds::F);
					std::string file_name = m_output_directory_name+"-222discretizationPointsFoundBdry_"+std::to_string(5*i+j)+".vtk";

					vtkWriterSing.write(file_name);
					SingularityPoint* geom_pnt;
					SingularityPoint::Slot* incoming_slot;
					math::Vector3d toDir = math::Vector3d(pointPaths[contSource][totalNumberOfSlots].first[pointPaths[contSource][totalNumberOfSlots].first.size()-2], pointPaths[contSource][totalNumberOfSlots].first.back());

					createGeometricSingularityPoint(pointPaths[contSource][totalNumberOfSlots].first.back(),      // the last point added
					                                toDir,      // the direction we come from
					                                to_cell_dim_sp[i][j], // the dim. of the cell
					                                to_cell_id_sp[i][j],  // the id of the cell
					                                geom_pnt,       // the created point
					                                incoming_slot); // and the slot

					current_slot->lineDeviation = streamlineDeviation[i][j];
					surf_line->addSingularityPoint(geom_pnt);

					current_slot->isFreeze = true;
					current_slot->isLaunched = true;
					if(pointPaths[contSource][totalNumberOfSlots].first.size()>0){
						for(unsigned int j7=0; j7<pointPaths[contSource][totalNumberOfSlots].first.size(); j7++) {
							surf_line->addDiscretizationPoint(pointPaths[contSource][totalNumberOfSlots].first[j7]);
						}
					}

					fixedVariables.push_back(make_pair(i*5+j, i*5+j));
					fixedVariablesOpt.resize(fixedVariables.size());
					addedAlready[i][j] = true;
					possibleTargetTriangles[i][j] = original_faces_number + 1;
				}
			}
		}
	}

		//Here we start the Dijkstra's approach;

	for (unsigned int i = 0; i < singPointNo; i++) {
		if(withGlobalComments)
			cout<<"i= "<<i<<endl;
		for (unsigned int j = 0; j < originalSlotsNoPerSing[i]; j++){
			if(withGlobalComments)
				cout<<"i= "<<i<<" , j= "<<j<<endl;

			if(possibleTargetTriangles[i][j]<original_faces_number+1){
				contSource = 5*i+j;
				gmds::TCellID source_id = targetFaceOpt[contSource];
				gmds::Face source = m_mesh->get<Face>(source_id);
				vector<unsigned int> tempTargetFaces;
				vector<gmds::math::Vector3d> tempTargetDir;
				gmds::math::Vector3d prevDir = to_dir[i][j];
				prevDir.normalize();
				vector<pair<double, unsigned int>> min_distance;
				vector<int> previous, previousBdry;
				int found;
				vector<vector<double>> face2FaceMatch;
				pair<gmds::math::Vector3d, gmds::math::Vector3d> prevDirCross;
				gmds::math::Vector3d prevCross;

				prevDirCross.first = prevDir;
				prevDirCross.second = prevDir;

				auto t0Opti = Clock::now();

				vector<TCellID> tempTargetFaceOpt(5);
				vector<gmds::math::Point> tempTargetPoints(5);
				for(int jk=0;jk<5;jk++){
					tempTargetFaceOpt[jk] = targetFaceOpt[5*i+jk];
					targetFaceOpt[5*i+jk] = original_faces_number + 1;

					tempTargetPoints[jk] = targetPoints[5*i+jk];
					targetPoints[5*i+jk] = zeroPoint;
				}
				lastVisTri = original_faces_number;
				if(withGlobalComments)
					cout<<"for source_id "<<source_id<<endl;

				getShortestPathBtwFacesOptimized(source_id,
					                              targetFaceOpt,
					                              targetPoints,
					                              min_distance,
					                              previous,
					                              found,
					                              face2FaceMatch,
					                              face2FaceDeviation,
					                              prevDirCross,
					                              alltargetDirCross,
					                              is_bdry_face,
					                              last_discr_point[i][j],
					                              endPoint,
					                              maxDist,
					                              to_cell_dim_sp[i][j],
					                              to_cell_id_sp[i][j],
					                              lastVisTri,
					                              contSource,
					                              totalNumberOfVariables,
					                              totalNumberOfSlots,
					                              faceNo2Cont,
					                              line_discretization,
					                              line_triangles,
					                              pointPaths,
					                              finalPaths,
					                              finalCenterPoints,
					                              traversedSPTris,
					                              distances,
					                              IllegalCross,
					                              isTraversedNodeCompVector,
					                              isTraversedFaceCompVector,
					                              isTraversedFaceCompVector_SegmentPathCode,
					                              isTraversedFaceCompVector_SegmentPathCont,
					                              contSourceToSingularity,
					                              singOrGeomFaces);
				if(withGlobalComments)
					cout<<"done getShortestPathBtwFacesOptimized between source "<<source_id<<" and all others "<<endl;

				for(int jk=0;jk<5;jk++){
					targetFaceOpt[5*i+jk] = tempTargetFaceOpt[jk];
					targetPoints[5*i+jk] = tempTargetPoints[jk];
				}

//				auto t1Opti = Clock::now();
//				cout << "total execution time Opimized "
//				<< std::chrono::duration_cast<std::chrono::milliseconds>(t1Opti - t0Opti).count()
//				<< " milliseconds" << std::endl;
			}
			if(possibleTargetTriangles[i][j]==original_faces_number+1){
				if(withGlobalComments)
					cout<<"possibleTargetTriangles["<<i<<"]["<<j<<"]==original_faces_number+1"<<endl;
				for(unsigned int l1=0;l1<fixedVariables.size(); l1++){
					if(withGlobalComments)
						cout<<"fixedVariables[l1].first "<<fixedVariables[l1].first<<" , "<<fixedVariables[l1].second<<endl;

					if(fixedVariables[l1].first==fixedVariables[l1].second){
						if(withGlobalComments)
							cout<<"fixed equal "<<fixedVariables[l1].first<<endl;
						if((i==(int)(fixedVariables[l1].first/maxNoSing))&&(j==fmod(fixedVariables[l1].first,maxNoSing))){
							fixedVariablesOpt[l1].first = fixedVariables[l1].first*totalNumberOfVariables + totalNumberOfSlots;
							fixedVariablesOpt[l1].second = fixedVariablesOpt[l1].first;
							if(withGlobalComments)
								cout<<"0 fixedVariables[l1].first "<<fixedVariables[l1].first<<" , "<<fixedVariables[l1].second<<" -> "<<fixedVariablesOpt[l1].first<<" "<<fixedVariablesOpt[l1].second<<endl;
						}
					}
					else{
						if((i==(int)(fixedVariables[l1].first/maxNoSing))&&(j==fmod(fixedVariables[l1].first,maxNoSing))){
							fixedVariablesOpt[l1].first = fixedVariables[l1].first*totalNumberOfVariables + fixedVariables[l1].second;
							fixedVariablesOpt[l1].second = fixedVariables[l1].second*totalNumberOfVariables + fixedVariables[l1].first;
							if(withGlobalComments)
								cout<<"1 fixedVariables[l1].first "<<fixedVariables[l1].first<<" , "<<fixedVariables[l1].second<<" -> "<<fixedVariablesOpt[l1].first<<" "<<fixedVariablesOpt[l1].second<<endl;
						}
						else{
							if((i==(int)(fixedVariables[l1].second/maxNoSing))&&(j==fmod(fixedVariables[l1].second,maxNoSing))){
								fixedVariablesOpt[l1].first = fixedVariables[l1].first*totalNumberOfVariables + fixedVariables[l1].second;
								fixedVariablesOpt[l1].second = fixedVariables[l1].second*totalNumberOfVariables + fixedVariables[l1].first;
								if(withGlobalComments)
									cout<<"2 fixedVariables[l1].first "<<fixedVariables[l1].first<<" , "<<fixedVariables[l1].second<<" -> "<<fixedVariablesOpt[l1].first<<" "<<fixedVariablesOpt[l1].second<<endl;
							}
						}
					}
				}
			}
		}
	}
	auto t3 = Clock::now();
	cout << "graph construction "
	<< std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t0).count()
	<< " milliseconds" << std::endl;
	cout << "travTri "<< timerTravTri<< " milliseconds" << std::endl;

	bool fixSPBdry = false;
	//==================================================================
	// fixSPBdry
	/* Description of this approach (fixSPBdry==true)
	* it seems that if two singulairties are very close to one another and their slots have similar directions
	* (situation especially found for meshes with many interior geometric slots), then the corresponding found
	* boundary paths could intersect (especially if the mesh is coarse) in an illegal manner - merge usually.
	* since in many cases we would like to have these paths in the final graph, we can try two approaches :
		* 1) remesh (less chances to be succesfull) or
		* 2) this approach - fixSPBdry
		* fixSPBdry - detect in the IllegalCross pairs of illegal paths that are both towards the boundary
		* for each such pair, choose the path with the highest deviation and recompute it
		* if even after the recomputation the path still intersects with its problematic pair,
	* recompute that one as well;
	*/
	//==================================================================
	if(fixSPBdry){
		if(withGlobalComments){
			cout<<"fixSPBdry ; ";
			cout<<"IllegalCross before ";
			cout<<"IllegalCross.size() "<<IllegalCross.size()<<endl;
			for(unsigned int t=0; t<IllegalCross.size(); t++){
				cout<<"IllegalCross btw ("<<(int)(IllegalCross[t].first/totalNumberOfVariables)<<" , "
				<<fmod(IllegalCross[t].first,totalNumberOfVariables)<<") and ("<<
				(int)(IllegalCross[t].second/totalNumberOfVariables)<<" , "<<fmod(IllegalCross[t].second,totalNumberOfVariables)<<
				") ; ["<<IllegalCross[t].first<<" , "<<IllegalCross[t].second<<"]"<<endl;
			}
		}

		set<std::pair<double,unsigned int> > conflictingSPBdry;
		vector<std::pair<unsigned int,unsigned int> > SPToFixAndToRecompute;
		vector<bool> fixedIsFixedVar;
		unsigned int tempCont = 0;

		while(tempCont<IllegalCross.size()){
			unsigned int problematicSP1_first = (int) (IllegalCross[tempCont].first/totalNumberOfVariables); //contSource
			unsigned int problematicSP1_second = fmod(IllegalCross[tempCont].first,totalNumberOfVariables); //faceNo2Cont[u_id]

			unsigned int problematicSP2_first = (int) (IllegalCross[tempCont].second/totalNumberOfVariables); //contSource
			unsigned int problematicSP2_second = fmod(IllegalCross[tempCont].second,totalNumberOfVariables); //faceNo2Cont[u_id]

			if((problematicSP1_second==totalNumberOfSlots)&&(problematicSP2_second==totalNumberOfSlots)){
				unsigned int spToFix;
				unsigned int spToRedo;

				if(distances[problematicSP1_first][totalNumberOfSlots]<distances[problematicSP2_first][totalNumberOfSlots]){
					spToFix = problematicSP1_first;
					spToRedo = problematicSP2_first;
				}
				else{
					spToFix = problematicSP2_first;
					spToRedo = problematicSP1_first;
				}
				if(addedAlready[(int) (spToFix/5)][fmod(spToFix,5)]){
					//this is part of fixedVariables
					fixedIsFixedVar.push_back(true);
				}
				else{
					if(addedAlready[(int) (spToRedo/5)][fmod(spToRedo,5)]){
						fixedIsFixedVar.push_back(true);
						unsigned int tempFix = spToFix;
						spToFix = spToRedo;
						spToRedo = tempFix;
					}
					else
						fixedIsFixedVar.push_back(false);
				}
				if(withGlobalComments){
					cout<<"spToFix "<<spToFix<<" , spToRedo "<<spToRedo<<endl;
					cout<<"(int) (spToFix/5) "<<(int) (spToFix/5)<<"; fmod(spToFix,5) "<<fmod(spToFix,5)<<endl;
				}
				SPToFixAndToRecompute.push_back(make_pair(spToFix, spToRedo));
			}
			tempCont++;
		}

		vector<bool> recomputeIllegalProblematicBdryPair(SPToFixAndToRecompute.size(), false);

		for(unsigned int i=0; i<SPToFixAndToRecompute.size(); i++){
			unsigned int tempCS = SPToFixAndToRecompute[i].second;
			for(unsigned int j=0; j<traversedSPTris[tempCS][totalNumberOfSlots].size(); j++){
				unsigned int tempFace = traversedSPTris[tempCS][totalNumberOfSlots][j];
				std::vector<unsigned int>::iterator itTemp;
				itTemp = std::find (isTraversedFaceCompVector[tempFace].first.begin(), isTraversedFaceCompVector[tempFace].first.end(), tempCS);
				if(itTemp!=isTraversedFaceCompVector[tempFace].first.end()){
					int indexTemp = std::distance(isTraversedFaceCompVector[tempFace].first.begin(), itTemp);
					isTraversedFaceCompVector[tempFace].first.erase(isTraversedFaceCompVector[tempFace].first.begin()+indexTemp);
					isTraversedFaceCompVector_SegmentPathCode[tempFace].first.erase(isTraversedFaceCompVector_SegmentPathCode[tempFace].first.begin()+indexTemp);
					isTraversedFaceCompVector_SegmentPathCont[tempFace].first.erase(isTraversedFaceCompVector_SegmentPathCont[tempFace].first.begin()+indexTemp);
				}

				itTemp = std::find (isTraversedFaceCompVector[tempFace].second.begin(), isTraversedFaceCompVector[tempFace].second.end(), tempCS);
				if(itTemp!=isTraversedFaceCompVector[tempFace].second.end()){
					int indexTemp = std::distance(isTraversedFaceCompVector[tempFace].second.begin(), itTemp);
					isTraversedFaceCompVector[tempFace].second.erase(isTraversedFaceCompVector[tempFace].second.begin()+indexTemp);
					isTraversedFaceCompVector_SegmentPathCode[tempFace].second.erase(isTraversedFaceCompVector_SegmentPathCode[tempFace].second.begin()+indexTemp);
					isTraversedFaceCompVector_SegmentPathCont[tempFace].second.erase(isTraversedFaceCompVector_SegmentPathCont[tempFace].second.begin()+indexTemp);
				}
			}
		}
		vector<bool> singularTriangles(original_faces_number, false);
		for(auto sing_it:m_singularities_3){
			singularTriangles[sing_it.id()] = true;
		}
		for(auto sing_it:m_singularities_5){
			singularTriangles[sing_it.id()] = true;
		}

		vector<bool> recomputedPath(totalNumberOfSlots,false);
		for(unsigned int k=0; k<SPToFixAndToRecompute.size(); k++){
			unsigned int tempCS = SPToFixAndToRecompute[k].second;
			if(withGlobalComments)
				cout<<"try to recompute "<<tempCS<<" and/or "<<SPToFixAndToRecompute[k].first<<endl;

			if(!recomputedPath[tempCS]){
				recomputedPath[tempCS] = true;
				unsigned int problematicSL1_A = (int) (tempCS/totalNumberOfVariables); //contSource
				unsigned int problematicSL1_B = fmod(tempCS,totalNumberOfVariables); //faceNo2Cont[u_id]

				unsigned int i = (int) (tempCS/5);
				unsigned int j = fmod(tempCS,5);

				vector<SingularityPoint::Slot*> pi_slots = singularity_points[i]->getSlots();
				SingularityPoint::Slot* current_slot = pi_slots[j];
				SingularityPoint       *to_sing_pnt  = 0;
				SingularityPoint::Slot *to_slot = 0;

				bool find_end_bdry = false;
				bool	end_on_free_slot = false;
				bool find_end = false;

				bool noBetterChoice = false;

				vector<math::Point> copy_line_discretization(line_discretization[i][j]);
				vector<TCellID>     copy_line_triangles(line_triangles[i][j]);
				int 				copy_to_cell_dim = to_cell_dim[i][j];
				gmds::TCellID		copy_to_cell_id = to_cell_id[i][j];
				double 			copy_streamlineDeviation = streamlineDeviation[i][j];
				math::Vector3d      copy_to_dir = to_dir[i][j];

				while((!find_end_bdry)&&(!end_on_free_slot)&&(!noBetterChoice)&&(!find_end)){
					if(withGlobalComments)
						cout<<"line_triangles[i][j].back() "<<line_triangles[i][j].back()<<endl;
					growLineRK4(current_slot->from_point,
					            current_slot,
					            to_sing_pnt,
					            to_slot,
					            copy_line_discretization.back(),
					            copy_to_dir,
					            copy_line_discretization,
					            copy_line_triangles,
					            copy_to_cell_dim,
					            copy_to_cell_id,
					            copy_streamlineDeviation,
					            searchStep,
					            find_end_bdry,
					            end_on_free_slot,
					            find_end);
					if(withGlobalComments){
						cout<<"find_end_bdry "<<find_end_bdry<<", end_on_free_slot "<<end_on_free_slot<<", find_end "<<find_end;
						cout<<" , noBetterChoice "<<noBetterChoice<<endl;
					}
					if(singularTriangles[line_triangles[i][j].back()])
						noBetterChoice = true;
				}

				if(find_end_bdry){
					line_discretization[i][j].clear();
					line_discretization[i][j].insert(line_discretization[i][j].end(), copy_line_discretization.begin(), copy_line_discretization.end());
					line_triangles[i][j].clear();
					line_triangles[i][j].insert(line_triangles[i][j].end(), copy_line_triangles.begin(), copy_line_triangles.end());
					to_cell_dim[i][j] = copy_to_cell_dim;
					to_cell_id[i][j] = copy_to_cell_id;
					streamlineDeviation[i][j] = copy_streamlineDeviation;
					to_dir[i][j] = copy_to_dir;
					if(withGlobalComments)
						cout<<"find_end_bdry"<<endl;
						/*
          	           	if(m_build_geometric_singularities){
                         //WARNING TODO check if we land on the boundary in the vicinity of a geometric point
							cout<<"m_build_geometric_singularities"<<endl;
							for(auto sing_it:m_singularities_3){
	  							if(sing_it.id()==line_triangles[i][j].back()){

									cout<<"sing_it.id()==line_triangles[i][j][k].back()= "<<sing_it.id()<<endl;
								}
							}
						}*/

					math::Vector3d firstDir = math::Vector3d(line_discretization[i][j][0], line_discretization[i][j][1]);
					to_dir[i][j] = firstDir;
					distances[tempCS][totalNumberOfSlots] = streamlineDeviation[i][j];
					finalPaths[tempCS][totalNumberOfSlots].clear();

					pointPaths[tempCS][totalNumberOfSlots].first.clear();
					pointPaths[tempCS][totalNumberOfSlots].first.insert(pointPaths[tempCS][totalNumberOfSlots].first.end(),line_discretization[i][j].begin(),line_discretization[i][j].end());

					finalPaths[tempCS][totalNumberOfSlots].clear();

					lastVisTri = original_faces_number;

					unsigned int firstPointsNo = pointPaths[tempCS][totalNumberOfSlots].first.size();
					if(withGlobalComments)
						cout<<"firstPointsNo "<<firstPointsNo<<", finalPaths[tempCS][totalNumberOfSlots].size() "<<finalPaths[tempCS][totalNumberOfSlots].size()<<endl;
					if(line_discretization[i][j].size()>=2){
						std::string file_name_bdry = "ShortestPathsBdryNew_"+std::to_string(tempCS)+".vtk";
						writeTestPoints(line_discretization[i][j], file_name_bdry);
					}
					unsigned int illegalProblematicBdryPair = SPToFixAndToRecompute[k].first;

					unsigned int contMatrix = tempCS*totalNumberOfVariables + totalNumberOfSlots;
					contSourceToSingularity[contMatrix].first = i;
					contSourceToSingularity[contMatrix].second = j;

					//remove all from IllegalCross
					unsigned int tempContMatrix = tempCS*totalNumberOfVariables + totalNumberOfSlots;
					unsigned int tempCont = 0;

					while(tempCont<IllegalCross.size()){
						if((IllegalCross[tempCont].first==tempContMatrix)||(IllegalCross[tempCont].second==tempContMatrix)){
							IllegalCross.erase(IllegalCross.begin()+tempCont);
							tempCont--;
						}
						tempCont++;
					}
					bool recomputeIt = false;
					bool reAddIt = false;
					if(withGlobalComments)
						cout<<"before endPoint["<<tempCS<<"] "<<endPoint[tempCS]<<endl;

					endPoint[tempCS] = line_discretization[i][j].back();
					if(withGlobalComments){
						cout<<"endPoint["<<tempCS<<"] "<<endPoint[tempCS]<<" ; tempCS "<<tempCS<<" , contMatrix "<<contMatrix<<endl;
					}
					finalPaths[tempCS][totalNumberOfSlots].clear();
					getIllegalCrossByTraversedTris(line_discretization,
					                               line_triangles,
					                               traversedSPTris[tempCS][totalNumberOfSlots],
					                               finalPaths,
					                               isTraversedNodeCompVector,
					                               isTraversedFaceCompVector,
					                               isTraversedFaceCompVector_SegmentPathCode,
					                               isTraversedFaceCompVector_SegmentPathCont,
					                               contMatrix,
					                               totalNumberOfVariables,
					                               tempCS,
					                               totalNumberOfSlots,
					                               IllegalCross,
					                               contSourceToSingularity,
					                               lastVisTri,
					                               endPoint,
					                               illegalProblematicBdryPair,
					                               recomputeIt,
					                               reAddIt,
					                               i,
					                               j);
					traversedSPTris[tempCS][totalNumberOfSlots].clear();
					traversedSPTris[tempCS][totalNumberOfSlots].insert(traversedSPTris[tempCS][totalNumberOfSlots].end(), line_triangles[i][j].begin(), line_triangles[i][j].end());

					recomputeIllegalProblematicBdryPair[k] = recomputeIt;
					if(fixedIsFixedVar[k])
						recomputeIllegalProblematicBdryPair[k] = false;

					to_cell_dim_sp[i][j] = to_cell_dim[i][j];
					to_cell_id_sp[i][j] = to_cell_id[i][j];
					if(visualizeTraversedTriangles){
						std::string file_name_trav_tri = "trav_tri_BDRYnew_"+std::to_string(tempCS)+".vtk";
						writeTestMeshTrianglesIds(line_triangles[i][j], file_name_trav_tri);
					}
				}
			}
		}
		for(unsigned int k=0; k<recomputeIllegalProblematicBdryPair.size(); k++){
			if(recomputeIllegalProblematicBdryPair[k]){
				unsigned int tempCS = SPToFixAndToRecompute[k].first;
				if(!recomputedPath[tempCS]){
					recomputedPath[tempCS] = true;
					for(unsigned int j=0; j<traversedSPTris[tempCS][totalNumberOfSlots].size(); j++){
						unsigned int tempFace = traversedSPTris[tempCS][totalNumberOfSlots][j];
						std::vector<unsigned int>::iterator itTemp;
						itTemp = std::find (isTraversedFaceCompVector[tempFace].first.begin(), isTraversedFaceCompVector[tempFace].first.end(), tempCS);
						if(itTemp!=isTraversedFaceCompVector[tempFace].first.end()){
							int indexTemp = std::distance(isTraversedFaceCompVector[tempFace].first.begin(), itTemp);
							isTraversedFaceCompVector[tempFace].first.erase(isTraversedFaceCompVector[tempFace].first.begin()+indexTemp);
							isTraversedFaceCompVector_SegmentPathCode[tempFace].first.erase(isTraversedFaceCompVector_SegmentPathCode[tempFace].first.begin()+indexTemp);
							isTraversedFaceCompVector_SegmentPathCont[tempFace].first.erase(isTraversedFaceCompVector_SegmentPathCont[tempFace].first.begin()+indexTemp);
						}
						itTemp = std::find (isTraversedFaceCompVector[tempFace].second.begin(), isTraversedFaceCompVector[tempFace].second.end(), tempCS);
						if(itTemp!=isTraversedFaceCompVector[tempFace].second.end()){
							int indexTemp = std::distance(isTraversedFaceCompVector[tempFace].second.begin(), itTemp);
							isTraversedFaceCompVector[tempFace].second.erase(isTraversedFaceCompVector[tempFace].second.begin()+indexTemp);
							isTraversedFaceCompVector_SegmentPathCode[tempFace].second.erase(isTraversedFaceCompVector_SegmentPathCode[tempFace].second.begin()+indexTemp);
							isTraversedFaceCompVector_SegmentPathCont[tempFace].second.erase(isTraversedFaceCompVector_SegmentPathCont[tempFace].second.begin()+indexTemp);
						}
					}

					unsigned int problematicSL1_A = (int) (tempCS/totalNumberOfVariables);
					unsigned int problematicSL1_B = fmod(tempCS,totalNumberOfVariables);

 					unsigned int i = (int) (tempCS/5);
					unsigned int j = fmod(tempCS,5);

					vector<SingularityPoint::Slot*> pi_slots = singularity_points[i]->getSlots();
					SingularityPoint::Slot*     current_slot = pi_slots[j];
					SingularityPoint           *to_sing_pnt  = 0;
					SingularityPoint::Slot          *to_slot = 0;

					bool    find_end_bdry = false;
					bool end_on_free_slot = false;
					bool         find_end = false;
					bool   noBetterChoice = false;

					vector<math::Point> copy_line_discretization(line_discretization[i][j]);
					vector<TCellID>     copy_line_triangles(line_triangles[i][j]);
					int                 copy_to_cell_dim = to_cell_dim[i][j];
					gmds::TCellID        copy_to_cell_id = to_cell_id[i][j];
					double      copy_streamlineDeviation = streamlineDeviation[i][j];
					math::Vector3d           copy_to_dir = to_dir[i][j];

					while((!find_end_bdry)&&(!end_on_free_slot)&&(!noBetterChoice)&&(!find_end)){
						growLineRK4(current_slot->from_point,
										current_slot,
										to_sing_pnt,
										to_slot,
										copy_line_discretization.back(),
										copy_to_dir,
										copy_line_discretization,
										copy_line_triangles,
										copy_to_cell_dim,
										copy_to_cell_id,
										copy_streamlineDeviation,
										searchStep,
										find_end_bdry,
										end_on_free_slot,
										find_end);
						if(withGlobalComments)
							cout<<"find_end_bdry "<<find_end_bdry<<", end_on_free_slot "<<end_on_free_slot<<endl;
					}

					if(find_end_bdry){
						line_discretization[i][j].clear();
						line_discretization[i][j].insert(line_discretization[i][j].end(), copy_line_discretization.begin(), copy_line_discretization.end());
						line_triangles[i][j].clear();
						line_triangles[i][j].insert(line_triangles[i][j].end(), copy_line_triangles.begin(), copy_line_triangles.end());
						to_cell_dim[i][j] = copy_to_cell_dim;
						to_cell_id[i][j] = copy_to_cell_id;
						streamlineDeviation[i][j] = copy_streamlineDeviation;
						to_dir[i][j] = copy_to_dir;
							/*
						if(m_build_geometric_singularities){
                                   	cout<<"m_build_geometric_singularities"<<endl;
								for(auto sing_it:m_singularities_3){
	  								if(sing_it.id()==line_triangles[i][j].back()){

										cout<<"sing_it.id()==line_triangles[i][j][k].back()= "<<sing_it.id()<<endl;
									}
								}
							}*/

						math::Vector3d firstDir = math::Vector3d(line_discretization[i][j][0], line_discretization[i][j][1]);
						to_dir[i][j] = firstDir;
						distances[tempCS][totalNumberOfSlots] = streamlineDeviation[i][j];
						finalPaths[tempCS][totalNumberOfSlots].clear();

						pointPaths[tempCS][totalNumberOfSlots].first.clear();
						pointPaths[tempCS][totalNumberOfSlots].first.insert(pointPaths[tempCS][totalNumberOfSlots].first.end(),line_discretization[i][j].begin(),line_discretization[i][j].end());

						finalPaths[tempCS][totalNumberOfSlots].clear();

						lastVisTri = original_faces_number;

						unsigned int firstPointsNo = pointPaths[tempCS][totalNumberOfSlots].first.size();

						if(line_discretization[i][j].size()>=2){
							std::string file_name_bdry = "ShortestPathsBdryNew_"+std::to_string(tempCS)+".vtk";
							writeTestPoints(line_discretization[i][j], file_name_bdry);
						}
						unsigned int illegalProblematicBdryPair = SPToFixAndToRecompute[k].first;

						unsigned int contMatrix = tempCS*totalNumberOfVariables + totalNumberOfSlots;
						contSourceToSingularity[contMatrix].first = i;
						contSourceToSingularity[contMatrix].second = j;

						unsigned int tempContMatrix = tempCS*totalNumberOfVariables + totalNumberOfSlots;
						unsigned int tempCont = 0;

						while(tempCont<IllegalCross.size()){
							if((IllegalCross[tempCont].first==tempContMatrix)||(IllegalCross[tempCont].second==tempContMatrix)){
								IllegalCross.erase(IllegalCross.begin()+tempCont);
								tempCont--;
							}
							tempCont++;
						}
						bool recomputeIt = false;
						bool reAddIt = true;
						endPoint[tempCS] = line_discretization[i][j].back();
						if(withGlobalComments)
							cout<<"tempCS "<<tempCS<<" , contMatrix "<<contMatrix<<endl;
						finalPaths[tempCS][totalNumberOfSlots].clear();
  						getIllegalCrossByTraversedTris(line_discretization,
						                               line_triangles,
						                               traversedSPTris[tempCS][totalNumberOfSlots],
						                               finalPaths,
						                               isTraversedNodeCompVector,
						                               isTraversedFaceCompVector,
						                               isTraversedFaceCompVector_SegmentPathCode,
						                               isTraversedFaceCompVector_SegmentPathCont,
						                               contMatrix,
						                               totalNumberOfVariables,
						                               tempCS,
						                               totalNumberOfSlots,
						                               IllegalCross,
						                               contSourceToSingularity,
						                               lastVisTri,
						                               endPoint,
						                               illegalProblematicBdryPair,
						                               recomputeIt,
						                               reAddIt,
						                               i,
						                               j);

						traversedSPTris[tempCS][totalNumberOfSlots].clear();
						traversedSPTris[tempCS][totalNumberOfSlots].insert(traversedSPTris[tempCS][totalNumberOfSlots].end(), line_triangles[i][j].begin(), line_triangles[i][j].end());
						recomputeIllegalProblematicBdryPair[k] = recomputeIt;
						to_cell_dim_sp[i][j] = to_cell_dim[i][j];
						to_cell_id_sp[i][j] = to_cell_id[i][j];
						if(visualizeTraversedTriangles){
							std::string file_name_trav_tri = "trav_tri_BDRYnew_"+std::to_string(tempCS)+".vtk";
							writeTestMeshTrianglesIds(line_triangles[i][j], file_name_trav_tri);
						}
					}
				}
			}
		}
	}

		//////////////////////////////////////////////////////////////////////////////////////////////
          //REMESH AND REDO
          //////////////////////////////////////////////////////////////////////////////////////////////

         /* if(remeshConflictingTriangles){
               //REMESH
               vector<gmds::TCellID> trianglesToRemesh;
               vector<bool> trianglesToRemeshBool(original_faces_number, false);
          	cout<<"IllegalCross "<<endl;
			cout<<"IllegalCross.size() "<<IllegalCross.size()<<endl;
			for(unsigned int t=0; t<IllegalCross.size(); t++){
	  			cout<<"IllegalCross btw ("<<(int)(IllegalCross[t].first/totalNumberOfVariables)<<" , "
	  				<<fmod(IllegalCross[t].first,totalNumberOfVariables)<<") and ("<<
	  				(int)(IllegalCross[t].second/totalNumberOfVariables)<<" , "<<fmod(IllegalCross[t].second,totalNumberOfVariables)<<
	  				") ; ["<<IllegalCross[t].first<<" , "<<IllegalCross[t].second<<"]"<<endl;
        		}

        		gmds::Mesh newLocalMesh(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
			gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

               //WARNING TODO bdry new edges and nodes (and triangles)? new bdry edge normals...
			//WARNING TODO also new triangle_centers and triangle_centers_cross
		  	vector<TCellID> local_index_to_m_id;
		  	gmds::Variable<gmds::math::Cross2D>* local_cross_field_2D = newLocalMesh.newVariable<math::Cross2D,GMDS_NODE>("cross_X");
        		vector<bool> remeshedSL(totalNumberOfSlots*(totalNumberOfSlots+1),false);

        		while(!IllegalCross.empty()){
               	unsigned int problematicSL1 = IllegalCross.begin()->first;
                    unsigned int problematicSL2 = IllegalCross.begin()->second;
                    IllegalCross.erase(IllegalCross.begin());

                    unsigned int problematicSL1_A = (int) (problematicSL1/totalNumberOfVariables); //contSource
                    unsigned int problematicSL1_B = fmod(problematicSL1,totalNumberOfVariables); //faceNo2Cont[u_id]

                    unsigned int problematicSL1_A_i = (int) (problematicSL1_A/5);
		       	unsigned int problematicSL1_A_j = fmod(problematicSL1_A,5);
                    unsigned int problematicSL1_B_i = (int) (problematicSL1_B/5);
		       	unsigned int problematicSL1_B_j = fmod(problematicSL1_B,5);
                    //vector<gmds::TCellID> trianglesToRemesh1;
                    //remesh it
				 //if(problematicSL1_B==totalNumberOfSlots)   //bdry

                    for(unsigned int t=0; t<traversedSPTris[problematicSL1_A][problematicSL1_B].size(); t++){
					if(!trianglesToRemeshBool[traversedSPTris[problematicSL1_A][problematicSL1_B][t]]){
				   		trianglesToRemesh.push_back(traversedSPTris[problematicSL1_A][problematicSL1_B][t]);
						trianglesToRemeshBool[traversedSPTris[problematicSL1_A][problematicSL1_B][t]] = true;
					}
				}

                    remeshedSL[problematicSL1] = true;

                    unsigned int problematicSL2_A = (int) (problematicSL1/totalNumberOfVariables);
                    unsigned int problematicSL2_B = fmod(problematicSL1,totalNumberOfVariables);
                    if(problematicSL2_B==totalNumberOfSlots){
                           //bdry
                    }
                    unsigned int problematicSL2_A_i = (int) (problematicSL1_A/5);
		       	unsigned int problematicSL2_A_j = fmod(problematicSL1_A,5);
                    unsigned int problematicSL2_B_i = (int) (problematicSL1_B/5);
		       	unsigned int problematicSL2_B_j = fmod(problematicSL1_B,5);
                    //vector<gmds::TCellID> trianglesToRemesh2;
                     for(unsigned int t=0; t<traversedSPTris[problematicSL2_A][problematicSL2_B].size(); t++){
					if(!trianglesToRemeshBool[traversedSPTris[problematicSL2_A][problematicSL2_B][t]]){
				   		trianglesToRemesh.push_back(traversedSPTris[problematicSL2_A][problematicSL2_B][t]);
						trianglesToRemeshBool[traversedSPTris[problematicSL2_A][problematicSL2_B][t]] = true;
					}
				}

                    remeshedSL[problematicSL2] = true;
		    	}
		    	//remesh it
		    	vector<gmds::math::Point> newTriangleCenters;
			vector<gmds::math::Cross2D> newTriangleCenterCrosses;
               vector<gmds::math::Vector3d> newBdryEdgeNormals;
               vector<gmds::math::Vector3d> newBdryNodeNormals;
               vector<gmds::TCellID> newLocalMesh_id_to_mesh_id_node;
               vector<gmds::TCellID> mesh_id_to_newLocalMesh_id_node;
               vector<gmds::TCellID> newLocalMesh_id_to_mesh_id_face;
               vector<gmds::TCellID> mesh_id_to_newLocalMesh_id_face;
			vector<bool> isCurveEdge;
     		vector<bool> isCurveNode;
               remeshTrianglesNewMesh(&newLocalMesh, local_cross_field_2D, trianglesToRemesh,
							   trianglesToRemeshBool, newTriangleCenters, newTriangleCenterCrosses,
							   newBdryEdgeNormals, newBdryNodeNormals, isCurveEdge, isCurveNode,
							   newLocalMesh_id_to_mesh_id_node, mesh_id_to_newLocalMesh_id_node,
							   newLocalMesh_id_to_mesh_id_face, mesh_id_to_newLocalMesh_id_face);

//unsigned int totalNumberOfVariables = totalNumberOfSlots + 1;
//distances[contSource][totalNumberOfSlots]   |   distances[contSource][faceNo2Cont[u_id]]
//pointPaths[contSource][totalNumberOfSlots].first     |     pointPaths[contSource][faceNo2Cont[u_id]].first & Paths[contSource][faceNo2Cont[u_id]].second
//finalPaths[contSource]         |     finalPaths[contSource][faceNo2Cont[u_id]]
//WARNING TODO  finalCenterPoints[contSource][5*t1_orig+t2_orig]  ... for bdry these are determined afterwards...
//unsigned int t1_orig = (int) (faceNo2Cont[u_id]/5);
//unsigned int t2_orig = fmod(faceNo2Cont[u_id],5);
//traversedSPTris[contSource][totalNumberOfSlots]    |     traversedSPTris[contSource][faceNo2Cont[u_id]]
//IllegalCross.push_back(make_pair(tempCS, contMatrix));    tempCS = isTraversedFaceCompVector[currentFace_id].first[l];
//contMatrix = contSource*totalNumberOfVariables + totalNumberOfSlots;     contSourceToSingularity[contMatrix].first = i_orig;        contSourceToSingularity[contMatrix].second = i_orig;
//contMatrix = contSource*totalNumberOfVariables + faceNo2Cont[u_id];	contSourceToSingularity[contMatrix].first = i_orig;  	contSourceToSingularity[contMatrix].second = t1_orig;
//m_mesh

               //REDO
               vector<pair<unsigned int, unsigned int>> IllegalCrossFinal;
			for (unsigned int i = 0; i < singPointNo; i++) {
		  		cout<<"i= "<<i<<endl;
				for (unsigned int j = 0; j < originalSlotsNoPerSing[i]; j++){	 //bef slotsNoPerSing
			  		cout<<"i= "<<i<<" , j= "<<j<<endl;
					cout<<"originalSlotsNoPerSing[i] "<<originalSlotsNoPerSing[i]<<endl;

					if(possibleTargetTriangles[i][j]<original_faces_number+1){

    						gmds::Face source = newLocalMesh.get<Face>(possibleTargetTriangles[i][j]);

						gmds::TCellID source_id = possibleTargetTriangles[i][j];
    						vector<unsigned int> tempTargetFaces;
						vector<gmds::math::Vector3d> tempTargetDir;
						gmds::math::Vector3d prevDir = to_dir[i][j];
						prevDir.normalize();
						vector<pair<double, unsigned int>> min_distance;
                         	vector<int> previous, previousBdry;
						int found;
						vector<vector<double>> face2FaceMatch;
						pair<gmds::math::Vector3d, gmds::math::Vector3d> prevDirCross;
						gmds::math::Vector3d prevCross;
						// prevCross = triangle_centers_cross[source_id].closestComponentVector(prevDir);

						prevDirCross.first = prevDir;
						prevDirCross.second = prevDir;//prevCross;
						//vector<pair<gmds::math::Vector3d, gmds::math::Vector3d>> targetDirCross(1);

                         	contSource = 5*i+j;
						cout<<"bdry"<<endl;
						bool detectedBdry = false;
						vector<math::Point> tempBdry;
						if(is_bdry_face[source_id]){

							math::Point intersectionPnt;
							double intersectionParam;
							math::Segment from_ray(triangle_centers[source_id], triangle_centers[source_id] + prevDirCross.second);
							//cout<<"from_ray "<<triangle_centers[source_id]<<" "<<triangle_centers[source_id]+prevDirCross.second<<endl;
			 				vector<gmds::Edge> currentEdges = source.get<gmds::Edge>();
			 				bool hasBdryEdge = false;
			 				for(unsigned int t7=0; t7<currentEdges.size(); t7++){
								if (m_mesh->isMarked(currentEdges[t7], m_mark_edges_on_curve)){
									hasBdryEdge = true;
									to_cell_dim_sp[i][j] = 1;
									to_cell_id_sp[i][j] = currentEdges[t7].id();
									vector<gmds::Node> currentNodes =  currentEdges[t7].get<Node>();
									math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
									if(from_ray.SecondMetIntersect2D(oppSeg, intersectionPnt, intersectionParam, temp_epsilon)){
										double tempAngleVal = 1.0 - fabs((prevDirCross.second).dot(bdry_edge_normals[currentEdges[t7].id()]));
										if(tempAngleVal<gmds::math::Constants::PIDIV4){
								  			detectedBdry = true;
					    						distances[contSource][totalNumberOfSlots] = tempAngleVal;
											tempBdry.push_back(intersectionPnt);
											break;
										}
									}
								}
			  				}
			  				if(!hasBdryEdge){
							// VERTEX-> angle with vert normal
					 		// cout<<"!hasBdryEdge"<<endl;
			  					vector<gmds::Node> currentNodes = source.get<gmds::Node>();
								for(unsigned int t7=0; t7<currentNodes.size(); t7++){
									if(m_mesh->isMarked(currentNodes[t7], m_mark_nodes_on_point) ||
           						 		m_mesh->isMarked(currentNodes[t7], m_mark_nodes_on_curve)){
							  			to_cell_dim_sp[i][j] = 0;
										to_cell_id_sp[i][j] = currentNodes[t7].id();
 										vector<gmds::Edge> currentEdges = currentNodes[t7].get<gmds::Edge>();
										for(unsigned int t9=0; t9<currentEdges.size(); t9++){
					  						if(m_mesh->isMarked(currentEdges[t9], m_mark_edges_on_curve)){
						  						vector<gmds::Node> edge_nodes = currentEdges[t9].get<gmds::Node>();
												//cout<<"edge_nodes "<<edge_nodes[0]<<" "<<edge_nodes[1]<<endl;
												math::Segment oppSeg(edge_nodes[0].getPoint(), edge_nodes[1].getPoint());
												if(from_ray.SecondMetIntersect2D(oppSeg, intersectionPnt, intersectionParam, temp_epsilon)){
										  			double tempAngleVal = 1.0 - fabs((prevDirCross.second).dot(bdry_edge_normals[currentEdges[t9].id()]));
													if(tempAngleVal<gmds::math::Constants::PIDIV4){
											  			vector<gmds::Face> adj_faces_temp = currentEdges[t9].get<gmds::Face>();
					    									line_triangles[i][j].push_back(adj_faces_temp[0].id());
											  			detectedBdry = true;
								  						distances[contSource][totalNumberOfSlots] = tempAngleVal;
														tempBdry.push_back(intersectionPnt);
														break;
							    						}
												}
											}
										}
									}
			  					}
			  				}
						}
						cout<<"detectedBdry "<<detectedBdry<<" contSource "<<contSource<<endl;
						if(detectedBdry){
							cout<<"!!!detectedBdry"<<endl;
							pointPaths[contSource][totalNumberOfSlots].first.insert(pointPaths[contSource][totalNumberOfSlots].first.end(),line_discretization[i][j].begin(),line_discretization[i][j].end());
							pointPaths[contSource][totalNumberOfSlots].first.insert(pointPaths[contSource][totalNumberOfSlots].first.end(),tempBdry.begin(),tempBdry.end());
							if(visualizeTraversedTriangles){
								std::string file_name_trav_tri = "trav_tri_BDRY_"+std::to_string(5*i+j)+".vtk";

								writeTestMeshTrianglesIds(line_triangles[i][j], file_name_trav_tri);
							}
							//add sing line to graph
							SurfaceSingularityLine* surf_line = m_graph.newSurfaceLine();
    							int sepNumberTmp = m_graph.getNbLines();
    							surf_line->setNumber(sepNumberTmp);
							//connect line to initial singularity point
							SingularityPoint* from_sing_pnt = singularity_points[i];
							std::vector<SingularityPoint::Slot*> pi_slots = from_sing_pnt->getSlots();
    							//SingularityPoint* from_sing_pnt = pi_slots[j]->from_point;
							SingularityPoint::Slot* current_slot = pi_slots[j];
    							surf_line->addSingularityPoint(from_sing_pnt);
    							//surf_line->addDiscretizationPoint(from_sing_pnt->getLocation());
							current_slot->line = surf_line;
                         		math::Vector3d firstDir = math::Vector3d(from_sing_pnt->getLocation(), line_discretization[i][j][1]);
							current_slot->line_direction = firstDir;
							current_slot->isLaunched = true;

							for(unsigned int k3=0; k3<line_triangles[i][j].size(); k3++) {
    								surf_line->addTraversedFace(line_triangles[i][j][k3]);
							}

   							SingularityPoint* geom_pnt;
        						SingularityPoint::Slot* incoming_slot;
							math::Vector3d toDir = math::Vector3d(pointPaths[contSource][totalNumberOfSlots].first[pointPaths[contSource][totalNumberOfSlots].first.size()-2], pointPaths[contSource][totalNumberOfSlots].first.back());

        						createGeometricSingularityPoint(pointPaths[contSource][totalNumberOfSlots].first.back(),      // the last point added
                                 							toDir,      // the direction we come from
                                  							to_cell_dim_sp[i][j], // the dim. of the cell
                                  							to_cell_id_sp[i][j],  // the id of the cell
                                  							geom_pnt,       // the created point
                                  							incoming_slot); // and the slot
							//WARNING didn't add the last line dev
                        			current_slot->lineDeviation = streamlineDeviation[i][j];
							//surf_line->addDiscretizationPoint(geom_pnt->getLocation());
							surf_line->addSingularityPoint(geom_pnt);

							current_slot->isFreeze = true;
							cout<<"i= "<<i<<" , j= "<<j<<" current_slot->isFreeze bdry"<<endl;
							current_slot->isLaunched = true;
			    				cout<<"!!!create point "<<pointPaths[contSource][totalNumberOfSlots].first.back()<<endl;
				  			if(pointPaths[contSource][totalNumberOfSlots].first.size()>0){
								for(unsigned int j7=0; j7<pointPaths[contSource][totalNumberOfSlots].first.size(); j7++) {
						  			surf_line->addDiscretizationPoint(pointPaths[contSource][totalNumberOfSlots].first[j7]);
								}
	   						}
	   						addedAlready[i][j] = true;
						}

					  	auto t0Opti = Clock::now();

 					  	if(!detectedBdry){
							vector<TCellID> tempTargetFaceOpt(5);
							vector<gmds::math::Point> tempTargetPoints(5);
							for(int jk=0;jk<5;jk++){
								tempTargetFaceOpt[jk] = targetFaceOpt[5*i+jk];
								targetFaceOpt[5*i+jk] = original_faces_number + 1;

								tempTargetPoints[jk] = targetPoints[5*i+jk];
								targetPoints[5*i+jk] = zeroPoint;
							}
							lastVisTri = original_faces_number;
							cout<<"for source_id "<<source_id<<endl;

							getShortestPathBtwFacesOptimized(source_id,
                                      					targetFaceOpt,
												targetPoints,
                                           				min_distance,
                                           				previous,
								   				found,
								 				face2FaceMatch,
								 				face2FaceDeviation,
								 				prevDirCross,
									     		alltargetDirCross,
												//forbiddenFacesJustSing,
											    	is_bdry_face,
												last_discr_point[i][j],
												endPoint,
												maxDist,
												to_cell_dim_sp[i][j],
									    			to_cell_id_sp[i][j],
												lastVisTri,
												contSource,
												totalNumberOfVariables,
											    	totalNumberOfSlots,
												faceNo2Cont,
												line_discretization,
											    	line_triangles,
											    	pointPaths,
												finalPaths,
												finalCenterPoints,
												traversedSPTris,
												distances,
												IllegalCross,
												isTraversedNodeCompVector,
											    	isTraversedFaceCompVector,
												isTraversedFaceCompVector_SegmentPathCode,
												isTraversedFaceCompVector_SegmentPathCont,
												contSourceToSingularity,
												singOrGeomFaces);
							cout<<"getShortestPathBtwFacesOptimized between source "<<endl;
							cout<<source_id<<" and all others "<<endl;

							for(int jk=0;jk<5;jk++){
								targetFaceOpt[5*i+jk] = tempTargetFaceOpt[jk];
								targetPoints[5*i+jk] = tempTargetPoints[jk];
							}

							auto t1Opti = Clock::now();
							cout << "total execution time Non Opimized "
         							<< std::chrono::duration_cast<std::chrono::milliseconds>(t1Opti - t0Opti).count()
          						<< " milliseconds" << std::endl;
				     	}
					}
					//else{ can't be else...could add some in if
						if(possibleTargetTriangles[i][j]==original_faces_number+1){
					  		cout<<"possibleTargetTriangles["<<i<<"]["<<j<<"]==original_faces_number+1"<<endl;
				    			for(unsigned int l1=0;l1<fixedVariables.size(); l1++){

						  		if(fixedVariables[l1].first==fixedVariables[l1].second){
							  		cout<<"fixed equal "<<fixedVariables[l1].first<<endl;
					 				if((i==(int)(fixedVariables[l1].first/maxNoSing))&&(j==fmod(fixedVariables[l1].first,maxNoSing))){

					   					fixedVariablesOpt[l1].first = fixedVariables[l1].first*totalNumberOfVariables + totalNumberOfSlots;// + totalNumberOfSlots;
										fixedVariablesOpt[l1].second = fixedVariablesOpt[l1].first;//contSource*totalNumberOfVariables + totalNumberOfSlots;
										cout<<"0 fixedVariables[l1].first "<<fixedVariables[l1].first<<" , "<<fixedVariables[l1].second<<" -> "<<fixedVariablesOpt[l1].first<<" "<<fixedVariablesOpt[l1].second<<endl;
									}
						  		}
						  		else{
						    			if((i==(int)(fixedVariables[l1].first/maxNoSing))&&(j==fmod(fixedVariables[l1].first,maxNoSing))){

					   					//fixedVariablesOpt[l1].first = contSource*totalNumberOfVariables + fixedVariables[l1].second;
									  	fixedVariablesOpt[l1].first = fixedVariables[l1].first*totalNumberOfVariables + fixedVariables[l1].second;
										fixedVariablesOpt[l1].second = fixedVariables[l1].second*totalNumberOfVariables + fixedVariables[l1].first;
										cout<<"1 fixedVariables[l1].first "<<fixedVariables[l1].first<<" , "<<fixedVariables[l1].second<<" -> "<<fixedVariablesOpt[l1].first<<" "<<fixedVariablesOpt[l1].second<<endl;
									}
									else{
							  			if((i==(int)(fixedVariables[l1].second/maxNoSing))&&(j==fmod(fixedVariables[l1].second,maxNoSing))){
					   						//fixedVariablesOpt[l1].first = contSource*totalNumberOfVariables + fixedVariables[l1].second;
										  	fixedVariablesOpt[l1].first = fixedVariables[l1].first*totalNumberOfVariables + fixedVariables[l1].second;
											fixedVariablesOpt[l1].second = fixedVariables[l1].second*totalNumberOfVariables + fixedVariables[l1].first;
											cout<<"2 fixedVariables[l1].first "<<fixedVariables[l1].first<<" , "<<fixedVariables[l1].second<<" -> "<<fixedVariablesOpt[l1].first<<" "<<fixedVariablesOpt[l1].second<<endl;
										}
						    			}
						  		}
				    			}
				    			//contSourceToSing.push_back(i);
							//contSource++;
				  		}
					//}
				}
			}
			auto t3 = Clock::now();
     		cout << "graph construction "
              		<< std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t0).count()
              		<< " milliseconds" << std::endl;
	 		cout << "travTri "<< timerTravTri<< " milliseconds" << std::endl;
		} */

	//////////////////////////////////////////////////////////////////////////////////////////////

	/*TODO sometimes (a,b) exists but not (b,a) - in this case add (b,a) as (a,b) in the
	* optimization problem; therefore also add them in the IllegalCross*/
	if(withGlobalComments){
		cout<<"IllegalCross after "<<endl;
		cout<<"IllegalCross.size() "<<IllegalCross.size()<<endl;
		for(unsigned int t=0; t<IllegalCross.size(); t++){
			cout<<"IllegalCross btw ("<<(int)(IllegalCross[t].first/totalNumberOfVariables)<<" , "
			<<fmod(IllegalCross[t].first,totalNumberOfVariables)<<") and ("<<
			(int)(IllegalCross[t].second/totalNumberOfVariables)<<" , "<<fmod(IllegalCross[t].second,totalNumberOfVariables)<<
			") ; ["<<IllegalCross[t].first<<" , "<<IllegalCross[t].second<<"]"<<endl;
		}
		cout<<"fixedVariables"<<endl;
	}
	for(unsigned int i=0; i<fixedVariables.size(); i++){
		if(withGlobalComments){
			cout<<fixedVariables[i].first<<" , "<<fixedVariables[i].second<<endl;
			cout<<"variable x "<<(int)fixedVariables[i].first/totalNumberOfVariables<<" , "<<fmod(fixedVariables[i].first,totalNumberOfVariables)<<endl;
		}
		if(fixedVariables[i].first==true){
			if(fixedVariables[i].first==fixedVariables[i].second)
				distances[(int)fixedVariables[i].first/totalNumberOfVariables][fmod(fixedVariables[i].first,totalNumberOfVariables)] = 0.001;
			else
				throw GMDSException("SingularityGraphBuilder2D::treat case where it already finds slot2slot before graph");
		}
	}

	vector<double> w(totalNumberOfSlots*totalNumberOfVariables);
	vector<bool> validSingAndSlotsLines(totalNumberOfSlots*totalNumberOfVariables, true);
	vector<bool> validSingAndSlotsCols(totalNumberOfSlots*totalNumberOfVariables, true);
	if(withGlobalComments)
		cout<<"w "<<endl;

	for(unsigned int i=0; i<distances.size(); i++){
		for(unsigned int j=0; j<distances[i].size(); j++){
			w[i*totalNumberOfVariables+j] = distances[i][j];
			if(withGlobalComments)
				cout<<"distances["<<i<<"]["<<j<<"] = "<<distances[i][j]<<" ";
		}
		if(withGlobalComments)
			cout<<endl;
	}

	if(withGlobalComments){
		cout<<"IllegalCross "<<endl;
		cout<<"!!!distances "<<endl;
		cout<<"   "<<"| ";
		for(unsigned int i=0; i<distances.size(); i++)
			cout<<i<<"  ";
		cout<<endl;
		for(unsigned int i=0; i<distances.size(); i++){
			cout<<i<<"   | ";
			for(unsigned int j=0; j<distances[i].size(); j++)
				cout<<distances[i][j]<<"  ";
			cout<<endl;
		}

		cout<<"distances as optimization "<<endl;
		unsigned int deleteCont = 0;
		for(unsigned int i=0; i<distances.size(); i++){
			for(unsigned int j=0; j<distances[i].size(); j++){
				cout<<distances[i][j]<<"z_"<<(i*totalNumberOfVariables+j+1)<<" ";
				deleteCont++;
				if(deleteCont==3){
					cout<<endl;
					deleteCont = 0;
				}
			}
		}
	}
	vector<pair<unsigned int, unsigned int>> shortest(distances.size());
	for(unsigned int i=0; i<distances.size(); i++){
		double minDist = 100000000.0;
		for(unsigned int j=0; j<distances[i].size()-1; j++){
			if(distances[i][j]<minDist){
				shortest[i].first = i;
				shortest[i].second = j;
				minDist = distances[i][j];
			}
		}
		if(withGlobalComments)
			cout<<shortest[i].first<<" "<<shortest[i].second<<endl;
	}
	if(withGlobalComments)
		cout<<"distances.size() "<<distances.size()<<endl;

	for(unsigned int i=0; i<distances.size(); i++){
		if(withGlobalComments)
			cout<<"i= "<<i<<endl;
		unsigned int firstPointsNo = pointPaths[shortest[i].first][shortest[i].second].first.size();
		unsigned int lastPointsNo = pointPaths[shortest[i].first][shortest[i].second].second.size();
		if(withGlobalComments)
			cout<<"firstPointsNo "<<firstPointsNo<<", lastPointsNo "<<lastPointsNo<<endl;
		vector<math::Point> centerPoints(firstPointsNo + lastPointsNo + finalPaths[shortest[i].first][shortest[i].second].size());
		if(finalPaths[shortest[i].first][shortest[i].second].size()>1){
			for(unsigned int j=0; j<firstPointsNo;j++)
				centerPoints[j] = pointPaths[shortest[i].first][shortest[i].second].first[j];
			for(unsigned int j=0; j<finalPaths[shortest[i].first][shortest[i].second].size(); j++){
				centerPoints[firstPointsNo+j] = triangle_centers[finalPaths[shortest[i].first][shortest[i].second][j]];
			}

			for(unsigned int j=0; j<lastPointsNo;j++)
				centerPoints[j+firstPointsNo+finalPaths[shortest[i].first][shortest[i].second].size()] = pointPaths[shortest[i].first][shortest[i].second].second[j];
		}
	}


	//////////////////////////////////////////////////////////////////////////////////////////////

	unsigned int varNumber = totalNumberOfSlots*totalNumberOfVariables;
	if(withGlobalComments){
		cout<<"!!!optimization step ; totalNumberOfSlots = "<<totalNumberOfSlots<<";"<<" , varNumber = "<<varNumber<<endl;
		cout<<"IllegalCross.size() = "<<IllegalCross.size()<<endl;
	}
	auto t1 = Clock::now();

	glp_prob *lp;
	int *ia, *ja;
	double *ar;

	lp = glp_create_prob();

	glp_set_prob_name(lp, "SingGraphBuildOpt");
	glp_set_obj_dir(lp, GLP_MIN);
	vector<bool> linesNotSum2(totalNumberOfSlots, false);
	vector<int> testMat(totalNumberOfSlots*totalNumberOfVariables, 33);
	//objective function //WARNING TODO if non launched slots -> their sum should be 0, not 2; or dont put condition
	glp_add_cols(lp, totalNumberOfSlots*totalNumberOfVariables);

	vector<bool> fixedZeroVars(totalNumberOfSlots*totalNumberOfVariables,false);

	for(unsigned int i=0; i< singPointNo; i++){
		for(unsigned int j=0; j< 5; j++){
			if(!possibleVariablesSlotsLaunched[5*i+j]){
				//entire lines fixed to 0
				for(unsigned int k=0; k<totalNumberOfVariables; k++){
					unsigned int tempVar = (5*i+j)*totalNumberOfVariables + k;
					if(!fixedZeroVars[tempVar]){
						linesNotSum2[5*i+j] = true;

						fixedZeroVars[tempVar] = true;
						glp_set_col_bnds(lp, tempVar+1, GLP_FX, 0.0, 0.0);
						glp_set_obj_coef(lp, tempVar+1, 0);
						glp_set_obj_coef(lp, tempVar+1, 0);
						testMat[tempVar] = 0;
					}
				}
				//columns fixed to 0
				for(unsigned int k=0; k<totalNumberOfSlots; k++){
					unsigned int tempVar = totalNumberOfVariables * k + (5*i+j);
					if(!fixedZeroVars[tempVar]){
						fixedZeroVars[tempVar] = true;
						glp_set_col_bnds(lp, tempVar+1, GLP_FX, 0.0, 0.0);
						glp_set_obj_coef(lp, tempVar+1, 0);
						glp_set_obj_coef(lp, tempVar+1, 0);
						testMat[tempVar] = 0;
					}
				}
			}
		}
	}

 	// CONDITION: slots corresponding to same singularity are not allowed to be connected in between themselves    ; however this condition should be removed

	for(unsigned int i=0; i< totalNumberOfSlots; i++){
		for(unsigned int j=0; j<5; j++){
			unsigned int tempVar = totalNumberOfVariables * i + (5* ((int)i/5)) + j;
			if(!fixedZeroVars[tempVar]){
				fixedZeroVars[tempVar] = true;
				glp_set_col_bnds(lp, tempVar+1, GLP_FX, 0.0, 0.0);
				glp_set_obj_coef(lp, tempVar+1, 0);
				glp_set_obj_coef(lp, tempVar+1, 0);
				testMat[tempVar] = 0;
			 }
		}
	}
	if(withGlobalComments){
		cout<<"before fixedVariables "<<endl;
		for(unsigned int i=0; i<fixedVariables.size(); i++){
			cout<<"fixedVariables["<<i<<"] "<<fixedVariables[i].first<<" (sing "<<(int)(fixedVariables[i].first/5)<<") slot (, "<<fmod(fixedVariables[i].first,5)<<")"<<endl;
			cout<<"towards "<<endl;
			cout<<"fixedVariables["<<i<<"] "<<fixedVariables[i].second<<" (sing "<<(int)(fixedVariables[i].second/5)<<") slot (, "<<fmod(fixedVariables[i].second,5)<<")"<<endl;
		}
	}
	vector<bool> setFixedVars(fixedVariablesOpt.size(),false);

	for(unsigned int i=0; i<totalNumberOfSlots*totalNumberOfVariables; i++){
		bool setVal = false;

		for(unsigned int j=0; j<fixedVariablesOpt.size(); j++){
			if((fixedVariablesOpt[j].first==i)||(fixedVariablesOpt[j].second==i)){
				setVal = true;
				if(!setFixedVars[j]){
					setFixedVars[j] = true;
					if(withGlobalComments){
						cout<<"sets fixedVariables true "<<fixedVariablesOpt[j].first<<" and "<<fixedVariablesOpt[j].second<<endl;
						cout<<"fixedVar.f slot a = (int) fixedVariablesOpt[j].first/(totalNumberOfSlots+1) "<<(int) fixedVariablesOpt[j].first/(totalNumberOfSlots+1)<<endl;
						cout<<"fixedVar.f slot b = std::fmod(fixedVariablesOpt[j].first,(totalNumberOfSlots+1)) "<<std::fmod(fixedVariablesOpt[j].first,(totalNumberOfSlots+1))<<endl;
						cout<<"fixedVar.s slot a = (int) fixedVariablesOpt[j].second/(totalNumberOfSlots+1) "<<(int) fixedVariablesOpt[j].second/(totalNumberOfSlots+1)<<endl;
						cout<<"fixedVar.s slot b = std::fmod(fixedVariablesOpt[j].second,(totalNumberOfSlots+1)) "<<std::fmod(fixedVariablesOpt[j].second,(totalNumberOfSlots+1))<<endl;
					}
					glp_set_col_bnds(lp, fixedVariablesOpt[j].first+1, GLP_FX, 1.0, 1.0);
					glp_set_obj_coef(lp, fixedVariablesOpt[j].first+1, 0.0001);
					testMat[fixedVariablesOpt[j].first] = 1;
				}
			}
		}

		if(!setVal){
			glp_set_col_kind(lp, i+1 , GLP_BV); //GLP_IV - integer, but in our case we want binary
			glp_set_obj_coef(lp, i+1, w[i]);
		}
	}

	glp_add_rows(lp, totalNumberOfSlots+IllegalCross.size());//!!!!!!!!!!!

	//first equality constraints: all valid slots must have exactly 1 singularity line associated
	for(unsigned int i=0; i<totalNumberOfSlots; i++){
		if(possibleVariablesSlotsLaunched[i]){
			glp_set_row_bnds(lp, i+1, GLP_FX, 1, 2);
		}
		else
			glp_set_row_bnds(lp, i+1, GLP_FX, 0, 0);
	}

	//second equality constraints: all IllegalCross pairs must have a maximum sum of 1 (at most 1 of them should be present in the final solution)
	for(unsigned int i=totalNumberOfSlots; i<totalNumberOfSlots+IllegalCross.size(); i++){
		glp_set_row_bnds(lp, i+1, GLP_LO, 0.0, 0);
		glp_set_row_bnds(lp, i+1, GLP_UP, 0.0, 1);
	}

	// dynamic allocation of the matrix where ia represents the row index, ja the col index and ar the entry
	ia = (int *) calloc(2 + (2*varNumber + IllegalCross.size()*2), sizeof(int));
	ja = (int *) calloc(2 + (2*varNumber + IllegalCross.size()*2), sizeof(int));
	ar = (double *) calloc(2 + (2*varNumber + IllegalCross.size()*2), sizeof(double));

	unsigned int count = 1;

	for (unsigned int i = 0; i < totalNumberOfSlots; i++) {
		for (unsigned int j = 0; j <totalNumberOfVariables; j++) {
			ia[count] = i + 1;
			ja[count] = totalNumberOfVariables*i+j+1;
			ar[count] = 1;
			count++;
			if((j!=totalNumberOfSlots)&&(i!=j)){
				ia[count] = j+1;
				ja[count] = totalNumberOfVariables*i+j+1;
				ar[count] = 1;
				count++;
			}
		}
	}

	if(withGlobalComments)
		cout<<"count after Slots "<<count<<endl;

	for (unsigned int i = 0; i < IllegalCross.size(); i++) {
		ia[count] = totalNumberOfVariables+i;
		ja[count] = IllegalCross[i].first+1;
		ar[count] = 1;
		count++;

		ia[count] = totalNumberOfVariables+i;
		ja[count] = IllegalCross[i].second+1;
		ar[count] = 1;
		count++;
	}
	if(withGlobalComments){
		cout<<"allocated "<<(varNumber + IllegalCross.size()*2)<<endl;
	}
	glp_load_matrix(lp, count-1, ia, ja, ar);

	// the next line disables the default terminal report
	glp_term_out(GLP_ON);
	glp_iocp glpParams;
	glp_init_iocp(&glpParams);
	glpParams.presolve = GLP_ON;

	glp_write_lp(lp, NULL, "checkMe0.txt");

	int glpErr = 0;
	glpErr = glp_intopt(lp, &glpParams);
	switch (glpErr) {
		case 0:
			cout << "GLP OK" << std::endl;
			break;
		default:
			std::cout << "pb solving in GLP." << std::endl;
			throw GMDSException("SingularityGraphBuilder2D::SingGraphBuildOpt pb solving in GLPK.");
			break;
	}
	glp_print_sol(lp, "SingGraphBuildOpt.txt");
	cout<<"status of MIP sol "<<glp_mip_status(lp)<<endl;
	for (unsigned int i = 0; i < varNumber; i++) {
		if(glp_mip_col_val(lp, i + 1)!=0){
			cout<<(int)(i/totalNumberOfVariables) <<" "<<fmod(i,totalNumberOfVariables)<< endl;
		}
	}

	gmds::Mesh meshSingPaths(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
	vector<SurfaceSingularityLine* > surf_lines = m_graph.getSurfaceLines();
	for(int i=0;i<surf_lines.size();i++){
		vector<math::Point> surf_line_discretization = surf_lines[i]->getDiscretizationPoints();
		gmds::Node mySing0 = meshSingPaths.newNode(surf_line_discretization[0].X(), surf_line_discretization[0].Y(), surf_line_discretization[0].Z());

		for(unsigned int j=1; j<surf_line_discretization.size(); j++) {
			gmds::Node mySing = meshSingPaths.newNode(surf_line_discretization[j].X(), surf_line_discretization[j].Y(), surf_line_discretization[j].Z());
			meshSingPaths.newTriangle(mySing0, mySing, mySing);
			mySing0 = mySing;
		}
	}
	if(withGlobalComments)
		cout<<"totalNumberOfVariables-1 "<<totalNumberOfVariables-1<<", totalNumberOfSlots "<<totalNumberOfSlots<<endl;

	unsigned int temp_i, temp_j;

	//file_name = "bdyLine";
	// writeOutput(file_name);

	for (unsigned int i = 0; i < varNumber; i++) {
		if(glp_mip_col_val(lp, i + 1)!=0){
			unsigned int i_local = (int)(i/totalNumberOfVariables);
			temp_i = (int)(i_local/5);
			temp_j = fmod(i_local,5);
			if(fmod(i,totalNumberOfVariables)==totalNumberOfSlots){
				if(!addedAlready[temp_i][temp_j]){
					if(withGlobalComments)
						cout<<"if(!addedAlready["<<temp_i<<"]["<<temp_j<<"]){"<<endl;
					SurfaceSingularityLine* surf_line = m_graph.newSurfaceLine();
					int sepNumberTmp = m_graph.getNbLines();
					surf_line->setNumber(sepNumberTmp);
					//connect line to initial singularity point
					SingularityPoint* from_sing_pnt = singularity_points[temp_i];
					std::vector<SingularityPoint::Slot*> pi_slots = from_sing_pnt->getSlots();
					SingularityPoint::Slot* current_slot = pi_slots[temp_j];
					surf_line->addSingularityPoint(from_sing_pnt);
					current_slot->line = surf_line;
					math::Vector3d firstDir = math::Vector3d(from_sing_pnt->getLocation(), line_discretization[temp_i][temp_j][1]);

					current_slot->line_direction = firstDir;
					current_slot->isLaunched = true;
					//line_triangles[temp_i][temp_j].back()==traversedSPTris[i_local][totalNumberOfSlots][0]
					//line_triangles[temp_i][temp_j][0] -> sing tri
					for(unsigned int k=1; k<line_triangles[temp_i][temp_j].size(); k++) {
						surf_line->addTraversedFace(line_triangles[temp_i][temp_j][k]);
					}
					for(unsigned int k=1; k<traversedSPTris[i_local][totalNumberOfSlots].size(); k++) {
						surf_line->addTraversedFace(traversedSPTris[i_local][totalNumberOfSlots][k]);
					}
					SingularityPoint* geom_pnt;
					SingularityPoint::Slot* incoming_slot;
					if(withGlobalComments)
						cout<<"endPoint[i_local] "<<endPoint[i_local]<<endl;
					math::Vector3d toDir = math::Vector3d(triangle_centers[finalPaths[i_local][totalNumberOfSlots].back()], endPoint[i_local]);

					createGeometricSingularityPoint(endPoint[i_local],      // the last point added
					                                toDir,      // the direction we come from
					                                to_cell_dim_sp[temp_i][temp_j], // the dim. of the cell
					                                to_cell_id_sp[temp_i][temp_j],  // the id of the cell
					                                geom_pnt,       // the created point
					                                incoming_slot); // and the slot
					current_slot->lineDeviation = distances[i_local][totalNumberOfSlots];

					surf_line->addSingularityPoint(geom_pnt);

					current_slot->isFreeze = true;
					if(withGlobalComments)
						cout<<"temp_i= "<<temp_i<<" , temp_j= "<<temp_j<<" current_slot->isFreeze bdry"<<endl;
					current_slot->isLaunched = true;
					unsigned int firstPointsNo = pointPaths[i_local][totalNumberOfSlots].first.size();
					unsigned int finalPathsNo = finalPaths[i_local][totalNumberOfSlots].size();

					if(finalPathsNo>=1){
						vector<math::Point> centerPoints(firstPointsNo + finalPathsNo);
						for(unsigned int j=0; j<firstPointsNo;j++){
							centerPoints[j] = pointPaths[i_local][totalNumberOfSlots].first[j];
						}

	 					for(unsigned int j=0; j<finalPathsNo; j++){
							centerPoints[j+firstPointsNo] = triangle_centers[finalPaths[i_local][totalNumberOfSlots][j]];
						}

						centerPoints.push_back(endPoint[i_local]);

						if(centerPoints.size()>=2){
							std::string file_name = "ShortestPathsBdryEnd_"+std::to_string(i_local)+".vtk";
							writeTestPoints(centerPoints, file_name);
						}
						for(unsigned int j=0; j<centerPoints.size(); j++) {
							surf_line->addDiscretizationPoint(centerPoints[j]);
						}
						if(visualizeTraversedTriangles){
							vector<gmds::TCellID> travTri = surf_line->getTraversedFaces();
							std::string file_name_trav_tri = "trav_tri_BDRY_"+std::to_string(i_local)+".vtk";
							writeTestMeshTrianglesIds(travTri, file_name_trav_tri);
						}
					}
					else{
						if(firstPointsNo>0){
							for(unsigned int j=0; j<firstPointsNo; j++) {
							   surf_line->addDiscretizationPoint(pointPaths[i_local][totalNumberOfSlots].first[j]);
							}
							if(pointPaths[i_local][totalNumberOfSlots].first.size()>=2){
								std::string file_name = "ShortestPathsBdryEnd_"+std::to_string(i_local)+".vtk";
								writeTestPoints(pointPaths[i_local][totalNumberOfSlots].first, file_name);
							}
						}
						if(visualizeTraversedTriangles){
							vector<gmds::TCellID> travTri = surf_line->getTraversedFaces();
							std::string file_name_trav_tri = "trav_tri_BDRY_"+std::to_string(i_local)+".vtk";
							writeTestMeshTrianglesIds(travTri, file_name_trav_tri);
						}
					}
				    //std::string file_name = "bdyLine";
					//writeOutput(file_name);
			 	}

			 	/////////////////////////////////////////////////////////

				unsigned int firstPointsNo = pointPaths[i_local][totalNumberOfSlots].first.size();
				unsigned int finalPathsNo = finalPaths[i_local][totalNumberOfSlots].size();

				if(finalPathsNo>=1){
					vector<math::Point> centerPoints(firstPointsNo + finalPathsNo);
					for(unsigned int j=0; j<firstPointsNo;j++){
						centerPoints[j] = pointPaths[i_local][totalNumberOfSlots].first[j];
					}

					for(unsigned int j=0; j<finalPathsNo; j++){
						centerPoints[j+firstPointsNo] = triangle_centers[finalPaths[i_local][totalNumberOfSlots][j]];
					}

					centerPoints.push_back(endPoint[i_local]);

					gmds::Node mySing0 = meshSingPaths.newNode(centerPoints[0].X(), centerPoints[0].Y(), centerPoints[0].Z());

					for(unsigned int j=1; j<centerPoints.size(); j++) {
						gmds::Node mySing = meshSingPaths.newNode(centerPoints[j].X(), centerPoints[j].Y(), centerPoints[j].Z());
						meshSingPaths.newTriangle(mySing0, mySing, mySing);
						mySing0 = mySing;
					}
				}
				else{

					if(firstPointsNo>0){

						gmds::Node mySing0 = meshSingPaths.newNode(pointPaths[i_local][totalNumberOfSlots].first[0].X(), pointPaths[i_local][totalNumberOfSlots].first[0].Y(), pointPaths[i_local][totalNumberOfSlots].first[0].Z());

						for(unsigned int j=1; j<firstPointsNo; j++) {
							gmds::Node mySing = meshSingPaths.newNode(pointPaths[i_local][totalNumberOfSlots].first[j].X(), pointPaths[i_local][totalNumberOfSlots].first[j].Y(), pointPaths[i_local][totalNumberOfSlots].first[j].Z());
							meshSingPaths.newTriangle(mySing0, mySing, mySing);
							mySing0 = mySing;
						}
					}
				}
			}
			else{// not boundary sing line
				if(!addedAlready[temp_i][temp_j]){
					if(finalCenterPoints[i_local][fmod(i,totalNumberOfVariables)].size()>1){
						SurfaceSingularityLine* surf_line = m_graph.newSurfaceLine();
						int sepNumberTmp = m_graph.getNbLines();
						surf_line->setNumber(sepNumberTmp);
						SingularityPoint* from_sing_pnt = singularity_points[temp_i];
						std::vector<SingularityPoint::Slot*> pi_slots = from_sing_pnt->getSlots();
						SingularityPoint::Slot* current_slot = pi_slots[temp_j];

						surf_line->addSingularityPoint(from_sing_pnt);

						current_slot->line = surf_line;

						math::Vector3d firstDir = math::Vector3d(from_sing_pnt->getLocation(), line_discretization[temp_i][temp_j][1]);

						current_slot->line_direction = firstDir;
						current_slot->isLaunched = true;
						current_slot->lineDeviation = streamlineDeviation[temp_i][temp_j];

						//line_triangles[temp_i][temp_j].back()==traversedSPTris[i_local][fmod(i,totalNumberOfVariables)][0]

						for(unsigned int k=1; k<line_triangles[temp_i][temp_j].size();k++) {
							surf_line->addTraversedFace(line_triangles[temp_i][temp_j][k]);
						}
						for(unsigned int k=1; k<traversedSPTris[i_local][fmod(i,totalNumberOfVariables)].size(); k++) {
							surf_line->addTraversedFace(traversedSPTris[i_local][fmod(i,totalNumberOfVariables)][k]);
						}

						unsigned int temp_ii, temp_jj;
						temp_ii = (int)((fmod(i,totalNumberOfVariables))/5);
						temp_jj = fmod(fmod(i,totalNumberOfVariables),5);

						SingularityPoint* to_sing_pnt = singularity_points[temp_ii];
						std::vector<SingularityPoint::Slot*> pi_to_slots = to_sing_pnt->getSlots();

						SingularityPoint::Slot* to_slot = pi_to_slots[temp_jj];

						to_slot->isLaunched = true;
						current_slot->isFreeze = true;
						to_slot->isFreeze = true;

						to_slot->line = surf_line;
						to_slot->line_direction = to_dir[temp_i][temp_j];

						to_slot->lineDeviation = streamlineDeviation[temp_i][temp_j];

						surf_line->addSingularityPoint(to_sing_pnt);
						for(unsigned int j=0; j<finalCenterPoints[i_local][fmod(i,totalNumberOfVariables)].size(); j++) {
							surf_line->addDiscretizationPoint(finalCenterPoints[i_local][fmod(i,totalNumberOfVariables)][j]);
						}

						if(finalCenterPoints[i_local][fmod(i,totalNumberOfVariables)].size()>=2){
							std::string file_name = "ShortestPathsEnd_"+std::to_string(i_local)+"-"+std::to_string(fmod(i,totalNumberOfVariables))+".vtk";
							writeTestPoints(finalCenterPoints[i_local][fmod(i,totalNumberOfVariables)], file_name);
						}
						if(visualizeTraversedTriangles){
							vector<gmds::TCellID> travTri = surf_line->getTraversedFaces();
							std::string file_name_trav_tri = "trav_tri_"+std::to_string(i_local)+"-"+std::to_string(fmod(i,totalNumberOfVariables))+".vtk";
							writeTestMeshTrianglesIds(travTri, file_name_trav_tri);
						}

					}

				}
				if(withGlobalComments)
					cout<<"finalCenterPoints[i_local][fmod(i,totalNumberOfVariables)].size() "<<finalCenterPoints[i_local][fmod(i,totalNumberOfVariables)].size()<<endl;
				if(finalCenterPoints[i_local][fmod(i,totalNumberOfVariables)].size()>1){
					gmds::Node mySing0 = meshSingPaths.newNode(finalCenterPoints[i_local][fmod(i,totalNumberOfVariables)][0].X(), finalCenterPoints[i_local][fmod(i,totalNumberOfVariables)][0].Y(), finalCenterPoints[i_local][fmod(i,totalNumberOfVariables)][0].Z());

					for(unsigned int j=1; j<finalCenterPoints[i_local][fmod(i,totalNumberOfVariables)].size(); j++) {
						gmds::Node mySing = meshSingPaths.newNode(finalCenterPoints[i_local][fmod(i,totalNumberOfVariables)][j].X(), finalCenterPoints[i_local][fmod(i,totalNumberOfVariables)][j].Y(), finalCenterPoints[i_local][fmod(i,totalNumberOfVariables)][j].Z());
						meshSingPaths.newTriangle(mySing0, mySing, mySing);
						mySing0 = mySing;
					}
				}
				else
					cout<<"error: Graph -> finalCenterPoints["<<i_local <<"]["<< fmod(i,totalNumberOfVariables)<<"] doesn't exist"<<endl;

			}
		}
	}

	gmds::IGMeshIOService ioServiceSingPaths(&meshSingPaths);
	gmds::VTKWriter vtkWriterSingPaths(&ioServiceSingPaths);
	vtkWriterSingPaths.setCellOptions(gmds::N|gmds::F);
	vtkWriterSingPaths.setDataOptions(gmds::N|gmds::F);
	std::string file_name_paths = m_output_directory_name+"-ShortestPaths_result.vtk";
	vtkWriterSingPaths.write(file_name_paths);
	if(withGlobalComments)
		cout<<"has written ShortestPaths_result "<<endl;

	/*vector<SurfaceSingularityLine*> surfSingLines = m_graph.getSurfaceLines();

	for(unsigned int i=0; i<surfSingLines.size(); i++){
		cout<<"i "<<i<<" getNo "<<surfSingLines[i]->getNumber()<<endl;
		vector<gmds::math::Point> discrPts = surfSingLines[i]->getDiscretizationPoints();
		if(discrPts.size()>0){
		  	gmds::Mesh meshSingx(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
			gmds::Node firstPt = meshSingx.newNode(discrPts[0].X(), discrPts[0].Y(), discrPts[0].Z());
			for(unsigned int j=1; j<discrPts.size(); j++){
	 			 gmds::Node nextPt = meshSingx.newNode(discrPts[j].X(), discrPts[j].Y(), discrPts[j].Z());
	  			meshSingx.newTriangle(firstPt, nextPt, nextPt);
				firstPt = nextPt;
			}

			gmds::IGMeshIOService ioServiceSingx(&meshSingx);
			gmds::VTKWriter vtkWriterSingx(&ioServiceSingx);
			vtkWriterSingx.setCellOptions(gmds::N|gmds::F);
			vtkWriterSingx.setDataOptions(gmds::N|gmds::F);

			std::stringstream vtk_fileSingx;
			vtk_fileSingx <<m_output_directory_name<<"-!!!surf_Sing_lines" <<i<<".vtk";
			vtkWriterSingx.write(vtk_fileSingx.str());
		}
	}
	*/

	// free everything
	free(ia);
	free(ja);
	free(ar);

	glp_delete_prob(lp);

	auto t2 = Clock::now();

	cout << "optimization "
		<< std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
		 << " milliseconds" << std::endl;

	if(withGlobalComments)
		cout<<"end createSingularityLinesShortestPaths"<<endl;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------

void SingularityGraphBuilder2D::getTraversedTrisFaceNeighbyVertsExhaustive(vector<vector<vector<gmds::math::Point>>>&                  line_discretization,
                                                                           vector<vector<vector<gmds::TCellID>>>&                      line_triangles,
                                                                           vector<vector<vector<gmds::TCellID>>>&                      finalPaths,
                                                                           bool&                                                       isTargetBdry,
                                                                           vector<gmds::TCellID>&                                      traversedTriangles,
                                                                           vector<pair<vector<unsigned int>, vector<unsigned int >>>&  isTraversedNodeCompVector,
                                                                           vector<pair<vector<unsigned int>, vector<unsigned int >>>&  isTraversedFaceCompVector,
                                                                           vector<pair<vector<unsigned int>, vector<unsigned int >>>&  isTraversedFaceCompVector_SegmentPathCode,
                                                                           vector<pair<vector<unsigned int>, vector<unsigned int >>>&  isTraversedFaceCompVector_SegmentPathCont,
                                                                           unsigned int&                                               contMatrix,
                                                                           unsigned int&                                               totalNumberOfVariables, //totalNumberOfSlots+1
                                                                           unsigned int&                                               contSource,
                                                                           unsigned int&                                               contTarget,
                                                                           vector<pair<unsigned int, unsigned int>>&                   IllegalCross,
                                                                           vector<pair<unsigned int, unsigned int>>&                   contMatrixToSingularity,
                                                                           gmds::TCellID&                                              lastVisTri,
                                                                           vector<gmds::math::Point>&                                  finalEndPoint,
                                                                           unsigned int&                                               i_orig,
                                                                           unsigned int&                                               j_orig,
                                                                           unsigned int&                                               t1_orig,
                                                                           unsigned int&                                               t2_orig)
{
	/* in getShortestPathBtwFacesOptimized we obtain traversed triangles that that have at least one vertex
	* common; therfore the list of traversedTriangles is not complete (we might have some missing
	* triangles in between)
	* This function detects the entire set of traversed triangles but also marks each visited triangle
	* with the cross component direction that the shortest path (encoded by contMatrix) follows
	* -> needed for the detection of IllegalCross*/

	//WARNING TODO problem: 45 degrees -> 2 lines might not be found as illegal
	if(withGlobalComments)
		cout<<"getTraversedTrisFaceNeighbyVertsExhaustive contMatrix "<<contMatrix<<endl;

	vector<unsigned int> currentPath(finalPaths[contSource][contTarget].begin(), finalPaths[contSource][contTarget].end());
	//vector<gmds::TCellID> traversedTriangles(traversedSPTris[contSource][contTarget].begin(), traversedSPTris[contSource][contTarget].end());
	vector<gmds::math::Point> discrPtsSource(line_discretization[i_orig][j_orig].begin(), line_discretization[i_orig][j_orig].end());
	vector<gmds::math::Point> discrPtsTarget;
	vector<gmds::TCellID> travTriSource(line_triangles[i_orig][j_orig].begin(), line_triangles[i_orig][j_orig].end());
	vector<gmds::TCellID> travTriTarget;
	//discrPtsSource.insert(discrPtsSource.end(), line_discretization[i_orig][j_orig].begin(), line_discretization[i_orig][j_orig].end());
	//travTriSource.insert(travTriSource.end(), line_triangles[i_orig][j_orig].begin(), line_triangles[i_orig][j_orig].end());

	if(!isTargetBdry){
		//std::reverse_copy(line_discretization[t1_orig][t2_orig].begin(), line_discretization[t1_orig][t2_orig].end(), std::back_inserter(discrPtsTarget));
		discrPtsTarget.insert(discrPtsTarget.end(), line_discretization[t1_orig][t2_orig].rbegin(), line_discretization[t1_orig][t2_orig].rend());

		//std::reverse_copy(line_discretization[t1_orig][t2_orig].begin(), line_discretization[t1_orig][t2_orig].end(), std::back_inserter(discrPtsTarget));
		travTriTarget.insert(travTriTarget.end(), line_triangles[t1_orig][t2_orig].rbegin(), line_triangles[t1_orig][t2_orig].rend());
	}

	bool withComments = false;
/*	if(contMatrix==-1){
		withComments = true;
	}
*/
	gmds::math::Point lastSlotPointTarget;
	gmds::TCellID lastSlotFaceTarget;

	if(withComments){
		cout<<"isTargetBdry "<<isTargetBdry<<endl;
		cout<<"travTriSource "<<endl;
		for(unsigned int i=0; i<travTriSource.size();i++)
			cout<<"travTriSource["<<i<<"] "<<travTriSource[i]<<endl;

		cout<<"discrPtsSource "<<endl;
		for(unsigned int i=0; i<discrPtsSource.size();i++)
			cout<<"discrPtsSource["<<i<<"] "<<discrPtsSource[i]<<endl;
		cout<<"currentPath "<<endl;
		for(unsigned int i=0; i<currentPath.size();i++)
			cout<<"currentPath["<<i<<"] "<<currentPath[i]<<endl;
	}

	if(!isTargetBdry){
		//lastSlotPointTarget = discrPtsTarget[0];   	//lastSlotFaceTarget = travTriTarget[0];
		if(withComments){
			for(unsigned int i=0; i<travTriTarget.size();i++)
				cout<<"travTriTarget["<<i<<"] "<<travTriTarget[i]<<endl;
			for(unsigned int i=0; i<discrPtsTarget.size();i++)
				cout<<"discrPtsTarget["<<i<<"] "<<discrPtsTarget[i]<<endl;
		}
	}

	gmds::TCellID lastSlotFaceSource = travTriSource.back();

	gmds::math::Point currentPoint = discrPtsSource[0];

	gmds::math::Point nextPoint;
	gmds::TCellID nextFace, prevFace = travTriSource[max(0,(int)travTriSource.size()-2)];

	gmds::math::Point intersectionPoint;
	double intersectionParam;

	vector<bool> visitedFaces(original_faces_number, false);
	vector<bool> visitedNodes(original_nodes_number, false);
	gmds::TCellID currentFace_id = travTriSource[0];
	gmds::Face currentFace;

	if(discrPtsSource.size()>2){
		if(withComments){
			cout<<"if(discrPtsSource.size()>2){"<<endl;
		}
		nextFace = currentFace_id;
		for(unsigned int i=1; i<discrPtsSource.size(); i++){
			math::Segment from_seg(discrPtsSource[i-1], discrPtsSource[i]);
			math::Vector3d tempDir(discrPtsSource[i-1], discrPtsSource[i]);
			for(unsigned int j=1; j<travTriSource.size(); j++){
				currentFace_id = nextFace;
				nextFace = travTriSource[j];

				bool foundNextFace = false;

				if(withComments){
					cout<<"currentFace_id "<<currentFace_id<<endl;
					cout<<"nextFace "<<nextFace<<endl;
					cout<<"prevFace "<<prevFace<<endl;
					cout<<"from_seg "<<discrPtsSource[i-1]<<" -> "<<discrPtsSource[i]<<endl;
				}
				vector<gmds::Edge> adj_edges;
				vector<gmds::Node> temp_nodes = (m_mesh->get<gmds::Face>(currentFace_id)).get<gmds::Node>();
				for(unsigned int j1=0; j1<3; j1++){
					vector<gmds::Edge> temp_edges = temp_nodes[j1].get<gmds::Edge>();
					adj_edges.insert(adj_edges.end(), temp_edges.begin(), temp_edges.end());
				}
				for(unsigned int j1=0; j1<adj_edges.size(); j1++){
					vector<gmds::Node> currentNodes =  adj_edges[j1].get<gmds::Node>();
					math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
					if(withComments){
						cout<<"edge number j = "<<j<<" between nodes "<<currentNodes[0].id()<<" & "<<currentNodes[1]<<endl;
						from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon);

						cout<<"from_seg.SecondMetIntersect2D "<<from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon)
						<<" ; intersectionParam "<<intersectionParam<<endl;
					}
					if(from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon)){

						vector<gmds::Face> edgeAdjFaces = adj_edges[j1].get<gmds::Face>();
						for(unsigned int k=0; k<edgeAdjFaces.size(); k++){
							if(!visitedFaces[edgeAdjFaces[k].id()]){
								currentFace_id = edgeAdjFaces[k].id();
								if(currentFace_id==nextFace)
									foundNextFace = true;

								math::Cross2D crossV = triangle_centers_cross[currentFace_id];
								unsigned int compVId = crossV.closestComponentVectorAsIndex(tempDir)%2;

								if(compVId==0){
									for(unsigned int j=0; j<isTraversedFaceCompVector[currentFace_id].first.size(); j++){
										unsigned int tempCS = isTraversedFaceCompVector[currentFace_id].first[j];
										if(tempCS!=contMatrix){
											if((contMatrixToSingularity[contMatrix].first!=contMatrixToSingularity[tempCS].first)&&
											((contMatrixToSingularity[contMatrix].second!=contMatrixToSingularity[tempCS].second)||
											(fmod(contMatrix, totalNumberOfVariables)==totalNumberOfVariables-1))){
												if((( (int)(contMatrix/totalNumberOfVariables)!=((int)(tempCS/totalNumberOfVariables)))&&
												( ( fmod(contMatrix, totalNumberOfVariables))!=(fmod(tempCS, totalNumberOfVariables)))&&
												( (int)(contMatrix/totalNumberOfVariables)!=(fmod(tempCS, totalNumberOfVariables)))&&
												( (fmod(contMatrix, totalNumberOfVariables))!=((int)(tempCS/totalNumberOfVariables)) ) )||
												(( ( fmod(contMatrix, totalNumberOfVariables))==(fmod(tempCS, totalNumberOfVariables)))&&(( fmod(contMatrix, totalNumberOfVariables))==totalNumberOfVariables-1))){
													bool toAdd = true;
													for(unsigned int k=0; k<IllegalCross.size(); k++){
														if(((contMatrix==IllegalCross[k].second)&&(tempCS==IllegalCross[k].first))||
														((tempCS==IllegalCross[k].second)&&(contMatrix==IllegalCross[k].first))){
															toAdd = false;
															break;
														}
													}
													if(toAdd){
														math::Point tempPt1, tempPt2;
														unsigned int pairContSource = (int)(tempCS/totalNumberOfVariables);
														unsigned int pairContTarget = fmod(tempCS,totalNumberOfVariables);
														unsigned int tempS_i = (int)(pairContSource/5);// source sing point
														unsigned int tempS_j = fmod(pairContSource,5);// source slot no
														unsigned int tempT_i = (int)(pairContTarget/5);// target sing point
														unsigned int tempT_j = fmod(pairContTarget,5);// target slot no
														if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==0){
															tempPt1 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]];
															tempPt2 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1];
														}
														else{
															if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==1){
																tempPt1 = line_discretization[tempS_i][tempS_j].back();
																if(finalPaths[pairContSource][pairContTarget].size()!=0)
																	tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][0]];
																else
																	tempPt2 = line_discretization[tempT_i][tempT_j].back();
															}
															else{
																if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==2){
																	tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]]];
																	tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1]];
																}
																else{
																	if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==3){
																		if(finalPaths[pairContSource][pairContTarget].size()!=0)
																			tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																		else
																		tempPt1 = line_discretization[tempS_i][tempS_j].back();
																		tempPt2 = line_discretization[tempT_i][tempT_j].back();// !!!it is the inverse
																	}
																	else{
																		if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==4){
																		 	tempPt1 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]];
																		 	tempPt2 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1];// !!!it is the inverse
																		}
																		else{
																			if(finalPaths[pairContSource][pairContTarget].size()!=0)
																				tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																			else
																				tempPt1 = line_discretization[tempS_i][tempS_j].back();
																			tempPt2 = finalEndPoint[pairContSource];
																		}
																	}
																}
															}
														}

														math::Segment segment1(tempPt1, tempPt2);
														if(from_seg.SecondMetIntersect2D(segment1, intersectionPoint, intersectionParam, temp_epsilon))
															IllegalCross.push_back(make_pair(tempCS, contMatrix));
														//perhaps here also check if their direction corresponds...however they folow the same cross component
													}
												}
											}
										}
									}
									isTraversedFaceCompVector[currentFace_id].first.push_back(contMatrix);
									isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first.push_back(0);
									isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first.push_back(i-1);
								}
								else{//compV ==1
									for(unsigned int j=0; j<isTraversedFaceCompVector[currentFace_id].second.size(); j++){
										unsigned int tempCS = isTraversedFaceCompVector[currentFace_id].second[j];
										if(tempCS!=contMatrix){
											if((contMatrixToSingularity[contMatrix].first!=contMatrixToSingularity[tempCS].first)&&
											((contMatrixToSingularity[contMatrix].second!=contMatrixToSingularity[tempCS].second)||
											(fmod(contMatrix, totalNumberOfVariables)==totalNumberOfVariables-1))){
												if((( (int)(contMatrix/totalNumberOfVariables)!=((int)(tempCS/totalNumberOfVariables)))&&
												( ( fmod(contMatrix, totalNumberOfVariables))!=(fmod(tempCS, totalNumberOfVariables)))&&
												( (int)(contMatrix/totalNumberOfVariables)!=(fmod(tempCS, totalNumberOfVariables)))&&
												( (fmod(contMatrix, totalNumberOfVariables))!=((int)(tempCS/totalNumberOfVariables)) ) )||
												(( ( fmod(contMatrix, totalNumberOfVariables))==(fmod(tempCS, totalNumberOfVariables)))&&(( fmod(contMatrix, totalNumberOfVariables))==totalNumberOfVariables-1))){
													bool toAdd = true;
													for(unsigned int k=0; k<IllegalCross.size(); k++){
														if(((contMatrix==IllegalCross[k].second)&&(tempCS==IllegalCross[k].first))||
														((tempCS==IllegalCross[k].second)&&(contMatrix==IllegalCross[k].first))){
															toAdd = false;
															break;
														}
													}
													if(toAdd){
														math::Point tempPt1, tempPt2;
														unsigned int pairContSource = (int)(tempCS/totalNumberOfVariables);
														unsigned int pairContTarget = fmod(tempCS,totalNumberOfVariables);
														unsigned int tempS_i = (int)(pairContSource/5);// source sing point
														unsigned int tempS_j = fmod(pairContSource,5);// source slot no
														unsigned int tempT_i = (int)(pairContTarget/5);// target sing point
														unsigned int tempT_j = fmod(pairContTarget,5);// target slot no
														if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==0){
															tempPt1 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]];
															tempPt2 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1];
														}
														else{
															if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==1){
																tempPt1 = line_discretization[tempS_i][tempS_j].back();
																if(finalPaths[pairContSource][pairContTarget].size()!=0)
																	tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][0]];
																else
																	tempPt2 = line_discretization[tempT_i][tempT_j].back();
															}
															else{
																if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==2){
																	tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]]];
																	tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1]];
																}
																else{
																	if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==3){
																		if(finalPaths[pairContSource][pairContTarget].size()!=0)
																			tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																		else
																			tempPt1 = line_discretization[tempS_i][tempS_j].back();
																		tempPt2 = line_discretization[tempT_i][tempT_j].back();// !!!it is the inverse
																	}
																	else{
																		if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==4){
																			tempPt1 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]];
																			tempPt2 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1];// !!!it is the inverse
																		}
																		else{
																			if(finalPaths[pairContSource][pairContTarget].size()!=0)
																				tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																			else
																				tempPt1 = line_discretization[tempS_i][tempS_j].back();
																			tempPt2 = finalEndPoint[pairContSource];
																		}
																	}
																}
															}
														}

														math::Segment segment1(tempPt1, tempPt2);
														if(from_seg.SecondMetIntersect2D(segment1, intersectionPoint, intersectionParam, temp_epsilon))
															IllegalCross.push_back(make_pair(tempCS, contMatrix));
															//perhaps here also check if their direction corresponds...however they folow the cross component

													}
												}
											}
										}
									}
									isTraversedFaceCompVector[currentFace_id].second.push_back(contMatrix);
									isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second.push_back(0);
									isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second.push_back(i-1);
								}

								visitedFaces[edgeAdjFaces[k].id()] = true;
								traversedTriangles.push_back(edgeAdjFaces[k].id());
							}
						}
					}
				}
				if((!foundNextFace)&&(withComments))
					cout<<"!foundNextFace0"<<endl;
			}
		}
	}

	nextPoint = discrPtsSource.back();
	if(withComments){
		cout<<"currentFace_id "<<currentFace_id<<endl;
		cout<<"currentPoint "<<currentPoint<<endl;
		cout<<"traversedTriangles.size() "<<traversedTriangles.size()<<endl;
		cout<<"prevFace "<<prevFace<<endl;
	}
	traversedTriangles.push_back(currentFace_id);
	visitedFaces[currentFace_id] = true;

	unsigned int currentPathSize = currentPath.size();
	//if(((fmod(contMatrix, totalNumberOfVariables))!=(totalNumberOfVariables-1))&&(currentPathSize>1))
	//	currentPathSize--;
	if(withComments){
		cout<<"lastSlotFaceSource "<<lastSlotFaceSource<<endl;
		cout<<"currentFace_id "<<currentFace_id<<endl;
		cout<<"!!!currentPathSize "<<currentPathSize<<endl;
	}
	nextFace = travTriSource.back();
	for(unsigned int i=0; i<currentPathSize; i++){
		currentFace_id = nextFace;
		currentPoint = nextPoint;
		nextFace = currentPath[i];
		nextPoint = triangle_centers[nextFace];

		math::Segment from_seg(currentPoint, nextPoint);
		math::Vector3d tempDir(currentPoint, nextPoint);

		if(withComments){
			cout<<"currentPath["<<i<<"] "<<currentPath[i]<<endl;
		}

		bool foundNextFace = false;

		if(withComments){
			cout<<"currentFace_id "<<currentFace_id<<endl;
			cout<<"nextFace "<<nextFace<<endl;
			cout<<"prevFace "<<prevFace<<endl;
			cout<<"from_seg "<<currentPoint<<" -> "<<nextPoint<<endl;
		}
		vector<gmds::Edge> adj_edges;
		vector<gmds::Node> temp_nodes = (m_mesh->get<gmds::Face>(currentFace_id)).get<gmds::Node>();
		for(unsigned int j1=0; j1<3; j1++){
			vector<gmds::Edge> temp_edges = temp_nodes[j1].get<gmds::Edge>();
			adj_edges.insert(adj_edges.end(), temp_edges.begin(), temp_edges.end());
		}
		for(unsigned int j1=0; j1<adj_edges.size(); j1++){
			vector<gmds::Node> currentNodes =  adj_edges[j1].get<gmds::Node>();
			math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
			if(withComments){
				cout<<"edge number j1 = "<<j1<<" between nodes "<<currentNodes[0].id()<<" & "<<currentNodes[1]<<endl;
				from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon);

				cout<<"from_seg.SecondMetIntersect2D "<<from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon)
				<<" ; intersectionParam "<<intersectionParam<<endl;

			}
			if(from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon)){

				vector<gmds::Face> edgeAdjFaces = adj_edges[j1].get<gmds::Face>();
				for(unsigned int k=0; k<edgeAdjFaces.size(); k++){
					if(!visitedFaces[edgeAdjFaces[k].id()]){
						currentFace_id = edgeAdjFaces[k].id();
						if(currentFace_id==nextFace)
							foundNextFace = true;
						math::Cross2D crossV = triangle_centers_cross[currentFace_id];
						unsigned int compVId = crossV.closestComponentVectorAsIndex(tempDir)%2;
						if(compVId==0){
							for(unsigned int j=0; j<isTraversedFaceCompVector[currentFace_id].first.size(); j++){
								unsigned int tempCS = isTraversedFaceCompVector[currentFace_id].first[j];
								if(tempCS!=contMatrix){
									if((contMatrixToSingularity[contMatrix].first!=contMatrixToSingularity[tempCS].first)&&
									((contMatrixToSingularity[contMatrix].second!=contMatrixToSingularity[tempCS].second)||
									(fmod(contMatrix, totalNumberOfVariables)==totalNumberOfVariables-1))){
										if((( (int)(contMatrix/totalNumberOfVariables)!=((int)(tempCS/totalNumberOfVariables)))&&
										( ( fmod(contMatrix, totalNumberOfVariables))!=(fmod(tempCS, totalNumberOfVariables)))&&
										( (int)(contMatrix/totalNumberOfVariables)!=(fmod(tempCS, totalNumberOfVariables)))&&
										( (fmod(contMatrix, totalNumberOfVariables))!=((int)(tempCS/totalNumberOfVariables)) ) )||
										(( ( fmod(contMatrix, totalNumberOfVariables))==(fmod(tempCS, totalNumberOfVariables)))&&(( fmod(contMatrix, totalNumberOfVariables))==totalNumberOfVariables-1))){
											bool toAdd = true;
											for(unsigned int k=0; k<IllegalCross.size(); k++){
												if(((contMatrix==IllegalCross[k].second)&&(tempCS==IllegalCross[k].first))||
												((tempCS==IllegalCross[k].second)&&(contMatrix==IllegalCross[k].first))){
													toAdd = false;
													break;
												}
											}
											if(toAdd){
												if(withComments){
													cout<<"j= "<<j<<endl;
													cout<<"isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first.size() "<<
													isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first.size()<<endl;
													cout<<"isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j] "<<
													isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]<<endl;
													cout<<"isTraversedFaceCompVector[currentFace_id].first[j] "<<
													isTraversedFaceCompVector[currentFace_id].first[j]<<endl;
												}
												math::Point tempPt1, tempPt2;
												unsigned int pairContSource = (int)(tempCS/totalNumberOfVariables);
												unsigned int pairContTarget = fmod(tempCS,totalNumberOfVariables);
												unsigned int tempS_i = (int)(pairContSource/5);// source sing point
												unsigned int tempS_j = fmod(pairContSource,5);// source slot no
												unsigned int tempT_i = (int)(pairContTarget/5);// target sing point
												unsigned int tempT_j = fmod(pairContTarget,5);// target slot no
												if(withComments){
													cout<<"tempS_i "<<tempS_i<<endl;
												 cout<<"tempS_j "<<tempS_j<<endl;
												 cout<<"tempT_i "<<tempT_i<<endl;
												 cout<<"tempT_j "<<tempT_j<<endl;

												}
												if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==0){
													tempPt1 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]];
													tempPt2 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1];
												}
												else{
													if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==1){
														tempPt1 = line_discretization[tempS_i][tempS_j].back();
														if(finalPaths[pairContSource][pairContTarget].size()!=0)
															tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][0]];
														else
															tempPt2 = line_discretization[tempT_i][tempT_j].back();
													}
													else{
														if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==2){
															tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]]];
															tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1]];
														}
														else{
															if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==3){

																if(finalPaths[pairContSource][pairContTarget].size()!=0)
																	tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																else
																	tempPt1 = line_discretization[tempS_i][tempS_j].back();

																tempPt2 = line_discretization[tempT_i][tempT_j].back();// !!!it is the inverse
															}
															else{
																if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==4){
																	tempPt1 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]];
																	tempPt2 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1];// !!!it is the inverse
																}
																else{
																	if(finalPaths[pairContSource][pairContTarget].size()!=0)
																		tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																	else
																		tempPt1 = line_discretization[tempS_i][tempS_j].back();
																	tempPt2 = finalEndPoint[pairContSource];
																}
															}
														}
													}
												}
												math::Segment segment1(tempPt1, tempPt2);
												from_seg.SecondMetIntersect2D(segment1, intersectionPoint, intersectionParam, temp_epsilon);
												if(from_seg.SecondMetIntersect2D(segment1, intersectionPoint, intersectionParam, temp_epsilon))
													IllegalCross.push_back(make_pair(tempCS, contMatrix));
													 //perhaps here also check if their direction corresponds...however they folow the cross component

											}
										}
									}
						 		}
							}
							isTraversedFaceCompVector[currentFace_id].first.push_back(contMatrix);
							if(i==0){
								isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first.push_back(1);
								isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first.push_back(0);
							}
							else{
								isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first.push_back(2);
								isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first.push_back(i-1);
							}
						}
						else{//compVId==1
							for(unsigned int j=0; j<isTraversedFaceCompVector[currentFace_id].second.size(); j++){
								unsigned int tempCS = isTraversedFaceCompVector[currentFace_id].second[j];
								if(tempCS!=contMatrix){
									if((contMatrixToSingularity[contMatrix].first!=contMatrixToSingularity[tempCS].first)&&
									((contMatrixToSingularity[contMatrix].second!=contMatrixToSingularity[tempCS].second)||
									(fmod(contMatrix, totalNumberOfVariables)==totalNumberOfVariables-1))){
										if((( (int)(contMatrix/totalNumberOfVariables)!=((int)(tempCS/totalNumberOfVariables)))&&
										( ( fmod(contMatrix, totalNumberOfVariables))!=(fmod(tempCS, totalNumberOfVariables)))&&
										( (int)(contMatrix/totalNumberOfVariables)!=(fmod(tempCS, totalNumberOfVariables)))&&
										( (fmod(contMatrix, totalNumberOfVariables))!=((int)(tempCS/totalNumberOfVariables)) ) )||
										(( ( fmod(contMatrix, totalNumberOfVariables))==(fmod(tempCS, totalNumberOfVariables)))&&(( fmod(contMatrix, totalNumberOfVariables))==totalNumberOfVariables-1))){
											bool toAdd = true;
											for(unsigned int k=0; k<IllegalCross.size(); k++){
												if(((contMatrix==IllegalCross[k].second)&&(tempCS==IllegalCross[k].first))||
												((tempCS==IllegalCross[k].second)&&(contMatrix==IllegalCross[k].first))){
													toAdd = false;
													break;
												}
											}
											if(toAdd){
												if(withComments){
													cout<<"j= "<<j<<endl;
													cout<<"isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second.size() "<<
													isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second.size()<<endl;
													cout<<"isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j] "<<
													isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]<<endl;
												}
												math::Point tempPt1, tempPt2;
												unsigned int pairContSource = (int)(tempCS/totalNumberOfVariables);
												unsigned int pairContTarget = fmod(tempCS,totalNumberOfVariables);
												unsigned int tempS_i = (int)(pairContSource/5);// source sing point
												unsigned int tempS_j = fmod(pairContSource,5);// source slot no
												unsigned int tempT_i = (int)(pairContTarget/5);// target sing point
												unsigned int tempT_j = fmod(pairContTarget,5);// target slot no

												if(withComments){
													cout<<"tempS_i "<<tempS_i<<endl;
													cout<<"tempS_j "<<tempS_j<<endl;
													cout<<"tempT_i "<<tempT_i<<endl;
													cout<<"tempT_j "<<tempT_j<<endl;
												}
												if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==0){
													tempPt1 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]];
													tempPt2 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1];
												}
												else{
													if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==1){
														tempPt1 = line_discretization[tempS_i][tempS_j].back();
														if(finalPaths[pairContSource][pairContTarget].size()!=0)
														 	tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][0]];
														else
														 	tempPt2 = line_discretization[tempT_i][tempT_j].back();
													}
													else{
														if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==2){
															tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]]];
															tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1]];
														}
														else{
															if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==3){
																if(finalPaths[pairContSource][pairContTarget].size()!=0)
																	tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																else
																	tempPt1 = line_discretization[tempS_i][tempS_j].back();

																tempPt2 = line_discretization[tempT_i][tempT_j].back();// !!!it is the inverse
															}
															else{
																if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==4){
																	tempPt1 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]];
																	tempPt2 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1];// !!!it is the inverse
																}
																else{
																	if(finalPaths[pairContSource][pairContTarget].size()!=0)
																		tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																	else
																		tempPt1 = line_discretization[tempS_i][tempS_j].back();
																	tempPt2 = finalEndPoint[pairContSource];
																}
															}
														}
													}
												}
												math::Segment segment1(tempPt1, tempPt2);
												if(from_seg.SecondMetIntersect2D(segment1, intersectionPoint, intersectionParam, temp_epsilon))
													IllegalCross.push_back(make_pair(tempCS, contMatrix));
												//perhaps here also check if their direction corresponds...however they folow the cross component
											}
										}
									}
								}
							}
							isTraversedFaceCompVector[currentFace_id].second.push_back(contMatrix);
							if(i==0){
								isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second.push_back(1);
								isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second.push_back(0);
							}
							else{
								isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second.push_back(2);
								isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second.push_back(i-1);
							}
						}

						visitedFaces[edgeAdjFaces[k].id()] = true;
						traversedTriangles.push_back(edgeAdjFaces[k].id());
					}
				}
			}
		}
		if((!foundNextFace)&&(withComments))
			cout<<"!foundNextFace1"<<endl;
	}

	//nextFace = currentPath[currentPathSize-1];
	currentPoint = nextPoint;//triangle_centers[nextFace];
	if(withComments){
		cout<<"before check if target bdry or not currentFace_id "<<currentFace_id<<endl;
		cout<<"currentPoint "<<currentPoint<<endl;
	}
	if(!isTargetBdry){

		if((!travTriTarget.empty())&&(!discrPtsTarget.empty())){
			if(withComments)
				cout<<"if((!travTriTarget.empty())&&(!discrPtsTarget.empty())){"<<endl;
			currentFace_id = nextFace;
			nextFace = travTriTarget[0];
			nextPoint = discrPtsTarget[0];
			math::Segment from_seg(currentPoint, nextPoint);
			math::Vector3d tempDir(currentPoint, nextPoint);

			bool foundNextFace = false;

			if(withComments){
				cout<<"currentFace_id "<<currentFace_id<<endl;
				cout<<"nextFace "<<nextFace<<endl;
				cout<<"prevFace "<<prevFace<<endl;
			}
			vector<gmds::Edge> adj_edges;
			vector<gmds::Node> temp_nodes = (m_mesh->get<gmds::Face>(currentFace_id)).get<gmds::Node>();
			for(unsigned int j1=0; j1<3; j1++){
				vector<gmds::Edge> temp_edges = temp_nodes[j1].get<gmds::Edge>();
				adj_edges.insert(adj_edges.end(), temp_edges.begin(), temp_edges.end());
			}
			for(unsigned int j1=0; j1<adj_edges.size(); j1++){
				vector<gmds::Node> currentNodes =  adj_edges[j1].get<gmds::Node>();
				math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
				if(withComments){
					cout<<"edge number j = "<<j1<<" between nodes "<<currentNodes[0].id()<<" & "<<currentNodes[1]<<endl;
					from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon);

					cout<<"from_seg.SecondMetIntersect2D "<<from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon)
					<<endl;
					cout<<"temp_epsilon "<<temp_epsilon<<endl;
					cout<<"intersectionParam "<<intersectionParam<<endl;
				}
				if(from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon)){

					//WARNING if we want with nodes, here should be inserted

					vector<gmds::Face> edgeAdjFaces = adj_edges[j1].get<gmds::Face>();
					for(unsigned int k=0; k<edgeAdjFaces.size(); k++){
						if(!visitedFaces[edgeAdjFaces[k].id()]){
							currentFace_id = edgeAdjFaces[k].id();
							if(currentFace_id==nextFace)
								foundNextFace = true;

							math::Cross2D crossV = triangle_centers_cross[currentFace_id];
							unsigned int compVId = crossV.closestComponentVectorAsIndex(tempDir)%2;

							if(compVId==0){
								for(unsigned int j=0; j<isTraversedFaceCompVector[currentFace_id].first.size(); j++){
									unsigned int tempCS = isTraversedFaceCompVector[currentFace_id].first[j];
									if(tempCS!=contMatrix){
										if((contMatrixToSingularity[contMatrix].first!=contMatrixToSingularity[tempCS].first)&&
										((contMatrixToSingularity[contMatrix].second!=contMatrixToSingularity[tempCS].second)||
										(fmod(contMatrix, totalNumberOfVariables)==totalNumberOfVariables-1))){
											if((( (int)(contMatrix/totalNumberOfVariables)!=((int)(tempCS/totalNumberOfVariables)))&&
											( ( fmod(contMatrix, totalNumberOfVariables))!=(fmod(tempCS, totalNumberOfVariables)))&&
											( (int)(contMatrix/totalNumberOfVariables)!=(fmod(tempCS, totalNumberOfVariables)))&&
											( (fmod(contMatrix, totalNumberOfVariables))!=((int)(tempCS/totalNumberOfVariables)) ) )||
											(( ( fmod(contMatrix, totalNumberOfVariables))==(fmod(tempCS, totalNumberOfVariables)))&&(( fmod(contMatrix, totalNumberOfVariables))==totalNumberOfVariables-1))){
												bool toAdd = true;
												for(unsigned int k=0; k<IllegalCross.size(); k++){
													if(((contMatrix==IllegalCross[k].second)&&(tempCS==IllegalCross[k].first))||
													((tempCS==IllegalCross[k].second)&&(contMatrix==IllegalCross[k].first))){
														toAdd = false;
														break;
													}
												}
												if(toAdd){
													math::Point tempPt1, tempPt2;
													unsigned int pairContSource = (int)(tempCS/totalNumberOfVariables);
													unsigned int pairContTarget = fmod(tempCS,totalNumberOfVariables);
													unsigned int tempS_i = (int)(pairContSource/5);// source sing point
													unsigned int tempS_j = fmod(pairContSource,5);// source slot no
													unsigned int tempT_i = (int)(pairContTarget/5);// target sing point
													unsigned int tempT_j = fmod(pairContTarget,5);// target slot no
													if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==0){
														tempPt1 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]];
														tempPt2 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1];
													}
													else{
														if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==1){
															tempPt1 = line_discretization[tempS_i][tempS_j].back();
															if(finalPaths[pairContSource][pairContTarget].size()!=0)
															 	tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][0]];
															else
															 	tempPt2 = line_discretization[tempT_i][tempT_j].back();
														}
														else{
															if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==2){
															 	tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]]];
															 	tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1]];
															}
															else{
																if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==3){

																	if(finalPaths[pairContSource][pairContTarget].size()!=0)
																		tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																	else
																		tempPt1 = line_discretization[tempS_i][tempS_j].back();
																	tempPt2 = line_discretization[tempT_i][tempT_j].back();// !!!it is the inverse
																}
																else{
																	if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==4){
																	 	tempPt1 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]];
																	 	tempPt2 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1];// !!!it is the inverse
																	}
																	else{
																		if(finalPaths[pairContSource][pairContTarget].size()!=0)
																			tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																		else
																			tempPt1 = line_discretization[tempS_i][tempS_j].back();
																		tempPt2 = finalEndPoint[pairContSource];
																	}
																}
															}
														}
													}

													math::Segment segment1(tempPt1, tempPt2);
													if(from_seg.SecondMetIntersect2D(segment1, intersectionPoint, intersectionParam, temp_epsilon))
														IllegalCross.push_back(make_pair(tempCS, contMatrix));
													//perhaps here also check if their direction corresponds...however they folow the cross component
												}
											}
										}
									}
								}
								isTraversedFaceCompVector[currentFace_id].first.push_back(contMatrix);
								isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first.push_back(3);
								isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first.push_back(0);
							}
							else{//compVId == 1
								for(unsigned int j=0; j<isTraversedFaceCompVector[currentFace_id].second.size(); j++){
									unsigned int tempCS = isTraversedFaceCompVector[currentFace_id].second[j];
									if(tempCS!=contMatrix){
										if((contMatrixToSingularity[contMatrix].first!=contMatrixToSingularity[tempCS].first)&&
										((contMatrixToSingularity[contMatrix].second!=contMatrixToSingularity[tempCS].second)||
										(fmod(contMatrix, totalNumberOfVariables)==totalNumberOfVariables-1))){
											if((( (int)(contMatrix/totalNumberOfVariables)!=((int)(tempCS/totalNumberOfVariables)))&&
											( ( fmod(contMatrix, totalNumberOfVariables))!=(fmod(tempCS, totalNumberOfVariables)))&&
											( (int)(contMatrix/totalNumberOfVariables)!=(fmod(tempCS, totalNumberOfVariables)))&&
											( (fmod(contMatrix, totalNumberOfVariables))!=((int)(tempCS/totalNumberOfVariables)) ) )||
											(( ( fmod(contMatrix, totalNumberOfVariables))==(fmod(tempCS, totalNumberOfVariables)))&&(( fmod(contMatrix, totalNumberOfVariables))==totalNumberOfVariables-1))){
												bool toAdd = true;
												for(unsigned int k=0; k<IllegalCross.size(); k++){
													if(((contMatrix==IllegalCross[k].second)&&(tempCS==IllegalCross[k].first))||
													((tempCS==IllegalCross[k].second)&&(contMatrix==IllegalCross[k].first))){
														toAdd = false;
														break;
													}
												}
												if(toAdd){
													math::Point tempPt1, tempPt2;
													unsigned int pairContSource = (int)(tempCS/totalNumberOfVariables);
													unsigned int pairContTarget = fmod(tempCS,totalNumberOfVariables);
													unsigned int tempS_i = (int)(pairContSource/5);// source sing point
													unsigned int tempS_j = fmod(pairContSource,5);// source slot no
													unsigned int tempT_i = (int)(pairContTarget/5);// target sing point
													unsigned int tempT_j = fmod(pairContTarget,5);// target slot no
													if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==0){
														tempPt1 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]];
														tempPt2 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1];
													}
													else{
														if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==1){
															tempPt1 = line_discretization[tempS_i][tempS_j].back();
															if(finalPaths[pairContSource][pairContTarget].size()!=0)
																tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][0]];
															else
																tempPt2 = line_discretization[tempT_i][tempT_j].back();
														}
														else{
															if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==2){
																tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]]];
																tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1]];
															}
															else{
																if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==3){
																	if(finalPaths[pairContSource][pairContTarget].size()!=0)
																		tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																	else
																		tempPt1 = line_discretization[tempS_i][tempS_j].back();
																	tempPt2 = line_discretization[tempT_i][tempT_j].back();// !!!it is the inverse
																}
																else{
																	if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==4){
																		tempPt1 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]];
																		tempPt2 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1];// !!!it is the inverse
																	}
																	else{
																		if(finalPaths[pairContSource][pairContTarget].size()!=0)
																			tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																		else
																			tempPt1 = line_discretization[tempS_i][tempS_j].back();
																		tempPt2 = finalEndPoint[pairContSource];
																	}
																}
															}
														}
													}

													math::Segment segment1(tempPt1, tempPt2);
													if(from_seg.SecondMetIntersect2D(segment1, intersectionPoint, intersectionParam, temp_epsilon))
														IllegalCross.push_back(make_pair(tempCS, contMatrix));
								 					//perhaps here also check if their direction corresponds...however they folow the cross component
												}
											}
										}
									}
								}
								isTraversedFaceCompVector[currentFace_id].second.push_back(contMatrix);
								isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second.push_back(3);
								isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second.push_back(0);
							}

							visitedFaces[edgeAdjFaces[k].id()] = true;
							traversedTriangles.push_back(edgeAdjFaces[k].id());
						}
					}
				}
			}
			if((!foundNextFace)&&(withComments))
				cout<<"!foundNextFace2"<<endl;

			currentPoint = nextPoint;
		}

		if(withComments){
			cout<<"if(!isTargetBdry){"<<endl;
			cout<<"discrPtsTarget.size() "<<discrPtsTarget.size()<<endl;
			for(unsigned int i=1; i<discrPtsTarget.size(); i++){
				cout<<"discrPtsTarget["<<i<<"] "<<discrPtsTarget[i]<<endl;
			}
		}
		nextFace = travTriTarget[0];
		if(discrPtsTarget.size()>2){

			for(unsigned int i=1; i<discrPtsTarget.size(); i++){

				if(withComments)
					cout<<"from_seg "<<discrPtsTarget[i-1]<<" -> "<<discrPtsTarget[i]<<endl;
				math::Segment from_seg(discrPtsTarget[i-1], discrPtsTarget[i]);
				math::Vector3d tempDir(discrPtsTarget[i-1], discrPtsTarget[i]);
				for(unsigned int j=1; j<travTriTarget.size(); j++){// first one is the sing face
					if(withComments)
						cout<<"travTriTarget[ "<<j<<"] "<<travTriTarget[j]<<endl;
					currentFace_id = nextFace;
					nextFace = travTriTarget[j];

					bool foundNextFace = false;

					if(withComments){
						cout<<"currentFace_id "<<currentFace_id<<endl;
						cout<<"nextFace "<<nextFace<<endl;
						cout<<"prevFace "<<prevFace<<endl;
					}
					vector<gmds::Edge> adj_edges;
					vector<gmds::Node> temp_nodes = (m_mesh->get<gmds::Face>(currentFace_id)).get<gmds::Node>();
					for(unsigned int j1=0; j1<3; j1++){
						vector<gmds::Edge> temp_edges = temp_nodes[j1].get<gmds::Edge>();
						adj_edges.insert(adj_edges.end(), temp_edges.begin(), temp_edges.end());
					}
					for(unsigned int j1=0; j1<adj_edges.size(); j1++){
						vector<gmds::Node> currentNodes =  adj_edges[j1].get<gmds::Node>();
						math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
						if(withComments){
							cout<<"edge number j = "<<j<<" between nodes "<<currentNodes[0].id()<<" & "<<currentNodes[1]<<endl;
							from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon);

							cout<<"from_seg.SecondMetIntersect2D "<<from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon)
							<<endl;
							cout<<"temp_epsilon "<<temp_epsilon<<endl;
							cout<<"intersectionParam "<<intersectionParam<<endl;

						}
						if(from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon)){
							vector<gmds::Face> edgeAdjFaces = adj_edges[j1].get<gmds::Face>();
							for(unsigned int k=0; k<edgeAdjFaces.size(); k++){
								if(!visitedFaces[edgeAdjFaces[k].id()]){
									currentFace_id = edgeAdjFaces[k].id();
									if(currentFace_id==nextFace)
										foundNextFace = true;
									if(nextFace!=travTriTarget[j1]){
										math::Cross2D crossV = triangle_centers_cross[currentFace_id];
										unsigned int compVId = crossV.closestComponentVectorAsIndex(tempDir)%2;

										if(compVId==0){
											for(unsigned int j=0; j<isTraversedFaceCompVector[currentFace_id].first.size(); j++){
												unsigned int tempCS = isTraversedFaceCompVector[currentFace_id].first[j];
												if(tempCS!=contMatrix){
													if((contMatrixToSingularity[contMatrix].first!=contMatrixToSingularity[tempCS].first)&&
													((contMatrixToSingularity[contMatrix].second!=contMatrixToSingularity[tempCS].second)||
													(fmod(contMatrix, totalNumberOfVariables)==totalNumberOfVariables-1))){
														if((( (int)(contMatrix/totalNumberOfVariables)!=((int)(tempCS/totalNumberOfVariables)))&&
														( ( fmod(contMatrix, totalNumberOfVariables))!=(fmod(tempCS, totalNumberOfVariables)))&&
														( (int)(contMatrix/totalNumberOfVariables)!=(fmod(tempCS, totalNumberOfVariables)))&&
														( (fmod(contMatrix, totalNumberOfVariables))!=((int)(tempCS/totalNumberOfVariables)) ) )||
														(( ( fmod(contMatrix, totalNumberOfVariables))==(fmod(tempCS, totalNumberOfVariables)))&&(( fmod(contMatrix, totalNumberOfVariables))==totalNumberOfVariables-1))){
															bool toAdd = true;
															for(unsigned int k=0; k<IllegalCross.size(); k++){
																if(((contMatrix==IllegalCross[k].second)&&(tempCS==IllegalCross[k].first))||
																((tempCS==IllegalCross[k].second)&&(contMatrix==IllegalCross[k].first))){
																	toAdd = false;
																	break;
																}
															}
															if(toAdd){
																math::Point tempPt1, tempPt2;
																unsigned int pairContSource = (int)(tempCS/totalNumberOfVariables);
																unsigned int pairContTarget = fmod(tempCS,totalNumberOfVariables);
																unsigned int tempS_i = (int)(pairContSource/5);// source sing point
																unsigned int tempS_j = fmod(pairContSource,5);// source slot no
																unsigned int tempT_i = (int)(pairContTarget/5);// target sing point
																unsigned int tempT_j = fmod(pairContTarget,5);// target slot no
																if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==0){
																	tempPt1 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]];
																	tempPt2 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1];
																}
																else{
																	if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==1){
																		tempPt1 = line_discretization[tempS_i][tempS_j].back();
																		if(finalPaths[pairContSource][pairContTarget].size()!=0)
																			tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][0]];
																		else
																			tempPt2 = line_discretization[tempT_i][tempT_j].back();
																	}
																	else{
																		if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==2){
																			tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]]];
																			tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1]];
																		}
																		else{
																			if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==3){
																				if(finalPaths[pairContSource][pairContTarget].size()!=0)
																					tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																				else
																					tempPt1 = line_discretization[tempS_i][tempS_j].back();
																				tempPt2 = line_discretization[tempT_i][tempT_j].back();// !!!it is the inverse
																			}
																			else{
																				if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==4){
																					tempPt1 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]];
																					tempPt2 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1];// !!!it is the inverse
																				}
																				else{
																					if(finalPaths[pairContSource][pairContTarget].size()!=0)
																						tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																					else
																						tempPt1 = line_discretization[tempS_i][tempS_j].back();
																				 tempPt2 = finalEndPoint[pairContSource];
																				}
																			}
																		}
																	}
																}

																math::Segment segment1(tempPt1, tempPt2);
																if(from_seg.SecondMetIntersect2D(segment1, intersectionPoint, intersectionParam, temp_epsilon))
																	IllegalCross.push_back(make_pair(tempCS, contMatrix));
													 			//perhaps here also check if their direction corresponds...however they folow the cross component

															}
														}
													}
												}
											}
											isTraversedFaceCompVector[currentFace_id].first.push_back(contMatrix);

											isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first.push_back(4);
											isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first.push_back(i-1);
										}
										else{
											for(unsigned int j=0; j<isTraversedFaceCompVector[currentFace_id].second.size(); j++){
												unsigned int tempCS = isTraversedFaceCompVector[currentFace_id].second[j];
												if(tempCS!=contMatrix){
												if((contMatrixToSingularity[contMatrix].first!=contMatrixToSingularity[tempCS].first)&&
													((contMatrixToSingularity[contMatrix].second!=contMatrixToSingularity[tempCS].second)||
													(fmod(contMatrix, totalNumberOfVariables)==totalNumberOfVariables-1))){
														if((( (int)(contMatrix/totalNumberOfVariables)!=((int)(tempCS/totalNumberOfVariables)))&&
														( ( fmod(contMatrix, totalNumberOfVariables))!=(fmod(tempCS, totalNumberOfVariables)))&&
														( (int)(contMatrix/totalNumberOfVariables)!=(fmod(tempCS, totalNumberOfVariables)))&&
														( (fmod(contMatrix, totalNumberOfVariables))!=((int)(tempCS/totalNumberOfVariables)) ) )||
														(( ( fmod(contMatrix, totalNumberOfVariables))==(fmod(tempCS, totalNumberOfVariables)))&&(( fmod(contMatrix, totalNumberOfVariables))==totalNumberOfVariables-1))){
															bool toAdd = true;
															for(unsigned int k=0; k<IllegalCross.size(); k++){
																if(((contMatrix==IllegalCross[k].second)&&(tempCS==IllegalCross[k].first))||
																((tempCS==IllegalCross[k].second)&&(contMatrix==IllegalCross[k].first))){
																	toAdd = false;
																	break;
																}
															}
															if(toAdd){
																math::Point tempPt1, tempPt2;
																unsigned int pairContSource = (int)(tempCS/totalNumberOfVariables);
																unsigned int pairContTarget = fmod(tempCS,totalNumberOfVariables);
																unsigned int tempS_i = (int)(pairContSource/5);// source sing point
																unsigned int tempS_j = fmod(pairContSource,5);// source slot no
																unsigned int tempT_i = (int)(pairContTarget/5);// target sing point
																unsigned int tempT_j = fmod(pairContTarget,5);// target slot no
																if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==0){
																	tempPt1 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]];
																	tempPt2 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1];
																}
																else{
																	if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==1){
																		tempPt1 = line_discretization[tempS_i][tempS_j].back();
																		if(finalPaths[pairContSource][pairContTarget].size()!=0)
																			tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][0]];
																		else
																			tempPt2 = line_discretization[tempT_i][tempT_j].back();
																	}
																	else{
																		if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==2){
																			tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]]];
																			tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1]];
																		}
																		else{
																			if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==3){
																				if(finalPaths[pairContSource][pairContTarget].size()!=0)
																					tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																				else
																					tempPt1 = line_discretization[tempS_i][tempS_j].back();
																				tempPt2 = line_discretization[tempT_i][tempT_j].back();// !!!it is the inverse
																			}
																			else{
																				if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==4){
																					tempPt1 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]];
																					tempPt2 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1];// !!!it is the inverse
																				}
																				else{
																					if(finalPaths[pairContSource][pairContTarget].size()!=0)
																					 	tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																					else
																						tempPt1 = line_discretization[tempS_i][tempS_j].back();
																				tempPt2 = finalEndPoint[pairContSource];
																				}
																			}
																		}
																	}
																}

																math::Segment segment1(tempPt1, tempPt2);
																if(from_seg.SecondMetIntersect2D(segment1, intersectionPoint, intersectionParam, temp_epsilon))
																	IllegalCross.push_back(make_pair(tempCS, contMatrix));
																//perhaps here also check if their direction corresponds...however they folow the cross component

															}
														}
								 					}
												}
											}
											isTraversedFaceCompVector[currentFace_id].second.push_back(contMatrix);
											isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second.push_back(4);
											isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second.push_back(i-1);
										}
									}
									visitedFaces[edgeAdjFaces[k].id()] = true;
									traversedTriangles.push_back(edgeAdjFaces[k].id());
								}
							}
						}
					}
					if((!foundNextFace)&&(withComments))
						cout<<"!foundNextFace3"<<endl;

				}
					//currentPoint = discrPtsTarget[i];
			}
		}
	}
	else{
		if(withComments)
			cout<<"isTargetBdry"<<endl;
		if(lastVisTri != original_faces_number){
			currentFace_id = nextFace;
			nextFace = lastVisTri;

			nextPoint = finalEndPoint[contSource];
			math::Segment from_seg(currentPoint, nextPoint);

			math::Vector3d tempDir(currentPoint, nextPoint);
			if(withComments){
				cout<<"if(lastVisTri != original_faces_number){ "<<endl;
				cout<<"currentFace_id "<<currentFace_id<<endl;
				cout<<"nextFace "<<nextFace<<endl;
				cout<<"prevFace "<<prevFace<<endl;
				cout<<"finalEndPoint[contSource] "<<finalEndPoint[contSource]<<endl;
			}
			bool foundNextFace = false;

			if(withComments){
				cout<<"currentFace_id "<<currentFace_id<<endl;
				cout<<"nextFace "<<nextFace<<endl;
				cout<<"prevFace "<<prevFace<<endl;
			}
			vector<gmds::Edge> adj_edges;
			vector<gmds::Node> temp_nodes = (m_mesh->get<gmds::Face>(currentFace_id)).get<gmds::Node>();
			for(unsigned int j1=0; j1<3; j1++){
				vector<gmds::Edge> temp_edges = temp_nodes[j1].get<gmds::Edge>();
				adj_edges.insert(adj_edges.end(), temp_edges.begin(), temp_edges.end());
			}
			for(unsigned int j1=0; j1<adj_edges.size(); j1++){
				vector<gmds::Node> currentNodes =  adj_edges[j1].get<gmds::Node>();
				math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
				if(withComments){
					cout<<"edge number j = "<<j1<<" between nodes "<<currentNodes[0].id()<<" & "<<currentNodes[1]<<endl;
					from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon);

					cout<<"from_seg.SecondMetIntersect2D "<<from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon)
					<<endl;
					cout<<"temp_epsilon "<<temp_epsilon<<endl;
					cout<<"intersectionParam "<<intersectionParam<<endl;
				}
				if(from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon)){

					//WARNING if we want with nodes, here should be inserted

					vector<gmds::Face> edgeAdjFaces = adj_edges[j1].get<gmds::Face>();
					for(unsigned int k=0; k<edgeAdjFaces.size(); k++){
						if(!visitedFaces[edgeAdjFaces[k].id()]){
							currentFace_id = edgeAdjFaces[k].id();
							if(currentFace_id==nextFace)
								foundNextFace = true;

							math::Cross2D crossV = triangle_centers_cross[currentFace_id];
							unsigned int compVId = crossV.closestComponentVectorAsIndex(tempDir)%2;

							if(compVId==0){
								for(unsigned int j=0; j<isTraversedFaceCompVector[currentFace_id].first.size(); j++){
									unsigned int tempCS = isTraversedFaceCompVector[currentFace_id].first[j];
									if(tempCS!=contMatrix){
									if((contMatrixToSingularity[contMatrix].first!=contMatrixToSingularity[tempCS].first)&&
										((contMatrixToSingularity[contMatrix].second!=contMatrixToSingularity[tempCS].second)||
										(fmod(contMatrix, totalNumberOfVariables)==totalNumberOfVariables-1))){
											if((( (int)(contMatrix/totalNumberOfVariables)!=((int)(tempCS/totalNumberOfVariables)))&&
											( ( fmod(contMatrix, totalNumberOfVariables))!=(fmod(tempCS, totalNumberOfVariables)))&&
											( (int)(contMatrix/totalNumberOfVariables)!=(fmod(tempCS, totalNumberOfVariables)))&&
											( (fmod(contMatrix, totalNumberOfVariables))!=((int)(tempCS/totalNumberOfVariables)) ) )||
											(( ( fmod(contMatrix, totalNumberOfVariables))==(fmod(tempCS, totalNumberOfVariables)))&&(( fmod(contMatrix, totalNumberOfVariables))==totalNumberOfVariables-1))){
												bool toAdd = true;
												for(unsigned int k=0; k<IllegalCross.size(); k++){
													if(((contMatrix==IllegalCross[k].second)&&(tempCS==IllegalCross[k].first))||
													((tempCS==IllegalCross[k].second)&&(contMatrix==IllegalCross[k].first))){
														toAdd = false;
														break;
													}
												}
												if(toAdd){
													math::Point tempPt1, tempPt2;
													unsigned int pairContSource = (int)(tempCS/totalNumberOfVariables);
													unsigned int pairContTarget = fmod(tempCS,totalNumberOfVariables);
													unsigned int tempS_i = (int)(pairContSource/5);// source sing point
													unsigned int tempS_j = fmod(pairContSource,5);// source slot no
													unsigned int tempT_i = (int)(pairContTarget/5);// target sing point
													unsigned int tempT_j = fmod(pairContTarget,5);// target slot no
													if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==0){
														tempPt1 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]];
														tempPt2 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1];
													}
													else{
														if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==1){
															tempPt1 = line_discretization[tempS_i][tempS_j].back();
															if(finalPaths[pairContSource][pairContTarget].size()!=0)
																tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][0]];
															else
																tempPt2 = line_discretization[tempT_i][tempT_j].back();
														}
														else{
															if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==2){
																tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]]];
																tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1]];
															}
															else{
																if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==3){
																	if(finalPaths[pairContSource][pairContTarget].size()!=0)
																		tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																	else
																		tempPt1 = line_discretization[tempS_i][tempS_j].back();
																	tempPt2 = line_discretization[tempT_i][tempT_j].back();// !!!it is the inverse
																}
																else{
																	if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==4){
																		tempPt1 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]];
																		tempPt2 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1];// !!!it is the inverse
																	}
																	else{
																		if(finalPaths[pairContSource][pairContTarget].size()!=0)
																			tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																		else
																			tempPt1 = line_discretization[tempS_i][tempS_j].back();
																		tempPt2 = finalEndPoint[pairContSource];
																	}
																}
															}
														}
													}

													math::Segment segment1(tempPt1, tempPt2);
													if(from_seg.SecondMetIntersect2D(segment1, intersectionPoint, intersectionParam, temp_epsilon))
														IllegalCross.push_back(make_pair(tempCS, contMatrix));
													//perhaps here also check if their direction corresponds...however they folow the cross component

												}
											}
							 			}
									}
								}
								isTraversedFaceCompVector[currentFace_id].first.push_back(contMatrix);
								isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first.push_back(5);
								isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first.push_back(0);
							}
							else{ //compVId==1
								for(unsigned int j=0; j<isTraversedFaceCompVector[currentFace_id].second.size(); j++){
									unsigned int tempCS = isTraversedFaceCompVector[currentFace_id].second[j];
									if(tempCS!=contMatrix){
										if((contMatrixToSingularity[contMatrix].first!=contMatrixToSingularity[tempCS].first)&&
										((contMatrixToSingularity[contMatrix].second!=contMatrixToSingularity[tempCS].second)||
										(fmod(contMatrix, totalNumberOfVariables)==totalNumberOfVariables-1))){
											if((( (int)(contMatrix/totalNumberOfVariables)!=((int)(tempCS/totalNumberOfVariables)))&&
											( ( fmod(contMatrix, totalNumberOfVariables))!=(fmod(tempCS, totalNumberOfVariables)))&&
											( (int)(contMatrix/totalNumberOfVariables)!=(fmod(tempCS, totalNumberOfVariables)))&&
											( (fmod(contMatrix, totalNumberOfVariables))!=((int)(tempCS/totalNumberOfVariables)) ) )||
											(( ( fmod(contMatrix, totalNumberOfVariables))==(fmod(tempCS, totalNumberOfVariables)))&&(( fmod(contMatrix, totalNumberOfVariables))==totalNumberOfVariables-1))){
												bool toAdd = true;
												for(unsigned int k=0; k<IllegalCross.size(); k++){
													if(((contMatrix==IllegalCross[k].second)&&(tempCS==IllegalCross[k].first))||
													((tempCS==IllegalCross[k].second)&&(contMatrix==IllegalCross[k].first))){
														toAdd = false;
														break;
													}
												}
												if(toAdd){
													math::Point tempPt1, tempPt2;
													unsigned int pairContSource = (int)(tempCS/totalNumberOfVariables);
													unsigned int pairContTarget = fmod(tempCS,totalNumberOfVariables);
													unsigned int tempS_i = (int)(pairContSource/5);// source sing point
													unsigned int tempS_j = fmod(pairContSource,5);// source slot no
													unsigned int tempT_i = (int)(pairContTarget/5);// target sing point
													unsigned int tempT_j = fmod(pairContTarget,5);// target slot no
													if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==0){
														tempPt1 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]];
														tempPt2 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1];
													}
													else{
														if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==1){
															tempPt1 = line_discretization[tempS_i][tempS_j].back();
															if(finalPaths[pairContSource][pairContTarget].size()!=0)
																tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][0]];
															else
																tempPt2 = line_discretization[tempT_i][tempT_j].back();
														}
														else{
															if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==2){
																tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]]];
																tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1]];
															}
															else{
																if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==3){
																	if(finalPaths[pairContSource][pairContTarget].size()!=0)
																		tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																	else
																		tempPt1 = line_discretization[tempS_i][tempS_j].back();
																	tempPt2 = line_discretization[tempT_i][tempT_j].back();// !!!it is the inverse
																}
																else{
																	if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==4){
																		tempPt1 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]];
																		tempPt2 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1];// !!!it is the inverse
																	}
																	else{
																		if(finalPaths[pairContSource][pairContTarget].size()!=0)
																			tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																		else
																			tempPt1 = line_discretization[tempS_i][tempS_j].back();
																		tempPt2 = finalEndPoint[pairContSource];
																	}
																}
															}
														}
													}

													math::Segment segment1(tempPt1, tempPt2);
													if(from_seg.SecondMetIntersect2D(segment1, intersectionPoint, intersectionParam, temp_epsilon))
														IllegalCross.push_back(make_pair(tempCS, contMatrix));
													//perhaps here also check if their direction corresponds...however they folow the cross component

												}
											}
										}
									}
								}
								isTraversedFaceCompVector[currentFace_id].second.push_back(contMatrix);
								isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second.push_back(5);
								isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second.push_back(0);
							}

							visitedFaces[edgeAdjFaces[k].id()] = true;
							traversedTriangles.push_back(edgeAdjFaces[k].id());
						}
					}
				}
			}
			if((!foundNextFace)&&(withComments))
				cout<<"!foundNextFace4"<<endl;
		}
	}
	if(withGlobalComments)
		cout<<"getTraversedTrisFaceNeighbyVertsExhaustive end"<<endl;
}

//-------------------------------------------------------------------------------------------------

void SingularityGraphBuilder2D::getShortestPathBtwFacesOptimized(gmds::TCellID&                                                              source,
                                                                 vector<unsigned int>&                                                       targets,
                                                                 vector<gmds::math::Point>&                                                  targetPoints,
                                                                 vector<pair<double, unsigned int>>&                                         min_distance,
                                                                 vector<int>&                                                                previous,
                                                                 int&                                                                        found,
                                                                 vector<vector<double>>&                                                     face2FaceMatch,
                                                                 vector<vector<unsigned int>>&                                               face2FaceDeviation,
                                                                 pair<gmds::math::Vector3d, gmds::math::Vector3d>&                           prevDirCrossGlobal,
                                                                 vector<pair<gmds::math::Vector3d, gmds::math::Vector3d>>&                   targetDirCross,
                                                                 vector<bool>&                                                               is_bdry_face,
                                                                 gmds::math::Point&                                                          startPoint,
                                                                 vector<gmds::math::Point>&                                                  finalEndPoint,
                                                                 double&                                                                     maxDist,
                                                                 int&                                                                        final_to_cell_dim,
                                                                 gmds::TCellID&                                                              final_to_cell_id,
                                                                 gmds::TCellID&                                                              finalLastVisTri,
                                                                 unsigned int&                                                               contSource,
                                                                 unsigned int&                                                               totalNumberOfVariables,
                                                                 unsigned int&                                                               totalNumberOfSlots,
                                                                 vector<unsigned int>&                                                       faceNo2Cont,
                                                                 vector<vector<vector<gmds::math::Point>>>&                                  line_discretization,
                                                                 vector<vector<vector<gmds::TCellID>>>&                                      line_triangles,
                                                                 vector<vector<pair<vector<gmds::math::Point>, vector<gmds::math::Point>>>>& pointPaths,
                                                                 vector<vector<vector<gmds::TCellID>>>&                                      finalPaths,
                                                                 vector<vector<vector<gmds::math::Point>>>&                                  finalCenterPoints,
                                                                 vector<vector<vector<gmds::TCellID>>>&                                      traversedSPTris,
                                                                 vector<vector<double>>&                                                     distances,
                                                                 vector<pair<unsigned int, unsigned int>>&                                   IllegalCross,
                                                                 vector<pair<vector<unsigned int>, vector<unsigned int >>>&                  isTraversedNodeCompVector,
                                                                 vector<pair<vector<unsigned int>, vector<unsigned int >>>&                  isTraversedFaceCompVector,
                                                                 vector<pair<vector<unsigned int>, vector<unsigned int >>>&                  isTraversedFaceCompVector_SegmentPathCode,
                                                                 vector<pair<vector<unsigned int>, vector<unsigned int >>>&                  isTraversedFaceCompVector_SegmentPathCont,
                                                                 vector<pair<unsigned int, unsigned int>>&                                   contSourceToSingularity,
                                                                 vector<bool>&                                                               singOrGeomFaces)
{


	/*Description: get the shortest path between a singular triangle \param[in] source
	* and all other singular triangles \param[in] targets by 'walking from face center to face center'*/

	/*
	std::string file_name_bdry_faces = "forbiddenFaces.vtk";
	writeTestMeshTrianglesIds(forbiddenFacesNow,file_name_bdry_faces);*/

	unsigned int totNoTargets = 0;
	vector<vector<unsigned int>> neighbouringTarget(original_faces_number, vector<unsigned int> (0));
	vector<unsigned int> targetIndex2Cont(targets.size(), 0);
	for(unsigned int i=0; i<targets.size(); i++){
		if(targets[i]<original_faces_number){
			targetIndex2Cont[i] = totNoTargets;
			totNoTargets++;
			neighbouringTarget[targets[i]].push_back(i);

			for(unsigned int j=0; j<face2Face_neighbours_by_verts[targets[i]].size(); j++){
				neighbouringTarget[face2Face_neighbours_by_verts[targets[i]][j].id()].push_back(i);

			}
		}
	}

	bool visualizeTraversedTriangles = false;
	bool visualizeAll = false;

	if(withGlobalComments){
		cout<<"getShortestPathBtwFacesOptimized"<<endl;
		cout<<"contSource "<<contSource<<endl;
		cout<<"startPoint "<<startPoint<<endl;
		cout<<"finalEndPoint[contSource] "<<finalEndPoint[contSource]<<endl;
		cout<<"source tri  "<<source<<endl;
	}
	gmds::math::Point endPoint;

	gmds::TCellID lastVisTri = finalLastVisTri;
	int to_cell_dim = final_to_cell_dim;
	gmds::TCellID to_cell_id = final_to_cell_id;

	vector<unsigned int> targetsTemp(targets);
	std::vector<unsigned int>::iterator it;
	gmds::Face v, u;
	gmds::TCellID u_id, v_id;
	unsigned int tempFacesNo = original_faces_number + 1 + totNoTargets;
	vector<pair<gmds::math::Vector3d, gmds::math::Vector3d>> tempPrevDirCross(tempFacesNo, prevDirCrossGlobal); // pair of previous direction and previous followed cross component

	vector<pair<gmds::math::Vector3d, gmds::math::Vector3d>> tempTargetDirCross(targetDirCross);

	bool withComments = false;
	bool withDetailedComments = false;

/*	if((contSource==-1)){
		withComments = true;
		withDetailedComments = true;
	}
*/
	 //consider first face as imaginary (as well as the target faces)
	min_distance.clear();
	min_distance.resize(tempFacesNo, make_pair(maxDist, 1));
	min_distance[original_faces_number].first = 0.0;
	previous.clear();
	previous.resize(tempFacesNo, tempFacesNo+1);
	vector<unsigned int> paths2;
	unsigned int i_orig = (int) (contSource/5);
	unsigned int j_orig = fmod(contSource,5);

	double turnPenalizationTerm = 0.1 * maxDist;

	bool isTargetBdry;

	bool reachedBdry = false;

	std::set<std::pair<double,unsigned int> > vertex_queue;

	math::Vector3d closestToPrevDir, closestToPrevCross;

	for(unsigned int i=0; i<tempFacesNo; i++){
		vertex_queue.insert(std::make_pair(min_distance[i].first,i));
	}

	vector<bool> visitedFaces(tempFacesNo, false);

	tempFacesNo++;
	//////////////////////////////////////////////////////////////////////////////
	u = m_mesh->get<gmds::Face>(source);
	u_id = u.id();
	if(withComments){
		cout<<"treat first point (possibleTargetTriangles) "<<u_id<<endl;
		cout<<"original_faces_number "<<original_faces_number<<endl;
	}
	double total_dist_u = vertex_queue.begin()->first;
	unsigned int visitedFacesNo = min_distance[vertex_queue.begin()->second].second;

	vertex_queue.erase(vertex_queue.begin());
	visitedFaces[original_faces_number] = true;

	math::Vector3d prevDir_u =  prevDirCrossGlobal.first;

	bool validNeighbour = true;
	bool toCheck = false;

	if(validNeighbour){

		math::Vector3d tri2tri(startPoint, triangle_centers[u_id]);
		tri2tri.normalize();

		math::Cross2D crossV = triangle_centers_cross[u_id];
		closestToPrevDir =  crossV.closestComponentVector(prevDir_u);
		closestToPrevCross =  crossV.closestComponentVector(prevDir_u);

		double distance_through_u;

		if(prevDir_u.angle(tri2tri)>gmds::math::Constants::PIDIV4)
			distance_through_u =  total_dist_u + prevDir_u.angle(tri2tri)
			                                   + prevDir_u.angle(closestToPrevCross)
			                                   + turnPenalizationTerm;
			                                   //+ (floor)((prevDir_u.angle(tri2tri))/gmds::math::Constants::PIDIV2)*turnPenalizationTerm;
			                                   //+ (floor)((prevDir_u.angle(closestToPrevCross))/gmds::math::Constants::PIDIV2)*turnPenalizationTerm;

		else
			distance_through_u =  total_dist_u + prevDir_u.angle(tri2tri)
			                                   + prevDir_u.angle(closestToPrevCross);
			                                   //+ (floor)((prevDir_u.angle(closestToPrevCross))/gmds::math::Constants::PIDIV2)*9.0;

		if(withComments){
			cout<<"u_id "<<u_id<<endl;
			cout<<"tri2tri "<<startPoint<<" -> "<<triangle_centers[u_id]<<endl;
			if(withDetailedComments){
				cout<<"prevDir_u "<<startPoint<<" - > "<<prevDir_u[0] + startPoint.X()<<" , "<<prevDir_u[1] + startPoint.Y()<<endl;
				cout<<"closestToPrevCross "<<triangle_centers[u_id]<<" - > "<<closestToPrevCross[0] + triangle_centers[u_id].X()<<" , "<<closestToPrevCross[1] + triangle_centers[u_id].Y()<<endl;
				cout<<"prevDir_u.angle(tri2tri) rad "<<prevDir_u.angle(tri2tri)<<endl;

				cout<<"prevDir_u.angle(tri2tri) deg "<<prevDir_u.angle(tri2tri)* 180 / math::Constants::PI<<endl;

				cout<<"closestToPrevCross.angle(tri2tri) "<<closestToPrevCross.angle(tri2tri)* 180 / math::Constants::PI<<endl;
				cout<<"prevDir_u.angle(closestToPrevCross) rad "<<prevDir_u.angle(closestToPrevCross)<<endl;
				cout<<"prevDir_u.angle(closestToPrevCross) deg "<<prevDir_u.angle(closestToPrevCross)* 180 / math::Constants::PI<<endl;
				cout<<"(floor)((prevDir_u.angle(tri2tri))/gmds::math::Constants::PIDIV2) "<<(floor)((prevDir_u.angle(tri2tri))/gmds::math::Constants::PIDIV2)<<endl;

				cout<<"(floor)((prevDir_u.angle(closestToPrevCross))/gmds::math::Constants::PIDIV2) "<<(floor)((prevDir_u.angle(closestToPrevCross))/gmds::math::Constants::PIDIV2)<<endl;
				cout<<"distance_through_u "<<distance_through_u<<endl;
				cout<<"distance_through_u/(visitedFacesNo + 1)  "<<distance_through_u/(visitedFacesNo + 1) <<endl;
			}
		}

		//WARNING division par visitedFacesNo not particularly ok; especially for very irregular meshes
		if (distance_through_u/(visitedFacesNo + 1) < min_distance[u_id].first/min_distance[u_id].second) {

			if(withComments){
				cout<<"if(distance_through_u/VFN<min_distance[u_id]/VFN"<<endl;
			}
			vertex_queue.erase(std::make_pair(min_distance[u_id].first,u_id));

			min_distance[u_id].first = distance_through_u;
			min_distance[u_id].second = visitedFacesNo + 1;
			previous[u_id] = original_faces_number;

			vertex_queue.insert(std::make_pair(min_distance[u_id].first,u_id));

			tempPrevDirCross[u_id].first = tri2tri;
			tempPrevDirCross[u_id].second = closestToPrevCross;
		}
	}

	for(unsigned int i=0; i<face2Face_neighbours_by_verts[u_id].size(); i++){
		v = face2Face_neighbours_by_verts[u_id][i];
		v_id = v.id();
		bool validNeighbour = true;
		if(m_mesh->isMarked(v, m_mark_faces_with_sing_point))
			validNeighbour = false;

		if(validNeighbour){

			math::Vector3d tri2tri(startPoint, triangle_centers[v_id]);
			tri2tri.normalize();
			math::Cross2D crossV = triangle_centers_cross[v_id];
			closestToPrevDir =  crossV.closestComponentVector(prevDir_u);
			closestToPrevCross =  crossV.closestComponentVector(prevDir_u);
			double distance_through_u;

			if(prevDir_u.angle(tri2tri)>gmds::math::Constants::PIDIV4)
				distance_through_u =  total_dist_u + prevDir_u.angle(tri2tri)
				                                   + prevDir_u.angle(closestToPrevCross)
				                                   + turnPenalizationTerm;
				                                   // + (floor)((prevDir_u.angle(closestToPrevCross))/gmds::math::Constants::PIDIV2)*turnPenalizationTerm;

			else
				distance_through_u =  total_dist_u + prevDir_u.angle(tri2tri)
				                                   + prevDir_u.angle(closestToPrevCross);
				                                   // + (floor)((prevDir_u.angle(closestToPrevCross))/gmds::math::Constants::PIDIV2)*9.0;

			if(withComments){
				cout<<"v_id "<<v_id<<endl;
				cout<<"tri2tri "<<startPoint<<" -> "<<triangle_centers[v_id]<<endl;
				cout<<"prevDir_u "<<startPoint<<" - > "<<prevDir_u[0] + startPoint.X()<<" , "<<prevDir_u[1] + startPoint.Y()<<endl;

				if(withDetailedComments){
					cout<<"closestToPrevDir "<<triangle_centers[v_id].X()<<" -> "<<closestToPrevDir[0] + triangle_centers[v_id].X()<<" , "<<closestToPrevDir[1] + triangle_centers[v_id].Y()<<endl;
					cout<<"prevDir_u.angle(tri2tri) rad "<<prevDir_u.angle(tri2tri)<<endl;
					cout<<"prevDir_u.angle(tri2tri) deg "<<prevDir_u.angle(tri2tri)* 180 / math::Constants::PI<<endl;
					cout<<"not in eq closestToPrevCross.angle(tri2tri) "<<closestToPrevCross.angle(tri2tri)* 180 / math::Constants::PI<<endl;
					cout<<"not in eq tri2tri.angle(closestToPrevCross) "<<tri2tri.angle(closestToPrevCross)* 180 / math::Constants::PI<<endl;
					cout<<"(floor)((prevDir_u.angle(closestToPrevCross))/gmds::math::Constants::PIDIV2) "<<(floor)((prevDir_u.angle(closestToPrevCross))/gmds::math::Constants::PIDIV2)<<endl;
					cout<<"distance_through_u "<<distance_through_u<<endl;
				}
			}

			//WARNING division par visitedFacesNo not particularly ok; especially for very irregular meshes
			if (distance_through_u/(visitedFacesNo + 1) < min_distance[v_id].first/min_distance[v_id].second) {
				if(withComments){
					cout<<"if(distance_through_u/VFN<min_distance[v_id]/VFN"<<endl;
				}

				vertex_queue.erase(std::make_pair(min_distance[v_id].first,v_id));

				min_distance[v_id].first = distance_through_u;
				min_distance[v_id].second = visitedFacesNo + 1;
				previous[v_id] = original_faces_number;

				vertex_queue.insert(std::make_pair(min_distance[v_id].first,v_id));

				tempPrevDirCross[v_id].first = tri2tri;
				tempPrevDirCross[v_id].second = closestToPrevCross;

			}
		}
	}


	while(!vertex_queue.empty()){//while(vertex_queue.size()>(1+totNoTargets)){
		u_id = vertex_queue.begin()->second;
		double total_dist_u = vertex_queue.begin()->first;

		unsigned int visitedFacesNo = min_distance[vertex_queue.begin()->second].second;

		vertex_queue.erase(vertex_queue.begin());

		if(u_id<original_faces_number){

			u = m_mesh->get<gmds::Face>(u_id);
			if(withComments){
				cout<<"-------------------------------------------------------------------"<<endl;
				cout<<"u_id "<<u_id<<"|"<<endl;
				cout<<"!!!vertex_queue"<<endl;
				cout<<"total_dist_u "<<total_dist_u<<endl;
			}

			prevDir_u =  tempPrevDirCross[u_id].first;
			math::Vector3d prevCross_u =  tempPrevDirCross[u_id].second;

			if(!neighbouringTarget[u_id].empty()){

				for(unsigned int z=0; z<neighbouringTarget[u_id].size(); z++){
					unsigned int targetIndex = neighbouringTarget[u_id][z];
					unsigned int actualTarget = targets[targetIndex];
					if(withComments){
						cout<<"targetIndex "<<targetIndex<<endl;
						cout<<"actualTarget "<<actualTarget<<endl;
						cout<<"targetPoints[targetIndex] "<<targetPoints[targetIndex]<<endl;
					}
					math::Vector3d tri2tri(triangle_centers[u_id], targetPoints[targetIndex]);
					tri2tri.normalize();
					double distance_through_u;
					math::Vector3d tempTDC =  targetDirCross[targetIndex].first;

					if(withComments){
						cout<<"ShortestPaths between "<<contSource<<" and "<<targetIndex<<endl;
						cout<<"tempTDC = targetPoints[targetIndex]-> "<<tempTDC[0]+targetPoints[targetIndex].X()<<" , "<<tempTDC[1]+targetPoints[targetIndex].Y()<<endl;
						cout<<"prevCross_u "<<triangle_centers[u_id]<<" -> "<<prevCross_u[0]+triangle_centers[u_id].X()<<" , "<<prevCross_u[1]+triangle_centers[u_id].Y()<<endl;
						cout<<"tri2tri "<<triangle_centers[u_id]<<" -> "<<targetPoints[targetIndex]+triangle_centers[u_id]<<endl;
						cout<<"tempTDC.angle(prevCross_u) "<<tempTDC.angle(prevCross_u)<<endl;
						cout<<"tri2tri.angle(prevCross_u) "<<tri2tri.angle(prevCross_u)<<endl;
					}
					if((tempTDC.angle(prevCross_u)>gmds::math::Constants::PIDIV4)
					||(tri2tri.angle(prevCross_u)>gmds::math::Constants::PIDIV4)){
						distance_through_u =  total_dist_u + tri2tri.angle(prevCross_u)
						                                   + tempTDC.angle(prevCross_u)
 						                                   + 3*turnPenalizationTerm;
						if(visualizeAll){
							if(((double)distance_through_u/(visitedFacesNo+1)<distances[contSource][targetIndex])){
								distances[contSource][targetIndex] = (double)distance_through_u/(visitedFacesNo+1);
								if(withComments){
									cout<<"if distance_through_u/VFN < distances[contSource][targetIndex] "<<distance_through_u<<endl;
									cout<<"tri2tri "<<triangle_centers[u_id]<<" -> "<<targetPoints[targetIndex]<<endl;
									cout<<"targetIndex2Cont["<<targetIndex<<"] "<<targetIndex2Cont[targetIndex]<<endl;
									cout<<"retraceShortestPath"<<endl;
								}

								retraceShortestPath(u_id, previous, paths2, tempFacesNo);
								reverse(paths2.begin(),paths2.end());
								if(withComments){
									cout<<"paths2 "<<endl;
									for(unsigned int k77=0;k77<paths2.size(); k77++)
										cout<<paths2[k77]<<endl;
									cout<<"!---!"<<endl;
								}

								finalPaths[contSource][targetIndex].clear();
								if(paths2.size()>1)
									finalPaths[contSource][targetIndex].insert(finalPaths[contSource][targetIndex].end(),
									                                    paths2.begin(),paths2.end()-1);
								min_distance[original_faces_number + 1 + targetIndex2Cont[targetIndex]].first = distance_through_u;
								min_distance[original_faces_number + 1 + targetIndex2Cont[targetIndex]].second = visitedFacesNo + 1;
								tempPrevDirCross[original_faces_number + 1 + targetIndex2Cont[targetIndex]].first = tri2tri;
								tempPrevDirCross[original_faces_number + 1 + targetIndex2Cont[targetIndex]].second = closestToPrevCross;

							}
						}
					}
					else{
						distance_through_u =  total_dist_u + tri2tri.angle(prevCross_u)
						                                   + tempTDC.angle(prevCross_u);
						                                   //+10*(tempTDC[targetIndex].first).angle(prevDir);

						if(withComments){
							cout<<"distance_through_u "<<distance_through_u<<" , (double)distance_through_u/(visitedFacesNo+1) "<<(double)distance_through_u/(visitedFacesNo+1)<<endl;
							cout<<"distances[contSource][targetIndex] "<<distances[contSource][targetIndex]<<endl;
						}
						if(((double)distance_through_u/(visitedFacesNo+1)<distances[contSource][targetIndex])){
							distances[contSource][targetIndex] = (double)distance_through_u/(visitedFacesNo+1);
							if(withComments){
								cout<<"if distance_through_u/VFN < distances[contSource][targetIndex] "<<distance_through_u<<endl;
								cout<<"tri2tri "<<triangle_centers[u_id]<<" -> "<<targetPoints[targetIndex]<<endl;
								cout<<"targetIndex2Cont["<<targetIndex<<"] "<<targetIndex2Cont[targetIndex]<<endl;
							}
							retraceShortestPath(u_id, previous, paths2, tempFacesNo);
							reverse(paths2.begin(),paths2.end());
							if(withComments){
								cout<<"paths2 "<<endl;
								for(unsigned int k77=0;k77<paths2.size(); k77++)
									cout<<paths2[k77]<<endl;
								cout<<"!---!"<<endl;
							}

							finalPaths[contSource][targetIndex].clear();
							if(paths2.size()>1)
								finalPaths[contSource][targetIndex].insert(finalPaths[contSource][targetIndex].end(),
								                                    paths2.begin(),paths2.end()-1);
							min_distance[original_faces_number + 1 + targetIndex2Cont[targetIndex]].first = distance_through_u;
							min_distance[original_faces_number + 1 + targetIndex2Cont[targetIndex]].second = visitedFacesNo + 1;
							tempPrevDirCross[original_faces_number + 1 + targetIndex2Cont[targetIndex]].first = tri2tri;
							tempPrevDirCross[original_faces_number + 1 + targetIndex2Cont[targetIndex]].second = closestToPrevCross;

						}
					}
				}
			}

			if(is_bdry_face[u_id]){ //(!reachedBdry)
				double bdry_dist;
				to_cell_dim = 2;
				to_cell_id = u_id;
				//math::Point intersectionPnt;

				double intersectionParam;
				math::Ray from_ray(triangle_centers[u_id], prevCross_u);
				if(withComments){
					cout<<"is_bdry_face[u_id]"<<endl;
					if(withDetailedComments){
						cout<<"previous["<<u_id<<"]= "<<previous[u_id]<<endl;
						cout<<"prevDir previous[u_id] "<<triangle_centers[previous[u_id]]<<" -> "<<
						triangle_centers[previous[u_id]].X()+tempPrevDirCross[previous[u_id]].first[0]<<" , "<<
						triangle_centers[previous[u_id]].Y()+tempPrevDirCross[previous[u_id]].first[1]<<endl;
						cout<<"prevCross previous[u_id] "<<triangle_centers[previous[u_id]]<<" -> "<<
						triangle_centers[previous[u_id]].X()+tempPrevDirCross[previous[u_id]].second[0]<<" , "<<
						triangle_centers[previous[u_id]].Y()+tempPrevDirCross[previous[u_id]].second[1]<<endl;
						cout<<"prevDir_u "<<triangle_centers[u_id]<<" -> "<<triangle_centers[u_id].X()+prevDir_u[0]<<" , "<<triangle_centers[u_id].Y()+prevDir_u[1]<<endl;

						cout<<"from_ray (prevCross_u) "<<triangle_centers[u_id]<<" -> "<<triangle_centers[u_id].X()+prevCross_u[0]<<" , "<<triangle_centers[u_id].Y()+prevCross_u[1]<<endl;
						cout<<"previous[previous["<<u_id<<"]]= "<<previous[previous[u_id]]<<endl;
						cout<<"prevDir previous[previous[u_id]] "<<triangle_centers[previous[previous[u_id]]]<<" -> "<<
						triangle_centers[previous[previous[u_id]]].X()+tempPrevDirCross[previous[previous[u_id]]].first[0]<<" , "<<
						triangle_centers[previous[previous[u_id]]].Y()+tempPrevDirCross[previous[previous[u_id]]].first[1]<<endl;
						cout<<"prevCross previous[previous[u_id]] "<<triangle_centers[previous[previous[u_id]]]<<" -> "<<
						triangle_centers[previous[previous[u_id]]].X()+tempPrevDirCross[previous[previous[u_id]]].second[0]<<" , "<<
						triangle_centers[previous[previous[u_id]]].Y()+tempPrevDirCross[previous[previous[u_id]]].second[1]<<endl;
					}
				}

				vector<gmds::Edge> currentEdges = u.get<gmds::Edge>();
				bool hasBdryEdge = false;
				reachedBdry = false;
				for(unsigned int i=0; i<currentEdges.size(); i++){
					if (m_mesh->isMarked(currentEdges[i], m_mark_edges_on_curve)){

						hasBdryEdge = true;
						to_cell_dim = 1;
						to_cell_id = currentEdges[i].id();

						vector<gmds::Node> currentNodes =  currentEdges[i].get<Node>();
						math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
						if(withDetailedComments){
						 cout<<" is bdry edge betinitConfusingBallsween  "<<currentNodes[0].id()<<" and "<<currentNodes[1].id()<<endl;
							from_ray.SecondMetIntersect2D(oppSeg, endPoint, intersectionParam, temp_epsilon);
							cout<<" intersectionParam "<<intersectionParam<<endl;
						}

						if(from_ray.SecondMetIntersect2D(oppSeg, endPoint, intersectionParam, temp_epsilon)){
							double tempAngleVal = 1.0 - fabs(prevCross_u.dot(bdry_edge_normals[currentEdges[i].id()]));
							if(withDetailedComments){
								cout<<"tempAngleVal "<<tempAngleVal<<" between from_ray and the bdry edge above"<<endl;
							}
							if(tempAngleVal<gmds::math::Constants::PIDIV4){
								bdry_dist =  total_dist_u + tempAngleVal;

								if(withComments){
									cout<<" u_id "<<u_id<<" , endPoint "<<endPoint<<endl;
									cout<<"tempAngleVal "<<tempAngleVal<<endl;
								}
								reachedBdry = true;
								lastVisTri = u_id;

							}
						}
					}
			 	}
				if(!hasBdryEdge){
			 		// VERTEX-> angle with vert normal
					if(withComments){
						cout<<"if(!hasBdryEdge){"<<endl;
					}
					vector<gmds::Node> currentNodes = u.get<gmds::Node>();
					for(unsigned int i=0; i<currentNodes.size(); i++){
						if(m_mesh->isMarked(currentNodes[i], m_mark_nodes_on_point) ||
						m_mesh->isMarked(currentNodes[i], m_mark_nodes_on_curve)){
							if(withDetailedComments){
								cout<<"node "<<currentNodes[i].id()<<" is bdry node"<<endl;
							}
							to_cell_dim = 0;
							to_cell_id = currentNodes[i].id();
 							vector<gmds::Edge> currentEdges = currentNodes[i].get<gmds::Edge>();
							for(unsigned int j=0; j<currentEdges.size(); j++){
								if(m_mesh->isMarked(currentEdges[j], m_mark_edges_on_curve)){
									vector<gmds::Node> edge_nodes = currentEdges[j].get<gmds::Node>();
									math::Segment oppSeg(edge_nodes[0].getPoint(), edge_nodes[1].getPoint());
									if(withDetailedComments){
										from_ray.SecondMetIntersect2D(oppSeg, endPoint, intersectionParam, temp_epsilon);
										cout<<" check edge between  "<<edge_nodes[0].id()<<" and "<<edge_nodes[1].id()<<endl;
										cout<<" intersectionParam "<<intersectionParam<<endl;
									}
									if(from_ray.SecondMetIntersect2D(oppSeg, endPoint, intersectionParam, temp_epsilon)){
										double tempAngleVal = 1.0 - fabs(prevCross_u.dot(bdry_edge_normals[currentEdges[j].id()]));
										to_cell_dim = 1;
										to_cell_id = currentEdges[j].id();
										if(withDetailedComments){
											cout<<"tempAngleVal "<<tempAngleVal<<" between from_ray and the check edge above"<<endl;

										}
										if(tempAngleVal<gmds::math::Constants::PIDIV4){
											//if(((double)(total_dist_u + tempAngleVal)/visitedFacesNo<(double)min_distance[u_id].first/(min_distance[u_id].second))||(!visitedFaces[u_id])){

											bdry_dist =  total_dist_u + tempAngleVal;

											vector<gmds::Face> adj_faces1 = currentEdges[j].get<gmds::Face>();
											lastVisTri = adj_faces1[0].id();
											if(withComments){
												cout<<" u_id "<<u_id<<" , endPoint "<<endPoint<<" , tempAngleVal "<<tempAngleVal<<endl;

											}
											reachedBdry = true;
										}
									}
								}
							}
						}
					}
				}
				if(reachedBdry){
					if(withComments){
						cout<<"!!!reachedBdry"<<"; bdry_dist= "<<bdry_dist<<" , (double)bdry_dist/(visitedFacesNo) "<<(double)bdry_dist/(visitedFacesNo)<<endl;
						cout<<"before distances[contSource][totalNumberOfSlots] "<<distances[contSource][totalNumberOfSlots]<<endl;
					}
					if(distances[contSource][totalNumberOfSlots] > (double)bdry_dist/(visitedFacesNo)){
						retraceShortestPath(u_id, previous, paths2, tempFacesNo);
						if(withComments){
							cout<<"paths2 "<<endl;
							for(unsigned int k77=0;k77<paths2.size(); k77++)
								cout<<paths2[k77]<<endl;
						}
						bool toConsider = true;
						for(unsigned int tk=0; tk<paths2.size(); tk++){
							if(singOrGeomFaces[paths2[tk]]){
								toConsider = false;
									break;
							}
						}
						if(toConsider){
							//WARNING TODO 	CHECK IF this is a geom face and if so, check how close are endPoint and the actual geom point
							finalEndPoint[contSource] = endPoint;
							finalLastVisTri = lastVisTri;
							final_to_cell_dim = to_cell_dim;
							final_to_cell_id = to_cell_id;
							distances[contSource][totalNumberOfSlots] = (double)bdry_dist/(visitedFacesNo);
							if(withComments){
								cout<<"after distances[contSource][totalNumberOfSlots] "<<distances[contSource][totalNumberOfSlots]<<endl;
								cout<<"finalEndPoint[contSource] "<<finalEndPoint[contSource]<<endl;
							}
							pointPaths[contSource][totalNumberOfSlots].first.clear();
							pointPaths[contSource][totalNumberOfSlots].first.insert(pointPaths[contSource][totalNumberOfSlots].first.end(),line_discretization[i_orig][j_orig].begin(),line_discretization[i_orig][j_orig].end());

							reverse(paths2.begin(),paths2.end());
							finalPaths[contSource][totalNumberOfSlots].clear();
							if(paths2.size()>0)
          						finalPaths[contSource][totalNumberOfSlots].insert(finalPaths[contSource][totalNumberOfSlots].end(),
									                                      paths2.begin(),paths2.end());
							lastVisTri = original_faces_number;
						}
					}
				}
				else{
					if(withComments){
						cout<<"NOT reachedBdry"<<endl;
					}
				}
			}
			if(withComments){
				cout<<"!!!    min_distance["<<u_id<<"].first "<<min_distance[u_id].first<<endl;
			}
			visitedFaces[u_id] = true;
			toCheck = false;

			for(unsigned int i=0; i<face2Face_neighbours_by_verts[u_id].size(); i++){

				v = face2Face_neighbours_by_verts[u_id][i];
				v_id = v.id();

				bool validNeighbour = true;
				if(m_mesh->isMarked(v, m_mark_faces_with_sing_point))
					validNeighbour = false;

				if((!visitedFaces[v_id])&&(validNeighbour)){
					math::Vector3d tri2tri(triangle_centers[u_id], triangle_centers[v_id]);
					tri2tri.normalize();
					math::Cross2D crossV = triangle_centers_cross[v_id];
					closestToPrevDir =  crossV.closestComponentVector(prevDir_u);
					closestToPrevCross =  crossV.closestComponentVector(prevCross_u);

					double distance_through_u;

					if(prevCross_u.angle(tri2tri)>gmds::math::Constants::PIDIV4)//if(tri2tri.angle(prevCross_u)>=gmds::math::Constants::PIDIV4)
						distance_through_u =  total_dist_u + prevCross_u.angle(tri2tri)
						                                   + prevCross_u.angle(closestToPrevCross)
						                                   + turnPenalizationTerm;
						                                   //+ (floor)((prevDir_u.angle(closestToPrevCross))/gmds::math::Constants::PIDIV2)*turnPenalizationTerm;

					else{
						distance_through_u =  total_dist_u + prevCross_u.angle(tri2tri)
						                                   + prevCross_u.angle(closestToPrevCross);	//max 3.14*10
						                                   //+ (floor)((prevDir_u.angle(closestToPrevCross))/gmds::math::Constants::PIDIV2)*9.0;// max 2*90
						                                    //10 and 90 to overcome the problems introduced by the division by number of visitedFaces
					}
                         if(withComments){
						cout<<"tri2tri "<<triangle_centers[u_id]<<" -> "<<triangle_centers[v_id]<<endl;

						cout<<"prevDir_u "<<triangle_centers[u_id]<<" - > "<<prevDir_u[0] + triangle_centers[u_id].X()<<" , "<<prevDir_u[1] + triangle_centers[u_id].Y()<<endl;
						if(withDetailedComments){
							cout<<"prevCross_u ->"<<prevCross_u<<endl;
							cout<<"closestToPrevDir "<<triangle_centers[v_id].X()<<" -> "<<closestToPrevDir[0] + triangle_centers[v_id].X()<<" , "<<closestToPrevDir[1] + triangle_centers[v_id].Y()<<endl;

							cout<<"condition tri2tri.angle(prevCross_u) "<<tri2tri.angle(prevCross_u)* 180 / math::Constants::PI<<endl;

							cout<<"closestToPrevCross.angle(tri2tri) "<<closestToPrevCross.angle(tri2tri)* 180 / math::Constants::PI<<endl;
							cout<<"closestToPrevCross.angle(prevCross_u) "<<closestToPrevCross.angle(prevCross_u)* 180 / math::Constants::PI<<endl;

							cout<<"not in eq (floor)((prevCross_u.angle(closestToPrevCross))/gmds::math::Constants::PIDIV2) "<<(floor)((prevCross_u.angle(closestToPrevCross))/gmds::math::Constants::PIDIV2)<<endl;
							cout<<" in eq (floor)((tri2tri.angle(closestToPrevCross))/gmds::math::Constants::PIDIV2) "<<(floor)((tri2tri.angle(closestToPrevCross))/gmds::math::Constants::PIDIV2)<<endl;
							cout<<"distance_through_u "<<distance_through_u<<endl;
							if((floor)((closestToPrevDir.angle(closestToPrevCross))/gmds::math::Constants::PIDIV2)>0){
								cout<<"!!! CTPD!=CTPC; to look into"<<endl;
							}
						}
					}

               		//WARNING division par visitedFacesNo not particularly ok; especially for very irregular meshes
					//if (distance_through_u/(visitedFacesNo + 1) < min_distance[v_id].first/min_distance[v_id].second) {
					//if(min_distance[v_id].first==maxDist||((distance_through_u<4*turnPenalizationTerm)&&(distance_through_u/(visitedFacesNo + 1) < min_distance[v_id].first/min_distance[v_id].second))) {
					if(((min_distance[v_id].first==maxDist)&&(distance_through_u<maxDist))||((distance_through_u<4*turnPenalizationTerm)&&(distance_through_u/(visitedFacesNo + 1) < min_distance[v_id].first/min_distance[v_id].second))) {

						if(withComments){
							cout<<"if distance_through_u<4*turnPenalizationTerm "<<endl;
						}

						vertex_queue.erase(std::make_pair(min_distance[v_id].first,v_id));

						min_distance[v_id].first = distance_through_u;
						min_distance[v_id].second = visitedFacesNo + 1;

						vertex_queue.insert(std::make_pair(min_distance[v_id].first,v_id));

						previous[v_id] = u_id;

						tempPrevDirCross[v_id].first = tri2tri;
						tempPrevDirCross[v_id].second = closestToPrevCross;
					}
				}
			}
		}
	}

 	if(withComments){
 		cout<<" while((!vertex_queue.empty())&&(!reachedBdry)){"<<endl;
		cout<<"reachedBdry "<<reachedBdry<<endl;
		cout<<"vertex_queue.size() "<<vertex_queue.size()<<endl;
	}

	//visited all

	for(unsigned int i=0; i<targetsTemp.size();i++){
		u_id = targetsTemp[i];
		bool withCommentsLocal = false;
		if(u_id<original_faces_number){

			if(withCommentsLocal)
				cout<<"distances[contSource][i] "<<distances[contSource][i]<<endl;

			if((min_distance[original_faces_number + 1 + targetIndex2Cont[i]].first<=0.3*maxDist)||(visualizeAll)){
				//this means less than 3 turns
				if(withComments){
					cout<<"targetIndex2Cont["<<i<<"] "<<targetIndex2Cont[i]<<endl;
					cout<<"min_distance[original_faces_number + targetIndex2Cont[i]].first= "<<min_distance[original_faces_number + targetIndex2Cont[i]].first<<endl; //
				}
				/*if want to visualize all if(distances[contSource][i]<maxDist){	*/
				pointPaths[contSource][i].first.clear();
				pointPaths[contSource][i].first.insert(pointPaths[contSource][i].first.end(),line_discretization[i_orig][j_orig].begin(),line_discretization[i_orig][j_orig].end());

				unsigned int t1_orig = (int) (i/5);
				unsigned int t2_orig = fmod(i,5);

				reverse(line_discretization[t1_orig][t2_orig].begin(), line_discretization[t1_orig][t2_orig].end());
				pointPaths[contSource][i].second.insert(pointPaths[contSource][i].second.end(),line_discretization[t1_orig][t2_orig].begin(),line_discretization[t1_orig][t2_orig].end());

				if(withCommentsLocal){
					cout<<"pointPaths[contSource][i].second.size() "<<pointPaths[contSource][i].second.size()<<endl;
					cout<<"pointPaths[contSource][i].second.size() "<<pointPaths[contSource][i].second.size()<<endl;
					cout<<"finalPaths[contSource][i][tt]"<<endl;

					for(unsigned int tt=0; tt<finalPaths[contSource][i].size(); tt++){
						cout<<finalPaths[contSource][i][tt]<<endl;
					}
				}
				vector<math::Point> centerPoints;
				centerPoints.insert(centerPoints.end(),line_discretization[i_orig][j_orig].begin(),line_discretization[i_orig][j_orig].end());

				for(unsigned int tt=0; tt<finalPaths[contSource][i].size(); tt++){
					centerPoints.push_back(triangle_centers[finalPaths[contSource][i][tt]]);
				}

				bool passesThroughDifferentSing = false;
				for(unsigned int k77=0; k77<finalPaths[contSource][i].size(); k77++){
					if(singOrGeomFaces[finalPaths[contSource][i][k77]]){
						passesThroughDifferentSing = true;
					}
				}
				if(passesThroughDifferentSing)
					distances[contSource][i] = distances[contSource][i] + 0.5*turnPenalizationTerm;


				centerPoints.insert(centerPoints.end(),line_discretization[t1_orig][t2_orig].begin(),line_discretization[t1_orig][t2_orig].end());

				finalCenterPoints[contSource][5*t1_orig+t2_orig].insert(finalCenterPoints[contSource][5*t1_orig+t2_orig].end(), centerPoints.begin(), centerPoints.end());
				if(withComments)
					cout<<"distances[contSource][i] "<<distances[contSource][i]<<endl;
				if((visualizeShortestPaths)&&(centerPoints.size()>=2)){
					std::string file_name = "ShortestPaths_"+std::to_string(contSource)+"-"+std::to_string(i)+".vtk";
					writeTestPoints(centerPoints, file_name);
				}
				else{
					std::string file_name = "ShortestPaths!!!In_"+std::to_string(contSource)+"-"+std::to_string(i)+".vtk";
					writeTestPoints(centerPoints, file_name);
				}

				isTargetBdry = false;

				unsigned int contMatrix = contSource*totalNumberOfVariables + i;
				contSourceToSingularity[contMatrix].first = i_orig;
				contSourceToSingularity[contMatrix].second = t1_orig;

				reverse(line_triangles[t1_orig][t2_orig].begin(),line_triangles[t1_orig][t2_orig].end());

				gmds::TCellID tempLastVisTri = original_faces_number;
				/*if((contSource==1)&&(i==17)){
					cout<<"before"<<endl;
					cout<<"if((contSource==1)&&(i==17)){"<<endl;
					for(unsigned int k7=0; k7<traversedSPTris[contSource][i].size(); k7++)
						cout<<"traversedSPTris[contSource][i]["<<k7<<"]= "<<traversedSPTris[contSource][i][k7]<<endl;
						cout<<"-------------------------------------------------------------------"<<endl;
						for(unsigned int k7=0; k7<line_triangles[i_orig][j_orig].size(); k7++)
						    cout<<"line_triangles[i_orig][j_orig]"<<k7<<"]= "<<line_triangles[i_orig][j_orig][k7]<<endl;

						cout<<"-------------------------------------------------------------------"<<endl;
						for(unsigned int k7=0; k7<line_triangles[t1_orig][t2_orig].size(); k7++)
						    cout<<"line_triangles[t1_orig][t2_orig]"<<k7<<"]= "<<line_triangles[t1_orig][t2_orig][k7]<<endl;
				 }*/

				if(withDetailedComments){
					cout<<"finalPaths[contSource][i]"<<endl;
					for(unsigned int k77=0; k77<finalPaths[contSource][i].size(); k77++)
						cout<<finalPaths[contSource][i][k77]<<endl;

					cout<<"line_triangles[i_orig][j_orig]"<<endl;
					for(unsigned int k77=0; k77<line_triangles[i_orig][j_orig].size(); k77++)
						cout<<line_triangles[i_orig][j_orig][k77]<<endl;

					cout<<"line_discretization[i_orig][j_orig]"<<endl;
					for(unsigned int k77=0; k77<line_discretization[i_orig][j_orig].size(); k77++)
						cout<<line_discretization[i_orig][j_orig][k77]<<endl;

					cout<<"line_discretization[t1_orig][t2_orig]"<<endl;
					for(unsigned int k77=0; k77<line_discretization[t1_orig][t2_orig].size(); k77++)
						cout<<line_discretization[t1_orig][t2_orig][k77]<<endl;

					cout<<"contSource "<<contSource<<endl;
				}
				unsigned int contTarget = i;
				getTraversedTrisFaceNeighbyVertsExhaustive(line_discretization,
				                                           line_triangles,
				                                           finalPaths,//[contSource][i],
				                                           isTargetBdry,
				                                           traversedSPTris[contSource][i],
				                                           isTraversedNodeCompVector,
				                                           isTraversedFaceCompVector,
				                                           isTraversedFaceCompVector_SegmentPathCode,
				                                           isTraversedFaceCompVector_SegmentPathCont,
				                                           contMatrix,
				                                           totalNumberOfVariables,
				                                           contSource,
				                                           contTarget,
				                                           IllegalCross,
				                                           contSourceToSingularity,
				                                           tempLastVisTri,
				                                           finalEndPoint,
				                                           i_orig,
				                                           j_orig,
				                                           t1_orig,
				                                           t2_orig);
				if(withComments){
					cout<<"traversedSPTris[contSource][i]"<<endl;
					for(unsigned int k77=0; k77<traversedSPTris[contSource][i].size(); k77++)
						cout<<traversedSPTris[contSource][i][k77]<<endl;
				}
				reverse(line_triangles[t1_orig][t2_orig].begin(),line_triangles[t1_orig][t2_orig].end());
				reverse(line_discretization[t1_orig][t2_orig].begin(),line_discretization[t1_orig][t2_orig].end());

				if(visualizeTraversedTriangles){
					std::string file_name_trav_tri = "trav_tri_"+std::to_string(contSource)+"-"+std::to_string(i)+".vtk";
					writeTestMeshTrianglesIds(traversedSPTris[contSource][i], file_name_trav_tri);
				}
			}
		}
	}



	unsigned int firstPointsNo = pointPaths[contSource][totalNumberOfSlots].first.size();
	unsigned int finalPathsNo = finalPaths[contSource][totalNumberOfSlots].size();
	if(withComments){
		cout<<"firstPointsNo "<<firstPointsNo<<endl;
		cout<<"finalPathsNo "<<finalPathsNo<<endl;
	}
	if(finalPathsNo>=1){
		bool passesThroughDifferentSing = false;
		for(unsigned int j=0; j<finalPathsNo; j++){
			if(singOrGeomFaces[finalPaths[contSource][totalNumberOfSlots][j]]){
				passesThroughDifferentSing = true;
			}
		}
		if(passesThroughDifferentSing)
			distances[contSource][totalNumberOfSlots] = distances[contSource][totalNumberOfSlots] + 0.5*turnPenalizationTerm;

		vector<gmds::math::Point> centerPoints(firstPointsNo + finalPathsNo);
		//finalCenterPoints[contSource][totalNumberOfSlots].clear();
		//finalCenterPoints[contSource][totalNumberOfSlots].resize(firstPointsNo + finalPathsNo);
		for(unsigned int j=0; j<firstPointsNo;j++){
			centerPoints[j] = pointPaths[contSource][totalNumberOfSlots].first[j];
		}

		for(unsigned int j=0; j<finalPathsNo; j++){
			centerPoints[j+firstPointsNo] = triangle_centers[finalPaths[contSource][totalNumberOfSlots][j]];
		}

		centerPoints.push_back(finalEndPoint[contSource]);
		if(withGlobalComments)
			cout<<"finalEndPoint[contSource] "<<finalEndPoint[contSource]<<endl;
		if(centerPoints.size()>=2){
			std::string file_name_bdry = "ShortestPathsBdry_"+std::to_string(contSource)+".vtk";
			writeTestPoints(centerPoints, file_name_bdry);
			cout<<"has written "<<file_name_bdry<<endl;

		}
	}
	else{
		if(firstPointsNo>=2){
			//vector<gmds::math::Point> centerPoints(firstPointsNo + finalPathsNo);

			//finalCenterPoints[contSource][totalNumberOfSlots].clear();
			//finalCenterPoints[contSource][totalNumberOfSlots].insert(finalCenterPoints[contSource][totalNumberOfSlots].end(), pointPaths[contSource][totalNumberOfSlots].first.begin(), pointPaths[contSource][totalNumberOfSlots].first.end());
			//finalCenterPoints[contSource][totalNumberOfSlots].insert(finalCenterPoints[contSource][totalNumberOfSlots].end(), centerPoints.begin(), centerPoints.end());
			std::string file_name_bdry = "ShortestPathsBdry_"+std::to_string(contSource)+".vtk";
			writeTestPoints(pointPaths[contSource][totalNumberOfSlots].first, file_name_bdry);
			cout<<"has written "<<file_name_bdry<<endl;
		}

	}

	if(pointPaths[contSource][totalNumberOfSlots].first.size()>0){
		isTargetBdry = true;
		unsigned int contMatrix = contSource*totalNumberOfVariables + totalNumberOfSlots;
		contSourceToSingularity[contMatrix].first = i_orig;
		contSourceToSingularity[contMatrix].second = i_orig;

		if(withDetailedComments){
			cout<<"finalPaths[contSource][totalNumberOfSlots]"<<endl;
			for(unsigned int k77=0; k77<finalPathsNo; k77++)
				cout<<finalPaths[contSource][totalNumberOfSlots][k77]<<endl;

			cout<<"line_discretization[i_orig][j_orig]"<<endl;
			for(unsigned int k77=0; k77<line_discretization[i_orig][j_orig].size(); k77++)
			cout<<line_discretization[i_orig][j_orig][k77]<<endl;

			cout<<"contSource "<<contSource<<endl;
		}
		getTraversedTrisFaceNeighbyVertsExhaustive(line_discretization,
		                                           line_triangles,
		                                           finalPaths,
		                                           isTargetBdry,
		                                           traversedSPTris[contSource][totalNumberOfSlots],
		                                           isTraversedNodeCompVector,
		                                           isTraversedFaceCompVector,
		                                           isTraversedFaceCompVector_SegmentPathCode,
		                                           isTraversedFaceCompVector_SegmentPathCont,
		                                           contMatrix,
		                                           totalNumberOfVariables,
		                                           contSource,
		                                           totalNumberOfSlots,
		                                           IllegalCross,
		                                           contSourceToSingularity,
		                                           finalLastVisTri,
		                                           finalEndPoint,
		                                           i_orig,
		                                           j_orig,
		                                           i_orig,
		                                           j_orig);
		if(withDetailedComments){
			for(unsigned int k77=0; k77<traversedSPTris[contSource][totalNumberOfSlots].size(); k77++)
				cout<<traversedSPTris[contSource][totalNumberOfSlots][k77]<<endl;
		}

		finalLastVisTri = original_faces_number;
		if(visualizeTraversedTriangles){
			std::string file_name_trav_tri = "trav_tri_BDRY_"+std::to_string(contSource)+".vtk";
			writeTestMeshTrianglesIds(traversedSPTris[contSource][totalNumberOfSlots], file_name_trav_tri);
			if(withGlobalComments)
				cout<<"has written "<<file_name_trav_tri<<endl;

		}
	}
	if(withGlobalComments)
		cout<<"getShortestPathBtwFacesOptimized end"<<endl;

}

/*-----------------------------------------------------------------------------------------------------*/


//-------------------------------------------------------------------------------------------------
/*void SingularityGraphBuilder2D::getShortestPathBtwFacesOptimizedNewLocalMesh(gmds::TCellID &source,
                                            			vector<unsigned int> &targets,
                                           			vector<pair<double, unsigned int>> &min_distance,
                                          			vector<int> &previous,
								   			int &found,
											vector<vector<double>>& face2FaceMatch,
											vector<vector<unsigned int>>& face2FaceDeviation,
											pair<gmds::math::Vector3d, gmds::math::Vector3d>& prevDirCrossGlobal,
											vector<pair<gmds::math::Vector3d, gmds::math::Vector3d>>& targetDirCross,
											bool& searchForTargetBdry,
											vector<bool>& is_bdry_face,
											gmds::math::Point& startPoint,
											gmds::math::Point& finalEndPoint,
											double& maxDist,
											int& final_to_cell_dim,
											gmds::TCellID& final_to_cell_id,
											gmds::TCellID& finalLastVisTri,
											unsigned int& contSource,
											unsigned int& totalNumberOfVariables,
											unsigned int& totalNumberOfSlots,
											vector<unsigned int>& faceNo2Cont,
											vector<vector<vector<gmds::math::Point>>>& line_discretization,
											vector<vector<vector<gmds::TCellID>>>& line_triangles,
											vector<vector<pair<vector<gmds::math::Point>, vector<gmds::math::Point>>>>& pointPaths,
											vector<vector<vector<gmds::TCellID>>>& finalPaths,
											vector<vector<vector<gmds::math::Point>>>& finalCenterPoints,
											vector<vector<vector<gmds::TCellID>>>& traversedSPTris,
											vector<vector<double>>& distances,
											vector<pair<unsigned int, unsigned int>>& IllegalCross,
											vector<pair<vector<unsigned int>, vector<unsigned int >>>& isTraversedFaceCompVector,
											vector<pair<unsigned int, unsigned int>>& contSourceToSingularity,
											vector<bool>& singOrGeomFaces,
                                                       gmds::Mesh* newLocalMesh,
                                                       vector<gmds::TCellID>& newLocalMesh_id_to_mesh_id_node,
                                                       vector<gmds::TCellID>& mesh_id_to_newLocalMesh_id_node,
                                                       vector<gmds::TCellID>& newLocalMesh_id_to_mesh_id_face,
                                                       vector<gmds::TCellID>& mesh_id_to_newLocalMesh_id_face,
                                                       gmds::Variable<gmds::math::Cross2D>* newMesh_cross_field_2D,
                                                       vector<bool>& trianglesToRemeshBool,
											vector<gmds::math::Vector3d>& newBdryEdgeNormals,
                                                       vector<gmds::math::Vector3d>& newBdryNodeNormals,
											vector<bool>& isCurveEdge,
     										vector<bool>& isCurveNode,
											vector<gmds::math::Point>& newTriangleCenters,
                                                       vector<gmds::math::Cross2D>& newTriangleCenterCrosses)
{

  	// vector<vector<gmds::Face>> face2Face_neighbours_by_verts; - for each face , for each of its Nodes store adjacent edges
 	// vector<vector<double>> face2FaceMatch; - for each face, for each of its outgoing faces (face2Face_neighbours_by_verts - in the exact same order) - store the face 2 face matching

 // vector<vector<unsigned int>> face2FaceDeviation; - for each face, for each of its outgoing faces (face2Face_neighbours_by_verts - in the exact same order) - store the face 2 face deviation (of cross field associated to its center)
   	//TODO color tris according to distance in order to visualize
	//startPoint is in triangle source
  	//WARNING always by edge...otherwise in surface with holes (or non-convex bdry) it can try to pass through the hole

	//WARNING TODO also need newLocalMesh_id_to_mesh_id_face; take care of target faces!!!
 	bool visualizeTraversedTriangles = false;
 	bool visualizeAll = false; // including the lines for which the distance is greater than maxDist
 	unsigned int newFacesNumber = newLocalMesh->getNbFaces();
     gmds::TCellID newSource = mesh_id_to_newLocalMesh_id_face[source];
     vector<unsigned int> targetsTemp(targets.size());
     //WARNING TODO make sure the targetfaces are not!!! remeshed
	unsigned int contNonTarget = 0;
     for(unsigned int i=0; i<targets.size(); i++){
         	if(targets[i]>=original_faces_number){
			contNonTarget++;
              	targetsTemp[i] = newFacesNumber;
         	}
         	else
              	targetsTemp[i] = mesh_id_to_newLocalMesh_id_face[targets[i]];
     }

     cout<<"getShortestPathBtwFacesOptimizedNewLocalMesh"<<endl;

	cout<<"startPoint "<<startPoint<<endl;
     gmds::math::Point endPoint;
	int to_cell_dim = final_to_cell_dim;
	gmds::TCellID to_cell_id = final_to_cell_id;
	gmds::TCellID lastVisTri = finalLastVisTri;
     cout<<"finalEndPoint "<<finalEndPoint<<endl;

     cout<<"source tri  "<<source<<endl;

     vector<pair<gmds::math::Vector3d, gmds::math::Vector3d>> tempTargetDirCross(targetDirCross);

	bool withComments = false;
	unsigned int mycont = 0;

  	min_distance.clear();
  	min_distance.resize(newFacesNumber, make_pair(maxDist, 1));

  	min_distance[newSource] = make_pair(0.0,1);
     previous.clear();
 	previous.resize(newFacesNumber, -1);
	vector<unsigned int> paths2;
	unsigned int i_orig = (int) (contSource/5);
	unsigned int j_orig = fmod(contSource,5);

     double turnPenalizationTerm = 0.1 * maxDist;

	bool isTargetBdry;

	bool reachedBdry = false;

     std::set<std::pair<double,unsigned int> > vertex_queue;

	math::Vector3d closestToPrevDir, closestToPrevCross, prevDir;

	vector<pair<gmds::math::Vector3d, gmds::math::Vector3d>> tempPrevDirCross(newFacesNumber, prevDirCrossGlobal); // pair of previous direction and previous followed cross component

	for(unsigned int i=0; i<newFacesNumber; i++){
	    vertex_queue.insert(std::make_pair(min_distance[i].first,i));
     }

  	std::vector<unsigned int>::iterator it;
  	gmds::Face v, u;
	gmds::TCellID u_id, v_id;

	vector<bool> visitedFaces(newFacesNumber, false);

	cout<<"treat first point"<<endl;
	u = newLocalMesh->get<gmds::Face>(vertex_queue.begin()->second);
	u_id = u.id();
	double total_dist = vertex_queue.begin()->first;
	unsigned int visitedFacesNo = min_distance[vertex_queue.begin()->second].second;

	vertex_queue.erase(vertex_queue.begin());
	visitedFaces[u_id] = true;

	it = std::find (targetsTemp.begin(), targetsTemp.end(), u_id);
	math::Vector3d prevDirCross =  tempPrevDirCross[u_id].second;
	if(it!=targetsTemp.end()){

		found = std::distance(targetsTemp.begin(), it );
		prevDir =  tempPrevDirCross[u_id].first;

		if(withComments){
		 	cout<<"if(it!=targetsTemp.end()){"<<endl;
		 	cout<<"u_id "<<u_id<<endl;
		  	cout<<"prevDirCross "<<newTriangleCenters[u_id]<<" - > "<<prevDirCross[0] + newTriangleCenters[u_id].X()<<" , "<<prevDirCross[1] + newTriangleCenters[u_id].Y()<<endl;
			cout<<"tempTargetDirCross "<<" - > "<<newTriangleCenters[u_id].X() - tempTargetDirCross[it - targetsTemp.begin()].second[0]<<" , "<<newTriangleCenters[u_id].Y() - tempTargetDirCross[it - targetsTemp.begin()].second[1]<<newTriangleCenters[u_id]<<endl;
			cout<<"tri2tri "<<newTriangleCenters[previous[u_id]]<<" - > "<<newTriangleCenters[previous[u_id]] + newTriangleCenters[u_id]<<endl;
			math::Vector3d tri2tritest(newTriangleCenters[previous[u_id]], newTriangleCenters[u_id]);
			tri2tritest.normalize();
			cout<<"(tempTargetDirCross[it - targetsTemp.begin()].second).angle(tri2tritest) "<<((tempTargetDirCross[it - targetsTemp.begin()].second).angle(tri2tritest)* 180 / math::Constants::PI)<<endl;
			//if angle prevDir prevDirCRoss close to 45 deg - let it
			cout<<"(tempTargetDirCross[it - targetsTemp.begin()].second).angle(prevDirCross) "<<(tempTargetDirCross[it - targetsTemp.begin()].second).angle(prevDirCross)<<endl;
		}
		if((tempTargetDirCross[it - targetsTemp.begin()].second).angle(prevDirCross)>=gmds::math::Constants::PIDIV4)
              	min_distance[u_id].first =  total_dist + turnPenalizationTerm + 1.0;
		else{
		  	min_distance[u_id].first =  total_dist + 10*(tempTargetDirCross[it - targetsTemp.begin()].second).angle(prevDirCross);
			//return u_id;
		}
		/////////////////////////////////////////////////
		if((visualizeAll)||(min_distance[u_id].first<maxDist)){
		  	unsigned int original_u_id = newLocalMesh_id_to_mesh_id_face[u_id];
			distances[contSource][faceNo2Cont[original_u_id]] = (double)min_distance[u_id].first/(min_distance[u_id].second+1);
			pointPaths[contSource][faceNo2Cont[original_u_id]].first.insert(pointPaths[contSource][faceNo2Cont[original_u_id]].first.end(),line_discretization[i_orig][j_orig].begin(),line_discretization[i_orig][j_orig].end());

			retraceShortestPath(u_id, previous, paths2);
			reverse(paths2.begin(),paths2.end());
			finalPaths[contSource][faceNo2Cont[original_u_id]].clear();
			if(paths2.size()>1)
          		finalPaths[contSource][faceNo2Cont[original_u_id]].insert(finalPaths[contSource][faceNo2Cont[original_u_id]].end(),
                                    paths2.begin(),paths2.end()-1);
			unsigned int t1_orig = (int) (faceNo2Cont[original_u_id]/5);
			unsigned int t2_orig = fmod(faceNo2Cont[original_u_id],5);

			reverse(line_discretization[t1_orig][t2_orig].begin(), line_discretization[t1_orig][t2_orig].end());
          	pointPaths[contSource][faceNo2Cont[original_u_id]].second.insert(pointPaths[contSource][faceNo2Cont[original_u_id]].second.end(),line_discretization[t1_orig][t2_orig].begin(),line_discretization[t1_orig][t2_orig].end());

			vector<math::Point> centerPoints;
			centerPoints.insert(centerPoints.end(),line_discretization[i_orig][j_orig].begin(),line_discretization[i_orig][j_orig].end());

			for(unsigned int tt=0; tt<finalPaths[contSource][faceNo2Cont[original_u_id]].size(); tt++){
	   			centerPoints.push_back(triangle_centers[finalPaths[contSource][faceNo2Cont[original_u_id]][tt]]);
			}
			centerPoints.insert(centerPoints.end(),line_discretization[t1_orig][t2_orig].begin(),line_discretization[t1_orig][t2_orig].end());

			finalCenterPoints[contSource][5*t1_orig+t2_orig].insert(finalCenterPoints[contSource][5*t1_orig+t2_orig].end(), centerPoints.begin(), centerPoints.end());
			std::string file_name = "ShortestPaths_"+std::to_string(contSource)+"-"+std::to_string(faceNo2Cont[original_u_id])+".vtk";
			writeTestPoints(centerPoints, file_name);
			if(withComments)
				cout<<"has written "<<file_name<<endl;

			isTargetBdry = false;

			unsigned int contMatrix = contSource*totalNumberOfVariables + faceNo2Cont[original_u_id];
			contSourceToSingularity[contMatrix].first = i_orig;
			contSourceToSingularity[contMatrix].second = t1_orig;

			reverse(line_triangles[t1_orig][t2_orig].begin(),line_triangles[t1_orig][t2_orig].end());

			if(min_distance[u_id].first<maxDist){
				getTraversedTrisFaceNeighbyVerts(line_discretization[i_orig][j_orig],
										line_triangles[i_orig][j_orig],
									 	paths2,
									 	isTargetBdry,
									 	line_discretization[t1_orig][t2_orig],
								 		line_triangles[t1_orig][t2_orig],
										traversedSPTris[contSource][faceNo2Cont[original_u_id]],
										isTraversedFaceCompVector,
									 	contMatrix,
									 	totalNumberOfVariables,
										IllegalCross,
									 	contSourceToSingularity,
									 	finalLastVisTri,
									 	endPoint);
			}
			reverse(line_triangles[t1_orig][t2_orig].begin(),line_triangles[t1_orig][t2_orig].end());
			reverse(line_discretization[t1_orig][t2_orig].begin(),line_discretization[t1_orig][t2_orig].end());
			if(visualizeTraversedTriangles){
				std::string file_name_trav_tri = "trav_tri_"+std::to_string(contSource)+"-"+std::to_string(faceNo2Cont[original_u_id])+".vtk";
				writeTestMeshTrianglesIds(traversedSPTris[contSource][faceNo2Cont[original_u_id]], file_name_trav_tri);
			}
			/////////////////////////////////////////////////////////////////////
			int index = std::distance(targetsTemp.begin(), it);
  			targetsTemp.erase(targetsTemp.begin()+index);
			tempTargetDirCross.erase(tempTargetDirCross.begin()+index);
		}
	}
	else{//WARNING what if target face is also a bdry face? treat the case
		if((is_bdry_face[u_id])&&(searchForTargetBdry)){

			// math::Point intersectionPnt;
			to_cell_dim = 2;
			to_cell_id = u_id;
			double intersectionParam;
			math::Ray from_ray(startPoint, startPoint + prevDirCross);
			if(withComments){
				cout<<"u_id "<<u_id<<endl;
			  	cout<<"from_ray "<<startPoint<<" -> "<<startPoint + prevDirCross<<endl;
			  	cout<<"prevDirCross "<<newTriangleCenters[u_id]<<" - > "<<prevDirCross[0] + newTriangleCenters[u_id].X()<<" , "<<prevDirCross[1] + triangle_centers[u_id].Y()<<endl;
			}

			vector<gmds::Edge> currentEdges = u.get<gmds::Edge>();
			bool hasBdryEdge = false;
			for(unsigned int i=0; i<currentEdges.size(); i++){
				if (m_mesh->isMarked(currentEdges[i], m_mark_edges_on_curve)){
					hasBdryEdge = true;
					to_cell_dim = 1;
					to_cell_id = currentEdges[i].id();
					vector<gmds::Node> currentNodes =  currentEdges[i].get<Node>();
					math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
					if(from_ray.SecondMetIntersect2D(oppSeg, endPoint, intersectionParam, temp_epsilon)){
				  		double tempAngleVal = 1.0 - fabs(prevDirCross.dot(newBdryEdgeNormals[currentEdges[i].id()]));
						if(tempAngleVal<gmds::math::Constants::PIDIV4){
					 		min_distance[u_id].first =  total_dist + tempAngleVal;
							reachedBdry = true;
						}
						if(withComments)
							cout<<"endPoint "<<endPoint<<endl;
					}
				}
			}
			if(!hasBdryEdge){
				// VERTEX-> angle with vert normal
				vector<gmds::Node> currentNodes = u.get<gmds::Node>();
				for(unsigned int i=0; i<currentNodes.size(); i++){
					if(m_mesh->isMarked(currentNodes[i], m_mark_nodes_on_point) ||
           		   	m_mesh->isMarked(currentNodes[i], m_mark_nodes_on_curve)){
				  		to_cell_dim = 0;
				 		to_cell_id = currentNodes[i].id();
 				  		vector<gmds::Edge> currentEdges = currentNodes[i].get<gmds::Edge>();
						for(unsigned int j=0; j<currentEdges.size(); j++){
					  		if(m_mesh->isMarked(currentEdges[j], m_mark_edges_on_curve)){
						  		vector<gmds::Node> edge_nodes = currentEdges[j].get<gmds::Node>();
								math::Segment oppSeg(edge_nodes[0].getPoint(), edge_nodes[1].getPoint());
								if(from_ray.SecondMetIntersect2D(oppSeg, endPoint, intersectionParam, temp_epsilon)){
							  		to_cell_dim = 1;
									to_cell_id = currentEdges[j].id();
				  					double tempAngleVal = 1.0 - fabs(prevDirCross.dot(bdry_edge_normals[currentEdges[j].id()]));

					 				if(tempAngleVal<gmds::math::Constants::PIDIV4){
								  		min_distance[u_id].first =  total_dist + tempAngleVal;
										vector<gmds::Face> adj_faces1 = currentEdges[j].get<gmds::Face>();

										reachedBdry = true;
										lastVisTri = adj_faces1[0].id();
							    		}
							    		if(withComments)
							    			cout<<"endPoint "<<endPoint<<endl;
								}
							}
						}
					}
				}
			}

    			/////////////////////////////////////////////////
    			if(reachedBdry){
		  		if((min_distance[u_id].first<maxDist)){
				  	retraceShortestPath( u_id, previous, paths2);
					bool toConsider = true;
					for(unsigned int tk=0; tk<paths2.size(); tk++){
					  	if(singOrGeomFaces[paths2[tk]])
					    		toConsider = false;
				    	}
				    	if(toConsider){
					  	finalEndPoint = endPoint;
					  	finalLastVisTri = lastVisTri;
					  	final_to_cell_dim = to_cell_dim;
						final_to_cell_id = to_cell_id;

						distances[contSource][totalNumberOfSlots] = (double)min_distance[u_id].first/(min_distance[u_id].second);
						pointPaths[contSource][totalNumberOfSlots].first.clear();
						pointPaths[contSource][totalNumberOfSlots].first.insert(pointPaths[contSource][totalNumberOfSlots].first.end(),line_discretization[i_orig][j_orig].begin(),line_discretization[i_orig][j_orig].end()-1);


						reverse(paths2.begin(),paths2.end());
						finalPaths[contSource][totalNumberOfSlots].clear();
						if(paths2.size()>1)
          					finalPaths[contSource][totalNumberOfSlots].insert(finalPaths[contSource][totalNumberOfSlots].end(),
                                   	 paths2.begin(),paths2.end());


					}

				}
    			}
		}
	}
		/////////////

	bool toCheck = false;
	if(is_bdry_face[u_id])
	 	toCheck = true;

	for(unsigned int i=0; i<face2Face_neighbours_by_verts[u_id].size(); i++){
		v = face2Face_neighbours_by_verts[u_id][i];
		v_id = v.id();

		bool validNeighbour = true;
		if(toCheck==true){
	  		if(is_bdry_face[v_id]){
				math::Point intersectionPnt;
				double intersectionParam;
				math::Ray from_ray(triangle_centers[u_id], triangle_centers[v_id]);
				vector<gmds::Edge> currentEdges = u.get<gmds::Edge>();
				for(unsigned int j=0; j<3; j++){
					if (m_mesh->isMarked(currentEdges[j], m_mark_edges_on_curve)){
			  			vector<gmds::Node> currentNodes =  currentEdges[j].get<Node>();
						math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
						if(from_ray.SecondMetIntersect2D(oppSeg, intersectionPnt, intersectionParam, temp_epsilon)){
				  			if((intersectionParam>0.000000001)||(intersectionParam<0.999999999)){
					  			validNeighbour = false;
							}
						}
					}
				}
	 		}
		}

		if((!visitedFaces[v_id])&&(validNeighbour)){

			math::Vector3d tri2tri(startPoint, triangle_centers[v_id]);
			tri2tri.normalize();
			math::Cross2D crossV = triangle_centers_cross[v_id];
			prevDir =  tempPrevDirCross[u_id].first;
			closestToPrevDir =  crossV.closestComponentVector(prevDir);
			closestToPrevCross =  crossV.closestComponentVector(tempPrevDirCross[u_id].second);

               double distance_through_u;

			if(withComments){
				cout<<"withComments"<<endl;
			  	cout<<"u_id "<<u_id<<endl;
			   	cout<<"v_id "<<v_id<<endl;
			  	cout<<"tri2tri "<<startPoint<<" -> "<<triangle_centers[v_id]<<endl;

			  	cout<<"prevDirCross "<<triangle_centers[u_id]<<" - > "<<tempPrevDirCross[u_id].second[0] + triangle_centers[u_id].X()<<" , "<<tempPrevDirCross[u_id].second[1] + triangle_centers[u_id].Y()<<endl;
				cout<<"closestToPrevCross "<<triangle_centers[v_id]<<" - > "<<closestToPrevCross[0] + triangle_centers[v_id].X()<<" , "<<closestToPrevCross[1] + triangle_centers[v_id].Y()<<endl;
				cout<<"tri2tri.angle(tempPrevDirCross[u_id].second) "<<tri2tri.angle(tempPrevDirCross[u_id].second)* 180 / math::Constants::PI<<endl;
			}

			if(tri2tri.angle(tempPrevDirCross[u_id].second)>=gmds::math::Constants::PIDIV4)
               	distance_through_u =  total_dist + maxDist + 1.0;
 			else
				distance_through_u =  total_dist + //(tempPrevDirCross[u_id].second).angle(tri2tri)
			    							   + closestToPrevCross.angle(tri2tri)
										   + (floor)((tempPrevDirCross[u_id].second.angle(closestToPrevCross))/gmds::math::Constants::PIDIV2)*10.0;


			if (distance_through_u/(visitedFacesNo+1) < min_distance[v_id].first/min_distance[v_id].second) {
			//if (distance_through_u < min_distance[v_id].first) {

				vertex_queue.erase(std::make_pair(min_distance[v_id].first,v_id));

				min_distance[v_id].first = distance_through_u;
				min_distance[v_id].second = visitedFacesNo + 1;
				previous[v_id] = u_id;
				vertex_queue.insert(std::make_pair(min_distance[v_id].first,v_id));
				tempPrevDirCross[v_id].first = tri2tri;
				tempPrevDirCross[v_id].second = closestToPrevCross;

			}

		}

	}
	///////////////////////////////////////////////////////////////////////////////////////////
	if(withComments)
		cout<<"while ((!vertex_queue.empty())&&(!targetsTemp.empty())){"<<endl;


	while ((!vertex_queue.empty())&&(targetsTemp.size()>contNonTarget)){

		u = m_mesh->get<gmds::Face>(vertex_queue.begin()->second);
		u_id = u.id();

		double total_dist = vertex_queue.begin()->first;
		unsigned int visitedFacesNo = min_distance[vertex_queue.begin()->second].second;

		vertex_queue.erase(vertex_queue.begin());
		visitedFaces[u_id] = true;

		it = std::find (targetsTemp.begin(), targetsTemp.end(), u_id);
		math::Vector3d prevDirCross =  tempPrevDirCross[u_id].second;
		if(it!=targetsTemp.end()){

			found = std::distance( targetsTemp.begin(), it );
			prevDir =  tempPrevDirCross[u_id].first;

			if(withComments){
			   	cout<<"if(it!=targetsTemp.end()){"<<endl;
				cout<<"u_id "<<u_id<<endl;
				cout<<"prevDirCross "<<triangle_centers[u_id].X()<<" , "<<triangle_centers[u_id].Y()<<" - > "<<prevDirCross[0] + triangle_centers[u_id].X()<<" , "<<prevDirCross[1] + triangle_centers[u_id].Y()<<endl;

				cout<<"is_bdry_face[u_id] "<<is_bdry_face[u_id]<<endl;
				cout<<"tempTargetDirCross "<<" - > "<<triangle_centers[targetsTemp[0]].X() - tempTargetDirCross[it - targetsTemp.begin()].second[0]<<" , "<<triangle_centers[targetsTemp[0]].Y() - tempTargetDirCross[it - targetsTemp.begin()].second[1]<<triangle_centers[targetsTemp[0]]<<endl;
				cout<<"tri2tri "<<triangle_centers[previous[u_id]]<<" - > "<<triangle_centers[u_id]<<endl;
				math::Vector3d tri2tritest(triangle_centers[previous[u_id]], triangle_centers[u_id]);
				tri2tritest.normalize();
				cout<<"tempTargetDirCross.angle(tri2tritest) "<<((tempTargetDirCross[it - targetsTemp.begin()].second).angle(tri2tritest)* 180 / math::Constants::PI)<<endl;

				cout<<"(tempTargetDirCross[it - targetsTemp.begin()].second).angle(prevDirCross) "<<(tempTargetDirCross[it - targetsTemp.begin()].second).angle(prevDirCross)<<endl;
			}

			if((tempTargetDirCross[it - targetsTemp.begin()].second).angle(prevDirCross)>=gmds::math::Constants::PIDIV4)
               	min_distance[u_id].first =  total_dist + turnPenalizationTerm + 1.0;
			else{
			  	min_distance[u_id].first =  total_dist + 10*(tempTargetDirCross[it - targetsTemp.begin()].second).angle(prevDirCross);
				 							  //+10*(tempTargetDirCross[it - targetsTemp.begin()].first).angle(prevDir);
			}
			/////////////////////////////////////////////////
			if((visualizeAll)||(min_distance[u_id].first<maxDist)){
				distances[contSource][faceNo2Cont[u_id]] = (double)min_distance[u_id].first/(min_distance[u_id].second+1);
				pointPaths[contSource][faceNo2Cont[u_id]].first.insert(pointPaths[contSource][faceNo2Cont[u_id]].first.end(),line_discretization[i_orig][j_orig].begin(),line_discretization[i_orig][j_orig].end());

				retraceShortestPath(u_id, previous, paths2);
				reverse(paths2.begin(),paths2.end());
				finalPaths[contSource][faceNo2Cont[u_id]].clear();
				if(paths2.size()>1)
          			finalPaths[contSource][faceNo2Cont[u_id]].insert(finalPaths[contSource][faceNo2Cont[u_id]].end(),
                                    paths2.begin(),paths2.end()-1);

				unsigned int t1_orig = (int) (faceNo2Cont[u_id]/5);
				unsigned int t2_orig = fmod(faceNo2Cont[u_id],5);

				reverse(line_discretization[t1_orig][t2_orig].begin(), line_discretization[t1_orig][t2_orig].end());
          		pointPaths[contSource][faceNo2Cont[u_id]].second.insert(pointPaths[contSource][faceNo2Cont[u_id]].second.end(),line_discretization[t1_orig][t2_orig].begin(),line_discretization[t1_orig][t2_orig].end());

				vector<math::Point> centerPoints;
				centerPoints.insert(centerPoints.end(),line_discretization[i_orig][j_orig].begin(),line_discretization[i_orig][j_orig].end());

				for(unsigned int tt=0; tt<finalPaths[contSource][faceNo2Cont[u_id]].size(); tt++){
	   				centerPoints.push_back(triangle_centers[finalPaths[contSource][faceNo2Cont[u_id]][tt]]);
				}
				centerPoints.insert(centerPoints.end(),line_discretization[t1_orig][t2_orig].begin(),line_discretization[t1_orig][t2_orig].end());

				finalCenterPoints[contSource][5*t1_orig+t2_orig].insert(finalCenterPoints[contSource][5*t1_orig+t2_orig].end(), centerPoints.begin(), centerPoints.end());
				std::string file_name = "ShortestPaths_"+std::to_string(contSource)+"-"+std::to_string(faceNo2Cont[u_id])+".vtk";
				writeTestPoints(centerPoints, file_name);
				cout<<"has written "<<file_name<<endl;

				isTargetBdry = false;

				unsigned int contMatrix = contSource*totalNumberOfVariables + faceNo2Cont[u_id];
				contSourceToSingularity[contMatrix].first = i_orig;
				contSourceToSingularity[contMatrix].second = t1_orig;

				reverse(line_triangles[t1_orig][t2_orig].begin(),line_triangles[t1_orig][t2_orig].end());

				if(min_distance[u_id].first<maxDist){
			  		finalLastVisTri = newFacesNumber;

					getTraversedTrisFaceNeighbyVerts(line_discretization[i_orig][j_orig],
										line_triangles[i_orig][j_orig],
									 	paths2,
									 	isTargetBdry,
									 	line_discretization[t1_orig][t2_orig],
								 		line_triangles[t1_orig][t2_orig],
										traversedSPTris[contSource][faceNo2Cont[u_id]],
										isTraversedFaceCompVector,
									 	contMatrix,
									 	totalNumberOfVariables,
										IllegalCross,
									 	contSourceToSingularity,
									 	finalLastVisTri,
									 	endPoint);
				}

				reverse(line_triangles[t1_orig][t2_orig].begin(),line_triangles[t1_orig][t2_orig].end());
				reverse(line_discretization[t1_orig][t2_orig].begin(),line_discretization[t1_orig][t2_orig].end());

				if(visualizeTraversedTriangles){
					std::string file_name_trav_tri = "trav_tri_"+std::to_string(contSource)+"-"+std::to_string(faceNo2Cont[u_id])+".vtk";
					writeTestMeshTrianglesIds(traversedSPTris[contSource][faceNo2Cont[u_id]], file_name_trav_tri);
				}

				/////////////////////////////////////////////////////////////////////

				int index = std::distance(targetsTemp.begin(), it);

  				targetsTemp.erase(targetsTemp.begin()+index);

				tempTargetDirCross.erase(tempTargetDirCross.begin()+index);
			}

		}

		if((is_bdry_face[u_id])&&(searchForTargetBdry)){ //(!reachedBdry)

			to_cell_dim = 2;
			to_cell_id = u_id;
			//math::Point intersectionPnt;

			double intersectionParam;
			math::Ray from_ray(triangle_centers[u_id], prevDirCross);
			if(withComments){
			  	cout<<"if((!reachedBdry)&&(is_bdry_face[u_id])&&(searchForTargetBdry)){"<<endl;
				cout<<"target is bdry"<<endl;
				cout<<"from_ray "<<triangle_centers[u_id]<<" -> "<<triangle_centers[u_id].X()+prevDirCross[0]<<" , "<<triangle_centers[u_id].Y()+prevDirCross[1]<<endl;
			}

			 vector<gmds::Edge> currentEdges = u.get<gmds::Edge>();
			 bool hasBdryEdge = false;
			 reachedBdry = false;
			 for(unsigned int i=0; i<currentEdges.size(); i++){

				if (m_mesh->isMarked(currentEdges[i], m_mark_edges_on_curve)){
					hasBdryEdge = true;
					to_cell_dim = 1;
					to_cell_id = currentEdges[i].id();

					vector<gmds::Node> currentNodes =  currentEdges[i].get<Node>();
					math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());


					if(from_ray.SecondMetIntersect2D(oppSeg, endPoint, intersectionParam, temp_epsilon)){
					  	double tempAngleVal = 1.0 - fabs(prevDirCross.dot(bdry_edge_normals[currentEdges[i].id()]));

						if(tempAngleVal<gmds::math::Constants::PIDIV4){
						 	min_distance[u_id].first =  total_dist + tempAngleVal;

							if(withComments)
								cout<<"1 u_id "<<u_id<<" , endPoint "<<endPoint<<endl;
							reachedBdry = true;
						}
					}
				}
			 }
			 if(!hasBdryEdge){
			 	// VERTEX-> angle with vert normal
			    	if(withComments){
					cout<<"if(!hasBdryEdge){"<<endl;
				}
			  	vector<gmds::Node> currentNodes = u.get<gmds::Node>();
				for(unsigned int i=0; i<currentNodes.size(); i++){
					if(m_mesh->isMarked(currentNodes[i], m_mark_nodes_on_point) ||
           			   m_mesh->isMarked(currentNodes[i], m_mark_nodes_on_curve)){
					  	to_cell_dim = 0;
					 	to_cell_id = currentNodes[i].id();
 					  	vector<gmds::Edge> currentEdges = currentNodes[i].get<gmds::Edge>();
						for(unsigned int j=0; j<currentEdges.size(); j++){
						  	if(m_mesh->isMarked(currentEdges[j], m_mark_edges_on_curve)){
							  	vector<gmds::Node> edge_nodes = currentEdges[j].get<gmds::Node>();
								math::Segment oppSeg(edge_nodes[0].getPoint(), edge_nodes[1].getPoint());

								if(from_ray.SecondMetIntersect2D(oppSeg, endPoint, intersectionParam, temp_epsilon)){
					  				double tempAngleVal = 1.0 - fabs(prevDirCross.dot(bdry_edge_normals[currentEdges[j].id()]));
									to_cell_dim = 1;
					 				to_cell_id = currentEdges[j].id();
									//<0.5
								 	if(tempAngleVal<gmds::math::Constants::PIDIV4){
									  	min_distance[u_id].first =  total_dist + tempAngleVal;
										vector<gmds::Face> adj_faces1 = currentEdges[j].get<gmds::Face>();
					     				lastVisTri = adj_faces1[0].id();
										reachedBdry = true;
								    }
								}
							}
						}
					}
			  	}
			 }
			 if(reachedBdry){
			   	if(distances[contSource][totalNumberOfSlots] > (double)min_distance[u_id].first/(min_distance[u_id].second)){
					retraceShortestPath( u_id, previous, paths2);
				  	bool toConsider = true;
					for(unsigned int tk=0; tk<paths2.size(); tk++){
					  	if(singOrGeomFaces[paths2[tk]])
					    		toConsider = false;
				    	}
				    	if(toConsider){
				  		finalEndPoint = endPoint;
					  	finalLastVisTri = lastVisTri;
					  	final_to_cell_dim = to_cell_dim;
						final_to_cell_id = to_cell_id;
				  		distances[contSource][totalNumberOfSlots] = (double)min_distance[u_id].first/(min_distance[u_id].second);
						pointPaths[contSource][totalNumberOfSlots].first.clear();
						pointPaths[contSource][totalNumberOfSlots].first.insert(pointPaths[contSource][totalNumberOfSlots].first.end(),line_discretization[i_orig][j_orig].begin(),line_discretization[i_orig][j_orig].end()-1);

						reverse(paths2.begin(),paths2.end());
						finalPaths[contSource][totalNumberOfSlots].clear();
						if(paths2.size()>1)
          					finalPaths[contSource][totalNumberOfSlots].insert(finalPaths[contSource][totalNumberOfSlots].end(),
                         	           paths2.begin(),paths2.end());
					}
				}
			 }
		}

		toCheck = false;
		if(is_bdry_face[u_id])
		  	toCheck = true;

		for(unsigned int i=0; i<face2Face_neighbours_by_verts[u_id].size(); i++){
			v = face2Face_neighbours_by_verts[u_id][i];
			v_id = v.id();
			//cout<<"v_id "<<v_id<<endl;
			bool validNeighbour = true;
			if(toCheck==true){
		  		if(is_bdry_face[v_id]){
					math::Point intersectionPnt;
					double intersectionParam;
					math::Ray from_ray(triangle_centers[u_id], triangle_centers[v_id]);
					vector<gmds::Edge> currentEdges = u.get<gmds::Edge>();
					for(unsigned int j=0; j<3; j++){
						if (m_mesh->isMarked(currentEdges[j], m_mark_edges_on_curve)){
				  			vector<gmds::Node> currentNodes =  currentEdges[j].get<Node>();
							math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
							if(from_ray.SecondMetIntersect2D(oppSeg, intersectionPnt, intersectionParam, temp_epsilon)){
					  			if((intersectionParam>0.000000001)||(intersectionParam<0.999999999)){
						  			validNeighbour = false;
								}
							}
						}
					}
		 		}
			}
			if((!visitedFaces[v_id])&&(validNeighbour)){

				math::Vector3d tri2tri(triangle_centers[u_id], triangle_centers[v_id]);
				tri2tri.normalize();
				//cout<<"valid v_id "<<v_id<<endl;
				math::Cross2D crossV = triangle_centers_cross[v_id];
				prevDir =  tempPrevDirCross[u_id].first;
	 			closestToPrevDir =  crossV.closestComponentVector(prevDir);
				closestToPrevCross =  crossV.closestComponentVector(tempPrevDirCross[u_id].second);

				//closestToPrevCross =  crossV.closestComponentVector(tri2tri);

               	double distance_through_u;

				if(tri2tri.angle(tempPrevDirCross[u_id].second)>=gmds::math::Constants::PIDIV4)
                  		distance_through_u =  total_dist + turnPenalizationTerm + 1.0;
 				else
			  		distance_through_u =  total_dist + //(tempPrevDirCross[u_id].second).angle(tri2tri)
				    							   + closestToPrevCross.angle(tri2tri)
											   + (floor)((tempPrevDirCross[u_id].second.angle(closestToPrevCross))/gmds::math::Constants::PIDIV2)*10.0;

			    	if(withComments){
				  	cout<<"withComments"<<endl;
			  		cout<<"u_id "<<u_id<<endl;
			   		cout<<"v_id "<<v_id<<endl;
			  		cout<<"tri2tri "<<triangle_centers[u_id]<<" -> "<<triangle_centers[v_id]<<endl;

			  		cout<<"prevDirCross "<<triangle_centers[u_id]<<" - > "<<tempPrevDirCross[u_id].second[0] + triangle_centers[u_id].X()<<" , "<<tempPrevDirCross[u_id].second[1] + triangle_centers[u_id].Y()<<endl;
					cout<<"closestToPrevCross "<<triangle_centers[v_id]<<" - > "<<closestToPrevCross[0] + triangle_centers[v_id].X()<<" , "<<closestToPrevCross[1] + triangle_centers[v_id].Y()<<endl;
					cout<<"tri2tri.angle(tempPrevDirCross[u_id].second) "<<tri2tri.angle(tempPrevDirCross[u_id].second)* 180 / math::Constants::PI<<endl;
					cout<<"distance_through_u "<<distance_through_u<<endl;
				}


               	//WARNING division par visitedFacesNo not particularly ok; especially for very irregular meshes
				if (distance_through_u/(visitedFacesNo+1) < min_distance[v_id].first/min_distance[v_id].second) {
				//if (distance_through_u < min_distance[v_id].first) {
			  		if(withComments)
			    			cout<<"if distance_through_u<"<<endl;
					vertex_queue.erase(std::make_pair(min_distance[v_id].first,v_id));

					min_distance[v_id].first = distance_through_u;
					min_distance[v_id].second = visitedFacesNo + 1;
					previous[v_id] = u_id;
					vertex_queue.insert(std::make_pair(min_distance[v_id].first,v_id));

					tempPrevDirCross[v_id].first = tri2tri;
					tempPrevDirCross[v_id].second = closestToPrevCross;
				}

			}

		}
	}

 	if(withComments){
 		cout<<" while((!vertex_queue.empty())&&(!reachedBdry)&&(searchForTargetBdry)){"<<endl;
		cout<<"reachedBdry "<<reachedBdry<<endl;
		cout<<"searchForTargetBdry "<<searchForTargetBdry<<endl;
		cout<<"vertex_queue.size() "<<vertex_queue.size()<<endl;
	}
 	while((!vertex_queue.empty())&&(searchForTargetBdry)){//&&(!reachedBdry)

		u = m_mesh->get<gmds::Face>(vertex_queue.begin()->second);
		u_id = u.id();

		double total_dist = vertex_queue.begin()->first;
		unsigned int visitedFacesNo = min_distance[vertex_queue.begin()->second].second;

		vertex_queue.erase(vertex_queue.begin());
		visitedFaces[u_id] = true;

		math::Vector3d prevDirCross =  tempPrevDirCross[u_id].second;

		if(is_bdry_face[u_id]){//target is bdry
			to_cell_dim = 2;
			to_cell_id = u_id;
			//math::Point intersectionPnt;

			double intersectionParam;
			math::Ray from_ray(triangle_centers[u_id], prevDirCross);
			if(withComments){
			  	cout<<"u_id "<<u_id<<endl;
				cout<<"target is bdry"<<endl;
				cout<<"from_ray "<<triangle_centers[u_id]<<" -> "<<triangle_centers[u_id].X()+prevDirCross[0]<<" , "<<triangle_centers[u_id].Y()+prevDirCross[1]<<endl;
			}

			 vector<gmds::Edge> currentEdges = u.get<gmds::Edge>();
			 bool hasBdryEdge = false;
			 reachedBdry = false;
			 for(unsigned int i=0; i<currentEdges.size(); i++){

				if (m_mesh->isMarked(currentEdges[i], m_mark_edges_on_curve)){
					hasBdryEdge = true;
					to_cell_dim = 1;
					to_cell_id = currentEdges[i].id();

					vector<gmds::Node> currentNodes =  currentEdges[i].get<Node>();
					math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());


					if(from_ray.SecondMetIntersect2D(oppSeg, endPoint, intersectionParam, temp_epsilon)){
					  	double tempAngleVal = 1.0 - fabs(prevDirCross.dot(bdry_edge_normals[currentEdges[i].id()]));

						if(tempAngleVal<gmds::math::Constants::PIDIV4){
						 	min_distance[u_id].first =  total_dist + tempAngleVal;
							if(withComments)
								cout<<"1 u_id "<<u_id<<" , endPoint "<<endPoint<<endl;
							reachedBdry = true;
						}
					}
				}
			 }
			 if(!hasBdryEdge){
			 	// VERTEX-> angle with vert normal
			    	if(withComments){
					cout<<"if(!hasBdryEdge){"<<endl;
				}
			  	vector<gmds::Node> currentNodes = u.get<gmds::Node>();
				for(unsigned int i=0; i<currentNodes.size(); i++){
					if(m_mesh->isMarked(currentNodes[i], m_mark_nodes_on_point) ||
           			   m_mesh->isMarked(currentNodes[i], m_mark_nodes_on_curve)){
					  	to_cell_dim = 0;
					 	to_cell_id = currentNodes[i].id();
 					  	vector<gmds::Edge> currentEdges = currentNodes[i].get<gmds::Edge>();
						for(unsigned int j=0; j<currentEdges.size(); j++){
						  	if(m_mesh->isMarked(currentEdges[j], m_mark_edges_on_curve)){
							  	vector<gmds::Node> edge_nodes = currentEdges[j].get<gmds::Node>();
								math::Segment oppSeg(edge_nodes[0].getPoint(), edge_nodes[1].getPoint());

								if(from_ray.SecondMetIntersect2D(oppSeg, endPoint, intersectionParam, temp_epsilon)){
					  				double tempAngleVal = 1.0 - fabs(prevDirCross.dot(bdry_edge_normals[currentEdges[j].id()]));
									to_cell_dim = 1;
					 				to_cell_id = currentEdges[j].id();
									//<0.5
								 	if(tempAngleVal<gmds::math::Constants::PIDIV4){
									  	min_distance[u_id].first =  total_dist + tempAngleVal;
										vector<gmds::Face> adj_faces1 = currentEdges[j].get<gmds::Face>();
					     				lastVisTri = adj_faces1[0].id();
										reachedBdry = true;
								    }
								}
							}
						}
					}
			  	}
			 }
			 if(reachedBdry){
				if(distances[contSource][totalNumberOfSlots] > (double)min_distance[u_id].first/(min_distance[u_id].second)){
					retraceShortestPath( u_id, previous, paths2);
				  	bool toConsider = true;
					for(unsigned int tk=0; tk<paths2.size(); tk++){
					  	if(singOrGeomFaces[paths2[tk]])
					    		toConsider = false;
				    	}
				    	if(toConsider){
				  		finalEndPoint = endPoint;
					  	finalLastVisTri = lastVisTri;
					  	final_to_cell_dim = to_cell_dim;
						final_to_cell_id = to_cell_id;
				  		distances[contSource][totalNumberOfSlots] = (double)min_distance[u_id].first/(min_distance[u_id].second);
						pointPaths[contSource][totalNumberOfSlots].first.clear();
						pointPaths[contSource][totalNumberOfSlots].first.insert(pointPaths[contSource][totalNumberOfSlots].first.end(),line_discretization[i_orig][j_orig].begin(),line_discretization[i_orig][j_orig].end()-1);

						reverse(paths2.begin(),paths2.end());
						finalPaths[contSource][totalNumberOfSlots].clear();
						if(paths2.size()>1)
          					finalPaths[contSource][totalNumberOfSlots].insert(finalPaths[contSource][totalNumberOfSlots].end(),
                                   	 paths2.begin(),paths2.end());
				    }
			 	}
			 }
		}

		toCheck = false;
		if(is_bdry_face[u_id])
		  	toCheck = true;

		for(unsigned int i=0; i<face2Face_neighbours_by_verts[u_id].size(); i++){
			v = face2Face_neighbours_by_verts[u_id][i];
			v_id = v.id();
			bool validNeighbour = true;
			if(toCheck==true){
		  		if(is_bdry_face[v_id]){
					math::Point intersectionPnt;
					double intersectionParam;
					math::Ray from_ray(triangle_centers[u_id], triangle_centers[v_id]);
					vector<gmds::Edge> currentEdges = u.get<gmds::Edge>();
					for(unsigned int j=0; j<3; j++){
						if (m_mesh->isMarked(currentEdges[j], m_mark_edges_on_curve)){
				  			vector<gmds::Node> currentNodes =  currentEdges[j].get<Node>();
							math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
							if(from_ray.SecondMetIntersect2D(oppSeg, intersectionPnt, intersectionParam, temp_epsilon)){
					  			if((intersectionParam>0.000000001)||(intersectionParam<0.999999999)){
						  			validNeighbour = false;
								}
							}
						}
					}
		 		}
			}
			if((!visitedFaces[v_id])&&(validNeighbour)){

				math::Vector3d tri2tri(triangle_centers[u_id], triangle_centers[v_id]);
				tri2tri.normalize();
				//cout<<"valid v_id "<<v_id<<endl;
				math::Cross2D crossV = triangle_centers_cross[v_id];
				prevDir =  tempPrevDirCross[u_id].first;
	 			closestToPrevDir =  crossV.closestComponentVector(prevDir);
				closestToPrevCross =  crossV.closestComponentVector(tempPrevDirCross[u_id].second);

				//closestToPrevCross =  crossV.closestComponentVector(tri2tri);

               	double distance_through_u;

				if(tri2tri.angle(tempPrevDirCross[u_id].second)>=gmds::math::Constants::PIDIV4)
                  		distance_through_u =  total_dist + turnPenalizationTerm + 1.0;
 				else
			  		distance_through_u =  total_dist + //(tempPrevDirCross[u_id].second).angle(tri2tri)
				    							   + closestToPrevCross.angle(tri2tri)
											   + (floor)((tempPrevDirCross[u_id].second.angle(closestToPrevCross))/gmds::math::Constants::PIDIV2)*10.0;

			    	if(withComments){
				  	cout<<"withComments"<<endl;
			  		cout<<"u_id "<<u_id<<endl;
			   		cout<<"v_id "<<v_id<<endl;
			  		cout<<"tri2tri "<<triangle_centers[u_id]<<" -> "<<triangle_centers[v_id]<<endl;

			  		cout<<"prevDirCross "<<triangle_centers[u_id]<<" - > "<<tempPrevDirCross[u_id].second[0] + triangle_centers[u_id].X()<<" , "<<tempPrevDirCross[u_id].second[1] + triangle_centers[u_id].Y()<<endl;
					cout<<"closestToPrevCross "<<triangle_centers[v_id]<<" - > "<<closestToPrevCross[0] + triangle_centers[v_id].X()<<" , "<<closestToPrevCross[1] + triangle_centers[v_id].Y()<<endl;
					cout<<"tri2tri.angle(tempPrevDirCross[u_id].second) "<<tri2tri.angle(tempPrevDirCross[u_id].second)* 180 / math::Constants::PI<<endl;
					cout<<"distance_through_u "<<distance_through_u<<endl;
				}


               	//WARNING division par visitedFacesNo not particularly ok; especially for very irregular meshes
				if (distance_through_u/(visitedFacesNo+1) < min_distance[v_id].first/min_distance[v_id].second) {
				//if (distance_through_u < min_distance[v_id].first) {
			  		if(withComments)
			    			cout<<"if distance_through_u<"<<endl;
					vertex_queue.erase(std::make_pair(min_distance[v_id].first,v_id));

					min_distance[v_id].first = distance_through_u;
					min_distance[v_id].second = visitedFacesNo + 1;
					previous[v_id] = u_id;
					vertex_queue.insert(std::make_pair(min_distance[v_id].first,v_id));

					tempPrevDirCross[v_id].first = tri2tri;
					tempPrevDirCross[v_id].second = closestToPrevCross;
				}

			}

		}
	}





	if(finalPaths[contSource][totalNumberOfSlots].size()>=1){

		unsigned int firstPointsNo = pointPaths[contSource][totalNumberOfSlots].first.size();
		finalCenterPoints[contSource][totalNumberOfSlots].clear();
		finalCenterPoints[contSource][totalNumberOfSlots].resize(firstPointsNo + finalPaths[contSource][totalNumberOfSlots].size());
		for(unsigned int j=0; j<firstPointsNo;j++){
			finalCenterPoints[contSource][totalNumberOfSlots][j] = pointPaths[contSource][totalNumberOfSlots].first[j];
		}
 		for(unsigned int j=0; j<finalPaths[contSource][totalNumberOfSlots].size(); j++){
			finalCenterPoints[contSource][totalNumberOfSlots][j+firstPointsNo] = triangle_centers[finalPaths[contSource][totalNumberOfSlots][j]];
		}

		finalCenterPoints[contSource][totalNumberOfSlots].push_back(endPoint);
		std::string file_name_bdry = "ShortestPathsBdry_"+std::to_string(contSource)+".vtk";
		writeTestPoints(finalCenterPoints[contSource][totalNumberOfSlots], file_name_bdry);

	}
	else{
  		if(pointPaths[contSource][totalNumberOfSlots].first.size()>0){
		  	finalCenterPoints[contSource][totalNumberOfSlots].clear();
			std::string file_name_bdry = "ShortestPathsBdry_"+std::to_string(contSource)+".vtk";
			finalCenterPoints[contSource][totalNumberOfSlots].insert(finalCenterPoints[contSource][totalNumberOfSlots].end(), pointPaths[contSource][totalNumberOfSlots].first.begin(), pointPaths[contSource][totalNumberOfSlots].first.end());
			writeTestPoints(finalCenterPoints[contSource][totalNumberOfSlots], file_name_bdry);

		}

	}
	if(pointPaths[contSource][totalNumberOfSlots].first.size()>0){
		isTargetBdry = true;
		unsigned int contMatrix = contSource*totalNumberOfVariables + totalNumberOfSlots;
		contSourceToSingularity[contMatrix].first = i_orig;
		contSourceToSingularity[contMatrix].second = i_orig;

		getTraversedTrisFaceNeighbyVerts(line_discretization[i_orig][j_orig],
										line_triangles[i_orig][j_orig],
									 	finalPaths[contSource][totalNumberOfSlots],
									 	isTargetBdry,
									 	line_discretization[i_orig][j_orig],
								 		line_triangles[i_orig][j_orig],
										traversedSPTris[contSource][totalNumberOfSlots],
										isTraversedFaceCompVector,
									 	contMatrix,
									 	totalNumberOfVariables,
										IllegalCross,
									 	contSourceToSingularity,
									 	finalLastVisTri,
									 	finalEndPoint);

		finalLastVisTri = newFacesNumber;
		if(visualizeTraversedTriangles){
			std::string file_name_trav_tri = "trav_tri_BDRY_"+std::to_string(contSource)+".vtk";
			writeTestMeshTrianglesIds(traversedSPTris[contSource][totalNumberOfSlots], file_name_trav_tri);
		}
	}

 	cout<<"getShortestPathBtwFacesOptimizedNewLocalMesh end"<<endl;
	//WARNING TODO singOrGeomFaces - exclude targets!!!!!!!!!!!!!!!!!!!!!!!!!!1

}*/

/*-----------------------------------------------------------------------------------------------------*/


void SingularityGraphBuilder2D::getIllegalCrossByTraversedTris(vector<vector<vector<gmds::math::Point>>>&                    line_discretization,
                                                               vector<vector<vector<gmds::TCellID>>>&                        line_triangles,
                                                               vector<gmds::TCellID>&                                        traversedTriangles,
                                                               vector<vector<vector<gmds::TCellID>>>&                        finalPaths,
                                                               vector<pair<vector<unsigned int>, vector<unsigned int >>>&    isTraversedNodeCompVector,
                                                               vector<pair<vector<unsigned int>, vector<unsigned int >>>&    isTraversedFaceCompVector,
                                                               vector<pair<vector<unsigned int>, vector<unsigned int >>>&    isTraversedFaceCompVector_SegmentPathCode,
                                                               vector<pair<vector<unsigned int>, vector<unsigned int >>>&    isTraversedFaceCompVector_SegmentPathCont,
                                                               unsigned int&                                                 contMatrix,
                                                               unsigned int&                                                 totalNumberOfVariables, //totalNumberOfSlots+1
                                                               unsigned int&                                                 contSource,
                                                               unsigned int&                                                 contTarget,
                                                               vector<pair<unsigned int, unsigned int>>&                     IllegalCross,
                                                               vector<pair<unsigned int, unsigned int>>&                     contMatrixToSingularity,
                                                               gmds::TCellID&                                                lastVisTri,
                                                               vector<gmds::math::Point>&                                    finalEndPoint,
                                                               unsigned int&                                                 illegalProblematicBdryPair,
                                                               bool&                                                         recomputeIllegalProblematicBdryPair,
                                                               bool&                                                         reAddIt,
                                                               unsigned int&                                                 i_orig,
                                                               unsigned int&                                                 j_orig)
 {

	/* after each detection of ShortestPaths, this function detects the entire set of triangles traveresed and as well it detects with which component vector of the cross vector this path is allligned
	and it stores it;
	isTraversedFaceCompVector[face_id].first = path_id => path_id passes through face_id alligned with the first direction (in 2D case this one is horizontal); .second => alligned with the other direction*/
	if(withGlobalComments){
		cout<<"getIllegalCrossByTraversedTris contMatrix "<<contMatrix<<endl;
		cout<<" -------------------------------------------- "<<endl;
	}
	gmds::math::Point endPoint = finalEndPoint[contSource];
	vector<gmds::math::Point> discrPtsSource(line_discretization[i_orig][j_orig].begin(), line_discretization[i_orig][j_orig].end());
	vector<gmds::TCellID> travTriSource(line_triangles[i_orig][j_orig].begin(), line_triangles[i_orig][j_orig].end());
	bool withComments = false;
	//if(contMatrix==-1)
	//	withComments = true;

	if(withComments){
		cout<<"travTriSource "<<endl;
		for(unsigned int i=0; i<travTriSource.size();i++)
			cout<<"travTriSource["<<i<<"] "<<travTriSource[i]<<endl;

		cout<<"discrPtsSource "<<endl;
		for(unsigned int i=0; i<discrPtsSource.size();i++)
			cout<<"discrPtsSource["<<i<<"] "<<discrPtsSource[i]<<endl;
	}
	/*
	 if(illegalProblematicBdryPair==tempCS){*

	  recomputeIllegalProblematicBdryPair = true;
	}
	if((!recomputeIllegalProblematicBdryPair)||((recomputeIllegalProblematicBdryPair)&&(reAddIt))){
  */

	gmds::TCellID lastSlotFaceSource = travTriSource.back();

	gmds::math::Point currentPoint = discrPtsSource[0];

	gmds::math::Point nextPoint;
	gmds::TCellID nextFace, prevFace = travTriSource[max(0,(int)travTriSource.size()-2)];

	gmds::math::Point intersectionPoint;
	double intersectionParam;

	vector<bool> visitedFaces(original_faces_number, false);
	vector<bool> visitedNodes(original_nodes_number, false);
	gmds::TCellID currentFace_id = travTriSource[0];
	gmds::Face currentFace;

	if(discrPtsSource.size()>2){
		if(withComments){
			cout<<"if(discrPtsSource.size()>2){"<<endl;
		}
		nextFace = currentFace_id;
		for(unsigned int i=1; i<discrPtsSource.size(); i++){
			cout<<"i= "<<i<<endl;
			math::Segment from_seg(discrPtsSource[i-1], discrPtsSource[i]);
			math::Vector3d tempDir(discrPtsSource[i-1], discrPtsSource[i]);
			for(unsigned int j=1; j<travTriSource.size(); j++){
				currentFace_id = nextFace;
				nextFace = travTriSource[j];

				bool foundNextFace = false;

				if(withComments){
					cout<<"currentFace_id "<<currentFace_id<<endl;
					cout<<"nextFace "<<nextFace<<endl;
					cout<<"prevFace "<<prevFace<<endl;
					cout<<"from_seg "<<discrPtsSource[i-1]<<" -> "<<discrPtsSource[i]<<endl;
				}
				vector<gmds::Edge> adj_edges;
				vector<gmds::Node> temp_nodes = (m_mesh->get<gmds::Face>(currentFace_id)).get<gmds::Node>();
				for(unsigned int j1=0; j1<3; j1++){
					vector<gmds::Edge> temp_edges = temp_nodes[j1].get<gmds::Edge>();
					adj_edges.insert(adj_edges.end(), temp_edges.begin(), temp_edges.end());
				}
				for(unsigned int j1=0; j1<adj_edges.size(); j1++){
					vector<gmds::Node> currentNodes =  adj_edges[j1].get<gmds::Node>();
					math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
					if(withComments){
						cout<<"edge number j = "<<j<<" between nodes "<<currentNodes[0].id()<<" & "<<currentNodes[1]<<endl;
						from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon);

						cout<<"from_seg.SecondMetIntersect2D "<<from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon)
						<<" ; intersectionParam "<<intersectionParam<<endl;
					}
					if(from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon)){

						//WARNING if we want with nodes, here should be inserted

						vector<gmds::Face> edgeAdjFaces = adj_edges[j1].get<gmds::Face>();
						for(unsigned int k=0; k<edgeAdjFaces.size(); k++){
							if(!visitedFaces[edgeAdjFaces[k].id()]){
								currentFace_id = edgeAdjFaces[k].id();
								if(currentFace_id==nextFace)
									foundNextFace = true;

								math::Cross2D crossV = triangle_centers_cross[currentFace_id];
								unsigned int compVId = crossV.closestComponentVectorAsIndex(tempDir)%2;

								if(compVId==0){
									for(unsigned int j=0; j<isTraversedFaceCompVector[currentFace_id].first.size(); j++){
										unsigned int tempCS = isTraversedFaceCompVector[currentFace_id].first[j];
										if(tempCS!=contMatrix){
											if(illegalProblematicBdryPair==tempCS){
												recomputeIllegalProblematicBdryPair = true;
											}
											if((!recomputeIllegalProblematicBdryPair)||((recomputeIllegalProblematicBdryPair)&&(reAddIt))){
												if((contMatrixToSingularity[contMatrix].first!=contMatrixToSingularity[tempCS].first)&&
												((contMatrixToSingularity[contMatrix].second!=contMatrixToSingularity[tempCS].second)||
												(fmod(contMatrix, totalNumberOfVariables)==totalNumberOfVariables-1))){
													if((( (int)(contMatrix/totalNumberOfVariables)!=((int)(tempCS/totalNumberOfVariables)))&&
													( ( fmod(contMatrix, totalNumberOfVariables))!=(fmod(tempCS, totalNumberOfVariables)))&&
													( (int)(contMatrix/totalNumberOfVariables)!=(fmod(tempCS, totalNumberOfVariables)))&&
													( (fmod(contMatrix, totalNumberOfVariables))!=((int)(tempCS/totalNumberOfVariables)) ) )||
													(( ( fmod(contMatrix, totalNumberOfVariables))==(fmod(tempCS, totalNumberOfVariables)))&&(( fmod(contMatrix, totalNumberOfVariables))==totalNumberOfVariables-1))){
														bool toAdd = true;
														for(unsigned int k=0; k<IllegalCross.size(); k++){
															if(((contMatrix==IllegalCross[k].second)&&(tempCS==IllegalCross[k].first))||
															((tempCS==IllegalCross[k].second)&&(contMatrix==IllegalCross[k].first))){
																toAdd = false;
																break;
															}
														}
														if(toAdd){
															math::Point tempPt1, tempPt2;
															unsigned int pairContSource = (int)(tempCS/totalNumberOfVariables);
															unsigned int pairContTarget = fmod(tempCS,totalNumberOfVariables);
															unsigned int tempS_i = (int)(pairContSource/5);// source sing point
															unsigned int tempS_j = fmod(pairContSource,5);// source slot no
															unsigned int tempT_i = (int)(pairContTarget/5);// target sing point
															unsigned int tempT_j = fmod(pairContTarget,5);// target slot no
															if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==0){
																tempPt1 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]];
																tempPt2 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1];
															}
															else{
																if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==1){
																	tempPt1 = line_discretization[tempS_i][tempS_j].back();
																	if(finalPaths[pairContSource][pairContTarget].size()!=0)
																		tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][0]];
																	else
																		tempPt2 = line_discretization[tempT_i][tempT_j].back();
																}
																else{
																	if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==2){
																		tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]]];
																		tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1]];
																	}
																	else{
																		if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==3){
																			if(finalPaths[pairContSource][pairContTarget].size()!=0)
																				tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																			else
																				tempPt1 = line_discretization[tempS_i][tempS_j].back();
																			tempPt2 = line_discretization[tempT_i][tempT_j].back();// !!!it is the inverse
																		}
																		else{
																			if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==4){
																				tempPt1 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]];
																				tempPt2 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1];// !!!it is the inverse
																			}
																			else{
																				if(finalPaths[pairContSource][pairContTarget].size()!=0)
																					tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																				else
																					tempPt1 = line_discretization[tempS_i][tempS_j].back();
																				tempPt2 = finalEndPoint[pairContSource];
																			}
																		}
																	}
																}
															}

															math::Segment segment1(tempPt1, tempPt2);
															if(from_seg.SecondMetIntersect2D(segment1, intersectionPoint, intersectionParam, temp_epsilon))
																IllegalCross.push_back(make_pair(tempCS, contMatrix));
														}
													}
												}
											}
										}
									}
									isTraversedFaceCompVector[currentFace_id].first.push_back(contMatrix);
									isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first.push_back(0);
									isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first.push_back(i-1);
								}
								else{
									for(unsigned int j=0; j<isTraversedFaceCompVector[currentFace_id].second.size(); j++){
										unsigned int tempCS = isTraversedFaceCompVector[currentFace_id].second[j];
										if(tempCS!=contMatrix){
											if(illegalProblematicBdryPair==tempCS){
												recomputeIllegalProblematicBdryPair = true;
											}
											if((!recomputeIllegalProblematicBdryPair)||((recomputeIllegalProblematicBdryPair)&&(reAddIt))){
												if((contMatrixToSingularity[contMatrix].first!=contMatrixToSingularity[tempCS].first)&&
												((contMatrixToSingularity[contMatrix].second!=contMatrixToSingularity[tempCS].second)||
												(fmod(contMatrix, totalNumberOfVariables)==totalNumberOfVariables-1))){
													if((( (int)(contMatrix/totalNumberOfVariables)!=((int)(tempCS/totalNumberOfVariables)))&&
													( ( fmod(contMatrix, totalNumberOfVariables))!=(fmod(tempCS, totalNumberOfVariables)))&&
													( (int)(contMatrix/totalNumberOfVariables)!=(fmod(tempCS, totalNumberOfVariables)))&&
													( (fmod(contMatrix, totalNumberOfVariables))!=((int)(tempCS/totalNumberOfVariables)) ) )||
													(( ( fmod(contMatrix, totalNumberOfVariables))==(fmod(tempCS, totalNumberOfVariables)))&&(( fmod(contMatrix, totalNumberOfVariables))==totalNumberOfVariables-1))){
														bool toAdd = true;
														for(unsigned int k=0; k<IllegalCross.size(); k++){
															if(((contMatrix==IllegalCross[k].second)&&(tempCS==IllegalCross[k].first))||
															((tempCS==IllegalCross[k].second)&&(contMatrix==IllegalCross[k].first))){
																toAdd = false;
																break;
															}
														}
														if(toAdd){
															math::Point tempPt1, tempPt2;
															unsigned int pairContSource = (int)(tempCS/totalNumberOfVariables);
															unsigned int pairContTarget = fmod(tempCS,totalNumberOfVariables);
															unsigned int tempS_i = (int)(pairContSource/5);// source sing point
															unsigned int tempS_j = fmod(pairContSource,5);// source slot no
															unsigned int tempT_i = (int)(pairContTarget/5);// target sing point
															unsigned int tempT_j = fmod(pairContTarget,5);// target slot no
															if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==0){
																tempPt1 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]];
																tempPt2 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1];
															}
															else{
																if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==1){
																	tempPt1 = line_discretization[tempS_i][tempS_j].back();
																	if(finalPaths[pairContSource][pairContTarget].size()!=0)
																		tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][0]];
																	else
																		tempPt2 = line_discretization[tempT_i][tempT_j].back();
																}
																else{
																	if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==2){
																		tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]]];
																		tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1]];
																	}
																	else{
																		if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==3){
																			if(finalPaths[pairContSource][pairContTarget].size()!=0)
																				tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																			else
																				tempPt1 = line_discretization[tempS_i][tempS_j].back();
																			tempPt2 = line_discretization[tempT_i][tempT_j].back();// !!!it is the inverse
																		}
																		else{
																			if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==4){
																				tempPt1 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]];
																				tempPt2 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1];// !!!it is the inverse
																			}
																			else{
																				if(finalPaths[pairContSource][pairContTarget].size()!=0)
																					tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																				else
																					tempPt1 = line_discretization[tempS_i][tempS_j].back();
																				tempPt2 = finalEndPoint[pairContSource];
																			}
																		}
																	}
																}
															}

															math::Segment segment1(tempPt1, tempPt2);
															if(from_seg.SecondMetIntersect2D(segment1, intersectionPoint, intersectionParam, temp_epsilon))
																IllegalCross.push_back(make_pair(tempCS, contMatrix));

														}
													}
												}
											}
										}
									}
									isTraversedFaceCompVector[currentFace_id].second.push_back(contMatrix);
									isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second.push_back(0);
									isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second.push_back(i-1);
								}

								visitedFaces[edgeAdjFaces[k].id()] = true;
								traversedTriangles.push_back(edgeAdjFaces[k].id());
							}
						}
					}
				}
				if((!foundNextFace)&&(withComments))
					cout<<"!foundNextFace0"<<endl;
			}
		}
	}

	nextPoint = discrPtsSource.back();
	if(withComments){
		cout<<"currentFace_id "<<currentFace_id<<endl;
		cout<<"currentPoint "<<currentPoint<<endl;
		cout<<"traversedTriangles.size() "<<traversedTriangles.size()<<endl;
		cout<<"prevFace "<<prevFace<<endl;
	}
	traversedTriangles.push_back(currentFace_id);
	visitedFaces[currentFace_id] = true;

	currentPoint = nextPoint;
	if(lastVisTri != original_faces_number){
		currentFace_id = nextFace;
		nextFace = lastVisTri;

		nextPoint = finalEndPoint[contSource];
		math::Segment from_seg(currentPoint, nextPoint);

		math::Vector3d tempDir(currentPoint, nextPoint);
		if(withComments){
			cout<<"if(lastVisTri != original_faces_number){ "<<endl;
			cout<<"currentFace_id "<<currentFace_id<<endl;
			cout<<"nextFace "<<nextFace<<endl;
			cout<<"prevFace "<<prevFace<<endl;
			cout<<"finalEndPoint[contSource] "<<finalEndPoint[contSource]<<endl;
		}
		bool foundNextFace = false;

		if(withComments){
			cout<<"currentFace_id "<<currentFace_id<<endl;
			cout<<"nextFace "<<nextFace<<endl;
			cout<<"prevFace "<<prevFace<<endl;
		}
		vector<gmds::Edge> adj_edges;
		vector<gmds::Node> temp_nodes = (m_mesh->get<gmds::Face>(currentFace_id)).get<gmds::Node>();
		for(unsigned int j1=0; j1<3; j1++){
			vector<gmds::Edge> temp_edges = temp_nodes[j1].get<gmds::Edge>();
			adj_edges.insert(adj_edges.end(), temp_edges.begin(), temp_edges.end());
		}
		for(unsigned int j1=0; j1<adj_edges.size(); j1++){
			vector<gmds::Node> currentNodes =  adj_edges[j1].get<gmds::Node>();
			math::Segment oppSeg(currentNodes[0].getPoint(), currentNodes[1].getPoint());
			if(withComments){
				cout<<"edge number j = "<<j1<<" between nodes "<<currentNodes[0].id()<<" & "<<currentNodes[1]<<endl;
				from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon);

				cout<<"from_seg.SecondMetIntersect2D "<<from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon)
				<<endl;
				cout<<"temp_epsilon "<<temp_epsilon<<endl;
				cout<<"intersectionParam "<<intersectionParam<<endl;
			}
			if(from_seg.SecondMetIntersect2D(oppSeg, intersectionPoint, intersectionParam, temp_epsilon)){

				//WARNING if we want with nodes, here should be inserted

				vector<gmds::Face> edgeAdjFaces = adj_edges[j1].get<gmds::Face>();
				for(unsigned int k=0; k<edgeAdjFaces.size(); k++){
					if(!visitedFaces[edgeAdjFaces[k].id()]){
						currentFace_id = edgeAdjFaces[k].id();
						if(currentFace_id==nextFace)
							foundNextFace = true;

						math::Cross2D crossV = triangle_centers_cross[currentFace_id];
						unsigned int compVId = crossV.closestComponentVectorAsIndex(tempDir)%2;

						if(compVId==0){
							for(unsigned int j=0; j<isTraversedFaceCompVector[currentFace_id].first.size(); j++){
								unsigned int tempCS = isTraversedFaceCompVector[currentFace_id].first[j];
								if(tempCS!=contMatrix){
									if(illegalProblematicBdryPair==tempCS)
										recomputeIllegalProblematicBdryPair = true;

									if((!recomputeIllegalProblematicBdryPair)||((recomputeIllegalProblematicBdryPair)&&(reAddIt))){
										if((contMatrixToSingularity[contMatrix].first!=contMatrixToSingularity[tempCS].first)&&
										((contMatrixToSingularity[contMatrix].second!=contMatrixToSingularity[tempCS].second)||
										(fmod(contMatrix, totalNumberOfVariables)==totalNumberOfVariables-1))){
											if((( (int)(contMatrix/totalNumberOfVariables)!=((int)(tempCS/totalNumberOfVariables)))&&
											( ( fmod(contMatrix, totalNumberOfVariables))!=(fmod(tempCS, totalNumberOfVariables)))&&
											( (int)(contMatrix/totalNumberOfVariables)!=(fmod(tempCS, totalNumberOfVariables)))&&
											( (fmod(contMatrix, totalNumberOfVariables))!=((int)(tempCS/totalNumberOfVariables)) ) )||
											(( ( fmod(contMatrix, totalNumberOfVariables))==(fmod(tempCS, totalNumberOfVariables)))&&(( fmod(contMatrix, totalNumberOfVariables))==totalNumberOfVariables-1))){
												bool toAdd = true;
												for(unsigned int k=0; k<IllegalCross.size(); k++){
													if(((contMatrix==IllegalCross[k].second)&&(tempCS==IllegalCross[k].first))||
													((tempCS==IllegalCross[k].second)&&(contMatrix==IllegalCross[k].first))){
														toAdd = false;
														break;
													}
												}
												if(toAdd){
													math::Point tempPt1, tempPt2;
													unsigned int pairContSource = (int)(tempCS/totalNumberOfVariables);
													unsigned int pairContTarget = fmod(tempCS,totalNumberOfVariables);
													unsigned int tempS_i = (int)(pairContSource/5);// source sing point
													unsigned int tempS_j = fmod(pairContSource,5);// source slot no
													unsigned int tempT_i = (int)(pairContTarget/5);// target sing point
													unsigned int tempT_j = fmod(pairContTarget,5);// target slot no
													if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==0){
														tempPt1 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]];
														tempPt2 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1];
													}
													else{
														if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==1){
															tempPt1 = line_discretization[tempS_i][tempS_j].back();
															if(finalPaths[pairContSource][pairContTarget].size()!=0)
																tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][0]];
															else
																tempPt2 = line_discretization[tempT_i][tempT_j].back();
														}
														else{
															if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==2){
																tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]]];
																tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1]];
															}
															else{
																if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==3){
																	if(finalPaths[pairContSource][pairContTarget].size()!=0)
																		tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																	else
																		tempPt1 = line_discretization[tempS_i][tempS_j].back();
																	tempPt2 = line_discretization[tempT_i][tempT_j].back();// !!!it is the inverse
																}
																else{
																	if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first[j]==4){
																		tempPt1 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]];
																		tempPt2 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first[j]+1];// !!!it is the inverse
																	}
																	else{
																		if(finalPaths[pairContSource][pairContTarget].size()!=0)
																			tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																		else
																			tempPt1 = line_discretization[tempS_i][tempS_j].back();
																		tempPt2 = finalEndPoint[pairContSource];
																	}
																}
															}
														}
													}

													math::Segment segment1(tempPt1, tempPt2);
													if(from_seg.SecondMetIntersect2D(segment1, intersectionPoint, intersectionParam, temp_epsilon))
														IllegalCross.push_back(make_pair(tempCS, contMatrix));
												}
											}
										}
									}
								}
							}
							isTraversedFaceCompVector[currentFace_id].first.push_back(contMatrix);
							isTraversedFaceCompVector_SegmentPathCode[currentFace_id].first.push_back(5);
							isTraversedFaceCompVector_SegmentPathCont[currentFace_id].first.push_back(0);
						}
						else{
							for(unsigned int j=0; j<isTraversedFaceCompVector[currentFace_id].second.size(); j++){
								unsigned int tempCS = isTraversedFaceCompVector[currentFace_id].second[j];
								if(tempCS!=contMatrix){
									if(illegalProblematicBdryPair==tempCS)
										recomputeIllegalProblematicBdryPair = true;

									if((!recomputeIllegalProblematicBdryPair)||((recomputeIllegalProblematicBdryPair)&&(reAddIt))){

										if((contMatrixToSingularity[contMatrix].first!=contMatrixToSingularity[tempCS].first)&&
										((contMatrixToSingularity[contMatrix].second!=contMatrixToSingularity[tempCS].second)||
										(fmod(contMatrix, totalNumberOfVariables)==totalNumberOfVariables-1))){
											if((( (int)(contMatrix/totalNumberOfVariables)!=((int)(tempCS/totalNumberOfVariables)))&&
											( ( fmod(contMatrix, totalNumberOfVariables))!=(fmod(tempCS, totalNumberOfVariables)))&&
											( (int)(contMatrix/totalNumberOfVariables)!=(fmod(tempCS, totalNumberOfVariables)))&&
											( (fmod(contMatrix, totalNumberOfVariables))!=((int)(tempCS/totalNumberOfVariables)) ) )||
											(( ( fmod(contMatrix, totalNumberOfVariables))==(fmod(tempCS, totalNumberOfVariables)))&&(( fmod(contMatrix, totalNumberOfVariables))==totalNumberOfVariables-1))){
												bool toAdd = true;
												for(unsigned int k=0; k<IllegalCross.size(); k++){
													if(((contMatrix==IllegalCross[k].second)&&(tempCS==IllegalCross[k].first))||
													((tempCS==IllegalCross[k].second)&&(contMatrix==IllegalCross[k].first))){
														toAdd = false;
														break;
													}
												}
												if(toAdd){
													math::Point tempPt1, tempPt2;
													unsigned int pairContSource = (int)(tempCS/totalNumberOfVariables);
													unsigned int pairContTarget = fmod(tempCS,totalNumberOfVariables);
													unsigned int tempS_i = (int)(pairContSource/5);// source sing point
													unsigned int tempS_j = fmod(pairContSource,5);// source slot no
													unsigned int tempT_i = (int)(pairContTarget/5);// target sing point
													unsigned int tempT_j = fmod(pairContTarget,5);// target slot no
													if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==0){
														tempPt1 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]];
														tempPt2 = line_discretization[tempS_i][tempS_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1];
													}
													else{
														if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==1){
															tempPt1 = line_discretization[tempS_i][tempS_j].back();
															if(finalPaths[pairContSource][pairContTarget].size()!=0)
																tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][0]];
															else
																tempPt2 = line_discretization[tempT_i][tempT_j].back();
														}
														else{
															if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==2){
																tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]]];
																tempPt2 = triangle_centers[finalPaths[pairContSource][pairContTarget][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1]];
															}
															else{
																if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==3){
																	if(finalPaths[pairContSource][pairContTarget].size()!=0)
																		tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																	else
																		tempPt1 = line_discretization[tempS_i][tempS_j].back();
																	tempPt2 = line_discretization[tempT_i][tempT_j].back();// !!!it is the inverse
																}
																else{
																	if(isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second[j]==4){
																		tempPt1 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]];
																		tempPt2 = line_discretization[tempT_i][tempT_j][isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second[j]+1];// !!!it is the inverse
																	}
																	else{
																		if(finalPaths[pairContSource][pairContTarget].size()!=0)
																			tempPt1 = triangle_centers[finalPaths[pairContSource][pairContTarget].back()];
																		else
																			tempPt1 = line_discretization[tempS_i][tempS_j].back();
																		tempPt2 = finalEndPoint[pairContSource];
																	}
																}
															}
														}
													}

													math::Segment segment1(tempPt1, tempPt2);
													if(from_seg.SecondMetIntersect2D(segment1, intersectionPoint, intersectionParam, temp_epsilon))
														IllegalCross.push_back(make_pair(tempCS, contMatrix));
												}
											}
										}
									}
								}
							}
							isTraversedFaceCompVector[currentFace_id].second.push_back(contMatrix);
							isTraversedFaceCompVector_SegmentPathCode[currentFace_id].second.push_back(5);
							isTraversedFaceCompVector_SegmentPathCont[currentFace_id].second.push_back(0);
						}

						visitedFaces[edgeAdjFaces[k].id()] = true;
						traversedTriangles.push_back(edgeAdjFaces[k].id());
					}
				}
			}
		}
		if((!foundNextFace)&&(withComments))
			cout<<"!foundNextFace4"<<endl;
	}


	if(withGlobalComments)
		cout<<"getIllegalCrossByTraversedTris end"<<endl;
}
