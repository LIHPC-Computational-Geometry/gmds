//
// Created by bourmaudp on 20/12/22.
//
#include <gmds/blocking_Quality/mainBis1.h>
using namespace gmds;

void test(){
	std::cout<<"TEST"<<std::endl;
}

void fonctionMain1(){
	std::cout<<"TEST"<<std::endl;

	/*
	Mesh mImprint (	MeshModel(DIM3|F|E|N|F2N|N2F|E2N));
	IGMeshIOService ioServiceRead(&mImprint);
	VTKReader vtkReader(&ioServiceRead);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	//vtkReader.read("/home/bourmaudp/Documents/GMDS/gmds/test_samples/HolesInSquare0.vtk");
	vtkReader.read("/home/bourmaudp/Documents/GMDS/gmds/test_samples/Cube.vtk");
	 */

	/*=========================================================*/
	Mesh simpleCube(MeshModel(DIM3|F|E|N|F2N|N2F|E2N|N2E|E2F|F2E|E2R));


	Variable<int> *boundaryEdge = simpleCube.newVariable<int,gmds::GMDS_EDGE>("boundaryEdge");

	Node n0 = simpleCube.newNode(math::Point(0,0,0));
	Node n1 = simpleCube.newNode(math::Point(1,0,0));
	Node n2 = simpleCube.newNode(math::Point(1,0,-1));
	Node n3 = simpleCube.newNode(math::Point(0,0,-1));
	Node n4 = simpleCube.newNode(math::Point(0,1,0));
	Node n5 = simpleCube.newNode(math::Point(1,1,0));
	Node n6 = simpleCube.newNode(math::Point(1,1,-1));
	Node n7 = simpleCube.newNode(math::Point(0,1,-1));


	Edge e0 = simpleCube.newEdge(n0,n1);
	simpleCube.get<Node>(n0.id()).add(e0);
	simpleCube.get<Node>(n1.id()).add(e0);
	Edge e1 = simpleCube.newEdge(n1,n2);
	simpleCube.get<Node>(n1.id()).add(e1);
	simpleCube.get<Node>(n2.id()).add(e1);
	Edge e2 = simpleCube.newEdge(n2,n3);
	simpleCube.get<Node>(n2.id()).add(e2);
	simpleCube.get<Node>(n3.id()).add(e2);
	Edge e3 =simpleCube.newEdge(n3,n0);
	simpleCube.get<Node>(n3.id()).add(e3);
	simpleCube.get<Node>(n0.id()).add(e3);

	Edge e4 = simpleCube.newEdge(n4,n5);
	simpleCube.get<Node>(n4.id()).add(e4);
	simpleCube.get<Node>(n5.id()).add(e4);
	Edge e5 = simpleCube.newEdge(n5,n6);
	simpleCube.get<Node>(n5.id()).add(e5);
	simpleCube.get<Node>(n6.id()).add(e5);
	Edge e6 = simpleCube.newEdge(n6,n7);
	simpleCube.get<Node>(n6.id()).add(e6);
	simpleCube.get<Node>(n7.id()).add(e6);
	Edge e7 = simpleCube.newEdge(n7,n4);
	simpleCube.get<Node>(n7.id()).add(e7);
	simpleCube.get<Node>(n4.id()).add(e7);

	Edge e8 = simpleCube.newEdge(n0,n4);
	simpleCube.get<Node>(n0.id()).add(e8);
	simpleCube.get<Node>(n4.id()).add(e8);
	Edge e9 = simpleCube.newEdge(n1,n5);
	simpleCube.get<Node>(n1.id()).add(e9);
	simpleCube.get<Node>(n5.id()).add(e9);
	Edge e10 = simpleCube.newEdge(n2,n6);
	simpleCube.get<Node>(n2.id()).add(e10);
	simpleCube.get<Node>(n6.id()).add(e10);
	Edge e11 = simpleCube.newEdge(n3,n7);
	simpleCube.get<Node>(n3.id()).add(e11);
	simpleCube.get<Node>(n7.id()).add(e11);


	//Set value boundaryEdge
	boundaryEdge->set(e0.id(),1);
	boundaryEdge->set(e1.id(),1);
	boundaryEdge->set(e2.id(),1);
	boundaryEdge->set(e3.id(),1);
	boundaryEdge->set(e4.id(),1);
	boundaryEdge->set(e5.id(),1);
	boundaryEdge->set(e6.id(),1);
	boundaryEdge->set(e7.id(),1);
	boundaryEdge->set(e8.id(),1);
	boundaryEdge->set(e9.id(),1);
	boundaryEdge->set(e10.id(),1);
	boundaryEdge->set(e11.id(),1);


	Face f0 = simpleCube.newQuad(n0,n1,n2,n3);
	simpleCube.get<Edge>(e0.id()).add(f0);
	simpleCube.get<Edge>(e1.id()).add(f0);
	simpleCube.get<Edge>(e2.id()).add(f0);
	simpleCube.get<Edge>(e3.id()).add(f0);
	Face f1 = simpleCube.newQuad(n4,n5,n6,n7);
	simpleCube.get<Edge>(e4.id()).add(f1);
	simpleCube.get<Edge>(e5.id()).add(f1);
	simpleCube.get<Edge>(e6.id()).add(f1);
	simpleCube.get<Edge>(e7.id()).add(f1);
	Face f2 = simpleCube.newQuad(n1,n2,n6,n5);
	simpleCube.get<Edge>(e1.id()).add(f2);
	simpleCube.get<Edge>(e5.id()).add(f2);
	simpleCube.get<Edge>(e10.id()).add(f2);
	simpleCube.get<Edge>(e9.id()).add(f2);
	Face f3 = simpleCube.newQuad(n0,n3,n7,n4);
	simpleCube.get<Edge>(e2.id()).add(f3);
	simpleCube.get<Edge>(e6.id()).add(f3);
	simpleCube.get<Edge>(e10.id()).add(f3);
	simpleCube.get<Edge>(e11.id()).add(f3);
	Face f4 = simpleCube.newQuad(n0,n1,n5,n4);
	simpleCube.get<Edge>(e3.id()).add(f4);
	simpleCube.get<Edge>(e7.id()).add(f4);
	simpleCube.get<Edge>(e11.id()).add(f4);
	simpleCube.get<Edge>(e8.id()).add(f4);
	Face f5 = simpleCube.newQuad(n3,n2,n6,n7);
	simpleCube.get<Edge>(e0.id()).add(f5);
	simpleCube.get<Edge>(e4.id()).add(f5);
	simpleCube.get<Edge>(e8.id()).add(f5);
	simpleCube.get<Edge>(e9.id()).add(f5);


	Region r1 = simpleCube.newHex(n0,n1,n2,n3,n4,n5,n6,n7);


	simpleCube.get<Edge>(e0.id()).add(r1);
	simpleCube.get<Edge>(e1.id()).add(r1);
	simpleCube.get<Edge>(e2.id()).add(r1);
	simpleCube.get<Edge>(e3.id()).add(r1);
	simpleCube.get<Edge>(e4.id()).add(r1);
	simpleCube.get<Edge>(e5.id()).add(r1);
	simpleCube.get<Edge>(e6.id()).add(r1);
	simpleCube.get<Edge>(e7.id()).add(r1);
	simpleCube.get<Edge>(e8.id()).add(r1);
	simpleCube.get<Edge>(e9.id()).add(r1);
	simpleCube.get<Edge>(e10.id()).add(r1);
	simpleCube.get<Edge>(e11.id()).add(r1);



	/*=========================================================*/

	Mesh deformedCube(MeshModel(DIM3|F|E|N|F2N|N2F|R2N));
	Node dn0 = deformedCube.newNode(math::Point(0,0,0));
	Node dn1 = deformedCube.newNode(math::Point(1,0,0));
	Node dn2 = deformedCube.newNode(math::Point(1,0,-1));
	Node dn3 = deformedCube.newNode(math::Point(0,0,-1));
	Node dn4 = deformedCube.newNode(math::Point(0,1,0));
	Node dn5 = deformedCube.newNode(math::Point(2,1,0));
	Node dn6 = deformedCube.newNode(math::Point(1,1,-1));
	Node dn7 = deformedCube.newNode(math::Point(0,1,-1));

	deformedCube.newQuad(dn0,dn1,dn2,dn3);
	deformedCube.newQuad(dn4,dn5,dn6,dn7);
	deformedCube.newQuad(dn1,dn2,dn6,dn5);
	deformedCube.newQuad(dn0,dn3,dn7,dn4);
	deformedCube.newQuad(dn0,dn1,dn5,dn4);
	deformedCube.newQuad(dn3,dn2,dn6,dn7);

	deformedCube.newHex(dn0,dn1,dn2,dn3,dn4,dn5,dn6,dn7);

	/*=========================================================*/

	Mesh simpleRectangle(MeshModel(DIM3|F|E|N|F2N|N2F|R2N));
	Node nr0 = simpleRectangle.newNode(math::Point(0,0,0));
	Node nr1 = simpleRectangle.newNode(math::Point(3,0,0));
	Node nr2 = simpleRectangle.newNode(math::Point(3,0,-1));
	Node nr3 = simpleRectangle.newNode(math::Point(0,0,-1));
	Node nr4 = simpleRectangle.newNode(math::Point(0,1,0));
	Node nr5 = simpleRectangle.newNode(math::Point(3,1,0));
	Node nr6 = simpleRectangle.newNode(math::Point(3,1,-1));
	Node nr7 = simpleRectangle.newNode(math::Point(0,1,-1));

	simpleRectangle.newQuad(nr0,nr1,nr2,nr3);
	simpleRectangle.newQuad(nr4,nr5,nr6,nr7);
	simpleRectangle.newQuad(nr1,nr2,nr6,nr5);
	simpleRectangle.newQuad(nr0,nr3,nr7,nr4);
	simpleRectangle.newQuad(nr0,nr1,nr5,nr4);
	simpleRectangle.newQuad(nr3,nr2,nr6,nr7);

	simpleRectangle.newHex(nr0,nr1,nr2,nr3,nr4,nr5,nr6,nr7);

	/*=========================================================*/
	// RECTANGLE SIMPLE COPY

	Mesh simpleRectangleCopy(MeshModel(DIM3|F|E|N|F2N|N2F|R2N));
	Node nr0B = simpleRectangleCopy.newNode(math::Point(0,0,0));
	Node nr1B = simpleRectangleCopy.newNode(math::Point(3,0,0));
	Node nr2B = simpleRectangleCopy.newNode(math::Point(3,0,-1));
	Node nr3B = simpleRectangleCopy.newNode(math::Point(0,0,-1));
	Node nr4B = simpleRectangleCopy.newNode(math::Point(0,1,0));
	Node nr5B = simpleRectangleCopy.newNode(math::Point(3,1,0));
	Node nr6B = simpleRectangleCopy.newNode(math::Point(3,1,-1));
	Node nr7B = simpleRectangleCopy.newNode(math::Point(0,1,-1));

	simpleRectangleCopy.newQuad(nr0B,nr1B,nr2B,nr3B);
	simpleRectangleCopy.newQuad(nr4B,nr5B,nr6B,nr7B);
	simpleRectangleCopy.newQuad(nr1B,nr2B,nr6B,nr5B);
	simpleRectangleCopy.newQuad(nr0B,nr3B,nr7B,nr4B);
	simpleRectangleCopy.newQuad(nr0B,nr1B,nr5B,nr4B);
	simpleRectangleCopy.newQuad(nr3B,nr2B,nr6B,nr7B);

	simpleRectangle.newHex(nr0B,nr1B,nr2B,nr3B,nr4B,nr5B,nr6B,nr7B);


	/*=========================================================*/
	Mesh doubleRectangle(MeshModel(DIM3|F|E|N|F2N|N2F|R2N));
	Node dnr0 = doubleRectangle.newNode(math::Point(0,0,0));
	Node dnr1 = doubleRectangle.newNode(math::Point(3,0,0));
	Node dnr2 = doubleRectangle.newNode(math::Point(3,0,-1));
	Node dnr3 = doubleRectangle.newNode(math::Point(0,0,-1));
	Node dnr4 = doubleRectangle.newNode(math::Point(0,1,0));
	Node dnr5 = doubleRectangle.newNode(math::Point(3,1,0));
	Node dnr6 = doubleRectangle.newNode(math::Point(3,1,-1));
	Node dnr7 = doubleRectangle.newNode(math::Point(0,1,-1));
	Node dnr8 = doubleRectangle.newNode(math::Point(6,0,0));
	Node dnr9 = doubleRectangle.newNode(math::Point(6,0,-1));
	Node dnr10 = doubleRectangle.newNode(math::Point(6,1,0));
	Node dnr11 = doubleRectangle.newNode(math::Point(6,1,-1));


	doubleRectangle.newQuad(dnr0,dnr1,dnr2,dnr3);
	doubleRectangle.newQuad(dnr4,dnr5,dnr6,dnr7);
	doubleRectangle.newQuad(dnr1,dnr2,dnr6,dnr5);
	doubleRectangle.newQuad(dnr0,dnr3,dnr7,dnr4);
	doubleRectangle.newQuad(dnr0,dnr1,dnr5,dnr4);
	doubleRectangle.newQuad(dnr3,dnr2,dnr6,dnr7);

	doubleRectangle.newQuad(dnr1,dnr8,dnr9,dnr2);
	doubleRectangle.newQuad(dnr1,dnr8,dnr10,dnr5);
	doubleRectangle.newQuad(dnr5,dnr10,dnr11,dnr6);
	doubleRectangle.newQuad(dnr2,dnr9,dnr11,dnr6);
	doubleRectangle.newQuad(dnr8,dnr9,dnr11,dnr10);

	doubleRectangle.newHex(dnr0,dnr1,dnr2,dnr3,dnr4,dnr5,dnr6,dnr7);
	doubleRectangle.newHex(dnr1,dnr8,dnr9,dnr2,dnr5,dnr10,dnr11,dnr6);

	/*=========================================================*/
	// CUBE SIMPLE POUR LE RUBIKSCUBE DE TAILLE 2x2x2
	Mesh cubeRubiksCube(MeshModel(DIM3|F|R|E|N|R2F|R2E|R2N|F2N|N2F|E2N|N2E|E2F|F2E|E2R));
	Node crn0 = cubeRubiksCube.newNode(math::Point(0,0,0));
	Node crn1 = cubeRubiksCube.newNode(math::Point(4,0,0));
	Node crn2 = cubeRubiksCube.newNode(math::Point(4,0,4));
	Node crn3 = cubeRubiksCube.newNode(math::Point(0,0,4));
	Node crn4 = cubeRubiksCube.newNode(math::Point(0,4,0));
	Node crn5 = cubeRubiksCube.newNode(math::Point(4,4,0));
	Node crn6 = cubeRubiksCube.newNode(math::Point(4,4,4));
	Node crn7 = cubeRubiksCube.newNode(math::Point(0,4,4));

	Region crr1 = cubeRubiksCube.newHex(crn3,crn2,crn1,crn0,crn7,crn6,crn5,crn4);

	MeshDoctor cubeDoctor(&cubeRubiksCube);
	cubeDoctor.buildFacesAndR2F();
	cubeDoctor.buildEdgesAndX2E();
	cubeDoctor.updateUpwardConnectivity();

	Variable<int> *boundaryEdgeCubeRubiksCube = cubeRubiksCube.newVariable<int,gmds::GMDS_EDGE>("boundaryEdge");

	for (auto e : cubeRubiksCube.edges()){
		//std::cout<<"in edge : "<<e<<std::endl;
		//std::cout<<"Nb regions : "<<rubiksCube.get<Edge>(e).nbRegions()<<std::endl;
		if(cubeRubiksCube.get<Edge>(e).nbRegions()== 1){
			//std::cout<<"dans if "<<std::endl;
			boundaryEdgeCubeRubiksCube->set(e,1);
		}
		else{
			//std::cout<<"Dans else"<<std::endl;
			boundaryEdgeCubeRubiksCube->set(e,0);
		}
	}


	/*=========================================================*/
	// FORME GRILLE 2x2x2

	Mesh rubiksCube(MeshModel(DIM3|F|R|E|N|R2F|R2E|R2N|F2N|N2F|E2N|N2E|E2F|F2E|E2R));
	GridBuilder gridBuilder(&rubiksCube,3);
	double step = 1;
	int nb = 5;
	int AXnb = nb;
	int AYnb = nb;
	int AZnb = nb;

	gridBuilder.execute(AXnb,step,AYnb,step,AZnb,step);
	MeshDoctor doctor(&rubiksCube);
	doctor.buildFacesAndR2F();
	doctor.buildEdgesAndX2E();
	doctor.updateUpwardConnectivity();

	Variable<int> *boundaryEdgeRubiksCube = rubiksCube.newVariable<int,gmds::GMDS_EDGE>("boundaryEdge");

	/*	for(auto r : rubiksCube.regions()){
	      Region region = rubiksCube.get<Region>(r);
	      std::cout<<"region "<< r <<std::endl;
	      std::vector<TCellID> edges;
	      region.delegateGetEdgeIDs(edges);
	      for(auto e : edges){
	         rubiksCube.get<Edge>(e).add(region);
	      }
	   }*/

	for (auto e : rubiksCube.edges()){
		//std::cout<<"in edge : "<<e<<std::endl;
		//std::cout<<"Nb regions : "<<rubiksCube.get<Edge>(e).nbRegions()<<std::endl;
		if(rubiksCube.get<Edge>(e).nbRegions()== 1){
			//std::cout<<"dans if "<<std::endl;
			boundaryEdgeRubiksCube->set(e,1);
		}
		else{
			//std::cout<<"Dans else"<<std::endl;
			boundaryEdgeRubiksCube->set(e,0);
		}
	}





	double qualityValue = 0;

	//HexQuality simpleCube
	//quality::HexQuality he = quality::HexQuality::build(n0.point(),n1.point(),n2.point(),n3.point(),n4.point(),n5.point(),n6.point(),n7.point());

	//HexQuality deformedCube
	//quality::HexQuality he = quality::HexQuality::build(dn0.point(),dn1.point(),dn2.point(),dn3.point(),dn4.point(),dn5.point(),dn6.point(),dn7.point());

	//HexQuality simpleRectangle
	//quality::HexQuality he = quality::HexQuality::build(nr0.point(),nr1.point(),nr2.point(),nr3.point(),nr4.point(),nr5.point(),nr6.point(),nr7.point());

	//HexQuality bigRectangle
	//quality::HexQuality he = quality::HexQuality::build(dnr0.point(),dnr1.point(),dnr2.point(),dnr3.point(),dnr4.point(),dnr5.point(),dnr6.point(),dnr7.point());

	//HexQuality 2 rectangles -> bigRectangle
	//quality::HexQuality he = quality::HexQuality::build(dnr0.point(),dnr1.point(),dnr2.point(),dnr3.point(),dnr4.point(),dnr5.point(),dnr6.point(),dnr7.point());
	quality::HexQuality he2 = quality::HexQuality::build(dnr1.point(),dnr8.point(),dnr9.point(),dnr2.point(),dnr5.point(),dnr10.point(),dnr11.point(),dnr6.point());





	/*=============================					TEST FONCTIONS 				==============================================*/

	std::cout<<"test Fonction :"<<std::endl;
	//blockQuality(doubleRectangle.get<Region>(0));
	//blockQuality(doubleRectangle.get<Region>(1));
	//allBlocksQuality(&doubleRectangle);
	//layeringNodes(&simpleRectangle, &simpleRectangleCopy);
	//std::cout<<"ratio global : \n "<<ratioGlobalEdgeValence(&rubiksCube)<<std::endl;
	//std::cout<<"nb Blocs "<<blockingSimplicity(&simpleCube);
	std::cout<<"Test blockingQuality : " <<blockingQuality(&rubiksCube,&cubeRubiksCube);




	/*=======================================================================================================================*/
	/*
	std::cout<<"Quality bloc 2"<<std::endl;
	std::cout<<"Le volume "<<he2.volume()<<std::endl;
	std::cout<<"La diagonal "<<he2.diagonal()<<std::endl;
	std::cout<<"Le edgeRatio "<<he2.edgeRatio()<<std::endl;
	std::cout<<"Le maxi edgeRatio "<<he2.maximumEdgeRatio()<<std::endl;
	std::cout<<"Le Jacobien "<<he2.jacobian()<<std::endl;
	std::cout<<"Le Scaled Jacobien "<<he2.scaledJacobian()<<std::endl;
	std::cout<<"Le shape "<<he2.shape()<<std::endl;
	std::cout<<"Le maxiAF "<<he2.maximumAspectFrobenius()<<std::endl;
	std::cout<<"Le meanAF "<<he2.meanAspectFrobenius()<<std::endl;
	std::cout<<"Le oddy "<<he2.oddy()<<std::endl;
	std::cout<<"Le shear "<<he2.shear()<<std::endl;
	std::cout<<"Le skew "<<he2.skew()<<std::endl;
	std::cout<<"Le stretch "<<he2.stretch()<<std::endl;
	std::cout<<"Le taper "<<he2.taper()<<std::endl;
	 */
	/*
	gmds::MeshDoctor doc_imprint(&simpleCube);
	doc_imprint.updateUpwardConnectivity();
	doc_imprint.buildBoundaryCells();

	std::cout<<"Le volume apres doctor "<<he.volume()<<std::endl;
	*/


	gmds::IGMeshIOService ioService_write(&rubiksCube);
	gmds::VTKWriter vtkWriter(&ioService_write);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("rubiksCube.vtk");

	gmds::IGMeshIOService ioService_writeBis(&cubeRubiksCube);
	gmds::VTKWriter vtkWriterBis(&ioService_writeBis);
	vtkWriterBis.setCellOptions(gmds::N|gmds::F);
	vtkWriterBis.setDataOptions(gmds::N|gmds::F);
	vtkWriterBis.write("cubeRubiksCube.vtk");
}

void fonctionMain2(){

	// TEST BLOCKMESHER

	std::cout << "============== Blocker ================" << std::endl;
	//==================================================================
	// PARAMETERS' PARSING
	//==================================================================
	std::string file_geom, file_mesh, file_mesh_out, file_block_out;
	int nb_curve_smooth_iterations=0;
	int nb_surface_smooth_iterations=0;
	int nb_volume_smooth_iterations=0;
	int edge_discretization=0;
	/*
	if (argc != 9) {
	   std::cout << "Require four paramaters : \n";
	   std::cout << "  - [IN ] tetrahedral mesh (.vtk) that describes the geometry, \n";
	   std::cout << "  - [IN ] mesh (.vtk) that describes the block structure to refine and smoothed, \n";
	   std::cout << "  - [IN ] edge discretization of each block edge,\n";
	   std::cout << "  - [IN ] number of curve smoothing iterations,\n";
	   std::cout << "  - [IN ] number of surface smoothing iterations,\n";
	   std::cout << "  - [IN ] number of volume smoothing iterations,\n";
	   std::cout << "  - [OUT] the final blocks (.vtk),\n";
	   std::cout << "  - [OUT] the final mesh (.vtk). \n"<< std::endl;
	   throw gmds::GMDSException("Wrong number of parameters");
	}

	file_geom = std::string(argv[1]);
	file_mesh = std::string(argv[2]);
	nb_curve_smooth_iterations = atoi(std::string(argv[3]).c_str());
	nb_surface_smooth_iterations = atoi(std::string(argv[4]).c_str());
	nb_volume_smooth_iterations = atoi(std::string(argv[5]).c_str());
	edge_discretization = atoi(std::string(argv[6]).c_str());
	file_block_out = std::string(argv[7]);
	file_mesh_out = std::string(argv[8]);
	std::cout << "Parameters " << std::endl;
	std::cout << "  - Geometry file: " << file_geom << std::endl;
	std::cout << "  - Mesh file    : " << file_mesh << std::endl;
	std::cout << "  - Nb curve iterations: " << nb_curve_smooth_iterations << std::endl;
	std::cout << "  - Nb surface iterations: " << nb_surface_smooth_iterations << std::endl;
	std::cout << "  - Nb volume iterations: " << nb_volume_smooth_iterations << std::endl;
	std::cout << "  - Edge discretization: " << edge_discretization << std::endl;
	std::cout << "  - Output block file  : " << file_block_out << std::endl;
	std::cout << "  - Output mesh file  : " << file_mesh_out << std::endl;
	std::cout << "=======================================" << std::endl;*/

	/*=========================================================*/
	// CUBE SIMPLE POUR LE RUBIKSCUBE DE TAILLE 2x2x2
	Mesh cubeRubiksCube(MeshModel(DIM3 | R | F | E | N |
	                              R2N | R2F | R2E |
	                              F2N | F2R | F2E |
	                              E2F | E2N | E2R| N2E | N2R | N2F));


	Node crn0 = cubeRubiksCube.newNode(math::Point(0,0,0));
	Node crn1 = cubeRubiksCube.newNode(math::Point(4,0,0));
	Node crn2 = cubeRubiksCube.newNode(math::Point(4,0,4));
	Node crn3 = cubeRubiksCube.newNode(math::Point(0,0,4));
	Node crn4 = cubeRubiksCube.newNode(math::Point(0,4,0));
	Node crn5 = cubeRubiksCube.newNode(math::Point(4,4,0));
	Node crn6 = cubeRubiksCube.newNode(math::Point(4,4,4));
	Node crn7 = cubeRubiksCube.newNode(math::Point(0,4,4));

	Region crr1 = cubeRubiksCube.newHex(crn3,crn2,crn1,crn0,crn7,crn6,crn5,crn4);

	MeshDoctor cubeDoctor(&cubeRubiksCube);
	cubeDoctor.buildFacesAndR2F();
	cubeDoctor.buildEdgesAndX2E();
	cubeDoctor.updateUpwardConnectivity();

	Variable<int> *boundaryEdgeCubeRubiksCube = cubeRubiksCube.newVariable<int,gmds::GMDS_EDGE>("boundaryEdge");

	for (auto e : cubeRubiksCube.edges()){
		//std::cout<<"in edge : "<<e<<std::endl;
		//std::cout<<"Nb regions : "<<rubiksCube.get<Edge>(e).nbRegions()<<std::endl;
		if(cubeRubiksCube.get<Edge>(e).nbRegions()== 1){
			//std::cout<<"dans if "<<std::endl;
			boundaryEdgeCubeRubiksCube->set(e,1);
		}
		else{
			//std::cout<<"Dans else"<<std::endl;
			boundaryEdgeCubeRubiksCube->set(e,0);
		}
	}


	/*=========================================================*/
	// FORME GRILLE 2x2x2

	Mesh rubiksCube(MeshModel(DIM3 | R | F | E | N |
	                          R2N | R2F | R2E | F2N | F2R | F2E | E2F | E2N |E2R| N2E | N2R | N2F));


	GridBuilder gridBuilder(&rubiksCube,3);
	double step = 1;
	int nb = 5;
	int AXnb = nb;
	int AYnb = nb;
	int AZnb = nb;

	gridBuilder.execute(AXnb,step,AYnb,step,AZnb,step);
	MeshDoctor doctor(&rubiksCube);
	doctor.buildFacesAndR2F();
	doctor.buildEdgesAndX2E();
	doctor.updateUpwardConnectivity();

	Variable<int> *boundaryEdgeRubiksCube = rubiksCube.newVariable<int,gmds::GMDS_EDGE>("boundaryEdge");

	/*	for(auto r : rubiksCube.regions()){
	      Region region = rubiksCube.get<Region>(r);
	      std::cout<<"region "<< r <<std::endl;
	      std::vector<TCellID> edges;
	      region.delegateGetEdgeIDs(edges);
	      for(auto e : edges){
	         rubiksCube.get<Edge>(e).add(region);
	      }
	   }*/

	for (auto e : rubiksCube.edges()){
		//std::cout<<"in edge : "<<e<<std::endl;
		//std::cout<<"Nb regions : "<<rubiksCube.get<Edge>(e).nbRegions()<<std::endl;
		if(rubiksCube.get<Edge>(e).nbRegions()== 1){
			//std::cout<<"dans if "<<std::endl;
			boundaryEdgeRubiksCube->set(e,1);
		}
		else{
			//std::cout<<"Dans else"<<std::endl;
			boundaryEdgeRubiksCube->set(e,0);
		}
	}

	// DEBUT LINK !
	std::cout<<"> Start link"<<std::endl;
	cad::FACManager geom_manager;
	geom_manager.initFrom3DMesh(&cubeRubiksCube);

	//==================================================================
	// MARK ALL THE BOUNDARY CELL OF THE INIT MESH
	//==================================================================

	std::cout<<"> Start block boundary retrieval"<<std::endl;

	BoundaryOperator op(&rubiksCube);
	auto mark_node_NAN = rubiksCube.newMark<Node>();
	auto mark_node_on_pnt = rubiksCube.newMark<Node>();
	auto mark_node_on_crv = rubiksCube.newMark<Node>();
	auto mark_node_on_srf = rubiksCube.newMark<Node>();
	auto mark_edge_on_crv = rubiksCube.newMark<Edge>();
	auto mark_edge_on_srf = rubiksCube.newMark<Edge>();
	auto mark_face_on_srf = rubiksCube.newMark<Face>();

	op.markCellOnGeometry(mark_face_on_srf,
	                      mark_edge_on_srf,
	                      mark_node_on_srf,
	                      mark_edge_on_crv,
	                      mark_node_on_crv,
	                      mark_node_on_pnt,
	                      mark_node_NAN);

	IGMeshIOService ioService_geom(&cubeRubiksCube);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("surf_input_classification.vtk");


	std::cout<<"> Start block->geometry classification"<<std::endl;

	cad::GeomMeshLinker linker(&rubiksCube, &geom_manager);
	std::vector<cad::GeomSurface*> surfaces;
	std::vector<cad::GeomCurve*> curves;
	std::vector<cad::GeomPoint*> points;

	geom_manager.getSurfaces(surfaces);
	geom_manager.getCurves(curves);
	geom_manager.getPoints(points);

	//==================================================================
	//First, we classify each face

	for(auto f_id: rubiksCube.faces()){
		Face f= rubiksCube.get<Face>(f_id);
		if(rubiksCube.isMarked(f, mark_face_on_srf) ){
			//we've got a boundary face
			math::Point p = f.center();
			math::Vector3d v = f.normal();
			//we ensure to get the outward normal that goes from the inside block
			//towards the outside
			//as f is a boundary face, it is connected to only one block
			Region r = f.get<Region>()[0];
			math::Vector3d v2inside= r.center()-p;
			if(v.dot(v2inside)>0){
				//v points inside
				v = v.opp();
			}

			v.normalize();
			std::cout<<"FACE TO PROJECT "<<f_id<<" with center point  "<<p<<" and normal direction "<<v<<std::endl;
			double min_dist = 100000;
			int min_entity_dim=-1;
			int min_entity_id = -1;
			std::map<TCellID , double> surf_dist;
			std::map<TCellID , double> surf_dot;
			bool found_surf=false;
			for(auto s:surfaces){
				math::Point closest_pnt=p;
				s->project(closest_pnt);
				double dist = p.distance2(closest_pnt);
				math::Vector v_proj=closest_pnt-p;
				v_proj.normalize();
				surf_dot[s->id()]=(v.dot(v_proj));
				surf_dist[s->id()]=dist;

				if(dist<min_dist) {
					//either the points are the same or we check the vector orientation
					if (dist<1e-4 || fabs(v.dot(v_proj)) > 0.2 ) {
						bool keep_projection = true;

						if (v.dot(v_proj)<0){
							//we have to check that we do not cross the volume. It can occur for thin curved
							//shape (like cylinder with a hole)

							//We take other points closed to the face corners and we project all of them, if they go to
							// the same surf, we keep it. Otherwise, we will put another surface later.
							keep_projection=false;
							std::vector<Node> corners = f.get<Node>();
							std::vector<math::Point> corner_pnts;
							for(auto c:corners){
								corner_pnts.push_back(0.1*p+0.9* c.point());
							}
							std::set<int> corner_surf;
							for(auto cp:corner_pnts) {
								double cp_min_dist = 100000;
								auto cp_min_surf_id = -1;
								for (auto s: surfaces) {
									math::Point closest_pnt = cp;
									s->project(closest_pnt);

									double local_dist = cp.distance2(closest_pnt);
									if (local_dist < cp_min_dist) {
										cp_min_dist = local_dist;
										cp_min_surf_id = s->id();
									}
								}
								corner_surf.insert(cp_min_surf_id);
							}
							if(corner_surf.size()==1 && (*corner_surf.begin()==s->id())){
								//ok we project on it
								keep_projection= true;
							}

						}
						if(keep_projection) {
							min_dist = dist;
							min_entity_dim = 2;
							min_entity_id = s->id();
							found_surf = true;
						}
					}
				}
			}
			if(!found_surf){
				min_dist = 100000;
				min_entity_dim=-1;
				min_entity_id = -1;
				//means we have a block face in a concave area
				for(auto s:surfaces){
					math::Point closest_pnt=p;
					s->project(closest_pnt);

					double dist = p.distance2(closest_pnt);
					math::Vector v_proj=closest_pnt-p;
					v_proj.normalize();
					surf_dot[s->id()]=(v.dot(v_proj));
					surf_dist[s->id()]=dist;
					if(dist<min_dist) {
						min_dist = dist;
						min_entity_dim = 2;
						min_entity_id = s->id();
						found_surf=true;

					}
				}
			}
			if (min_entity_dim==2){
				std::cout<<"==> Link face "<<f_id<<" on surf "<<min_entity_id<<" with distance:" <<min_dist<<std::endl;
				linker.linkFaceToSurface(f_id,min_entity_id);
			}
			else{
				std::cout<<"\t ====> Link error for classifying a face"<<std::endl;
				throw GMDSException("Link error for classifying a face");
			}

		}
	}


	//==================================================================
	//Second, we classify each edge
	for(auto e_id: rubiksCube.edges()){
		Edge e= rubiksCube.get<Edge>(e_id);
		if(rubiksCube.isMarked(e, mark_edge_on_crv) ||
		    rubiksCube.isMarked(e, mark_edge_on_srf) ){
			//we've got a boundary edge, now we get the 2 boundary faces
			// around e
			std::vector<Node> e_nodes = e.get<Node>();
			std::vector<TCellID> adj_face_id = e.getIDs<Face>();
			std::vector<TCellID> adj_bnd_face_ids;
			for(auto f_id:adj_face_id){
				if(rubiksCube.isMarked<Face>(f_id, mark_face_on_srf)){
					adj_bnd_face_ids.push_back(f_id);
				}
			}
			if(adj_bnd_face_ids.size()!=2){
				std::cout<<"ERROR: One boundary edge is adjacent to more "
				          <<"than 2 boundary faces ("
				          <<"face "<<e_id<<")"<<std::endl;
			}
			auto surf0_id = linker.getGeomId<Face>(adj_bnd_face_ids[0]);
			auto surf1_id = linker.getGeomId<Face>(adj_bnd_face_ids[1]);

			if(surf0_id==surf1_id){
				//edge embedded in the surface
				linker.linkEdgeToSurface(e_id,surf1_id);
			}
			else{
				//it must be an edge classified on curves
				// Warning, like for the face, it is not necessary the geometrically closest to be connected to.
				//We check which curve is bounded by surf0_id and surf1_id
				bool found = false;

				std::vector<cad::GeomCurve*> candidate_curves;
				for(auto ci:curves) {
					std::vector<cad::GeomSurface *> c_surfs = ci->surfaces();
					if (c_surfs.size() == 2) {
						if ((c_surfs[0]->id() == surf0_id && c_surfs[1]->id() == surf1_id) ||
						    (c_surfs[1]->id() == surf0_id && c_surfs[0]->id() == surf1_id)) {
							std::cout<<"Curve "<<ci->id()<<" adj to surfaces "<<c_surfs[0]->id()<<" and "<<c_surfs[1]->id()<<std::endl;
							found = true;
							candidate_curves.push_back(ci);
						}
					}
				}
				if(found) {
					if(candidate_curves.size()==1) {
						linker.linkEdgeToCurve(e_id, candidate_curves[0]->id());
					}
					else{
						//we project on each of curve and keep the closest one
						double min_dist = 1e6;
						cad::GeomCurve* selected_curve = NULL;
						math::Point center_edge = e.center();
						for(auto ci:candidate_curves) {
							math::Point pi = center_edge;
							ci->project(pi);
							auto di = pi.distance2(center_edge);
							if(di<min_dist){
								min_dist=di;
								selected_curve = ci;
							}
						}

						linker.linkEdgeToCurve(e_id, selected_curve->id());

					}
				}
				else{
					throw GMDSException("BlockMesher error: impossible to link a block edge onto the geometry");
				}
			}
		}
	}




	//==================================================================
	//we classify each node
	auto on_pnt=0, on_curve=0, on_surf=0;
	for(auto n_id: rubiksCube.nodes()){
		Node n= rubiksCube.get<Node>(n_id);
		if(rubiksCube.isMarked(n, mark_node_on_pnt) ||
		    rubiksCube.isMarked(n, mark_node_on_crv) ||
		    rubiksCube.isMarked(n, mark_node_on_srf) ){
			//we've got a boundary node
			//As face and edges are already classified, we used it
			// A node is on a surface if all its adjacent boundary faces are on the same
			// surface. If they are on two it is on a curve potentially
			std::vector<TCellID> adj_faces = n.getIDs<Face>();
			std::vector<TCellID> adj_bnd_faces;
			for(auto id_face:adj_faces){
				if(rubiksCube.isMarked<Face>(id_face, mark_face_on_srf)) {
					adj_bnd_faces.push_back(id_face);
				}
			}
			std::set<int> bnd_surfaces;

			for (auto  id_face:adj_bnd_faces){
				bnd_surfaces.insert(linker.getGeomId<Face>(id_face));
			}
			int min_entity_dim=-1;
			int min_entity_id = -1;
			double min_dist = 100000;

			if(bnd_surfaces.size()==1){
				//on a single surface
				min_entity_dim=2;
				min_entity_id=*(bnd_surfaces.begin());
			}
			else if(bnd_surfaces.size()==2){
				//on a curve or a point that is connected to a single curve
				math::Point node_loc = n.point();
				for(auto p:points){
					math::Point closest_pnt=node_loc;
					double dist = node_loc.distance(p->point());
					if(dist<min_dist){
						min_dist =dist;
						min_entity_dim=0;
						min_entity_id=p->id();
					}
				}
				for(auto c:curves){
					math::Point closest_pnt=node_loc;
					c->project(closest_pnt);
					double dist = node_loc.distance(closest_pnt);
					if(dist<min_dist){
						//WARNING: Take care of this trick that is not good at all but mandatory to be staying on the
						// curve and not on the surface
						if(dist<1e-4)
							dist=0;
						min_dist =dist;
						min_entity_dim=1;
						min_entity_id=c->id();
					}
				}

			} else{
				//On a point!!
				math::Point node_loc = n.point();
				for(auto p:points){
					math::Point closest_pnt=node_loc;
					double dist = node_loc.distance(p->point());
					if(dist<min_dist){
						min_dist =dist;
						min_entity_dim=0;
						min_entity_id=p->id();
					}
				}
			}

			if(min_entity_dim==0){
				on_pnt++;
				linker.linkNodeToPoint(n_id,min_entity_id);
				std::cout<<"Node "<<n_id<<" is on point "<<min_entity_id<<std::endl;
			}
			else if (min_entity_dim==1){
				on_curve++;
				linker.linkNodeToCurve(n_id,min_entity_id);
			}
			else if (min_entity_dim==2){
				on_surf++;
				linker.linkNodeToSurface(n_id,min_entity_id);
			}
			else{
				throw GMDSException("Link error for classifying a node");
			}

		}
	}
	std::cout<<"  info [node classified on points "<<on_pnt;
	std::cout<<", on curves "<<on_curve;
	std::cout<<", on surfs "<<on_surf<<"]"<<std::endl;

	linker.writeVTKDebugMesh("linker_debug.vtk");

	//==================================================================
	// PERFORM THE BLOCK SMOOTHING NOW
	//==================================================================
	std::cout<<"> Start block smoothing"<<std::endl;
	smoothy::LaplacianSmoother smoother(&linker);
	if(!smoother.isValid())
	{
		std::cout<<"INVALID MODEL"<<std::endl;
		exit(1);
	}
	std::cout<<"  - start  smoothing ("<<nb_curve_smooth_iterations<<" iterations)"<<std::endl;
	smoother.smoothCurves(nb_curve_smooth_iterations);
	smoother.smoothSurfaces(nb_surface_smooth_iterations);
	// std::cout<<"  - start volume smoothing ("<<nb_volume_smooth_iterations<<" iterations)"<<std::endl;
	// smoother.smoothVolumes(nb_volume_smooth_iterations);

	//==================================================================
	// COMPUTE BLOCK QUALITY
	//==================================================================
	Variable<double>* var_quality = rubiksCube.newVariable<double,gmds::GMDS_REGION>("GMDS_Scaled_Jacobian");
	for(auto r_id:rubiksCube.regions()){
		Region r = rubiksCube.get<Region>(r_id);
		var_quality->value(r_id)=fabs(r.computeScaledJacobian());
	}
	std::cout << "> Write smoothed block mesh in: " << file_mesh_out << std::endl;
	IGMeshIOService ioService3(&rubiksCube);
	VTKWriter vtkWriter(&ioService3);
	vtkWriter.setCellOptions(N|R);
	vtkWriter.setDataOptions(N|R);
	vtkWriter.write(file_block_out);


	//==================================================================
	// NOW WE CREATE THE FINAL MESH
	//==================================================================
	//The first stage consists in performing a transfinite interpolation of blocks
	BlockMesher bm(&rubiksCube,&linker);
	std::cout<<"> Start mesh generation via transfinite interpolation"<<std::endl;
	bm.execute(edge_discretization);
	//All boundary node, edge and curves are classified.
	std::map<TCellID, Cell::Data> mesh_node_info = bm.mesh_node_classification();
	std::map<TCellID, Cell::Data> mesh_edge_info = bm.mesh_edge_classification();
	std::map<TCellID, Cell::Data> mesh_face_info = bm.mesh_face_classification();

	std::cout<<"> Start block smoothing"<<std::endl;
	bm.mesh()->changeModel(MeshModel(DIM3|R|N|E|F|R2N|F2N|E2N|N2F|N2E|N2R));
	MeshDoctor doc3(bm.mesh());
	doc3.updateUpwardConnectivity();
	cad::GeomMeshLinker linker_final(bm.mesh(),&geom_manager);
	for(auto info:mesh_node_info){
		Node n = bm.mesh()->get<Node>(info.first);
		if(info.second.dim==0) {
			linker_final.linkNodeToPoint(n.id(), info.second.id);
		}
		else if(info.second.dim==1) {
			linker_final.linkNodeToCurve(n.id(), info.second.id);
		}
		else if(info.second.dim==2)
			linker_final.linkNodeToSurface(n.id(),info.second.id);
	}
	for(auto info:mesh_edge_info){
		Edge e = bm.mesh()->get<Edge>(info.first);
		if(info.second.dim==1) {
			linker_final.linkEdgeToCurve(e.id(), info.second.id);
			std::cout<<"Edge "<<e.getIDs<Node>()[0]<<"-"<<e.getIDs<Node>()[1]<<" on curve "<<info.second.id<<std::endl;
		}
	}
	for(auto info:mesh_face_info){
		Face f = bm.mesh()->get<Face>(info.first);
		if(info.second.dim==2)
			linker_final.linkFaceToSurface(f.id(),info.second.id);
	}
	smoothy::LaplacianSmoother smoother_final(&linker_final);
	if(!smoother_final.isValid())
	{
		std::cout<<"INVALID MODEL"<<std::endl;
		exit(1);
	}
	std::cout<<"  - start smoothing ("<<nb_curve_smooth_iterations<<" iterations)"<<std::endl;
	smoother_final.smoothCurves(nb_curve_smooth_iterations);
	smoother_final.smoothSurfaces(nb_surface_smooth_iterations);
	std::cout<<"> Start final mesh writing"<<std::endl;
	IGMeshIOService ioService_mesh(bm.mesh());
	VTKWriter vtkWriter_mesh(&ioService_mesh);
	vtkWriter_mesh.setCellOptions(N|F);
	vtkWriter_mesh.setDataOptions(N|F);
	vtkWriter_mesh.write("file_mesh_out.vtk");

	VTKWriter vtkWriter_curves(&ioService_mesh);
	vtkWriter_curves.setCellOptions(N|E);
	vtkWriter_curves.setDataOptions(N|E);
	vtkWriter_curves.write("mesh_curves.vtk");

	std::cout << "======== Task done by blocker =========" << std::endl;
}