//
// Created by Paul Bourmaud on 23/03/2022.
//
#include <cstdlib>
#include <iostream>
#include <ctime>

#include "gmds/paul/Grid.h"
#include "gmds/ig/Mesh.h"
#include "gmds/igalgo/VolFracComputation.h"
#include "gmds/paul/Actions_Agent.h"
#include "gmds/paul/Tools.h"
#include "gmds/paul/Environment.h"
#include "gmds/paul/Q_Learning_Algo.h"
#include "gmds/paul/Politique.h"

#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKWriter.h"
#include "gmds/io/VTKReader.h"
#include "gmds/ig/MeshDoctor.h"

using namespace gmds;

int main(){


	srand (time(NULL));

	//Create the agent environment
	std::cout<<"=== Test Grid ==="<<std::endl;
	Grid g(2,4);
	std::cout<<g.getX()<<" - "<<g.getY()<<std::endl;

	//Create triangle reference for BuildGridAround
	Mesh m(MeshModel(DIM3|F|N|F2N|N2F));
	Node n0 = m.newNode(math::Point(0,0,0));
	Node n1 = m.newNode(math::Point(0,3,0));
	Node n2 = m.newNode(math::Point(3,3,0));
	/*Variable<double> *v = m.newVariable<double,GMDS_NODE>("val");
	v->set(n0.id(),10);
	v->set(n1.id(),40);*/
	m.newTriangle(n0,n1,n2);
	for(auto face_id:m.faces()){
		Face f = m.get<Face>(face_id);
		std::vector<Node> f_nodes = f.get<Node>();
		for(auto n:f_nodes){
			std::cout<<n<<std::endl;
		}
	}

	std::cout<<"Creation triangle target"<<std::endl;
	//Create triangle target
	Mesh mT(MeshModel(DIM3|F|N|F2N|N2F));
	Node n0b = mT.newNode(math::Point(0,0,0));
	Node n1b = mT.newNode(math::Point(0,3,0));
	Node n2b = mT.newNode(math::Point(3,3,0));
	Variable<double> *vb = mT.newVariable<double,GMDS_NODE>("valTarget");
	vb->set(n0b.id(),10);
	vb->set(n1b.id(),30);
	mT.newTriangle(n0b,n1b,n2b);
	Variable<int> *titi = mT.newVariable<int,gmds::GMDS_FACE>("titi");
	titi->set(0,1);

	for (auto face_id:mT.faces()){
		Face fT=mT.get<Face>(face_id);
		std::vector<Node> f_nodes = fT.get<Node>();
		for (auto n : f_nodes){
			std::cout<<n.id()<<std::endl;
			vb->set(n.id(),40); //attribution val aux noeuds
		}
	}

	//test VolFrac function create by nicolas legoff
	std::cout<<"=== Valeur Fraction ==="<<std::endl;


	Mesh mImprint (	MeshModel(DIM2|F|E|N|F2N|N2F|E2N));
	IGMeshIOService ioServiceRead(&mImprint);
	VTKReader vtkReader(&ioServiceRead);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read("/home/bourmaudp/Documents/Stage_CEA_2022/Repo_GMDS/gmds/test_samples/HolesInSquare0.vtk");



/*
	gmds::MeshDoctor doc_imprint(&mImprint);
	doc_imprint.updateUpwardConnectivity();
	doc_imprint.buildBoundaryCells();
*/
	Mesh mGridAround(MeshModel (DIM2|F|N|F2N|N2F));
	IGMeshIOService ioServiceRead2(&mGridAround);
	VTKReader vtkReader2(&ioServiceRead2);
	vtkReader2.setCellOptions(gmds::N|gmds::F);
	vtkReader2.read("/home/bourmaudp/Documents/Stage_CEA_2022/Repo_GMDS/gmds/test_samples/HolesInSquare0.vtk");
	GridBuilderAround gridAround(&mGridAround,&mImprint,2);
	gridAround.executeGrid2D(5);


	Actions actionb(&gridAround);
	Tools toolsB(&gridAround);

	Mesh mGridAroundCopy(mGridAround.getModel());
	GridBuilderAround gridAroundCopy(&mGridAroundCopy,&mImprint,2);
	toolsB.cloneMesh(gridAround,gridAroundCopy);

	Actions actionCopy(&gridAroundCopy);
	Tools copyTools(&gridAroundCopy);



	//actionCopy.executeDeleteFace(gridAroundCopy.m_mesh.get<Face>(9).id());
	//actionCopy.executeCutFace(gridAroundCopy.m_mesh.get<Face>(2),0);
	//actionCopy.executeCutFace(gridAroundCopy.m_mesh.get<Face>(16),1);
	//auto edgeCopy = copyTools.getHorizontalEdge(gridAroundCopy.m_mesh.get<Face>(2));
	//std::cout<<"Edge Horiz Face 2 : "<<edgeCopy.front()<<" "<<edgeCopy.back()<<std::endl;
	//std::cout<<"Edge Horiz Face 2 : "<<edgeCopy.size()<<std::endl;

	gmds::IGMeshIOService ioService_write_Copy(&mGridAroundCopy);
	gmds::VTKWriter vtkWriterGba_Copy(&ioService_write_Copy);
	vtkWriterGba_Copy.setCellOptions(gmds::N|gmds::F);
	vtkWriterGba_Copy.setDataOptions(gmds::N|gmds::F);
	vtkWriterGba_Copy.write("VolFracComputation_test2DCopy.vtk");




	// =============================== TEST NOUVELLES ACTIONS SEULEMENT SUR LES FACES =================================
	/*
	actionb.executeCutFace(mGridAround.get<Face>(9),0);
	actionb.executeCutFace(mGridAround.get<Face>(20),1);
	actionb.executeGlidNodeFace(mGridAround.get<Face>(18));
	actionb.executeCutFace(mGridAround.get<Face>(18),1);
	actionb.executeGlidNodeFace(mGridAround.get<Face>(20));*/

	//actionb.executeGlideMaxNodeFace(mGridAround.get<Face>(12));
	//actionb.executeGlideMinNodeFace(gridAround.m_mesh.get<Face>(12));


	MeshDoctor doc(&mGridAround);
	//Variable<double>* volFrac = mGridAround.newVariable<double,GMDS_FACE>("volFrac");
	Variable<double>* IoU = gridAround.meshTarget.newVariable<double,GMDS_FACE>("IoU");

	toolsB.intersectionTargetWithGrid(&gridAround.m_mesh,&gridAround.meshTarget,IoU);


/*
	//coupe
	actionb.executeCutEdge(mGridAround.get<Node>(21),mGridAround.get<Node>(20));
	actionb.executeCutEdge(mGridAround.get<Node>(5),mGridAround.get<Node>(10));
	actionb.executeCutEdge(mGridAround.get<Node>(4),mGridAround.get<Node>(9));
	actionb.executeCutEdge(mGridAround.get<Node>(4),mGridAround.get<Node>(3));
	actionb.executeCutEdge(mGridAround.get<Node>(19),mGridAround.get<Node>(24));
	actionb.executeCutEdge(mGridAround.get<Node>(22),mGridAround.get<Node>(21));
	actionb.executeCutEdge(mGridAround.get<Node>(4),mGridAround.get<Node>(41));
	actionb.executeCutEdge(mGridAround.get<Node>(3),mGridAround.get<Node>(2));
	actionb.executeCutEdge(mGridAround.get<Node>(14),mGridAround.get<Node>(19));
	actionb.executeCutEdge(mGridAround.get<Node>(40),mGridAround.get<Node>(8));
	actionb.executeCutEdge(mGridAround.get<Node>(8),mGridAround.get<Node>(34));
	actionb.executeCutEdge(mGridAround.get<Node>(72),mGridAround.get<Node>(2));
	actionb.executeCutEdge(mGridAround.get<Node>(2),mGridAround.get<Node>(56));

	//Glide

	actionb.executeGlideNode(mGridAround.get<Node>(130),&mImprint);
	actionb.executeGlideNode(mGridAround.get<Node>(61),&mImprint);
	actionb.executeGlideNode(mGridAround.get<Node>(84),&mImprint);
	actionb.executeGlideNode(mGridAround.get<Node>(83),&mImprint);
	actionb.executeGlideNode(mGridAround.get<Node>(52),&mImprint);
	actionb.executeGlideNode(mGridAround.get<Node>(11),&mImprint);
	actionb.executeGlideNode(mGridAround.get<Node>(27),&mImprint);

	//trou
	actionb.executeGlideNode(mGridAround.get<Node>(75),&mImprint);
	actionb.executeGlideNode(mGridAround.get<Node>(95),&mImprint);
	actionb.executeGlideNode(mGridAround.get<Node>(111),&mImprint);
	actionb.executeGlideNode(mGridAround.get<Node>(94),&mImprint);
	actionb.executeGlideNode(mGridAround.get<Node>(123),&mImprint);
	actionb.executeGlideNode(mGridAround.get<Node>(93),&mImprint);
	actionb.executeGlideNode(mGridAround.get<Node>(58),&mImprint);
	actionb.executeGlideNode(mGridAround.get<Node>(102),&mImprint);
	actionb.executeGlideNode(mGridAround.get<Node>(125),&mImprint);
	actionb.executeGlideNode(mGridAround.get<Node>(103),&mImprint);
	actionb.executeGlideNode(mGridAround.get<Node>(113),&mImprint);
	actionb.executeGlideNode(mGridAround.get<Node>(104),&mImprint);*/

	//Supp face
	//trou
	/*actionb.executeDeleteFace(0);
	actionb.executeDeleteFace(96);
	actionb.executeDeleteFace(8);
	actionb.executeDeleteFace(89);
	actionb.executeDeleteFace(80);
	actionb.executeDeleteFace(40);
	actionb.executeDeleteFace(81);
	actionb.executeDeleteFace(97);

	//rebord droit
	actionb.executeDeleteFace(26);
	actionb.executeDeleteFace(10);
	actionb.executeDeleteFace(27);
	actionb.executeDeleteFace(53);
	actionb.executeDeleteFace(55);
	actionb.executeDeleteFace(118);
	actionb.executeDeleteFace(120);
	actionb.executeDeleteFace(119);
	actionb.executeDeleteFace(15);
	actionb.executeDeleteFace(11);
	actionb.executeDeleteFace(42);
	actionb.executeDeleteFace(43);
	actionb.executeDeleteFace(28);
	actionb.executeDeleteFace(29);
	actionb.executeDeleteFace(10);*/

	volfraccomputation_2d(&mGridAround,&mImprint,mGridAround.getVariable<double,GMDS_FACE>("volFrac"));

	Environment environment(&gridAround,&gridAround.meshTarget,&actionb);
	//environment.executeAction(environment.faceSelect(),1);

	Politique policy(&environment);






	//std::cout<<"interval : "<<policy.getInterval(1)<<std::endl;



	environment.localIoU(mGridAround.get<Face>(6));

	//===================================== TEST ACTIONS CUT FACE ET DIRECTION ===============================
	/*actionb.executeCutFace(gridAround.m_mesh.get<Face>(6),1);
	actionb.executeCutFace(gridAround.m_mesh.get<Face>(16),0);
	actionb.executeCutFace(gridAround.m_mesh.get<Face>(2),1);
	actionb.executeCutFace(gridAround.m_mesh.get<Face>(31),1);*/


	environment.localIoU(mGridAround.get<Face>(6));
	std::cout<<"Reward : "<<environment.reward(mGridAround.get<Face>(6))<<std::endl;



	//0.745944


	//+===================================== TEST BOUNDARY EXTRACTOR =======================================
	/*
	Mesh boundaryMesh(MeshModel (DIM2|E|N|N2E|E2N));
	BoundaryExtractor2D boundary_extractor(&mImprint,&boundaryMesh);
	std::vector<Node> listBoundaryNodes={};
	std::map<TCellID ,TCellID > aNodeMap;
	std::map<TCellID ,TCellID > aEdgeMap;
	std::map<TCellID ,TCellID > aNodeMapInv;
	std::map<TCellID ,TCellID > aEdgeMapInv;
	boundary_extractor.setMappings(&aNodeMap,&aEdgeMap,&aNodeMapInv,&aEdgeMapInv);
	boundary_extractor.execute();

	std::cout<<"$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
	if(aNodeMap.empty()){
		std::cout<<"LA MAP EST VIDE"<<std::endl;
	}
	for (auto n : aNodeMap){
		std::cout<<"N First : "<<n.first<< " " << "N Second : "<<n.second<<std::endl;

	}*/

	//=======================================================================================================================




	gmds::IGMeshIOService ioService_write(&mGridAround);
	gmds::VTKWriter vtkWriterGba(&ioService_write);
	vtkWriterGba.setCellOptions(gmds::N|gmds::F);
	vtkWriterGba.setDataOptions(gmds::N|gmds::F);
	vtkWriterGba.write("VolFracComputation_test2D.vtk");

	gmds::IGMeshIOService ioService_write_target(&gridAround.meshTarget);
	gmds::VTKWriter vtkWriterGbaTarget(&ioService_write_target);
	vtkWriterGbaTarget.setCellOptions(gmds::N|gmds::F);
	vtkWriterGbaTarget.setDataOptions(gmds::N|gmds::F);
	vtkWriterGbaTarget.write("VolFracComputation_Target.vtk");


	//Create the grid around the target (triangle in this case)
	std::cout<<"=== Test Grid Builder Around ==="<<std::endl;

	Mesh meshTarget (	MeshModel(DIM3|F|N|F2N|N2F));
	IGMeshIOService ioServiceRead1(&meshTarget);
	VTKReader vtkReader1(&ioServiceRead1);
	vtkReader1.setCellOptions(gmds::N|gmds::F);
	vtkReader1.read("/home/bourmaudp/Documents/Stage_CEA_2022/Repo_GMDS/gmds/test_samples/HolesInSquare0.vtk");

	GridBuilderAround gba(&meshTarget,&mImprint,2);
	int nbNodes = 5;
	gba.executeGrid2D(nbNodes);
	Actions action(&gba);
	//Variable<int> *activate = m.newVariable<int,GMDS_FACE>("exist");

	for(auto face_id:gba.m_mesh.faces()){
		Face f =gba.m_mesh.get<Face>(face_id);
		std::vector<Node>f_nodes = f.get<Node>();
		//std::cout<<"Avant bool activate"<<std::endl;
		//std::cout<<activate<<std::endl;
		//std::cout<<"Apres bool activate"<<std::endl;
		for(auto n:f_nodes){
			//std::cout<<n<<std::endl;
			//std::cout<<"print de f.id()"<<std::endl;
			//std::cout<<f.id()<<std::endl;
			//std::cout<<"================================================="<<std::endl;
			//gba.activate->set(f.id(), 1);
		}
	}

	/*
	for(int id=0;id <9;id++) {
		action.executeDeleteFace(id);
	}*/

	Tools tool(&gba);

	/*int i1=4;
	int i2=9;
	tool.checkVertical(tool.g_grid.m_mesh.get<Node>(i1),tool.g_grid.m_mesh.get<Node>(i2));
	*/

	//======================================== TEST CUT NODE GLIDE =============================================
/*
	gba.m_mesh.get<Node>(4).X()=0.5;
	gba.m_mesh.get<Node>(4).Y()=4.5;

	action.executeCutEdge(gba.m_mesh.get<Node>(4),gba.m_mesh.get<Node>(9));
*/

	//======================================== TEST MULTI CUT  =============================================
	/*
	action.executeCutEdge(tool.g_grid.m_mesh.get<Node>(18),tool.g_grid.m_mesh.get<Node>(13));
	action.executeCutEdge(tool.g_grid.m_mesh.get<Node>(1),tool.g_grid.m_mesh.get<Node>(2));
	action.executeCutEdge(tool.g_grid.m_mesh.get<Node>(3),tool.g_grid.m_mesh.get<Node>(2));
	action.executeCutEdge(tool.g_grid.m_mesh.get<Node>(14),tool.g_grid.m_mesh.get<Node>(29));
	action.executeCutEdge(tool.g_grid.m_mesh.get<Node>(19),tool.g_grid.m_mesh.get<Node>(29));
	action.executeCutEdge(tool.g_grid.m_mesh.get<Node>(4),tool.g_grid.m_mesh.get<Node>(3));
	*/

	//====================== TEST =====================
		/*std::srand(std::time(nullptr));
		for (int i = 0; i < 1000; i++) {
			int i1 = rand() % 82;
			std::cout<<i1<<std::endl;
			int i2 = rand() % 82;
			std::cout<<i2<<std::endl;
			Node node1 = gba.m_mesh.get<Node>(i1);
			Node node2 = gba.m_mesh.get<Node>(i2);
			action.executeCutEdge(node1, node2);
		}

		for (int i; i < 18; i++) {
			int idFace = rand() % 50;
			action.executeDeleteFace(idFace);
		}*/

	// exemple lecture fichier vtk
	   /*Mesh meshTarget (	MeshModel(DIM3|F|N|F2N|N2F));
	   IGMeshIOService ioServiceRead(&meshTarget);
	   VTKReader vtkReader(&ioServiceRead);
	   vtkReader.setCellOptions(gmds::N|gmds::R);
	   vtkReader.read("/home/bourmaudp/Documents/Stage_CEA_2022/Repo_GMDS/gmds/test_samples/A.vtk");
*/





// Save Triangle Generation
	IGMeshIOService ioService(&mT);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("lulu.vtk");


// Save GridBuilder Generation
	IGMeshIOService ioService2(&meshTarget);
	VTKWriter vtkWriter2(&ioService2);
	vtkWriter2.setCellOptions(gmds::N|gmds::F);
	vtkWriter2.setDataOptions(gmds::N|gmds::F);
	vtkWriter2.write("lili.vtk");


	//executeTrainQlearning(environment);

	exit(3);

}
