//
// Created by Paul Bourmaud on 23/03/2022.
//
#include <iostream>
#include "gmds/paul/Grid.h"
#include "gmds/ig/Mesh.h"
#include "gmds/igalgo/VolFracComputation.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKWriter.h"
#include "gmds/paul/Actions_Agent.h"
#include "gmds/paul/Tools.h"

using namespace gmds;

int main(){

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

	std::cout<<"================== Apres le test ========================="<<std::endl;
	for(auto face_id:mT.faces()){
		Face fT = mT.get<Face>(face_id);
		std::vector<Node> f_nodes = fT.get<Node>();
		//std::cout<<"Viens de print le node"<<std::endl;
		for(auto n:f_nodes){
			//std::cout<<"Avant Print boucle for auto 2"<<std::endl;
			std::cout<<n<<std::endl;
			//std::cout<<"Apres Print boucle for auto 2"<<std::endl;
		}
	}
	//test VolFrac function create by nicolas legoff (pas possible pour le moment parce que pas bon type
	// de maillage en entrée, besoin maillage tétraédrique je crois)
	/*std::cout<<"=== Valeur Fraction ==="<<std::endl;
	Variable<double> *volFrac = mT.newVariable<double,GMDS_FACE>("Valeur Fraction");
	gmds::volfraccomputation_2d(&mT,&m,volFrac);
	std::cout<<volFrac<<std::endl;*/

	//Create the grid around the target (triangle in this case)
	std::cout<<"=== Test Grid Builder Around ==="<<std::endl;
	GridBuilderAround gba(&m,2);
	gba.executeGrid2D(5,3,5,3);
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

	for(int id=0;id <2;id++) {
		action.executeDeleteFace(id);
	}
	/*for (auto face_id : gba.m_mesh.faces()) {
		std::cout << "Valeur Activate :" << action.getValueActivateFace(face_id)<<"\nFace_ID :"<<face_id << std::endl;
	}*/

	Tools tool(&gba);
	/*for (auto nodes_id: gba.m_mesh.nodes()){
		tool.getListFacesOfNode(nodes_id);
	}*/
	/*for (int i=0;i < 20;i++) {
		int i1alea,i2alea,faceAlea;
		i1alea = rand() % 16;
		i2alea = rand() % 16;
		faceAlea=rand() % 9;
		std::cout<<"Noeud 1 :"<<i1alea<<"\nNoeud 2 :"<<i2alea<<"\nFace :"<<tool.getIdOneCommonFace(i1alea,i2alea)<<std::endl;
		tool.checkExistEdge(i1alea,i2alea,tool.getIdOneCommonFace(i1alea,i2alea));
	}*/
	//tool.getOtherNodes(2,7);
	//tool.getListSecondNodesChain(11,16);
	tool.creatMiddleNode(tool.g_grid.m_mesh.get<Node>(1),tool.g_grid.m_mesh.get<Node>(2));

// Save Triangle Generation
	IGMeshIOService ioService(&mT);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("tata.vtk");


// Save GridBuilder Generation
	IGMeshIOService ioService2(&m);
	VTKWriter vtkWriter2(&ioService2);
	vtkWriter2.setCellOptions(gmds::N|gmds::F);
	vtkWriter2.setDataOptions(gmds::N|gmds::F);
	vtkWriter2.write("toto.vtk");

	exit(3);

}