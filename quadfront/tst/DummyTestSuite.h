#include <gtest/gtest.h>
#include <gmds/quadfront/Quadfront.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

TEST(DummyTestClass, aaa)
{

	ASSERT_EQ(0, 0);
	quadfront::Quadfront quad("/home/pagea/Documents/Travail/data/exemple.vtk");

	quad.initialise_nodeBoundary();

	for (auto n : quad.m_edgeBoundary){
		quad.edgeRecovery(n.first);
		}


	// TEST NODE

	/*
	for (auto n : quad.m_mesh->nodes()) {
		Node node = quad.m_mesh->get<Node>(n);
		std::tuple<int, Edge, Edge> map = quad.get_nodeBoundary(node);
		std::cout<<"Ici : "<<std::get<0>(map)<<" ; "<<std::get<1>(map).id()<<" ; "<<std::get<2>(map).id()<<std::endl;
	}

	for (auto i : quad.m_nodeBoundary){
		std::cout<<"\n"<<std::endl;
		std::cout<<"Node : "<<i.first<<std::endl;
		std::cout<<"Frontière : "<<std::get<0>(i.second)<<std::endl;
		std::cout<<"Edge 1 : "<<std::get<1>(i.second).id()<<std::endl;
		std::cout<<"Edge 2 : "<<std::get<2>(i.second).id()<<"\n"<<std::endl;
	}
	 */

/*
	// TEST FACE
	for (auto f : quad.m_mesh->faces()) {
		Face face = quad.m_mesh->get<Face>(f);
		std::map<Node, std::pair<Edge,Edge>> ex = quad.get_edgeInFace(face);
	}
*/
}

TEST(DummyTestClass, swap){
	ASSERT_EQ(0, 0);

	quadfront::Quadfront quad("/home/pagea/Documents/Travail/data/swap.vtk");

	std::cout<<"================= Avant ================="<<std::endl;

	Node Anode;
	Node Anode_opposit;

	for (auto& n : quad.m_mesh->edges()){
		Edge edge = quad.m_mesh->get<Edge>(n);
		std::cout<<"Noeud 1 : "<<edge.get<Node>()[0]<<" Noeud 2 : "<<edge.get<Node>()[1]<<std::endl;
	}

	for (auto& e : quad.m_nodeBoundary){
		if (e.first.X()==0 and e.first.Y()==0 and e.first.Z()==0){
			Anode = e.first;
		}
		if (e.first.X()==1 and e.first.Y()==1 and e.first.Z()==0){
			Anode_opposit = e.first;
		}
	}

	Face Face_anode = Anode.get<Face>()[0];
	Face Face_opposit = Anode_opposit.get<Face>()[0];

	std::map<Node, std::pair<Edge, Edge>> node2edge = quad.get_edgeInFace(Face_anode);
	std::map<Node, std::pair<Edge, Edge>> node2edge_opposit = quad.get_edgeInFace(Face_opposit);

	Edge edge_inter = quad.get_opposit2Node(Anode, Face_anode);

	quad.swap(edge_inter, Anode, Anode_opposit, node2edge_opposit, std::get<0>(node2edge[Anode]), std::get<1>(node2edge[Anode]));

	gmds::IGMeshIOService ioService(quad.m_mesh);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("/home/pagea/Documents/Travail/data/toto_link.vtk");

	std::cout<<"================= Après ================="<<std::endl;

	for (auto& n : quad.m_mesh->edges()){
		Edge edge = quad.m_mesh->get<Edge>(n);
		std::cout<<"Noeud 1 : "<<edge.get<Node>()[0]<<" Noeud 2 : "<<edge.get<Node>()[1]<<std::endl;
	}
}

TEST(DummyTestClass, split){
	ASSERT_EQ(0, 0);

	quadfront::Quadfront quad("/home/pagea/Documents/Travail/data/swap.vtk");

	std::cout<<"================= Avant ================="<<std::endl;

	Node Anode;
	Node Anode_opposit;


	for (auto& n : quad.m_mesh->edges()){
		Edge edge = quad.m_mesh->get<Edge>(n);
		std::cout<<"Edge : "<<edge.id()<<std::endl;
		std::cout<<edge.get<Node>()[0]<<" : "<<edge.get<Node>()[1]<<std::endl;
	}

	for (auto& e : quad.m_nodeBoundary){
		if (e.first.X()==0 and e.first.Y()==0 and e.first.Z()==0){
			Anode = e.first;
		}
		if (e.first.X()==1 and e.first.Y()==1 and e.first.Z()==0){
			Anode_opposit = e.first;
		}
	}

	Face Face_anode = Anode.get<Face>()[0];
	Face Face_opposit = Anode_opposit.get<Face>()[0];

	std::map<Node, std::pair<Edge, Edge>> node2edge = quad.get_edgeInFace(Face_anode);
	std::map<Node, std::pair<Edge, Edge>> node2edge_opposit = quad.get_edgeInFace(Face_opposit);

	Edge edge_inter = quad.get_opposit2Node(Anode, Face_anode);

	Node node_inter = quad.m_mesh->newNode(0.5, 0.5, 0);

	std::cout<<"================= Après ================="<<std::endl;

	quad.split(node_inter, edge_inter, Anode, Anode_opposit, node2edge_opposit, std::get<0>(node2edge[Anode]), std::get<1>(node2edge[Anode]));

	gmds::IGMeshIOService ioService(quad.m_mesh);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("/home/pagea/Documents/Travail/data/toto_link.vtk");

	std::cout<<"================= Après ================="<<std::endl;

	for (auto& n : quad.m_mesh->edges()){
		Edge edge = quad.m_mesh->get<Edge>(n);
		std::cout<<"Edge : "<<edge.id()<<std::endl;
		std::cout<<"Noeud 1 : "<<edge.get<Node>()[0]<<" Noeud 2 : "<<edge.get<Node>()[1]<<std::endl;
	}


	std::cout<<"\n"<<std::endl;

	for (auto& n : quad.m_mesh->faces()){
		Face face = quad.m_mesh->get<Face>(n);
		std::cout<<face<<std::endl;
		std::vector<Edge> edge = face.get<Edge>();
		for (auto& e : edge){
			std::cout<<e.id()<<std::endl;
			std::cout<<e.get<Node>()[0]<<" : "<<e.get<Node>()[1]<<std::endl;
		}
		std::cout<<"\n"<<std::endl;
	}
}

