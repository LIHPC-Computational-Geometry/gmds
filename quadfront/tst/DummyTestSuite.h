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

		std::cout<<"---------- edge_boundary --------"<<std::endl;
		std::cout<<n.first.id()<<" : "<<n.first.get<Node>()[0]<<" ; "<<n.first.get<Node>()[1]<<std::endl;
		std::cout<<"\n"<<std::endl;
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

TEST(DummyTestClass, intersection){

	quadfront::Quadfront quad("/home/pagea/Documents/Travail/data/intersection.vtk");

	Node node = quad.m_mesh->get<Node>(9);
	std::cout<<node<<std::endl;

	Node AnodeRecovery0 = quad.m_mesh->get<Node>(1);
	std::cout<<AnodeRecovery0<<std::endl;
	Node AnodeRecovery1 = quad.m_mesh->get<Node>(3);
	std::cout<<AnodeRecovery1<<std::endl;

	std::vector<Edge> vEdgeIntersect = quad.interestionS(AnodeRecovery0, AnodeRecovery1);

	gmds::IGMeshIOService ioService(quad.m_mesh);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F|gmds::E);
	vtkWriter.setDataOptions(gmds::N|gmds::F|gmds::E);
	vtkWriter.write("/home/pagea/Documents/Travail/data/toto_link.vtk");

	std::cout<<"intersection : "<<vEdgeIntersect.size()<<std::endl;
}

TEST(DummyTestClass, edgeRecovery) {
	quadfront::Quadfront quad("/home/pagea/Documents/Travail/data/intersection.vtk");

	quad.initialise_nodeBoundary();

	Edge edge0 = quad.m_mesh->get<Edge>(0);
	std::cout<<" edge id 0 : "<<edge0.get<Node>()[0]<<" ; "<<edge0.get<Node>()[1]<<std::endl;
	std::cout<<"\n=== Voisin de Edge 0 ===\n "<<std::endl;
	for (auto& e : edge0.get<Face>()){
		std::cout<<e<<std::endl;
	}


	Edge edge = quad.m_edgeBoundary.begin()->first;

	std::cout<<"edge m_edge boundary "<<edge.id()<<" : "<< edge.get<Node>()[0]<<" ; "<<edge.get<Node>()[1]<<std::endl;


	quad.edgeRecovery(edge);

	gmds::IGMeshIOService ioService(quad.m_mesh);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F|gmds::E);
	vtkWriter.setDataOptions(gmds::N|gmds::F|gmds::E);
	vtkWriter.write("/home/pagea/Documents/Travail/data/toto_link.vtk");

}

TEST(DummyTestClass, superpose) {
	quadfront::Quadfront quad("/home/pagea/Documents/Travail/data/swap.vtk");

	Edge edge0 = quad.m_mesh->get<Edge>(3);
	std::cout<<" edge0 : "<<edge0.get<Node>()[0]<<" ; "<<edge0.get<Node>()[1]<<std::endl;
	Edge edge1 = quad.m_mesh->get<Edge>(0);
	std::cout<<" edge1 : "<<edge1.get<Node>()[0]<<" ; "<<edge1.get<Node>()[1]<<std::endl;
	Edge edge2 = quad.m_mesh->get<Edge>(1);
	std::cout<<" edge2 : "<<edge2.get<Node>()[0]<<" ; "<<edge2.get<Node>()[1]<<std::endl;
	Edge edge3 = quad.m_mesh->get<Edge>(4);
	std::cout<<" edge3 : "<<edge3.get<Node>()[0]<<" ; "<<edge3.get<Node>()[1]<<std::endl;
	Node node0 = quad.m_mesh->get<Node>(2);


	std::vector<Face> vFace = quad.createQuad(node0, edge0, edge1, edge3, edge2);
	std::cout<<" size : "<<vFace.size()<<std::endl;


}

