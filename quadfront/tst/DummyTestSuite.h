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


	/*
	for (auto n : quad.m_edgeBoundary){

		std::cout<<"---------- edge_boundary --------"<<std::endl;
		std::cout<<n.first.id()<<" : "<<n.first.get<Node>()[0]<<" ; "<<n.first.get<Node>()[1]<<std::endl;
		std::cout<<"\n"<<std::endl;
		quad.edgeRecovery(n.first);
		}
	 */


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

	for (auto& e : quad.m_nodeFront){
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

	//quad.swap(edge_inter, Anode, Anode_opposit, node2edge_opposit, std::get<0>(node2edge[Anode]), std::get<1>(node2edge[Anode]));


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

	for (auto& e : quad.m_nodeFront){
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

	//quad.split(node_inter, edge_inter, Anode, Anode_opposit, node2edge_opposit, std::get<0>(node2edge[Anode]), std::get<1>(node2edge[Anode]));

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

	Node AnodeRecovery0 = quad.m_mesh->get<Node>(9);
	std::cout<<AnodeRecovery0<<std::endl;
	Node AnodeRecovery1 = quad.m_mesh->get<Node>(5);
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

	Edge edge0 = quad.m_mesh->get<Edge>(0);
	std::cout<<" edge id 0 : "<<edge0.get<Node>()[0]<<" ; "<<edge0.get<Node>()[1]<<std::endl;
	std::cout<<"\n=== Voisin de Edge 0 ===\n "<<std::endl;
	for (auto& e : edge0.get<Face>()){
		std::cout<<e<<std::endl;
	}


	//Edge edge = quad.m_edgeBoundary.begin()->first;

	//std::cout<<"edge m_edge boundary "<<edge.id()<<" : "<< edge.get<Node>()[0]<<" ; "<<edge.get<Node>()[1]<<std::endl;


	//quad.edgeRecovery(edge);

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
	std::cout<<" node0 : "<<node0<<std::endl;
	Node node1 = quad.m_mesh->get<Node>(3);
	std::cout<<" node1 : "<<node1<<std::endl;

	std::vector<Face> vFace;
	std::vector<Face> vFace0 = node0.get<Face>();
	std::vector<Face> vFace1 = node1.get<Face>();
	std::set_intersection(vFace0.begin(), vFace0.end(), vFace1.begin(), vFace1.end(), std::back_inserter(vFace));
	for (auto& e : vFace){
		std::cout<<e<<std::endl;
	}
}


TEST(DummyTestClass, createQuad){
	quadfront::Quadfront quad("/home/pagea/Documents/Travail/data/swap.vtk");


	Edge e0 = quad.m_mesh->get<Edge>(0);
	Edge e1 = quad.m_mesh->get<Edge>(1);
	Edge e2 = quad.m_mesh->get<Edge>(3);
	Edge e3 = quad.m_mesh->get<Edge>(4);
	std::cout<<" e0 : "<<e0.get<Node>()[0]<<" ; "<<e0.get<Node>()[1]<<std::endl;
	std::cout<<" e1 : "<<e1.get<Node>()[0]<<" ; "<<e1.get<Node>()[1]<<std::endl;
	std::cout<<" e2 : "<<e2.get<Node>()[0]<<" ; "<<e2.get<Node>()[1]<<std::endl;
	std::cout<<" e3 : "<<e3.get<Node>()[0]<<" ; "<<e3.get<Node>()[1]<<std::endl;

	Node node0 = quad.m_mesh->get<Node>(0);
	std::cout<<" node0 : "<<node0<<"\n"<<std::endl;

	for(auto& f : quad.m_mesh->faces()){
		std::cout<<quad.m_mesh->get<Face>(f)<<std::endl;
	}

	quad.createQuad(node0, e0, e1, e2, e3);

	for(auto& f : quad.m_mesh->faces()){
		Face face = quad.m_mesh->get<Face>(f);
		std::cout<<face<<std::endl;
		for (auto& e : face.get<Edge>()){
			std::cout<<" edge : "<<e.get<Node>()[0]<<" ; "<<e.get<Node>()[1]<<std::endl;
		}
	}


	gmds::IGMeshIOService ioService(quad.m_mesh);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F|gmds::E);
	vtkWriter.setDataOptions(gmds::N|gmds::F|gmds::E);
	vtkWriter.write("/home/pagea/Documents/Travail/data/toto_link.vtk");
}



TEST(DummyTestClass, smooth1)
{
	quadfront::Quadfront mesh("/home/pagea/Documents/Travail/data/cercle2.vtk");

	mesh.initialise_Boundary();

	Variable<int> *front = mesh.m_mesh->newVariable<int, GMDS_EDGE>("front");
	Variable<int> *quad_var = mesh.m_mesh->newVariable<int, GMDS_FACE>("quad");

	std::cout << mesh.m_nodeFront.size() << std::endl;

	int pr = 1;

	while (pr <= 9) {
		std::vector<std::vector<Edge>> listBoundary;
		for (auto &l : mesh.m_listFront) {
			std::vector<Face> listQuad;
			std::vector<Edge> listEdge;
			std::vector<Edge> listEdge_bis;
			for (auto &e : l) {
				if (e.get<Face>()[0].nbEdges() == 3 or e.get<Face>()[1].nbEdges() == 3) {
					Face quad = mesh.qmorphAlgo(e);
					listQuad.push_back(quad);
					quad_var->set(quad.id(), pr);

					Edge edge_up = mesh.get_oppositEdge_Quad(e, quad);

					mesh.interiorSmoothing(edge_up.get<Node>()[0]);
					mesh.interiorSmoothing(edge_up.get<Node>()[1]);
					mesh.triangleSmoothing(edge_up);
				}
			}

			// Define the new front
			for (auto &f : listQuad) {
				for (auto &ef : f.get<Edge>()) {
					if (ef.get<Face>().size() == 2) {
						if (ef.get<Face>()[0].nbEdges() == 3 or ef.get<Face>()[1].nbEdges() == 3) {
							front->set(ef.id(), pr * 10);
							listEdge_bis.push_back(ef);

							mesh.interiorSmoothing(ef.get<Node>()[0]);
							mesh.interiorSmoothing(ef.get<Node>()[1]);

							Edge edge_base = mesh.get_oppositEdge_Quad(ef, f);

							if (edge_base.get<Face>().size() == 2) {
								mesh.interiorSmoothing(edge_base.get<Node>()[0]);
								mesh.interiorSmoothing(edge_base.get<Node>()[1]);
							}

							mesh.interiorSmoothing(ef.get<Node>()[0]);
							mesh.interiorSmoothing(ef.get<Node>()[1]);
						}
					}
				}
			}

			for (auto &ea : listEdge_bis) {
				mesh.triangleSmoothing(ea);
			}

			listBoundary.push_back(listEdge_bis);
		}

		mesh.m_listFront = listBoundary;
		pr += 1;

		mesh.initialise_Boundary();

		gmds::IGMeshIOService ioService(mesh.m_mesh);
		gmds::VTKWriter vtkWriter(&ioService);
		vtkWriter.setCellOptions(gmds::N | gmds::F | gmds::E);
		vtkWriter.setDataOptions(gmds::N | gmds::F | gmds::E);
		vtkWriter.write("/home/pagea/Documents/Travail/data/toto_link2.vtk");
	}
}


TEST(DummyTestClass, exemple){
	quadfront::Quadfront mesh("/home/pagea/Documents/Travail/data/petitangle.vtk");

	mesh.initialise_Boundary();

	std::vector<Face> listQuad_init;
	for(auto& i : mesh.m_nodeBoundary){
		Node node = i.first;
		mesh.initialiseQuad(node, listQuad_init);
	}




	Variable<int> *front = mesh.m_mesh->newVariable<int, GMDS_EDGE>("front");
	Variable<int> *quad_var = mesh.m_mesh->newVariable<int, GMDS_FACE>("quad");

	std::cout<<mesh.m_nodeFront.size()<<std::endl;

	int pr = 1;

	while (pr <= 1) {
		std::vector<std::vector<Edge>> listBoundary;
		for (auto &l : mesh.m_listFront) {
			std::vector<Face> listQuad;
			std::vector<Edge> listEdge;
			std::vector<Edge> listEdge_bis;
			for (auto &e : l) {
				if ((e.get<Face>().size() == 2 and (e.get<Face>()[0].nbEdges() == 3 or e.get<Face>()[1].nbEdges() == 3)) or
				    e.get<Face>().size() == 1 and (e.get<Face>()[0].nbEdges() == 3)){
					Face quad = mesh.qmorphAlgo(e);
					listQuad.push_back(quad);
					quad_var->set(quad.id(), pr);


					Edge edge_up = mesh.get_oppositEdge_Quad(e, quad);

					if(mesh.m_nodeBoundary.find(edge_up.get<Node>()[0]) == mesh.m_nodeBoundary.end()){
						mesh.interiorSmoothing(edge_up.get<Node>()[0]);
					}
					if(mesh.m_nodeBoundary.find(edge_up.get<Node>()[1]) == mesh.m_nodeBoundary.end()){
						mesh.interiorSmoothing(edge_up.get<Node>()[1]);
					}
					mesh.triangleSmoothing(edge_up);
				}
				for(auto& n : e.get<Node>()){
					for(auto& f : n.get<Face>()){
						if(f.nbEdges() == 4 and std::find(listQuad.begin(), listQuad.end(), f) == listQuad.end()){
							listQuad.push_back(f);
							quad_var->set(f.id(), pr);
						}
					}
				}
				gmds::IGMeshIOService ioService(mesh.m_mesh);
				gmds::VTKWriter vtkWriter(&ioService);
				vtkWriter.setCellOptions(gmds::N | gmds::F | gmds::E);
				vtkWriter.setDataOptions(gmds::N | gmds::F | gmds::E);
				vtkWriter.write("/home/pagea/Documents/Travail/data/toto_link3.vtk");
			}

			// Define the new front
			for (auto &f : listQuad) {
				for (auto &ef : f.get<Edge>()) {
					if (ef.get<Face>().size() == 2) {
						if (ef.get<Face>()[0].nbEdges() == 3 or ef.get<Face>()[1].nbEdges() == 3) {
							front->set(ef.id(), pr * 10);
							listEdge_bis.push_back(ef);
						}
					}
				}
			}

			if (listEdge_bis.size() != 0) {
				std::vector<double> stat = mesh.statFront(listEdge_bis);

				for (auto &ea : listEdge_bis) {
					mesh.blSmoothing(ea.get<Node>()[0], stat[0], stat[1]);
					mesh.blSmoothing(ea.get<Node>()[1], stat[0], stat[1]);

					Face quad = ea.get<Face>()[0].nbEdges() == 4 ? ea.get<Face>()[0] : ea.get<Face>()[1];

					Edge edge_base = mesh.get_oppositEdge_Quad(ea, quad);

					if (edge_base.get<Face>().size() == 2) {
						mesh.interiorSmoothing(edge_base.get<Node>()[0]);
						mesh.interiorSmoothing(edge_base.get<Node>()[1]);
					}
				}
			}

			for(auto& i : mesh.m_mesh->nodes()){
				Node node = mesh.m_mesh->get<Node>(i);
				if (mesh.m_nodeFront.find(node) == mesh.m_nodeFront.end() and
				    mesh.m_nodeBoundary.find(node) == mesh.m_nodeBoundary.end()) {
					mesh.interiorSmoothing(node);
				}
			}

			for(auto &eb : listEdge_bis){
				mesh.triangleSmoothing(eb);
			}

			listBoundary.push_back(listEdge_bis);
		}

		mesh.m_listFront = listBoundary;
		pr += 1;

		mesh.initialise_Boundary();

		gmds::IGMeshIOService ioService(mesh.m_mesh);
		gmds::VTKWriter vtkWriter(&ioService);
		vtkWriter.setCellOptions(gmds::N | gmds::F | gmds::E);
		vtkWriter.setDataOptions(gmds::N | gmds::F | gmds::E);
		vtkWriter.write("/home/pagea/Documents/Travail/data/toto_link3.vtk");
	}


	for (auto k = 0 ; k != 10 ; k++) {
		for (auto &i : mesh.m_mesh->nodes()) {
			Node node = mesh.m_mesh->get<Node>(i);
			if (mesh.m_nodeBoundary.find(node) == mesh.m_nodeBoundary.end() and mesh.m_nodeFront.find(node) == mesh.m_nodeFront.end()) {
				mesh.interiorSmoothing(node);
			}
		}
	}

	std::cout<<"\nsuivant\n"<<std::endl;
	for(auto& n : mesh.m_nodeFront) {
		Node node = n.first;
		std::cout << node << std::endl;
	}

	std::cout<<"\nseam\n"<<std::endl;
	std::vector<Node> listNode;

	for(auto& n : mesh.m_nodeFront){
		std::cout<<"listNode : "<<listNode.size()<<std::endl;
		if (std::find(listNode.begin(), listNode.end(), n.first) == listNode.end()){
			Node node = n.first;
			std::cout << node << std::endl;
			mesh.seam(node, listNode);
		}
		gmds::IGMeshIOService ioService(mesh.m_mesh);
		gmds::VTKWriter vtkWriter(&ioService);
		vtkWriter.setCellOptions(gmds::N | gmds::F | gmds::E);
		vtkWriter.setDataOptions(gmds::N | gmds::F | gmds::E);
		vtkWriter.write("/home/pagea/Documents/Travail/data/toto_link3.vtk");
	}
}

TEST(DummyTestClass, stat){
	quadfront::Quadfront mesh("/home/pagea/Documents/Travail/data/petitangle.vtk");


	Node node = mesh.m_mesh->get<Node>(3);
	std::vector<Face> listFace;
	std::vector<Edge> listEdge;
	mesh.petitAngle(node.get<Face>()[0], listFace, listEdge);
	for(auto& f : listFace){
		std::cout<<f<<std::endl;
	}
	std::cout<<"\n"<<std::endl;
	for(auto& e : listEdge){
		std::cout<<e.id()<<" : "<<e.get<Node>()[0]<<" ; "<<e.get<Node>()[1]<<std::endl;
	}



	gmds::IGMeshIOService ioService(mesh.m_mesh);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N | gmds::F | gmds::E);
	vtkWriter.setDataOptions(gmds::N | gmds::F | gmds::E);
	vtkWriter.write("/home/pagea/Documents/Travail/data/toto_link3.vtk");
}

TEST(DummyTestClass, petitangle){
	quadfront::Quadfront mesh("/home/pagea/Documents/Travail/data/petitangle.vtk");

	mesh.initialise_Boundary();
	std::vector<Face> listQuad;

	for(auto& i : mesh.m_nodeBoundary){
		Node node = i.first;
		mesh.initialiseQuad(node, listQuad);
	}


	gmds::IGMeshIOService ioService(mesh.m_mesh);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N | gmds::F | gmds::E);
	vtkWriter.setDataOptions(gmds::N | gmds::F | gmds::E);
	vtkWriter.write("/home/pagea/Documents/Travail/data/toto_link3.vtk");
}

TEST(DummyTestClass, revers){
	quadfront::Quadfront mesh("/home/pagea/Documents/Travail/data/petitangle.vtk");


	std::vector<Face> listQuad;

	for(auto& i : mesh.m_nodeBoundary) {
		Node node = i.first;
		std::cout<<node<<std::endl;
		mesh.initialiseQuad(node, listQuad);
		gmds::IGMeshIOService ioService(mesh.m_mesh);
		gmds::VTKWriter vtkWriter(&ioService);
		vtkWriter.setCellOptions(gmds::N | gmds::F | gmds::E);
		vtkWriter.setDataOptions(gmds::N | gmds::F | gmds::E);
		vtkWriter.write("/home/pagea/Documents/Travail/data/toto_link3.vtk");
	}

	std::cout<<listQuad.size()<<std::endl;
}
