/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
TEST(EdgeTestSuite, testMeshEdge)
{
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::E | gmds::N | gmds::E2N));
	gmds::Node n0 = m.newNode(0, 0, 0);
	gmds::Node n1 = m.newNode(1, 0, 0);
	gmds::Node n2 = m.newNode(0, 1, 0);
	gmds::Edge e0 = m.newEdge(n0,n1);
	gmds::Edge e1 = m.newEdge(n0,n2);

	gmds::Edge ea = m.get<gmds::Edge> (e0.id());

	ASSERT_TRUE(e0 == ea);
	ASSERT_FALSE(e0 == e1);

	ASSERT_FALSE(e0 != ea);
	ASSERT_TRUE(e0 != e1);
}
/*----------------------------------------------------------------------------*/
TEST(EdgeTestSuite, testAccessorIssues) {

	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E));
	gmds::Node n0 = m.newNode(0, 0, 0);
	gmds::Node n1 = m.newNode(1, 0, 0);
	gmds::Node n2 = m.newNode(2, 0, 0);
	gmds::Edge e0 = m.newEdge(n0,n1);

//	ASSERT_ANY_THROW(e0.getIDs<gmds::Node>());
//	ASSERT_ANY_THROW(e0.getIDs<gmds::Edge>());
//	ASSERT_ANY_THROW(e0.getIDs<gmds::Face>());
//	ASSERT_ANY_THROW(e0.getIDs<gmds::Region>());
//   ASSERT_ANY_THROW(e0.getAllIDs<gmds::Node>());
//	ASSERT_ANY_THROW(e0.getAllIDs<gmds::Edge>());
//	ASSERT_ANY_THROW(e0.getAllIDs<gmds::Face>());
//	ASSERT_ANY_THROW(e0.getAllIDs<gmds::Region>());
//	ASSERT_ANY_THROW(e0.get<gmds::Node>());
//	ASSERT_ANY_THROW(e0.get<gmds::Edge>());
//	ASSERT_ANY_THROW(e0.get<gmds::Face>());
//	ASSERT_ANY_THROW(e0.get<gmds::Region>());
//	ASSERT_ANY_THROW(e0.getAll<gmds::Node>());
//	ASSERT_ANY_THROW(e0.getAll<gmds::Edge>());
//	ASSERT_ANY_THROW(e0.getAll<gmds::Face>());
//	ASSERT_ANY_THROW(e0.getAll<gmds::Region>());

	gmds::MeshModel mod(gmds::DIM2 | gmds::F | gmds::N | gmds::E | gmds::E2N);
	m.changeModel(mod);

	std::vector<gmds::Node> nodes = {n0, n1};
	e0.set<gmds::Node> (nodes);

	std::vector<gmds::TCellID> ids = e0.getIDs<gmds::Node>();
	ASSERT_EQ(ids.size(), 2);
	auto it = std::find(ids.begin(), ids.end(), n0.id());
	ASSERT_TRUE(it != ids.end());
	it = std::find(ids.begin(), ids.end(), n1.id());
	ASSERT_TRUE(it != ids.end());
	it = std::find(ids.begin(), ids.end(), n2.id());
	ASSERT_TRUE(it == ids.end());

	e0.replace<gmds::Node>(n0,n2);
	ids = e0.getIDs<gmds::Node>();
	it = std::find(ids.begin(), ids.end(), n0.id());
	ASSERT_TRUE(it == ids.end());
	it = std::find(ids.begin(), ids.end(), n1.id());
	ASSERT_TRUE(it != ids.end());
	it = std::find(ids.begin(), ids.end(), n2.id());
	ASSERT_TRUE(it != ids.end());

	e0.remove(n1);
	ids = e0.getIDs<gmds::Node>();
	ASSERT_EQ(ids.size(), 1);
}
/*----------------------------------------------------------------------------*/
TEST(EdgeTestSuite, testAccessorAll) {
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::N | gmds::F | gmds::E |
	                             gmds::F2N | gmds::R2N| gmds::E2N |
	                             gmds::R2F | gmds::F2R | gmds::R2E ));

	gmds::Node n0 = m.newNode(0, 0, 0);
	gmds::Node n1 = m.newNode(0, 1, 0);
	gmds::Node n2 = m.newNode(1, 1, 0);
	gmds::Node n3 = m.newNode(0, 1, 1);
	gmds::Node n4 = m.newNode(0, 1, 2);


	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.buildFacesAndR2F();
	doc.updateUpwardConnectivity();

	gmds::Edge e0 = m.get<gmds::Edge> (0);
//	ASSERT_EQ(e0.nbFaces(), 1);
}
/*----------------------------------------------------------------------------*/