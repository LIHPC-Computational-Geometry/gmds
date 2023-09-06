/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
TEST(MeshDocClass, buildN2F)
{
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::F2N | gmds::N2F));

	gmds::Node n0 = m.newNode(0,0,0);
    gmds::Node n1 = m.newNode(1,1,0);
    gmds::Node n2 = m.newNode(0,1,0);
    gmds::Node n3 = m.newNode(1,0,0);

    m.newTriangle(n0,n1,n3);
    m.newTriangle(n0,n3,n2);
    m.newQuad(n0,n1,n2,n3);


    ASSERT_EQ(n0.get<gmds::Face>().size(),0);
    ASSERT_EQ(n1.get<gmds::Face>().size(),0);
    ASSERT_EQ(n2.get<gmds::Face>().size(),0);
    ASSERT_EQ(n3.get<gmds::Face>().size(),0);

    gmds::MeshDoctor doc(&m);

    doc.updateUpwardConnectivity();


    ASSERT_EQ(n0.get<gmds::Face>().size(),3);
    ASSERT_EQ(n1.get<gmds::Face>().size(),2);
    ASSERT_EQ(n2.get<gmds::Face>().size(),2);
    ASSERT_EQ(n3.get<gmds::Face>().size(),3);
}

TEST(MeshDocClass, buildE)
{
	gmds::Mesh m(gmds::MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                             F2E | E2F | R2E | E2R | N2R | N2F | N2E));

	gmds::Node n0 = m.newNode(0,0,1);
	gmds::Node n1 = m.newNode(1,0,1);
	gmds::Node n2 = m.newNode(1,1,1);
	gmds::Node n3 = m.newNode(0,1,1);
	gmds::Node n4 = m.newNode(0,0,1);
	gmds::Node n5 = m.newNode(1,0,1);
	gmds::Node n6 = m.newNode(1,1,1);
	gmds::Node n7 = m.newNode(0,1,1);

	m.newTriangle(n0,n1,n2);
	m.newTriangle(n0,n2,n3);
	m.newQuad(n0,n1,n2,n3);

	ASSERT_EQ(m.getNbEdges(), 0);
	gmds::MeshDoctor doc(&m);
	doc.buildE();
	ASSERT_EQ(m.getNbEdges(), 0);
	m.newEdge(n0,n7);
	ASSERT_EQ(m.getNbEdges(), 1);
	doc.buildEfromR();
	ASSERT_EQ(m.getNbEdges(), 1);
	doc.buildEfromF();
	ASSERT_EQ(m.getNbEdges(), 6);

	m.newHex(n0,n1,n2,n3,n4,n5,n6,n7);
	m.newPyramid(n0,n3,n7,n4,n1);
	m.newPrism3(n0,n1,n2,n4,n5,n6);
	doc.buildEfromR();
	ASSERT_EQ(m.getNbEdges(), 18);

}
/*----------------------------------------------------------------------------*/
