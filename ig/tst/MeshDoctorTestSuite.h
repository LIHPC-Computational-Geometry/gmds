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

/*----------------------------------------------------------------------------*/
