/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/io/LimaWriter.h>
#include <gmds/io/LimaReader.h>
#include <gtest/gtest.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(LimaWriterTestSuite, testLimaWriter1)
{
    // WE WRITE
    Mesh m(MeshModel(DIM3|F|N|F2N));

    Node n0 = m.newNode(0,0,0);
    Node n1 = m.newNode(1,0,0);
    Node n2 = m.newNode(1,1,0);
    Node n3 = m.newNode(0,1,0);
    Node n4 = m.newNode(2,0,0);

    Face f0 = m.newTriangle(n1,n2,n4);
    Face f1 = m.newQuad(n0,n1,n2,n3);

	 CellGroup<Node>* node_grp = m.newGroup<Node>("surface1");
	 node_grp->add(n1);
	 node_grp->add(n2);
	 node_grp->add(n3);
	 CellGroup<Face>* face_grp = m.newGroup<Face>("surface1");
	 face_grp->add(f0);

	 LimaWriter writer(m);
    writer.write("test_lima.mli2", F|N);

    // WE READ
    Mesh m2(MeshModel(DIM3|F|N|F2N));
	 LimaReader reader(m2);
    reader.read("test_lima.mli2",F|N);

    ASSERT_EQ(m2.getNbNodes(),5);
    ASSERT_EQ(m2.getNbFaces(),2);
	 ASSERT_EQ(m2.getNbGroups<Face>(),1);
	 CellGroup<Node>* node_grp2 = m2.getGroup<Node>(0);
	 ASSERT_EQ(node_grp2->name(),"surface1");
	 CellGroup<Face>* face_grp2 = m2.getGroup<Face>(0);
	 ASSERT_EQ(face_grp2->name(),"surface1");
}