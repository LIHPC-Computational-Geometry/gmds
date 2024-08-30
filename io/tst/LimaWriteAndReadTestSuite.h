/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/io/LimaReader.h>
#include <gmds/io/LimaWriter.h>
#include <gmds/io/LimaWriterStreamMli2.h>
#include <gtest/gtest.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(LimaWriterTestSuite, testLimaWriterStreamMli2)
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

	 LimaWriterStreamMli2 writer(m);
    writer.write("test_lima0.mli2", F|N);

    // WE READ
    Mesh m2(MeshModel(DIM3|F|N|F2N));
	 LimaReader reader(m2);
    reader.read("test_lima0.mli2",F|N);

    ASSERT_EQ(m2.getNbNodes(),5);
    ASSERT_EQ(m2.getNbFaces(),2);
	 ASSERT_EQ(m2.getNbGroups<Face>(),1);
	 CellGroup<Node>* node_grp2 = m2.getGroup<Node>(0);
	 ASSERT_EQ(node_grp2->name(),"surface1");
	 CellGroup<Face>* face_grp2 = m2.getGroup<Face>(0);
	 ASSERT_EQ(face_grp2->name(),"surface1");
}
/*----------------------------------------------------------------------------*/
TEST(LimaWriterTestSuite, testLimaWriterStreamMli2_1)
{
	 // WE WRITE
	 Mesh m(MeshModel(DIM3|R|F|E|N|R2N|F2N|E2N));

	 Node n0 = m.newNode(0,0,0);
	 Node n1 = m.newNode(1,0,0);
	 Node n2 = m.newNode(1,1,0);
	 Node n3 = m.newNode(0,1,0);
	 Node n4 = m.newNode(2,0,0);

	 Edge e0 = m.newEdge(n1, n2);
	 Edge e1 = m.newEdge(n0, n3);

	 Face f0 = m.newTriangle(n1,n2,n4);
	 Face f1 = m.newQuad(n0,n1,n2,n3);

	 Region r0 = m.newTet(n0,n1,n2,n3);
	 Region r1 = m.newTet(n0,n1,n2,n4);


	 CellGroup<Node>* node_grp = m.newGroup<Node>("cloud1");
	 node_grp->add(n1);
	 node_grp->add(n2);
	 node_grp->add(n3);
	 CellGroup<Edge>* edge_grp = m.newGroup<Edge>("line1");
	 edge_grp->add(e1);
	 CellGroup<Face>* face_grp = m.newGroup<Face>("surface1");
	 face_grp->add(f0);
	 CellGroup<Region>* region_grp = m.newGroup<Region>("vol1");
	 region_grp->add(r0);

	 LimaWriterStreamMli2 writer(m);
	 writer.write("test_lima1.mli2", R|F|E|N);

	 // WE READ
	 Mesh m2(MeshModel(DIM3|R|F|E|N|R2N|F2N|E2N));
	 LimaReader reader(m2);
	 reader.read("test_lima1.mli2",R|F|E|N);

	 ASSERT_EQ(m2.getNbNodes(),5);
	 ASSERT_EQ(m2.getNbEdges(),2);
	 ASSERT_EQ(m2.getNbFaces(),2);
	 ASSERT_EQ(m2.getNbRegions(),2);
	 ASSERT_EQ(m2.getNbGroups<Face>(),1);
	 CellGroup<Node>* node_grp2 = m2.getGroup<Node>(0);
	 ASSERT_EQ(node_grp2->name(),"cloud1");
	 CellGroup<Edge>* edge_grp2 = m2.getGroup<Edge>(0);
	 ASSERT_EQ(edge_grp2->name(),"line1");
	 CellGroup<Face>* face_grp2 = m2.getGroup<Face>(0);
	 ASSERT_EQ(face_grp2->name(),"surface1");
	 CellGroup<Region>* region_grp2 = m2.getGroup<Region>(0);
	 ASSERT_EQ(region_grp2->name(),"vol1");
}
/*----------------------------------------------------------------------------*/
TEST(LimaWriterTestSuite, testLimaWriterUnf)
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
	 writer.write("test_lima1.unf", F|N);

	 // WE READ
	 Mesh m2(MeshModel(DIM3|F|N|F2N));
	 LimaReader reader(m2);
	 reader.read("test_lima1.unf",F|N);

	 ASSERT_EQ(m2.getNbNodes(),5);
	 ASSERT_EQ(m2.getNbFaces(),2);
	 ASSERT_EQ(m2.getNbGroups<Face>(),1);
	 CellGroup<Node>* node_grp2 = m2.getGroup<Node>(0);
	 ASSERT_EQ(node_grp2->name(),"surface1");
	 CellGroup<Face>* face_grp2 = m2.getGroup<Face>(0);
	 ASSERT_EQ(face_grp2->name(),"surface1");
}
/*----------------------------------------------------------------------------*/
TEST(LimaWriterTestSuite, testLimaWriterUnf_1)
{
	 // WE WRITE
	 Mesh m(MeshModel(DIM3|R|F|E|N|R2N|F2N|E2N));

	 Node n0 = m.newNode(0,0,0);
	 Node n1 = m.newNode(1,0,0);
	 Node n2 = m.newNode(1,1,0);
	 Node n3 = m.newNode(0,1,0);
	 Node n4 = m.newNode(2,0,0);

	 Edge e0 = m.newEdge(n1, n2);
	 Edge e1 = m.newEdge(n0, n3);

	 Face f0 = m.newTriangle(n1,n2,n4);
	 Face f1 = m.newQuad(n0,n1,n2,n3);

	 Region r0 = m.newTet(n0,n1,n2,n3);
	 Region r1 = m.newTet(n0,n1,n2,n4);


	 CellGroup<Node>* node_grp = m.newGroup<Node>("cloud1");
	 node_grp->add(n1);
	 node_grp->add(n2);
	 node_grp->add(n3);
	 CellGroup<Edge>* edge_grp = m.newGroup<Edge>("line1");
	 edge_grp->add(e1);
	 CellGroup<Face>* face_grp = m.newGroup<Face>("surface1");
	 face_grp->add(f0);
	 CellGroup<Region>* region_grp = m.newGroup<Region>("vol1");
	 region_grp->add(r0);

	 LimaWriter writer(m);
	 writer.write("test_lima2.unf", R|F|E|N);

	 // WE READ
	 Mesh m2(MeshModel(DIM3|R|F|E|N|R2N|F2N|E2N));
	 LimaReader reader(m2);
	 reader.read("test_lima2.unf",R|F|E|N);

	 ASSERT_EQ(m2.getNbNodes(),5);
	 ASSERT_EQ(m2.getNbEdges(),2);
	 ASSERT_EQ(m2.getNbFaces(),2);
	 ASSERT_EQ(m2.getNbRegions(),2);
	 ASSERT_EQ(m2.getNbGroups<Face>(),1);
	 CellGroup<Node>* node_grp2 = m2.getGroup<Node>(0);
	 ASSERT_EQ(node_grp2->name(),"cloud1");
	 CellGroup<Edge>* edge_grp2 = m2.getGroup<Edge>(0);
	 ASSERT_EQ(edge_grp2->name(),"line1");
	 CellGroup<Face>* face_grp2 = m2.getGroup<Face>(0);
	 ASSERT_EQ(face_grp2->name(),"surface1");
	 CellGroup<Region>* region_grp2 = m2.getGroup<Region>(0);
	 ASSERT_EQ(region_grp2->name(),"vol1");
}
/*----------------------------------------------------------------------------*/
TEST(LimaWriterTestSuite, testLimaWriter)
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
	 writer.write("test_lima2.mli2", F|N);

	 // WE READ
	 Mesh m2(MeshModel(DIM3|F|N|F2N));
	 LimaReader reader(m2);
	 reader.read("test_lima2.mli2",F|N);

	 ASSERT_EQ(m2.getNbNodes(),5);
	 ASSERT_EQ(m2.getNbFaces(),2);
	 ASSERT_EQ(m2.getNbGroups<Face>(),1);
	 CellGroup<Node>* node_grp2 = m2.getGroup<Node>(0);
	 ASSERT_EQ(node_grp2->name(),"surface1");
	 CellGroup<Face>* face_grp2 = m2.getGroup<Face>(0);
	 ASSERT_EQ(face_grp2->name(),"surface1");
}
/*----------------------------------------------------------------------------*/