/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
TEST(BooleanMarkMesh, testNodeMarks)
{
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::F2N));

	gmds::Node n0 = m.newNode(0,0,0);
	gmds::Node n1 = m.newNode(1,1,1);
	gmds::Node n2 = m.newNode(0,1,1);
	gmds::Node n3 = m.newNode(1,0,1);

	m.newTriangle(n0,n1,n3);
	m.newTriangle(n0,n3,n2);
	m.newQuad(n0,n1,n2,n3);

	int mark = m.newMark<gmds::Node>();

	ASSERT_EQ(gmds::Marks32(), m.getMarks(n0));
	ASSERT_EQ(gmds::Marks32(), m.getMarks(n1));
	ASSERT_EQ(gmds::Marks32(), m.getMarks(n2));
	ASSERT_EQ(gmds::Marks32(), m.getMarks(n3));

	int nb_marked_nodes=0;
	for(auto i : m.nodes()){
		if(m.isMarked<gmds::Node>(i,mark))
			nb_marked_nodes++;
	}
	ASSERT_EQ(0,nb_marked_nodes);

	m.mark(n0,mark);
	m.mark<gmds::Node>(n1.id(),mark);
	m.mark<gmds::Node>(n2.id(),mark);
	nb_marked_nodes=0;
	for(auto i : m.nodes()){
		if(m.isMarked(m.get<gmds::Node>(i),mark))
			nb_marked_nodes++;
	}
	ASSERT_EQ(3,nb_marked_nodes);

	gmds::Marks32 marks32;
	marks32.set(mark);
	ASSERT_EQ(marks32, m.getMarks(n0));
	ASSERT_EQ(gmds::Marks32(), m.getMarks(n3));

	m.negateMaskMark<gmds::Node>(mark);

	nb_marked_nodes=0;
	for(auto i : m.nodes()){
		if(m.isMarked(m.get<gmds::Node>(i),mark))
			nb_marked_nodes++;
	}
	ASSERT_EQ(1,nb_marked_nodes);

	ASSERT_EQ(marks32, m.getMarks(n0));
	ASSERT_EQ(gmds::Marks32(), m.getMarks(n3));

	m.unmarkAll<gmds::Node>(mark);

	nb_marked_nodes=0;
	for(auto i : m.nodes()){
		if(m.isMarked(m.get<gmds::Node>(i),mark))
			nb_marked_nodes++;
	}

	ASSERT_EQ(0,nb_marked_nodes);

	m.freeMark<gmds::Node>(mark);
}
/*----------------------------------------------------------------------------*/
TEST(BooleanMarkMesh, testFaceMarks)
{
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::F2N));

	gmds::Node n0 = m.newNode(0,0,0);
	gmds::Node n1 = m.newNode(1,1,1);
	gmds::Node n2 = m.newNode(0,1,1);
	gmds::Node n3 = m.newNode(1,0,1);

	gmds::Face f1 = m.newTriangle(n0,n1,n3);
	gmds::Face f2 = m.newTriangle(n0,n3,n2);
	gmds::Face f3 = m.newQuad(n0,n1,n2,n3);

	int mark = m.newMark<gmds::Face>();

	ASSERT_EQ(gmds::Marks32(), m.getMarks(f1));
	ASSERT_EQ(gmds::Marks32(), m.getMarks(f2));
	ASSERT_EQ(gmds::Marks32(), m.getMarks(f3));

	int nb_marked_faces=0;
	for(auto i : m.faces()){
		if(m.isMarked<gmds::Face>(i,mark))
			nb_marked_faces++;
	}
	ASSERT_EQ(0,nb_marked_faces);

	m.mark(f1,mark);
	m.mark<gmds::Face>(f2.id(),mark);

	nb_marked_faces=0;
	for(auto i : m.faces()){
		if(m.isMarked(m.get<gmds::Face>(i),mark))
			nb_marked_faces++;
	}
	ASSERT_EQ(2,nb_marked_faces);

	gmds::Marks32 marks32;
	marks32.set(mark);
	ASSERT_EQ(marks32, m.getMarks(f1));
	ASSERT_EQ(gmds::Marks32(), m.getMarks(f3));

	m.negateMaskMark<gmds::Face>(mark);

	nb_marked_faces=0;
	for(auto i : m.faces()){
		if(m.isMarked(m.get<gmds::Face>(i),mark))
			nb_marked_faces++;
	}
	ASSERT_EQ(1,nb_marked_faces);

	ASSERT_EQ(marks32, m.getMarks(f1));
	ASSERT_EQ(gmds::Marks32(), m.getMarks(f3));

	m.unmarkAll<gmds::Face>(mark);

	nb_marked_faces=0;
	for(auto i : m.faces()){
		if(m.isMarked(m.get<gmds::Face>(i),mark))
			nb_marked_faces++;
	}

	ASSERT_EQ(0,nb_marked_faces);

	m.freeMark<gmds::Face>(mark);
}

/*----------------------------------------------------------------------------*/
TEST(BooleanMarkMesh, testRegionMarks)
{
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::R|gmds::N|gmds::R2N));

	gmds::Node n0 = m.newNode(0, 0, 0);
	gmds::Node n1 = m.newNode(0, 1, 0);
	gmds::Node n2 = m.newNode(1, 1, 0);
	gmds::Node n3 = m.newNode(0, 1, 0);

	gmds::Node n4 = m.newNode(0, 0, 1);
	gmds::Node n5 = m.newNode(0, 1, 1);
	gmds::Node n6 = m.newNode(1, 1, 1);
	gmds::Node n7 = m.newNode(0, 1, 1);


	gmds::Region t0 =m.newTet(n0,n1,n3,n4);
	gmds::Region t1 =m.newTet(n1,n2,n3,n6);
	gmds::Region t2 =m.newTet(n4,n5,n6,n1);

	int mark = m.newMark<gmds::Region>();

	ASSERT_EQ(gmds::Marks32(), m.getMarks(t0));
	ASSERT_EQ(gmds::Marks32(), m.getMarks(t1));
	ASSERT_EQ(gmds::Marks32(), m.getMarks(t2));

	int nb_marked=0;
	for(auto i : m.regions()){
		if(m.isMarked<gmds::Region>(i,mark))
			nb_marked++;
	}
	ASSERT_EQ(0,nb_marked);

	m.mark(t0,mark);
	m.mark<gmds::Region>(t1.id(),mark);

	nb_marked=0;
	for(auto i : m.regions()){
		if(m.isMarked(m.get<gmds::Region>(i),mark))
			nb_marked++;
	}
	ASSERT_EQ(2,nb_marked);

	gmds::Marks32 marks32;
	marks32.set(mark);
	ASSERT_EQ(marks32, m.getMarks(t0));
	ASSERT_EQ(gmds::Marks32(), m.getMarks(t2));

	m.negateMaskMark<gmds::Region>(mark);

	nb_marked=0;
	for(auto i : m.regions()){
		if(m.isMarked(m.get<gmds::Region>(i),mark))
			nb_marked++;
	}
	ASSERT_EQ(1,nb_marked);

	ASSERT_EQ(marks32, m.getMarks(t0));
	ASSERT_EQ(gmds::Marks32(), m.getMarks(t2));

	m.unmarkAll<gmds::Region>(mark);

	nb_marked=0;
	for(auto i : m.regions()){
		if(m.isMarked(m.get<gmds::Region>(i),mark))
			nb_marked++;
	}

	ASSERT_EQ(0,nb_marked);

	m.freeMark<gmds::Region>(mark);
}

/*----------------------------------------------------------------------------*/
TEST(BooleanMarkMesh, testEdgeMarks)
{
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|
										  gmds::N|gmds::E|
										  gmds::F2N|gmds::E2N|
										  gmds::F2E));

	gmds::Node n0 = m.newNode(0,0,0);
	gmds::Node n1 = m.newNode(1,0,0);
	gmds::Node n2 = m.newNode(0,1,0);
	gmds::Node n3 = m.newNode(1,1,0);
	gmds::Node n4 = m.newNode(2,0,0);
	gmds::Node n5 = m.newNode(2,1,0);

	m.newTriangle(n0,n1,n3);
	m.newTriangle(n1,n3,n2);
	m.newQuad(n1,n4,n5,n3);

	gmds::MeshDoctor doc(&m);
	doc.buildE();
	ASSERT_EQ(8,m.getNbEdges());
	int mark = m.newMark<gmds::Edge>();

	int nb_marked_edges=0;
	for(auto i : m.edges()){
		if(m.isMarked<gmds::Edge>(i,mark)) {
			nb_marked_edges++;
		}
		ASSERT_EQ(gmds::Marks32(), m.getMarks(m.get<gmds::Edge>(i)));
	}
	ASSERT_EQ(0,nb_marked_edges);


	m.mark(m.get<gmds::Edge>(1),mark);
	m.mark<gmds::Edge>(2,mark);

	nb_marked_edges=0;
	for(auto i : m.edges()){
		if(m.isMarked(m.get<gmds::Edge>(i),mark))
			nb_marked_edges++;
	}
	ASSERT_EQ(2,nb_marked_edges);

	gmds::Marks32 marks32;
	marks32.set(mark);
	ASSERT_EQ(marks32, m.getMarks(m.get<gmds::Edge>(1)));
	ASSERT_EQ(gmds::Marks32(), m.getMarks(m.get<gmds::Edge>(3)));

	m.negateMaskMark<gmds::Edge>(mark);

	nb_marked_edges=0;
	for(auto i : m.edges()){
		if(m.isMarked(m.get<gmds::Edge>(i),mark))
			nb_marked_edges++;
	}
	ASSERT_EQ(6,nb_marked_edges);

	ASSERT_EQ(marks32, m.getMarks(m.get<gmds::Edge>(1)));
	ASSERT_EQ(gmds::Marks32(), m.getMarks(m.get<gmds::Edge>(3)));

	m.unmarkAll<gmds::Edge>(mark);

	nb_marked_edges=0;
	for(auto i : m.edges()){
		if(m.isMarked(m.get<gmds::Edge>(i),mark))
			nb_marked_edges++;
	}

	ASSERT_EQ(0,nb_marked_edges);

	m.freeMark<gmds::Edge>(mark);
}

/*----------------------------------------------------------------------------*/
TEST(BooleanMarkMesh, testTooMuchNodeMarks){
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|
										  gmds::N|gmds::E|
										  gmds::F2N|gmds::E2N));
	EXPECT_ANY_THROW(
		int nb=32;
		for(auto i=0; i<nb;i++){
			m.newMark<gmds::Node>();
		}
	);

}

/*----------------------------------------------------------------------------*/
TEST(BooleanMarkMesh, testTooMuchEdgeMarks){
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|
										  gmds::N|gmds::E|
										  gmds::F2N|gmds::E2N));
	EXPECT_ANY_THROW(
		int nb=32;
		for(auto i=0; i<nb;i++){
			m.newMark<gmds::Edge>();
		}
	);

}

/*----------------------------------------------------------------------------*/
TEST(BooleanMarkMesh, testTooMuchFaceMarks){
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|
										  gmds::N|gmds::E|
										  gmds::F2N|gmds::E2N));
	EXPECT_ANY_THROW(
		int nb=32;
		for(auto i=0; i<nb;i++){
			m.newMark<gmds::Face>();
		}
	);

}

/*----------------------------------------------------------------------------*/
TEST(BooleanMarkMesh, testTooMuchRegionMarks){
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|
										  gmds::N|gmds::E|
										  gmds::F2N|gmds::E2N));
	EXPECT_ANY_THROW(
		int nb=32;
		for(auto i=0; i<nb;i++){
			m.newMark<gmds::Region>();
		}
	);

}
/*----------------------------------------------------------------------------*/
TEST(BooleanMarkMesh, testTDoubleFreeMarkNode){
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::R|gmds::F|
	                             gmds::N|gmds::E|gmds::R2N|
	                             gmds::F2N|gmds::E2N));
	int mark1 = m.newMark<gmds::Node>();
	m.freeMark<gmds::Node>(mark1);
	m.freeMark<gmds::Node>(mark1);
	int mark2 = m.newMark<gmds::Node>();
	int mark3 = m.newMark<gmds::Node>();
	ASSERT_EQ(mark2,0);
	ASSERT_EQ(mark3,1);
}
/*----------------------------------------------------------------------------*/
TEST(BooleanMarkMesh, testTDoubleFreeMarkEdge){
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::R|gmds::F|
	                             gmds::N|gmds::E|gmds::R2N|
	                             gmds::F2N|gmds::E2N));
	int mark1 = m.newMark<gmds::Edge>();
	m.freeMark<gmds::Edge>(mark1);
	m.freeMark<gmds::Edge>(mark1);
	int mark2 = m.newMark<gmds::Edge>();
	int mark3 = m.newMark<gmds::Edge>();
	ASSERT_EQ(mark2,0);
	ASSERT_EQ(mark3,1);
}
/*----------------------------------------------------------------------------*/
TEST(BooleanMarkMesh, testTDoubleFreeMarkFace){
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::R|gmds::F|
	                             gmds::N|gmds::E|gmds::R2N|
	                             gmds::F2N|gmds::E2N));
	int mark1 = m.newMark<gmds::Face>();
	m.freeMark<gmds::Face>(mark1);
	m.freeMark<gmds::Face>(mark1);
	int mark2 = m.newMark<gmds::Face>();
	int mark3 = m.newMark<gmds::Face>();
	ASSERT_EQ(mark2,0);
	ASSERT_EQ(mark3,1);
}
/*----------------------------------------------------------------------------*/
TEST(BooleanMarkMesh, testTDoubleFreeMarkRegion){
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::R|gmds::F|
	                             gmds::N|gmds::E|gmds::R2N|
	                             gmds::F2N|gmds::E2N));
	int mark1 = m.newMark<gmds::Region>();
	m.freeMark<gmds::Region>(mark1);
	m.freeMark<gmds::Region>(mark1);
	int mark2 = m.newMark<gmds::Region>();
	int mark3 = m.newMark<gmds::Region>();
	ASSERT_EQ(mark2,0);
	ASSERT_EQ(mark3,1);
}