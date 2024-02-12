/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/EdgeContainer.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/Node.h>
#include <gmds/ig/Edge.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
TEST(EdgeContainer, testEdgeContainerConstructor)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::E));
    gmds::EdgeContainer edgeContainer(&mesh);

    ASSERT_TRUE(edgeContainer.has(0));
    ASSERT_FALSE(edgeContainer.has(1));
}
/*----------------------------------------------------------------------------*/
TEST(EdgeContainer, testAddAndGetNodesData)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E));
    gmds::Node n1 = mesh.newNode(0, 0, 0);
    gmds::Node n2 = mesh.newNode(1, 0, 0);

    gmds::EdgeContainer edgeContainer(&mesh);
    gmds::Edge edge = edgeContainer.add(n1.id(), n2.id());

    TInt nbNodes;
    edgeContainer.getNodesData(edge.id(), nbNodes);

    ASSERT_EQ(2, nbNodes);
}
/*----------------------------------------------------------------------------*/
TEST(EdgeContainer, testAddAndGetEdgesData)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E));
    gmds::Node n1 = mesh.newNode(0, 0, 0);
    gmds::Node n2 = mesh.newNode(1, 0, 0);
    gmds::Node n3 = mesh.newNode(1, 1, 0);

    gmds::EdgeContainer edgeContainer(&mesh);
    gmds::Edge e1 = mesh.newEdge(n1, n2);
    gmds::Edge e2 = mesh.newEdge(n2, n3);
    gmds::Edge e3 = mesh.newEdge(n1, n3);

    gmds::Edge edge = edgeContainer.add(n1.id(), n2.id());
    edgeContainer.add(n2.id(), n3.id());
    edgeContainer.add(n1.id(), n3.id());

    TInt nbEdges;
    edgeContainer.getEdgesData(edge.id(), nbEdges);

    ASSERT_EQ(2, nbEdges);
}
/*----------------------------------------------------------------------------*/
TEST(EdgeContainer, testAddAndGetFacesData)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E | gmds::F));
    gmds::Node n1 = mesh.newNode(0, 0, 0);
    gmds::Node n2 = mesh.newNode(1, 0, 0);
    gmds::Node n3 = mesh.newNode(1, 1, 0);

    gmds::FaceContainer faceContainer(&mesh);
    gmds::Face f1 = mesh.newFace(n1, n2, n3);

    gmds::EdgeContainer edgeContainer(&mesh);
    gmds::Edge edge = edgeContainer.add(n1.id(), n2.id());

    TInt nbFaces;
    edgeContainer.getFacesData(edge.id(), nbFaces);

    ASSERT_EQ(1, nbFaces);
}
/*----------------------------------------------------------------------------*/
TEST(EdgeContainer, testAddAndGetRegionsData)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E | gmds::R));
    gmds::Node n1 = mesh.newNode(0, 0, 0);
    gmds::Node n2 = mesh.newNode(1, 0, 0);
    gmds::Node n3 = mesh.newNode(1, 1, 0);

    gmds::RegionContainer regionContainer(&mesh);
    gmds::Region r1 = mesh.newRegion(n1, n2, n3);

    gmds::EdgeContainer edgeContainer(&mesh);
    gmds::Edge edge = edgeContainer.add(n1.id(), n2.id());

    TInt nbRegions;
    edgeContainer.getRegionsData(edge.id(), nbRegions);

    ASSERT_EQ(1, nbRegions);
}
/*----------------------------------------------------------------------------*/
TEST(EdgeContainer, testAddAndBuildEdge)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E));
    gmds::Node n1 = mesh.newNode(0, 0, 0);
    gmds::Node n2 = mesh.newNode(1, 0, 0);

    gmds::EdgeContainer edgeContainer(&mesh);
    gmds::Edge edge = edgeContainer.add(n1.id(), n2.id());

    gmds::Edge builtEdge = edgeContainer.buildEdge(edge.id());

    ASSERT_EQ(edge.id(), builtEdge.id());
}
/*----------------------------------------------------------------------------*/
TEST(EdgeContainer, testHas)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E));
    gmds::EdgeContainer edgeContainer(&mesh);

    ASSERT_FALSE(edgeContainer.has(0));

    gmds::Edge edge = edgeContainer.add(1, 2);
    ASSERT_TRUE(edgeContainer.has(edge.id()));
}
/*----------------------------------------------------------------------------*/
TEST(EdgeContainer, testClear)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E));
    gmds::EdgeContainer edgeContainer(&mesh);

    gmds::Edge edge1 = edgeContainer.add(1, 2);
    gmds::Edge edge2 = edgeContainer.add(3, 4);

    ASSERT_TRUE(edgeContainer.has(edge1.id()));
    ASSERT_TRUE(edgeContainer.has(edge2.id()));

    edgeContainer.clear();

    ASSERT_FALSE(edgeContainer.has(edge1.id()));
    ASSERT_FALSE(edgeContainer.has(edge2.id()));
}
/*----------------------------------------------------------------------------*/
TEST(EdgeContainer, testResize)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E));
    gmds::EdgeContainer edgeContainer(&mesh);

    gmds::Edge edge1 = edgeContainer.add(1, 2);
    gmds::Edge edge2 = edgeContainer.add(3, 4);

    ASSERT_TRUE(edgeContainer.has(edge1.id()));
    ASSERT_TRUE(edgeContainer.has(edge2.id()));

    edgeContainer.resize(1);

    ASSERT_TRUE(edgeContainer.has(edge1.id()));
    ASSERT_FALSE(edgeContainer.has(edge2.id()));
}
