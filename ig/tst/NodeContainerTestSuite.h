/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/NodeContainer.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/Node.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
TEST(NodeContainer, testNodeContainerConstructor)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N));
    gmds::NodeContainer nodeContainer(&mesh);

    ASSERT_TRUE(nodeContainer.has(0));
    ASSERT_FALSE(nodeContainer.has(1));
}
/*----------------------------------------------------------------------------*/
TEST(NodeContainer, testAddAndGetNodesData)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E));
    gmds::NodeContainer nodeContainer(&mesh);
    gmds::Node n1 = nodeContainer.add(0, 0, 0);

    TInt nbNodes;
    nodeContainer.getNodesData(n1.id(), nbNodes);

    ASSERT_EQ(0, nbNodes);
}
/*----------------------------------------------------------------------------*/
TEST(NodeContainer, testAddAndGetEdgesData)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E));
    gmds::NodeContainer nodeContainer(&mesh);
    gmds::Node n1 = nodeContainer.add(0, 0, 0);
    gmds::Node n2 = nodeContainer.add(1, 0, 0);

    gmds::Edge e1 = mesh.newEdge(n1, n2);

    TInt nbEdges;
    nodeContainer.getEdgesData(n1.id(), nbEdges);

    ASSERT_EQ(1, nbEdges);
}
/*----------------------------------------------------------------------------*/
TEST(NodeContainer, testAddAndGetFacesData)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E | gmds::F));
    gmds::NodeContainer nodeContainer(&mesh);
    gmds::Node n1 = nodeContainer.add(0, 0, 0);
    gmds::Node n2 = nodeContainer.add(1, 0, 0);
    gmds::Node n3 = nodeContainer.add(1, 1, 0);

    gmds::Face f1 = mesh.newFace(n1, n2, n3);

    TInt nbFaces;
    nodeContainer.getFacesData(n1.id(), nbFaces);

    ASSERT_EQ(1, nbFaces);
}
/*----------------------------------------------------------------------------*/
TEST(NodeContainer, testAddAndGetRegionsData)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E | gmds::R));
    gmds::NodeContainer nodeContainer(&mesh);
    gmds::Node n1 = nodeContainer.add(0, 0, 0);
    gmds::Node n2 = nodeContainer.add(1, 0, 0);
    gmds::Node n3 = nodeContainer.add(1, 1, 0);

    gmds::Region r1 = mesh.newRegion(n1, n2, n3);

    TInt nbRegions;
    nodeContainer.getRegionsData(n1.id(), nbRegions);

    ASSERT_EQ(1, nbRegions);
}
/*----------------------------------------------------------------------------*/
TEST(NodeContainer, testAddAndBuildNode)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N));
    gmds::NodeContainer nodeContainer(&mesh);
    gmds::Node n1 = nodeContainer.add(0, 0, 0);

    gmds::Node builtNode = nodeContainer.buildNode(n1.id());

    ASSERT_EQ(n1.id(), builtNode.id());
}
/*----------------------------------------------------------------------------*/
TEST(NodeContainer, testHas)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N));
    gmds::NodeContainer nodeContainer(&mesh);

    ASSERT_FALSE(nodeContainer.has(0));

    gmds::Node n1 = nodeContainer.add(1, 2, 3);
    ASSERT_TRUE(nodeContainer.has(n1.id()));
}
/*----------------------------------------------------------------------------*/
TEST(NodeContainer, testClear)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N));
    gmds::NodeContainer nodeContainer(&mesh);

    gmds::Node n1 = nodeContainer.add(1, 2, 3);
    gmds::Node n2 = nodeContainer.add(4, 5, 6);

    ASSERT_TRUE(nodeContainer.has(n1.id()));
    ASSERT_TRUE(nodeContainer.has(n2.id()));

    nodeContainer.clear();

    ASSERT_FALSE(nodeContainer.has(n1.id()));
    ASSERT_FALSE(nodeContainer.has(n2.id()));
}
/*----------------------------------------------------------------------------*/
TEST(NodeContainer, testResize)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N));
    gmds::NodeContainer nodeContainer(&mesh);

    gmds::Node n1 = nodeContainer.add(1, 2, 3);
    gmds::Node n2 = nodeContainer.add(4, 5, 6);

    ASSERT_TRUE(nodeContainer.has(n1.id()));
    ASSERT_TRUE(nodeContainer.has(n2.id()));

    nodeContainer.resize(1);

    ASSERT_TRUE(nodeContainer.has(n1.id()));
    ASSERT_FALSE(nodeContainer.has(n2.id()));
}
/*----------------------------------------------------------------------------*/
TEST(NodeContainer, testUpdate)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N));
    gmds::NodeContainer nodeContainer(&mesh);

    gmds::Node n1 = nodeContainer.add(1, 2, 3);
    gmds::Node n2 = nodeContainer.add(4, 5, 6);

    ASSERT_TRUE(nodeContainer.has(n1.id()));
    ASSERT_TRUE(nodeContainer.has(n2.id()));

    nodeContainer.update();

    ASSERT_TRUE(nodeContainer.has(n1.id()));
    ASSERT_TRUE(nodeContainer.has(n2.id()));
}
/*----------------------------------------------------------------------------*/
TEST(NodeContainer, testSerializeAndUnserialize)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N));
    gmds::NodeContainer nodeContainer(&mesh);

    gmds::Node n1 = nodeContainer.add(1, 2, 3);
    gmds::Node n2 = nodeContainer.add(4, 5, 6);

    std::ostringstream oss;
    nodeContainer.serialize(oss);

    gmds::NodeContainer newContainer(&mesh);
    std::istringstream iss(oss.str());
    newContainer.unserialize(iss);

    ASSERT_TRUE(newContainer.has(n1.id()));
    ASSERT_TRUE(newContainer.has(n2.id()));
}