/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Node.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/EdgeContainer.h>
#include <gmds/ig/Edge.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
TEST(Node, testNodeConstructor)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N));
    gmds::Node node(&mesh, 1, 0.0, 0.0, 0.0);

    ASSERT_EQ(1, node.id());
    ASSERT_EQ(0.0, node.X());
    ASSERT_EQ(0.0, node.Y());
    ASSERT_EQ(0.0, node.Z());
}
/*----------------------------------------------------------------------------*/
TEST(Node, testNodeConstructorWithPoint)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N));
    gmds::math::Point pt(1.0, 2.0, 3.0);
    gmds::Node node(&mesh, 1, pt);

    ASSERT_EQ(1, node.id());
    ASSERT_EQ(1.0, node.X());
    ASSERT_EQ(2.0, node.Y());
    ASSERT_EQ(3.0, node.Z());
}
/*----------------------------------------------------------------------------*/
TEST(Node, testNodeCopyConstructor)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N));
    gmds::Node node1(&mesh, 1, 0.0, 0.0, 0.0);
    gmds::Node node2 = node1;

    ASSERT_EQ(node1.id(), node2.id());
    ASSERT_EQ(node1.X(), node2.X());
    ASSERT_EQ(node1.Y(), node2.Y());
    ASSERT_EQ(node1.Z(), node2.Z());
}
/*----------------------------------------------------------------------------*/
TEST(Node, testNodeEqualityOperator)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N));
    gmds::Node node1(&mesh, 1, 0.0, 0.0, 0.0);
    gmds::Node node2(&mesh, 2, 0.0, 0.0, 0.0);
    gmds::Node node3(&mesh, 1, 1.0, 0.0, 0.0);

    ASSERT_TRUE(node1 == node1);
    ASSERT_FALSE(node1 == node2);
    ASSERT_FALSE(node1 == node3);
}
/*----------------------------------------------------------------------------*/
TEST(Node, testNodeInequalityOperator)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N));
    gmds::Node node1(&mesh, 1, 0.0, 0.0, 0.0);
    gmds::Node node2(&mesh, 2, 0.0, 0.0, 0.0);
    gmds::Node node3(&mesh, 1, 1.0, 0.0, 0.0);

    ASSERT_FALSE(node1 != node1);
    ASSERT_TRUE(node1 != node2);
    ASSERT_TRUE(node1 != node3);
}
/*----------------------------------------------------------------------------*/
TEST(Node, testNodeCenter)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N));
    gmds::Node node(&mesh, 1, 1.0, 2.0, 3.0);

    gmds::math::Point center = node.center();

    ASSERT_EQ(1.0, center.X());
    ASSERT_EQ(2.0, center.Y());
    ASSERT_EQ(3.0, center.Z());
}
/*----------------------------------------------------------------------------*/
TEST(Node, testNodeNbNodes)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E));
    gmds::Node node(&mesh, 1, 0.0, 0.0, 0.0);
    mesh.newNode(2, 0.0, 1.0, 0.0);
    mesh.newNode(3, 1.0, 0.0, 0.0);
    mesh.newNode(4, 1.0, 1.0, 0.0);

    TInt nbNodes = node.nbNodes();

    ASSERT_EQ(3, nbNodes);
}
/*----------------------------------------------------------------------------*/
TEST(Node, testNodeNbEdges)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E));
    gmds::Node node(&mesh, 1, 0.0, 0.0, 0.0);
    gmds::Edge e1 = mesh.newEdge(2, 3);
    gmds::Edge e2 = mesh.newEdge(3, 4);
    gmds::Edge e3 = mesh.newEdge(2, 4);

    TInt nbEdges = node.nbEdges();

    ASSERT_EQ(3, nbEdges);
}
/*----------------------------------------------------------------------------*/
TEST(Node, testNodeNbFaces)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E | gmds::F));
    gmds::Node node(&mesh, 1, 0.0, 0.0, 0.0);
    mesh.newNode(2, 0.0, 1.0, 0.0);
    mesh.newNode(3, 1.0, 0.0, 0.0);
    gmds::Face f = mesh.newFace(1, 2, 3);

    TInt nbFaces = node.nbFaces();

    ASSERT_EQ(1, nbFaces);
}
/*----------------------------------------------------------------------------*/
TEST(Node, testNodeNbRegions)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E | gmds::R));
    gmds::Node node(&mesh, 1, 0.0, 0.0, 0.0);
    mesh.newNode(2, 0.0, 1.0, 0.0);
    mesh.newNode(3, 1.0, 0.0, 0.0);
    gmds::Region r = mesh.newRegion(1, 2, 3);

    TInt nbRegions = node.nbRegions();

    ASSERT_EQ(1, nbRegions);
}
/*----------------------------------------------------------------------------*/
TEST(Node, testNodeDelegateGetNodes)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E));
    gmds::Node node(&mesh, 1, 0.0, 0.0, 0.0);
    mesh.newNode(2, 0.0, 1.0, 0.0);
    mesh.newNode(3, 1.0, 0.0, 0.0);
    mesh.newNode(4, 1.0, 1.0, 0.0);

    std::vector<gmds::Node> nodes;
    node.delegateGet(nodes);

    ASSERT_EQ(3, nodes.size());
    ASSERT_EQ(2, nodes[0].id());
    ASSERT_EQ(3, nodes[1].id());
    ASSERT_EQ(4, nodes[2].id());
}
/*----------------------------------------------------------------------------*/
TEST(Node, testNodeDelegateGetEdges)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E));
    gmds::Node node(&mesh, 1, 0.0, 0.0, 0.0);
    gmds::Edge e1 = mesh.newEdge(2, 3);
    gmds::Edge e2 = mesh.newEdge(3, 4);
    gmds::Edge e3 = mesh.newEdge(2, 4);

    std::vector<gmds::Edge> edges;
    node.delegateGet(edges);

    ASSERT_EQ(3, edges.size());
    ASSERT_EQ(e1.id(), edges[0].id());
    ASSERT_EQ(e2.id(), edges[1].id());
    ASSERT_EQ(e3.id(), edges[2].id());
}
/*----------------------------------------------------------------------------*/
TEST(Node, testNodeDelegateGetFaces)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E | gmds::F));
    gmds::Node node(&mesh, 1, 0.0, 0.0, 0.0);
    mesh.newNode(2, 0.0, 1.0, 0.0);
    mesh.newNode(3, 1.0, 0.0, 0.0);
    gmds::Face f = mesh.newFace(1, 2, 3);

    std::vector<gmds::Face> faces;
    node.delegateGet(faces);

    ASSERT_EQ(1, faces.size());
    ASSERT_EQ(f.id(), faces[0].id());
}
/*----------------------------------------------------------------------------*/
TEST(Node, testNodeDelegateGetRegions)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E | gmds::R));
    gmds::Node node(&mesh, 1, 0.0, 0.0, 0.0);
    mesh.newNode(2, 0.0, 1.0, 0.0);
    mesh.newNode(3, 1.0, 0.0, 0.0);
    gmds::Region r = mesh.newRegion(1, 2, 3);

    std::vector<gmds::Region> regions;
    node.delegateGet(regions);

    ASSERT_EQ(1, regions.size());
    ASSERT_EQ(r.id(), regions[0].id());
}
/*----------------------------------------------------------------------------*/
TEST(Node, testNodeDelegateGetNodeIDs)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E));
    gmds::Node node(&mesh, 1, 0.0, 0.0, 0.0);
    mesh.newNode(2, 0.0, 1.0, 0.0);
    mesh.newNode(3, 1.0, 0.0, 0.0);
    mesh.newNode(4, 1.0, 1.0, 0.0);

    std::vector<TCellID> nodeIDs;
    node.delegateGetNodeIDs(nodeIDs);

    ASSERT_EQ(3, nodeIDs.size());
    ASSERT_EQ(2, nodeIDs[0]);
    ASSERT_EQ(3, nodeIDs[1]);
    ASSERT_EQ(4, nodeIDs[2]);
}
/*----------------------------------------------------------------------------*/
TEST(Node, testNodeDelegateGetEdgeIDs)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E));
    gmds::Node node(&mesh, 1, 0.0, 0.0, 0.0);
    gmds::Edge e1 = mesh.newEdge(2, 3);
    gmds::Edge e2 = mesh.newEdge(3, 4);
    gmds::Edge e3 = mesh.newEdge(2, 4);

    std::vector<TCellID> edgeIDs;
    node.delegateGetEdgeIDs(edgeIDs);

    ASSERT_EQ(3, edgeIDs.size());
    ASSERT_EQ(e1.id(), edgeIDs[0]);
    ASSERT_EQ(e2.id(), edgeIDs[1]);
    ASSERT_EQ(e3.id(), edgeIDs[2]);
}