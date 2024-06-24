/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Cell.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/Node.h>
#include <gmds/ig/Edge.h>
#include <gmds/ig/Face.h>
#include <gmds/ig/Region.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
TEST(Cell, testCellConstructorAndGetters)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::R));
    gmds::Cell cell(&mesh, gmds::GMDS_NODE, 1);

    ASSERT_EQ(1, cell.id());
    ASSERT_EQ(gmds::GMDS_NODE, cell.type());
}
/*----------------------------------------------------------------------------*/
TEST(Cell, testComputeBoundingBoxForNodes)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N));
    gmds::Cell node(&mesh, gmds::GMDS_NODE, 1);

    EXPECT_THROW(node.computeBoundingBox(nullptr, nullptr), gmds::GMDSException);
}
/*----------------------------------------------------------------------------*/
TEST(Cell, testComputeBoundingBoxForElements)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::R));
    gmds::Node n1 = mesh.newNode(0, 0, 0);
    gmds::Node n2 = mesh.newNode(1, 1, 1);
    gmds::Node n3 = mesh.newNode(0, 1, 1);
    gmds::Node n4 = mesh.newNode(1, 0, 1);

    gmds::Cell element(&mesh, gmds::GMDS_QUAD, 1);
    element.set<std::vector<gmds::Node>>({n1, n2, n3, n4});

    TCoord minXYZ[3], maxXYZ[3];
    element.computeBoundingBox(minXYZ, maxXYZ);

    ASSERT_EQ(0, minXYZ[0]);
    ASSERT_EQ(0, minXYZ[1]);
    ASSERT_EQ(0, minXYZ[2]);
    ASSERT_EQ(1, maxXYZ[0]);
    ASSERT_EQ(1, maxXYZ[1]);
    ASSERT_EQ(1, maxXYZ[2]);
}
/*----------------------------------------------------------------------------*/
TEST(Cell, testSerializeAndUnserializeCellData)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::R));
    gmds::Cell cell(&mesh, gmds::GMDS_QUAD, 1);

    std::stringstream ss;
    cell.serializeCellData(ss);

    gmds::Cell newCell(&mesh, gmds::GMDS_NODE, 2);
    newCell.unserializeCellData(ss);

    ASSERT_EQ(gmds::GMDS_QUAD, newCell.type());
    ASSERT_EQ(1, newCell.id());
}
/*----------------------------------------------------------------------------*/
TEST(Cell, testGetIDsAndSetForNodes)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N));
    gmds::Node n1 = mesh.newNode(0, 0, 0);
    gmds::Node n2 = mesh.newNode(1, 1, 1);
    gmds::Node n3 = mesh.newNode(0, 1, 1);

    gmds::Cell node(&mesh, gmds::GMDS_NODE, 1);
    node.set<std::vector<gmds::Node>>({n1, n2, n3});

    std::vector<gmds::TCellID> ids;
    node.getIDs<gmds::Node>(ids);

    ASSERT_EQ(3, ids.size());
    ASSERT_EQ(0, ids[0]);
    ASSERT_EQ(1, ids[1]);
    ASSERT_EQ(2, ids[2]);
}
/*----------------------------------------------------------------------------*/
TEST(Cell, testAddRemoveAndReplaceForEdges)
{
    gmds::Mesh mesh(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E));
    gmds::Node n1 = mesh.newNode(0, 0, 0);
    gmds::Node n2 = mesh.newNode(1, 0, 0);
    gmds::Node n3 = mesh.newNode(1, 1, 0);
    gmds::Node n4 = mesh.newNode(0, 1, 0);

    gmds::Edge e1 = mesh.newEdge(n1, n2);
    gmds::Edge e2 = mesh.newEdge(n2, n3);

    gmds::Cell element(&mesh, gmds::GMDS_QUAD, 1);
    element.set<std::vector<gmds::Edge>>({e1, e2});

    std::vector<gmds::TCellID> edgeIDs;
    element.getIDs<gmds::Edge>(edgeIDs);

    ASSERT_EQ(2, edgeIDs.size());
    ASSERT_EQ(0, edgeIDs[0]);
    ASSERT_EQ(1, edgeIDs[1]);

    element.remove<gmds::Edge>(1);
    edgeIDs.clear();
    element.getIDs<gmds::Edge>(edgeIDs);

    ASSERT_EQ(1, edgeIDs.size());
    ASSERT_EQ(0, edgeIDs[0]);

    element.add<gmds::Edge>(2);
    edgeIDs.clear();
    element.getIDs<gmds::Edge>(edgeIDs);

    ASSERT_EQ(2, edgeIDs.size());
    ASSERT_EQ(0, edgeIDs[0]);
    ASSERT_EQ(2, edgeIDs[1]);

    element.replace<gmds::Edge>(0, 3);
    edgeIDs.clear();
    element.getIDs<gmds::Edge>(edgeIDs);

    ASSERT_EQ(2, edgeIDs.size());
    ASSERT_EQ(3, edgeIDs[0]);
    ASSERT_EQ(2, edgeIDs[1]);
}