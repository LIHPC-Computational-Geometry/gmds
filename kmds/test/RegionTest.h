/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Google Test headers
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// GMDS headers
#include <KM/DS/Mesh.h>
/*----------------------------------------------------------------------------*/
// using namespace kmds;
/*----------------------------------------------------------------------------*/
class RegionTest : public ::testing::Test
{
protected:
    RegionTest()
    {
        ;
    }
    virtual ~RegionTest()
    {
        ;
    }

    static void
    SetUpTestCase()
    {
        // Kokkos::Serial::initialize();
        // Kokkos::Threads::initialize();
        Kokkos::InitArguments kargs;
        kargs.num_threads = 3;
//        int num_threads = 4;
//        int use_numa = 1;
//        int use_core = 1;
//        Kokkos::OpenMP::initialize(num_threads, use_numa, use_core);
        Kokkos::initialize(kargs);
    }

    static void
    TearDownTestCase()
    {
        // Kokkos::Serial::finalize();
        // Kokkos::Threads::finalize();
//        Kokkos::OpenMP::finalize();
        Kokkos::finalize();
    }
};
/*----------------------------------------------------------------------------*/
TEST_F(RegionTest, RegionCreation_hexpyrtetprism)
{
    kmds::Mesh m;

    m.updateNodeCapacity(13);
    m.updateRegionCapacity(4);

    kmds::TCellID nodes[13];
    nodes[0] = m.newNode(0., 0., 0.);
    nodes[1] = m.newNode(1., 0., 0.);
    nodes[2] = m.newNode(1., 1., 0.);
    nodes[3] = m.newNode(0., 1., 0.);
    nodes[4] = m.newNode(0., 0., 1.);
    nodes[5] = m.newNode(1., 0., 1.);
    nodes[6] = m.newNode(1., 1., 1.);
    nodes[7] = m.newNode(0., 1., 1.);

    nodes[8] = m.newNode(0.5, 0.5, 2.);
    nodes[9] = m.newNode(1.5, 0.5, 1.5);

    nodes[10] = m.newNode(0., 1., 1.);
    nodes[11] = m.newNode(0., 1., 1.);
    nodes[12] = m.newNode(0., 1., 1.);


    kmds::TCellID hexId = m.newHexahedron(
            nodes[0],
            nodes[1],
            nodes[2],
            nodes[3],
            nodes[4],
            nodes[5],
            nodes[6],
            nodes[7]
    );

    kmds::TCellID pyrId = m.newPyramid(
            nodes[4],
            nodes[5],
            nodes[6],
            nodes[7],
            nodes[8]
    );

    kmds::TCellID tetId = m.newTetrahedron(
            nodes[5],
            nodes[6],
            nodes[8],
            nodes[9]
    );

    kmds::TCellID prismId = m.newPrism3(
            nodes[7],
            nodes[4],
            nodes[8],
            nodes[10],
            nodes[11],
            nodes[12]
    );

    EXPECT_EQ(4, m.getNbRegions());
    EXPECT_EQ(1, m.getNbHexahedra());
    EXPECT_EQ(1, m.getNbTetrahedra());
    EXPECT_EQ(1, m.getNbPyramids());
    EXPECT_EQ(1, m.getNbPrism3s());

    kmds::Region r = m.getRegion(hexId);
    Kokkos::View<kmds::TCellID *> ids;
    r.nodeIds(ids);

    EXPECT_EQ(kmds::KMDS_HEX, r.computeType());
    EXPECT_EQ(8, r.getNbNodes());
    EXPECT_EQ(8, ids.extent(0));
    EXPECT_EQ(8, ids.size());

    EXPECT_EQ(nodes[0], ids(0));
    EXPECT_EQ(nodes[1], ids(1));
    EXPECT_EQ(nodes[2], ids(2));
    EXPECT_EQ(nodes[3], ids(3));
    EXPECT_EQ(nodes[4], ids(4));
    EXPECT_EQ(nodes[5], ids(5));
    EXPECT_EQ(nodes[6], ids(6));
    EXPECT_EQ(nodes[7], ids(7));

    r = m.getRegion(pyrId);
    r.nodeIds(ids);

    EXPECT_EQ(kmds::KMDS_PYRAMID, r.computeType());
    EXPECT_EQ(5, r.getNbNodes());
    EXPECT_EQ(5, ids.extent(0));
    EXPECT_EQ(5, ids.size());

    EXPECT_EQ(nodes[4], ids(0));
    EXPECT_EQ(nodes[5], ids(1));
    EXPECT_EQ(nodes[6], ids(2));
    EXPECT_EQ(nodes[7], ids(3));
    EXPECT_EQ(nodes[8], ids(4));

    r = m.getRegion(tetId);
    r.nodeIds(ids);

    EXPECT_EQ(kmds::KMDS_TETRA, r.computeType());
    EXPECT_EQ(4, r.getNbNodes());
    EXPECT_EQ(4, ids.extent(0));
    EXPECT_EQ(4, ids.size());

    EXPECT_EQ(nodes[5], ids(0));
    EXPECT_EQ(nodes[6], ids(1));
    EXPECT_EQ(nodes[8], ids(2));
    EXPECT_EQ(nodes[9], ids(3));

    r = m.getRegion(prismId);
    r.nodeIds(ids);

    EXPECT_EQ(kmds::KMDS_PRISM3, r.computeType());
    EXPECT_EQ(6, r.getNbNodes());
    EXPECT_EQ(6, ids.extent(0));
    EXPECT_EQ(6, ids.size());

    EXPECT_EQ(nodes[7], ids(0));
    EXPECT_EQ(nodes[4], ids(1));
    EXPECT_EQ(nodes[8], ids(2));
    EXPECT_EQ(nodes[10], ids(3));
    EXPECT_EQ(nodes[11], ids(4));
    EXPECT_EQ(nodes[12], ids(5));

}
/*----------------------------------------------------------------------------*/