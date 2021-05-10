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
class FakeTypesTest : public ::testing::Test
{
protected:
    FakeTypesTest()
    {
        ;
    }
    virtual ~FakeTypesTest()
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
TEST_F(FakeTypesTest, hexpyrtetprism)
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

    kmds::Region r = m.getRegion(hexId);

    std::vector<kmds::FakeFace> fakeFaces = r.getFakeFaces();
    EXPECT_EQ(6, fakeFaces.size());

    for(auto ff: fakeFaces) {
        std::vector<kmds::FakeEdge> fakeEdges = ff.getFakeEdges();
        EXPECT_EQ(4, fakeEdges.size());
    }

    r = m.getRegion(pyrId);
    fakeFaces = r.getFakeFaces();
    EXPECT_EQ(5, fakeFaces.size());

    r = m.getRegion(tetId);
    fakeFaces = r.getFakeFaces();
    EXPECT_EQ(4, fakeFaces.size());

    r = m.getRegion(prismId);
    fakeFaces = r.getFakeFaces();
    EXPECT_EQ(5, fakeFaces.size());
}
/*----------------------------------------------------------------------------*/