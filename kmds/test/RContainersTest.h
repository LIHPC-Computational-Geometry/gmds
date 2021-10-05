/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Google Test headers
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// GMDS headers
#include <KM/DS/RContainer.h>
/*----------------------------------------------------------------------------*/
// using namespace kmds;
/*----------------------------------------------------------------------------*/
class RContainersTest : public ::testing::Test
{
 protected:
        RContainersTest()
        {
                ;
        }
        virtual ~RContainersTest()
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
//            int num_threads = 4;
//            int use_numa = 1;
//            int use_core = 1;
//            Kokkos::OpenMP::initialize(num_threads, use_numa, use_core);
            Kokkos::initialize(kargs);
    }

    static void
    TearDownTestCase()
    {
            // Kokkos::Serial::finalize();
            // Kokkos::Threads::finalize();
//            Kokkos::OpenMP::finalize();
            Kokkos::finalize();
    }
};
/*----------------------------------------------------------------------------*/
TEST_F(RContainersTest, init)
{
    kmds::RContainer rc(nullptr);

    EXPECT_EQ(0,rc.nbCells());
    EXPECT_EQ(0,rc.nbTetrahedra());
    EXPECT_EQ(0,rc.nbHexahedra());
    EXPECT_EQ(0,rc.nbPyramids());
}
/*----------------------------------------------------------------------------*/
TEST_F(RContainersTest, addCells)
{
    kmds::RContainer rc(nullptr);

    kmds::TCellID id0 = rc.addHexahedron();

    EXPECT_EQ(1,rc.nbCells());
    EXPECT_EQ(0,rc.nbTetrahedra());
    EXPECT_EQ(1,rc.nbHexahedra());
    EXPECT_EQ(0,rc.nbPyramids());

    kmds::TCellID id1 = rc.addHexahedron();

    EXPECT_EQ(2,rc.nbCells());
    EXPECT_EQ(2,rc.nbHexahedra());

    rc.remove(id0);

    EXPECT_EQ(1,rc.nbCells());
    EXPECT_EQ(0,rc.nbTetrahedra());
    EXPECT_EQ(1,rc.nbHexahedra());
    EXPECT_EQ(0,rc.nbPyramids());
}
/*----------------------------------------------------------------------------*/
TEST_F(RContainersTest, addCellsMultiple)
{
    kmds::RContainer rc(nullptr);

    kmds::TCellID id0 = rc.addHexahedra(10);

    EXPECT_EQ(10,rc.nbCells());
    EXPECT_EQ(0,rc.nbTetrahedra());
    EXPECT_EQ(10,rc.nbHexahedra());
    EXPECT_EQ(0,rc.nbPyramids());
}
/*----------------------------------------------------------------------------*/