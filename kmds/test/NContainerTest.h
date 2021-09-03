/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Google Test headers
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// GMDS headers
#include <KM/DS/NContainer.h>
/*----------------------------------------------------------------------------*/
// using namespace kmds;
/*----------------------------------------------------------------------------*/
class NContainerTest : public ::testing::Test
{
protected:
    NContainerTest()
    {
        ;
    }
    virtual ~NContainerTest()
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
TEST_F(NContainerTest, init)
{
    kmds::NContainer nc(nullptr);

    EXPECT_EQ(0,nc.nbCells());
}
/*----------------------------------------------------------------------------*/
TEST_F(NContainerTest, addCells)
{
    kmds::NContainer nc(nullptr);

    kmds::TCellID id0 = nc.add();

    EXPECT_EQ(1,nc.nbCells());

    kmds::TCellID id1 = nc.add();

    EXPECT_EQ(2,nc.nbCells());

    nc.remove(id0);

    EXPECT_EQ(1,nc.nbCells());
}
/*----------------------------------------------------------------------------*/
TEST_F(NContainerTest, addCellsMultiple)
{
    kmds::NContainer nc(nullptr);

    kmds::TCellID id0 = nc.add(10);

    EXPECT_EQ(10,nc.nbCells());
}
/*----------------------------------------------------------------------------*/