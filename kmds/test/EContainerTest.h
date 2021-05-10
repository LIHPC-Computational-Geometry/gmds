/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Google Test headers
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// GMDS headers
#include <KM/DS/EContainer.h>
/*----------------------------------------------------------------------------*/
// using namespace kmds;
/*----------------------------------------------------------------------------*/
class EContainerTest : public ::testing::Test
{
protected:
    EContainerTest()
    {
        ;
    }
    virtual ~EContainerTest()
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
TEST_F(EContainerTest, init)
{
    kmds::EContainer ec(nullptr);

    EXPECT_EQ(0,ec.nbCells());
}
/*----------------------------------------------------------------------------*/
TEST_F(EContainerTest, addCells)
{
    kmds::EContainer ec(nullptr);

    kmds::TCellID id0 = ec.add();

    EXPECT_EQ(1,ec.nbCells());

    kmds::TCellID id1 = ec.add();

    EXPECT_EQ(2,ec.nbCells());

    ec.remove(id0);

    EXPECT_EQ(1,ec.nbCells());
}
/*----------------------------------------------------------------------------*/
TEST_F(EContainerTest, addCellsMultiple)
{
    kmds::EContainer ec(nullptr);

    kmds::TCellID id0 = ec.add(10);

    EXPECT_EQ(10,ec.nbCells());
}
/*----------------------------------------------------------------------------*/