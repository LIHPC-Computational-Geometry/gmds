/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Google Test headers
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// GMDS headers
#include <KM/DS/FContainer.h>
/*----------------------------------------------------------------------------*/
// using namespace kmds;
/*----------------------------------------------------------------------------*/
class FContainersTest : public ::testing::Test
{
 protected:
        FContainersTest()
        {
                ;
        }
        virtual ~FContainersTest()
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
TEST_F(FContainersTest, init)
{
    kmds::FContainer fc(nullptr);

    EXPECT_EQ(0,fc.nbCells());
    EXPECT_EQ(0,fc.nbTriangles());
    EXPECT_EQ(0,fc.nbQuads());
    EXPECT_EQ(0,fc.nbPentagons());
}
/*----------------------------------------------------------------------------*/
TEST_F(FContainersTest, addCells)
{
    kmds::FContainer fc(nullptr);

    kmds::TCellID id0 = fc.addQuad();

    EXPECT_EQ(1,fc.nbCells());
    EXPECT_EQ(0,fc.nbTriangles());
    EXPECT_EQ(1,fc.nbQuads());
    EXPECT_EQ(0,fc.nbPentagons());

    kmds::TCellID id1 = fc.addQuad();

    EXPECT_EQ(2,fc.nbCells());
    EXPECT_EQ(2,fc.nbQuads());

    fc.remove(id0);

    EXPECT_EQ(1,fc.nbCells());
    EXPECT_EQ(0,fc.nbTriangles());
    EXPECT_EQ(1,fc.nbQuads());
    EXPECT_EQ(0,fc.nbPentagons());
}
/*----------------------------------------------------------------------------*/
TEST_F(FContainersTest, addCellsMultiple)
{
    kmds::FContainer fc(nullptr);

    kmds::TCellID id0 = fc.addQuads(10);

    EXPECT_EQ(10,fc.nbCells());
    EXPECT_EQ(0,fc.nbTriangles());
    EXPECT_EQ(10,fc.nbQuads());
    EXPECT_EQ(0,fc.nbPentagons());


}
/*----------------------------------------------------------------------------*/