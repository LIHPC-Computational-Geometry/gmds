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
    }

    static void
    TearDownTestCase()
    {
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