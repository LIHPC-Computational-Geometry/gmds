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
    }

    static void
    TearDownTestCase()
    {
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