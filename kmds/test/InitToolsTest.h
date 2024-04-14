/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Google Test headers
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Kokkos headers
//#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// GMDS headers
#include <KM/DS/Mesh.h>
#include <KM/Utils/InitTools.h>
/*----------------------------------------------------------------------------*/
// using namespace kmds;
/*----------------------------------------------------------------------------*/
class InitToolsTest : public ::testing::Test
{
 protected:
        InitToolsTest()
        {
                ;
        }
        virtual ~InitToolsTest()
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
TEST_F(InitToolsTest, init_2D)
{
        kmds::Mesh m;

        double xyz_min[3] = {0., 0., 0};
        double xyz_max[3] = {1., 1., 1};

        kmds::InitTools_createGrid_2D(&m, xyz_min, xyz_max, 3, 4);

        EXPECT_EQ(20, m.getNbNodes());
        EXPECT_EQ(12, m.getNbFaces());
        EXPECT_EQ(0, m.getNbRegions());
}
/*----------------------------------------------------------------------------*/
TEST_F(InitToolsTest, init_3D)
{
    kmds::Mesh m;

    double xyz_min[3] = {0., 0., 0};
    double xyz_max[3] = {1., 0.7, 2.};

    kmds::InitTools_createGrid_3D(&m, xyz_min, xyz_max, 3, 3, 2);

    EXPECT_EQ(48, m.getNbNodes());
    EXPECT_EQ(0, m.getNbFaces());
    EXPECT_EQ(18, m.getNbRegions());
}
/*----------------------------------------------------------------------------*/