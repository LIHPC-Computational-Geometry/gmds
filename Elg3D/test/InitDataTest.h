/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <KM/DS/Mesh.h>
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/InitData.h"
#include "ELG3D/DATACMPNT/FracPres.h"
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
class InitDataTest : public ::testing::Test
{
 protected:
        InitDataTest()
        {
                ;
        }
        virtual ~InitDataTest()
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
TEST_F(InitDataTest, dummy)
{
        EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(InitDataTest, createData_3x3_2D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;

    elg3d::initData_3x3_2D(&mesh, &fp);

    EXPECT_EQ(9, mesh.getNbFaces());
    EXPECT_EQ(9, mesh.getNbQuads());
    EXPECT_EQ(16, mesh.getNbNodes());

    EXPECT_EQ(3, fp.getNbMaterials());
    EXPECT_EQ(0., fp.getFracPres(1,0));
    EXPECT_EQ(0.2, fp.getFracPres(1,1));
}
/*----------------------------------------------------------------------------*/
TEST_F(InitDataTest, createData_3x3x3_3D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_3x3x3_3D(&mesh, &fp);

    EXPECT_EQ(27, mesh.getNbRegions());
    EXPECT_EQ(27, mesh.getNbHexahedra());
    EXPECT_EQ(64, mesh.getNbNodes());

    EXPECT_EQ(3, fp.getNbMaterials());
    EXPECT_EQ(0., fp.getFracPres(1,0));
    EXPECT_EQ(0.2, fp.getFracPres(1,1));
}
/*----------------------------------------------------------------------------*/