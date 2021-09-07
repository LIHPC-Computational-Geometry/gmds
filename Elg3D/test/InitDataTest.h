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
                // Kokkos::Serial::initialize();
                // Kokkos::Threads::initialize();
                Kokkos::InitArguments kargs;
                kargs.num_threads = 3;
//                int num_threads = 4;
//                int use_numa = 1;
//                int use_core = 1;
//                Kokkos::OpenMP::initialize(num_threads, use_numa, use_core);
                Kokkos::initialize(kargs);
        }

        static void
        TearDownTestCase()
        {
                // Kokkos::Serial::finalize();
                // Kokkos::Threads::finalize();
//                Kokkos::OpenMP::finalize();
                Kokkos::finalize();
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