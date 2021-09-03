/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/AssignCells.h"
#include "ELG3D/ALGOCMPNT/ExtractGeomModel.h"
#include "ELG3D/ALGOCMPNT/InitData.h"
#include "ELG3D/ALGOCMPNT/Tools.h"
#include "ELG3D/DATACMPNT/FracPres.h"
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
class ExtractGeomTest : public ::testing::Test
{
 protected:
        ExtractGeomTest()
        {
                ;
        }
        virtual ~ExtractGeomTest()
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
//            int num_threads = 3;
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
TEST_F(ExtractGeomTest, dummy)
{
        EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(ExtractGeomTest, sample_3x3_2D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_3x3_2D(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getFaceCapacity());

    elg3d::assignCellsMajorityCriteria_2D(&mesh, &fp, &ma);

    elg3d::Tools_write_2D(&mesh, &fp, &ma, "ExtractGeomTest_sample_3x3_2D");


    kmds::Connectivity* c_E2F = mesh.createConnectivity(kmds::E2F);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildEandE2F_2D_variant_1();

    gmds::cad::FACManager AGeomModel;
    elg3d::extractGeomModel_extract_2D(&mesh, &ma, c_E2F, &AGeomModel);

//    EXPECT_EQ(4, AGeomModel.getNbPoints());
//    EXPECT_EQ(3, AGeomModel.getNbCurves());
}
/*----------------------------------------------------------------------------*/
TEST_F(ExtractGeomTest, sample_3x3x3_3D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_3x3x3_3D(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getFaceCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);

    elg3d::Tools_write_3D(&mesh, &fp, &ma, "ExtractGeomTest_sample_3x3x3_3D");


    kmds::Connectivity* c_F2R = mesh.createConnectivity(kmds::F2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildFandF2R();

    kmds::Variable <std::uintptr_t> *geomassoc = mesh.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr),
                                                                                                   kmds::KMDS_NODE, "geomassoc");
    gmds::cad::FACManager geomModel;

    elg3d::extractGeomModel_extract_3D(&mesh, &ma, c_F2R, &geomModel, geomassoc);

    EXPECT_EQ(3, geomModel.getNbSurfaces());
    EXPECT_EQ(4, geomModel.getNbCurves());
    // this is a multipoint
    EXPECT_EQ(1, geomModel.getNbPoints());
}
/*----------------------------------------------------------------------------*/