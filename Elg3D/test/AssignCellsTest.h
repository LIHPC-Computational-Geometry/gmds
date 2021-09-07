/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/AssignCells.h"
#include "ELG3D/ALGOCMPNT/InitData.h"
#include "ELG3D/ALGOCMPNT/ManifoldDetection.h"
#include "ELG3D/ALGOCMPNT/MaterialInterfaces.h"
#include "ELG3D/ALGOCMPNT/Tools.h"
#include "ELG3D/DATACMPNT/FracPres.h"
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
class AssignCellsTest : public ::testing::Test
{
 protected:
        AssignCellsTest()
        {
                ;
        }
        virtual ~AssignCellsTest()
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
TEST_F(AssignCellsTest, dummy)
{
        EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(AssignCellsTest, majorityCriteria_3x3_2D)
{
        kmds::Mesh mesh;
        elg3d::FracPres fp;
        elg3d::initData_3x3_2D(&mesh, &fp);

        elg3d::MaterialAssignment ma(mesh.getFaceCapacity());

        elg3d::assignCellsMajorityCriteria_2D(&mesh, &fp, &ma);

        EXPECT_EQ(ma.getMaterial(0), 0);
        EXPECT_EQ(ma.getMaterial(1), 0);
        EXPECT_EQ(ma.getMaterial(3), 0);
        EXPECT_EQ(ma.getMaterial(4), 0);
        EXPECT_EQ(ma.getMaterial(2), 1);
        EXPECT_EQ(ma.getMaterial(5), 1);
        EXPECT_EQ(ma.getMaterial(6), 2);
        EXPECT_EQ(ma.getMaterial(7), 2);
        EXPECT_EQ(ma.getMaterial(8), 2);
}
/*----------------------------------------------------------------------------*/
TEST_F(AssignCellsTest, majorityCriteria_2x2_2D_nonmanifold)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_2x2_2D_nonmanifold(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getFaceCapacity());

    elg3d::assignCellsMajorityCriteria_2D(&mesh, &fp, &ma);

    EXPECT_EQ(ma.getMaterial(0), 0);
    EXPECT_EQ(ma.getMaterial(1), 1);
    EXPECT_EQ(ma.getMaterial(2), 1);
    EXPECT_EQ(ma.getMaterial(3), 0);
}
/*----------------------------------------------------------------------------*/
TEST_F(AssignCellsTest, majorityCriteria_3x3x3_3D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_3x3x3_3D(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);

    EXPECT_EQ(ma.getMaterial(0), 0);
    EXPECT_EQ(ma.getMaterial(3), 0);
    EXPECT_EQ(ma.getMaterial(6), 0);
    EXPECT_EQ(ma.getMaterial(1), 0);
    EXPECT_EQ(ma.getMaterial(4), 0);
    EXPECT_EQ(ma.getMaterial(7), 0);
    EXPECT_EQ(ma.getMaterial(9), 0);
    EXPECT_EQ(ma.getMaterial(12), 0);
    EXPECT_EQ(ma.getMaterial(15), 0);
    EXPECT_EQ(ma.getMaterial(10), 0);
    EXPECT_EQ(ma.getMaterial(13), 0);
    EXPECT_EQ(ma.getMaterial(16), 0);
}
/*----------------------------------------------------------------------------*/
TEST_F(AssignCellsTest, majorityCriteria_2x1x2_3D_nonmanifold)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_2x1x2_3D_nonmanifold(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);

    EXPECT_EQ(ma.getMaterial(0), 0);
    EXPECT_EQ(ma.getMaterial(1), 1);
    EXPECT_EQ(ma.getMaterial(2), 1);
    EXPECT_EQ(ma.getMaterial(3), 0);
}
/*----------------------------------------------------------------------------*/
TEST_F(AssignCellsTest, assignCells_2x1x2_3D_nonmanifold)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_2x1x2_3D_nonmanifold(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);

    EXPECT_EQ(ma.getMaterial(0), 0);
    EXPECT_EQ(ma.getMaterial(1), 1);
    EXPECT_EQ(ma.getMaterial(2), 1);
    EXPECT_EQ(ma.getMaterial(3), 0);

    // connectivities
    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2R();

    elg3d::assignCellsCorrection_3D(&mesh, c_N2R, &fp, &ma);

    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);

    kmds::GrowingView<kmds::TCellID> nodesManifold("NODES_MANIFOLD", nodesInterfaces.getNbElems());
    elg3d::ManifoldDetection_getNonManifoldNodes_3D(&nodesInterfaces, &mesh, c_N2R, &ma, &nodesManifold);

    EXPECT_EQ(0, nodesManifold.getNbElems());
}
/*----------------------------------------------------------------------------*/
TEST_F(AssignCellsTest, assignCells_I_2D_nonmanifold)
{
    const int nb_i = 30;

    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_I_2D_nonmanifold(&mesh, &fp, nb_i);

    elg3d::MaterialAssignment ma(mesh.getFaceCapacity());

    elg3d::assignCellsMajorityCriteria_2D(&mesh, &fp, &ma);

//    EXPECT_EQ(ma.getMaterial(0), 0);
//    EXPECT_EQ(ma.getMaterial(1), 1);
//    EXPECT_EQ(ma.getMaterial(2), 1);
//    EXPECT_EQ(ma.getMaterial(3), 0);


    elg3d::Tools_write_2D(&mesh, &fp, &ma, "AssignCellsTest_assignCells_I_2D_majority");


    // connectivities
    kmds::Connectivity* c_N2F = mesh.createConnectivity(kmds::N2F);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2F();

    elg3d::assignCellsCorrection_2D(&mesh, c_N2F, &fp, &ma);

    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2F, &ma, &nodesInterfaces);

    kmds::GrowingView<kmds::TCellID> nodesManifold("NODES_MANIFOLD", nodesInterfaces.getNbElems());
    elg3d::ManifoldDetection_getNonManifoldNodes_2D(&nodesInterfaces, &mesh, c_N2F, &ma, &nodesManifold);

    elg3d::Tools_write_2D(&mesh, &fp, &ma, "AssignCellsTest_assignCells_I_2D_correction");

    EXPECT_EQ(0, nodesManifold.getNbElems());
}
/*----------------------------------------------------------------------------*/
TEST_F(AssignCellsTest, assignCells_IxJxI_3D_nonmanifold)
{
    const int nb_i = 30;
    const int nb_j = 10;

    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_IxJxI_3D_nonmanifold(&mesh, &fp, nb_i, nb_j);

    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);

//    EXPECT_EQ(ma.getMaterial(0), 0);
//    EXPECT_EQ(ma.getMaterial(1), 1);
//    EXPECT_EQ(ma.getMaterial(2), 1);
//    EXPECT_EQ(ma.getMaterial(3), 0);


    elg3d::Tools_write_3D(&mesh, &fp, &ma, "AssignCellsTest_assignCells_IxJxI_3D_majority");


    // connectivities
    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2R();

    elg3d::assignCellsCorrection_3D(&mesh, c_N2R, &fp, &ma);

    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);

    kmds::GrowingView<kmds::TCellID> nodesManifold("NODES_MANIFOLD", nodesInterfaces.getNbElems());
    elg3d::ManifoldDetection_getNonManifoldNodes_3D(&nodesInterfaces, &mesh, c_N2R, &ma, &nodesManifold);

    elg3d::Tools_write_3D(&mesh, &fp, &ma, "AssignCellsTest_assignCells_IxJxI_3D_correction");

    EXPECT_EQ(0, nodesManifold.getNbElems());
}
/*----------------------------------------------------------------------------*/
TEST_F(AssignCellsTest, assignCells_rainbow_3D_nonmanifold)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_rainbow_3D_nonmanifold(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);

    // connectivities
    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2R();

    elg3d::assignCellsCorrection_3D(&mesh, c_N2R, &fp, &ma);

    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);

    kmds::GrowingView<kmds::TCellID> nodesManifold("NODES_MANIFOLD", nodesInterfaces.getNbElems());
    elg3d::ManifoldDetection_getNonManifoldNodes_3D(&nodesInterfaces, &mesh, c_N2R, &ma, &nodesManifold);

    EXPECT_EQ(0, nodesManifold.getNbElems());
}
/*----------------------------------------------------------------------------*/
TEST_F(AssignCellsTest, assignCells_2x2x2_3D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_2x2x2_3D(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);


    elg3d::Tools_write_3D(&mesh, &fp, &ma, "AssignCellsTest_assignCells_2x2x2_3D_majority");
}
/*----------------------------------------------------------------------------*/