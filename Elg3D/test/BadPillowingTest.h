/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#include <set>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/AssignCells.h"
#include "ELG3D/ALGOCMPNT/BadPillowDetection.h"
#include "ELG3D/ALGOCMPNT/InitData.h"
#include "ELG3D/ALGOCMPNT/ManifoldDetection.h"
#include "ELG3D/ALGOCMPNT/MaterialInterfaces.h"
#include "ELG3D/DATACMPNT/FracPres.h"
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
class BadPillowingTest : public ::testing::Test
{
protected:
    BadPillowingTest()
    {
        ;
    }
    virtual ~BadPillowingTest()
    {
        ;
    }

    static void
    SetUpTestCase()
    {
//        // Kokkos::Serial::initialize();
//        // Kokkos::Threads::initialize();
//        Kokkos::InitArguments kargs;
//        kargs.num_threads = 1;
////        int num_threads = 3;
////        int use_numa = 1;
////        int use_core = 1;
////        Kokkos::OpenMP::initialize(num_threads, use_numa, use_core);
//        Kokkos::initialize(kargs);
    }

    static void
    TearDownTestCase()
    {
    }
};
/*----------------------------------------------------------------------------*/
TEST_F(BadPillowingTest, dummy)
{
    EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(BadPillowingTest, 2x2x2_active_3D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_2x2x2_badpillow_3D(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);

    EXPECT_EQ(ma.getMaterial(0), 0);
    EXPECT_EQ(ma.getMaterial(1), 1);
    EXPECT_EQ(ma.getMaterial(2), 0);
    EXPECT_EQ(ma.getMaterial(3), 0);
    EXPECT_EQ(ma.getMaterial(4), 1);
    EXPECT_EQ(ma.getMaterial(5), 1);
    EXPECT_EQ(ma.getMaterial(6), 1);
    EXPECT_EQ(ma.getMaterial(7), 0);

    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2R();

    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);

    kmds::GrowingView<kmds::TCellID> nodesManifold("NODES_MANIFOLD", nodesInterfaces.getNbElems());
    elg3d::ManifoldDetection_getNonManifoldNodes_3D(&nodesInterfaces, &mesh, c_N2R, &ma, &nodesManifold);

    elg3d::Tools_write_3D(&mesh, &fp, &ma, "poyop_before");

    EXPECT_EQ(1, nodesManifold.getNbElems());
    EXPECT_EQ(13, nodesManifold.get(0));

    elg3d::assignCellsCorrection_3D(&mesh, c_N2R, &fp, &ma);

    elg3d::Tools_write_3D(&mesh, &fp, &ma, "poyop_after");

    kmds::GrowingView<kmds::TCellID> nodesInterfaces_after("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces_after);

    kmds::GrowingView<kmds::TCellID> nodesManifold_after("NODES_MANIFOLD", nodesInterfaces_after.getNbElems());
    elg3d::ManifoldDetection_getNonManifoldNodes_3D(&nodesInterfaces_after, &mesh, c_N2R, &ma, &nodesManifold_after);

    EXPECT_EQ(0, nodesManifold_after.getNbElems());
    EXPECT_EQ(ma.getMaterial(7), 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(BadPillowingTest, 2x2x2_ortho_3D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_2x2x2_badpillow_3D(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);

    EXPECT_EQ(ma.getMaterial(0), 0);
    EXPECT_EQ(ma.getMaterial(1), 1);
    EXPECT_EQ(ma.getMaterial(2), 0);
    EXPECT_EQ(ma.getMaterial(3), 0);
    EXPECT_EQ(ma.getMaterial(4), 1);
    EXPECT_EQ(ma.getMaterial(5), 1);
    EXPECT_EQ(ma.getMaterial(6), 1);
    EXPECT_EQ(ma.getMaterial(7), 0);

    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2R();

    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);

    kmds::GrowingView<kmds::TCellID> nodesManifold("NODES_MANIFOLD", nodesInterfaces.getNbElems());
    elg3d::ManifoldDetection_getNonManifoldNodes_3D(&nodesInterfaces, &mesh, c_N2R, &ma, &nodesManifold);

    elg3d::Tools_write_3D(&mesh, &fp, &ma, "poyop_before");

    EXPECT_EQ(1, nodesManifold.getNbElems());
    EXPECT_EQ(13, nodesManifold.get(0));

//    elg3d::assignCellsCorrection_3D(&mesh, c_N2R, &fp, &ma);
//
//    elg3d::Tools_write_3D(&mesh, &fp, &ma, "poyop_after");
//
//    kmds::GrowingView<kmds::TCellID> nodesInterfaces_after("NODES_ON_INTERFACES", mesh.getNbNodes());
//    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces_after);
//
//    kmds::GrowingView<kmds::TCellID> nodesManifold_after("NODES_MANIFOLD", nodesInterfaces_after.getNbElems());
//    elg3d::ManifoldDetection_getNonManifoldNodes_3D(&nodesInterfaces_after, &mesh, c_N2R, &ma, &nodesManifold_after);
//
//    EXPECT_EQ(0, nodesManifold_after.getNbElems());
//    EXPECT_EQ(ma.getMaterial(7), 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(BadPillowingTest, 2x2x2_good_3D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_2x2x2_badpillow_3D(&mesh, &fp);

    mesh.setNodeLocation(22, 2, 1, 1.5);
    mesh.setNodeLocation(25, 2, 2, 1.5);

    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);

    EXPECT_EQ(ma.getMaterial(0), 0);
    EXPECT_EQ(ma.getMaterial(1), 1);
    EXPECT_EQ(ma.getMaterial(2), 0);
    EXPECT_EQ(ma.getMaterial(3), 0);
    EXPECT_EQ(ma.getMaterial(4), 1);
    EXPECT_EQ(ma.getMaterial(5), 1);
    EXPECT_EQ(ma.getMaterial(6), 1);
    EXPECT_EQ(ma.getMaterial(7), 0);

    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2R();

    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);

    kmds::GrowingView<kmds::TCellID> nodesManifold("NODES_MANIFOLD", nodesInterfaces.getNbElems());
    elg3d::ManifoldDetection_getNonManifoldNodes_3D(&nodesInterfaces, &mesh, c_N2R, &ma, &nodesManifold);

    elg3d::Tools_write_3D(&mesh, &fp, &ma, "poyop_before");

    EXPECT_EQ(0, nodesManifold.getNbElems());
}
/*----------------------------------------------------------------------------*/
TEST_F(BadPillowingTest, 2x2x2_bad_3D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_2x2x2_badpillow_3D(&mesh, &fp);

    mesh.setNodeLocation(22, 2, 1, 0.5);
    mesh.setNodeLocation(25, 2, 2, 0.5);

    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);

    EXPECT_EQ(ma.getMaterial(0), 0);
    EXPECT_EQ(ma.getMaterial(1), 1);
    EXPECT_EQ(ma.getMaterial(2), 0);
    EXPECT_EQ(ma.getMaterial(3), 0);
    EXPECT_EQ(ma.getMaterial(4), 1);
    EXPECT_EQ(ma.getMaterial(5), 1);
    EXPECT_EQ(ma.getMaterial(6), 1);
    EXPECT_EQ(ma.getMaterial(7), 0);

    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2R();

    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);

    kmds::GrowingView<kmds::TCellID> nodesManifold("NODES_MANIFOLD", nodesInterfaces.getNbElems());
    elg3d::ManifoldDetection_getNonManifoldNodes_3D(&nodesInterfaces, &mesh, c_N2R, &ma, &nodesManifold);

    elg3d::Tools_write_3D(&mesh, &fp, &ma, "poyop_before");

    EXPECT_EQ(1, nodesManifold.getNbElems());
    EXPECT_EQ(13, nodesManifold.get(0));

//    elg3d::assignCellsCorrection_3D(&mesh, c_N2R, &fp, &ma);
//
//    elg3d::Tools_write_3D(&mesh, &fp, &ma, "poyop_after");
//
//    kmds::GrowingView<kmds::TCellID> nodesInterfaces_after("NODES_ON_INTERFACES", mesh.getNbNodes());
//    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces_after);
//
//    kmds::GrowingView<kmds::TCellID> nodesManifold_after("NODES_MANIFOLD", nodesInterfaces_after.getNbElems());
//    elg3d::ManifoldDetection_getNonManifoldNodes_3D(&nodesInterfaces_after, &mesh, c_N2R, &ma, &nodesManifold_after);
//
//    EXPECT_EQ(0, nodesManifold_after.getNbElems());
}
/*----------------------------------------------------------------------------*/