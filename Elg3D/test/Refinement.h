/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#include <set>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/Refinement.h"
/*----------------------------------------------------------------------------*/
#include <KM/IO/VTKWriter.h>
#include <KM/Utils/InitTools.h>
/*----------------------------------------------------------------------------*/
class RefinementTest : public ::testing::Test
{
protected:
    RefinementTest()
    {
        ;
    }
    virtual ~RefinementTest()
    {
        ;
    }

    static void
    SetUpTestCase()
    {
        // Kokkos::Serial::initialize();
        // Kokkos::Threads::initialize();
        Kokkos::InitArguments kargs;
        kargs.num_threads = 1;
//        int num_threads = 3;
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
TEST_F(RefinementTest, dummy)
{
    EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(RefinementTest, 1x1_corner_2D)
{
    kmds::Mesh mesh;
    double xyz_min[3] = {0., 0., 0};
    double xyz_max[3] = {1., 1., 0};

    kmds::InitTools_createGrid_2D(&mesh, xyz_min, xyz_max, 1, 1);

    EXPECT_TRUE(elg3d::Refinement_validityCheck_2D(&mesh));

    kmds::GrowingView <kmds::TCellID> nodeIDs2Refine("NODES", mesh.getNbNodes());
    nodeIDs2Refine.push_back(0);

    kmds::ConnectivityHelper ch(&mesh);
    kmds::Connectivity *c_N2F = mesh.createConnectivity(kmds::N2F);
    ch.buildN2F();
    kmds::Connectivity *c_F2F_byN = mesh.createConnectivity(kmds::F2F_byN);
    ch.buildF2F_byN(c_N2F);

    kmds::Connectivity* c_E2F = mesh.createConnectivity(kmds::E2F);
    ch.buildEandE2F_2D_variant_1();

    elg3d::Refinement_refine_2D(&nodeIDs2Refine, false, &mesh, c_F2F_byN, c_E2F);

    EXPECT_EQ(3, mesh.getNbFaces());

    kmds::VTKWriter<kmds::Mesh> w(mesh);
    w.write("RefinementTest_1x1_corner_2D", kmds::F);
}
/*----------------------------------------------------------------------------*/
TEST_F(RefinementTest, 1x1_edge_2D)
{
    kmds::Mesh mesh;
    double xyz_min[3] = {0., 0., 0};
    double xyz_max[3] = {1., 1., 0};

    kmds::InitTools_createGrid_2D(&mesh, xyz_min, xyz_max, 1, 1);

    EXPECT_TRUE(elg3d::Refinement_validityCheck_2D(&mesh));

    kmds::GrowingView <kmds::TCellID> nodeIDs2Refine("NODES", mesh.getNbNodes());
    nodeIDs2Refine.push_back(0);
    nodeIDs2Refine.push_back(1);

    kmds::ConnectivityHelper ch(&mesh);
    kmds::Connectivity *c_N2F = mesh.createConnectivity(kmds::N2F);
    ch.buildN2F();
    kmds::Connectivity *c_F2F_byN = mesh.createConnectivity(kmds::F2F_byN);
    ch.buildF2F_byN(c_N2F);

    kmds::Connectivity* c_E2F = mesh.createConnectivity(kmds::E2F);
    ch.buildEandE2F_2D_variant_1();

    elg3d::Refinement_refine_2D(&nodeIDs2Refine, false, &mesh, c_F2F_byN, c_E2F);

    EXPECT_EQ(7, mesh.getNbFaces());

    kmds::VTKWriter<kmds::Mesh> w(mesh);
    w.write("RefinementTest_1x1_edge_2D", kmds::F);
}
/*----------------------------------------------------------------------------*/
TEST_F(RefinementTest, 1x1_edge_right_2D)
{
    kmds::Mesh mesh;
    double xyz_min[3] = {0., 0., 0};
    double xyz_max[3] = {1., 1., 0};

    kmds::InitTools_createGrid_2D(&mesh, xyz_min, xyz_max, 1, 1);

    EXPECT_TRUE(elg3d::Refinement_validityCheck_2D(&mesh));

    kmds::GrowingView <kmds::TCellID> nodeIDs2Refine("NODES", mesh.getNbNodes());
    nodeIDs2Refine.push_back(2);
    nodeIDs2Refine.push_back(3);

    kmds::ConnectivityHelper ch(&mesh);
    kmds::Connectivity *c_N2F = mesh.createConnectivity(kmds::N2F);
    ch.buildN2F();
    kmds::Connectivity *c_F2F_byN = mesh.createConnectivity(kmds::F2F_byN);
    ch.buildF2F_byN(c_N2F);

    kmds::Connectivity* c_E2F = mesh.createConnectivity(kmds::E2F);
    ch.buildEandE2F_2D_variant_1();

    elg3d::Refinement_refine_2D(&nodeIDs2Refine, false, &mesh, c_F2F_byN, c_E2F);

    EXPECT_EQ(7, mesh.getNbFaces());

    kmds::VTKWriter<kmds::Mesh> w(mesh);
    w.write("RefinementTest_1x1_edge_right_2D", kmds::F);
}
/*----------------------------------------------------------------------------*/
TEST_F(RefinementTest, 1x1_diag_2D)
{
    kmds::Mesh mesh;
    double xyz_min[3] = {0., 0., 0};
    double xyz_max[3] = {1., 1., 0};

    kmds::InitTools_createGrid_2D(&mesh, xyz_min, xyz_max, 1, 1);

    EXPECT_TRUE(elg3d::Refinement_validityCheck_2D(&mesh));

    kmds::GrowingView <kmds::TCellID> nodeIDs2Refine("NODES", mesh.getNbNodes());
    nodeIDs2Refine.push_back(0);
    nodeIDs2Refine.push_back(3);

    kmds::ConnectivityHelper ch(&mesh);
    kmds::Connectivity *c_N2F = mesh.createConnectivity(kmds::N2F);
    ch.buildN2F();
    kmds::Connectivity *c_F2F_byN = mesh.createConnectivity(kmds::F2F_byN);
    ch.buildF2F_byN(c_N2F);

    kmds::Connectivity* c_E2F = mesh.createConnectivity(kmds::E2F);
    ch.buildEandE2F_2D_variant_1();

    elg3d::Refinement_refine_2D(&nodeIDs2Refine, false, &mesh, c_F2F_byN, c_E2F);

    EXPECT_EQ(7, mesh.getNbFaces());

    kmds::VTKWriter<kmds::Mesh> w(mesh);
    w.write("RefinementTest_1x1_diag_2D", kmds::F);
}
/*----------------------------------------------------------------------------*/
TEST_F(RefinementTest, 1x1_three_2D)
{
    kmds::Mesh mesh;
    double xyz_min[3] = {0., 0., 0};
    double xyz_max[3] = {1., 1., 0};

    kmds::InitTools_createGrid_2D(&mesh, xyz_min, xyz_max, 1, 1);

    EXPECT_TRUE(elg3d::Refinement_validityCheck_2D(&mesh));

    kmds::GrowingView <kmds::TCellID> nodeIDs2Refine("NODES", mesh.getNbNodes());
    nodeIDs2Refine.push_back(1);
    nodeIDs2Refine.push_back(2);
    nodeIDs2Refine.push_back(3);

    kmds::ConnectivityHelper ch(&mesh);
    kmds::Connectivity *c_N2F = mesh.createConnectivity(kmds::N2F);
    ch.buildN2F();
    kmds::Connectivity *c_F2F_byN = mesh.createConnectivity(kmds::F2F_byN);
    ch.buildF2F_byN(c_N2F);

    kmds::Connectivity* c_E2F = mesh.createConnectivity(kmds::E2F);
    ch.buildEandE2F_2D_variant_1();

    elg3d::Refinement_refine_2D(&nodeIDs2Refine, false, &mesh, c_F2F_byN, c_E2F);

    EXPECT_EQ(8, mesh.getNbFaces());

    kmds::VTKWriter<kmds::Mesh> w(mesh);
    w.write("RefinementTest_1x1_three_2D", kmds::F);
}
/*----------------------------------------------------------------------------*/
TEST_F(RefinementTest, 1x1_all_2D)
{
    kmds::Mesh mesh;
    double xyz_min[3] = {0., 0., 0};
    double xyz_max[3] = {1., 1., 0};

    kmds::InitTools_createGrid_2D(&mesh, xyz_min, xyz_max, 1, 1);

    EXPECT_TRUE(elg3d::Refinement_validityCheck_2D(&mesh));

    kmds::GrowingView <kmds::TCellID> nodeIDs2Refine("NODES", mesh.getNbNodes());
    nodeIDs2Refine.push_back(0);

    kmds::ConnectivityHelper ch(&mesh);
    kmds::Connectivity *c_N2F = mesh.createConnectivity(kmds::N2F);
    ch.buildN2F();
    kmds::Connectivity *c_F2F_byN = mesh.createConnectivity(kmds::F2F_byN);
    ch.buildF2F_byN(c_N2F);

    kmds::Connectivity* c_E2F = mesh.createConnectivity(kmds::E2F);
    ch.buildEandE2F_2D_variant_1();

    elg3d::Refinement_refine_2D(&nodeIDs2Refine, true, &mesh, c_F2F_byN, c_E2F);

    EXPECT_EQ(9, mesh.getNbFaces());

    kmds::VTKWriter<kmds::Mesh> w(mesh);
    w.write("RefinementTest_1x1_all_2D", kmds::F);
}
/*----------------------------------------------------------------------------*/
TEST_F(RefinementTest, 3x3_2D)
{
    kmds::Mesh mesh;
    double xyz_min[3] = {0., 0., 0};
    double xyz_max[3] = {1., 1., 0};

    kmds::InitTools_createGrid_2D(&mesh, xyz_min, xyz_max, 3, 3);

    EXPECT_TRUE(elg3d::Refinement_validityCheck_2D(&mesh));

    kmds::GrowingView <kmds::TCellID> cellIDs2Refine("CELLS", mesh.getNbFaces());
    cellIDs2Refine.push_back(0);

    kmds::ConnectivityHelper ch(&mesh);
    kmds::Connectivity *c_N2F = mesh.createConnectivity(kmds::N2F);
    ch.buildN2F();
    kmds::Connectivity *c_F2F_byN = mesh.createConnectivity(kmds::F2F_byN);
    ch.buildF2F_byN(c_N2F);

    kmds::Connectivity* c_E2F = mesh.createConnectivity(kmds::E2F);
    ch.buildEandE2F_2D_variant_1();

    elg3d::Refinement_refine_2D(&cellIDs2Refine, true, &mesh, c_F2F_byN, c_E2F);

    EXPECT_EQ(5+9+7+7+3, mesh.getNbFaces());

    kmds::VTKWriter<kmds::Mesh> w(mesh);
    w.write("RefinementTest_3x3_2D", kmds::F);
}
/*----------------------------------------------------------------------------*/
TEST_F(RefinementTest, 1x1x1_corner_3D)
{
    kmds::Mesh mesh;
    double xyz_min[3] = {0., 0., 0};
    double xyz_max[3] = {3., 3., 3};

    kmds::InitTools_createGrid_3D(&mesh, xyz_min, xyz_max, 1, 1, 1);

    EXPECT_TRUE(elg3d::Refinement_validityCheck_3D(&mesh));

    kmds::GrowingView <kmds::TCellID> cellIDs2Refine("CELLS", mesh.getNbRegions());
    cellIDs2Refine.push_back(0);

    kmds::ConnectivityHelper ch(&mesh);
    kmds::Connectivity *c_N2R = mesh.createConnectivity(kmds::N2R);
    ch.buildN2R();
    kmds::Connectivity *c_R2R_byN = mesh.createConnectivity(kmds::R2R_byN);
    ch.buildF2F_byN(c_N2R);

    kmds::Connectivity* c_E2R = mesh.createConnectivity(kmds::E2R);
    ch.buildEandE2R();

    kmds::Connectivity* c_F2R = mesh.createConnectivity(kmds::F2R);
    ch.buildFandF2R();

    elg3d::Refinement_refine_3D(&cellIDs2Refine, false, &mesh, c_R2R_byN, c_E2R, c_F2R);

    EXPECT_EQ(4, mesh.getNbRegions());

    kmds::VTKWriter<kmds::Mesh> w(mesh);
    w.write("RefinementTest_1x1x1_corner_3D", kmds::R);
}
/*----------------------------------------------------------------------------*/
TEST_F(RefinementTest, 1x1x1_edge_3D)
{
    kmds::Mesh mesh;
    double xyz_min[3] = {0., 0., 0};
    double xyz_max[3] = {3., 3., 3};

    kmds::InitTools_createGrid_3D(&mesh, xyz_min, xyz_max, 1, 1, 1);

    EXPECT_TRUE(elg3d::Refinement_validityCheck_3D(&mesh));

    kmds::GrowingView <kmds::TCellID> cellIDs2Refine("CELLS", mesh.getNbRegions());
    cellIDs2Refine.push_back(0);
    cellIDs2Refine.push_back(1);

    kmds::ConnectivityHelper ch(&mesh);
    kmds::Connectivity *c_N2R = mesh.createConnectivity(kmds::N2R);
    ch.buildN2R();
    kmds::Connectivity *c_R2R_byN = mesh.createConnectivity(kmds::R2R_byN);
    ch.buildF2F_byN(c_N2R);

    kmds::Connectivity* c_E2R = mesh.createConnectivity(kmds::E2R);
    ch.buildEandE2R();

    kmds::Connectivity* c_F2R = mesh.createConnectivity(kmds::F2R);
    ch.buildFandF2R();

    elg3d::Refinement_refine_3D(&cellIDs2Refine, false, &mesh, c_R2R_byN, c_E2R, c_F2R);

    EXPECT_EQ(11, mesh.getNbRegions());

    kmds::VTKWriter<kmds::Mesh> w(mesh);
    w.write("RefinementTest_1x1x1_edge_3D", kmds::R);
}
/*----------------------------------------------------------------------------*/
TEST_F(RefinementTest, 1x1x1_face_3D)
{
    kmds::Mesh mesh;
    double xyz_min[3] = {0., 0., 0};
    double xyz_max[3] = {3., 3., 3};

    kmds::InitTools_createGrid_3D(&mesh, xyz_min, xyz_max, 1, 1, 1);

    EXPECT_TRUE(elg3d::Refinement_validityCheck_3D(&mesh));

    kmds::GrowingView <kmds::TCellID> cellIDs2Refine("CELLS", mesh.getNbRegions());
    cellIDs2Refine.push_back(0);
    cellIDs2Refine.push_back(1);
    cellIDs2Refine.push_back(2);
    cellIDs2Refine.push_back(3);


    kmds::ConnectivityHelper ch(&mesh);
    kmds::Connectivity *c_N2R = mesh.createConnectivity(kmds::N2R);
    ch.buildN2R();
    kmds::Connectivity *c_R2R_byN = mesh.createConnectivity(kmds::R2R_byN);
    ch.buildF2F_byN(c_N2R);

    kmds::Connectivity* c_E2R = mesh.createConnectivity(kmds::E2R);
    ch.buildEandE2R();

    kmds::Connectivity* c_F2R = mesh.createConnectivity(kmds::F2R);
    ch.buildFandF2R();

    elg3d::Refinement_refine_3D(&cellIDs2Refine, false, &mesh, c_R2R_byN, c_E2R, c_F2R);

    EXPECT_EQ(22, mesh.getNbRegions());

    kmds::VTKWriter<kmds::Mesh> w(mesh);
    w.write("RefinementTest_1x1x1_face_3D", kmds::R);
}
/*----------------------------------------------------------------------------*/
TEST_F(RefinementTest, 1x1x1_all_3D)
{
    kmds::Mesh mesh;
    double xyz_min[3] = {0., 0., 0};
    double xyz_max[3] = {3., 3., 3};

    kmds::InitTools_createGrid_3D(&mesh, xyz_min, xyz_max, 1, 1, 1);

    EXPECT_TRUE(elg3d::Refinement_validityCheck_3D(&mesh));

    kmds::GrowingView <kmds::TCellID> cellIDs2Refine("CELLS", mesh.getNbRegions());
    cellIDs2Refine.push_back(0);

    kmds::ConnectivityHelper ch(&mesh);
    kmds::Connectivity *c_N2R = mesh.createConnectivity(kmds::N2R);
    ch.buildN2R();
    kmds::Connectivity *c_R2R_byN = mesh.createConnectivity(kmds::R2R_byN);
    ch.buildF2F_byN(c_N2R);

    kmds::Connectivity* c_E2R = mesh.createConnectivity(kmds::E2R);
    ch.buildEandE2R();

    kmds::Connectivity* c_F2R = mesh.createConnectivity(kmds::F2R);
    ch.buildFandF2R();

    elg3d::Refinement_refine_3D(&cellIDs2Refine, true, &mesh, c_R2R_byN, c_E2R, c_F2R);

    EXPECT_EQ(27, mesh.getNbRegions());

    kmds::VTKWriter<kmds::Mesh> w(mesh);
    w.write("RefinementTest_1x1x1_all_3D", kmds::R);
}
/*----------------------------------------------------------------------------*/
TEST_F(RefinementTest, 3x3x3_3D)
{
    kmds::Mesh mesh;
    double xyz_min[3] = {0., 0., 0};
    double xyz_max[3] = {3., 3., 3};

    kmds::InitTools_createGrid_3D(&mesh, xyz_min, xyz_max, 3, 3, 3);

    EXPECT_TRUE(elg3d::Refinement_validityCheck_3D(&mesh));

    kmds::GrowingView <kmds::TCellID> cellIDs2Refine("CELLS", mesh.getNbRegions());
    cellIDs2Refine.push_back(0);
    cellIDs2Refine.push_back(4);


    kmds::ConnectivityHelper ch(&mesh);
    kmds::Connectivity *c_N2R = mesh.createConnectivity(kmds::N2R);
    ch.buildN2R();
    kmds::Connectivity *c_R2R_byN = mesh.createConnectivity(kmds::R2R_byN);
    ch.buildF2F_byN(c_N2R);

    kmds::Connectivity* c_E2R = mesh.createConnectivity(kmds::E2R);
    ch.buildEandE2R();

    kmds::Connectivity* c_F2R = mesh.createConnectivity(kmds::F2R);
    ch.buildFandF2R();

    elg3d::Refinement_refine_3D(&cellIDs2Refine, true, &mesh, c_R2R_byN, c_E2R, c_F2R);

    kmds::VTKWriter<kmds::Mesh> w(mesh);
    w.write("RefinementTest_3x3x3_3D", kmds::R);
}

/*----------------------------------------------------------------------------*/