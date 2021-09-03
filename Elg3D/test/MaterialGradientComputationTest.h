/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/AssignCells.h"
#include "ELG3D/ALGOCMPNT/InitData.h"
#include "ELG3D/ALGOCMPNT/InterfaceNodesPos.h"
#include "ELG3D/ALGOCMPNT/ManifoldDetection.h"
#include "ELG3D/ALGOCMPNT/MaterialGradientComputation.h"
#include "ELG3D/ALGOCMPNT/MaterialInterfaces.h"
#include "ELG3D/ALGOCMPNT/Tools.h"
#include "ELG3D/DATACMPNT/FracPres.h"
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
class MaterialGradientComputationTest : public ::testing::Test
{
protected:
    MaterialGradientComputationTest()
    {
        ;
    }
    virtual ~MaterialGradientComputationTest()
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
TEST_F(MaterialGradientComputationTest, dummy)
{
    EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(MaterialGradientComputationTest, majorityCriteria_3x3_2D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_3x3_2D(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_2D(&mesh, &fp, &ma);


    // store cells midpoints
    kmds::TSize nbCells = mesh.getNbFaces();
    kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", nbCells);
    mesh.getFaceIDs_dummy(&cellIDs);
    kmds::Variable<gmds::math::Point>* varMidpoints = mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_FACE, "midpoint");
    elg3d::Tools_computeMidPointVar_2D(&cellIDs, &mesh, varMidpoints);


    // gradient computation
    kmds::Connectivity* c_N2F = mesh.createConnectivity(kmds::N2F);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2F();
    kmds::Connectivity* c_F2F_byN = mesh.createConnectivity(kmds::F2F_byN);
    ch.buildF2F_byN(c_N2F);

    kmds::Variable<elg3d::MaterialGradientComputation_midcellGradients>* varMidcellgrads =
            mesh.createVariable<elg3d::MaterialGradientComputation_midcellGradients>(elg3d::MaterialGradientComputation_MIDCELLGRADIENTS_NULL, kmds::KMDS_FACE, "midcellgrads");
    elg3d::MaterialGradientComputation_leastsquare_2D(&cellIDs, c_F2F_byN, &fp, &ma, varMidpoints, varMidcellgrads);

    // write data
    elg3d::Tools_write_2D(&mesh, &fp, &ma, varMidcellgrads, "MaterialGradientComputationTest_majorityCriteria_3x3_2D");
}
/*----------------------------------------------------------------------------*/
TEST_F(MaterialGradientComputationTest, majorityCriteria_3x3x3_3D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_3x3x3_3D(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);


    // store cells midpoints
    kmds::TSize nbCells = mesh.getNbRegions();
    kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", nbCells);
    mesh.getRegionIDs_dummy(&cellIDs);
    kmds::Variable<gmds::math::Point>* varMidpoints = mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_REGION, "midpoint");
    elg3d::Tools_computeMidPointVar_3D(&cellIDs, &mesh, varMidpoints);


    // gradient computation
    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2R();
    kmds::Connectivity* c_R2R_byN = mesh.createConnectivity(kmds::R2R_byN);
    ch.buildR2R_byN(c_N2R);

    kmds::Variable<elg3d::MaterialGradientComputation_midcellGradients>* varMidcellgrads =
            mesh.createVariable<elg3d::MaterialGradientComputation_midcellGradients>(elg3d::MaterialGradientComputation_MIDCELLGRADIENTS_NULL, kmds::KMDS_REGION, "midcellgrads");
    elg3d::MaterialGradientComputation_leastsquare_3D(&cellIDs, c_R2R_byN, &fp, &ma, varMidpoints, varMidcellgrads);

    // write data
    elg3d::Tools_write_3D(&mesh, &fp, &ma, varMidcellgrads, "MaterialGradientComputationTest_majorityCriteria_3x3x3_3D");
}
/*----------------------------------------------------------------------------*/
TEST_F(MaterialGradientComputationTest, assignCells_IxJxI_3D_nonmanifold)
{
//    const int nb_i = 30;
//    const int nb_j = 10;
//
//    kmds::Mesh mesh;
//    elg3d::FracPres fp;
//    elg3d::initData_IxJxI_3D_nonmanifold(&mesh, &fp, nb_i, nb_j);
//
//    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());
//
//    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);
//
//    // connectivities
//    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
//    kmds::ConnectivityHelper ch(&mesh);
//    ch.buildN2R();
//
//    assignCellsCorrection_3D(&mesh, c_N2R, &fp, &ma);
//
//    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
//    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);
//
//    kmds::GrowingView<kmds::TCellID> nodesManifold("NODES_MANIFOLD", nodesInterfaces.getNbElems());
//    elg3d::ManifoldDetection_getNonManifoldNodes_3D(&nodesInterfaces, &mesh, c_N2R, &ma, &nodesManifold);

}
/*----------------------------------------------------------------------------*/