/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#include <set>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/AssignCells.h"
#include "ELG3D/ALGOCMPNT/InitData.h"
#include "ELG3D/ALGOCMPNT/ManifoldDetection.h"
#include "ELG3D/ALGOCMPNT/MaterialInterfaces.h"
#include "ELG3D/DATACMPNT/FracPres.h"
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
class ManifoldDetectionTest : public ::testing::Test
{
protected:
    ManifoldDetectionTest()
    {
        ;
    }
    virtual ~ManifoldDetectionTest()
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
TEST_F(ManifoldDetectionTest, dummy)
{
    EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(ManifoldDetectionTest, 3x3_2D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_3x3_2D(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getFaceCapacity());

    elg3d::assignCellsMajorityCriteria_2D(&mesh, &fp, &ma);

    kmds::Connectivity* c_N2F = mesh.createConnectivity(kmds::N2F);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2F();

    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2F, &ma, &nodesInterfaces);


}
/*----------------------------------------------------------------------------*/
TEST_F(ManifoldDetectionTest, 2x2_2D_nonmanifold)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_2x2_2D_nonmanifold(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getFaceCapacity());

    elg3d::assignCellsMajorityCriteria_2D(&mesh, &fp, &ma);

    kmds::Connectivity* c_N2F = mesh.createConnectivity(kmds::N2F);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2F();

    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2F, &ma, &nodesInterfaces);

    EXPECT_EQ(5, nodesInterfaces.getNbElems());

}
/*----------------------------------------------------------------------------*/
TEST_F(ManifoldDetectionTest, 3x3x3_3D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_3x3x3_3D(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);

    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2R();

    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);

    kmds::GrowingView<kmds::TCellID> nodesManifold("NODES_MANIFOLD", nodesInterfaces.getNbElems());
    elg3d::ManifoldDetection_getNonManifoldNodes_3D(&nodesInterfaces, &mesh, c_N2R, &ma, &nodesManifold);

    EXPECT_EQ(0, nodesManifold.getNbElems());

    kmds::GrowingView<kmds::TCellID> nodesInterfaces_matA("NODES_ON_INTERFACES_MATA", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(ma.getMaterialID("matA"), &mesh, c_N2R, &ma, &nodesInterfaces_matA);

    kmds::GrowingView<kmds::TCellID> nodesManifold_matA("NODES_MANIFOLD_MATA", nodesInterfaces.getNbElems());
    elg3d::ManifoldDetection_getNonManifoldNodes_3D(&nodesInterfaces_matA, ma.getMaterialID("matA"), &mesh, c_N2R, &ma, &nodesManifold_matA);

    EXPECT_EQ(0, nodesManifold_matA.getNbElems());

    kmds::GrowingView<kmds::TCellID> nodesInterfaces_matB("NODES_ON_INTERFACES_MATB", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(ma.getMaterialID("matB"), &mesh, c_N2R, &ma, &nodesInterfaces_matB);

    kmds::GrowingView<kmds::TCellID> nodesManifold_matB("NODES_MANIFOLD_MATB", nodesInterfaces.getNbElems());
    elg3d::ManifoldDetection_getNonManifoldNodes_3D(&nodesInterfaces_matB, ma.getMaterialID("matB"), &mesh, c_N2R, &ma, &nodesManifold_matB);

    EXPECT_EQ(0, nodesManifold_matB.getNbElems());

    kmds::GrowingView<kmds::TCellID> nodesInterfaces_matC("NODES_ON_INTERFACES_MATC", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(ma.getMaterialID("matC"), &mesh, c_N2R, &ma, &nodesInterfaces_matC);

    kmds::GrowingView<kmds::TCellID> nodesManifold_matC("NODES_MANIFOLD_MATC", nodesInterfaces.getNbElems());
    elg3d::ManifoldDetection_getNonManifoldNodes_3D(&nodesInterfaces_matC, ma.getMaterialID("matC"), &mesh, c_N2R, &ma, &nodesManifold_matC);

    EXPECT_EQ(0, nodesManifold_matC.getNbElems());
}
/*----------------------------------------------------------------------------*/
TEST_F(ManifoldDetectionTest, 2x1x2_3D_nonmanifold)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_2x1x2_3D_nonmanifold(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);

    kmds::Connectivity *c_N2R = mesh.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2R();

    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);

    kmds::GrowingView<kmds::TCellID> nodesManifold("NODES_MANIFOLD", nodesInterfaces.getNbElems());
    elg3d::ManifoldDetection_getNonManifoldNodes_3D(&nodesInterfaces, &mesh, c_N2R, &ma, &nodesManifold);

    EXPECT_EQ(2, nodesManifold.getNbElems());

    std::set<kmds::TCellID> nodesManifoldSet;
    for(int i=0; i<nodesManifold.getNbElems(); i++) {
        nodesManifoldSet.insert(nodesManifold.get(i));
    }

    EXPECT_TRUE(nodesManifoldSet.find(7) != nodesManifoldSet.end());
    EXPECT_TRUE(nodesManifoldSet.find(10) != nodesManifoldSet.end());

    kmds::GrowingView<kmds::TCellID> nodesInterfaces_matA("NODES_ON_INTERFACES_MATA", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(ma.getMaterialID("matA"), &mesh, c_N2R, &ma, &nodesInterfaces_matA);

    kmds::GrowingView<kmds::TCellID> nodesManifold_matA("NODES_MANIFOLD_MATA", nodesInterfaces.getNbElems());
    elg3d::ManifoldDetection_getNonManifoldNodes_3D(&nodesInterfaces_matA, ma.getMaterialID("matA"),&mesh, c_N2R, &ma, &nodesManifold_matA);

    EXPECT_EQ(2, nodesManifold_matA.getNbElems());

    kmds::GrowingView<kmds::TCellID> nodesInterfaces_matB("NODES_ON_INTERFACES_MATB", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(ma.getMaterialID("matB"), &mesh, c_N2R, &ma, &nodesInterfaces_matB);

    kmds::GrowingView<kmds::TCellID> nodesManifold_matB("NODES_MANIFOLD_MATB", nodesInterfaces.getNbElems());
    elg3d::ManifoldDetection_getNonManifoldNodes_3D(&nodesInterfaces_matB, ma.getMaterialID("matB"),&mesh, c_N2R, &ma, &nodesManifold_matB);

    EXPECT_EQ(2, nodesManifold_matB.getNbElems());
}
/*----------------------------------------------------------------------------*/
TEST_F(ManifoldDetectionTest, 2x1x2_3D_computeCost)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_2x1x2_3D_nonmanifold(&mesh, &fp);

// compute and store cells volume
    kmds::GrowingView<kmds::TCellID> cellIDAccessor("CELLS", mesh.getNbRegions());
    mesh.getRegionIDs_dummy(&cellIDAccessor);
    kmds::Variable<kmds::TCoord>* v = mesh.createVariable<kmds::TCoord>(-HUGE_VALF, kmds::KMDS_REGION, "cellVolume");
    for(int icell=0; icell<cellIDAccessor.getNbElems(); icell++) {
        kmds::TCellID id = cellIDAccessor.get(icell);
        (*v)[id] = mesh.getRegion(id).surfvol();
    }

    std::vector<kmds::TCellID> cellIDs = {0,1,2,3};

    std::vector<int> cells2mat_0110 = {0,1,1,0};
    double cost_0110 = elg3d::ManifoldDetection_ComputeCost(cellIDs, cells2mat_0110, &fp, v);
    EXPECT_DOUBLE_EQ(0.6, cost_0110);

    std::vector<int> cells2mat_1010 = {1,0,1,0};
    double cost_1010 = elg3d::ManifoldDetection_ComputeCost(cellIDs, cells2mat_1010, &fp, v);
    EXPECT_DOUBLE_EQ(2.2, cost_1010);

    std::vector<int> cells2mat_0101 = {0,1,0,1};
    double cost_0101= elg3d::ManifoldDetection_ComputeCost(cellIDs, cells2mat_0101, &fp, v);
    EXPECT_DOUBLE_EQ(1.8, cost_0101);

    std::vector<int> cells2mat_0000 = {0,0,0,0};
    double cost_0000 = elg3d::ManifoldDetection_ComputeCost(cellIDs, cells2mat_0000, &fp, v);
    EXPECT_DOUBLE_EQ(1.6, cost_0000);
}
/*----------------------------------------------------------------------------*/
TEST_F(ManifoldDetectionTest, 2x1x2_3D_exploreConfigs)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_2x1x2_3D_nonmanifold(&mesh, &fp);

// compute and store cells volume
    kmds::GrowingView<kmds::TCellID> cellIDAccessor("CELLS", mesh.getNbRegions());
    mesh.getRegionIDs_dummy(&cellIDAccessor);
    kmds::Variable<kmds::TCoord>* v = mesh.createVariable<kmds::TCoord>(-HUGE_VALF, kmds::KMDS_REGION, "cellVolume");
    for(int icell=0; icell<cellIDAccessor.getNbElems(); icell++) {
        kmds::TCellID id = cellIDAccessor.get(icell);
        (*v)[id] = mesh.getRegion(id).surfvol();
    }

    std::vector<kmds::TCellID> cellIDs = {0,1,2,3};
    std::vector<int> cells2mat_0110;

    Kokkos::UnorderedMap<kmds::TCellID, elg3d::ManifoldDetection_configStorage> storageMap(mesh.getNbNodes());

    elg3d::ManifoldDetection_SolveNonManifoldAtOneNode_3D(&mesh, 7, cellIDs, &fp, v, cells2mat_0110, &storageMap);


}
/*----------------------------------------------------------------------------*/