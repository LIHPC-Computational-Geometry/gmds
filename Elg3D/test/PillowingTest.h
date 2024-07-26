/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <KM/IO/VTKWriter.h>
#include <KM/Utils/InitTools.h>
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/AssignCells.h"
#include "ELG3D/ALGOCMPNT/InitData.h"
#include "ELG3D/ALGOCMPNT/InterfaceNodesPos.h"
#include "ELG3D/ALGOCMPNT/ManifoldDetection.h"
#include "ELG3D/ALGOCMPNT/MaterialGradientComputation.h"
#include "ELG3D/ALGOCMPNT/MaterialInterfaces.h"
#include "ELG3D/ALGOCMPNT/MoveToNewPos.h"
#include "ELG3D/ALGOCMPNT/Pillow.h"
#include "ELG3D/ALGOCMPNT/SmartLaplacian.h"
#include "ELG3D/ALGOCMPNT/Tools.h"
#include "ELG3D/DATACMPNT/FracPres.h"
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
#include "ELG3D/DATACMPNT/Parameters.h"
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
class PillowingTest : public ::testing::Test
{
protected:
    PillowingTest()
    {
        ;
    }
    virtual ~PillowingTest()
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
////            int num_threads = 3;
////            int use_numa = 1;
////            int use_core = 1;
////            Kokkos::OpenMP::initialize(num_threads, use_numa, use_core);
//        Kokkos::initialize(kargs);
    }

    static void
    TearDownTestCase()
    {
    }
};
/*----------------------------------------------------------------------------*/
TEST_F(PillowingTest, dummy)
{
    EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(PillowingTest, majorityCriteria_3x3_2D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_3x3_2D(&mesh, &fp);


    kmds::Variable<std::uintptr_t>* v = mesh.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr), kmds::KMDS_NODE, "geomEntity");
    gmds::cad::FACManager* facGeomManager = new gmds::cad::FACManager;

    elg3d::BoundingBoxGeomAssociation_init_2D(&mesh, facGeomManager, v);


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


    // compute interface nodes new position
    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2F, &ma, &nodesInterfaces);

    ch.buildEandE2F();
    kmds::Connectivity* c_E2F = mesh.getConnectivity(kmds::E2F);
    kmds::Connectivity* c_N2E = mesh.createConnectivity(kmds::N2E);
    ch.buildN2E();

    kmds::Variable<gmds::math::Point>* varNewPos =
            mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_NODE, "newPos");

    elg3d::InterfaceNodesPos_computeNodesNewPos_2D(&nodesInterfaces,
                                                   &mesh,
                                                   c_N2E,
                                                   c_E2F,
                                                   &ma,
                                                   &fp,
                                                   varMidpoints,
                                                   varMidcellgrads,
                                                   varNewPos);


//    std::cout<<0<<" "<<(*varNewPos)[0]<<std::endl;
    std::cout<<2<<" "<<(*varNewPos)[2]<<std::endl;
    std::cout<<8<<" "<<(*varNewPos)[8]<<std::endl;
//    std::cout<<18<<" "<<(*varNewPos)[18]<<std::endl;


    elg3d::moveToNewPos_basicMove(&nodesInterfaces,
                                  &mesh,
                                  varNewPos,
                                  v);

    kmds::GrowingView<kmds::TCellID> facesInterfaces("EDGES_ON_INTERFACES", mesh.getNbEdges());
    elg3d::MaterialInterfaces_getEdgeOnInterfaces(&mesh, c_E2F, &ma, &facesInterfaces);


    kmds::VTKWriter<kmds::Mesh> w(mesh);
    w.write("PillowingTest_before_3x3_2D", kmds::F);

    elg3d::pillow_2D(&facesInterfaces,
                     &nodesInterfaces,
                     &mesh,
                     c_E2F,
                     c_N2F,
                     &ma,
                     v);

    w.write("PillowingTest_after_3x3_2D", kmds::F);


    mesh.deleteConnectivity(kmds::N2F);
    c_N2F = mesh.createConnectivity(kmds::N2F);
    ch.buildN2F();

    ch.buildN2N(kmds::F);
    kmds::Connectivity* c_N2N = mesh.getConnectivity(kmds::N2N);

    elg3d::smartLaplacian_2D(10,
                             &mesh,
                             c_N2F,
                             c_N2N,
                             v,
                             &nodesInterfaces);

    w.write("PillowingTest_smoothed_3x3x2D", kmds::F);

}
/*----------------------------------------------------------------------------*/
TEST_F(PillowingTest, majorityCriteria_3x3x3_3D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_3x3x3_3D(&mesh, &fp);


    kmds::Variable<std::uintptr_t>* v = mesh.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr), kmds::KMDS_NODE, "geomEntity");
    gmds::cad::FACManager* facGeomManager = new gmds::cad::FACManager;

    elg3d::BoundingBoxGeomAssociation_init_3D(&mesh, facGeomManager, v);


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


    // compute interface nodes new position
    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);

    ch.buildFandF2R();
    kmds::Connectivity* c_F2R = mesh.getConnectivity(kmds::F2R);
    kmds::Connectivity* c_N2F = mesh.createConnectivity(kmds::N2F);
    ch.buildN2F();

    kmds::Variable<gmds::math::Point>* varNewPos =
            mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_NODE, "newPos");

    elg3d::InterfaceNodesPos_computeNodesNewPos_3D(&nodesInterfaces,
                                                   &mesh,
                                                   c_N2F,
                                                   c_F2R,
                                                   &ma,
                                                   &fp,
                                                   varMidpoints,
                                                   varMidcellgrads,
                                                   varNewPos);


    std::cout<<0<<" "<<(*varNewPos)[0]<<std::endl;
    std::cout<<2<<" "<<(*varNewPos)[2]<<std::endl;
    std::cout<<17<<" "<<(*varNewPos)[17]<<std::endl;
    std::cout<<18<<" "<<(*varNewPos)[18]<<std::endl;


    elg3d::moveToNewPos_basicMove(&nodesInterfaces,
                                  &mesh,
                                  varNewPos,
                                  v);


    kmds::GrowingView<kmds::TCellID> facesInterfaces("FACES_ON_INTERFACES", mesh.getNbFaces());
    elg3d::MaterialInterfaces_getFaceOnInterfaces(&mesh, c_F2R, &ma, &facesInterfaces);


    kmds::VTKWriter<kmds::Mesh> w(mesh);
    w.write("PillowingTest_before_3x3x3_3D", kmds::F|kmds::R);

    elg3d::pillow_3D(&facesInterfaces,
                     &nodesInterfaces,
                     &mesh,
                     c_F2R,
                     c_N2R,
                     &ma,
                     v);

    w.write("PillowingTest_after_3x3x3_3D", kmds::F|kmds::R);


//    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
    mesh.deleteConnectivity(kmds::N2R);
    c_N2R = mesh.createConnectivity(kmds::N2R);
    ch.buildN2R();

    ch.buildN2N(kmds::R);
    kmds::Connectivity* c_N2N = mesh.getConnectivity(kmds::N2N);

    elg3d::smartLaplacian_3D(10,
                             &mesh,
                             c_N2R,
                             c_N2N,
                             v,
                             &nodesInterfaces);

    w.write("PillowingTest_smoothed_3x3x3_3D", kmds::R);



}
/*----------------------------------------------------------------------------*/
TEST_F(PillowingTest, assignCells_IxJxI_3D_nonmanifold)
{
    const int nb_i = 30;
    const int nb_j = 10;

    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_IxJxI_3D_nonmanifold(&mesh, &fp, nb_i, nb_j);


    kmds::Variable<std::uintptr_t>* v = mesh.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr), kmds::KMDS_NODE, "geomEntity");
    gmds::cad::FACManager* facGeomManager = new gmds::cad::FACManager;

    elg3d::BoundingBoxGeomAssociation_init_3D(&mesh, facGeomManager, v);


    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);

    // connectivities
    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2R();

    assignCellsCorrection_3D(&mesh, c_N2R, &fp, &ma);


    // store cells midpoints
    kmds::TSize nbCells = mesh.getNbRegions();
    kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", nbCells);
    mesh.getRegionIDs_dummy(&cellIDs);
    kmds::Variable<gmds::math::Point>* varMidpoints = mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_REGION, "midpoint");
    elg3d::Tools_computeMidPointVar_3D(&cellIDs, &mesh, varMidpoints);


    // gradient computation
    kmds::Connectivity* c_R2R_byN = mesh.createConnectivity(kmds::R2R_byN);
    ch.buildR2R_byN(c_N2R);

    kmds::Variable<elg3d::MaterialGradientComputation_midcellGradients>* varMidcellgrads =
            mesh.createVariable<elg3d::MaterialGradientComputation_midcellGradients>(elg3d::MaterialGradientComputation_MIDCELLGRADIENTS_NULL, kmds::KMDS_REGION, "midcellgrads");
    elg3d::MaterialGradientComputation_leastsquare_3D(&cellIDs, c_R2R_byN, &fp, &ma, varMidpoints, varMidcellgrads);


    // compute interface nodes new position
    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);

    ch.buildFandF2R();
    kmds::Connectivity* c_F2R = mesh.getConnectivity(kmds::F2R);
    kmds::Connectivity* c_N2F = mesh.createConnectivity(kmds::N2F);
    ch.buildN2F();

    kmds::Variable<gmds::math::Point>* varNewPos =
            mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_NODE, "newPos");

    elg3d::InterfaceNodesPos_computeNodesNewPos_3D(&nodesInterfaces,
                                                   &mesh,
                                                   c_N2F,
                                                   c_F2R,
                                                   &ma,
                                                   &fp,
                                                   varMidpoints,
                                                   varMidcellgrads,
                                                   varNewPos);


    std::cout<<0<<" "<<(*varNewPos)[0]<<std::endl;
    std::cout<<2<<" "<<(*varNewPos)[2]<<std::endl;
    std::cout<<17<<" "<<(*varNewPos)[17]<<std::endl;
    std::cout<<18<<" "<<(*varNewPos)[18]<<std::endl;


    elg3d::moveToNewPos_basicMove(&nodesInterfaces,
                                  &mesh,
                                  varNewPos,
                                  v);

    kmds::GrowingView<kmds::TCellID> facesInterfaces("FACES_ON_INTERFACES", mesh.getNbFaces());
    elg3d::MaterialInterfaces_getFaceOnInterfaces(&mesh, c_F2R, &ma, &facesInterfaces);


    kmds::VTKWriter<kmds::Mesh> w(mesh);
    w.write("PillowingTest_before_IxJxI_3D", kmds::F|kmds::R);

    elg3d::pillow_3D(&facesInterfaces,
                     &nodesInterfaces,
                     &mesh,
                     c_F2R,
                     c_N2R,
                     &ma,
                     v);

    w.write("PillowingTest_after_IxJxI_3D", kmds::F|kmds::R);


//    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
    mesh.deleteConnectivity(kmds::N2R);
    c_N2R = mesh.createConnectivity(kmds::N2R);
    ch.buildN2R();

    ch.buildN2N(kmds::R);
    kmds::Connectivity* c_N2N = mesh.getConnectivity(kmds::N2N);

    elg3d::smartLaplacian_3D(10,
                             &mesh,
                             c_N2R,
                             c_N2N,
                             v,
                             &nodesInterfaces);

    std::cout<<"smartLaplacian_3D after exit"<<std::endl;
    w.write("PillowingTest_smoothed_IxJxI_3D", kmds::R);
}
/*----------------------------------------------------------------------------*/
TEST_F(PillowingTest, internalFormat_3x3x3_3D) {
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_internalFormat(&mesh, &fp, "Elg3D/test/Samples/internalFormat_3x3x3_3D.txt");


    kmds::Variable<std::uintptr_t> *v = mesh.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr),
                                                                            kmds::KMDS_NODE, "geomEntity");
    gmds::cad::FACManager *facGeomManager = new gmds::cad::FACManager;

    elg3d::BoundingBoxGeomAssociation_init_3D(&mesh, facGeomManager, v);


    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);

    elg3d::Tools_write_lite_3D(&mesh, &ma, "PillowingTest_internalFormat_3x3x3_3D");

    // connectivities
    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2R();

//    assignCellsCorrection_3D(&mesh, c_N2R, &fp, &ma);


    // store cells midpoints
    kmds::TSize nbCells = mesh.getNbRegions();
    kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", nbCells);
    mesh.getRegionIDs_dummy(&cellIDs);
    kmds::Variable<gmds::math::Point>* varMidpoints = mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_REGION, "midpoint");
    elg3d::Tools_computeMidPointVar_3D(&cellIDs, &mesh, varMidpoints);


    // gradient computation
    kmds::Connectivity* c_R2R_byN = mesh.createConnectivity(kmds::R2R_byN);
    ch.buildR2R_byN(c_N2R);

    kmds::Variable<elg3d::MaterialGradientComputation_midcellGradients>* varMidcellgrads =
            mesh.createVariable<elg3d::MaterialGradientComputation_midcellGradients>(elg3d::MaterialGradientComputation_MIDCELLGRADIENTS_NULL, kmds::KMDS_REGION, "midcellgrads");
    elg3d::MaterialGradientComputation_leastsquare_3D(&cellIDs, c_R2R_byN, &fp, &ma, varMidpoints, varMidcellgrads);


    // compute interface nodes new position
    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);

    ch.buildFandF2R();
    kmds::Connectivity* c_F2R = mesh.getConnectivity(kmds::F2R);
    kmds::Connectivity* c_N2F = mesh.createConnectivity(kmds::N2F);
    ch.buildN2F();

    kmds::Variable<gmds::math::Point>* varNewPos =
            mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_NODE, "newPos");

    elg3d::InterfaceNodesPos_computeNodesNewPos_3D(&nodesInterfaces,
                                                   &mesh,
                                                   c_N2F,
                                                   c_F2R,
                                                   &ma,
                                                   &fp,
                                                   varMidpoints,
                                                   varMidcellgrads,
                                                   varNewPos);



    elg3d::moveToNewPos_basicMove(&nodesInterfaces,
                                  &mesh,
                                  varNewPos,
                                  v);

    kmds::GrowingView<kmds::TCellID> facesInterfaces("FACES_ON_INTERFACES", mesh.getNbFaces());
    elg3d::MaterialInterfaces_getFaceOnInterfaces(&mesh, c_F2R, &ma, &facesInterfaces);

    elg3d::Tools_write_interfaces_3D(&facesInterfaces,
                                     &nodesInterfaces,
                                     &mesh,
                                     c_F2R,
                                     &ma,
                                     "PillowingTest_internalFormat_interface_before_move_3x3x3_3D");


    //w.write("PillowingTest_before_IxJxI_3D", kmds::F|kmds::R);
    elg3d::Tools_write_lite_3D(&mesh, &ma, "PillowingTest_internalFormat_before_pillow_3x3x3_3D");

    elg3d::pillow_3D(&facesInterfaces,
                     &nodesInterfaces,
                     &mesh,
                     c_F2R,
                     c_N2R,
                     &ma,
                     v);

    elg3d::Tools_write_lite_3D(&mesh, &ma, "PillowingTest_internalFormat_after_pillow_3x3x3_3D");
//
////    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
//    mesh.deleteConnectivity(kmds::N2R);
//    c_N2R = mesh.createConnectivity(kmds::N2R);
//    ch.buildN2R();
//
//    ch.buildN2N(kmds::R);
//    kmds::Connectivity* c_N2N = mesh.getConnectivity(kmds::N2N);
//
//    elg3d::smartLaplacian_3D(10,
//                             &mesh,
//                             c_N2R,
//                             c_N2N,
//                             v,
//                             &nodesInterfaces);
//
//    std::cout<<"smartLaplacian_3D after exit"<<std::endl;
//    elg3d::Tools_write_lite_3D(&mesh, &ma, "PillowingTest_internalFormat_after_smooth_3x3x3_3D");
}
/*----------------------------------------------------------------------------*/
TEST_F(PillowingTest, execute_2D)
{
    kmds::Mesh mesh;

    const kmds::TCoord minXYZ[3] = {0.,0.,0.};
    const kmds::TCoord maxXYZ[3] = {5.,5.,0.};

    kmds::InitTools_createGrid_2D(&mesh, minXYZ, maxXYZ, 5, 5);



    kmds::Variable<std::uintptr_t>* v = mesh.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr), kmds::KMDS_NODE, "geomEntity");
    gmds::cad::FACManager* facGeomManager = new gmds::cad::FACManager;

    elg3d::BoundingBoxGeomAssociation_init_2D(&mesh, facGeomManager, v);


    elg3d::MaterialAssignment ma(mesh.getFaceCapacity());
    ma.createMaterial("matA");
    ma.createMaterial("matB");
    ma.setMaterial(0,0);
    ma.setMaterial(0,1);
    ma.setMaterial(0,2);
    ma.setMaterial(0,3);
    ma.setMaterial(1,4);
    ma.setMaterial(0,5);
    ma.setMaterial(0,6);
    ma.setMaterial(0,7);
    ma.setMaterial(0,8);
    ma.setMaterial(1,9);
    ma.setMaterial(1,10);
    ma.setMaterial(1,11);
    ma.setMaterial(1,12);
    ma.setMaterial(0,13);
    ma.setMaterial(1,14);
    ma.setMaterial(1,15);
    ma.setMaterial(1,16);
    ma.setMaterial(1,17);
    ma.setMaterial(0,18);
    ma.setMaterial(1,19);
    ma.setMaterial(1,20);
    ma.setMaterial(1,21);
    ma.setMaterial(1,22);
    ma.setMaterial(1,23);
    ma.setMaterial(1,24);

    elg3d::Tools_write_lite_2D(&mesh, &ma, "test_pillow_execute_ma");

    kmds::Variable<gmds::math::Point>* varNewPos =
            mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_NODE, "newPos");
    kmds::Variable<gmds::math::Point>* varRefPos =
            mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_NODE, "refPos");


    elg3d::Tools_storeNodePos_xD(&mesh, varRefPos);
    elg3d::Tools_storeNodePos_xD(&mesh, varNewPos);

    (*varNewPos)[28] = gmds::math::Point(4.2, 4. , 0.);
    (*varNewPos)[27] = gmds::math::Point(3.8, 3.6, 0.);
    (*varNewPos)[21] = gmds::math::Point(3. , 3.3, 0.);

    (*varNewPos)[14] = gmds::math::Point(1.5, 2.4, 0.);

    (*varNewPos)[13] = gmds::math::Point(1.5, 1., 0.);
    (*varNewPos)[12] = gmds::math::Point(1.5, 0. , 0.);

    elg3d::Tools_loadNodePos_xD(&mesh, varNewPos);

    kmds::Variable <double> *varCellQuality =
            mesh.createVariable<double>(-HUGE_VALF, kmds::KMDS_FACE, "varCellQuality");
    elg3d::Tools_computeScaledJacobian_2D(&mesh,varCellQuality);

    elg3d::Tools_write_lite_2D(&mesh, &ma, "test_pillow_execute_newPos");

    elg3d::Tools_loadNodePos_xD(&mesh, varRefPos);


    elg3d::Tools_write_lite_2D(&mesh, &ma, "test_pillow_execute_poyop");

    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    nodesInterfaces.push_back(12);
    nodesInterfaces.push_back(13);
    nodesInterfaces.push_back(14);
    nodesInterfaces.push_back(15);
    nodesInterfaces.push_back(21);
    nodesInterfaces.push_back(27);
    nodesInterfaces.push_back(28);
    nodesInterfaces.push_back(22);
    nodesInterfaces.push_back(16);
    nodesInterfaces.push_back(10);
    nodesInterfaces.push_back(4);



    kmds::Connectivity* c_N2F = mesh.createConnectivity(kmds::N2F);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2F();


    kmds::Variable <int> *varNbMoves =
            mesh.createVariable<int>(0, kmds::KMDS_NODE, "varNbMoves");

    elg3d::moveToNewPos_noBadMove_2D(16,
                                     elg3d::Parameters::quality_threshold,
                                     &nodesInterfaces,
                                     &mesh,
                                     c_N2F,
                                     varNewPos,
                                     v,
                                     varNbMoves);

    elg3d::Tools_computeScaledJacobian_2D(&mesh,varCellQuality);
    elg3d::Tools_write_lite_2D(&mesh, &ma, "test_pillow_execute_newPos_nobadmove");

    kmds::Variable <double> *varDist =
            mesh.createVariable<double>(0, kmds::KMDS_NODE, "varDist");
    elg3d::Tools_computeDistance_xD(&nodesInterfaces, &mesh, varNewPos, varDist);

    kmds::Connectivity* c_E2F = mesh.createConnectivity(kmds::E2F);
    ch.buildEandE2F();
    kmds::Connectivity* c_N2E = mesh.createConnectivity(kmds::N2E);
    ch.buildN2E();


    elg3d::Parameters::pillow_node_seed_range = 3;
//    elg3d::pillow_execute_2D(
//            &nodesInterfaces,
//            &mesh,
//            c_E2F,
//            c_N2E,
//            c_N2F,
//            &ma,
//            v,
//            varDist);
    elg3d::pillow_execute_propagate_2D(
            &nodesInterfaces,
            &mesh,
            c_E2F,
            c_N2E,
            c_N2F,
            &ma,
            v,
            varDist);
    elg3d::Tools_computeScaledJacobian_2D(&mesh,varCellQuality);

    elg3d::Tools_write_lite_2D(&mesh, &ma, "test_pillow_execute_localized");


    mesh.deleteConnectivity(kmds::N2F);
    c_N2F = mesh.createConnectivity(kmds::N2F);
    ch.buildN2F();

    ch.buildN2N(kmds::F);
    kmds::Connectivity* c_N2N = mesh.getConnectivity(kmds::N2N);

    elg3d::smartLaplacian_2D(10,
                             &mesh,
                             c_N2F,
                             c_N2N,
                             v,
                             &nodesInterfaces);

    double minscaledJacobian = elg3d::Tools_computeScaledJacobian_2D(&mesh,varCellQuality);
    elg3d::Tools_write_lite_2D(&mesh, &ma, "test_pillow_execute_smooth");

    elg3d::moveToNewPos_noBadMove_2D(16,
                                     elg3d::Parameters::quality_threshold,
                                     &nodesInterfaces,
                                     &mesh,
                                     c_N2F,
                                     varNewPos,
                                     v,
                                     varNbMoves);

    elg3d::Tools_computeScaledJacobian_2D(&mesh,varCellQuality);
    elg3d::Tools_write_lite_2D(&mesh, &ma, "test_pillow_execute_poyop");

    EXPECT_GT(minscaledJacobian, elg3d::Parameters::quality_threshold);
}

/*----------------------------------------------------------------------------*/
TEST_F(PillowingTest, execute_propagate_3D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_internalFormat(&mesh, &fp, "Elg3D/test/Samples/internalFormat_3x3x3_3D.txt");


    kmds::Variable<std::uintptr_t> *v = mesh.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr),
                                                                            kmds::KMDS_NODE, "geomEntity");
    gmds::cad::FACManager *facGeomManager = new gmds::cad::FACManager;

    elg3d::BoundingBoxGeomAssociation_init_3D(&mesh, facGeomManager, v);


    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);

    elg3d::Tools_write_lite_3D(&mesh, &ma, "PillowingTest_internalFormat_3x3x3_3D");

    // connectivities
    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2R();

//    assignCellsCorrection_3D(&mesh, c_N2R, &fp, &ma);


    // store cells midpoints
    kmds::TSize nbCells = mesh.getNbRegions();
    kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", nbCells);
    mesh.getRegionIDs_dummy(&cellIDs);
    kmds::Variable<gmds::math::Point>* varMidpoints = mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_REGION, "midpoint");
    elg3d::Tools_computeMidPointVar_3D(&cellIDs, &mesh, varMidpoints);


    // gradient computation
    kmds::Connectivity* c_R2R_byN = mesh.createConnectivity(kmds::R2R_byN);
    ch.buildR2R_byN(c_N2R);

    kmds::Variable<elg3d::MaterialGradientComputation_midcellGradients>* varMidcellgrads =
            mesh.createVariable<elg3d::MaterialGradientComputation_midcellGradients>(elg3d::MaterialGradientComputation_MIDCELLGRADIENTS_NULL, kmds::KMDS_REGION, "midcellgrads");
    elg3d::MaterialGradientComputation_leastsquare_3D(&cellIDs, c_R2R_byN, &fp, &ma, varMidpoints, varMidcellgrads);


    // compute interface nodes new position
    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);

    ch.buildFandF2R();
    kmds::Connectivity* c_F2R = mesh.getConnectivity(kmds::F2R);
    kmds::Connectivity* c_N2F = mesh.createConnectivity(kmds::N2F);
    ch.buildN2F();

    kmds::Variable<gmds::math::Point>* varNewPos =
            mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_NODE, "newPos");
    kmds::Variable<gmds::math::Point>* varRefPos =
            mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_NODE, "refPos");

    elg3d::Tools_storeNodePos_xD(&mesh, varRefPos);
    elg3d::Tools_storeNodePos_xD(&mesh, varNewPos);


    elg3d::InterfaceNodesPos_computeNodesNewPos_3D(&nodesInterfaces,
                                                   &mesh,
                                                   c_N2F,
                                                   c_F2R,
                                                   &ma,
                                                   &fp,
                                                   varMidpoints,
                                                   varMidcellgrads,
                                                   varNewPos);

    elg3d::moveToNewPos_basicMove(&nodesInterfaces,
                                  &mesh,
                                  varNewPos,
                                  v);

    elg3d::Tools_storeNodePos_xD(&mesh, varNewPos);
    elg3d::Tools_loadNodePos_xD(&mesh, varRefPos);

    kmds::Variable <int> *varNbMoves =
            mesh.createVariable<int>(0, kmds::KMDS_NODE, "varNbMoves");

    elg3d::moveToNewPos_noBadMove_3D(16,
                                     elg3d::Parameters::quality_threshold,
                                     &nodesInterfaces,
                                     &mesh,
                                     c_N2R,
                                     varNewPos,
                                     v,
                                     varNbMoves);

    // testing the noBadMove
    kmds::Variable <double> *varCellQuality =
            mesh.createVariable<double>(-HUGE_VALF, kmds::KMDS_REGION, "varCellQuality");
    elg3d::Tools_computeScaledJacobian_3D(&mesh,varCellQuality);

    kmds::Variable <double> *varDist =
            mesh.createVariable<double>(-HUGE_VALF, kmds::KMDS_NODE, "varDist");
    elg3d::Tools_computeDistance_xD(&nodesInterfaces, &mesh, varNewPos, varDist);


    elg3d::Tools_write_lite_3D(&mesh, &ma, "AssignCellsTest_assignCells_I_3D_afternoBadMove");





    kmds::GrowingView<kmds::TCellID> facesInterfaces("FACES_ON_INTERFACES", mesh.getNbFaces());
    elg3d::MaterialInterfaces_getFaceOnInterfaces(&mesh, c_F2R, &ma, &facesInterfaces);

    elg3d::Tools_write_interfaces_3D(&facesInterfaces,
                                     &nodesInterfaces,
                                     &mesh,
                                     c_F2R,
                                     &ma,
                                     "PillowingTest_internalFormat_interface_before_move_3x3x3_3D");


    //w.write("PillowingTest_before_IxJxI_3D", kmds::F|kmds::R);
    elg3d::Tools_write_lite_3D(&mesh, &ma, "PillowingTest_internalFormat_before_pillow_3x3x3_3D");

//    elg3d::pillow_3D(&facesInterfaces,
//                     &nodesInterfaces,
//                     &mesh,
//                     c_F2R,
//                     c_N2R,
//                     &ma,
//                     v);

    (*varDist)[33] = 1.;

    elg3d::pillow_execute_propagate_3D(&nodesInterfaces,
                                       &mesh,
                                       c_F2R,
                                       c_N2F,
                                       c_N2R,
                                       &ma,
                                       v,
                                       varDist);

    double minscaledJacobian = elg3d::Tools_computeScaledJacobian_3D(&mesh,varCellQuality);
    elg3d::Tools_computeDistance_xD(&nodesInterfaces, &mesh, varNewPos, varDist);

    elg3d::Tools_write_lite_3D(&mesh, &ma, "PillowingTest_internalFormat_after_pillow_3x3x3_3D");




    elg3d::Tools_write_lite_3D(&mesh, &ma, "test_pillow_execute_smooth");

    EXPECT_GT(minscaledJacobian, elg3d::Parameters::quality_threshold);
}

/*----------------------------------------------------------------------------*/