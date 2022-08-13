/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    main_exodusreader2d.cpp
 *  \author  legoff
 *  \date    03/19/2019
 */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <Eigen/Eigen>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/cadfac/FACManager.h>
//#include <GMDS/IG/IG.h>

#include <KM/DS/Mesh.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include <ELG3D/ALGOCMPNT/AssignCells.h>
#include <ELG3D/ALGOCMPNT/BoundingBoxGeomAssociation.h>
#include <ELG3D/ALGOCMPNT/FracPresEnforcement.h>
#include <ELG3D/ALGOCMPNT/InitData.h>
#include <ELG3D/ALGOCMPNT/InterfaceNodesPos.h>
#include <ELG3D/ALGOCMPNT/ManifoldDetection.h>
#include <ELG3D/ALGOCMPNT/MaterialGradientComputation.h>
#include <ELG3D/ALGOCMPNT/MaterialInterfaces.h>
#include <ELG3D/ALGOCMPNT/MoveToNewPos.h>
#include <ELG3D/ALGOCMPNT/Pillow.h>
#include <ELG3D/ALGOCMPNT/SmartLaplacian.h>
#include <ELG3D/ALGOCMPNT/Tools.h>
#include <ELG3D/DATACMPNT/FracPres.h>
#include <ELG3D/DATACMPNT/MaterialAssignment.h>
/*----------------------------------------------------------------------------*/

int main(int argc, char* argv[])
{
    int num_threads = -1;

    std::istringstream iss( argv[1] );
    int val;

    if (iss >> val)
    {
        num_threads = val;
    } else {
        std::cerr<<"could not convert number of threads argument."<<std::endl;
        exit(-1);
    }

    bool writeOutputFiles = false;
    std::istringstream iss_2(argv[2]);
    iss_2 >> writeOutputFiles;

    bool isFromFile = false;
    std::istringstream iss_3(argv[3]);
    iss_3 >> isFromFile;

    std::string case_filename;
    bool fill_void;
    if (isFromFile) {
        std::istringstream iss_4(argv[4]);
        iss_4 >> case_filename;
        std::istringstream iss_5(argv[5]);
        iss_5 >> fill_void;
    }


//    bool isWriteLite = true;
    bool isWriteLite = false;

    // necessary in case eigen is called by multiple threads
    Eigen::initParallel();

    Kokkos::InitArguments kargs;
    kargs.num_threads = num_threads;
    Kokkos::initialize(kargs);

    std::cout<<"num_threads "<<num_threads<<std::endl;


    Kokkos::Timer timer_tot;
    timer_tot.reset();
    Kokkos::Timer timer;
    timer.reset();


    kmds::Mesh mesh;
    elg3d::FracPres fp;

//    const double xyz_min[3] = {0., 0., 0.};
//    const double xyz_max[3] = {24., 13., 0.};
//    const double xyz_min[3] = {-2., -2., 0.};
//    const double xyz_max[3] = {12., 12., 0.};
//    const int ni = 12;
//    const int nj = 12;
//    const double xyz_min[3] = {-2., -6, 0.};
//    const double xyz_max[3] = {6., 2., 0.};
//    const int ni = 40;
//    const int nj = 40;
    const double xyz_min[3] = {0., 0., 0.};
    const double xyz_max[3] = {7., 7., 7.};
    const int ni = 12;
    const int nj = 12;

    elg3d::initData_fromExodus_2D(&mesh, &fp, case_filename, xyz_min, xyz_max, ni, nj, fill_void);


    // Enforce material selection
    elg3d::FracPres fp_enforce;
    {
        kmds::TSize nbCells = mesh.getNbFaces();
        kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", nbCells);
        mesh.getFaceIDs_dummy(&cellIDs);

        std::vector<std::string> list_of_mat;
        list_of_mat.push_back("block_1");
//        list_of_mat.push_back("block_2");

        elg3d::FracPresEnforcement_maintain_xD(&cellIDs,
                                               &fp,
                                               &fp_enforce,
                                               list_of_mat);
    }


////    elg3d::initData_internalFormat(&mesh, &fp, "Elg3D/test/Samples/internalFormat_3x3x3_3D.txt");
//    elg3d::initData_internalFormat(&mesh, &fp, "Elg3D/test/Samples/internalFormat_5x5_2D.txt");
////    elg3d::initData_I_2D_nonmanifold(&mesh, &fp, 2000);
//


    kmds::Variable<std::uintptr_t> *v = mesh.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr),
                                                                            kmds::KMDS_NODE, "geomEntity");
    gmds::cad::FACManager *facGeomManager = new gmds::cad::FACManager;

    elg3d::BoundingBoxGeomAssociation_init_2D(&mesh, facGeomManager, v);


    elg3d::MaterialAssignment ma(mesh.getFaceCapacity());

    timer.reset();
//    elg3d::assignCellsMajorityCriteria_2D(&mesh, &fp, &ma);
    elg3d::assignCellsMajorityCriteria_2D(&mesh, &fp_enforce, &ma);
    std::cout<<"timer assignCellsMajorityCriteria "<<timer.seconds()<<std::endl;

    if (writeOutputFiles) {
        if(isWriteLite) {
            elg3d::Tools_write_lite_2D(&mesh, &ma, "AssignCellsTest_assignCells_I_2D_majority");
        } else {
            elg3d::Tools_write_2D(&mesh, &fp, &ma, "AssignCellsTest_assignCells_I_2D_majority");
        }

    }

    // connectivities
    timer.reset();
    kmds::Connectivity* c_N2F = mesh.createConnectivity(kmds::N2F);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2F();
    std::cout<<"timer buildN2R_firsttime "<<timer.seconds()<<std::endl;


    timer.reset();
    elg3d::assignCellsCorrection_2D(&mesh, c_N2F, &fp, &ma);
    std::cout<<"timer assignCellsCorrection "<<timer.seconds()<<std::endl;

    if (writeOutputFiles) {
        if(isWriteLite) {
            elg3d::Tools_write_lite_2D(&mesh, &ma, "AssignCellsTest_assignCells_I_2D_correction");
        } else {
            elg3d::Tools_write_2D(&mesh, &fp, &ma, "AssignCellsTest_assignCells_I_2D_correction");
        }
    }


    // store cells midpoints
    timer.reset();
    kmds::TSize nbCells = mesh.getNbFaces();
    kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", nbCells);
    mesh.getFaceIDs_dummy(&cellIDs);
    kmds::Variable<gmds::math::Point>* varMidpoints = mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_FACE, "midpoint");
    elg3d::Tools_computeMidPointVar_2D(&cellIDs, &mesh, varMidpoints);
    std::cout<<"timer midpoints "<<timer.seconds()<<std::endl;

    // gradient computation
    timer.reset();
    kmds::Connectivity* c_F2F_byN = mesh.createConnectivity(kmds::F2F_byN);
    ch.buildF2F_byN(c_N2F);
    std::cout<<"timer buildF2F_byN "<<timer.seconds()<<std::endl;

    timer.reset();
    kmds::Variable<elg3d::MaterialGradientComputation_midcellGradients>* varMidcellgrads =
            mesh.createVariable<elg3d::MaterialGradientComputation_midcellGradients>(elg3d::MaterialGradientComputation_MIDCELLGRADIENTS_NULL, kmds::KMDS_FACE, "midcellgrads");
    elg3d::MaterialGradientComputation_leastsquare_2D(&cellIDs, c_F2F_byN, &fp, &ma, varMidpoints, varMidcellgrads);
    std::cout<<"timer GradientComputation "<<timer.seconds()<<std::endl;

    if (writeOutputFiles) {
        if(isWriteLite) {

        } else {
            elg3d::Tools_write_2D(&mesh, &fp, &ma, varMidcellgrads, "AssignCellsTest_gradient_I_2D");
        }
    }

    // compute interface nodes new position
    timer.reset();
    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2F, &ma, &nodesInterfaces);
    std::cout<<"timer getNodeOnInterfaces "<<timer.seconds()<<std::endl;

    timer.reset();
    ch.buildEandE2F();
    std::cout<<"timer buildEandE2F "<<timer.seconds()<<std::endl;
    timer.reset();
    kmds::Connectivity* c_E2F = mesh.getConnectivity(kmds::E2F);
    kmds::Connectivity* c_N2E = mesh.createConnectivity(kmds::N2E);
    ch.buildN2E();
    std::cout<<"timer buildN2E "<<timer.seconds()<<std::endl;

    kmds::Variable<gmds::math::Point>* varNewPos =
            mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_NODE, "newPos");

    timer.reset();
    elg3d::InterfaceNodesPos_computeNodesNewPos_2D(&nodesInterfaces,
                                                   &mesh,
                                                   c_N2E,
                                                   c_E2F,
                                                   &ma,
                                                   &fp,
                                                   varMidpoints,
                                                   varMidcellgrads,
                                                   varNewPos);
    std::cout<<"timer computeNodesNewPos "<<timer.seconds()<<std::endl;


    timer.reset();
    elg3d::moveToNewPos_basicMove(&nodesInterfaces,
                                  &mesh,
                                  varNewPos,
                                  v);
    std::cout<<"timer basicMove "<<timer.seconds()<<std::endl;


    timer.reset();
    kmds::GrowingView<kmds::TCellID> facesInterfaces("FACES_ON_INTERFACES", mesh.getNbEdges());
    elg3d::MaterialInterfaces_getEdgeOnInterfaces(&mesh, c_E2F, &ma, &facesInterfaces);
    std::cout<<"timer getEdgeOnInterfaces "<<timer.seconds()<<std::endl;


//    elg3d::Tools_write_interfaces_2D(&facesInterfaces,
//                                     &nodesInterfaces,
//                                     &mesh,
//                                     c_E2F,
//                                     &ma,
//                                     "PillowingTest_internalFormat_interface_before_move_3x3x3_3D");

    if (writeOutputFiles) {
        if(isWriteLite) {
            elg3d::Tools_write_lite_2D(&mesh, &ma, "AssignCellsTest_assignCells_I_2D_before_pillow");
        } else {
            elg3d::Tools_write_2D(&mesh, &fp, &ma, "AssignCellsTest_assignCells_I_2D_before_pillow");
        }
    }


    timer.reset();
    elg3d::pillow_2D(&facesInterfaces,
                     &nodesInterfaces,
                     &mesh,
                     c_E2F,
                     c_N2F,
                     &ma,
                     v);
    std::cout<<"timer pillow "<<timer.seconds()<<std::endl;


    if (writeOutputFiles) {
        if(isWriteLite) {
            elg3d::Tools_write_lite_2D(&mesh, &ma, "AssignCellsTest_assignCells_I_2D_after_pillow");
        } else {
            elg3d::Tools_write_2D(&mesh, &fp, &ma, "AssignCellsTest_assignCells_I_2D_after_pillow");
        }
    }


//    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
    mesh.deleteConnectivity(kmds::N2F);

    timer.reset();
    c_N2F = mesh.createConnectivity(kmds::N2F);
    ch.buildN2F();
    std::cout<<"timer buildN2F_secondtime "<<timer.seconds()<<std::endl;

    timer.reset();
    ch.buildN2N(kmds::F);
    std::cout<<"timer buildN2N "<<timer.seconds()<<std::endl;
    kmds::Connectivity* c_N2N = mesh.getConnectivity(kmds::N2N);


    timer.reset();
    elg3d::smartLaplacian_2D(10,
                             &mesh,
                             c_N2F,
                             c_N2N,
                             v,
                             &nodesInterfaces);
    std::cout<<"timer smartLaplacian "<<timer.seconds()<<std::endl;


//    std::cout<<"smartLaplacian_3D after exit"<<std::endl;

    if (writeOutputFiles) {
        if(isWriteLite) {
            elg3d::Tools_write_lite_2D(&mesh, &ma, "AssignCellsTest_assignCells_I_2D_after_smooth");
        } else {
            elg3d::Tools_write_2D(&mesh, &fp, &ma, "AssignCellsTest_assignCells_I_2D_after_smooth");
        }
    }


    std::cout<<"timer_tot "<<num_threads<<" "<<timer_tot.seconds()<<std::endl;

    Kokkos::finalize();
    return 0;
}