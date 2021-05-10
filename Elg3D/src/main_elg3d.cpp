/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    main_elg3d.cpp
 *  \author  legoff
 *  \date    07/04/2018
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
#include <gmds/cad/FACManager.h>
#include <KM/DS/Mesh.h>

#include <gmds/ig/Mesh.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/IGMeshIOService.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include <ELG3D/ALGOCMPNT/AssignCells.h>
#include <ELG3D/ALGOCMPNT/BoundingBoxGeomAssociation.h>
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
#include <ELG3D/ALGOCMPNT/ExtractGeomModel.h>

#include <ELG3D/ALGOCMPNT/FracPresEnforcement.h>

/*----------------------------------------------------------------------------*/

int main(int argc, char *argv[]) {
    int num_threads = -1;

    // nbthreads writeoutput fromfile filename deducevoid

    std::istringstream iss(argv[1]);
    int val;

    if (iss >> val) {
        num_threads = val;
    } else {
        std::cerr << "could not convert number of threads argument." << std::endl;
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

    // necessary in case eigen is called by multiple threads
    Eigen::initParallel();

    Kokkos::InitArguments kargs;
    kargs.num_threads = num_threads;
    Kokkos::initialize(kargs);

    std::cout << "num_threads " << num_threads << std::endl;

    Kokkos::Timer timer_tot;
    timer_tot.reset();

    Kokkos::Timer timer;
    timer.reset();


    kmds::Mesh mesh;
    elg3d::FracPres fp;

//    if (isFromFile) {
////        elg3d::Tools_read_fracpres_vtk_3D(&mesh, &fp, case_filename, fill_void);
//        elg3d::initData_internalFormat(&mesh, &fp, case_filename);
//    } else {
////        elg3d::initData_IxJxI_3D_nonmanifold(&mesh, &fp, 10, 10);
    const double xyz_min[3] = {-15., -15., -15.};
    const double xyz_max[3] = {15., 15., 15.};
//        elg3d::initData_fromTXT_3D(&mesh, &fp, "volfrac.txt", xyz_min, xyz_max, 30, 30, 30, 2);
//std::cout<<"case_filename "<<case_filename<<std::endl;
//elg3d::FracPres fp_tmp;
    elg3d::initData_fromExodus_3D(&mesh, &fp, case_filename, xyz_min, xyz_max, 10, 10, 10, 1);

//    kmds::GrowingView <kmds::TCellID> cellIDs_tmp("CELLS", mesh.getNbRegions());
//    mesh.getRegionIDs_dummy(&cellIDs_tmp);
    std::vector<std::string> names;
    names.push_back("block_0");
//    names.push_back();
//    elg3d::FracPresEnforcement_maintain_xD(&cellIDs_tmp, &fp_tmp, &fp, names);
//    std::map<int, std::string> ml_tmp = fp_tmp.getMaterialList();
//    std::map<int, std::string> ml = fp.getMaterialList();
//    for(auto m: ml_tmp) {
//        std::cout<<"mat_tmp "<<m.first<<" "<<m.second<<std::endl;
//    }
//    for(auto m: ml) {
//        std::cout<<"mat "<<m.first<<" "<<m.second<<std::endl;
//    }
//    }


//    elg3d::initData_internalFormat(&mesh, &fp, "Elg3D/test/Samples/internalFormat_3x3x3_3D.txt");
//    elg3d::initData_IxJxI_3D_nonmanifold(&mesh, &fp, 100, 100);
//

    kmds::Variable<std::uintptr_t> *v = mesh.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr),
                                                                            kmds::KMDS_NODE, "geomEntity");
    gmds::cad::FACManager *facGeomManager = new gmds::cad::FACManager;

    elg3d::BoundingBoxGeomAssociation_init_3D(&mesh, facGeomManager, v);


    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    timer.reset();
//    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);
    elg3d::assignCellsMajorityCriteria_maintain_3D(&mesh, &fp, &ma, names);
    std::cout << "timer assignCellsMajorityCriteria " << timer.seconds() << std::endl;

//    elg3d::Tools_write_lite_3D(&mesh, &ma, "PillowingTest_internalFormat_3x3x3_3D");
    if (writeOutputFiles) {
        elg3d::Tools_write_3D(&mesh, &fp, &ma, "elg3d_assignCellsMajorityCriteria");
    }

    // connectivities
    timer.reset();
    kmds::Connectivity *c_N2R = mesh.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2R();
    std::cout << "timer buildN2R_firsttime " << timer.seconds() << std::endl;

    timer.reset();
//    elg3d::assignCellsCorrection_3D(&mesh, c_N2R, &fp, &ma);
    std::cout << "timer assignCellsCorrection " << timer.seconds() << std::endl;

    if (writeOutputFiles) {
        elg3d::Tools_write_3D(&mesh, &fp, &ma, "elg3d_assignCellsCorrection");
    }

    // store cells midpoints
    timer.reset();
    kmds::TSize nbCells = mesh.getNbRegions();
    kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", nbCells);
    mesh.getRegionIDs_dummy(&cellIDs);
    kmds::Variable<gmds::math::Point> *varMidpoints = mesh.createVariable<gmds::math::Point>(
            gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_REGION, "midpoint");
    elg3d::Tools_computeMidPointVar_3D(&cellIDs, &mesh, varMidpoints);
    std::cout << "timer midpoints " << timer.seconds() << std::endl;

    // gradient computation
    timer.reset();
    kmds::Connectivity *c_R2R_byN = mesh.createConnectivity(kmds::R2R_byN);
    ch.buildR2R_byN(c_N2R);
    std::cout << "timer buildR2R_byN " << timer.seconds() << std::endl;

    timer.reset();
    kmds::Variable<elg3d::MaterialGradientComputation_midcellGradients> *varMidcellgrads =
            mesh.createVariable<elg3d::MaterialGradientComputation_midcellGradients>(
                    elg3d::MaterialGradientComputation_MIDCELLGRADIENTS_NULL, kmds::KMDS_REGION, "midcellgrads");
    elg3d::MaterialGradientComputation_leastsquare_3D(&cellIDs, c_R2R_byN, &fp, &ma, varMidpoints, varMidcellgrads);
    std::cout << "timer GradientComputation " << timer.seconds() << std::endl;


    // compute interface nodes new position
    timer.reset();
    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);
    std::cout << "timer getNodeOnInterfaces " << timer.seconds() << std::endl;

    timer.reset();
    ch.buildFandF2R_variant_0();
//    ch.buildFandF2R();
    std::cout << "timer buildFandF2R " << timer.seconds() << std::endl;
    timer.reset();
    kmds::Connectivity *c_F2R = mesh.getConnectivity(kmds::F2R);
    kmds::Connectivity *c_N2F = mesh.createConnectivity(kmds::N2F);
    ch.buildN2F();
    std::cout << "timer buildN2F " << timer.seconds() << std::endl;

    kmds::Variable<gmds::math::Point> *varNewPos =
            mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF),
                                                   kmds::KMDS_NODE, "newPos");

    timer.reset();
    elg3d::InterfaceNodesPos_computeNodesNewPos_3D(&nodesInterfaces,
                                                   &mesh,
                                                   c_N2F,
                                                   c_F2R,
                                                   &ma,
                                                   &fp,
                                                   varMidpoints,
                                                   varMidcellgrads,
                                                   varNewPos);
    std::cout << "timer computeNodesNewPos " << timer.seconds() << std::endl;


    timer.reset();
    elg3d::moveToNewPos_basicMove(&nodesInterfaces,
                                  &mesh,
                                  varNewPos,
                                  v);
    std::cout << "timer basicMove " << timer.seconds() << std::endl;


    {
// write interface geom model
        kmds::Variable<std::uintptr_t> *geomassoc = mesh.createVariable<std::uintptr_t>(
                reinterpret_cast<std::uintptr_t>(nullptr),
                kmds::KMDS_NODE, "geomassoc");

        gmds::cad::FACManager geomModel;
        elg3d::extractGeomModel_extract_3D(&mesh, &ma, c_F2R, &geomModel, geomassoc);
        gmds::Mesh &meshGeom = geomModel.getMeshView();
        gmds::IGMeshIOService ioService(&meshGeom);
        gmds::VTKWriter vtkWriter(&ioService);
        vtkWriter.setCellOptions(gmds::N | gmds::F);
        vtkWriter.setDataOptions(gmds::N | gmds::F);
        vtkWriter.write("geommodel.vtk");

//    geomModel.write_oneSurfPerFile("onesurf_");
        geomModel.write_surfaces("interfaces.vtk");

        exit(0);
    }

    timer.reset();
    kmds::GrowingView<kmds::TCellID> facesInterfaces("FACES_ON_INTERFACES", mesh.getNbFaces());
    elg3d::MaterialInterfaces_getFaceOnInterfaces(&mesh, c_F2R, &ma, &facesInterfaces);
    std::cout << "timer getFaceOnInterfaces " << timer.seconds() << std::endl;

    if (writeOutputFiles) {
        elg3d::Tools_write_interfaces_3D(&facesInterfaces,
                                         &nodesInterfaces,
                                         &mesh,
                                         c_F2R,
                                         &ma,
                                         "PillowingTest_internalFormat_interface_before_move_3x3x3_3D");


        elg3d::Tools_write_lite_3D(&mesh, &ma, "PillowingTest_internalFormat_before_pillow_3x3x3_3D");
    }

    timer.reset();
    elg3d::pillow_3D(&facesInterfaces,
                     &nodesInterfaces,
                     &mesh,
                     c_F2R,
                     c_N2R,
                     &ma,
                     v);
    std::cout << "timer pillow " << timer.seconds() << std::endl;

    if (writeOutputFiles) {
        elg3d::Tools_write_lite_3D(&mesh, &ma, "PillowingTest_internalFormat_after_pillow_3x3x3_3D");
    }

//    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
    mesh.deleteConnectivity(kmds::N2R);

    timer.reset();
    c_N2R = mesh.createConnectivity(kmds::N2R);
    ch.buildN2R();
    std::cout << "timer buildN2R_secondtime " << timer.seconds() << std::endl;

    timer.reset();
    ch.buildN2N(kmds::R);
    std::cout << "timer buildN2N " << timer.seconds() << std::endl;
    kmds::Connectivity *c_N2N = mesh.getConnectivity(kmds::N2N);


    timer.reset();
    elg3d::smartLaplacian_3D(10,
                             &mesh,
                             c_N2R,
                             c_N2N,
                             v,
                             &nodesInterfaces);
    std::cout << "timer smartLaplacian " << timer.seconds() << std::endl;


    std::cout << "smartLaplacian_3D after exit" << std::endl;
    if (writeOutputFiles) {
        elg3d::Tools_write_lite_3D(&mesh, &ma, "PillowingTest_internalFormat_after_smooth_3x3x3_3D");
    }


    // geometric model extraction
    timer.reset();
    mesh.deleteConnectivity(kmds::F2R);
    mesh.removeAllFaces();
    c_F2R = mesh.createConnectivity(kmds::F2R);
    ch.buildFandF2R_variant_0();

    kmds::Variable<std::uintptr_t> *geomassoc = mesh.createVariable<std::uintptr_t>(
            reinterpret_cast<std::uintptr_t>(nullptr),
            kmds::KMDS_NODE, "geomassoc");

    gmds::cad::FACManager geomModel;
    elg3d::extractGeomModel_extract_3D(&mesh, &ma, c_F2R, &geomModel, geomassoc);
    gmds::Mesh &meshGeom = geomModel.getMeshView();
    gmds::IGMeshIOService ioService(&meshGeom);
    gmds::VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N | gmds::F);
    vtkWriter.setDataOptions(gmds::N | gmds::F);
    vtkWriter.write("geommodel.vtk");

//    gmds::VTKWriter<gmds::Mesh> vtkwriter(meshGeom);
//    vtkwriter.write("geommodel", gmds::N|gmds::F);

//    geomModel.write_oneSurfPerFile("onesurf_");
//    geomModel.write_surfaces("interfaces");


    std::cout << "geommodel.getNbPoints   " << geomModel.getNbPoints() << std::endl;
    std::cout << "geommodel.getNbCurves   " << geomModel.getNbCurves() << std::endl;
    std::cout << "geommodel.getNbSurfaces " << geomModel.getNbSurfaces() << std::endl;
    std::cout << "timer geomextract " << timer.seconds() << std::endl;


    std::cout << "timer_tot " << num_threads << " " << timer_tot.seconds() << std::endl;

    Kokkos::finalize();
    return 0;
}