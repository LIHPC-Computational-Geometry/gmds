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

#include <KM/IO/VTKWriter.h>

#include <GMDS/IO/VTKReader.h>
#include <GMDS/IO/VTKWriter.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include <ELG3D/ALGOCMPNT/AssignCells.h>
#include <ELG3D/ALGOCMPNT/BoundingBoxGeomAssociation.h>
#include <ELG3D/ALGOCMPNT/ExtractGeomModel.h>
#include <ELG3D/ALGOCMPNT/InitData.h>
#include <ELG3D/ALGOCMPNT/InterfaceNodesPos.h>
#include <ELG3D/ALGOCMPNT/ManifoldDetection.h>
#include <ELG3D/ALGOCMPNT/MaterialGradientComputation.h>
#include <ELG3D/ALGOCMPNT/MaterialInterfaces.h>
#include <ELG3D/ALGOCMPNT/MeshExtractor.h>
#include <ELG3D/ALGOCMPNT/MoveToNewPos.h>
#include <ELG3D/ALGOCMPNT/Pillow.h>
#include <ELG3D/ALGOCMPNT/SmartLaplacian.h>
#include <ELG3D/ALGOCMPNT/Tools.h>
#include <ELG3D/DATACMPNT/FracPres.h>
#include <ELG3D/DATACMPNT/FacetedSurfaceGeomServices.h>
#include <ELG3D/DATACMPNT/MaterialAssignment.h>
#include <ELG3D/DATACMPNT/Parameters.h>
/*----------------------------------------------------------------------------*/

int main(int argc, char* argv[]) {
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

    std::cout<<"num_threads "<<num_threads<<std::endl;

    Kokkos::Timer timer_tot;
    timer_tot.reset();

    Kokkos::Timer timer;
    timer.reset();


    kmds::Mesh mesh;
    elg3d::FracPres fp;

    if (isFromFile) {
        elg3d::Tools_read_fracpres_vtk_3D(&mesh, &fp, case_filename, fill_void);
//        elg3d::Tools_read_fracpres_vtk_3D(&mesh, &fp, "", fill_void);
//        elg3d::initData_TEClike_3D(&mesh, &fp, "/home/legoffn/travail/conference/imr_2019/data/simulation_volfrac/cas_gautier/Output_50_3D.dat");
    } else {
//        elg3d::initData_internalFormat(&mesh, &fp, "Elg3D/test/Samples/internalFormat_3x3x3_3D.txt");
        const double xyz_min[3] = {0., 0., 0.};
        const double xyz_max[3] = {16., 8., 8.};
        elg3d::initData_fromTXT_3D(&mesh, &fp, "volfrac.txt", xyz_min, xyz_max, 120, 60, 60, 2);
//        elg3d::initData_IxJxI_3D_nonmanifold(&mesh, &fp, 10, 10);
    }

////    elg3d::initData_internalFormat(&mesh, &fp, "Elg3D/test/Samples/internalFormat_3x3x3_3D.txt");
////    elg3d::initData_IxJxI_3D_nonmanifold(&mesh, &fp, 100, 100);
////

//    const double xyz_min[3] = {0.033583, 0.015658, 0.033182};
//    const double xyz_max[3] = {7.966427, 7.948502, 7.966026};
//    const int ni = 26;
//    const int nj = 26;
//    const int nk = 26;
//
//    elg3d::initData_fromExodus_3D(&mesh, &fp, "/home/legoffn/travail/conference/imr_2019/data/simulation_volfrac/gridfluid/cubit_files/tetmesh.e", xyz_min, xyz_max, ni, nj, nk, fill_void);
//    elg3d::initData_fromExodus_vf_3D(0, &mesh, &fp, "/home/legoffn/travail/conference/imr_2019/data/simulation_volfrac/gridfluid/Dam_Break/cubit_000039/volfrac.e.1.0");

    elg3d::Parameters::pillow_node_seed_range = 3;
    elg3d::Parameters::quality_threshold = 0.3;
    elg3d::Parameters::output_prefix = "./";

    elg3d::Parameters::badpillowing_detection_ON = false;


    kmds::Variable <std::uintptr_t> *v = mesh.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr),
                                                                             kmds::KMDS_NODE, "geomEntity");
    gmds::cad::FACManager *facGeomManager = new gmds::cad::FACManager;

    elg3d::BoundingBoxGeomAssociation_init_3D(&mesh, facGeomManager, v);


    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    timer.reset();
    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);
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
    elg3d::assignCellsCorrection_3D(&mesh, c_N2R, &fp, &ma);
    std::cout << "timer assignCellsCorrection " << timer.seconds() << std::endl;

    if (writeOutputFiles) {
        elg3d::Tools_write_3D(&mesh, &fp, &ma, "elg3d_assignCellsCorrection");
    }

    // store cells midpoints
    timer.reset();
    kmds::TSize nbCells = mesh.getNbRegions();
    kmds::GrowingView <kmds::TCellID> cellIDs("CELLS", nbCells);
    mesh.getRegionIDs_dummy(&cellIDs);
    kmds::Variable <gmds::math::Point> *varMidpoints = mesh.createVariable<gmds::math::Point>(
            gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_REGION, "midpoint");
    elg3d::Tools_computeMidPointVar_3D(&cellIDs, &mesh, varMidpoints);
    std::cout << "timer midpoints " << timer.seconds() << std::endl;

    // gradient computation
    timer.reset();
    kmds::Connectivity *c_R2R_byN = mesh.createConnectivity(kmds::R2R_byN);
    ch.buildR2R_byN(c_N2R);
    std::cout << "timer buildR2R_byN " << timer.seconds() << std::endl ;

    timer.reset();
    kmds::Variable <elg3d::MaterialGradientComputation_midcellGradients> *varMidcellgrads =
            mesh.createVariable<elg3d::MaterialGradientComputation_midcellGradients>(
                    elg3d::MaterialGradientComputation_MIDCELLGRADIENTS_NULL, kmds::KMDS_REGION, "midcellgrads");
    elg3d::MaterialGradientComputation_leastsquare_3D(&cellIDs, c_R2R_byN, &fp, &ma, varMidpoints, varMidcellgrads);
    std::cout << "timer GradientComputation " << timer.seconds() << std::endl;


    // compute interface nodes new position
    timer.reset();
    kmds::GrowingView <kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);
    std::cout << "timer getNodeOnInterfaces " << timer.seconds() << std::endl;

    timer.reset();
    ch.buildFandF2R_variant_0();
    std::cout << "timer buildFandF2R " << timer.seconds() << std::endl;
    timer.reset();
    kmds::Connectivity *c_F2R = mesh.getConnectivity(kmds::F2R);
    kmds::Connectivity *c_N2F = mesh.createConnectivity(kmds::N2F);
    ch.buildN2F();
    std::cout << "timer buildN2F " << timer.seconds() << std::endl;

    kmds::Variable <gmds::math::Point> *varNewPos =
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


//    timer.reset();

    // we store the current node positions
    kmds::Variable <gmds::math::Point> *varRefPos =
            mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF),
                                                   kmds::KMDS_NODE, "refPos");

    elg3d::Tools_storeNodePos_xD(&mesh, varRefPos);



//    elg3d::moveToNewPos_basicMove(&nodesInterfaces,
//                                  &mesh,
//                                  varNewPos,
//                                  v);

    ////////////////////////////////////////

//// INTERFACE projection and smooth
//    gmds::cad::FACManager geomModel_interface;
//    gmds::Mesh &gmdsmesh_interface = geomModel_interface.getMeshView();
//    gmds::VTKReader<gmds::Mesh> vtkreader(gmdsmesh_interface);
//    vtkreader.read("geommodel_pixels_FACES.vtk");
//    exit(-1);
    {
        kmds::Mesh mesh_interface;

        kmds::GrowingView<kmds::TCellID> facesInterfaces("FACES_ON_INTERFACES", mesh.getNbFaces());
        elg3d::MaterialInterfaces_getFaceOnInterfaces(&mesh, c_F2R, &ma, &facesInterfaces);

        kmds::GrowingView<kmds::TCellID> nodesIDs_old2new("NODESIDS_old2new", nodesInterfaces.getNbElems());
        kmds::GrowingView<kmds::TCellID> facesIDs_old2new("FACESIDS_old2new", facesInterfaces.getNbElems());

        elg3d::MeshExtractor_extract_NandF(&mesh,
                                           &mesh_interface,
                                           &nodesInterfaces,
                                           &facesInterfaces,
                                           &nodesIDs_old2new,
                                           &facesIDs_old2new
        );


        kmds::VTKWriter<kmds::Mesh> w(mesh_interface);
        w.write("mesh_interface", kmds::F);

        kmds::ConnectivityHelper ch_interface(&mesh_interface);
        kmds::Connectivity *c_N2F_interface = mesh_interface.createConnectivity(kmds::N2F);
        ch_interface.buildN2F();

        kmds::Connectivity *c_N2N_interface = mesh_interface.createConnectivity(kmds::N2N);
        ch_interface.buildN2N(kmds::F);
//    kmds::Connectivity *c_N2N_interface = mesh_interface.getConnectivity(kmds::N2N);

        gmds::cad::FACManager geomModel_interface;
        gmds::Mesh &gmdsmesh_interface = geomModel_interface.getMeshView();
        gmds::VTKReader<gmds::Mesh> vtkreader(gmdsmesh_interface);
//        vtkreader.read("../../temp/GMDS/imr_2019/pixels/asteroid/geommodel_FACES.vtk");
        vtkreader.read("geommodel_pixels_FACES.vtk");

        std::cout << "geomModel_interface nbNodes " << gmdsmesh_interface.getNbNodes() << " nbFaces "
                  << gmdsmesh_interface.getNbFaces() << std::endl;

        gmds::cad::FACSurface *surf = dynamic_cast<gmds::cad::FACSurface *> (geomModel_interface.newSurface());

        std::vector<gmds::Face> faces_to_add;
        for (int i = 0; i < gmdsmesh_interface.getNbFaces(); i++) {
            faces_to_add.push_back(gmdsmesh_interface.get<gmds::Face>(i));
        }

        surf->setMeshFaces(faces_to_add);

        std::vector<gmds::cad::FACSurface *> surfaces;
        surfaces.push_back(surf);
        elg3d::FacetedSurfaceGeomServices geomservices_surf;
        geomservices_surf.buildAABBSurfacesTriangulationTrees(surfaces);

        elg3d::FacetedCurveGeomServices geomservices_curv;

//        elg3d::smartLaplacian_interface_fromF_3D(30,
//                                                 &mesh_interface,
//                                                 c_N2F_interface,
//                                                 c_N2N_interface,
//                                                 surf,
//                                                 &geomservices);


        kmds::Variable <std::uintptr_t> *v_interface_surf = mesh_interface.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(surf),
                                                                                                     kmds::KMDS_NODE, "geomEntity_surf");

        kmds::Variable <std::uintptr_t> *v_interface_bb = mesh_interface.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr),
                                                                                                        kmds::KMDS_NODE, "geomEntity_bb");
        const double xyz_min[3] = {0., 0., 0.};
        const double xyz_max[3] = {16., 8., 8.};
        elg3d::BoundingBoxGeomAssociation_assoc_3D(&mesh_interface, xyz_min, xyz_max, facGeomManager, v_interface_bb);

        elg3d::smartLaplacian_interface_fromF_doubleGeom_3D(30,
                                                            &mesh_interface,
                                                            c_N2F_interface,
                                                            c_N2N_interface,
                                                            v_interface_surf,
                                                            v_interface_bb,
                                                            &geomservices_surf,
                                                            &geomservices_curv);

        w.write("mesh_interface_projected", kmds::F);

        for(int i=0; i<nodesInterfaces.getNbElems(); i++) {
            (*varNewPos)[nodesInterfaces.get(i)] = mesh_interface.getNodeLocation(nodesIDs_old2new.get(i));
        }

    }


//    std::cout << "timer basicMove " << timer.seconds() << std::endl;

    // we extract the geom model
//    gmds::cad::FACManager AGeomModel;
//    elg3d::extractGeomModel_extract_3D(&mesh, &ma, c_F2R, &AGeomModel);
//    gmds::Mesh &meshGeom = AGeomModel.getMeshView();
//    gmds::VTKWriter<gmds::Mesh> vtkwriter(meshGeom);
//    vtkwriter.write("geommodel", gmds::N|gmds::F);
//    elg3d::FacetedSurfaceGeomServices geomServices;



    // we store the current node positions, which is the target node position projected on the
    // the geom association
    kmds::Variable <gmds::math::Point> *varnewPosWithProj =
            mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF),
                                                   kmds::KMDS_NODE, "newPosWithProj");
    elg3d::Tools_storeNodePos_xD(&mesh, varnewPosWithProj);

    // we rollback the node displacement
    elg3d::Tools_loadNodePos_xD(&mesh, varRefPos);


    ////////////////////////////////////////
//// testing the noBadMove
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

    elg3d::Parameters::quality_threshold = 0.2;

    // testing the noBadMove
    kmds::Variable <double> *varCellQuality =
            mesh.createVariable<double>(-HUGE_VALF, kmds::KMDS_REGION, "varCellQuality");
    elg3d::Tools_computeScaledJacobian_3D(&mesh,varCellQuality);

    kmds::Variable <double> *varDist =
            mesh.createVariable<double>(-HUGE_VALF, kmds::KMDS_NODE, "varDist");
    double sumdist_nobadmove = elg3d::Tools_computeDistance_xD(&nodesInterfaces, &mesh, varNewPos, varDist);
    std::cout<<"sumdist_nobadmove "<<sumdist_nobadmove<<std::endl;


    if (writeOutputFiles) {
            elg3d::Tools_write_lite_3D(&mesh, &ma, "AssignCellsTest_assignCells_I_3D_afternoBadMove");
    }




//    // we extract the geom model
//    gmds::cad::FACManager AGeomModel_nobadmove;
//    elg3d::extractGeomModel_extract_3D(&mesh, &ma, c_F2R, &AGeomModel_nobadmove);
//    std::vector<gmds::cad::GeomSurface*> surfaces;
//    AGeomModel_nobadmove.getSurfaces(surfaces);
//    std::cout<<"AGeomModel_nobadmove nbSurf "<<surfaces.size()<<std::endl;
//    for(int isurf=0; isurf<surfaces.size(); isurf++) {
//        std::cout<<surfaces[isurf]->getName()<<std::endl;
//    }
//    gmds::cad::GeomSurface* surf = surfaces[1];
//
//    kmds::Variable <std::uintptr_t> *varGeomAssoc_smoothing = mesh.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr),
//                                                                                                  kmds::KMDS_NODE, "geomAssoc_smoothing");
//    Kokkos::parallel_for(nodesInterfaces.getNbElems(),
//                         KOKKOS_LAMBDA(const int i) {
//
//                             kmds::TCellID nid = nodesInterfaces.get(i);
//                             (*varGeomAssoc_smoothing)[nid] = reinterpret_cast<std::uintptr_t> (surf);
//
//                         });


    // smooth the mesh
    // allow all the nodes to move
//    ch.buildN2N(kmds::R);
//    std::cout << "timer buildN2N " << timer.seconds() << std::endl;
//    kmds::Connectivity *c_N2N = mesh.getConnectivity(kmds::N2N);

//    kmds::GrowingView <kmds::TCellID> nodesFixed("FIXED_NODES", mesh.getNbNodes());
//
//    elg3d::smartLaplacian_interface_3D(10,
//                                       &mesh,
//                                       c_N2R,
//                                       c_N2N,
//                                       v,
//                                       varGeomAssoc_smoothing,
//                                       &nodesInterfaces);
//    elg3d::Tools_computeScaledJacobian_3D(&mesh,varCellQuality);
//    elg3d::Tools_computeDistance_xD(&nodesInterfaces, &mesh, varnewPosWithProj, varDist);
//
//
//
//    if (writeOutputFiles) {
//        elg3d::Tools_write_lite_3D(&mesh, &ma, "AssignCellsTest_assignCells_I_3D_interface_smoothing");
//    }
//
//    mesh.deleteConnectivity(kmds::N2N);
//    exit(-1);

//    mesh.createConnectivity(kmds::N2N);
//    ch.buildN2N(kmds::R);
//    std::cout << "timer buildN2N " << timer.seconds() << std::endl;
//    kmds::Connectivity *c_N2N = mesh.getConnectivity(kmds::N2N);
//    c_N2N = mesh.getConnectivity(kmds::N2N);
//
//
//    timer.reset();
//    elg3d::smartLaplacian_3D(3,
//                             &mesh,
//                             c_N2R,
//                             c_N2N,
//                             v,
//                             &nodesInterfaces);
//    mesh.deleteConnectivity(kmds::N2N);
//
//    if (writeOutputFiles) {
//        elg3d::Tools_write_lite_3D(&mesh, &ma, "AssignCellsTest_assignCells_I_3D_beforepillow_smoothing");
//    }
////////////////////////////////////////



//    timer.reset();
//    elg3d::moveToNewPos_basicMove(&nodesInterfaces,
//                                  &mesh,
//                                  varNewPos,
//                                  v);
//    std::cout<<"timer basicMove "<<timer.seconds()<<std::endl;
//
//    elg3d::Tools_computeScaledJacobian_3D(&mesh,varCellQuality);
//
//
//    timer.reset();
//    kmds::GrowingView <kmds::TCellID> facesInterfaces("FACES_ON_INTERFACES", mesh.getNbFaces());
//    elg3d::MaterialInterfaces_getFaceOnInterfaces(&mesh, c_F2R, &ma, &facesInterfaces);
//    std::cout << "timer getFaceOnInterfaces " << timer.seconds() << std::endl;
//
//    if (writeOutputFiles) {
//        elg3d::Tools_write_interfaces_3D(&facesInterfaces,
//                                         &nodesInterfaces,
//                                         &mesh,
//                                         c_F2R,
//                                         &ma,
//                                         "PillowingTest_internalFormat_interface_before_move_3x3x3_3D");
//
//
//        elg3d::Tools_write_lite_3D(&mesh, &ma, "PillowingTest_internalFormat_before_pillow_3x3x3_3D");
//    }

    timer.reset();
//    elg3d::Tools_loadNodePos_xD(&mesh, varRefPos);
//    elg3d::pillow_3D(&facesInterfaces,
//                     &nodesInterfaces,
//                     &mesh,
//                     c_F2R,
//                     c_N2R,
//                     &ma,
//                     v);
//    elg3d::Parameters::quality_threshold = 0.2;
    elg3d::pillow_execute_propagate_3D(&nodesInterfaces,
                                       &mesh,
                                       c_F2R,
                                       c_N2F,
                                       c_N2R,
                                       &ma,
                                       v,
                                       varDist);

    std::cout << "timer pillow " << timer.seconds() << std::endl;

    elg3d::Tools_computeScaledJacobian_3D(&mesh,varCellQuality);
    elg3d::Tools_computeDistance_xD(&nodesInterfaces, &mesh, varnewPosWithProj, varDist);

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
//  mesh.createConnectivity(kmds::N2N);
    ch.buildN2N(kmds::R);
    std::cout << "timer buildN2N " << timer.seconds() << std::endl;
    kmds::Connectivity *c_N2N = mesh.getConnectivity(kmds::N2N);
//    c_N2N = mesh.getConnectivity(kmds::N2N);

    timer.reset();
    elg3d::smartLaplacian_3D(10,
                             &mesh,
                             c_N2R,
                             c_N2N,
                             v,
                             &nodesInterfaces);
    std::cout << "timer smartLaplacian " << timer.seconds() << std::endl;

    elg3d::Tools_computeScaledJacobian_3D(&mesh,varCellQuality);
    elg3d::Tools_computeDistance_xD(&nodesInterfaces, &mesh, varnewPosWithProj, varDist);

    std::cout << "smartLaplacian_3D after exit" << std::endl;
    if (writeOutputFiles) {
        elg3d::Tools_write_lite_3D(&mesh, &ma, "PillowingTest_internalFormat_after_smooth_3x3x3_3D");
    }

    elg3d::moveToNewPos_noBadMove_3D(16,
                                     elg3d::Parameters::quality_threshold,
                                     &nodesInterfaces,
                                     &mesh,
                                     c_N2R,
                                     varNewPos,
                                     v,
                                     varNbMoves);
    elg3d::Tools_computeScaledJacobian_3D(&mesh,varCellQuality);
    double sumdist_final = elg3d::Tools_computeDistance_xD(&nodesInterfaces, &mesh, varNewPos, varDist);
    std::cout<<"sumdist_final "<<sumdist_final<<std::endl;
    if (writeOutputFiles) {
        elg3d::Tools_write_lite_3D(&mesh, &ma, "PillowingTest_internalFormat_poyop_3x3x3_3D");
    }

    ////////////////////////////////////////////
// GETME
//    elg3d::Tools_callGETME_3D(&mesh, &nodesInterfaces, v);
//
//    elg3d::Tools_computeScaledJacobian_3D(&mesh,varCellQuality);
//
//    if (writeOutputFiles) {
//        elg3d::Tools_write_lite_3D(&mesh, &ma, "PillowingTest_internalFormat_after_GETMe_3x3x3_3D");
//    }
////////////////////////////////////////////

    std::cout<<"timer_tot "<<num_threads<<" "<<timer_tot.seconds()<<std::endl;

    Kokkos::finalize();
    return 0;
}