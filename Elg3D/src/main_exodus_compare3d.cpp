/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    main_pixelVF3d.cpp
 *  \author  legoff
 *  \date    11/21/2018
 */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <KM/IO/VTKWriter.h>

#include <gmds/cad/FACManager.h>
#include <KM/DS/Mesh.h>
#include <KM/Utils/InitTools.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include <ELG3D/ALGOCMPNT/AssignCells.h>
#include <ELG3D/ALGOCMPNT/BoundingBoxGeomAssociation.h>
#include <ELG3D/ALGOCMPNT/InitData.h>
#include <ELG3D/ALGOCMPNT/InterfaceNodesPos.h>
#include <ELG3D/ALGOCMPNT/InterfaceNodesPosSmoothVF.h>
#include <ELG3D/ALGOCMPNT/MaterialInterfaces.h>
#include <ELG3D/ALGOCMPNT/MeshExtractor.h>
#include <ELG3D/ALGOCMPNT/MoveToNewPos.h>
#include <ELG3D/ALGOCMPNT/Pillow.h>
#include <ELG3D/ALGOCMPNT/PixelsRepartitioning.h>
#include <ELG3D/ALGOCMPNT/SmartLaplacian.h>
#include <ELG3D/ALGOCMPNT/Tools.h>
#include <ELG3D/DATACMPNT/FracPres.h>
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
    if (isFromFile) {
        std::istringstream iss_4(argv[4]);
        iss_4 >> case_filename;
    }

    int chosen_method;
    std::istringstream iss_5(argv[5]);
    iss_5 >> chosen_method;

    Kokkos::InitArguments kargs;
    kargs.num_threads = num_threads;
    Kokkos::initialize(kargs);

    std::cout<<"num_threads "<<num_threads<<std::endl;

    Kokkos::Timer timer_tot;
    timer_tot.reset();

    Kokkos::Timer timer;
    timer.reset();

    elg3d::Parameters::dimension = 3;
    elg3d::Parameters::celltype = kmds::KMDS_REGION;


    kmds::Mesh mesh_source;
    elg3d::FracPres fp_source;

    const double xyz_min[3] = {-2., -2., -12};
    const double xyz_max[3] = {12., 12., 2};
    const int ni = 12;
    const int nj = 12;
    const int nk = 12;

//    const double xyz_min[3] = {0., 0., 0};
//    const double xyz_max[3] = {10., 10., 10};
//    const int ni = 10;
//    const int nj = 10;
//    const int nk = 10;


    bool fill_void = true;

    elg3d::initData_fromExodus_3D(&mesh_source, &fp_source, case_filename, xyz_min, xyz_max, ni, nj, nk, fill_void);
//    elg3d::initData_internalFormat(&mesh_source, &fp_source, case_filename);

    kmds::VTKWriter<kmds::Mesh> w(mesh_source);
    w.write("glpk_3x3x3_source", kmds::R);


    kmds::Mesh mesh_target;
    elg3d::MaterialAssignment ma_target;

    kmds::Connectivity* c_N2R = mesh_source.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh_source);
    ch.buildN2R();

    kmds::Connectivity* c_E2R = mesh_source.createConnectivity(kmds::E2R);
    ch.buildEandE2R();

    kmds::Connectivity* c_F2R = mesh_source.createConnectivity(kmds::F2R);
    ch.buildFandF2R();

    kmds::Variable<kmds::TCellID> *varSubCells2oldCells =
            mesh_target.createVariable<kmds::TCellID>(kmds::NullID, kmds::KMDS_REGION, "subCells2oldCells");
    kmds::Variable<kmds::TCellID> *varOldCells2firstSubCells =
            mesh_source.createVariable<kmds::TCellID>(kmds::NullID, kmds::KMDS_REGION, "oldCells2firstSubCells");


    kmds::Variable<bool> *varMixedCells_source =
            mesh_source.createVariable<bool>(false, kmds::KMDS_REGION, "mixedCells_source");
    kmds::Variable<bool> *varMixedCells_target =
            mesh_target.createVariable<bool>(false, kmds::KMDS_REGION, "mixedCells_target");

    kmds::Variable<double> *varSurfvol_source =
            mesh_source.createVariable<double>(false, kmds::KMDS_REGION, "surfvol_source");
    kmds::Variable<double> *varSurfvol_target =
            mesh_target.createVariable<double>(false, kmds::KMDS_REGION, "surfvol_target");


//    elg3d::InterfaceNodesPosSmoothVF_buildSubMesh_3D(&mesh_source,
//                                                     &mesh_target,
//                                                     c_N2R,
//                                                     c_E2R,
//                                                     c_F2R,
//                                                     &fp_source,
//                                                     &ma_target,
//                                                     varSubCells2oldCells,
//                                                     varOldCells2firstSubCells
//    );
    elg3d::InterfaceNodesPosSmoothVF_buildSubMesh_3D(&mesh_source,
                                                     &mesh_target,
                                                     c_N2R,
                                                     c_E2R,
                                                     c_F2R,
                                                     &fp_source,
                                                     &ma_target,
                                                     varSubCells2oldCells,
                                                     varOldCells2firstSubCells,
                                                     varMixedCells_source,
                                                     varMixedCells_target,
                                                     varSurfvol_source,
                                                     varSurfvol_target
    );

    elg3d::Tools_write_lite_3D(&mesh_target, &ma_target, "glpk_3x3x3_before");
//
    kmds::Connectivity* c_N2R_target = mesh_target.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch_target(&mesh_target);
    ch_target.buildN2R();
    kmds::Connectivity* c_R2R_byN_target = mesh_target.createConnectivity(kmds::R2R_byN);
    ch_target.buildR2R_byN(c_N2R_target);

    elg3d::FracPres fp_target;


    // build graph
    kmds::Variable<kmds::TCellID> *varCells2vertices =
            mesh_target.createVariable<kmds::TCellID>(kmds::NullID, kmds::KMDS_REGION, "cell2vertices");

    kmds::GrowingView<kmds::TCellID> cellsIDs_target("CELLS", mesh_target.getNbRegions());
    mesh_target.getRegionIDs(&cellsIDs_target);

    kmds::GrowingView<kmds::TCellID> selectionCells("CELLS", mesh_target.getNbRegions());

    elg3d::InterfaceNodesPosSmoothVF_selectCells_xD(&cellsIDs_target,
                                                    c_R2R_byN_target,
                                                    varMixedCells_target,
                                                    varCells2vertices,
                                                    &selectionCells);

    kmds::Graph graph("CELLS_SUBMESH_GRAPH", selectionCells.getNbElems(), 40);
    elg3d::InterfaceNodesPosSmoothVF_buildGraph_xD(&selectionCells,
                                                   c_R2R_byN_target,
                                                   varCells2vertices,
                                                   &graph);


    kmds::GrowingView<kmds::TCellID> cellsIDs_source("CELLS", mesh_source.getNbRegions());
    mesh_source.getRegionIDs(&cellsIDs_source);


    const int nbSubPixels = elg3d::InterfaceNodesPosSmoothVF_NBSUB3;


    switch (chosen_method) {
        case 0 :

            elg3d::InterfaceNodesPosSmoothVF_assignHeuristic_xD(&cellsIDs_source,
                                                                &cellsIDs_target,
                                                                &fp_source,
                                                                &ma_target,
                                                                varSubCells2oldCells,
                                                                varOldCells2firstSubCells,
                                                                varMixedCells_source,
                                                                varCells2vertices,
                                                                varSurfvol_source,
                                                                varSurfvol_target,
                                                                &graph,
                                                                nbSubPixels
            );

            break;
        case 1 :

            elg3d::InterfaceNodesPosSmoothVF_assignMIP_xD(&cellsIDs_source,
                                                          &cellsIDs_target,
                                                          &fp_source,
                                                          &ma_target,
                                                          varSubCells2oldCells,
                                                          varOldCells2firstSubCells,
                                                          varMixedCells_source,
                                                          varCells2vertices,
                                                          &graph,
                                                          nbSubPixels
            );

            break;
        case 2 :

            elg3d::InterfaceNodesPosSmoothVF_assignLP_xD(&cellsIDs_source,
                                                         &cellsIDs_target,
                                                         &fp_source,
                                                         &fp_target,
                                                         &ma_target,
                                                         varSubCells2oldCells,
                                                         varOldCells2firstSubCells,
                                                         varMixedCells_source,
                                                         varCells2vertices,
                                                         &graph,
                                                         nbSubPixels
            );

            break;
        case 3 : {
            kmds::Variable<gmds::math::Point> *varMidPoints_target =
                    mesh_target.createVariable<gmds::math::Point>(kmds::NullID, kmds::KMDS_REGION, "midpoints_target");

            elg3d::Tools_computeMidPointVar_3D(&cellsIDs_target,
                                               &mesh_target,
                                               varMidPoints_target);

            elg3d::InterfaceNodesPosSmoothVF_assignGraphCut_xD(&cellsIDs_source,
                                                               &cellsIDs_target,
                                                               &fp_source,
                                                               &ma_target,
                                                               varSubCells2oldCells,
                                                               varOldCells2firstSubCells,
                                                               varMixedCells_source,
                                                               varCells2vertices,
                                                               &graph,
                                                               varMidPoints_target,
                                                               nbSubPixels
            );

            mesh_target.deleteVariable(kmds::KMDS_REGION, "midpoints_target");
        }
            break;
        case 4 : {
            kmds::Variable<gmds::math::Point> *varMidPoints_target =
                    mesh_target.createVariable<gmds::math::Point>(kmds::NullID, kmds::KMDS_REGION, "midpoints_target");

            elg3d::Tools_computeMidPointVar_3D(&cellsIDs_target,
                                               &mesh_target,
                                               varMidPoints_target);

            elg3d::InterfaceNodesPosSmoothVF_assignSimulatedAnnealing_xD(&cellsIDs_source,
                                                                         &cellsIDs_target,
                                                                         &fp_source,
                                                                         &ma_target,
                                                                         varSubCells2oldCells,
                                                                         varOldCells2firstSubCells,
                                                                         varMixedCells_source,
                                                                         varCells2vertices,
                                                                         &graph,
                                                                         varMidPoints_target,
                                                                         nbSubPixels
            );

            mesh_target.deleteVariable(kmds::KMDS_REGION, "midpoints_target");
        }
            break;

        case 5 :

            elg3d::InterfaceNodesPosSmoothVF_assignHeuristic_withVolumeCheck_xD(&cellsIDs_source,
                                                                                &cellsIDs_target,
                                                                                &fp_source,
                                                                                &ma_target,
                                                                                varSubCells2oldCells,
                                                                                varOldCells2firstSubCells,
                                                                                varMixedCells_source,
                                                                                varCells2vertices,
                                                                                varSurfvol_source,
                                                                                varSurfvol_target,
                                                                                &graph,
                                                                                nbSubPixels
            );
            break;

        case 6 : {
            kmds::Variable<gmds::math::Point> *varMidPoints_target =
                    mesh_target.createVariable<gmds::math::Point>(kmds::NullID, kmds::KMDS_REGION, "midpoints_target");

            elg3d::Tools_computeMidPointVar_3D(&cellsIDs_target,
                                               &mesh_target,
                                               varMidPoints_target);

            elg3d::InterfaceNodesPosSmoothVF_assignHeuristic_withVolumeCheck_horizontal_xD(&cellsIDs_source,
                                                                                           &cellsIDs_target,
                                                                                           &fp_source,
                                                                                           &ma_target,
                                                                                           varSubCells2oldCells,
                                                                                           varOldCells2firstSubCells,
                                                                                           varMixedCells_source,
                                                                                           varCells2vertices,
                                                                                           varSurfvol_source,
                                                                                           varSurfvol_target,
                                                                                           &graph,
                                                                                           varMidPoints_target,
                                                                                           nbSubPixels
            );

            mesh_target.deleteVariable(kmds::KMDS_REGION, "midpoints_target");
        }
            break;

        default:
            break;
    }

    elg3d::InterfaceNodesPosSmoothVF_summary_xD(&cellsIDs_source,
                                                &cellsIDs_target,
                                                &mesh_source,
                                                &mesh_target,
                                                &fp_source,
                                                &ma_target,
                                                varSubCells2oldCells,
                                                varOldCells2firstSubCells,
                                                varMixedCells_source,
                                                varCells2vertices,
                                                &graph,
                                                varSurfvol_source,
                                                varSurfvol_target,
                                                nbSubPixels
    );


    elg3d::Tools_write_3D(&mesh_target, &fp_target, &ma_target, "glpk_3x3x3_after");


    std::cout<<"timer_tot "<<num_threads<<" "<<timer_tot.seconds()<<std::endl;

    //==========================================================================
    // output display pixel model
    {
        kmds::Connectivity* c_F2R_target = mesh_target.createConnectivity(kmds::F2R);
        ch_target.buildFandF2R_variant_0();

        kmds::GrowingView <kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh_target.getNbNodes());
        elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh_target, c_N2R_target, &ma_target, &nodesInterfaces);

        kmds::Mesh mesh_pixels_interface;


        kmds::GrowingView<kmds::TCellID> facesInterfaces("FACES_ON_INTERFACES", mesh_target.getNbFaces());
        elg3d::MaterialInterfaces_getFaceOnInterfaces(&mesh_target, c_F2R_target, &ma_target, &facesInterfaces);

        kmds::GrowingView<kmds::TCellID> nodesIDs_old2new("NODESIDS_old2new", nodesInterfaces.getNbElems());
        kmds::GrowingView<kmds::TCellID> facesIDs_old2new("FACESIDS_old2new", facesInterfaces.getNbElems());

        elg3d::MeshExtractor_extract_NandF(&mesh_target,
                                           &mesh_pixels_interface,
                                           &nodesInterfaces,
                                           &facesInterfaces,
                                           &nodesIDs_old2new,
                                           &facesIDs_old2new
        );

        // index the faces by pair of materials
        kmds::Variable <int> *varInterfaces_int = mesh_pixels_interface.createVariable<int>(-1,kmds::KMDS_FACE, "varInterfaces_int");

        const kmds::TCellID nbFaces_interface = facesInterfaces.getNbElems();
//        Kokkos::parallel_for(nbFaces_interface,
//                             KOKKOS_LAMBDA(const int i)
        for(int i=0; i<nbFaces_interface; i++)
        {
            const kmds::TCellID fid = facesInterfaces.get(i);
            const kmds::TCellID fid_interface = facesIDs_old2new.get(i);

            Kokkos::View<kmds::TCellID*> rids;
            c_F2R_target->get(fid, rids);

            // there should always be 2 regions
            const int mat0 = ma_target.getMaterial(rids[0]);
            const int mat1 = ma_target.getMaterial(rids[1]);

            std::pair<int, int> mat_pair (std::min(mat0, mat1), std::max(mat0, mat1));

            (*varInterfaces_int)[fid_interface] = elg3d::Parameters::mat_pair_2_int.at(mat_pair);
        }
//        );


        kmds::VTKWriter<kmds::Mesh> w(mesh_pixels_interface);
        w.write("mesh_pixels_interface", kmds::F);

//        exit(-1);
    }

    {
        kmds::Variable<gmds::math::Point> *varMidPoints_target =
                mesh_target.createVariable<gmds::math::Point>(kmds::NullID, kmds::KMDS_REGION, "midpoints_target");

        elg3d::Tools_computeMidPointVar_3D(&cellsIDs_target,
                                           &mesh_target,
                                           varMidPoints_target);


        // gradient computation
        timer.reset();

        kmds::Variable<gmds::math::Point> *varMidPoints_source =
                mesh_source.createVariable<gmds::math::Point>(kmds::NullID, kmds::KMDS_REGION, "midpoints_source");

        elg3d::Tools_computeMidPointVar_3D(&cellsIDs_source,
                                           &mesh_source,
                                           varMidPoints_source);

        kmds::Connectivity* c_R2R_byN_source = mesh_source.createConnectivity(kmds::R2R_byN);
        ch.buildR2R_byN(c_N2R);
        std::cout<<"timer buildR2R_byN "<<timer.seconds()<<std::endl;

        timer.reset();
        kmds::Variable<elg3d::MaterialGradientComputation_midcellGradients>* varMidcellgrads_source =
                mesh_source.createVariable<elg3d::MaterialGradientComputation_midcellGradients>(elg3d::MaterialGradientComputation_MIDCELLGRADIENTS_NULL, kmds::KMDS_REGION, "midcellgrads_source");
        elg3d::MaterialGradientComputation_leastsquare_allmat_3D(&cellsIDs_source, c_R2R_byN_source, &fp_source, varMidPoints_source, varMidcellgrads_source);
        std::cout<<"timer GradientComputation "<<timer.seconds()<<std::endl;


        elg3d::Tools_write_3D(&mesh_source, &fp_source, varMidcellgrads_source, "AssignCellsTest_gradient_I_3D");


        elg3d::PixelsRepartitioning_KL_xD(&cellsIDs_source,
                                          &cellsIDs_target,
                                          &fp_source,
                                          &ma_target,
                                          varSubCells2oldCells,
                                          varOldCells2firstSubCells,
                                          varMixedCells_source,
                                          varCells2vertices,
                                          varSurfvol_source,
                                          varSurfvol_target,
                                          &graph,
                                          varMidPoints_target,
                                          nbSubPixels,
                                          &mesh_target,
                                          &fp_target
        );
//        elg3d::PixelsRepartitioning_FM_grad_xD(&cellsIDs_source,
//                                               &cellsIDs_target,
//                                               &fp_source,
//                                               &ma_target,
//                                               varSubCells2oldCells,
//                                               varOldCells2firstSubCells,
//                                               varMixedCells_source,
//                                               varCells2vertices,
//                                               varSurfvol_source,
//                                               varSurfvol_target,
//                                               &graph,
//                                               varMidPoints_target,
//                                               nbSubPixels,
//                                               &mesh_target,
//                                               &fp_target,
//                                               varMidcellgrads_source
//        );

        elg3d::InterfaceNodesPosSmoothVF_summary_xD(&cellsIDs_source,
                                                    &cellsIDs_target,
                                                    &mesh_source,
                                                    &mesh_target,
                                                    &fp_source,
                                                    &ma_target,
                                                    varSubCells2oldCells,
                                                    varOldCells2firstSubCells,
                                                    varMixedCells_source,
                                                    varCells2vertices,
                                                    &graph,
                                                    varSurfvol_source,
                                                    varSurfvol_target,
                                                    nbSubPixels
        );

        mesh_target.deleteVariable(kmds::KMDS_REGION, "midpoints_target");
        mesh_source.deleteVariable(kmds::KMDS_REGION, "midcellgrads_source");
    }

    elg3d::Tools_write_3D(&mesh_target, &fp_target, &ma_target, "glpk_3x3x3_kl");
//    exit(0);


    // start base algo
    {
        kmds::Mesh mesh;
        kmds::InitTools_createGrid_3D(&mesh, xyz_min, xyz_max, ni * 3, nj * 3, nk*3);
//        const double xyz_min[3] = {0., 0., 0.};
//        const double xyz_max[3] = {3., 3., 3.};
//        kmds::InitTools_createGrid_3D(&mesh, xyz_min, xyz_max, 9, 9, 9);

        elg3d::FracPres fp;

        elg3d::Tools_compute_fracpres_source_and_submesh_3D(&mesh_source,
                                                            &fp_source,
                                                            varOldCells2firstSubCells,
                                                            &mesh_target,
                                                            &ma_target,
                                                            &mesh,
                                                            &fp);

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

//        exit(-1);

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
        std::cout << "timer buildR2R_byN " << timer.seconds() << std::endl;

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


        timer.reset();
        elg3d::moveToNewPos_basicMove(&nodesInterfaces,
                                      &mesh,
                                      varNewPos,
                                      v);
        std::cout << "timer basicMove " << timer.seconds() << std::endl;


        timer.reset();
        kmds::GrowingView <kmds::TCellID> facesInterfaces("FACES_ON_INTERFACES", mesh.getNbFaces());
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
        elg3d::smartLaplacian_3D(2,
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

        std::cout<<"timer_tot "<<num_threads<<" "<<timer_tot.seconds()<<std::endl;


        //==========================================================================
        // output display geom model
        {
//            kmds::Connectivity* c_F2R = mesh.createConnectivity(kmds::F2R);
//            ch.buildFandF2R_variant_0();

            kmds::GrowingView <kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
            elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);

            kmds::Mesh mesh_pixels_interface;


            kmds::GrowingView<kmds::TCellID> facesInterfaces("FACES_ON_INTERFACES", mesh.getNbFaces());
            elg3d::MaterialInterfaces_getFaceOnInterfaces(&mesh, c_F2R, &ma, &facesInterfaces);

            kmds::GrowingView<kmds::TCellID> nodesIDs_old2new("NODESIDS_old2new", nodesInterfaces.getNbElems());
            kmds::GrowingView<kmds::TCellID> facesIDs_old2new("FACESIDS_old2new", facesInterfaces.getNbElems());

            elg3d::MeshExtractor_extract_NandF(&mesh,
                                               &mesh_pixels_interface,
                                               &nodesInterfaces,
                                               &facesInterfaces,
                                               &nodesIDs_old2new,
                                               &facesIDs_old2new
            );

            // index the faces by pair of materials
            kmds::Variable <int> *varInterfaces_int = mesh_pixels_interface.createVariable<int>(-1,kmds::KMDS_FACE, "varInterfaces_int");

            const kmds::TCellID nbFaces_interface = facesInterfaces.getNbElems();
//        Kokkos::parallel_for(nbFaces_interface,
//                             KOKKOS_LAMBDA(const int i)
            for(int i=0; i<nbFaces_interface; i++)
            {
                const kmds::TCellID fid = facesInterfaces.get(i);
                const kmds::TCellID fid_interface = facesIDs_old2new.get(i);

                Kokkos::View<kmds::TCellID*> rids;
                c_F2R->get(fid, rids);

                // there should always be 2 regions
                const int mat0 = ma.getMaterial(rids[0]);
                const int mat1 = ma.getMaterial(rids[1]);

                std::pair<int, int> mat_pair (std::min(mat0, mat1), std::max(mat0, mat1));

                (*varInterfaces_int)[fid_interface] = elg3d::Parameters::mat_pair_2_int.at(mat_pair);
            }
//        );


            kmds::VTKWriter<kmds::Mesh> w(mesh_pixels_interface);
            w.write("mesh_refined_interface", kmds::F);

//        exit(-1);
        }
    }



    Kokkos::finalize();
    return 0;
}