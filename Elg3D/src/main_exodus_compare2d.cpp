/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    main_pixelVF2d.cpp
 *  \author  legoff
 *  \date    11/16/2018
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
#include <KM/Utils/Graph.h>
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
#include <ELG3D/ALGOCMPNT/MoveToNewPos.h>
#include <ELG3D/ALGOCMPNT/Pillow.h>
#include <ELG3D/ALGOCMPNT/PixelsRepartitioning.h>
#include <ELG3D/ALGOCMPNT/Refinement.h>
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

    elg3d::Parameters::dimension = 2;
    elg3d::Parameters::celltype = kmds::KMDS_FACE;


    kmds::Mesh mesh_source;
    elg3d::FracPres fp_source;

//    if (isFromFile) {
//        elg3d::initData_internalFormat(&mesh_source, &fp_source, case_filename);
//        elg3d::initData_unstructFormat(&mesh_source, &fp_source, case_filename);
//    } else {
//        std::cerr<<"needs to be read from a file."<<std::endl;
//    }


//    const double xyz_min[3] = {-2., -2., 0.};
//    const double xyz_max[3] = {12., 12., 0.};
//    const int ni = 12;
//    const int nj = 12;

//    const double xyz_min[3] = {-2., -2., 0.};
//    const double xyz_max[3] = {22., 12., 0.};
//    const int ni = 30;
//    const int nj = 15;

//    const double xyz_min[3] = {0., 0., 0.};
//    const double xyz_max[3] = {10., 10., 0.};
//    const int ni = 10;
//    const int nj = 10;

    const double xyz_min[3] = {0., 0., 0.};
    const double xyz_max[3] = {7., 7., 7.};
    const int ni = 12;
    const int nj = 12;

    bool fill_void = true;
    elg3d::initData_fromExodus_2D(&mesh_source, &fp_source, case_filename, xyz_min, xyz_max, ni, nj, fill_void);


    kmds::VTKWriter<kmds::Mesh> w(mesh_source);
    w.write("glpk_3x2_source", kmds::F);

    kmds::Mesh mesh_target;
    elg3d::MaterialAssignment ma_target;

    kmds::Connectivity* c_N2F = mesh_source.createConnectivity(kmds::N2F);
    kmds::ConnectivityHelper ch(&mesh_source);
    ch.buildN2F();

    kmds::Connectivity* c_E2F = mesh_source.createConnectivity(kmds::E2F);
    ch.buildEandE2F_2D_variant_1();


    kmds::Variable<kmds::TCellID> *varSubCells2oldCells =
            mesh_target.createVariable<kmds::TCellID>(kmds::NullID, kmds::KMDS_FACE, "subCells2oldCells");
    kmds::Variable<kmds::TCellID> *varOldCells2firstSubCells =
            mesh_source.createVariable<kmds::TCellID>(kmds::NullID, kmds::KMDS_FACE, "oldCells2firstSubCells");

    kmds::Variable<bool> *varMixedCells_source =
            mesh_source.createVariable<bool>(false, kmds::KMDS_FACE, "mixedCells_source");
    kmds::Variable<bool> *varMixedCells_target =
            mesh_target.createVariable<bool>(false, kmds::KMDS_FACE, "mixedCells_target");

    kmds::Variable<double> *varSurfvol_source =
            mesh_source.createVariable<double>(false, kmds::KMDS_FACE, "surfvol_source");
    kmds::Variable<double> *varSurfvol_target =
            mesh_target.createVariable<double>(false, kmds::KMDS_FACE, "surfvol_target");


    elg3d::InterfaceNodesPosSmoothVF_buildSubMesh_2D(&mesh_source,
                                                     &mesh_target,
                                                     c_N2F,
                                                     c_E2F,
                                                     &fp_source,
                                                     &ma_target,
                                                     varSubCells2oldCells,
                                                     varOldCells2firstSubCells,
                                                     varMixedCells_source,
                                                     varMixedCells_target,
                                                     varSurfvol_source,
                                                     varSurfvol_target
    );

    elg3d::Tools_write_lite_2D(&mesh_target, &ma_target, "glpk_3x2_before");

    kmds::Connectivity* c_N2F_target = mesh_target.createConnectivity(kmds::N2F);
    kmds::ConnectivityHelper ch_target(&mesh_target);
    ch_target.buildN2F();
    kmds::Connectivity* c_F2F_byN_target = mesh_target.createConnectivity(kmds::F2F_byN);
    ch_target.buildF2F_byN(c_N2F_target);

    elg3d::FracPres fp_target;


    // build graph
    kmds::Variable<kmds::TCellID> *varCells2vertices =
            mesh_target.createVariable<kmds::TCellID>(kmds::NullID, kmds::KMDS_FACE, "cell2vertices");

    kmds::GrowingView<kmds::TCellID> cellsIDs_target("CELLS", mesh_target.getNbFaces());
    mesh_target.getFaceIDs(&cellsIDs_target);

    kmds::GrowingView<kmds::TCellID> selectionCells("CELLS", mesh_target.getNbFaces());

    elg3d::InterfaceNodesPosSmoothVF_selectCells_xD(&cellsIDs_target,
                                                    c_F2F_byN_target,
                                                    varMixedCells_target,
                                                    varCells2vertices,
                                                    &selectionCells);

    kmds::Graph graph("CELLS_SUBMESH_GRAPH", selectionCells.getNbElems(), 20);
    elg3d::InterfaceNodesPosSmoothVF_buildGraph_xD(&selectionCells,
                                                   c_F2F_byN_target,
                                                   varCells2vertices,
                                                   &graph);


    kmds::GrowingView<kmds::TCellID> cellsIDs_source("CELLS", mesh_source.getNbFaces());
    mesh_source.getFaceIDs(&cellsIDs_source);


    const int nbSubPixels = elg3d::InterfaceNodesPosSmoothVF_NBSUB2;


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
                    mesh_target.createVariable<gmds::math::Point>(kmds::NullID, kmds::KMDS_FACE, "midpoints_target");

            elg3d::Tools_computeMidPointVar_2D(&cellsIDs_target,
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

            mesh_target.deleteVariable(kmds::KMDS_FACE, "midpoints_target");
        }
            break;
        case 4 : {
            kmds::Variable<gmds::math::Point> *varMidPoints_target =
                    mesh_target.createVariable<gmds::math::Point>(kmds::NullID, kmds::KMDS_FACE, "midpoints_target");

            elg3d::Tools_computeMidPointVar_2D(&cellsIDs_target,
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

            mesh_target.deleteVariable(kmds::KMDS_FACE, "midpoints_target");
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
                    mesh_target.createVariable<gmds::math::Point>(kmds::NullID, kmds::KMDS_FACE, "midpoints_target");

            elg3d::Tools_computeMidPointVar_2D(&cellsIDs_target,
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

            mesh_target.deleteVariable(kmds::KMDS_FACE, "midpoints_target");
        }
            break;

        case 7 : {
            kmds::Variable<gmds::math::Point> *varMidPoints_target =
                    mesh_target.createVariable<gmds::math::Point>(kmds::NullID, kmds::KMDS_FACE, "midpoints_target");

            elg3d::Tools_computeMidPointVar_2D(&cellsIDs_target,
                                               &mesh_target,
                                               varMidPoints_target);

            elg3d::InterfaceNodesPosSmoothVF_assignSimulatedAnnealing_horizontal_xD(&cellsIDs_source,
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

            mesh_target.deleteVariable(kmds::KMDS_FACE, "midpoints_target");
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



//    elg3d::Tools_write_lite_2D(&mesh_target, &ma_target, "glpk_5x5_after");
    elg3d::Tools_write_2D(&mesh_target, &fp_target, &ma_target, "glpk_3x2_after");


    std::cout<<"timer_tot "<<num_threads<<" "<<timer_tot.seconds()<<std::endl;

    {
        kmds::Variable<gmds::math::Point> *varMidPoints_target =
                mesh_target.createVariable<gmds::math::Point>(kmds::NullID, kmds::KMDS_FACE, "midpoints_target");

        elg3d::Tools_computeMidPointVar_2D(&cellsIDs_target,
                                           &mesh_target,
                                           varMidPoints_target);


        // gradient computation
        timer.reset();

        kmds::Variable<gmds::math::Point> *varMidPoints_source =
                mesh_source.createVariable<gmds::math::Point>(kmds::NullID, kmds::KMDS_FACE, "midpoints_source");

        elg3d::Tools_computeMidPointVar_2D(&cellsIDs_source,
                                           &mesh_source,
                                           varMidPoints_source);

        kmds::Connectivity* c_F2F_byN_source = mesh_source.createConnectivity(kmds::F2F_byN);
        ch.buildF2F_byN(c_N2F);
        std::cout<<"timer buildF2F_byN "<<timer.seconds()<<std::endl;

        timer.reset();
        kmds::Variable<elg3d::MaterialGradientComputation_midcellGradients>* varMidcellgrads_source =
                mesh_source.createVariable<elg3d::MaterialGradientComputation_midcellGradients>(elg3d::MaterialGradientComputation_MIDCELLGRADIENTS_NULL, kmds::KMDS_FACE, "midcellgrads_source");
        elg3d::MaterialGradientComputation_leastsquare_allmat_2D(&cellsIDs_source, c_F2F_byN_source, &fp_source, varMidPoints_source, varMidcellgrads_source);
        std::cout<<"timer GradientComputation "<<timer.seconds()<<std::endl;


        elg3d::Tools_write_2D(&mesh_source, &fp_source, varMidcellgrads_source, "source_gradient_2D");


//        elg3d::PixelsRepartitioning_KL_xD(&cellsIDs_source,
//                                          &cellsIDs_target,
//                                          &fp_source,
//                                          &ma_target,
//                                          varSubCells2oldCells,
//                                          varOldCells2firstSubCells,
//                                          varMixedCells_source,
//                                          varCells2vertices,
//                                          varSurfvol_source,
//                                          varSurfvol_target,
//                                          &graph,
//                                          varMidPoints_target,
//                                          nbSubPixels,
//                                          &mesh_target,
//                                          &fp_target
//        );
        elg3d::PixelsRepartitioning_KL_grad_xD(&cellsIDs_source,
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
                                               &fp_target,
                                               varMidcellgrads_source
        );

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

        mesh_target.deleteVariable(kmds::KMDS_FACE, "midpoints_target");
        mesh_source.deleteVariable(kmds::KMDS_FACE, "midcellgrads_source");
    }

    elg3d::Tools_write_2D(&mesh_target, &fp_target, &ma_target, "glpk_3x2_kl");
//    exit(0);


    // start base algo
    {
        kmds::Mesh mesh;
//        kmds::InitTools_createGrid_2D(&mesh, xyz_min, xyz_max, ni*3, nj*3);
//        const double xyz_min[3] = {0., 0., 0.};
//        const double xyz_max[3] = {3., 3., 0.};
//        const int ni = 3;
//        const int nj = 3;
//        kmds::InitTools_createGrid_2D(&mesh, xyz_min, xyz_max, ni*3, nj*3);
        kmds::InitTools_createGrid_2D(&mesh, xyz_min, xyz_max, ni, nj);
        {
            kmds::GrowingView <kmds::TCellID> cellIDs2Refine("CELLS", mesh.getNbFaces());
            cellIDs2Refine.push_back(7);
            cellIDs2Refine.push_back(8);
            cellIDs2Refine.push_back(9);
            cellIDs2Refine.push_back(19);
            cellIDs2Refine.push_back(20);
            cellIDs2Refine.push_back(21);
            cellIDs2Refine.push_back(31);
            cellIDs2Refine.push_back(32);
            cellIDs2Refine.push_back(33);
            cellIDs2Refine.push_back(43);
            cellIDs2Refine.push_back(44);
            cellIDs2Refine.push_back(45);
            cellIDs2Refine.push_back(42);
            cellIDs2Refine.push_back(54);
            cellIDs2Refine.push_back(55);
            cellIDs2Refine.push_back(56);
            cellIDs2Refine.push_back(66);
            cellIDs2Refine.push_back(67);
            cellIDs2Refine.push_back(68);
            cellIDs2Refine.push_back(65);
            cellIDs2Refine.push_back(77);
            cellIDs2Refine.push_back(78);
            cellIDs2Refine.push_back(79);
            cellIDs2Refine.push_back(89);
            cellIDs2Refine.push_back(90);
            cellIDs2Refine.push_back(91);
            cellIDs2Refine.push_back(75);
            cellIDs2Refine.push_back(76);
            cellIDs2Refine.push_back(87);
            cellIDs2Refine.push_back(88);
            cellIDs2Refine.push_back(99);
            cellIDs2Refine.push_back(100);
            cellIDs2Refine.push_back(101);
            cellIDs2Refine.push_back(84);
            cellIDs2Refine.push_back(85);
            cellIDs2Refine.push_back(86);
            cellIDs2Refine.push_back(96);
            cellIDs2Refine.push_back(97);
            cellIDs2Refine.push_back(98);
            cellIDs2Refine.push_back(108);
            cellIDs2Refine.push_back(109);
            cellIDs2Refine.push_back(110);
            cellIDs2Refine.push_back(111);

            kmds::ConnectivityHelper ch(&mesh);
            kmds::Connectivity *c_N2F = mesh.createConnectivity(kmds::N2F);
            ch.buildN2F();
            kmds::Connectivity *c_F2F_byN = mesh.createConnectivity(kmds::F2F_byN);
            ch.buildF2F_byN(c_N2F);

            kmds::Connectivity* c_E2F = mesh.createConnectivity(kmds::E2F);
            ch.buildEandE2F_2D_variant_1();

            elg3d::Refinement_refine_2D(&cellIDs2Refine, true, &mesh, c_F2F_byN, c_E2F);

            mesh.deleteConnectivity(kmds::N2F);
            mesh.deleteConnectivity(kmds::F2F_byN);
            mesh.deleteConnectivity(kmds::E2F);
        }

        elg3d::FracPres fp;

        elg3d::Tools_compute_fracpres_source_and_submesh_2D(&mesh_source,
                                                            &fp_source,
                                                            varOldCells2firstSubCells,
                                                            &mesh_target,
                                                            &ma_target,
                                                            &mesh,
                                                            &fp);

        kmds::Variable<std::uintptr_t> *v = mesh.createVariable<std::uintptr_t>(
                reinterpret_cast<std::uintptr_t>(nullptr),
                kmds::KMDS_NODE, "geomEntity");
        gmds::cad::FACManager *facGeomManager = new gmds::cad::FACManager;

        elg3d::BoundingBoxGeomAssociation_init_2D(&mesh, facGeomManager, v);

        elg3d::MaterialAssignment ma(mesh.getFaceCapacity());

        timer.reset();
        elg3d::assignCellsMajorityCriteria_2D(&mesh, &fp, &ma);
        elg3d::Tools_write_2D(&mesh, &fp, &ma, "AssignCellsTest_assignCells_I_2D_majority");

        // connectivities
        timer.reset();
        kmds::Connectivity* c_N2F = mesh.createConnectivity(kmds::N2F);
        kmds::ConnectivityHelper ch(&mesh);
        ch.buildN2F();
        std::cout<<"timer buildN2R_firsttime "<<timer.seconds()<<std::endl;


        timer.reset();
        elg3d::assignCellsCorrection_2D(&mesh, c_N2F, &fp, &ma);
        std::cout<<"timer assignCellsCorrection "<<timer.seconds()<<std::endl;

        elg3d::Tools_write_2D(&mesh, &fp, &ma, "AssignCellsTest_assignCells_I_2D_correction");

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


                elg3d::Tools_write_2D(&mesh, &fp, &ma, varMidcellgrads, "AssignCellsTest_gradient_I_2D");


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

                  elg3d::Tools_write_2D(&mesh, &fp, &ma, "AssignCellsTest_assignCells_I_2D_before_pillow");



        timer.reset();
        elg3d::pillow_2D(&facesInterfaces,
                         &nodesInterfaces,
                         &mesh,
                         c_E2F,
                         c_N2F,
                         &ma,
                         v);
        std::cout<<"timer pillow "<<timer.seconds()<<std::endl;



                elg3d::Tools_write_2D(&mesh, &fp, &ma, "AssignCellsTest_assignCells_I_2D_after_pillow");


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

                elg3d::Tools_write_2D(&mesh, &fp, &ma, "AssignCellsTest_assignCells_I_2D_after_smooth");

    }

    Kokkos::finalize();
    return 0;
}