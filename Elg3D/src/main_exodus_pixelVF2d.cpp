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
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include <ELG3D/ALGOCMPNT/InitData.h>
#include <ELG3D/ALGOCMPNT/InterfaceNodesPosSmoothVF.h>
#include <ELG3D/ALGOCMPNT/Tools.h>
#include <ELG3D/DATACMPNT/FracPres.h>
#include <ELG3D/DATACMPNT/MaterialAssignment.h>
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


    kmds::Mesh mesh_source;
    elg3d::FracPres fp_source;

//    if (isFromFile) {
////        elg3d::initData_internalFormat(&mesh_source, &fp_source, case_filename);
//        elg3d::initData_unstructFormat(&mesh_source, &fp_source, case_filename);
//    } else {
//        std::cerr<<"needs to be read from a file."<<std::endl;
//    }

//    const double xyz_min[3] = {0., 0., 0.};
//    const double xyz_max[3] = {24., 13., 0.};
    const double xyz_min[3] = {-2., -2., 0.};
    const double xyz_max[3] = {12., 12., 0.};
    const int ni = 12;
    const int nj = 12;
//    const int ni = 34;
//    const int nj = 23;

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

    Kokkos::finalize();
    return 0;
}