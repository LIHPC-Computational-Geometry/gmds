/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * This software is a computer program whose purpose is to provide a set of
 * functionnalities to represent and handle any type of meshes (2D, 3D,
 * triangles, tetrahedra, quad, hexa, polygons, polyhedra, etc.) and write
 * meshing algorithms. So it gathers many mathematical objects like points,
 * segment, quaternions, etc. and basic algorithms useful to build more evolved
 * ones.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and, more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
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

#include <gmds/ig/Mesh.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/IGMeshIOService.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include <ELG3D/ALGOCMPNT/ExtractGeomModel.h>
#include <ELG3D/ALGOCMPNT/InitData.h>
#include <ELG3D/ALGOCMPNT/InterfaceNodesPosSmoothVF.h>
#include <ELG3D/ALGOCMPNT/PixelsRepartitioning.h>
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

    if (isFromFile) {
        elg3d::initData_internalFormat(&mesh_source, &fp_source, case_filename);
//        elg3d::initData_unstructFormat(&mesh_source, &fp_source, case_filename);
//        elg3d::Tools_read_fracpres_vtk_2D(&mesh_source, &fp_source, case_filename, 1);
//        elg3d::initData_TEClike_2D(&mesh_source, &fp_source, "/home/legoffn/travail/conference/imr_2019/data/simulation_volfrac/cas_gautier/Output_50.dat");
    } else {
        std::cerr<<"needs to be read from a file."<<std::endl;
    }


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

//            elg3d::InterfaceNodesPosSmoothVF_assignHeuristic_withVolumeCheck_xD(&cellsIDs_source,
//                                                                                &cellsIDs_target,
//                                                                                &fp_source,
//                                                                                &ma_target,
//                                                                                varSubCells2oldCells,
//                                                                                varOldCells2firstSubCells,
//                                                                                varMixedCells_source,
//                                                                                varCells2vertices,
//                                                                                varSurfvol_source,
//                                                                                varSurfvol_target,
//                                                                                &graph,
//                                                                                nbSubPixels
//            );
            elg3d::InterfaceNodesPosSmoothVF_assignHeuristic_withVolumeCheck_print_xD(
                    &mesh_target,&cellsIDs_source,
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

        elg3d::PixelsRepartitioning_FM_xD(&cellsIDs_source,
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
    }

    elg3d::Tools_write_2D(&mesh_target, &fp_target, &ma_target, "glpk_3x2_kl");
    exit(0);


    // geometric model extraction
    kmds::Connectivity* c_E2F_target = mesh_target.createConnectivity(kmds::E2F);
    ch_target.buildEandE2F_2D_variant_0();

    gmds::cad::FACManager AGeomModel;
    elg3d::extractGeomModel_extract_2D(&mesh_target, &ma_target, c_E2F_target, &AGeomModel);
    gmds::Mesh &meshGeom = AGeomModel.getMeshView();
    gmds::IGMeshIOService ioService(&meshGeom);
    gmds::VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::F);
    vtkWriter.setDataOptions(gmds::N|gmds::F);
    vtkWriter.write("geommodel");




    Kokkos::finalize();
    return 0;
}