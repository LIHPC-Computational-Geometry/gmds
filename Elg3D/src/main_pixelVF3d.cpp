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

#include <gmds/cadfac/FACManager.h>
#include <KM/DS/Mesh.h>

#include <gmds/ig/Mesh.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/IGMeshIOService.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include <ELG3D/ALGOCMPNT/ExtractGeomModel.h>
#include <ELG3D/ALGOCMPNT/InitData.h>
#include <ELG3D/ALGOCMPNT/InterfaceNodesPosSmoothVF.h>
#include <ELG3D/ALGOCMPNT/Tools.h>
#include <ELG3D/DATACMPNT/FracPres.h>
#include <ELG3D/DATACMPNT/MaterialAssignment.h>


#include <ELG3D/ALGOCMPNT/AssignCells.h>
#include <ELG3D/ALGOCMPNT/BoundingBoxGeomAssociation.h>
#include <ELG3D/ALGOCMPNT/InterfaceNodesPos.h>
#include <ELG3D/ALGOCMPNT/ManifoldDetection.h>
#include <ELG3D/ALGOCMPNT/MaterialGradientComputation.h>
#include <ELG3D/ALGOCMPNT/MaterialInterfaces.h>
#include <ELG3D/ALGOCMPNT/MeshExtractor.h>
#include <ELG3D/ALGOCMPNT/MoveToNewPos.h>
#include <ELG3D/ALGOCMPNT/Pillow.h>
#include <ELG3D/ALGOCMPNT/SmartLaplacian.h>
#include <ELG3D/DATACMPNT/FacetedCurveGeomServices.h>
#include <ELG3D/DATACMPNT/FacetedSurfaceGeomServices.h>
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

    if (isFromFile) {
        elg3d::initData_internalFormat(&mesh_source, &fp_source, case_filename);
//        elg3d::Tools_read_fracpres_vtk_3D(&mesh_source, &fp_source, case_filename, 1);

    } else {
        const double xyz_min[3] = {0., 0., 0.};
        const double xyz_max[3] = {16., 8., 8.};
        elg3d::initData_fromTXT_3D(&mesh_source, &fp_source, "volfrac.txt", xyz_min, xyz_max, 120, 60, 60, 2);
//        std::cerr<<"needs to be read from a file."<<std::endl;
    }

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

//    exit(-1);

    // geometric model extraction
    kmds::Connectivity* c_F2R_target = mesh_target.createConnectivity(kmds::F2R);
    ch_target.buildFandF2R_variant_0();

    kmds::Variable <std::uintptr_t> *geomassoc_pixels = mesh_target.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr),
                                                                                             kmds::KMDS_NODE, "geomassoc_pixels");

    gmds::cad::FACManager geomModel_pixels;
    elg3d::extractGeomModel_extract_3D(&mesh_target, &ma_target, c_F2R_target, &geomModel_pixels, geomassoc_pixels);
    gmds::Mesh &meshGeom_pixels = geomModel_pixels.getMeshView();

    gmds::IGMeshIOService ioService(&meshGeom_pixels);
    gmds::VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::F);
    vtkWriter.setDataOptions(gmds::N|gmds::F);
    vtkWriter.write("geommodel_pixels");

//    geomModel_pixels.write_oneSurfPerFile("poyop_");


    std::cout<<"geommodel_pixels.getNbPoints   "<<geomModel_pixels.getNbPoints()<<std::endl;
    std::cout<<"geommodel_pixels.getNbCurves   "<<geomModel_pixels.getNbCurves()<<std::endl;
    std::cout<<"geommodel_pixels.getNbSurfaces "<<geomModel_pixels.getNbSurfaces()<<std::endl;
    exit(-1);

    //==========================================================================
    // output display pixel model
    {
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

        exit(-1);
    }


    //==========================================================================



    kmds::Mesh mesh_bis;
    elg3d::FracPres fp_bis;

    if (isFromFile) {
//        elg3d::initData_internalFormat(&mesh_bis, &fp_bis, case_filename);
        elg3d::Tools_read_fracpres_vtk_3D(&mesh_bis, &fp_bis, case_filename, 1);
//        const double xyz_min[3] = {0., 0., 0.};
//        const double xyz_max[3] = {8., 8., 8.};
//        elg3d::initData_fromTXT_3D(&mesh_source, &fp_source, "volfrac.txt", xyz_min, xyz_max, 60, 60, 60, 2);
    } else {
        std::cerr<<"needs to be read from a file."<<std::endl;
    }

    elg3d::Parameters::pillow_node_seed_range = 3;
    elg3d::Parameters::quality_threshold = 0.3;
    elg3d::Parameters::output_prefix = "./";


    kmds::Variable <std::uintptr_t> *v_boundingbox = mesh_bis.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr),
                                                                             kmds::KMDS_NODE, "geomEntity");
    gmds::cad::FACManager *facGeomManager_bis = new gmds::cad::FACManager;

    elg3d::BoundingBoxGeomAssociation_init_3D(&mesh_bis, facGeomManager_bis, v_boundingbox);


    elg3d::MaterialAssignment ma_bis(mesh_bis.getRegionCapacity());

    timer.reset();
    elg3d::assignCellsMajorityCriteria_3D(&mesh_bis, &fp_bis, &ma_bis);

    if (writeOutputFiles) {
        elg3d::Tools_write_3D(&mesh_bis, &fp_bis, &ma_bis, "elg3d_assignCellsMajorityCriteria");
    }

    // geometric model extraction
    kmds::ConnectivityHelper ch_bis(&mesh_bis);

    kmds::Connectivity* c_F2R_bis = mesh_bis.createConnectivity(kmds::F2R);
    ch_bis.buildFandF2R_variant_0();

    kmds::Variable <std::uintptr_t> *geomassoc_bis = mesh_bis.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr),
                                                                                                   kmds::KMDS_NODE, "geomassoc_bis");

    gmds::cad::FACManager geomModel_bis;
    elg3d::extractGeomModel_extract_3D(&mesh_bis, &ma_bis, c_F2R_bis, &geomModel_bis, geomassoc_bis);
    gmds::Mesh &meshGeom_bis = geomModel_bis.getMeshView();
//    gmds::VTKWriter<gmds::Mesh> vtkwriter_bis(meshGeom_bis);
//    vtkwriter_bis.write("geommodel_bis", gmds::N|gmds::F);

    gmds::IGMeshIOService ioService_bis(&meshGeom_bis);
    gmds::VTKWriter vtkWriter_bis(&ioService_bis);
    vtkWriter_bis.setCellOptions(gmds::N|gmds::F);
    vtkWriter_bis.setDataOptions(gmds::N|gmds::F);
    vtkWriter_bis.write("geommodel_bis");

//    {
//        std::uintptr_t geom_ptr = (*geomassoc_bis)[22];
//        <gmds::cad::GeomEntity *geomEntity = reinterpret_cast<gmds::cad::GeomEntity *> (geom_ptr);
//        std::cout<<"node 22 dimension "<<geomEntity->getDim()<<std::endl;
//        std::cout<<"geomEntity_name "<<geomEntity->getName()<<std::endl;
//        gmds::cad::FACSurface *surf = reinterpret_cast<gmds::cad::FACSurface *> (geom_ptr);
//        std::cout<<"surf_name "<<surf->getName()<<std::endl;
//    }
//    {
//        std::uintptr_t curve_ptr = (*geomassoc_bis)[2];
//        <gmds::cad::GeomEntity *geomEntity = reinterpret_cast<gmds::cad::GeomEntity *> (curve_ptr);
//        std::cout<<"node 2 dimension "<<geomEntity->getDim()<<std::endl;
//        std::cout<<"geomEntity_name "<<geomEntity->getName()<<std::endl;
//        gmds::cad::FACCurve* curv = reinterpret_cast<gmds::cad::FACCurve*> (curve_ptr);
//        std::cout<<"curve_name "<<curv->getName()<<std::endl;
//    }


    // working on the geom models
    std::vector<gmds::cad::GeomSurface*> surfaces_pixels;
    std::vector<gmds::cad::GeomCurve*> curves_pixels;
    std::vector<gmds::cad::GeomPoint*> vertices_pixels;
    geomModel_pixels.getSurfaces(surfaces_pixels);
    geomModel_pixels.getCurves(curves_pixels);
    geomModel_pixels.getPoints(vertices_pixels);

    std::vector<gmds::cad::GeomSurface*> surfaces_bis;
    std::vector<gmds::cad::GeomCurve*> curves_bis;
    std::vector<gmds::cad::GeomPoint*> vertices_bis;
    geomModel_bis.getSurfaces(surfaces_bis);
    geomModel_bis.getCurves(curves_bis);
    geomModel_bis.getPoints(vertices_bis);


    // preparing the geom services for fast projection
    std::vector<gmds::cad::FACSurface*> fsurfaces_pixels;
    for(auto surf: surfaces_pixels) {
        fsurfaces_pixels.push_back(dynamic_cast<gmds::cad::FACSurface*> (surf));
    }
    std::vector<gmds::cad::FACCurve*> fcurves_pixels;
    for(auto curv: curves_pixels) {
        fcurves_pixels.push_back(dynamic_cast<gmds::cad::FACCurve*> (curv));
    }

    elg3d::FacetedSurfaceGeomServices surfaceGeomServices_pixels;
    surfaceGeomServices_pixels.buildAABBSurfacesTriangulationTrees(fsurfaces_pixels);
    elg3d::FacetedCurveGeomServices curveGeomServices_pixels;
    curveGeomServices_pixels.buildAABBCurvesTriangulationTrees(fcurves_pixels);


    // test projection
//    {
//        gmds::math::Point pt_0(1.5, 1.5, 1.5);
//        surfaceGeomServices_pixels.project(surfaces_pixels[0], pt_0);
//        std::cout<<"proj on surf "<<surfaces_pixels[0]->getName()<<" "<<pt_0<<std::endl;
//
//        gmds::math::Point pt_1(1.5, 1.5, 1.5);
//        surfaceGeomServices_pixels.project(surfaces_pixels[1], pt_1);
//        std::cout<<"proj on surf "<<surfaces_pixels[1]->getName()<<" "<<pt_1<<std::endl;
//
//        gmds::math::Point pt_2(1.5, 1.5, 1.5);
//        surfaceGeomServices_pixels.project(surfaces_pixels[2], pt_2);
//        std::cout<<"proj on surf "<<surfaces_pixels[2]->getName()<<" "<<pt_2<<std::endl;
//
//        gmds::math::Point pt_4(1, 1, 1.5);
//        curveGeomServices_pixels.project(curves_pixels[0], pt_4);
//        std::cout<<"proj on curv "<<curves_pixels[0]->getName()<<" "<<pt_4<<std::endl;
//
//        gmds::math::Point pt_5(1, 1, 1.5);
//        curveGeomServices_pixels.project(curves_pixels[1], pt_5);
//        std::cout<<"proj on curv "<<curves_pixels[1]->getName()<<" "<<pt_5<<std::endl;
//
//        gmds::math::Point pt_6(1, 1, 1.5);
//        curveGeomServices_pixels.project(curves_pixels[2], pt_6);
//        std::cout<<"proj on curv "<<curves_pixels[2]->getName()<<" "<<pt_6<<std::endl;
//
//        gmds::math::Point pt_7(1, 1, 1.5);
//        curveGeomServices_pixels.project(curves_pixels[3], pt_7);
//        std::cout<<"proj on curv "<<curves_pixels[3]->getName()<<" "<<pt_7<<std::endl;
//
//        gmds::math::Point pt_8(1, 1, 1.5);
//        vertices_pixels[0]->project(pt_8);
//        std::cout<<"proj on vert "<<vertices_pixels[0]->getName()<<" "<<pt_8<<std::endl;
//    }


    // mapping between the grid and the pixels geom models
    std::map<std::uintptr_t, std::uintptr_t> geom_bis2pixels;
    {

        for(auto surf_bis: surfaces_bis) {

            bool found = false;

            for(auto surf_pixels: surfaces_pixels) {

                if(surf_bis->name() == surf_pixels->name()) {
                    geom_bis2pixels[reinterpret_cast<std::uintptr_t> (surf_bis)] = reinterpret_cast<std::uintptr_t> (surf_pixels);
                    found = true;
                }
            }

            if(!found) {
                std::cout<<"SURF MAPPING NOT FOUND"<<std::endl;
                geom_bis2pixels[reinterpret_cast<std::uintptr_t> (surf_bis)] = reinterpret_cast<std::uintptr_t> (nullptr);
            }
        }

        for(auto curv_bis: curves_bis) {

            bool found = false;

            for(auto curv_pixels: curves_pixels) {

                if(curv_bis->name() == curv_pixels->name()) {
                    geom_bis2pixels[reinterpret_cast<std::uintptr_t> (curv_bis)] = reinterpret_cast<std::uintptr_t> (curv_pixels);
                    found = true;
                }
            }

            if(!found) {
                std::cout<<"CURVE MAPPING NOT FOUND "<<curv_bis->name()<<std::endl;
                geom_bis2pixels[reinterpret_cast<std::uintptr_t> (curv_bis)] = reinterpret_cast<std::uintptr_t> (nullptr);
            }
        }

        for(auto vert_bis: vertices_bis) {

            bool found = false;

            for(auto vert_pixels: vertices_pixels) {

                if(vert_bis->name() == vert_pixels->name()) {
                    geom_bis2pixels[reinterpret_cast<std::uintptr_t> (vert_bis)] = reinterpret_cast<std::uintptr_t> (vert_pixels);
                    found = true;
                }
            }

            if(!found) {
                geom_bis2pixels[reinterpret_cast<std::uintptr_t> (vert_bis)] = reinterpret_cast<std::uintptr_t> (nullptr);
            }
        }

    }


    // interface projection and smooth
    {
        kmds::Connectivity *c_N2R_bis = mesh_bis.createConnectivity(kmds::N2R);
        ch_bis.buildN2R();

        kmds::GrowingView <kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh_bis.getNbNodes());
        elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh_bis, c_N2R_bis, &ma_bis, &nodesInterfaces);



        kmds::Mesh mesh_interface;


        kmds::GrowingView<kmds::TCellID> facesInterfaces("FACES_ON_INTERFACES", mesh_bis.getNbFaces());
        elg3d::MaterialInterfaces_getFaceOnInterfaces(&mesh_bis, c_F2R_bis, &ma_bis, &facesInterfaces);

        kmds::GrowingView<kmds::TCellID> nodesIDs_old2new("NODESIDS_old2new", nodesInterfaces.getNbElems());
        kmds::GrowingView<kmds::TCellID> facesIDs_old2new("FACESIDS_old2new", facesInterfaces.getNbElems());

        elg3d::MeshExtractor_extract_NandF(&mesh_bis,
                                           &mesh_interface,
                                           &nodesInterfaces,
                                           &facesInterfaces,
                                           &nodesIDs_old2new,
                                           &facesIDs_old2new
        );

        // index the faces by pair of materials
        kmds::Variable <int> *varInterfaces_int = mesh_interface.createVariable<int>(-1,kmds::KMDS_FACE, "varInterfaces_int");

        const kmds::TCellID nbFaces_interface = facesInterfaces.getNbElems();
//        Kokkos::parallel_for(nbFaces_interface,
//                             KOKKOS_LAMBDA(const int i)
for(int i=0; i<nbFaces_interface; i++)
                             {
                                 const kmds::TCellID fid = facesInterfaces.get(i);
                                 const kmds::TCellID fid_interface = facesIDs_old2new.get(i);

                                 Kokkos::View<kmds::TCellID*> rids;
                                 c_F2R_bis->get(fid, rids);

                                 // there should always be 2 regions
                                 const int mat0 = ma_bis.getMaterial(rids[0]);
                                 const int mat1 = ma_bis.getMaterial(rids[1]);

                                 std::pair<int, int> mat_pair (std::min(mat0, mat1), std::max(mat0, mat1));

                                 (*varInterfaces_int)[fid_interface] = elg3d::Parameters::mat_pair_2_int.at(mat_pair);
                             }
//        );


        kmds::VTKWriter<kmds::Mesh> w(mesh_interface);
        w.write("mesh_interface", kmds::F);

        kmds::ConnectivityHelper ch_interface(&mesh_interface);
        kmds::Connectivity *c_N2F_interface = mesh_interface.createConnectivity(kmds::N2F);
        ch_interface.buildN2F();

        kmds::Connectivity *c_N2N_interface = mesh_interface.createConnectivity(kmds::N2N);
        ch_interface.buildN2N(kmds::F);


        // geom assoc for the interface mesh
        kmds::Variable <std::uintptr_t> *geomassoc_interface_pixels = mesh_interface.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr),
                                                                                                 kmds::KMDS_NODE, "geomassoc_interface_pixels");
        kmds::Variable <std::uintptr_t> *geomassoc_interface_boundingbox = mesh_interface.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr),
                                                                                                                    kmds::KMDS_NODE, "geomassoc_interface_boundingbox");

        const kmds::TCellID nbNodes_interface = nodesInterfaces.getNbElems();
        Kokkos::parallel_for(nbNodes_interface,
                             KOKKOS_LAMBDA(const int i)
                             {
                                 const kmds::TCellID nid = nodesInterfaces.get(i);
                                 const kmds::TCellID nid_interface = nodesIDs_old2new.get(i);

                                 (*geomassoc_interface_pixels)[nid_interface] = (*geomassoc_bis)[nid];
                                 (*geomassoc_interface_boundingbox)[nid_interface] = (*v_boundingbox)[nid];
                             }
        );


        // DEBUG
//        {
//            Kokkos::View<kmds::TCellID *> nids_0;
//            c_N2N_interface->get(0, nids_0);
//
//            std::cout<<"before N2N_interface 0 "<<nids_0.size()<<" "<<nids_0[0]<<" "<<nids_0[1]<<std::endl;
//            const <gmds::cad::GeomEntity* geomEntity = reinterpret_cast<gmds::cad::GeomEntity *> ((*geomassoc_interface_pixels)[0]);
//            std::cout<<"geomEntity->getDim() "<<geomEntity->getDim()<<std::endl;
//        }

        // remove nodes from neighborhood to take into account curve and vertex association
        const kmds::TCellID nbNodes_interface_tmp = mesh_interface.getNbNodes();

        kmds::GrowingView<kmds::TCellID> nodesInterfaces_tmp("NODES", nbNodes_interface_tmp);
        mesh_interface.getNodeIDs(&nodesInterfaces_tmp);

        Kokkos::parallel_for(nbNodes_interface_tmp,
                             KOKKOS_LAMBDA(const int i)
                             {
                                 const kmds::TCellID nid = nodesInterfaces_tmp.get(i);

                                 const gmds::cad::GeomEntity* geomEntity = reinterpret_cast<gmds::cad::GeomEntity *> ((*geomassoc_interface_pixels)[nid]);

                                 // if it is on a surface we keep all the nodes as being adjacent
                                 if(geomEntity->dim() != 2) {

                                     // if a vertex keep no neighbors
                                     if(geomEntity->dim() == 0) {

                                         kmds::TCellID nids[1];
                                         c_N2N_interface->set(nid, nids, 0);

                                     } else {

                                         Kokkos::View<kmds::TCellID*> nids;
                                         c_N2N_interface->get(nid, nids);

                                         kmds::TCellID nids_n[2];
                                         int index = 0;

                                         for(int i_n=0; i_n<nids.size(); i_n++) {

                                             const gmds::cad::GeomEntity* geomEntity_n = reinterpret_cast<gmds::cad::GeomEntity *> ((*geomassoc_interface_pixels)[nids[i_n]]);

                                             if((geomEntity->name() == geomEntity_n->name()) || (geomEntity_n->dim() == 0)) {
                                                 nids_n[index] = nids[i_n];
                                                 index++;
                                             }
                                         }

                                         if(index > 2) {
                                             std::cerr<<"adjusting neighbors failed "<<std::endl;
                                         }

                                         c_N2N_interface->set(nid, nids_n, 2);
                                     }



                                 }

                             }
        );


        // map geom assoc for the interface mesh to the pixels geom assoc
        Kokkos::parallel_for(nbNodes_interface_tmp,
                             KOKKOS_LAMBDA(const int i)
                             {
                                 const kmds::TCellID nid = nodesInterfaces_tmp.get(i);

                                 (*geomassoc_interface_pixels)[nid] = geom_bis2pixels.at((*geomassoc_interface_pixels)[nid]);
                             }
        );


        kmds::VTKWriter<kmds::Mesh> w_tmp(mesh_interface);
        w_tmp.write("mesh_interface_tmp", kmds::F);


        // DEBUG
        {
//            Kokkos::View<kmds::TCellID *> nids_0;
//            c_N2N_interface->get(1372, nids_0);
//            std::cout<<"N2N_interface 1372 "<<nids_0.size()<<" "<<nids_0[0]<<" "<<nids_0[1]<<std::endl;

//            Kokkos::View<kmds::TCellID *> nids_4;
//            c_N2N_interface->get(4, nids_4);
//            std::cout<<"N2N_interface 4 "<<nids_4.size()<<" "<<nids_4[0]<<" "<<nids_4[1]<<std::endl;
//
//            Kokkos::View<kmds::TCellID *> nids_10;
//            c_N2N_interface->get(10, nids_10);
//            std::cout<<"N2N_interface 10 "<<nids_10.size()<<std::endl;
//
//            Kokkos::View<kmds::TCellID *> nids_5;
//            c_N2N_interface->get(5, nids_5);
//            std::cout<<"N2N_interface 5 "<<nids_5.size()<<" "<<nids_5[0]<<" "<<nids_5[1]<<" "<<nids_5[2]<<" "<<nids_5[3]<<std::endl;
        }


        // projection


//        Kokkos::parallel_for(nbNodes_interface_tmp,
//                             KOKKOS_LAMBDA(const int i) {

        // TODO : cannot be called in parallel because gts is not threadsafe
        for(int i=0; i<nbNodes_interface_tmp; i++) {

                                 const kmds::TCellID nid = nodesInterfaces_tmp.get(i);

                                 gmds::math::Point pt = mesh_interface.getNodeLocation(nid);

                                 const std::uintptr_t ge_pixels_ptr =  (*geomassoc_interface_pixels)[nid];
                                 if(ge_pixels_ptr == reinterpret_cast<std::uintptr_t>(nullptr)) {

                                     // default projection if no info on specific geom entity
                                     surfaceGeomServices_pixels.project(pt);
                                 } else {

                                     gmds::cad::GeomEntity *ge_pixels = reinterpret_cast<gmds::cad::GeomEntity *> (ge_pixels_ptr);

                                     switch(ge_pixels->dim()) {
                                         case 0 : {
                                             gmds::cad::FACPoint *vert = reinterpret_cast<gmds::cad::FACPoint *> (ge_pixels);
                                             vert->project(pt);
                                             break;
                                         }
                                         case 1 : {
                                             gmds::cad::FACCurve *curv = reinterpret_cast<gmds::cad::FACCurve *> (ge_pixels);
                                             curveGeomServices_pixels.project(curv, pt);
                                             break;
                                         }
                                         case 2 : {
                                             gmds::cad::FACSurface *surf = reinterpret_cast<gmds::cad::FACSurface *> (ge_pixels);
                                             surfaceGeomServices_pixels.project(surf, pt);
                                             break;
                                         }
                                         default:
                                             exit(-1);
                                     }
                                 }

                                 // computational domain geom association
                                 const std::uintptr_t ge_bb_ptr =  (*geomassoc_interface_boundingbox)[nid];
                                 if(ge_bb_ptr != reinterpret_cast<std::uintptr_t>(nullptr)) {
                                     gmds::cad::GeomEntity *ge_bb = reinterpret_cast<gmds::cad::GeomEntity *> (ge_bb_ptr);

                                     ge_bb->project(pt);
                                 }

            mesh_interface.setNodeLocation(nid, pt);
                             }
//        );

        kmds::VTKWriter<kmds::Mesh> w_proj(mesh_interface);
        w_proj.write("mesh_interface_proj", kmds::F);


        // smooth
        elg3d::smartLaplacian_interface_fromF_doubleGeom_3D(30,
                                                            &mesh_interface,
                                                            c_N2F_interface,
                                                            c_N2N_interface,
                                                            geomassoc_interface_pixels,
                                                            geomassoc_interface_boundingbox,
                                                            &surfaceGeomServices_pixels,
                                                            &curveGeomServices_pixels);


        w_tmp.write("mesh_interface_aftersmooth", kmds::F);
    }

    //==========================================================================


    Kokkos::finalize();
    return 0;
}