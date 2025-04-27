//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/frame3d/PointGenerator.h>
#include <gmds/frame3d/PGPComputing.h>
#include <gmds/frame3d/PGPStructPointExtraction.h>
#include <gmds/frame3d/PointConnectionBuilder.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <gmds/io/VTKWriter.h>
#include <gmds/utils/Log.h>
#include <unit_test_config.h>
#include <chrono>
using namespace std::chrono;


using namespace gmds;


/*----------------------------------------------------------------------------*/
void setup_point_generation_tests(const std::string& AVtkFile,
                                  Mesh* AMesh, ParamsGlobal* AParamGlobal,
                                  ParamsFrameField* AParamFF,
                                  ParamsMark* AParamMark,
                                  std::map<gmds::TCellID, gmds::math::Vector3d>& ABndNormals)
{
    IGMeshIOService ioService(AMesh);
    VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::R);
    vtkReader.read(AVtkFile);


    MeshDoctor doctor(AMesh);
    doctor.buildFacesAndR2F();
    doctor.buildEdgesAndX2E();
    doctor.updateUpwardConnectivity();

    AParamGlobal->algo_choice=ParamsGlobal::STABLE_HEX_GENERATION;
    AParamGlobal->start_from=ParamsGlobal::FF_GEN;
    AParamGlobal->stop_at=ParamsGlobal::FF_SMOOTH;
    AParamGlobal->with_debug_files=true;
    AParamGlobal->with_quad_bnd_constraint=false;
    AParamGlobal->mesh_file=AVtkFile;
    AParamGlobal->output_dir=".";

    AParamFF->solver_type=ParamsFrameField::OPENNL;
    AParamFF->with_cotangent_weights=false;
    AParamFF->with_smoothing=false;
    AParamFF->smoothing_algo=ParamsFrameField::RAY;
    AParamFF->smoothing_nb_iter=10;
    AParamFF->smoothing_epsilon=1e-6;
    AParamFF->with_mesh_adaptation=false;
    AParamFF->premeshing_sing_lines=false;

    AParamMark->mark_node_on_surf = AMesh->newMark<Node>();
    AParamMark->mark_node_on_curv = AMesh->newMark<Node>();
    AParamMark->mark_node_on_pnt  = AMesh->newMark<Node>();
    AParamMark->mark_node_isolated= AMesh->newMark<Node>();

    AParamMark->mark_node_frame   = AMesh->newMark<Node>();
    AParamMark->mark_node_hard    = AMesh->newMark<Node>();

    AParamMark->mark_edge_on_surf = AMesh->newMark<Edge>();
    AParamMark->mark_edge_on_curv = AMesh->newMark<Edge>();

    AParamMark->mark_face_on_surf = AMesh->newMark<Face>();

    BoundaryOperator boundaryOp(AMesh);
    if (!boundaryOp.isValid()) {
        std::cout << "Invalid model for boundary operations" << std::endl;
        throw GMDSException("Invalid model for boundary operations");
    }

//==================================================================
// Mark boundary cells
    boundaryOp.markCellOnGeometry(AParamMark->mark_face_on_surf,
                                  AParamMark->mark_edge_on_surf,
                                  AParamMark->mark_node_on_surf,
                                  AParamMark->mark_edge_on_curv,
                                  AParamMark->mark_node_on_curv,
                                  AParamMark->mark_node_on_pnt,
                                  AParamMark->mark_node_isolated);

//==================================================================
// COLOR SURFACE AND CURVE NODES AND ASSIGN BND NORMALS

// Color names used to name the mesh variables are defined in
// the boundary operator class
    boundaryOp.colorEdges(AParamMark->mark_edge_on_curv,
                          AParamMark->mark_node_on_pnt);

    boundaryOp.colorNodes(AParamMark->mark_node_on_pnt);

    for (auto n_id :AMesh->nodes()) {
        Node n = AMesh->get<Node>(n_id);
        if (AMesh->isMarked(n, AParamMark->mark_node_on_surf)) {
        }
    }
    for (auto n_id : AMesh->nodes()) {
        Node n = AMesh->get<Node>(n_id);
        if (AMesh->isMarked(n, AParamMark->mark_node_on_surf)) {
            math::Vector nv= boundaryOp.getOutputNormalOfABoundaryNode(n);
            ABndNormals[n_id]=math::Vector3d({nv.X(), nv.Y(), nv.Z()});
        }
    }

    // now we color nodes on curves and surfaces
    Variable<int>* color_f = AMesh->getVariable<int, GMDS_FACE>("BND_SURFACE_COLOR");
    Variable<int>* color_c = AMesh->getVariable<int, GMDS_EDGE>("BND_CURVE_COLOR");

    Variable<int>* color_nf = AMesh->getOrCreateVariable<int,GMDS_NODE>("BND_SURFACE_COLOR");
    Variable<int>* color_nc = AMesh->getOrCreateVariable<int,GMDS_NODE>("BND_CURVE_COLOR");

    for (auto f_id: AMesh->faces()) {
        Face f = AMesh->get<Face>(f_id);

        // only faces on surface are of interest
        if (!AMesh->isMarked(f, AParamMark->mark_face_on_surf))
            continue;

        std::vector<Node> f_nodes = f.get<Node>();
        for (auto ni : f_nodes) {
            if (AMesh->isMarked(ni, AParamMark->mark_node_on_surf) &&
                !AMesh->isMarked(ni, AParamMark->mark_node_on_curv) &&
                !AMesh->isMarked(ni, AParamMark->mark_node_on_pnt)) {
                (*color_nf)[ni.id()] = (*color_f)[f.id()];
            }
        }
    }
    for (auto e_id:AMesh->edges()){
        Edge e = AMesh->get<Edge>(e_id);
        // only edges on surface are of interest
        if (!AMesh->isMarked(e, AParamMark->mark_edge_on_curv))
            continue;

        std::vector<Node> e_nodes = e.get<Node>();
        for (auto ni : e_nodes) {
            if (AMesh->isMarked(ni, AParamMark->mark_node_on_curv) &&
                !AMesh->isMarked(ni, AParamMark->mark_node_on_pnt)) {
                (*color_nc)[ni.id()] = (*color_c)[e.id()];
            }
        }
    }

    Variable<double>* cot_w = AMesh->newVariable<double,GMDS_EDGE>("cot_weight");
    // All weights equal to 1
    for (auto e_id:AMesh->edges()){

        (*cot_w)[e_id] = 1;
    }

}
/*----------------------------------------------------------------------------*/
void tear_down_point_generation_tests(Mesh* AMesh,
                                      ParamsMark* AParamsMark) {
    AMesh->unmarkAll<Node>(AParamsMark->mark_node_on_surf );
    AMesh->unmarkAll<Node>(AParamsMark->mark_node_on_curv );
    AMesh->unmarkAll<Node>(AParamsMark->mark_node_on_pnt  );
    AMesh->unmarkAll<Node>(AParamsMark->mark_node_isolated);
    AMesh->unmarkAll<Node>(AParamsMark->mark_node_frame   );
    AMesh->unmarkAll<Node>(AParamsMark->mark_node_hard    );
    AMesh->unmarkAll<Edge>(AParamsMark->mark_edge_on_surf );
    AMesh->unmarkAll<Edge>(AParamsMark->mark_edge_on_curv );
    AMesh->unmarkAll<Face>(AParamsMark->mark_face_on_surf );

    AMesh->freeMark<Node>(AParamsMark->mark_node_on_surf );
    AMesh->freeMark<Node>(AParamsMark->mark_node_on_curv );
    AMesh->freeMark<Node>(AParamsMark->mark_node_on_pnt  );
    AMesh->freeMark<Node>(AParamsMark->mark_node_isolated);
    AMesh->freeMark<Node>(AParamsMark->mark_node_frame   );
    AMesh->freeMark<Node>(AParamsMark->mark_node_hard    );
    AMesh->freeMark<Edge>(AParamsMark->mark_edge_on_surf );
    AMesh->freeMark<Edge>(AParamsMark->mark_edge_on_curv );
    AMesh->freeMark<Face>(AParamsMark->mark_face_on_surf );

}
/*----------------------------------------------------------------------------*/
TEST(PointGeneratorTestSuite, test_point_generator)
{
    Log::mng().clear();
    LogStream ls(&std::cout);
    Log::mng().addStream(ls);
    Mesh m(MeshModel(DIM3 | R | F | E | N |
    R2N | F2N | E2N | R2F | F2R |
    F2E | E2F | R2E | N2R | N2F | N2E));
    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/B0.vtk";

    ParamsGlobal pg;
    ParamsFrameField pf;
    ParamsMark pm;
    std::map<gmds::TCellID, gmds::math::Vector3d> bnd_normals;

    setup_point_generation_tests(vtk_file,&m,&pg,&pf, &pm, bnd_normals);


    OpenNLFieldSolverStrategy* solver = new OpenNLFieldSolverStrategy();
    FieldGenerator ffg(solver, &m, pg, pf, pm);
    ffg.execute();

    Variable<math::Chart>* ch = m.getVariable<math::Chart, GMDS_NODE>( "SHChart" );
    Variable<math::AxisAngleRotation>* axis_angle = m.newVariable<math::AxisAngleRotation, GMDS_NODE>("rotation_field");
    for(auto n_id:m.nodes()){
        math::AxisAngleRotation aar(ch->value(n_id));
        axis_angle->set(n_id, aar);
    }
    //==================================================================
    //We call the point generation algorithm
    PointGenerator ptg (&m,
                        pg,
                        bnd_normals,
                        pm,
                        1,
                        0.35);

    ptg.execute();

    std::cout<<"NB POINTS = "<<ptg.points().size()<<std::endl;

    MeshModel model(DIM3 | R  | N | R2N );
    Mesh point_mesh(model);
    for(auto pi:ptg.points()){
        Node ni = point_mesh.newNode(pi);
        point_mesh.newTet(ni,ni,ni,ni);
    }


    IGMeshIOService ioService2(&point_mesh);
    VTKWriter writer2(&ioService2);
    writer2.setCellOptions(gmds::N|gmds::R);
    writer2.write("GeneratedPoints.vtk");
    ASSERT_EQ(ptg.points().size(),372);
    tear_down_point_generation_tests(&m, &pm);
}
/*----------------------------------------------------------------------------*/
TEST(PointGeneratorTestSuite, test_hex_extraction)
{
    Log::mng().clear();
    LogStream ls(&std::cout);
    Log::mng().addStream(ls);
    Mesh m(MeshModel(DIM3 | R | F | E | N |
    R2N | F2N | E2N | R2F | F2R |
    F2E | E2F | R2E | N2R | N2F | N2E));
    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/B0.vtk";

    ParamsGlobal pg;
    ParamsFrameField pf;
    ParamsMark pm;
    std::map<gmds::TCellID, gmds::math::Vector3d> bnd_normals;

    setup_point_generation_tests(vtk_file,&m,&pg,&pf, &pm, bnd_normals);


    OpenNLFieldSolverStrategy* solver = new OpenNLFieldSolverStrategy();
    FieldGenerator ffg(solver, &m, pg, pf, pm);
    ffg.execute();

    Variable<math::Chart>* ch = m.getVariable<math::Chart, GMDS_NODE>( "SHChart" );
    Variable<math::AxisAngleRotation>* axis_angle = m.newVariable<math::AxisAngleRotation, GMDS_NODE>("rotation_field");
    for(auto n_id:m.nodes()){
        math::AxisAngleRotation aar(ch->value(n_id));
        axis_angle->set(n_id, aar);
    }
    //==================================================================
    //We call the point generation algorithm
    PointGenerator ptg (&m,
                        pg,
                        bnd_normals,
                        pm,
                        1,
                        0.35);

    ptg.execute();

    std::cout<<"NB POINTS = "<<ptg.points().size()<<std::endl;

    PointConnectionBuilder pcb(&m,
                               ptg.points(),
                               ptg.charts(),
                               ptg.pointMeshData(),
                               ptg.pointTypes(),
                               ptg.pointClassification(),
                               ptg.pointCurveNumbering(),
                               ptg.pointSurfaceNumbering(),
                               ptg.pointSurfaceNormal());
    pcb.setDebugInfo(true);
    auto start = high_resolution_clock::now();
    pcb.execute();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "Hex structure: "<<duration.count() << std::endl;
    std::vector<std::vector<int> > hexes;
    pcb.getHexes(hexes);
    std::cout<<"NB HEXES = "<<hexes.size()<<std::endl;

    ASSERT_EQ(200, hexes.size());
}
/*----------------------------------------------------------------------------*/
TEST(PointGeneratorTestSuite, test_PGP)
{
    Log::mng().clear();
    LogStream ls(&std::cout);
    Log::mng().addStream(ls);
    Mesh m(MeshModel(DIM3 | R | F | E | N |
                     R2N | F2N | E2N | R2F | F2R |
                     F2E | E2F | R2E | N2R | N2F | N2E));
    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/B0.vtk";

    ParamsGlobal pg;
    ParamsFrameField pf;
    ParamsMark pm;
    std::map<gmds::TCellID, gmds::math::Vector3d> bnd_normals;

    setup_point_generation_tests(vtk_file,&m,&pg,&pf, &pm, bnd_normals);


    OpenNLFieldSolverStrategy* solver = new OpenNLFieldSolverStrategy();
    FieldGenerator ffg(solver, &m, pg, pf, pm);
    ffg.execute();

    Variable<math::Chart>* ch = m.getVariable<math::Chart, GMDS_NODE>( "SHChart" );
    Variable<math::AxisAngleRotation>* axis_angle = m.newVariable<math::AxisAngleRotation, GMDS_NODE>("rotation_field");
    for(auto n_id:m.nodes()){
        math::AxisAngleRotation aar(ch->value(n_id));
        axis_angle->set(n_id, aar);
    }
    //==================================================================
    //We call the point generation algorithm
    PGPComputing ptg (&m,
                        pg,
                        pm,
                        1,
                        0.35);
    ptg.execute();

    PGPStructPointExtraction pte(&m,pg,bnd_normals,pm,ptg.getUI(),1,0.35);
    pte.execute();
    std::cout<<"NB POINTS = "<<pte.points().size()<<std::endl;
    //    writePoints(ptg->points());

    MeshModel model(DIM3 | R  | N | R2N );
    Mesh point_mesh(model);
    for(auto pi:pte.points()){
        Node ni = point_mesh.newNode(pi);
        point_mesh.newTet(ni,ni,ni,ni);
    }


    IGMeshIOService ioService2(&point_mesh);
    VTKWriter writer2(&ioService2);
    writer2.setCellOptions(gmds::N|gmds::R);
    writer2.write("GeneratedPoints.vtk");
    ASSERT_EQ(pte.points().size(),579);
    tear_down_point_generation_tests(&m, &pm);
}