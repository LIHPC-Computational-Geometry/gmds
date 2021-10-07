//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/frame3d/FieldGenerator.h>
#include <gmds/frame3d/OpenNLFieldSolverStrategy.h>
#include <gmds/frame3d/SingularityLineHelper.h>
#include <gmds/math/Cross2D.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <gmds/io/VTKWriter.h>
#include <gmds/utils/Log.h>
#include <unit_test_config.h>
using namespace gmds;
/*----------------------------------------------------------------------------*/
void setup_frame3d_tests(const std::string& AVtkFile,
                      Mesh* AMesh, ParamsGlobal* AParamGlobal,
                      ParamsFrameField* AParamFF,
                      ParamsMark* AParamMark)
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

// now we color nodes on curves and surfaces
    Variable<int>* color_f = AMesh->getVariable<int, GMDS_FACE>("BND_SURFACE_COLOR");
    Variable<int>* color_c = AMesh->getVariable<int, GMDS_EDGE>("BND_CURVE_COLOR");

    Variable<int>* color_nf = AMesh->getOrCreateVariable<int, GMDS_NODE>("BND_SURFACE_COLOR");
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
void tear_down_frame3d_tests(Mesh* AMesh,
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
TEST(Frame3dTestClass, test_B0)
{
    Log::mng().clear();
    Mesh m(MeshModel(DIM3 | R | F | E | N |
                     R2N | F2N | E2N | R2F | F2R |
                     F2E | E2F | R2E | N2R | N2F | N2E));
    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/B0.vtk";

    ParamsGlobal pg;
    ParamsFrameField pf;
    ParamsMark pm;

    setup_frame3d_tests(vtk_file,&m,&pg,&pf, &pm);

    OpenNLFieldSolverStrategy* solver = new OpenNLFieldSolverStrategy();
    FieldGenerator ffg(solver, &m, pg, pf, pm);
    ffg.execute();

    ASSERT_EQ(ffg.getSingularTetIds().size(),129);
    tear_down_frame3d_tests(&m, &pm);

}
/*----------------------------------------------------------------------------*/
TEST(Frame3dTestClass, test_B0_free_boundary)
{
    Log::mng().clear();
    Mesh m(MeshModel(DIM3 | R | F | E | N |
                     R2N | F2N | E2N | R2F | F2R |
                     F2E | E2F | R2E | N2R | N2F | N2E));
    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/B0.vtk";

    ParamsGlobal pg;
    ParamsFrameField pf;
    ParamsMark pm;

    setup_frame3d_tests(vtk_file,&m,&pg,&pf, &pm);

    Variable<int>* free_bnd = m.newVariable<int,GMDS_NODE>("free_bnd_nodes");
    //We relax the whole boundary
    Variable<int>* color_nf = m.getVariable<int, GMDS_NODE>("BND_SURFACE_COLOR");
    Variable<int>* color_nc = m.getVariable<int,GMDS_NODE>("BND_CURVE_COLOR");

    for(auto n_id:m.nodes()){
        int curv_color =color_nc->value(n_id);
        int surf_color = color_nf->value(n_id);
        if(surf_color==1 || curv_color==1 || curv_color==6 || curv_color==13 || curv_color==14) {
            math::Point p = m.get<Node>(n_id).point();
            if(p.X()<=0)
                free_bnd->set(n_id, 1);
        }
        else
            free_bnd->set(n_id, 0);
    }

    OpenNLFieldSolverStrategy* solver = new OpenNLFieldSolverStrategy();
    FieldGenerator ffg(solver, &m, pg, pf, pm);
    ffg.relaxBoundary(free_bnd);
    ffg.execute();

    IGMeshIOService ioService(&m);
    VTKWriter writer(&ioService);
    writer.setCellOptions(gmds::N|gmds::R);
    writer.setDataOptions(gmds::N|gmds::R);
    //only one singularity line for this case
    ASSERT_EQ(ffg.getSingularTetIds().size(),64);

    tear_down_frame3d_tests(&m, &pm);
}
/*----------------------------------------------------------------------------*/
TEST(Frame3dTestClass, test_B48_free_boundary)
{
    Log::mng().clear();
    Mesh m(MeshModel(DIM3 | R | F | E | N |
                     R2N | F2N | E2N | R2F | F2R |
                     F2E | E2F | R2E | N2R | N2F | N2E));
    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/B48.vtk";

    ParamsGlobal pg;
    ParamsFrameField pf;
    ParamsMark pm;

    setup_frame3d_tests(vtk_file,&m,&pg,&pf, &pm);

    Variable<int>* free_bnd = m.newVariable<int,GMDS_NODE>("free_bnd_nodes");
    //We relax a part of the boundary

    double tol = 0.1;
    double x_min=0.-tol;
    double y_min=-0.5-tol;
    double z_min=-5-tol;
    double x_max=1+tol;
    double y_max=0.5+tol;
    double z_max=0.0+tol;
    for(auto n_id:m.nodes()){
        free_bnd->set(n_id, 0);
        if(m.isMarked<Node>(n_id,pm.mark_node_on_surf) ||
           m.isMarked<Node>(n_id,pm.mark_node_on_curv)||
           m.isMarked<Node>(n_id,pm.mark_node_on_pnt)) {
            math::Point p = m.get<Node>(n_id).point();
            if(p.X()>x_min && p.X()<x_max &&
               p.Y()>y_min && p.Y()<y_max &&
               p.Z()>z_min && p.Z()<z_max)
                free_bnd->set(n_id, 1);
        }

    }
    x_min=-5.-tol;
    y_min=-0.5-tol;
    z_min=0-tol;
    x_max=1+tol;
    y_max=1+tol;
    z_max=1+tol;
    for(auto n_id:m.nodes()){
        if(m.isMarked<Node>(n_id,pm.mark_node_on_surf) ||
           m.isMarked<Node>(n_id,pm.mark_node_on_curv)||
           m.isMarked<Node>(n_id,pm.mark_node_on_pnt)) {
            math::Point p = m.get<Node>(n_id).point();
            if(p.X()>x_min && p.X()<x_max &&
               p.Y()>y_min && p.Y()<y_max &&
               p.Z()>z_min && p.Z()<z_max)
                free_bnd->set(n_id, 1);
        }

    }


    OpenNLFieldSolverStrategy* solver = new OpenNLFieldSolverStrategy();
    FieldGenerator ffg(solver, &m, pg, pf, pm);
    ffg.relaxBoundary(free_bnd);
    ffg.execute();

    IGMeshIOService ioService(&m);
    VTKWriter writer(&ioService);
    writer.setCellOptions(gmds::N|gmds::R);
    writer.setDataOptions(gmds::N|gmds::R);
    //only one singularity line for this case
    ASSERT_EQ(ffg.getSingularTetIds().size(),77);

    tear_down_frame3d_tests(&m, &pm);
}
/*----------------------------------------------------------------------------*/
TEST(Frame3dTestClass, test2)
{
    ASSERT_TRUE(true);
//    // WE WRITE
//   Mesh m(MeshModel(DIM3 | R | F | E | N |
//                     R2N | F2N | E2N | R2F | F2R |
//                     F2E | E2F | R2E | N2R | N2F | N2E));
//
//    std::string vtk_file = ("test_samples/S3.vtk");
//
//
//    IGMeshIOService ioService(&m);
//    VTKReader vtkReader(&ioService);
//    vtkReader.setCellOptions(gmds::N|gmds::R);
//    vtkReader.read(vtk_file);
//
//
//    MeshDoctor doctor(&m);
//    doctor.buildFacesAndR2F();
//    doctor.buildEdgesAndX2E();
//    doctor.updateUpwardConnectivity();
//
//
//    ParamsGlobal pg;
//    pg.algo_choice=ParamsGlobal::STABLE_HEX_GENERATION;
//    pg.start_from=ParamsGlobal::FF_GEN;
//    pg.stop_at=ParamsGlobal::FF_SMOOTH;
//    pg.with_debug_files=true;
//    pg.with_quad_bnd_constraint=false;
//    pg.mesh_file=vtk_file;
//    pg.output_dir=".";
//
//    ParamsFrameField pf;
//    pf.solver_type=ParamsFrameField::OPENNL;
//    pf.with_cotangent_weights=false;
//    pf.with_smoothing=false;
//    pf.smoothing_algo=ParamsFrameField::RAY;
//    pf.smoothing_nb_iter=10;
//    pf.smoothing_epsilon=1e-6;
//    pf.with_mesh_adaptation=false;
//    pf.premeshing_sing_lines=false;
//
//    ParamsMark pm;
//    pm.mark_node_on_surf = m.newMark<Node>();
//    pm.mark_node_on_curv = m.newMark<Node>();
//    pm.mark_node_on_pnt  = m.newMark<Node>();
//    pm.mark_node_isolated= m.newMark<Node>();
//
//    pm.mark_node_frame   = m.newMark<Node>();
//    pm.mark_node_hard    = m.newMark<Node>();
//
//    pm.mark_edge_on_surf = m.newMark<Edge>();
//    pm.mark_edge_on_curv = m.newMark<Edge>();
//
//    pm.mark_face_on_surf = m.newMark<Face>();
//
//    BoundaryOperator boundaryOp(&m);
//
//    if (!boundaryOp.isValid()) {
//        std::cout << "Invalid model for boundary operations" << std::endl;
//        throw GMDSException("Invalid model for boundary operations");
//    }
//
//    //==================================================================
//    // Mark boundary cells
//    boundaryOp.markCellOnGeometry(pm.mark_face_on_surf, pm.mark_edge_on_surf,
//                                  pm.mark_node_on_surf,
//                                  pm.mark_edge_on_curv, pm.mark_node_on_curv,
//                                  pm.mark_node_on_pnt,
//                                  pm.mark_node_isolated);
//
//    //==================================================================
//    // COLOR SURFACE AND CURVE NODES AND ASSIGN BND NORMALS
//
//    // Color names used to name the mesh variables are defined in
//    // the boundary operator class
//    boundaryOp.colorEdges(pm.mark_edge_on_curv, pm.mark_node_on_pnt);
//
//    boundaryOp.colorNodes(pm.mark_node_on_pnt);
//
//    for (auto n_id : m.nodes()) {
//        Node n = m.get<Node>(n_id);
//        if (m.isMarked(n, pm.mark_node_on_surf)) {
//        }
//    }
//
//    // now we color nodes on curves and surfaces
//    Variable<int>* color_f = m.getVariable<int, GMDS_FACE>("BND_SURFACE_COLOR");
//    Variable<int>* color_c = m.getVariable<int, GMDS_EDGE>("BND_CURVE_COLOR");
//
//    Variable<int>* color_nf = 0;
//    try {
//        color_nf = m.newVariable<int, GMDS_NODE>("BND_SURFACE_COLOR");
//    } catch (GMDSException& e) {
//        color_nf = m.getVariable<int,GMDS_NODE>("BND_SURFACE_COLOR");
//    }
//    Variable<int>* color_nc = 0;
//    try {
//        color_nc = m.newVariable<int,GMDS_NODE>("BND_CURVE_COLOR");
//    } catch (GMDSException& e) {
//        color_nc = m.getVariable<int,GMDS_NODE>("BND_CURVE_COLOR");
//    }
//
//    for (auto f_id: m.faces()) {
//        Face f = m.get<Face>(f_id);
//
//        // only faces on surface are of interest
//        if (!m.isMarked(f, pm.mark_face_on_surf))
//            continue;
//
//        std::vector<Node> f_nodes = f.get<Node>();
//        for (auto ni : f_nodes) {
//            if (m.isMarked(ni, pm.mark_node_on_surf) &&
//                !m.isMarked(ni, pm.mark_node_on_curv) &&
//                !m.isMarked(ni, pm.mark_node_on_pnt)) {
//                (*color_nf)[ni.id()] = (*color_f)[f.id()];
//            }
//        }
//    }
//    for (auto e_id:m.edges()){
//        Edge e = m.get<Edge>(e_id);
//        // only edges on surface are of interest
//        if (!m.isMarked(e, pm.mark_edge_on_curv))
//            continue;
//
//        std::vector<Node> e_nodes = e.get<Node>();
//        for (auto ni : e_nodes) {
//            if (m.isMarked(ni, pm.mark_node_on_curv) &&
//                !m.isMarked(ni, pm.mark_node_on_pnt)) {
//                (*color_nc)[ni.id()] = (*color_c)[e.id()];
//            }
//        }
//    }
//
//
//    Variable<double>* cot_w = m.newVariable<double,GMDS_EDGE>("cot_weight");
//    // All weights equal to 1
//    for (auto e_id:m.edges()){
//
//        (*cot_w)[e_id] = 1;
//    }
//
//    OpenNLFieldSolverStrategy* solver = new OpenNLFieldSolverStrategy();
//    FieldGenerator ffg(solver, &m, pg, pf, pm);
//    ffg.execute();
//
//
//    ASSERT_EQ(ffg.getSingularTetIds().size(),366);
//
//    VTKWriter writer(&ioService);
//    writer.setCellOptions(gmds::N|gmds::R);
//    writer.setDataOptions(gmds::N|gmds::R);
//    writer.write("S3_out.vtk");
//
//
//
//
//
//    m.unmarkAll<Node>(pm.mark_node_on_surf );
//    m.unmarkAll<Node>(pm.mark_node_on_curv );
//    m.unmarkAll<Node>(pm.mark_node_on_pnt  );
//    m.unmarkAll<Node>(pm.mark_node_isolated);
//    m.unmarkAll<Node>(pm.mark_node_frame   );
//    m.unmarkAll<Node>(pm.mark_node_hard    );
//    m.unmarkAll<Edge>(pm.mark_edge_on_surf );
//    m.unmarkAll<Edge>(pm.mark_edge_on_curv );
//    m.unmarkAll<Face>(pm.mark_face_on_surf );
//
//    m.freeMark<Node>(pm.mark_node_on_surf );
//    m.freeMark<Node>(pm.mark_node_on_curv );
//    m.freeMark<Node>(pm.mark_node_on_pnt  );
//    m.freeMark<Node>(pm.mark_node_isolated);
//    m.freeMark<Node>(pm.mark_node_frame   );
//    m.freeMark<Node>(pm.mark_node_hard    );
//    m.freeMark<Edge>(pm.mark_edge_on_surf );
//    m.freeMark<Edge>(pm.mark_edge_on_curv );
//    m.freeMark<Face>(pm.mark_face_on_surf );

}
/*----------------------------------------------------------------------------*/
TEST(Frame3dTestClass, test_modif_FF)
{
//// WE WRITE
//    Mesh m(MeshModel(DIM3 | R | F | E | N |
//                     R2N | F2N | E2N | R2F | F2R |
//                     F2E | E2F | R2E | N2R | N2F | N2E));
//
// //   std::string vtk_file = ("test_samples/sweep_ff.vtk");
//   std::string vtk_file = ("test.vtk");
//
//    IGMeshIOService ioService(&m);
//    VTKReader vtkReader(&ioService);
//    vtkReader.setCellOptions(gmds::N|gmds::R);
// //   vtkReader.setDataOptions(gmds::N|gmds::R);
//    vtkReader.read(vtk_file);
//
//
////===================================================
//   /* Variable<math::Chart>* v_ch =m.newVariable<math::Chart, GMDS_NODE>( "SHChart" );
//    Variable<math::Vector3d>* v_X =m.getVariable<math::Vector3d, GMDS_NODE>("FF_X_POS");
//    Variable<math::Vector3d>* v_Y =m.getVariable<math::Vector3d, GMDS_NODE>("FF_Y_POS");
//    Variable<math::Vector3d>* v_Z =m.getVariable<math::Vector3d, GMDS_NODE>("FF_Z_POS");
//
//    for(auto n_id:m.nodes()){
//        math::Chart ci(v_X->value(n_id),v_Y->value(n_id),v_Z->value(n_id));
//        v_ch->set(n_id,ci);
//    }
//*/
//    MeshDoctor doctor(&m);
//    doctor.buildFacesAndR2F();
//    doctor.buildEdgesAndX2E();
//    doctor.updateUpwardConnectivity();
//
//
//
//    ParamsGlobal pg;
//    pg.algo_choice=ParamsGlobal::STABLE_HEX_GENERATION;
//    pg.start_from=ParamsGlobal::FF_GEN;
//    pg.stop_at=ParamsGlobal::FF_SMOOTH;
//    pg.with_debug_files=true;
//    pg.with_quad_bnd_constraint=false;
//    pg.mesh_file=vtk_file;
//    pg.output_dir=".";
//
//    ParamsFrameField pf;
//    pf.solver_type=ParamsFrameField::OPENNL;
//    pf.with_cotangent_weights=false;
//    pf.with_smoothing=false;
//    pf.smoothing_algo=ParamsFrameField::RAY;
//    pf.smoothing_nb_iter=0;
//    pf.smoothing_epsilon=1e-6;
//    pf.with_mesh_adaptation=false;
//    pf.premeshing_sing_lines=false;
//
//    ParamsMark pm;
//    pm.mark_node_on_surf = m.newMark<Node>();
//    pm.mark_node_on_curv = m.newMark<Node>();
//    pm.mark_node_on_pnt  = m.newMark<Node>();
//    pm.mark_node_isolated= m.newMark<Node>();
//
//    pm.mark_node_frame   = m.newMark<Node>();
//    pm.mark_node_hard    = m.newMark<Node>();
//
//    pm.mark_edge_on_surf = m.newMark<Edge>();
//    pm.mark_edge_on_curv = m.newMark<Edge>();
//
//    pm.mark_face_on_surf = m.newMark<Face>();
//
//    BoundaryOperator boundaryOp(&m);
//
//    if (!boundaryOp.isValid()) {
//        std::cout << "Invalid model for boundary operations" << std::endl;
//        throw GMDSException("Invalid model for boundary operations");
//    }
//
////==================================================================
//// Mark boundary cells
//    boundaryOp.markCellOnGeometry(pm.mark_face_on_surf, pm.mark_edge_on_surf,
//                                  pm.mark_node_on_surf,
//                                  pm.mark_edge_on_curv, pm.mark_node_on_curv,
//                                  pm.mark_node_on_pnt,
//                                  pm.mark_node_isolated);
//
////==================================================================
//// COLOR SURFACE AND CURVE NODES AND ASSIGN BND NORMALS
//
//// Color names used to name the mesh variables are defined in
//// the boundary operator class
//    boundaryOp.colorEdges(pm.mark_edge_on_curv, pm.mark_node_on_pnt);
//
//    boundaryOp.colorNodes(pm.mark_node_on_pnt);
//
//    for (auto n_id : m.nodes()) {
//        Node n = m.get<Node>(n_id);
//        if (m.isMarked(n, pm.mark_node_on_surf)) {
//        }
//    }
//
//// now we color nodes on curves and surfaces
//    Variable<int>* color_f = m.getVariable<int, GMDS_FACE>("BND_SURFACE_COLOR");
//    Variable<int>* color_c = m.getVariable<int, GMDS_EDGE>("BND_CURVE_COLOR");
//
//    Variable<int>* color_nf = 0;
//    try {
//        color_nf = m.newVariable<int, GMDS_NODE>("BND_SURFACE_COLOR");
//    } catch (GMDSException& e) {
//        color_nf = m.getVariable<int,GMDS_NODE>("BND_SURFACE_COLOR");
//    }
//    Variable<int>* color_nc = 0;
//    try {
//        color_nc = m.newVariable<int,GMDS_NODE>("BND_CURVE_COLOR");
//    } catch (GMDSException& e) {
//        color_nc = m.getVariable<int,GMDS_NODE>("BND_CURVE_COLOR");
//    }
//
//    for (auto f_id: m.faces()) {
//        Face f = m.get<Face>(f_id);
//
//// only faces on surface are of interest
//        if (!m.isMarked(f, pm.mark_face_on_surf))
//            continue;
//
//        std::vector<Node> f_nodes = f.get<Node>();
//        for (auto ni : f_nodes) {
//            if (m.isMarked(ni, pm.mark_node_on_surf) &&
//                !m.isMarked(ni, pm.mark_node_on_curv) &&
//                !m.isMarked(ni, pm.mark_node_on_pnt)) {
//                (*color_nf)[ni.id()] = (*color_f)[f.id()];
//            }
//        }
//    }
//    for (auto e_id:m.edges()){
//        Edge e = m.get<Edge>(e_id);
//// only edges on surface are of interest
//        if (!m.isMarked(e, pm.mark_edge_on_curv))
//            continue;
//
//        std::vector<Node> e_nodes = e.get<Node>();
//        for (auto ni : e_nodes) {
//            if (m.isMarked(ni, pm.mark_node_on_curv) &&
//                !m.isMarked(ni, pm.mark_node_on_pnt)) {
//                (*color_nc)[ni.id()] = (*color_c)[e.id()];
//            }
//        }
//    }
//
//
//    Variable<double>* cot_w = m.newVariable<double,GMDS_EDGE>("cot_weight");
//// All weights equal to 1
//    for (auto e_id:m.edges()){
//
//        (*cot_w)[e_id] = 1;
//    }
//
//    OpenNLFieldSolverStrategy* solver = new OpenNLFieldSolverStrategy();
//    FieldGenerator ffg(solver, &m, pg, pf, pm);
//    ffg.execute();
//
//
//   // ASSERT_EQ(ffg.getSingularTetIds().size(),56);//47);
//
//    frame3d::SingularityLineHelper sm(&m,pg,pm);
//    sm.execute();
//
//    Variable<math::Chart>* v_chart = m.getVariable<math::Chart   , GMDS_NODE>("propag_chart");
//    Variable<math::Vector3d>* v_cx = m.newVariable<math::Vector3d, GMDS_NODE>("chart_X");
//    Variable<math::Vector3d>* v_cy = m.newVariable<math::Vector3d, GMDS_NODE>("chart_Y");
//    Variable<math::Vector3d>* v_cz = m.newVariable<math::Vector3d, GMDS_NODE>("chart_Z");
//
//    for(auto n: m.nodes()){
//        v_cx->set(n,(*v_chart)[n].X());
//        v_cy->set(n,(*v_chart)[n].Y());
//        v_cz->set(n,(*v_chart)[n].Z());
//    }
//
//
//    //AND WE START AGAIN THE FRAME FIELD GENERATION WITH EXTRA CONSTRAINT
//    Variable<math::Chart>* v_chart_constraint = m.getVariable<math::Chart,GMDS_NODE>("chart_constraint");
//    Variable<int>* v_has_constraint = m.getVariable<int,GMDS_NODE>("has_constraint");
//
//    Variable<math::Chart>* chart_field    = m.getVariable<math::Chart, GMDS_NODE>( "SHChart" );
//    Variable<int>* tet_to_update = m.newVariable<int,GMDS_REGION>("to_update");
//
//    std::set<TCellID> nids_to_update;
//
//    sm.updateSmoothingData(nids_to_update,
//                           chart_field,
//                           v_has_constraint,
//                           tet_to_update);
//
//    VTKWriter writer(&ioService);
//    writer.setCellOptions(gmds::N|gmds::R|gmds::F);
//    writer.setDataOptions(gmds::N|gmds::R|gmds::F);
//    writer.write("sweep_ff_out.vtk");
//
//
//
//    sm.initFrameSmoothing(nids_to_update,
//                          v_has_constraint,
//                          tet_to_update);
//    //We update the
//    sm.updateSmoothingData(nids_to_update,
//                           chart_field,
//                           v_has_constraint,
//                           tet_to_update);
//    VTKWriter writer2(&ioService);
//    writer2.setCellOptions(gmds::N|gmds::R|gmds::F);
//    writer2.setDataOptions(gmds::N|gmds::R|gmds::F);
//    writer2.write("sweep_ff_out_last.vtk");
//
//    for(auto t:m.regions()){
//        if((*tet_to_update)[t]==1){
//            Region r = m.get<Region>(t);
//            std::vector<TCellID> node_ids=r.getIDs<Node>();
//            for(auto i:node_ids){
//                math::Chart ci = (*chart_field)[i];
//                (*v_chart_constraint)[i]=ci;
//                (*v_has_constraint)[i]=1;
//            }
//        }
//    }
//
//    pf.smoothing_nb_iter=200;
//
//      OpenNLFieldSolverStrategy* solver2 = new OpenNLFieldSolverStrategy();
//      FieldGenerator ffg2(solver2, &m, pg, pf, pm);
//
//      ffg2.addHardConstraints(v_has_constraint, v_chart_constraint);
//      ffg2.execute();
//
//
//    m.unmarkAll<Node>(pm.mark_node_on_surf );
//    m.unmarkAll<Node>(pm.mark_node_on_curv );
//    m.unmarkAll<Node>(pm.mark_node_on_pnt  );
//    m.unmarkAll<Node>(pm.mark_node_isolated);
//    m.unmarkAll<Node>(pm.mark_node_frame   );
//    m.unmarkAll<Node>(pm.mark_node_hard    );
//    m.unmarkAll<Edge>(pm.mark_edge_on_surf );
//    m.unmarkAll<Edge>(pm.mark_edge_on_curv );
//    m.unmarkAll<Face>(pm.mark_face_on_surf );
//
//    m.freeMark<Node>(pm.mark_node_on_surf );
//    m.freeMark<Node>(pm.mark_node_on_curv );
//    m.freeMark<Node>(pm.mark_node_on_pnt  );
//    m.freeMark<Node>(pm.mark_node_isolated);
//    m.freeMark<Node>(pm.mark_node_frame   );
//    m.freeMark<Node>(pm.mark_node_hard    );
//    m.freeMark<Edge>(pm.mark_edge_on_surf );
//    m.freeMark<Edge>(pm.mark_edge_on_curv );
//    m.freeMark<Face>(pm.mark_face_on_surf );

}
