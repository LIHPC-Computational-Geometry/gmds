/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 4/1/19.
//
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/frame3d/FieldGenerator.h>
#include <gmds/frame3d/OpenNLFieldSolverStrategy.h>
#include <gmds/frame3d/PointGenerator.h>
#include <gmds/frame3d/PointConnectionBuilder.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <gmds/io/VTKWriter.h>

using namespace gmds;
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    std::cout << "==== 3D Frame Field Generator====" << std::endl;
    Mesh m(MeshModel(DIM3 | R | F | E | N |
                     R2N | F2N | E2N | R2F | F2R |
                     F2E | E2F | R2E | N2R | N2F | N2E));

    //==================================================================
    // MESH READING
    //==================================================================
    std::cout << "Reading " << std::endl;
    std::string fIn, fOut;
    fOut = fIn;
    double spacing = 1;
    if (argc < 3) {
        std::cout << "Require two paramaters : input tetrahedral mesh (.vtk) and expected size (1 by default)"
                  << std::endl;
        throw gmds::GMDSException("Wrong parameters");
    }
    if (argc >= 3) {
        fIn = std::string(argv[1]);
        fOut = ".";
    }
    if (argc == 3) {
        fIn = std::string(argv[1]);
        spacing = atof(std::string(argv[2]).c_str());
     //   fOut = ".";
        std::cout << "INPUT FILE: " << fIn << std::endl;
        std::cout << "OUTPUT DIR: " << fOut << std::endl;
    }


    IGMeshIOService ioService(&m);
    VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N | gmds::R);
    vtkReader.read(fIn);


    MeshDoctor doctor(&m);
    doctor.buildFacesAndR2F();
    doctor.buildEdgesAndX2E();
    doctor.updateUpwardConnectivity();


    ParamsGlobal pg;
    pg.algo_choice=ParamsGlobal::STABLE_HEX_GENERATION;
    pg.start_from=ParamsGlobal::FF_GEN;
    pg.stop_at=ParamsGlobal::FF_SMOOTH;
    pg.with_debug_files=true;
    pg.with_quad_bnd_constraint=false;
    pg.mesh_file=fIn;
    pg.output_dir=fOut;

    ParamsFrameField pf;
    pf.solver_type=ParamsFrameField::OPENNL;
    pf.with_cotangent_weights=false;
    pf.with_smoothing=false;
    pf.smoothing_algo=ParamsFrameField::RAY;
    pf.smoothing_nb_iter=10;
    pf.smoothing_epsilon=1e-6;
    pf.with_mesh_adaptation=false;
    pf.premeshing_sing_lines=false;

    ParamsMark pm;
    pm.mark_node_on_surf = m.newMark<Node>();
    pm.mark_node_on_curv = m.newMark<Node>();
    pm.mark_node_on_pnt  = m.newMark<Node>();
    pm.mark_node_isolated= m.newMark<Node>();

    pm.mark_node_frame   = m.newMark<Node>();
    pm.mark_node_hard    = m.newMark<Node>();

    pm.mark_edge_on_surf = m.newMark<Edge>();
    pm.mark_edge_on_curv = m.newMark<Edge>();

    pm.mark_face_on_surf = m.newMark<Face>();

    BoundaryOperator boundaryOp(&m);

    if (!boundaryOp.isValid()) {
        std::cout << "Invalid model for boundary operations" << std::endl;
        throw GMDSException("Invalid model for boundary operations");
    }

//==================================================================
// Mark boundary cells
    boundaryOp.markCellOnGeometry(pm.mark_face_on_surf, pm.mark_edge_on_surf,
                                  pm.mark_node_on_surf,
                                  pm.mark_edge_on_curv, pm.mark_node_on_curv,
                                  pm.mark_node_on_pnt,
                                  pm.mark_node_isolated);

//==================================================================
// COLOR SURFACE AND CURVE NODES AND ASSIGN BND NORMALS

// Color names used to name the mesh variables are defined in
// the boundary operator class
    boundaryOp.colorEdges(pm.mark_edge_on_curv, pm.mark_node_on_pnt);

    boundaryOp.colorNodes(pm.mark_node_on_pnt);

    std::map<gmds::TCellID, gmds::math::Vector3d> bnd_normals;

    for (auto n_id : m.nodes()) {
        Node n = m.get<Node>(n_id);
        if (m.isMarked(n, pm.mark_node_on_surf)) {
        math::Vector3d nv = boundaryOp.getOutputNormalOfABoundaryNode(n);
        bnd_normals[n.id()] = math::Vector3d({nv.X(), nv.Y(), nv.Z()});
        }
    }

// now we color nodes on curves and surfaces
    Variable<int>* color_f = m.getVariable<int, GMDS_FACE>("BND_SURFACE_COLOR");
    Variable<int>* color_c = m.getVariable<int, GMDS_EDGE>("BND_CURVE_COLOR");

    Variable<int>* color_nf = 0;
    try {
        color_nf = m.newVariable<int, GMDS_NODE>("BND_SURFACE_COLOR");
    } catch (GMDSException& e) {
        color_nf = m.getVariable<int,GMDS_NODE>("BND_SURFACE_COLOR");
    }
    Variable<int>* color_nc = 0;
    try {
        color_nc = m.newVariable<int,GMDS_NODE>("BND_CURVE_COLOR");
    } catch (GMDSException& e) {
        color_nc = m.getVariable<int,GMDS_NODE>("BND_CURVE_COLOR");
    }

    for (auto f_id: m.faces()) {
        Face f = m.get<Face>(f_id);

// only faces on surface are of interest
        if (!m.isMarked(f, pm.mark_face_on_surf))
            continue;

        std::vector<Node> f_nodes = f.get<Node>();
        for (auto ni : f_nodes) {
            if (m.isMarked(ni, pm.mark_node_on_surf) &&
                !m.isMarked(ni, pm.mark_node_on_curv) &&
                !m.isMarked(ni, pm.mark_node_on_pnt)) {
                (*color_nf)[ni.id()] = (*color_f)[f.id()];
            }
        }
    }
    for (auto e_id:m.edges()){
        Edge e = m.get<Edge>(e_id);
// only edges on surface are of interest
        if (!m.isMarked(e, pm.mark_edge_on_curv))
            continue;

        std::vector<Node> e_nodes = e.get<Node>();
        for (auto ni : e_nodes) {
            if (m.isMarked(ni, pm.mark_node_on_curv) &&
                !m.isMarked(ni, pm.mark_node_on_pnt)) {
                (*color_nc)[ni.id()] = (*color_c)[e.id()];
            }
        }
    }


    Variable<double>* cot_w = m.newVariable<double,GMDS_EDGE>("cot_weight");
// All weights equal to 1
    for (auto e_id:m.edges()){

        (*cot_w)[e_id] = 1;
    }

    OpenNLFieldSolverStrategy* solver = new OpenNLFieldSolverStrategy();
    FieldGenerator ffg(solver, &m, pg, pf, pm);
    ffg.execute();


    Variable<math::Chart>* ch = m.getVariable<math::Chart, GMDS_NODE>( "SHChart" );
    Variable<math::AxisAngleRotation>* axis_angle = m.newVariable<math::AxisAngleRotation, GMDS_NODE>("rotation_field");
    for(auto n_id:m.nodes()){
        math::AxisAngleRotation aar(ch->value(n_id));
        axis_angle->set(n_id, aar);
    }
    std::cout<<"========================================"<<std::endl;
    std::cout<<"> Start point generation "<<std::endl;
    //==================================================================
    //We call the point generation algorithm
    PointGenerator ptg (&m,
                        pg,
                        bnd_normals,
                        pm,
                        spacing,
                        0.35);

    ptg.execute();

    std::cout<<"NB POINTS = "<<ptg.points().size()<<std::endl;
    //  writePoints(ptg->points());

    //extraction of the edges
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
    pcb.execute();



    ////////
    MeshModel model(DIM3 | R  | F | N | R2N | F2N);
    Mesh point_mesh(model);

    std::vector<std::pair<int, int>> AEdges;
    pcb.getEdges(AEdges);

    std::vector<std::vector<int>> AHexes;
    pcb.getHexes(AHexes);



    Variable<int>* node_on_pnt = point_mesh.newVariable<int, GMDS_NODE>("BND_VERTEX_COLOR"  );
    Variable<int>* node_on_crv = point_mesh.newVariable<int, GMDS_NODE>("BND_CURVE_COLOR"  );
    Variable<int>* node_on_srf = point_mesh.newVariable<int, GMDS_NODE>("BND_SURFACE_COLOR");

    node_on_pnt->setValuesTo(-1);
    node_on_crv->setValuesTo(-1);
    node_on_srf->setValuesTo(-1);


    const std::vector<int>& pointClassification   = ptg.pointClassification();
    const std::vector<int>& pointCurveNumbering   = ptg.pointCurveNumbering();
    const std::vector<int>& pointSurfaceNumbering = ptg.pointSurfaceNumbering();

    unsigned int cpt = 0;
    for(auto pi:ptg.points()){

        if(pointClassification[cpt] == 0)
        {
          node_on_pnt->set(cpt, 0);
        }
        else if(pointClassification[cpt] == 1)
        {
          node_on_crv->set(cpt, pointCurveNumbering[cpt]);
        }
        else if(pointClassification[cpt] == 2)
        {
          node_on_srf->set(cpt, pointSurfaceNumbering[cpt]);
        }
        Node ni = point_mesh.newNode(pi);
        point_mesh.newTet(ni,ni,ni,ni);
        cpt++;
    }


    for(auto const e : AEdges)
    {
      point_mesh.newTriangle(e.first, e.second, e.second);
    }

    for(auto const h : AHexes)
    {
      math::Point coord0 = point_mesh.get<Node>(h[0]).point();
      math::Point coord1 = point_mesh.get<Node>(h[1]).point();
      math::Point coord2 = point_mesh.get<Node>(h[2]).point();
      math::Point coord3 = point_mesh.get<Node>(h[3]).point();
      math::Point coord4 = point_mesh.get<Node>(h[4]).point();
      math::Point coord5 = point_mesh.get<Node>(h[5]).point();
      math::Point coord6 = point_mesh.get<Node>(h[6]).point();
      math::Point coord7 = point_mesh.get<Node>(h[7]).point();

      math::Point coordCenterHex = 1.0 / 8.0 * (coord0 + coord1 + coord2 +
                                            coord3 + coord4 + coord5 +
                                            coord6 + coord7);
      math::Vector3d n = coordCenterHex - coord0;
      math::Vector3d v0 = coord1 - coord0;
      math::Vector3d v1 = coord2 - coord1;
      math::Vector3d v = v0.cross(v1);

      if(v.dot(n) > 0.0){

        point_mesh.newHex(h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7]);
      }
      else{
        point_mesh.newHex(h[0], h[3], h[2], h[1], h[4], h[7], h[6], h[5]);
      }

    }

    IGMeshIOService ioService2(&point_mesh);
    VTKWriter writer2(&ioService2);
    writer2.setCellOptions(gmds::N|gmds::R|gmds::F);
    writer2.setDataOptions(gmds::N|gmds::R|gmds::F);
    writer2.write("GeneratedPoints.vtk");


    m.unmarkAll<Node>(pm.mark_node_on_surf );
    m.unmarkAll<Node>(pm.mark_node_on_curv );
    m.unmarkAll<Node>(pm.mark_node_on_pnt  );
    m.unmarkAll<Node>(pm.mark_node_isolated);
    m.unmarkAll<Node>(pm.mark_node_frame   );
    m.unmarkAll<Node>(pm.mark_node_hard    );
    m.unmarkAll<Edge>(pm.mark_edge_on_surf );
    m.unmarkAll<Edge>(pm.mark_edge_on_curv );
    m.unmarkAll<Face>(pm.mark_face_on_surf );

    m.freeMark<Node>(pm.mark_node_on_surf );
    m.freeMark<Node>(pm.mark_node_on_curv );
    m.freeMark<Node>(pm.mark_node_on_pnt  );
    m.freeMark<Node>(pm.mark_node_isolated);
    m.freeMark<Node>(pm.mark_node_frame   );
    m.freeMark<Node>(pm.mark_node_hard    );
    m.freeMark<Edge>(pm.mark_edge_on_surf );
    m.freeMark<Edge>(pm.mark_edge_on_curv );
    m.freeMark<Face>(pm.mark_face_on_surf );

}
