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
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    std::cout << "==== Surface/Curve/Corner indexed ====" << std::endl;
    Mesh m(MeshModel(DIM3 | R | F | E | N |
                     R2N | F2N | E2N | R2F | F2R |
                     F2E | E2F | R2E | N2R | N2F | N2E));

    //==================================================================
    // MESH READING
    //==================================================================
    std::cout << "Reading " << std::endl;
    std::string fIn, fOut;
    double spacing = 1.0;
    if(argc != 3)
    {
        throw gmds::GMDSException("MISSING PARAMETER <file_name> <3D parameter>");
    }
    fIn = std::string(argv[1]);
    spacing = atof(std::string(argv[2]).c_str());

    if (fIn.find('.vtk') == std::string::npos) {
      throw gmds::GMDSException("NOT A .vtk FILE");
    }

    std::string extansion(".vtk");
    std::size_t position = fIn.find(extansion);
    fOut = fIn.substr(0,position) + "_POINT_GENERATED.vtk";
    std::cout << "INPUT FILE: " << fIn << std::endl;
    std::cout << "OUTPUT FILE: " << fOut << std::endl;


    IGMeshIOService ioService(&m);
    VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N | gmds::R);
    vtkReader.setDataOptions(gmds::N);
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
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
    std::map<gmds::TCellID, gmds::math::Vector3d> bnd_normals;
    for (auto n_id : m.nodes()) {
        Node n = m.get<Node>(n_id);
        if (m.isMarked(n, pm.mark_node_on_surf)) {
        math::Vector3d nv = boundaryOp.getOutputNormalOfABoundaryNode(n);
        std::cout << nv << std::endl;
        bnd_normals[n.id()] = math::Vector3d(nv.X(), nv.Y(), nv.Z());
        }
    }
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

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
        /*PointConnectionBuilder pcb(&m,
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
        std::cout << "pcb.execute(); done " << std::endl;*/


        ////////
        MeshModel model(DIM3 | R  | F | N | R2N | F2N);
        Mesh point_mesh(model);

        /*std::vector<std::pair<int, int>> AEdges;
        pcb.getEdges(AEdges);

        std::vector<std::vector<int>> AHexes;
        pcb.getHexes(AHexes);*/



        Variable<int>* node_on_pnt = point_mesh.newVariable<int, GMDS_NODE>("BND_VERTEX_COLOR"  );
        Variable<int>* node_on_crv = point_mesh.newVariable<int, GMDS_NODE>("BND_CURVE_COLOR"  );
        Variable<int>* node_on_srf = point_mesh.newVariable<int, GMDS_NODE>("BND_SURFACE_COLOR");

        node_on_pnt->setValuesTo(-1);
        node_on_crv->setValuesTo(-1);
        node_on_srf->setValuesTo(-1);


        const std::vector<int>& pointClassification   = ptg.pointClassification();
        const std::vector<int>& pointCurveNumbering   = ptg.pointCurveNumbering();
        const std::vector<int>& pointSurfaceNumbering = ptg.pointSurfaceNumbering();

        std::cout << "pointClassification.size() --> " << pointClassification.size() << std::endl;
        std::cout << "pointCurveNumbering.size() --> " << pointCurveNumbering.size() << std::endl;
        std::cout << "pointSurfaceNumbering.size() --> " << pointSurfaceNumbering.size() << std::endl;

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
        IGMeshIOService ioService3(&point_mesh);
        VTKWriter writer3(&ioService3);
        writer3.setCellOptions(gmds::N | gmds::R);
        writer3.setDataOptions(gmds::N | gmds::R);
        writer3.write("POINT_GENERATED.vtk");
        return 0 ;


        /*for(auto const e : AEdges)
        {
          point_mesh.newTriangle(e.first, e.second, e.second);
        }

        for(auto const h : AHexes)
        {
          point_mesh.newHex(h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7]);
        }

        IGMeshIOService ioService2(&point_mesh);
        VTKWriter writer2(&ioService2);
        writer2.setCellOptions(gmds::N|gmds::R|gmds::F);
        writer2.setDataOptions(gmds::N|gmds::R|gmds::F);
        writer2.write("CUBE_TRANSFORED_POINT_GENERATED.vtk");*/


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
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

}

/*----------------------------------------------------------------------------*/
