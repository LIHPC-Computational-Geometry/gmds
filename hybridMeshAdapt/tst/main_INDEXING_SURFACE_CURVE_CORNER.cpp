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
    if(argc != 2)
    {
        throw gmds::GMDSException("NO INPUT FILE");
    }
    fIn = std::string(argv[1]);
    if (fIn.find('.vtk') == std::string::npos) {
      throw gmds::GMDSException("NOT A .vtk FILE");
    }

    std::string extansion(".vtk");
    std::size_t position = fIn.find(extansion);
    fOut = fIn.substr(0,position) + "_INDEXED.vtk";
    std::cout << "INPUT FILE: " << fIn << std::endl;
    std::cout << "OUTPUT FILE: " << fOut << std::endl;


    IGMeshIOService ioService(&m);
    VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N | gmds::R);
    vtkReader.read(fIn);

    MeshDoctor doctor(&m);
    doctor.buildFacesAndR2F();
    doctor.buildEdgesAndX2E();
    doctor.updateUpwardConnectivity();

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
        bnd_normals[n.id()] = math::Vector3d(nv.X(), nv.Y(), nv.Z());
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

    IGMeshIOService ioService2(&m);
    VTKWriter writer2(&ioService2);
    writer2.setCellOptions(gmds::N|gmds::R);
    writer2.setDataOptions(gmds::N|gmds::R);
    writer2.write(fOut);
    return 0;

}

/*----------------------------------------------------------------------------*/
