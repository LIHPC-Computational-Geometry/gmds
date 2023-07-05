/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/cad/GeomPoint.h>
#include <gmds/cad/GeomCurve.h>
#include <gmds/cad/GeomSurface.h>
#include <gmds/smoothy/LaplacianSmoother.h>
#include <gmds/igalgo/BoundaryOperator.h>
/*----------------------------------------------------------------------------*/
#include "gmds/igalgo/BoundaryOperator2D.h"
#include <gmds/smoothy/EllipticSmoother2D.h>
#include <iostream>

/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    std::cout << "============== Morph 2D ================" << std::endl;

    //==================================================================
    // PARAMETERS' PARSING
    //==================================================================
    std::string file_input;
    if (argc != 2) {
        std::cout << "Require two parameters : \n";
        std::cout << "  - [IN ] input mesh (.vtk) that we must deform, \n";
        exit(0);
    }

    file_input = std::string(argv[1]);
    std::cout << "Parameters " << std::endl;
    std::cout << "  - Input file : " << file_input << std::endl;
    std::cout << "=======================================" << std::endl;

    //==================================================================
    // GEOMETRY READING
    //==================================================================
    std::cout<<"> Start geometry reading"<<std::endl;
    Mesh mesh_in(MeshModel(DIM3 | F | E | N |
                            F2N | F2E | E2F | E2N | N2E | N2F));


	 IGMeshIOService ioService(&mesh_in);
    VTKReader vtkReader(&ioService);
	 vtkReader.setCellOptions(N| F | E);
	// vtkReader.setDataOptions(N| F | E);
    vtkReader.read(file_input);
	 std::cout<<"Stats: ("<<mesh_in.getNbFaces()<<", "<<mesh_in.getNbEdges()<<")"<<std::endl;

	 MeshDoctor doc(&mesh_in);
	 doc.buildEdgesAndX2E();
	 doc.updateUpwardConnectivity();
	 doc.orient2DFaces();
	 for(auto f_id:mesh_in.faces()) {
		 Face f=mesh_in.get<Face>(f_id);
		 if (f.normal().dot(math::Vector3d({.0, .0, 1.0})) <= 0) {
			 std::vector<TCellID> ns = f.getIDs<Node>();
			 std::vector<TCellID> ns2(ns.size());
			 for (auto i = 0; i < ns.size(); i++)
				 ns2[ns.size() - 1 - i] = ns[i];
			 f.set<Node>(ns2);
		 }
	 }

	 //==================================================================
	 // PERFORM THE PERTURBATION
	 //==================================================================
	 auto mark_lock= mesh_in.newMark<Node>();

	 math::Point center({0,0,0});
	 double alpha = 2;

	 //==================================================================
	 // MARK ALL THE BOUNDARY CELL OF THE INIT MESH
	 //==================================================================
	 // we get all the nodes that are on the mesh boundary
	 BoundaryOperator2D op(&mesh_in);
	 auto mark_node_NAN = mesh_in.newMark<Node>();
	 auto mark_node_on_pnt = mesh_in.newMark<Node>();
	 auto mark_node_on_crv = mesh_in.newMark<Node>();
	 auto mark_edge_on_crv = mesh_in.newMark<Edge>();

	 op.markCellOnGeometry(mark_edge_on_crv, mark_node_on_crv, mark_node_on_pnt, mark_node_NAN);
	 auto nb_locked = 0;
	 for (auto n_id : mesh_in.nodes()) {
		 if (mesh_in.isMarked<Node>(n_id, mark_node_on_crv) || mesh_in.isMarked<Node>(n_id, mark_node_on_pnt)) {
			 nb_locked += 1;
			 mesh_in.mark<Node>(n_id, mark_lock);
			 Node n = mesh_in.get<Node>(n_id);
			 math::Point p= n.point();
			 n.setPoint(p+alpha*(p-center));
		 }
	 }
	/* std::set<TCellID> corners;
	 Variable<int>* var_lock = mesh_in.newVariable<int, gmds::GMDS_NODE>("lock");
	 for(auto e_id:mesh_in.edges()){
		 Edge e = mesh_in.get<Edge>(e_id);
		 auto e_nodes = e.getIDs<Node>();
		 corners.insert(e_nodes[0]);
		 corners.insert(e_nodes[1]);
	 }
	 for(auto c:corners){
		 Node n = mesh_in.get<Node>(c);
		 math::Point p= n.point();
		 n.setPoint(p+alpha*(p-center));
		 mesh_in.mark(n,mark_lock);
		// var_lock->set(c,1);
	 }*/
	 std::string file_output ="morph2d_perturb.vtk";
	 std::cout<<"> Write perturbed mesh in: "<<file_output<<std::endl;
	 VTKWriter vtkWriter(&ioService);
	 vtkWriter.setCellOptions(N|F);
	 vtkWriter.write(file_output);

	 //==================================================================
	 // PERFORM THE MESH SMOOTHING NOWEllipticSmoothingTestSuite.grid2D_smoother (1)
	 //==================================================================

	 smoothy::EllipticSmoother2D smoother2D(&mesh_in);
	 smoother2D.lock(mark_lock);
	 //smoother2D.setTheta(0);
	 smoother2D.execute();


	 file_output ="morph2d_result.vtk";
    std::cout<<"> Write output mesh in: "<<file_output<<std::endl;

    vtkWriter.setCellOptions(N|F);
    vtkWriter.write(file_output);
    std::cout << "======== Task done =========" << std::endl;
}