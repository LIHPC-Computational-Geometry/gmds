/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/smoothy/EllipticSmoothing.h>
#include <gmds/igalgo/GridBuilder.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <random>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::smoothy;
/*----------------------------------------------------------------------------*/
int
main(int argc, char *argv[])
{
	std::cout << "============== Invertible Mapper 2D ================" << std::endl;

	//==================================================================
	// PARAMETERS' PARSING
	//==================================================================
	std::string input_file, output_file;
	int nb_iterations = 0;
	if (argc != 3) {
		std::cout << "Require four paramaters : \n";
		std::cout << "  - [IN ] triangular mesh (.vtk) that gives the input mes, \n";
		std::cout << "  - [OUT] the smoothed mesh (.vtk). \n" << std::endl;
		throw gmds::GMDSException("Wrong number of parameters");
	}

	input_file = std::string(argv[1]);
	output_file = std::string(argv[2]);
	std::cout << "Parameters " << std::endl;
	std::cout << "  - Input  : " << input_file << std::endl;
	std::cout << "  - Output : " << output_file << std::endl;
	std::cout << "===================================" << std::endl;

	//==================================================================
	// INPUT READING
	//==================================================================
	std::cout << "> Start input reading" << std::endl;
	Mesh m(MeshModel(DIM3 | F | E | N | F2N | F2E | E2F | E2N | N2E | N2F));
	IGMeshIOService ioService(&m);
	/*VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(N | F);
	vtkReader.read(input_file);

*/
	GridBuilder gb(&m,2);
	gb.execute(5,1,5,1);

	MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();


	for(auto f_id:m.faces()){
		Face f = m.get<Face>(f_id);
		std::cout<<"Face "<<f_id<<f.normal()<<std::endl;
		if(f.normal().dot(math::Vector3d({.0,.0,1.0}))<=0){
			std::vector<TCellID> ns = f.getIDs<Node>();
			std::vector<TCellID> ns2(ns.size());
			for(auto i=0;i<ns.size();i++)
				ns2[ns.size()-1-i]=ns[i];
			f.set<Node>(ns2);
		}
	}
	//==================================================================
	// MARK ALL THE BOUNDARY CELL OF THE INIT MESH
	//==================================================================
	std::cout << "> Start mesh boundary retrieval" << std::endl;
	// we get all the nodes that are on the mesh boundary
	BoundaryOperator2D op(&m);
	auto mark_node_NAN = m.newMark<Node>();
	auto mark_node_on_pnt = m.newMark<Node>();
	auto mark_node_on_crv = m.newMark<Node>();
	auto mark_edge_on_crv = m.newMark<Edge>();

	op.markCellOnGeometry(mark_edge_on_crv, mark_node_on_crv, mark_node_on_pnt, mark_node_NAN);
	auto nb_locked = 0;
	auto mark_bnd_nodes = m.newMark<Node>();
	for (auto n_id : m.nodes()) {
		if (m.isMarked<Node>(n_id, mark_node_on_crv) || m.isMarked<Node>(n_id, mark_node_on_pnt)) {
			nb_locked += 1;
			m.mark<Node>(n_id, mark_bnd_nodes);
		}
	}
	std::cout<<"> Nb locked nodes: "<<nb_locked<<std::endl;
	//==================================================================
	// PERFORM THE PERTURBATION
	//==================================================================
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N | F);
	vtkWriter.setDataOptions(N | F);
	vtkWriter.write("init.vtk");
	std::cout << "> Start perturbation" << std::endl;
	constexpr int FLOAT_MIN = -100;
	constexpr int FLOAT_MAX = 108;
	std::random_device rd;
	std::default_random_engine eng(rd());
	std::uniform_real_distribution<float> distr(FLOAT_MIN, FLOAT_MAX);

	for(auto n_id:m.nodes()){
		if(!m.isMarked<Node>(n_id,mark_bnd_nodes)) {
			Node n = m.get<Node>(n_id);
			n.setXYZ(distr(eng), distr(eng), 0);
		}
	}
	vtkWriter.write("perturb.vtk");
	//==================================================================
	// PERFORM THE MESH SMOOTHING NOW
	//==================================================================
	std::cout << "> Start smoothing" << std::endl;


	std::vector<double> node_coords(2*m.getNbNodes());
	std::vector<bool> lock_nodes(2*m.getNbNodes());
	for(auto n_id:m.nodes()){
		math::Point pi = m.get<Node>(n_id).point();
		node_coords[2*n_id]=pi.X();
		node_coords[2*n_id+1]=pi.Y();
		lock_nodes[2*n_id]= (m.isMarked<Node>(n_id, mark_bnd_nodes));
		lock_nodes[2*n_id+1]= (m.isMarked<Node>(n_id, mark_bnd_nodes));
	}

	std::vector<std::array<int, 3>> triangles(4*m.getNbFaces());
	constexpr int quadrature[4][3]= {{0,1,3},{1,2,0},{2,3,1},{3,0,2}};
	for(auto f_id:m.faces()){
		std::vector<TCellID> nids = m.get<Face>(f_id).getIDs<Node>();
		for(auto i=0;i<4;i++) {
			triangles[4 * f_id + i] = {static_cast<int>(nids[quadrature[i][0]]),
			                           static_cast<int>(nids[quadrature[i][1]]),
			                           static_cast<int>(nids[quadrature[i][2]])};
		}
	}
	EllipticSmoothingMeshGlue2D tri(node_coords,triangles,lock_nodes);

	std::cout<<"==================="<<std::endl;
	for(auto c:node_coords)
		std::cout<<c<<" ";
	std::cout<<std::endl;
	std::cout<<"==================="<<std::endl;
	EllipticSmoothingOptions options = ESO_default2D;
	options.maxiter=10;//00;
	options.eps_from_theorem=true;
	options.theta = 1e-3;
	options.bfgs_threshold = 1e-9;
	options.debug=0;
	EllipticSmoothing2D smoother2D(tri, triangles.size(),options);
//	smoother2D.start_eps = 1E-4;
	bool res = smoother2D.execute();
	std::cout<<"Resultat: "<<res<<std::endl;
	std::cout<<"==================="<<std::endl;
	tri.get_verts(node_coords);
	for(auto c:node_coords)
		std::cout<<c<<" ";
	std::cout<<std::endl;
	std::cout<<"==================="<<std::endl;
for(auto n_id:m.nodes()){
	m.get<Node>(n_id).setXYZ(node_coords[2*n_id],node_coords[2*n_id+1],.0);
}
	std::cout << "> Write output mesh in: " << output_file << std::endl;
	vtkWriter.write(output_file);
	std::cout << "======== Task done by Smoothy =========" << std::endl;
}