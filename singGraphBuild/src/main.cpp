/*----------------------------------------------------------------------------*/
#include <chrono>
#include <iostream>
#include <string>
/*----------------------------------------------------------------------------*/
#include <gmds/singGraphBuild/SingGraphBuilder2DOriginal.h>
#include <gmds/singGraphBuild/SingGraphBuilder2DShortestPath.h>
#include <gmds/singGraphBuild/SingGraphBuilder2DSimultStartHeun.h>
#include <gmds/singGraphBuild/SingGraphBuilder2DSimultStartRK4.h>
#include <gmds/singGraphBuild/SingularityGraphBuilder2D.h>
/*----------------------------------------------------------------------------*/
#include <glpk.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
int
main(int argc, char *argv[])
{
	std::cout << "==== Singularity Graph Builder ====" << std::endl;
	MeshModel model(DIM3 | F | E | N | F2N | E2N | F2E | E2F | N2F | N2E | F2F);
	Mesh mesh(model);

	//==================================================================
	// MESH READING
	//==================================================================
	std::cout << "Reading " << std::endl;
	try {
		typedef std::chrono::high_resolution_clock Clock;

		auto t0 = Clock::now();

		if (argc != 4) throw gmds::GMDSException("Wrong parameters");

		std::string fIn, outDir;

		fIn = std::string(argv[1]);
		outDir = std::string(argv[2]);

		int method_index = std::stoi(std::string(argv[3]));
		SingularityGraphBuilder2D::Strategy strategy = SingularityGraphBuilder2D::original;
		std::cout << "Method index: " << method_index << std::endl;
		if (method_index == 0) {
			strategy = SingularityGraphBuilder2D::shortestPaths;
		}
		else if (method_index == 1) {
			strategy = SingularityGraphBuilder2D::simultaneousStartHeun;
		}
		else if (method_index == 2) {
			strategy = SingularityGraphBuilder2D::simultaneousStartRK4;
		}
		std::cout << "Start reading file " << fIn << std::endl;
		gmds::IGMeshIOService ioService(&mesh);
		gmds::VTKReader vtkReader(&ioService);
		vtkReader.setCellOptions(gmds::N | gmds::F);
		vtkReader.setDataOptions(gmds::N | gmds::F);
		vtkReader.read(fIn);

		//==================================================================
		// MESH TOPOLOGY PREPARATION
		//==================================================================
		gmds::MeshDoctor doctor(&mesh);
		doctor.buildEdgesAndX2E();
		doctor.updateUpwardConnectivity();
		doctor.orient2DFaces();

		//==================================================================
		// INIT MARKs FOR BOUNDARY NODES
		//==================================================================
		int edge_curve_mark = mesh.newMark<gmds::Edge>();
		int node_curve_mark = mesh.newMark<gmds::Node>();
		int node_point_mark = mesh.newMark<gmds::Node>();

		BoundaryOperator boundaryOp(&mesh);
		boundaryOp.markCellsOnCurves(edge_curve_mark, node_curve_mark);
		boundaryOp.markNodesOnPoint(edge_curve_mark, node_curve_mark, node_point_mark);

		//==================================================================
		// CROSS FIELD EXTRACTION FROM THE IMPORTED FILE
		//==================================================================

		Variable<math::Vector3d> *field_X = mesh.getVariable<math::Vector3d, GMDS_NODE>("cross_X");
		Variable<math::Cross2D> *field = mesh.newVariable<math::Cross2D, GMDS_NODE>("c");

		math::Vector3d OX(1, 0, 0);
		for (auto n_id : mesh.nodes()) {
			Node n = mesh.get<Node>(n_id);
			math::Vector3d vx = (*field_X)[n.id()];
			(*field)[n.id()] = math::Cross2D(4 * vx.angle(OX));
		}

		auto t1 = Clock::now();
		std::cout << "topology preparation " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " milliseconds" << std::endl;
		//==================================================================
		// SINGULARITY GRAPH EXTRACTION
		//==================================================================
		SingularityGraphBuilder2D *algo;
		const bool BuildGeomSing = true;
		// TODO : create factory method
		switch (strategy) {
		case SingularityGraphBuilder2D::original: algo = new SingGraphBuilder2DOriginal(&mesh, field, BuildGeomSing); break;
		case SingularityGraphBuilder2D::simultaneousStartHeun: algo = new SingGraphBuilder2DSimultStartHeun(&mesh, field, BuildGeomSing); break;
		case SingularityGraphBuilder2D::simultaneousStartRK4: algo = new SingGraphBuilder2DSimultStartRK4(&mesh, field, BuildGeomSing); break;
		case SingularityGraphBuilder2D::testPrescribedSing:
		case SingularityGraphBuilder2D::shortestPaths: algo = new SingGraphBuilder2DShortestPath(&mesh, field, BuildGeomSing); break;
		default: throw GMDSException("Wrong strategy value for singularity graph computation");
		}

		algo->setDebugPrefix(outDir + "/sgb");
		algo->initMarks(node_point_mark, node_curve_mark, edge_curve_mark);
		// User must choose the number of control Points for the Bezier curve computation (last step)
		// WARNING if the number chosen exceeds the original number of points for the curve, the latter will be chosen

		unsigned int number_of_control_points = 8;
		algo->execute(strategy, number_of_control_points);

		//==================================================================
		// MARKS CLEANING
		//==================================================================
		mesh.unmarkAll<Node>(node_curve_mark);
		mesh.unmarkAll<Node>(node_point_mark);
		mesh.unmarkAll<Edge>(edge_curve_mark);
		mesh.freeMark<Node>(node_curve_mark);
		mesh.freeMark<Node>(node_point_mark);
		mesh.freeMark<Edge>(edge_curve_mark);
	}
	catch (exception &e) {
		std::cout << e.what();
		std::cout << std::endl;
	}
	return 0;
}
/*----------------------------------------------------------------------------*/
