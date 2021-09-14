/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>

#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
/*----------------------------------------------------------------------------*/
#include <gmds/singGraphBuild/SingGraphBuilder2DOriginal.h>
#include <gmds/singGraphBuild/SingGraphBuilder2DShortestPath.h>
#include <gmds/singGraphBuild/SingGraphBuilder2DSimultStartHeun.h>
#include <gmds/singGraphBuild/SingGraphBuilder2DSimultStartRK4.h>
#include <gmds/singGraphBuild/SingularityGraphBuilder2D.h>
/*----------------------------------------------------------------------------*/
#include <chrono>
#include <iostream>
#include <math.h>
#include <vector>
#ifndef _WIN32
#	include <unistd.h>
#endif
typedef std::chrono::high_resolution_clock Clock;
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

TEST(SingGraphBuilder2DTest, halfDisk)
{
	auto t0 = Clock::now();

	MeshModel model(DIM3 | F | E | N | F2N | E2N | F2E | E2F | N2F | N2E | F2F);
	Mesh mesh(model);

	//==================================================================
	// MESH READING
	//==================================================================
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/singGraphBuilderTest.vtk";

	gmds::IGMeshIOService ioService(&mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::F);
	vtkReader.setDataOptions(gmds::N | gmds::F);
	vtkReader.read(vtk_file);

	//==================================================================
	// MESH TOPOLOGY PREPARATION
	//==================================================================
	MeshDoctor doctor(&mesh);
	doctor.buildEdgesAndX2E();
	doctor.updateUpwardConnectivity();
	doctor.orient2DFaces();

	EXPECT_TRUE(mesh.getNbFaces() != 0);
	EXPECT_TRUE(mesh.getNbEdges() != 0);
	EXPECT_TRUE(mesh.getNbNodes() != 0);

	//==================================================================
	// INIT MARKs FOR BOUNDARY NODES
	//==================================================================
	int edge_curve_mark = mesh.newMark<gmds::Edge>();
	int node_curve_mark = mesh.newMark<gmds::Node>();
	int node_point_mark = mesh.newMark<gmds::Node>();
	int node_forbiddenBdry_mark = mesh.newMark<gmds::Node>();

	BoundaryOperator boundaryOp(&mesh);
	boundaryOp.setSurfaceAngleDot(0.86);     // pi/6
	boundaryOp.markCellsOnCurves(edge_curve_mark, node_curve_mark);
	boundaryOp.markNodesOnPoint(edge_curve_mark, node_curve_mark, node_point_mark);

	//==================================================================
	// CROSS FIELD EXTRACTION FROM THE IMPORTED FILE
	//==================================================================
	Variable<math::Vector3d> *field_X = mesh.getVariable<math::Vector3d, gmds::GMDS_NODE>("cross_X");
	Variable<math::Cross2D> *field = mesh.newVariable<math::Cross2D, gmds::GMDS_NODE>("c");

	math::Vector3d OX(1, 0, 0);
	for (auto n_id : mesh.nodes()) {
		Node n = mesh.get<Node>(n_id);
		math::Vector3d vx = (*field_X)[n.id()];
		(*field)[n.id()] = math::Cross2D(4 * vx.angle(OX));
	}

	//==================================================================
	// SINGULARITY GRAPH EXTRACTION
	//==================================================================
	const bool BuildGeomSing = true;

	SingularityGraphBuilder2D::Strategy strat = SingularityGraphBuilder2D::shortestPaths;
	auto algo = SingGraphBuilder2DShortestPath(&mesh, field, BuildGeomSing);
	algo.setDebugPrefix("singGraphBuilderTest");
	algo.initMarks(node_point_mark, node_curve_mark, edge_curve_mark, node_forbiddenBdry_mark);
	algo.setQuadMeshSmoothingEnable(true);
	algo.setDebugFilesWritingEnable(true);

	auto t1 = Clock::now();
	std::cout << "Before execute() " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " milliseconds" << std::endl;
	algo.execute();
	auto t2 = Clock::now();
	std::cout << "execute ftc " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " milliseconds" << std::endl;
	std::cout << "TOTAL   " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t0).count() << " milliseconds" << std::endl;

	//==================================================================
	// MARKS CLEANING
	//==================================================================
	mesh.unmarkAll<Node>(node_curve_mark);
	mesh.unmarkAll<Node>(node_point_mark);
	mesh.unmarkAll<Edge>(edge_curve_mark);
	mesh.unmarkAll<Node>(node_forbiddenBdry_mark);
	mesh.freeMark<Node>(node_curve_mark);
	mesh.freeMark<Node>(node_point_mark);
	mesh.freeMark<Edge>(edge_curve_mark);
	mesh.freeMark<Node>(node_forbiddenBdry_mark);

	const ::SingularityGraph &g = algo.graph();
	ASSERT_EQ(g.getNbPoints(), 80);
	ASSERT_EQ(g.getNbLines(), 141);
	ASSERT_EQ(g.getNbSurfacePatches(), 58);
	ASSERT_EQ(g.getNLinkedEdgeGroup(), 26);
}
