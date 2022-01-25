//
// Created by rochec on 21/01/2022.
//

#include <gmds/claire/LevelSet.h>
#include <gmds/claire/LevelSetFromIntToOut.h>
#include <gmds/claire/GradientComputation2D.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gtest/gtest.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

TEST(GradientComputation2DTestClass, GradientComputation2D_Test1)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
	                             gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/Carre.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	//Get the boundary node ids
	BoundaryOperator2D bnd_op(&m);
	std::vector<TCellID> bnd_node_ids;
	bnd_op.getBoundaryNodes(bnd_node_ids);

	int markFrontNodes = m.newMark<gmds::Node>();
	for(auto id:bnd_node_ids){
		Node n = m.get<Node>(id);
		double coord_y = n.Y() ;
		if (coord_y == 0) {
			// For the square test case, the front to advance is the boundary where y=0
			m.mark<Node>(id,markFrontNodes);
		}
	}

	LevelSetNaif ls(&m, markFrontNodes);
	ls.execute();

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	GradientComputation2D grad2D(&m, m.getVariable<double,GMDS_NODE>("distance"));
	grad2D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("GradientComputation2D_Test1_Result.vtk");

	ASSERT_TRUE(true);
}


TEST(GradientComputation2DTestClass, GradientComputation2D_Test2)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
	                             gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/Carre_Quart_Cylindre_maxsize_0.1.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	// Noeud détruit car n'appartient pas au maillage ni à la géométrie.
	// Il apparaît à cause de la façon dont a été généré le cas test avec GMSH.
	m.deleteNode(5);

	//Get the boundary node ids
	BoundaryOperator2D bnd_op(&m);
	std::vector<TCellID> bnd_node_ids;
	bnd_op.getBoundaryNodes(bnd_node_ids);

	int markFrontNodes = m.newMark<gmds::Node>();
	for(auto id:bnd_node_ids){
		Node n = m.get<Node>(id);
		double coord_y = n.Y() ;
		double coord_x = n.X() ;
		if ( (sqrt( pow(coord_y,2) + pow(coord_x,2)) - 1) <= pow(10,-6) ) {
			// For this test case, the front to advance is the boundary where x²+y²=1
			m.mark<Node>(id,markFrontNodes);
			//std::cout << "Noeud marqué :" << id << std::endl;
		}
	}

	LevelSetNaif ls(&m, markFrontNodes);
	ls.execute();

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	GradientComputation2D grad2D(&m, m.getVariable<double,GMDS_NODE>("distance"));
	grad2D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("GradientComputation2D_Test2_Result.vtk");

	ASSERT_TRUE(true);
}


TEST(GradientComputation2DTestClass, GradientComputation2D_Test3)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
	                             gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	//std::string vtk_file = dir+"/Carre_Quart_Cylindre_maxsize_0.1.vtk";
	std::string vtk_file = dir+"/Carre_Quart_Cylindre.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	// Noeud détruit car n'appartient pas au maillage ni à la géométrie.
	// Il apparaît à cause de la façon dont a été généré le cas test avec GMSH.
	m.deleteNode(5);

	//Get the boundary node ids
	BoundaryOperator2D bnd_op(&m);
	std::vector<TCellID> bnd_node_ids;
	bnd_op.getBoundaryNodes(bnd_node_ids);

	int markFrontNodesInt = m.newMark<gmds::Node>();
	int markFrontNodesOut = m.newMark<gmds::Node>();
	for(auto id:bnd_node_ids){
		Node n = m.get<Node>(id);
		double coord_y = n.Y() ;
		double coord_x = n.X() ;
		if ( ( sqrt( pow(coord_x,2) + pow(coord_y,2)) - 1.0 ) <= pow(10,-6)) {
			// For this test case, the front to advance is the boundary where x²+y²=1
			m.mark<Node>(id,markFrontNodesInt);
		}
		else if (coord_x == -2 || coord_y == 2) {
			m.mark<Node>(id,markFrontNodesOut);
		}
	}

	LevelSetFromIntToOut lsCombined(&m, markFrontNodesInt, markFrontNodesOut);
	lsCombined.execute();

	m.unmarkAll<Node>(markFrontNodesInt);
	m.freeMark<Node>(markFrontNodesInt);
	m.unmarkAll<Node>(markFrontNodesOut);
	m.freeMark<Node>(markFrontNodesOut);

	GradientComputation2D grad2D(&m, m.getVariable<double,GMDS_NODE>("distance_combined"));
	grad2D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("GradientComputation2D_Test3_Result.vtk");

	ASSERT_TRUE(true);
}
