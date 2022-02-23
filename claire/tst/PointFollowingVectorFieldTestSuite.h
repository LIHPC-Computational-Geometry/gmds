//
// Created by rochec on 26/01/2022.
//

#include <gmds/claire/LevelSetCombined.h>
#include <gmds/claire/LevelSetEloi.h>
#include <gmds/claire/LevelSetExtended.h>
#include <gmds/claire/GradientComputation2D.h>
#include <gmds/claire/LeastSquaresGradientComputation.h>
#include <gmds/claire/AdvectedPointRK4_2D.h>
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


/*----------------------------------------------------------------------------*/
/*                 CAS TEST 2D CLASSE AdvectedPointRK4_2D                     */
/*----------------------------------------------------------------------------*/

TEST(PointFollowingVectorFieldTestClass, AdvectedPointRK4_2D_Test1)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
	                             gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/Carre.vtk";
	//std::string vtk_file = dir+"/Carre_maxsize_0.01.vtk";		// Same geo but refined mesh

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
	doc.orient2DFaces();

	//Get the boundary node ids
	BoundaryOperator2D bnd_op(&m);
	std::vector<TCellID> bnd_node_ids;
	bnd_op.getBoundaryNodes(bnd_node_ids);

	// Initialisation des marques sur le front à avancer
	int markFrontNodes = m.newMark<gmds::Node>();
	for(auto id:bnd_node_ids){
		Node n = m.get<Node>(id);
		double coord_y = n.Y() ;
		if (coord_y == 0) {
			// For the square test case, the front to advance is the boundary where y=0
			m.mark<Node>(id,markFrontNodes);
		}
	}

	// Calcul des Level Set
	Variable<double> *var_dist = m.newVariable<double,GMDS_NODE>("GMDS_Distance");
	LevelSetExtended ls(&m, markFrontNodes, var_dist);
	LevelSetExtended::STATUS result_ls = ls.execute();
	ASSERT_EQ(LevelSetExtended::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	// Calcul du gradient du champ de Level Set
	m.newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");
	LeastSquaresGradientComputation grad2D(&m, var_dist, m.getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	LeastSquaresGradientComputation::STATUS result_grad = grad2D.execute();
	ASSERT_EQ(LeastSquaresGradientComputation::SUCCESS, result_grad);

	// Placement du point P à la distance souhaitée suivant le champ de gradient
	math::Point M(0.5, 0.0, 0.0);
	double distance = 0.3;
	AdvectedPointRK4_2D advpoint(&m, M, distance, var_dist, m.getVariable<math::Vector3d ,GMDS_NODE>("GMDS_Gradient"));
	AdvectedPointRK4_2D::STATUS result = advpoint.execute();

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("AdvectedPointRK4_2D_Test1_Result.vtk");

	ASSERT_EQ(AdvectedPointRK4_2D::SUCCESS, result);
}

TEST(PointFollowingVectorFieldTestClass, AdvectedPointRK4_2D_Test2)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
	                             gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/Carre_Quart_Cylindre_maxsize_0.1.vtk";
	//std::string vtk_file = dir+"/Carre_Quart_Cylindre_maxsize_0.01.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	doc.orient2DFaces();

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
		if ( ( sqrt( pow(coord_x,2) + pow(coord_y,2)) - 1 ) <= pow(10,-6)) {
			// For this test case, the front to advance is the boundary where x²+y²=1
			m.mark<Node>(id,markFrontNodesInt);
		}
		else if (coord_x == -2 || coord_y == 2) {
			m.mark<Node>(id,markFrontNodesOut);
		}
	}

	m.newVariable<double,GMDS_NODE>("GMDS_Distance");
	m.newVariable<double,GMDS_NODE>("GMDS_Distance_Int");
	m.newVariable<double,GMDS_NODE>("GMDS_Distance_Out");
	LevelSetCombined lsCombined(&m, markFrontNodesInt, markFrontNodesOut,
	                            m.getVariable<double,GMDS_NODE>("GMDS_Distance"),
	                            m.getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
	                            m.getVariable<double,GMDS_NODE>("GMDS_Distance_Out"));
	LevelSetCombined::STATUS result_ls = lsCombined.execute();
	ASSERT_EQ(LevelSetCombined::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodesInt);
	m.freeMark<Node>(markFrontNodesInt);
	m.unmarkAll<Node>(markFrontNodesOut);
	m.freeMark<Node>(markFrontNodesOut);

	m.newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");
	LeastSquaresGradientComputation grad2D(&m, m.getVariable<double,GMDS_NODE>("GMDS_Distance"),
	   m.getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	LeastSquaresGradientComputation::STATUS result_grad = grad2D.execute();
	ASSERT_EQ(LeastSquaresGradientComputation::SUCCESS, result_grad);

	// Placement du point P à la distance souhaitée suivant le champ de gradient
	double ang(M_PI/4);
	math::Point M(-cos(ang), sin(ang), 0.0);
	double distance = 1.0;
	AdvectedPointRK4_2D advpoint(&m, M, distance, m.getVariable<double,GMDS_NODE>("GMDS_Distance"),
	   m.getVariable<math::Vector3d ,GMDS_NODE>("GMDS_Gradient"));
	AdvectedPointRK4_2D::STATUS result = advpoint.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("AdvectedPointRK4_2D_Test2_Result.vtk");

	ASSERT_EQ(AdvectedPointRK4_2D::SUCCESS, result);

}

TEST(PointFollowingVectorFieldTestClass, AdvectedPointRK4_2D_Test3)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
	                             gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	//std::string vtk_file = dir+"/Carre_dans_cercle.vtk";
	std::string vtk_file = dir+"/Carre_dans_cercle_maxsize_0.1.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	doc.orient2DFaces();

	// Noeud détruit car n'appartient pas au maillage ni à la géométrie.
	// Il apparaît à cause de la façon dont a été généré le cas test avec GMSH.
	//m.deleteNode(5);

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
		if ( abs( sqrt( pow(coord_x,2) + pow(coord_y,2)) - 4 ) <= pow(10,-6)) {
			// For this test case, the front to advance is the boundary where x²+y²=4
			m.mark<Node>(id,markFrontNodesOut);
		}
		else if ( (coord_x == -1 && coord_y <= 1 ) || (coord_y == 1 && coord_x >= -1 )) {
			m.mark<Node>(id,markFrontNodesInt);
		}
	}

	m.newVariable<double,GMDS_NODE>("GMDS_Distance");
	m.newVariable<double,GMDS_NODE>("GMDS_Distance_Int");
	m.newVariable<double,GMDS_NODE>("GMDS_Distance_Out");
	LevelSetCombined lsCombined(&m, markFrontNodesInt, markFrontNodesOut,
	                            m.getVariable<double,GMDS_NODE>("GMDS_Distance"),
	                            m.getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
	                            m.getVariable<double,GMDS_NODE>("GMDS_Distance_Out"));
	LevelSetCombined::STATUS result_ls = lsCombined.execute();
	ASSERT_EQ(LevelSetCombined::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodesInt);
	m.freeMark<Node>(markFrontNodesInt);
	m.unmarkAll<Node>(markFrontNodesOut);
	m.freeMark<Node>(markFrontNodesOut);

	m.newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");
	LeastSquaresGradientComputation grad2D(&m, m.getVariable<double,GMDS_NODE>("GMDS_Distance"), m.getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	LeastSquaresGradientComputation::STATUS result_grad = grad2D.execute();
	ASSERT_EQ(LeastSquaresGradientComputation::SUCCESS, result_grad);

	// Placement du point P à la distance souhaitée suivant le champ de gradient
	math::Point M(-1,1, 0.0);
	double distance = 1.0;
	AdvectedPointRK4_2D advpoint(&m, M, distance, m.getVariable<double,GMDS_NODE>("GMDS_Distance"), m.getVariable<math::Vector3d ,GMDS_NODE>("GMDS_Gradient"));
	AdvectedPointRK4_2D::STATUS result = advpoint.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("AdvectedPointRK4_2D_Test3_Result.vtk");

	ASSERT_EQ(AdvectedPointRK4_2D::SUCCESS, result);

}

TEST(PointFollowingVectorFieldTestClass, AdvectedPointRK4_2D_Test4)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
	                             gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	//std::string vtk_file = dir+"/Marche.vtk";
	std::string vtk_file = dir+"/Marche_maxsize_0.1.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	doc.orient2DFaces();

	// Noeud détruit car n'appartient pas au maillage ni à la géométrie.
	// Il apparaît à cause de la façon dont a été généré le cas test avec GMSH.
	//m.deleteNode(5);

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
		if ( abs( sqrt( pow(coord_x,2) + pow(coord_y,2)) - 4 ) <= pow(10,-6)) {
			// For this test case, the front to advance is the boundary where x²+y²=4
			m.mark<Node>(id,markFrontNodesOut);
		}
		else if ( (coord_x == -2 && coord_y <= 1 ) || (coord_y == 1 && coord_x >= -2 ) || (coord_x == -1 && coord_y <= 2 ) || (coord_y == 2 && coord_x >= -1 ) ) {
			m.mark<Node>(id,markFrontNodesInt);
		}
	}

	m.newVariable<double,GMDS_NODE>("GMDS_Distance");
	m.newVariable<double,GMDS_NODE>("GMDS_Distance_Int");
	m.newVariable<double,GMDS_NODE>("GMDS_Distance_Out");
	LevelSetCombined lsCombined(&m, markFrontNodesInt, markFrontNodesOut,
	                            m.getVariable<double,GMDS_NODE>("GMDS_Distance"),
	                            m.getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
	                            m.getVariable<double,GMDS_NODE>("GMDS_Distance_Out"));
	LevelSetCombined::STATUS result_ls = lsCombined.execute();
	ASSERT_EQ(LevelSetCombined::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodesInt);
	m.freeMark<Node>(markFrontNodesInt);
	m.unmarkAll<Node>(markFrontNodesOut);
	m.freeMark<Node>(markFrontNodesOut);

	m.newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");
	LeastSquaresGradientComputation grad2D(&m, m.getVariable<double,GMDS_NODE>("GMDS_Distance"),
	   m.getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	LeastSquaresGradientComputation::STATUS result_grad = grad2D.execute();
	ASSERT_EQ(LeastSquaresGradientComputation::SUCCESS, result_grad);

	// Placement du point P à la distance souhaitée suivant le champ de gradient
	math::Point M(-1,1, 0.0);
	double distance = 1.0;
	AdvectedPointRK4_2D advpoint(&m, M, distance, m.getVariable<double,GMDS_NODE>("GMDS_Distance"),
	   m.getVariable<math::Vector3d ,GMDS_NODE>("GMDS_Gradient"));
	AdvectedPointRK4_2D::STATUS result = advpoint.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("AdvectedPointRK4_2D_Test4_Result.vtk");

	ASSERT_EQ(AdvectedPointRK4_2D::SUCCESS, result);

}

/*----------------------------------------------------------------------------*/