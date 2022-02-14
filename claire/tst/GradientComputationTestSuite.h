//
// Created by rochec on 21/01/2022.
//

#include <gmds/claire/LevelSetExtended.h>
#include <gmds/claire/LevelSetEloi.h>
#include <gmds/claire/LevelSetCombined.h>
#include <gmds/claire/GradientComputation2D.h>
#include <gmds/claire/GradientComputation3D.h>
#include <gmds/claire/LeastSquaresGradientComputation.h>
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
/*                  CAS TEST CLASSE GradientComputation2D                     */
/*----------------------------------------------------------------------------*/

TEST(GradientComputationTestClass, GradientComputation2D_Test1)
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
	doc.orient2DFaces();

	math::Point P;

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

	Variable<double> *var_dist = m.newVariable<double,GMDS_NODE>("GMDS_Distance");
	LevelSetEloi ls(&m, markFrontNodes, var_dist);
	LevelSetEloi::STATUS result_ls = ls.execute();

	ASSERT_EQ(LevelSetEloi::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	GradientComputation2D grad2D(&m, m.getVariable<double,GMDS_NODE>("GMDS_Distance"));
	GradientComputation2D::STATUS result = grad2D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("GradientComputation2D_Test1_Result.vtk");

	ASSERT_EQ(GradientComputation2D::SUCCESS, result);
}

TEST(GradientComputationTestClass, GradientComputation2D_Test2)
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

	doc.orient2DFaces();

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

	Variable<double> *var_dist = m.newVariable<double,GMDS_NODE>("GMDS_Distance");
	LevelSetEloi ls(&m, markFrontNodes, var_dist);
	LevelSetEloi::STATUS result_ls = ls.execute();

	ASSERT_EQ(LevelSetEloi::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	GradientComputation2D grad2D(&m, m.getVariable<double,GMDS_NODE>("GMDS_Distance"));
	GradientComputation2D::STATUS result = grad2D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("GradientComputation2D_Test2_Result.vtk");

	ASSERT_EQ(GradientComputation2D::SUCCESS, result);
}

TEST(GradientComputationTestClass, GradientComputation2D_Test3)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
	                             gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/Carre_Quart_Cylindre_maxsize_0.1.vtk";
	//std::string vtk_file = dir+"/Carre_Quart_Cylindre.vtk";

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
		if ( ( sqrt( pow(coord_x,2) + pow(coord_y,2)) - 1.0 ) <= pow(10,-6)) {
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

	GradientComputation2D grad2D(&m, m.getVariable<double,GMDS_NODE>("GMDS_Distance"));
	GradientComputation2D::STATUS result = grad2D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("GradientComputation2D_Test3_Result.vtk");

	ASSERT_EQ(GradientComputation2D::SUCCESS, result);
}

TEST(GradientComputationTestClass, GradientComputation2D_Test4)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
	                             gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/Pont.vtk";

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

	int markFrontNodes = m.newMark<gmds::Node>();
	for(auto id:bnd_node_ids){
		Node n = m.get<Node>(id);
		double coord_y = n.Y() ;
		double coord_x = n.X() ;
		if ( coord_x == 0 && coord_y == 0 ) {
			// For this test case, the front to advance is the boundary where x²+y²=1
			m.mark<Node>(id,markFrontNodes);
		}
	}

	Variable<double> *var_dist = m.newVariable<double,GMDS_NODE>("GMDS_Distance");
	LevelSetEloi ls(&m, markFrontNodes, var_dist);
	LevelSetEloi::STATUS result_ls = ls.execute();

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	GradientComputation2D grad2D(&m, m.getVariable<double,GMDS_NODE>("GMDS_Distance"));
	GradientComputation2D::STATUS result = grad2D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("GradientComputation2D_Test4_Result.vtk");

	ASSERT_EQ(GradientComputation2D::SUCCESS, result);
}

/*----------------------------------------------------------------------------*/




/*----------------------------------------------------------------------------*/
/*                  CAS TEST CLASSE GradientComputation3D                     */
/*----------------------------------------------------------------------------*/

TEST(GradientComputationTestClass, GradientComputation3D_Test1)
{
	Mesh m(MeshModel(DIM3 | R | F | E | N |
	                 R2N | F2N | E2N | R2F | F2R |
	                 F2E | E2F | R2E | N2R | N2F | N2E));
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/Cube.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::R);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doctor(&m);
	doctor.buildFacesAndR2F();
	doctor.buildEdgesAndX2E();
	doctor.updateUpwardConnectivity();

	int markFrontNodes = m.newMark<gmds::Node>();

	// Initialisation de la marque pour noter quels fronts sont à avancer
	for(auto id:m.nodes()){
		Node n = m.get<Node>(id);
		double coord_x = n.X() ;
		if ( abs(coord_x) <= pow(10,-6) ) {
			// For this test case, the front to advance is the boundary where x=0
			m.mark<Node>(id,markFrontNodes);
		}
	}

	std::cout << "Fin de l'initialisation des marques" << std::endl ;

	Variable<double> *var_dist = m.newVariable<double,GMDS_NODE>("GMDS_Distance");
	LevelSetEloi ls(&m, markFrontNodes, var_dist);
	LevelSetEloi::STATUS result_ls = ls.execute();

	ASSERT_EQ(LevelSetEloi::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	GradientComputation3D grad3D(&m, m.getVariable<double,GMDS_NODE>("GMDS_Distance"));
	GradientComputation3D::STATUS result = grad3D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("GradientComputation3D_Test1_Result.vtk");

	ASSERT_EQ(GradientComputation3D::SUCCESS, result);
}

TEST(GradientComputationTestClass, GradientComputation3D_Test2)
{
	Mesh m(MeshModel(DIM3 | R | F | E | N |
	                 R2N | F2N | E2N | R2F | F2R |
	                 F2E | E2F | R2E | N2R | N2F | N2E));
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/B0.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::R);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doctor(&m);
	doctor.buildFacesAndR2F();
	doctor.buildEdgesAndX2E();
	doctor.updateUpwardConnectivity();

	// Initialisation des marques sur les bords internes et externes pour noter quels fronts sont à avancer
	int markFrontNodes = m.newMark<gmds::Node>();
	for(auto id:m.nodes()){
		Node n = m.get<Node>(id);
		double coord_y = n.Y() ;
		double coord_x = n.X() ;
		double rayon;
		rayon = sqrt( (pow(coord_x, 2) + pow(coord_y + 2.5, 2)) ) ;
		if ( (rayon - 2.5) < pow(10,-3)) {
			// For this test case, the front to advance is the boundary where x²+y²=2.5
			m.mark<Node>(id,markFrontNodes);
		}
	}


	Variable<double> *var_dist = m.newVariable<double,GMDS_NODE>("GMDS_Distance");
	LevelSetEloi ls(&m, markFrontNodes, var_dist);
	LevelSetEloi::STATUS result_ls = ls.execute();
	ASSERT_EQ(LevelSetEloi::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	GradientComputation3D grad3D(&m, m.getVariable<double,GMDS_NODE>("GMDS_Distance"));
	GradientComputation3D::STATUS result = grad3D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("GradientComputation3D_Test2_Result.vtk");

	ASSERT_EQ(GradientComputation3D::SUCCESS, result);
}

TEST(GradientComputationTestClass, GradientComputation3D_Test3)
{
	/*-----------------------------------------------------------*/
	// Cas test sur la géométrie B0. Les Level Sets sont calculés
	// en deux fois puis combinés de l'intérieur vers l'extérieur.
	/*-----------------------------------------------------------*/
	Mesh m(MeshModel(DIM3 | R | F | E | N |
	                 R2N | F2N | E2N | R2F | F2R |
	                 F2E | E2F | R2E | N2R | N2F | N2E));
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/B0.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::R);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doctor(&m);
	doctor.buildFacesAndR2F();
	doctor.buildEdgesAndX2E();
	doctor.updateUpwardConnectivity();

	int markFrontNodesInt = m.newMark<gmds::Node>();
	int markFrontNodesOut = m.newMark<gmds::Node>();
	for(auto id:m.nodes()){
		Node n = m.get<Node>(id);
		double coord_y = n.Y() ;
		double coord_x = n.X() ;
		double rayon = sqrt( (pow(coord_x, 2) + pow(coord_y + 2.5, 2)) ) ;
		if ( (rayon - 2.5) < pow(10,-6)) {
			// For this test case, the front to advance is the boundary where x²+y²=2.5
			m.mark<Node>(id,markFrontNodesInt);
		}
		else if ( abs(coord_x+5) < pow(10,-6) || abs(coord_x-5) < pow(10,-6) || abs(coord_y-2.5) < pow(10,-6) ) {
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
	LevelSetCombined::STATUS result_lscombined = lsCombined.execute();

	m.unmarkAll<Node>(markFrontNodesInt);
	m.freeMark<Node>(markFrontNodesInt);
	m.unmarkAll<Node>(markFrontNodesOut);
	m.freeMark<Node>(markFrontNodesOut);

	ASSERT_EQ(LevelSetCombined::SUCCESS, result_lscombined);

	GradientComputation3D grad3D(&m, m.getVariable<double,GMDS_NODE>("GMDS_Distance"));
	GradientComputation3D::STATUS result = grad3D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("GradientComputation3D_Test3_Result.vtk");

	ASSERT_EQ(GradientComputation3D::SUCCESS, result);
}

/*----------------------------------------------------------------------------*/




/*----------------------------------------------------------------------------*/
/*             CAS TEST 2D CLASSE LeastSquaresGradientComputation             */
/*----------------------------------------------------------------------------*/

TEST(GradientComputationTestClass, LeastSquaresGradientComputation_2D_Test1)
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

	Variable<double> *var_dist = m.newVariable<double,GMDS_NODE>("GMDS_Distance");
	LevelSetEloi ls(&m, markFrontNodes, var_dist);
	LevelSetEloi::STATUS result_ls = ls.execute();
	ASSERT_EQ(LevelSetEloi::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	m.newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");

	LeastSquaresGradientComputation grad2D(&m, var_dist, m.getVariable<math::Vector3d,GMDS_NODE>("GMDS_Gradient"));
	LeastSquaresGradientComputation::STATUS result = grad2D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("LeastSquaresGradientComputation_2D_Test1_Result.vtk");

	ASSERT_EQ(GradientComputation2D::SUCCESS, result);
}

TEST(GradientComputationTestClass, LeastSquaresGradientComputation_2D_Test2)
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

	doc.orient2DFaces();

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

	Variable<double> *var_dist = m.newVariable<double,GMDS_NODE>("GMDS_Distance");
	LevelSetExtended ls(&m, markFrontNodes, var_dist);
	LevelSetExtended::STATUS result_ls = ls.execute();
	ASSERT_EQ(LevelSetExtended::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	m.newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");

	LeastSquaresGradientComputation grad2D(&m, var_dist, m.getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	LeastSquaresGradientComputation::STATUS result = grad2D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("LeastSquaresGradientComputation_2D_Test2_Result.vtk");

	ASSERT_EQ(LeastSquaresGradientComputation::SUCCESS, result);
}

TEST(GradientComputationTestClass, LeastSquaresGradientComputation_2D_Test3)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
	                             gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/Pont.vtk";

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

	int markFrontNodes = m.newMark<gmds::Node>();
	for(auto id:bnd_node_ids){
		Node n = m.get<Node>(id);
		double coord_y = n.Y() ;
		double coord_x = n.X() ;
		if ( coord_x == 0 && coord_y == 0 ) {
			// For this test case, the front to advance is the boundary where x²+y²=1
			m.mark<Node>(id,markFrontNodes);
		}
	}

	Variable<double> *var_dist = m.newVariable<double,GMDS_NODE>("GMDS_Distance");
	LevelSetEloi ls(&m, markFrontNodes, var_dist);
	LevelSetEloi::STATUS result_ls = ls.execute();
	ASSERT_EQ(LevelSetEloi::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	m.newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");
	LeastSquaresGradientComputation grad2D(&m, m.getVariable<double,GMDS_NODE>("GMDS_Distance"), m.getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	LeastSquaresGradientComputation::STATUS result = grad2D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("LeastSquaresGradientComputation_2D_Test3_Result.vtk");

	ASSERT_EQ(LeastSquaresGradientComputation::SUCCESS, result);
}

/*----------------------------------------------------------------------------*/




/*----------------------------------------------------------------------------*/
/*             CAS TEST 3D CLASSE LeastSquaresGradientComputation             */
/*----------------------------------------------------------------------------*/

TEST(GradientComputationTestClass, LeastSquaresGradientComputation_3D_Test1)
{
	Mesh m(MeshModel(DIM3 | R | F | E | N |
	                 R2N | F2N | E2N | R2F | F2R |
	                 F2E | E2F | R2E | N2R | N2F | N2E));
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/Cube.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::R);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doctor(&m);
	doctor.buildFacesAndR2F();
	doctor.buildEdgesAndX2E();
	doctor.updateUpwardConnectivity();

	int markFrontNodes = m.newMark<gmds::Node>();

	// Initialisation de la marque pour noter quels fronts sont à avancer
	for(auto id:m.nodes()){
		Node n = m.get<Node>(id);
		double coord_x = n.X() ;
		if ( abs(coord_x) <= pow(10,-6) ) {
			// For this test case, the front to advance is the boundary where x=0
			m.mark<Node>(id,markFrontNodes);
		}
	}

	std::cout << "Fin de l'initialisation des marques" << std::endl ;

	Variable<double> *var_dist = m.newVariable<double,GMDS_NODE>("GMDS_Distance");
	LevelSetEloi ls(&m, markFrontNodes, var_dist);
	LevelSetEloi::STATUS result_ls = ls.execute();
	ASSERT_EQ(LevelSetEloi::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	m.newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");

	LeastSquaresGradientComputation grad3D(&m, var_dist, m.getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	LeastSquaresGradientComputation::STATUS result = grad3D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("LeastSquaresGradientComputation_3D_Test1_Result.vtk");

	ASSERT_EQ(LeastSquaresGradientComputation::SUCCESS, result);
}

TEST(GradientComputationTestClass, LeastSquaresGradientComputation_3D_Test2)
{
	Mesh m(MeshModel(DIM3 | R | F | E | N |
	                 R2N | F2N | E2N | R2F | F2R |
	                 F2E | E2F | R2E | N2R | N2F | N2E));
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/B0.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::R);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doctor(&m);
	doctor.buildFacesAndR2F();
	doctor.buildEdgesAndX2E();
	doctor.updateUpwardConnectivity();

	// Initialisation des marques sur les bords internes et externes pour noter quels fronts sont à avancer
	int markFrontNodes = m.newMark<gmds::Node>();
	for(auto id:m.nodes()){
		Node n = m.get<Node>(id);
		double coord_y = n.Y() ;
		double coord_x = n.X() ;
		double rayon;
		rayon = sqrt( (pow(coord_x, 2) + pow(coord_y + 2.5, 2)) ) ;
		if ( (rayon - 2.5) < pow(10,-3)) {
			// For this test case, the front to advance is the boundary where x²+y²=2.5
			m.mark<Node>(id,markFrontNodes);
		}
	}


	Variable<double> *var_dist = m.newVariable<double,GMDS_NODE>("GMDS_Distance");
	LevelSetEloi ls(&m, markFrontNodes, var_dist);
	LevelSetEloi::STATUS result_ls = ls.execute();
	ASSERT_EQ(LevelSetEloi::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	m.newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");

	LeastSquaresGradientComputation grad3D(&m, var_dist, m.getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	LeastSquaresGradientComputation::STATUS result = grad3D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("LeastSquaresGradientComputation_3D_Test2_Result.vtk");

	ASSERT_EQ(LeastSquaresGradientComputation::SUCCESS, result);
}

TEST(GradientComputationTestClass, LeastSquaresGradientComputation_3D_Test3)
{
	/* Cas test B0 (3D) avec Level Set calculé par la méthode LevelSetCombined
	 * 2 LS sont calculés, une de l'intérieur vers l'extérieur, une de l'extérieur
	 * vers l'intérieur puis elles sont combinées.
	 */
	Mesh m(MeshModel(DIM3 | R | F | E | N |
	                 R2N | F2N | E2N | R2F | F2R |
	                 F2E | E2F | R2E | N2R | N2F | N2E));
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/B0.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::R);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doctor(&m);
	doctor.buildFacesAndR2F();
	doctor.buildEdgesAndX2E();
	doctor.updateUpwardConnectivity();

	int markFrontNodesInt = m.newMark<gmds::Node>();
	int markFrontNodesOut = m.newMark<gmds::Node>();

	for(auto id:m.nodes()){
		Node n = m.get<Node>(id);
		double coord_y = n.Y() ;
		double coord_x = n.X() ;
		double rayon;
		rayon = sqrt( (pow(coord_x, 2) + pow(coord_y + 2.5, 2)) ) ;
		if ( (rayon - 2.5) < pow(10,-3)) {
			// For this test case, the front to advance is the boundary where x²+y²=2.5
			m.mark<Node>(id,markFrontNodesInt);
		}
		else if (coord_x == -5 || coord_x == 5 || coord_y == 2.5) {
			m.mark<Node>(id,markFrontNodesOut);
		}
	}

	std::cout << "Fin de l'initialisation des marques" << std::endl ;

	m.newVariable<double,GMDS_NODE>("GMDS_Distance");
	m.newVariable<double,GMDS_NODE>("GMDS_Distance_Int");
	m.newVariable<double,GMDS_NODE>("GMDS_Distance_Out");
	LevelSetCombined lsCombined(&m, markFrontNodesInt, markFrontNodesOut,
	                            m.getVariable<double,GMDS_NODE>("GMDS_Distance"),
	                            m.getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
	                            m.getVariable<double,GMDS_NODE>("GMDS_Distance_Out"));
	LevelSetCombined::STATUS result_ls = lsCombined.execute();

	m.unmarkAll<Node>(markFrontNodesInt);
	m.freeMark<Node>(markFrontNodesInt);
	m.unmarkAll<Node>(markFrontNodesOut);
	m.freeMark<Node>(markFrontNodesOut);

	m.newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");
	LeastSquaresGradientComputation grad3D(&m, m.getVariable<double,GMDS_NODE>("GMDS_Distance"),
	   m.getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	LeastSquaresGradientComputation::STATUS result = grad3D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("LeastSquaresGradientComputation_3D_Test3_Result.vtk");

	ASSERT_EQ(LeastSquaresGradientComputation::SUCCESS, result);

}

/*----------------------------------------------------------------------------*/
