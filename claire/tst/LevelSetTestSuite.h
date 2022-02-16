//
// Created by rochec on 13/01/2022.
//

#include <gmds/claire/AbstractLevelSet.h>
#include <gmds/claire/LevelSetEloi.h>
#include <gmds/claire/LevelSetExtended.h>
#include <gmds/claire/LevelSetCombined.h>
#include <gmds/claire/GradientComputation2D.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator.h>
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
/*                   CAS TEST 2D CLASSE LevelSetEloi                          */
/*----------------------------------------------------------------------------*/

TEST(LevelSetTestClass, LevelSetEloi_2D_Test1)
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
	LevelSetEloi::STATUS result = ls.execute();

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("LevelSetEloi_2D_Test1_Result.vtk");

	ASSERT_EQ(LevelSetEloi::SUCCESS, result);
}

TEST(LevelSetTestClass, LevelSetEloi_2D_Test2)
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

	// Source ponctuelle
	for(auto id:bnd_node_ids){
		Node n = m.get<Node>(id);
		double coord_x = n.X() ;
		double coord_y = n.Y() ;
		if (coord_y == 0 && coord_x == 0) {
			// For the square test case, the front to advance is the boundary where y=0
			m.mark<Node>(id,markFrontNodes);
		}
	}

	Variable<double> *var_dist = m.newVariable<double,GMDS_NODE>("GMDS_Distance");
	LevelSetEloi ls(&m, markFrontNodes, var_dist);
	LevelSetEloi::STATUS result = ls.execute();

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("LevelSetEloi_2D_Test2_Result.vtk");

	ASSERT_EQ(LevelSetEloi::SUCCESS, result);
}

TEST(LevelSetTestClass, LevelSetEloi_2D_Test3)
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
		if ( sqrt( pow(coord_y,2) + pow(coord_x,2)) == 1) {
			// For this test case, the front to advance is the boundary where x²+y²=1
			m.mark<Node>(id,markFrontNodes);
			//std::cout << "Noeud marqué :" << id << std::endl;
		}
	}

	Variable<double> *var_dist = m.newVariable<double,GMDS_NODE>("GMDS_Distance");
	LevelSetEloi ls(&m, markFrontNodes, var_dist);
	LevelSetEloi::STATUS result = ls.execute();

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("LevelSetEloi_2D_Test3_Result.vtk");

	ASSERT_EQ(LevelSetEloi::SUCCESS, result);
}

TEST(LevelSetTestClass, LevelSetEloi_2D_Test4)
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
		if ( coord_x == 0 || coord_y == 0 ) {
			// For this test case, the front to advance is the boundary where x²+y²=1
			m.mark<Node>(id,markFrontNodes);
		}
	}

	Variable<double> *var_dist = m.newVariable<double,GMDS_NODE>("GMDS_Distance");
	LevelSetEloi ls(&m, markFrontNodes, var_dist);
	LevelSetEloi::STATUS result = ls.execute();

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("LevelSetEloi_2D_Test4_Result.vtk");

	ASSERT_EQ(LevelSetEloi::SUCCESS, result);
}

TEST(LevelSetTestClass, LevelSetEloi_2D_Test5)
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
	LevelSetEloi::STATUS result = ls.execute();

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("LevelSetEloi_2D_Test5_Result.vtk");

	ASSERT_EQ(LevelSetEloi::SUCCESS, result);
}

/*----------------------------------------------------------------------------*/




/*----------------------------------------------------------------------------*/
/*                   CAS TEST 3D CLASSE LevelSetEloi                          */
/*----------------------------------------------------------------------------*/

TEST(LevelSetTestClass, LevelSetEloi_3D_Test1)
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

	int markFrontNodes = m.newMark<gmds::Node>();

	// Test avec une source ponctuelle
	m.mark<Node>(10,markFrontNodes);

	Variable<double> *var_dist = m.newVariable<double,GMDS_NODE>("GMDS_Distance");
	LevelSetEloi ls(&m, markFrontNodes, var_dist);
	LevelSetEloi::STATUS result = ls.execute();

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("LevelSetEloi_3D_Test1_Result.vtk");

	ASSERT_TRUE(true);

}

/*----------------------------------------------------------------------------*/




/*----------------------------------------------------------------------------*/
/*                   CAS TEST 2D CLASSE LevelSetExtended                      */
/*----------------------------------------------------------------------------*/

TEST(LevelSetTestClass, LevelSetExtended_2D_Test1)
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
	LevelSetExtended ls(&m, markFrontNodes, var_dist);
	LevelSetExtended::STATUS result = ls.execute();

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("LevelSetExtended_2D_Test1_Result.vtk");

	ASSERT_EQ(LevelSetExtended::SUCCESS, result);
}

/*----------------------------------------------------------------------------*/




/*----------------------------------------------------------------------------*/
/*                   CAS TEST 3D CLASSE LevelSetExtended                      */
/*----------------------------------------------------------------------------*/

TEST(LevelSetTestClass, LavelSetExtended_3D_Test1)
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

	// Initialisation de la marque pour noter quels fronts sont à avancer
	int markFrontNodes = m.newMark<gmds::Node>();
	for(auto id:m.nodes()){
		Node n = m.get<Node>(id);
		double coord_y = n.Y() ;
		if ( abs(coord_y) <= pow(10,-6) ) {
			// For this test case, the front to advance is the boundary where x=0
			m.mark<Node>(id,markFrontNodes);
		}
	}

	// Calcul des Level Set
	Variable<double> *var_dist = m.newVariable<double,GMDS_NODE>("GMDS_Distance");
	LevelSetExtended ls(&m, markFrontNodes, var_dist);
	LevelSetExtended::STATUS result = ls.execute();

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("LevelSetExtended_3D_Test1_Result.vtk");

	ASSERT_EQ(LevelSetExtended::SUCCESS, result);
}

/*----------------------------------------------------------------------------*/




/*----------------------------------------------------------------------------*/
/*                    CAS TEST 2D CLASSE LevelSetCombined                     */
/*----------------------------------------------------------------------------*/

TEST(LevelSetTestClass, LevelSetCombined_2D_Test1)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
	                             gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/Carre_maxsize_0.01.vtk";

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

	int markFrontNodesInt = m.newMark<gmds::Node>();
	int markFrontNodesOut = m.newMark<gmds::Node>();
	for(auto id:bnd_node_ids){
		Node n = m.get<Node>(id);
		double coord_y = n.Y() ;
		if (coord_y == 0) {
			// For the square test case, the front to advance is the boundary where y=0
			m.mark<Node>(id,markFrontNodesInt);
		}
		else if (coord_y == 1) {
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
	LevelSetCombined::STATUS result = lsCombined.execute();

	m.unmarkAll<Node>(markFrontNodesInt);
	m.freeMark<Node>(markFrontNodesInt);
	m.unmarkAll<Node>(markFrontNodesOut);
	m.freeMark<Node>(markFrontNodesOut);

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("LevelSetCombined_2D_Test1_Result.vtk");

	ASSERT_EQ(LevelSetCombined::SUCCESS, result);
}

TEST(LevelSetTestClass, LevelSetCombined_2D_Test2)
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
	LevelSetCombined::STATUS result = lsCombined.execute();

	m.unmarkAll<Node>(markFrontNodesInt);
	m.freeMark<Node>(markFrontNodesInt);
	m.unmarkAll<Node>(markFrontNodesOut);
	m.freeMark<Node>(markFrontNodesOut);

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("LevelSetCombined_2D_Test2_Result.vtk");

	ASSERT_EQ(LevelSetCombined::SUCCESS, result);

}

/*----------------------------------------------------------------------------*/




/*----------------------------------------------------------------------------*/
/*                    CAS TEST 3D CLASSE LevelSetCombined                     */
/*----------------------------------------------------------------------------*/

TEST(LevelSetTestClass, LevelSetCombined_3D_Test1)
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
	lsCombined.execute();

	m.unmarkAll<Node>(markFrontNodesInt);
	m.freeMark<Node>(markFrontNodesInt);
	m.unmarkAll<Node>(markFrontNodesOut);
	m.freeMark<Node>(markFrontNodesOut);

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("LevelSetCombined_3D_Test1_Result.vtk");

	ASSERT_TRUE(true);

}

/*----------------------------------------------------------------------------*/





/*----------------------------------------------------------------------------*/
/*        Etude convergence en maillage LevelSetEloi/Extended                 */
/*----------------------------------------------------------------------------*/


TEST(LevelSetTestClass, LevelSet_Cvg_2D_Test1)
{
	std::string dir(TEST_SAMPLES_DIR);

	std::vector<std::string> liste_fichiers_vtk;
	liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.5.vtk");
	liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.25.vtk");
	liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.15.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.1.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.075.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.05.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.025.vtk");	// Trop fin
	//liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.01.vtk"); // Trop fin

	// First, we create the file where we are going to store the info
	std::ofstream stream= std::ofstream("LevelSet_Etude_Convergence.table", std::ios::out);
	//set the numerical precision (number of digits)
	stream.precision(15);
	//Header indicating which type of file it is
	stream << "NbrNoeuds EloiL2 ExtendedL2 EloiMax ExtendedMax\n";

	for(int i=0;i<liste_fichiers_vtk.size();i++) {
		std::cout << "-----------------------------" << std::endl;
		std::cout << "-----------------------------" << std::endl;
		std::cout << "|         Maillage : " << i+1 << "      |" << std::endl;
		std::cout << "-----------------------------" << std::endl;

		// WE READ
		gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
		                             gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

		std::string vtk_file = liste_fichiers_vtk[i];

		gmds::IGMeshIOService ioService(&m);
		gmds::VTKReader vtkReader(&ioService);
		vtkReader.setCellOptions(gmds::N | gmds::F);
		vtkReader.read(vtk_file);

		gmds::MeshDoctor doc(&m);
		doc.buildEdgesAndX2E();
		doc.updateUpwardConnectivity();

		// Get the boundary node ids
		BoundaryOperator2D bnd_op(&m);
		std::vector<TCellID> bnd_node_ids;
		bnd_op.getBoundaryNodes(bnd_node_ids);

		double err = pow(10, -6);
		Variable<double> *var_dist_Eloi = m.newVariable<double, GMDS_NODE>("GMDS_Distance_Eloi");
		Variable<double> *var_dist_Extended = m.newVariable<double, GMDS_NODE>("GMDS_Distance_Extended");
		Variable<double> *var_dist_Exacte = m.newVariable<double, GMDS_NODE>("GMDS_Distance_Exacte");

		int markFrontNodes = m.newMark<gmds::Node>();
		for (auto id : bnd_node_ids) {
			Node n = m.get<Node>(id);
			double coord_y = n.Y();
			double coord_x = n.X();
			if (abs(sqrt(pow(coord_x, 2) + pow(coord_y, 2)) - 0.5) < err) {
				m.mark<Node>(id, markFrontNodes);
			}
		}
		std::cout << "-> Fin initialisation bords" << std::endl;

		LevelSetEloi ls_Eloi(&m, markFrontNodes, var_dist_Eloi);
		LevelSetEloi::STATUS result_Eloi = ls_Eloi.execute();
		std::cout << "-> Fin calcul LS Eloi" << std::endl;

		LevelSetExtended ls_Extended(&m, markFrontNodes, var_dist_Extended);
		LevelSetExtended::STATUS result_Extended = ls_Extended.execute();
		std::cout << "-> Fin calcul LS Extended" << std::endl;

		double err_nL1_Eloi(0);
		double err_nL1_Extended(0);

		double err_nL2_Eloi(0);
		double err_nL2_Extended(0);

		double err_max_Eloi(0);
		double err_max_Extended(0);

		double sum_dist_exacte(0);
		double rayon_int(0.5);
		int nbr_noeuds(0);
		for (auto n_id : m.nodes()) {
			Node n = m.get<Node>(n_id);
			math::Point P = n.point();
			double dist_exacte = sqrt(pow(P.X(), 2) + pow(P.Y(), 2)) - rayon_int;
			var_dist_Exacte->set(n_id, dist_exacte);
			sum_dist_exacte += pow(dist_exacte, 2);
			// Erreur en norme L2
			err_nL2_Eloi += pow(dist_exacte - var_dist_Eloi->value(n_id), 2);
			err_nL2_Extended += pow(dist_exacte - var_dist_Extended->value(n_id), 2);
			// Erreur max
			double err_loc_Eloi = abs(dist_exacte - var_dist_Eloi->value(n_id)) / abs(dist_exacte);
			double err_loc_Extended = abs(dist_exacte - var_dist_Extended->value(n_id)) / abs(dist_exacte);
			if (err_loc_Eloi > err_max_Eloi && (var_dist_Eloi->value(n_id) > err)) {
				err_max_Eloi = err_loc_Eloi;
			}
			if (err_loc_Extended > err_max_Extended && (var_dist_Extended->value(n_id) > err)) {
				err_max_Extended = err_loc_Extended;
			}
			nbr_noeuds++;
		}
		err_nL2_Eloi = sqrt(err_nL2_Eloi / sum_dist_exacte);
		err_nL2_Extended = sqrt(err_nL2_Extended / sum_dist_exacte);

		// Ecriture des erreurs dans le fichier
		stream << nbr_noeuds << " " <<  err_nL2_Eloi << " " << err_nL2_Extended << " " << err_max_Eloi << " " << err_max_Extended << "\n";

		m.unmarkAll<Node>(markFrontNodes);
		m.freeMark<Node>(markFrontNodes);

		gmds::VTKWriter vtkWriter(&ioService);
		vtkWriter.setCellOptions(gmds::N | gmds::F);
		vtkWriter.setDataOptions(gmds::N | gmds::F);
		vtkWriter.write("LevelSet_Cvg_2D_Test1_Result.vtk");

		ASSERT_EQ(LevelSetEloi::SUCCESS, result_Eloi);
		ASSERT_EQ(LevelSetExtended::SUCCESS, result_Extended);
	}

	stream.close();

}

TEST(LevelSetTestClass, LevelSet_Cvg_2D_Test2)
{
	std::string dir(TEST_SAMPLES_DIR);

	std::vector<std::string> liste_fichiers_vtk;
	liste_fichiers_vtk.push_back(dir+"/Aero/Poubelle/Carre_0.1.vtk");
	liste_fichiers_vtk.push_back(dir+"/Aero/Poubelle/Carre_0.05.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/Poubelle/Carre_0.01.vtk");
	liste_fichiers_vtk.push_back(dir+"/Aero/Poubelle/Carre_0.02.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/Poubelle/Carre_0.008.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/Poubelle/Carre_0.005.vtk");	// Trop fin

	// First, we create the file where we are going to store the info
	std::ofstream stream= std::ofstream("LevelSet_Etude_Convergence.table", std::ios::out);
	//set the numerical precision (number of digits)
	stream.precision(15);
	//Header indicating which type of file it is
	stream << "NbrNoeuds EloiL2 ExtendedL2 EloiMax ExtendedMax\n";

	for(int i=0;i<liste_fichiers_vtk.size();i++) {
		std::cout << "-----------------------------" << std::endl;
		std::cout << "-----------------------------" << std::endl;
		std::cout << "|         Maillage : " << i+1 << "      |" << std::endl;
		std::cout << "-----------------------------" << std::endl;

		// WE READ
		gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
		                             gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

		std::string vtk_file = liste_fichiers_vtk[i];

		gmds::IGMeshIOService ioService(&m);
		gmds::VTKReader vtkReader(&ioService);
		vtkReader.setCellOptions(gmds::N | gmds::F);
		vtkReader.read(vtk_file);

		gmds::MeshDoctor doc(&m);
		doc.buildEdgesAndX2E();
		doc.updateUpwardConnectivity();

		// Get the boundary node ids
		BoundaryOperator2D bnd_op(&m);
		std::vector<TCellID> bnd_node_ids;
		bnd_op.getBoundaryNodes(bnd_node_ids);

		double err = pow(10, -6);
		Variable<double> *var_dist_Eloi = m.newVariable<double, GMDS_NODE>("GMDS_Distance_Eloi");
		Variable<double> *var_dist_Extended = m.newVariable<double, GMDS_NODE>("GMDS_Distance_Extended");
		Variable<double> *var_dist_Exacte = m.newVariable<double, GMDS_NODE>("GMDS_Distance_Exacte");

		int markFrontNodes = m.newMark<gmds::Node>();
		for (auto id : bnd_node_ids) {
			Node n = m.get<Node>(id);
			double coord_y = n.Y();
			if (abs(coord_y) < err) {
				m.mark<Node>(id, markFrontNodes);
			}
		}
		std::cout << "-> Fin initialisation bords" << std::endl;

		LevelSetEloi ls_Eloi(&m, markFrontNodes, var_dist_Eloi);
		LevelSetEloi::STATUS result_Eloi = ls_Eloi.execute();
		std::cout << "-> Fin calcul LS Eloi" << std::endl;

		LevelSetExtended ls_Extended(&m, markFrontNodes, var_dist_Extended);
		LevelSetExtended::STATUS result_Extended = ls_Extended.execute();
		std::cout << "-> Fin calcul LS Extended" << std::endl;

		double err_nL1_Eloi(0);
		double err_nL1_Extended(0);

		double err_nL2_Eloi(0);
		double err_nL2_Extended(0);

		double err_max_Eloi(0);
		double err_max_Extended(0);

		double sum_dist_exacte(0);

		int nbr_noeuds(0);
		for (auto n_id : m.nodes()) {
			Node n = m.get<Node>(n_id);
			math::Point P = n.point();
			double dist_exacte = P.Y();
			var_dist_Exacte->set(n_id, dist_exacte);
			sum_dist_exacte += pow(dist_exacte, 2);
			// Erreur en norme L2
			err_nL2_Eloi += pow(dist_exacte - var_dist_Eloi->value(n_id), 2);
			err_nL2_Extended += pow(dist_exacte - var_dist_Extended->value(n_id), 2);
			// Erreur max
			double err_loc_Eloi = abs(dist_exacte - var_dist_Eloi->value(n_id)) / abs(dist_exacte);
			double err_loc_Extended = abs(dist_exacte - var_dist_Extended->value(n_id)) / abs(dist_exacte);
			if (err_loc_Eloi > err_max_Eloi && (dist_exacte > err)) {
				err_max_Eloi = err_loc_Eloi;
			}
			if (err_loc_Extended > err_max_Extended && (dist_exacte > err)) {
				err_max_Extended = err_loc_Extended;
			}
			nbr_noeuds++;
		}
		err_nL2_Eloi = sqrt(err_nL2_Eloi / sum_dist_exacte);
		err_nL2_Extended = sqrt(err_nL2_Extended / sum_dist_exacte);

		// Ecriture des erreurs dans le fichier
		stream << nbr_noeuds << " " <<  err_nL2_Eloi << " " << err_nL2_Extended << " " << err_max_Eloi << " " << err_max_Extended << "\n";

		std::cout << "Nbr de noeuds : " << nbr_noeuds << std::endl;

		m.unmarkAll<Node>(markFrontNodes);
		m.freeMark<Node>(markFrontNodes);

		gmds::VTKWriter vtkWriter(&ioService);
		vtkWriter.setCellOptions(gmds::N | gmds::F);
		vtkWriter.setDataOptions(gmds::N | gmds::F);
		vtkWriter.write("LevelSet_Cvg_2D_Test2_Result.vtk");

		ASSERT_EQ(LevelSetEloi::SUCCESS, result_Eloi);
		ASSERT_EQ(LevelSetExtended::SUCCESS, result_Extended);
	}

	stream.close();

}