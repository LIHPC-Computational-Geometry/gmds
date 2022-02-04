//
// Created by rochec on 26/01/2022.
//

#include <gmds/claire/LevelSet.h>
#include <gmds/claire/LevelSetCombined.h>
#include <gmds/claire/GradientComputation2D.h>
#include <gmds/claire/LeastSquaresGradientComputation.h>
#include <gmds/claire/PointFollowingVectorField2D.h>
#include <gmds/claire/PointFollowingVectorField3D.h>
#include <gmds/claire/PointFollowingVectorFieldOnNodes.h>
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
/*             CAS TEST 2D CLASSE PointFollowingVectorField2D                 */
/*----------------------------------------------------------------------------*/

TEST(PointFollowingVectorFieldTestClass, PointFollowingVectorField2D_Test1)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
	                             gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	//std::string vtk_file = dir+"/Carre.vtk";
	std::string vtk_file = dir+"/Carre_maxsize_0.01.vtk";		// Same geo but refined mesh

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
	LevelSetNaif ls(&m, markFrontNodes);
	LevelSetNaif::STATUS result_ls = ls.execute();
	ASSERT_EQ(LevelSetNaif::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	// Calcul du gradient du champ de Level Set
	GradientComputation2D grad2D(&m, m.getVariable<double,GMDS_NODE>("GMDS_Distance"));
	GradientComputation2D::STATUS result_grad = grad2D.execute();
	ASSERT_EQ(GradientComputation2D::SUCCESS, result_grad);

	// Placement du point P à la distance souhaitée suivant le champ de gradient
	math::Point M(0.5, 0.5, 0.0);
	double distance = 0.3;
	PointFollowingVectorField2D pfvf2D(&m, M, distance, m.getVariable<math::Vector3d ,GMDS_FACE>("GMDS_Gradient_2D"));
	PointFollowingVectorField2D::STATUS result = pfvf2D.execute();

	// Placement du point P à la distance souhaitée suivant le champ de gradient
	M.setXYZ(0.3, 0.3, 0.0);
	distance = 0.5;
	pfvf2D = PointFollowingVectorField2D(&m, M, distance, m.getVariable<math::Vector3d ,GMDS_FACE>("GMDS_Gradient_2D"));
	result = pfvf2D.execute();

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("PointFollowingVectorField2D_Test1_Result.vtk");

	ASSERT_EQ(PointFollowingVectorField2D::SUCCESS, result);
}

/*----------------------------------------------------------------------------*/




/*----------------------------------------------------------------------------*/
/*             CAS TEST 3D CLASSE PointFollowingVectorField3D                 */
/*----------------------------------------------------------------------------*/

TEST(PointFollowingVectorFieldTestClass, PointFollowingVectorField3D_Test1)
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
		double coord_y = n.Y() ;
		if ( abs(coord_y) <= pow(10,-6) ) {
			// For this test case, the front to advance is the boundary where x=0
			m.mark<Node>(id,markFrontNodes);
		}
	}

	std::cout << "Fin de l'initialisation des marques" << std::endl ;

	LevelSetNaif ls(&m, markFrontNodes);
	LevelSetNaif::STATUS result_ls = ls.execute();
	ASSERT_EQ(LevelSetNaif::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	GradientComputation3D grad3D(&m, m.getVariable<double,GMDS_NODE>("GMDS_Distance"));
	GradientComputation3D::STATUS result_grad = grad3D.execute();
	ASSERT_EQ(GradientComputation3D::SUCCESS, result_grad);

	// Placement du point P à la distance souhaitée suivant le champ de gradient
	math::Point M(0.1, 0.0, 0.3);
	double distance = 0.8;
	PointFollowingVectorField3D pfvf3D(&m, M, distance, m.getVariable<math::Vector3d ,GMDS_REGION>("GMDS_Gradient_3D"));
	PointFollowingVectorField3D::STATUS result = pfvf3D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("PointFollowingVectorField3D_Test1_Result.vtk");

	ASSERT_EQ(PointFollowingVectorField3D::SUCCESS, result);
}

/*----------------------------------------------------------------------------*/




/*----------------------------------------------------------------------------*/
/*             CAS TEST 2D CLASSE PointFollowingVectorFieldOnNodes            */
/*----------------------------------------------------------------------------*/

TEST(PointFollowingVectorFieldTestClass, PointFollowingVectorFieldOnNodes_2D_Test1)
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
	LevelSetNaif ls(&m, markFrontNodes);
	LevelSetNaif::STATUS result_ls = ls.execute();
	ASSERT_EQ(LevelSetNaif::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	// Calcul du gradient du champ de Level Set
	LeastSquaresGradientComputation grad2D(&m, m.getVariable<double,GMDS_NODE>("GMDS_Distance"));
	LeastSquaresGradientComputation::STATUS result_grad2D = grad2D.execute();
	ASSERT_EQ(LeastSquaresGradientComputation::SUCCESS, result_grad2D);

	// Placement du point P à la distance souhaitée suivant le champ de gradient
	math::Point M(0.5, 0.0, 0.0);
	double distance = 0.8;
	PointFollowingVectorFieldOnNodes pfvf2D(&m, M, distance, m.getVariable<double,GMDS_NODE>("GMDS_Distance"),
	   m.getVariable<math::Vector3d ,GMDS_NODE>("GMDS_Gradient"));
	PointFollowingVectorFieldOnNodes::STATUS result = pfvf2D.execute();

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("PointFollowingVectorFieldOnNodes_2D_Test1_Result.vtk");

	ASSERT_EQ(PointFollowingVectorFieldOnNodes::SUCCESS, result);
}

TEST(PointFollowingVectorFieldTestClass, PointFollowingVectorFieldOnNodes_2D_Test2)
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

	LevelSetCombined lsCombined(&m, markFrontNodesInt, markFrontNodesOut);
	LevelSetCombined::STATUS result_ls = lsCombined.execute();
	ASSERT_EQ(LevelSetCombined::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodesInt);
	m.freeMark<Node>(markFrontNodesInt);
	m.unmarkAll<Node>(markFrontNodesOut);
	m.freeMark<Node>(markFrontNodesOut);

	LeastSquaresGradientComputation grad2D(&m, m.getVariable<double,GMDS_NODE>("GMDS_Distance_Combined"));
	LeastSquaresGradientComputation::STATUS result_grad = grad2D.execute();
	ASSERT_EQ(LeastSquaresGradientComputation::SUCCESS, result_grad);

	// Placement du point P à la distance souhaitée suivant le champ de gradient
	double ang(M_PI/4);
	math::Point M(-cos(ang), sin(ang), 0.0);
	double distance = 1.0;
	PointFollowingVectorFieldOnNodes pfvf2D(&m, M, distance, m.getVariable<double,GMDS_NODE>("GMDS_Distance_Combined"),
	                                        m.getVariable<math::Vector3d ,GMDS_NODE>("GMDS_Gradient"));
	PointFollowingVectorFieldOnNodes::STATUS result = pfvf2D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("PointFollowingVectorFieldOnNodes_2D_Test2_Result.vtk");

	ASSERT_EQ(PointFollowingVectorFieldOnNodes::SUCCESS, result);

}

/*----------------------------------------------------------------------------*/




/*----------------------------------------------------------------------------*/
/*             CAS TEST 3D CLASSE PointFollowingVectorFieldOnNodes            */
/*----------------------------------------------------------------------------*/

TEST(PointFollowingVectorFieldTestClass, PointFollowingVectorFieldOnNodes_3D_Test1)
{
	// Cas test B0 (3D) avec Level Set calculé par la méthode LevelSetCombined
	// LS sont calculés, une de l'intérieur vers l'extérieur, une de l'extérieur
	//vers l'intérieur puis elles sont combinées.

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

	LevelSetCombined lsCombined(&m, markFrontNodesInt, markFrontNodesOut);
	LevelSetCombined::STATUS result_ls = lsCombined.execute();

	m.unmarkAll<Node>(markFrontNodesInt);
	m.freeMark<Node>(markFrontNodesInt);
	m.unmarkAll<Node>(markFrontNodesOut);
	m.freeMark<Node>(markFrontNodesOut);

	LeastSquaresGradientComputation grad3D(&m, m.getVariable<double,GMDS_NODE>("GMDS_Distance_Combined"));
	LeastSquaresGradientComputation::STATUS result_grad = grad3D.execute();
	ASSERT_EQ(LeastSquaresGradientComputation::SUCCESS, result_grad);

	// Placement du point P à la distance souhaitée suivant le champ de gradient
	math::Point M(2.5*cos(M_PI/4), -2.5 + 2.5*sin(M_PI/4), 0.0);
	double distance = 1.0;
	PointFollowingVectorFieldOnNodes pfvf2D(&m, M, distance, m.getVariable<double,GMDS_NODE>("GMDS_Distance_Combined"),
	   m.getVariable<math::Vector3d ,GMDS_NODE>("GMDS_Gradient"));
	PointFollowingVectorFieldOnNodes::STATUS result = pfvf2D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("PointFollowingVectorFieldOnNodes_3D_Test1_Result.vtk");

	ASSERT_EQ(PointFollowingVectorFieldOnNodes::SUCCESS, result);
}

/*
TEST(PointFollowingVectorFieldTestClass, PointFollowingVectorFieldOnNodes_3D_Test2)
{
	// Cas test bicone.

	// --- Lecture du maillage du bicone 3D ---
	Mesh m(MeshModel(DIM3 | R | F | E | N |
	                 R2N | F2N | E2N | R2F | F2R |
	                 F2E | E2F | R2E | N2R | N2F | N2E));
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/biconique.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::R);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doctor(&m);
	doctor.buildFacesAndR2F();
	doctor.buildEdgesAndX2E();
	doctor.updateUpwardConnectivity();
	// -------------------------------------------

	// --- Lecture du maillage de la peau du maillage ---
	Mesh m_peau_int(MeshModel(DIM3 | R | F | E | N |
	                 R2N | F2N | E2N | R2F | F2R |
	                 F2E | E2F | R2E | N2R | N2F | N2E));
	std::string dir_2(TEST_SAMPLES_DIR);
	std::string vtk_file_2 = dir+"/biconique_peau_int.vtk";

	gmds::IGMeshIOService ioService_2(&m_peau_int);
	gmds::VTKReader vtkReader_2(&ioService_2);
	vtkReader_2.setCellOptions(gmds::N|gmds::R);
	vtkReader_2.read(vtk_file_2);

	gmds::MeshDoctor doctor_2(&m_peau_int);
	doctor_2.buildFacesAndR2F();
	doctor_2.buildEdgesAndX2E();
	doctor_2.updateUpwardConnectivity();
	// ---------------------------------------------------

	int markFrontNodesInt = m.newMark<gmds::Node>();
	int markFrontNodesOut = m.newMark<gmds::Node>();

	for(auto id:m.nodes()){
		Node n = m.get<Node>(id);
		double coord_x = n.X() ;
		double coord_y = n.Y() ;
		double coord_z = n.Z() ;
		double rayon;
		rayon = sqrt( (pow(coord_x, 2) + pow(coord_y, 2) + pow(coord_z - 500.0, 2)) ) ;
		if ( (rayon - 1000.0) < pow(10,-6)) {
			// Pour ce cas test, le front extérieur est sur la sphère englobante
			m.mark<Node>(id,markFrontNodesOut);
		}
		// Marquage des noeuds du front intérieur
		// Comparaison avec le maillage de peau pour savoir quel noeud est dessus
		for (auto id_peau:m_peau_int.nodes()){
			Node n_peau = m_peau_int.get<Node>(id_peau);
			math::Vector3d Vec(n_peau.point()-n.point());
			if ( Vec.norm() < pow(10,-6) ){
				m.mark<Node>(id,markFrontNodesInt);
			}
		}
	}

	std::cout << "Fin de l'initialisation des marques" << std::endl ;

	LevelSetCombined lsCombined(&m, markFrontNodesInt, markFrontNodesOut);
	LevelSetCombined::STATUS result_ls = lsCombined.execute();
	ASSERT_EQ(LevelSetCombined::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodesInt);
	m.freeMark<Node>(markFrontNodesInt);
	m.unmarkAll<Node>(markFrontNodesOut);
	m.freeMark<Node>(markFrontNodesOut);

	LeastSquaresGradientComputation grad3D(&m, m.getVariable<double,GMDS_NODE>("GMDS_Distance_Combined"));
	LeastSquaresGradientComputation::STATUS result_grad = grad3D.execute();
	ASSERT_EQ(LeastSquaresGradientComputation::SUCCESS, result_grad);

	// Placement du point P à la distance souhaitée suivant le champ de gradient
	math::Point M(2.5*cos(M_PI/4), -2.5 + 2.5*sin(M_PI/4), 0.0);
	double distance = 1.0;
	PointFollowingVectorFieldOnNodes pfvf2D(&m, M, distance, m.getVariable<double,GMDS_NODE>("GMDS_Distance_Combined"),
	                                        m.getVariable<math::Vector3d ,GMDS_NODE>("GMDS_Gradient"));
	PointFollowingVectorFieldOnNodes::STATUS result = pfvf2D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("PointFollowingVectorFieldOnNodes_3D_Test1_Result.vtk");

	ASSERT_EQ(PointFollowingVectorFieldOnNodes::SUCCESS, result);
}
*/

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
	LevelSetNaif ls(&m, markFrontNodes);
	LevelSetNaif::STATUS result_ls = ls.execute();
	ASSERT_EQ(LevelSetNaif::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	// Calcul du gradient du champ de Level Set
	LeastSquaresGradientComputation grad2D(&m, m.getVariable<double,GMDS_NODE>("GMDS_Distance"));
	LeastSquaresGradientComputation::STATUS result_grad = grad2D.execute();
	ASSERT_EQ(LeastSquaresGradientComputation::SUCCESS, result_grad);

	// Placement du point P à la distance souhaitée suivant le champ de gradient
	math::Point M(0.5, 0.0, 0.0);
	double distance = 0.3;
	AdvectedPointRK4_2D advpoint(&m, M, distance, m.getVariable<double,GMDS_NODE>("GMDS_Distance"), m.getVariable<math::Vector3d ,GMDS_NODE>("GMDS_Gradient"));
	AdvectedPointRK4_2D::STATUS result = advpoint.execute();

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("PointFollowingVectorField2D_Test1_Result.vtk");

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

	LevelSetCombined lsCombined(&m, markFrontNodesInt, markFrontNodesOut);
	LevelSetCombined::STATUS result_ls = lsCombined.execute();
	ASSERT_EQ(LevelSetCombined::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodesInt);
	m.freeMark<Node>(markFrontNodesInt);
	m.unmarkAll<Node>(markFrontNodesOut);
	m.freeMark<Node>(markFrontNodesOut);

	LeastSquaresGradientComputation grad2D(&m, m.getVariable<double,GMDS_NODE>("GMDS_Distance_Combined"));
	LeastSquaresGradientComputation::STATUS result_grad = grad2D.execute();
	ASSERT_EQ(LeastSquaresGradientComputation::SUCCESS, result_grad);

	// Placement du point P à la distance souhaitée suivant le champ de gradient
	double ang(M_PI/4);
	math::Point M(-cos(ang), sin(ang), 0.0);
	double distance = 0.9;
	AdvectedPointRK4_2D advpoint(&m, M, distance, m.getVariable<double,GMDS_NODE>("GMDS_Distance_Combined"), m.getVariable<math::Vector3d ,GMDS_NODE>("GMDS_Gradient"));
	AdvectedPointRK4_2D::STATUS result = advpoint.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("AdvectedPointRK4_2D_Test2_Result.vtk");

	ASSERT_EQ(AdvectedPointRK4_2D::SUCCESS, result);

}

/*----------------------------------------------------------------------------*/