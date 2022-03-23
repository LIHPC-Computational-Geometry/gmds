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
/*                               TESTS UNITAIRES                              */
/*----------------------------------------------------------------------------*/

TEST(LevelSetTestClass, LevelSet_Test_Unit_Carre)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/Aero/Poubelle/Carre.vtk";

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

	double tol(pow(10,-6));

	int markFrontNodesInt = m.newMark<gmds::Node>();
	int markFrontNodesOut = m.newMark<gmds::Node>();
	for(auto id:bnd_node_ids){
		Node n = m.get<Node>(id);
		double coord_y = n.Y() ;
		if ( abs(coord_y) < tol) {
			// For the square test case, the front to advance is the boundary where y=0
			m.mark<Node>(id,markFrontNodesInt);
		}
		else if ( abs(coord_y-1.0) < tol) {
			m.mark<Node>(id,markFrontNodesOut);
		}
	}


	Variable<double> *var_dist = m.newVariable<double, GMDS_NODE>("GMDS_Distance");
	Variable<double> *var_dist_int = m.newVariable<double,GMDS_NODE>("GMDS_Distance_Int");
	Variable<double> *var_dist_out = m.newVariable<double,GMDS_NODE>("GMDS_Distance_Out");

	LevelSetEloi ls(&m, markFrontNodesInt, var_dist);
	LevelSetEloi::STATUS result = ls.execute();
	ASSERT_EQ(LevelSetEloi::SUCCESS, result);

	double eps = pow(10, -6);
	double eps_2 = pow(10, -3);

	{
		ASSERT_TRUE(abs(var_dist->value(0) - 0) < eps);
		ASSERT_TRUE(abs(var_dist->value(1) - 0) < eps);
		ASSERT_TRUE(abs(var_dist->value(2) - 1) < eps);
		ASSERT_TRUE(abs(var_dist->value(3) - 1) < eps);
		ASSERT_TRUE(abs(var_dist->value(4) - 0) < eps);
		ASSERT_TRUE(abs(var_dist->value(5) - 0) < eps);
		ASSERT_TRUE(abs(var_dist->value(6) - 0) < eps);
		ASSERT_TRUE(abs(var_dist->value(7) - 0) < eps);
		ASSERT_TRUE(abs(var_dist->value(8) - 0) < eps);
		ASSERT_TRUE(abs(var_dist->value(9) - 0) < eps);
		ASSERT_TRUE(abs(var_dist->value(10) - 0) < eps);
		ASSERT_TRUE(abs(var_dist->value(11) - 0.125) < eps);
		ASSERT_TRUE(abs(var_dist->value(12) - 0.25) < eps);
		ASSERT_TRUE(abs(var_dist->value(13) - 0.375) < eps);
		ASSERT_TRUE(abs(var_dist->value(14) - 0.5) < eps);
		ASSERT_TRUE(abs(var_dist->value(15) - 0.625) < eps);
		ASSERT_TRUE(abs(var_dist->value(16) - 0.75) < eps);
		ASSERT_TRUE(abs(var_dist->value(17) - 0.875) < eps);
		ASSERT_TRUE(abs(var_dist->value(18) - 1.03176) < eps_2);
		ASSERT_TRUE(abs(var_dist->value(19) - 1.06903) < eps_2);
		ASSERT_TRUE(abs(var_dist->value(20) - 1.0895) < eps);
		ASSERT_TRUE(abs(var_dist->value(21) - 1.10972) < eps_2);
		ASSERT_TRUE(abs(var_dist->value(22) - 1.08348) < eps_2);
		ASSERT_TRUE(abs(var_dist->value(23) - 1.06808) < eps_2);
		ASSERT_TRUE(abs(var_dist->value(24) - 1.02271) < eps_2);
		ASSERT_TRUE(abs(var_dist->value(25) - 0.875) < eps);
		ASSERT_TRUE(abs(var_dist->value(26) - 0.75) < eps);
		ASSERT_TRUE(abs(var_dist->value(27) - 0.625) < eps);
		ASSERT_TRUE(abs(var_dist->value(28) - 0.5) < eps);
		ASSERT_TRUE(abs(var_dist->value(29) - 0.375) < eps);
		ASSERT_TRUE(abs(var_dist->value(30) - 0.25) < eps);
		ASSERT_TRUE(abs(var_dist->value(31) - 0.125) < eps);
		ASSERT_TRUE(abs(var_dist->value(32) - 0.998903) < eps);
		ASSERT_TRUE(abs(var_dist->value(33) - 0.442778) < eps);
		ASSERT_TRUE(abs(var_dist->value(34) - 0.125) < eps);
		ASSERT_TRUE(abs(var_dist->value(35) - 0.579409) < eps);
		ASSERT_TRUE(abs(var_dist->value(36) - 0.703389) < eps);
		ASSERT_TRUE(abs(var_dist->value(37) - 0.964234) < eps);
		ASSERT_TRUE(abs(var_dist->value(38) - 0.121599) < eps);
		ASSERT_TRUE(abs(var_dist->value(39) - 0.338205) < eps);
		ASSERT_TRUE(abs(var_dist->value(40) - 0.821044) < eps);
		ASSERT_TRUE(abs(var_dist->value(41) - 0.945425) < eps);
		ASSERT_TRUE(abs(var_dist->value(42) - 0.125) < eps);
		ASSERT_TRUE(abs(var_dist->value(43) - 0.197234) < eps);
		ASSERT_TRUE(abs(var_dist->value(44) - 0.95887) < eps);
		ASSERT_TRUE(abs(var_dist->value(45) - 0.876705) < eps);
		ASSERT_TRUE(abs(var_dist->value(46) - 0.879217) < eps);
		ASSERT_TRUE(abs(var_dist->value(47) - 0.788033) < eps);
		ASSERT_TRUE(abs(var_dist->value(48) - 0.754264) < eps);
		ASSERT_TRUE(abs(var_dist->value(49) - 0.669756) < eps);
		ASSERT_TRUE(abs(var_dist->value(50) - 0.665699) < eps);
		ASSERT_TRUE(abs(var_dist->value(51) - 0.552643) < eps);
		ASSERT_TRUE(abs(var_dist->value(52) - 0.544019) < eps);
		ASSERT_TRUE(abs(var_dist->value(53) - 0.432382) < eps);
		ASSERT_TRUE(abs(var_dist->value(54) - 0.426745) < eps);
		ASSERT_TRUE(abs(var_dist->value(55) - 0.419651) < eps);
		ASSERT_TRUE(abs(var_dist->value(56) - 0.504019) < eps);
		ASSERT_TRUE(abs(var_dist->value(57) - 0.30741) < eps);
		ASSERT_TRUE(abs(var_dist->value(58) - 0.379157) < eps);
		ASSERT_TRUE(abs(var_dist->value(59) - 0.299818) < eps);
		ASSERT_TRUE(abs(var_dist->value(60) - 0.308527) < eps);
		ASSERT_TRUE(abs(var_dist->value(61) - 0.41078) < eps);
		ASSERT_TRUE(abs(var_dist->value(62) - 0.480459) < eps);
		ASSERT_TRUE(abs(var_dist->value(63) - 0.60867) < eps);
		ASSERT_TRUE(abs(var_dist->value(64) - 0.710588) < eps);
		ASSERT_TRUE(abs(var_dist->value(65) - 0.734787) < eps);
		ASSERT_TRUE(abs(var_dist->value(66) - 0.835822) < eps);
		ASSERT_TRUE(abs(var_dist->value(67) - 0.81513) < eps);
		ASSERT_TRUE(abs(var_dist->value(68) - 0.323446) < eps);
		ASSERT_TRUE(abs(var_dist->value(69) - 0.528492) < eps);
		ASSERT_TRUE(abs(var_dist->value(70) - 0.574415) < eps);
		ASSERT_TRUE(abs(var_dist->value(71) - 0.697732) < eps);
		ASSERT_TRUE(abs(var_dist->value(72) - 0.859325) < eps);
		ASSERT_TRUE(abs(var_dist->value(73) - 0.125) < eps);
		ASSERT_TRUE(abs(var_dist->value(74) - 0.125) < eps);
		ASSERT_TRUE(abs(var_dist->value(75) - 0.984645) < eps);
		ASSERT_TRUE(abs(var_dist->value(76) - 0.629131) < eps);
		ASSERT_TRUE(abs(var_dist->value(77) - 0.308233) < eps);
		ASSERT_TRUE(abs(var_dist->value(78) - 0.828026) < eps);
		ASSERT_TRUE(abs(var_dist->value(79) - 0.109618) < eps);
		ASSERT_TRUE(abs(var_dist->value(80) - 0.199993) < eps);
		ASSERT_TRUE(abs(var_dist->value(81) - 0.950542) < eps);
		ASSERT_TRUE(abs(var_dist->value(82) - 0.302962) < eps);
		ASSERT_TRUE(abs(var_dist->value(83) - 0.26747) < eps);
		ASSERT_TRUE(abs(var_dist->value(84) - 0.757383) < eps);
		ASSERT_TRUE(abs(var_dist->value(85) - 0.450133) < eps);
		ASSERT_TRUE(abs(var_dist->value(86) - 0.541353) < eps);
		ASSERT_TRUE(abs(var_dist->value(87) - 0.70725) < eps);
		ASSERT_TRUE(abs(var_dist->value(88) - 0.925271) < eps);
		ASSERT_TRUE(abs(var_dist->value(89) - 0.934315) < eps);
		ASSERT_TRUE(abs(var_dist->value(90) - 0.0974435) < eps);
		ASSERT_TRUE(abs(var_dist->value(91) - 0.0974435) < eps);
		ASSERT_TRUE(abs(var_dist->value(92) - 0.631916) < eps);
		ASSERT_TRUE(abs(var_dist->value(93) - 0.216616) < eps);
		ASSERT_TRUE(abs(var_dist->value(94) - 0.216616) < eps);
		ASSERT_TRUE(abs(var_dist->value(95) - 0.215794) < eps);
		ASSERT_TRUE(abs(var_dist->value(96) - 0.210062) < eps);
		ASSERT_TRUE(abs(var_dist->value(97) - 0.210517) < eps);
	}

	/*
	for(auto id:m.nodes()){
	   std::cout << "ASSERT_TRUE(abs(var_dist->value( " << id << ") -  " << var_dist->value(id) << ") < eps);" << std::endl;
	}
	 */


	LevelSetExtended lsExtended(&m, markFrontNodesInt, var_dist);
	LevelSetExtended::STATUS resultExtended = lsExtended.execute();
	ASSERT_EQ(LevelSetExtended::SUCCESS, resultExtended);

	{
		ASSERT_TRUE(abs(var_dist->value(0) - 0) < eps);
		ASSERT_TRUE(abs(var_dist->value(1) - 0) < eps);
		ASSERT_TRUE(abs(var_dist->value(2) - 1) < eps);
		ASSERT_TRUE(abs(var_dist->value(3) - 1) < eps);
		ASSERT_TRUE(abs(var_dist->value(4) - 0) < eps);
		ASSERT_TRUE(abs(var_dist->value(5) - 0) < eps);
		ASSERT_TRUE(abs(var_dist->value(6) - 0) < eps);
		ASSERT_TRUE(abs(var_dist->value(7) - 0) < eps);
		ASSERT_TRUE(abs(var_dist->value(8) - 0) < eps);
		ASSERT_TRUE(abs(var_dist->value(9) - 0) < eps);
		ASSERT_TRUE(abs(var_dist->value(10) - 0) < eps);
		ASSERT_TRUE(abs(var_dist->value(11) - 0.125) < eps);
		ASSERT_TRUE(abs(var_dist->value(12) - 0.25) < eps);
		ASSERT_TRUE(abs(var_dist->value(13) - 0.375) < eps);
		ASSERT_TRUE(abs(var_dist->value(14) - 0.5) < eps);
		ASSERT_TRUE(abs(var_dist->value(15) - 0.625) < eps);
		ASSERT_TRUE(abs(var_dist->value(16) - 0.75) < eps);
		ASSERT_TRUE(abs(var_dist->value(17) - 0.875) < eps);
		ASSERT_TRUE(abs(var_dist->value(18) - 1.00272) < eps_2);
		ASSERT_TRUE(abs(var_dist->value(19) - 1.0013) < eps_2);
		ASSERT_TRUE(abs(var_dist->value(20) - 1.00004) < eps_2);
		ASSERT_TRUE(abs(var_dist->value(21) - 1) < eps_2);
		ASSERT_TRUE(abs(var_dist->value(22) - 1.0004) < eps_2);
		ASSERT_TRUE(abs(var_dist->value(23) - 1.00393) < eps_2);
		ASSERT_TRUE(abs(var_dist->value(24) - 1.00161) < eps_2);
		ASSERT_TRUE(abs(var_dist->value(25) - 0.875) < eps);
		ASSERT_TRUE(abs(var_dist->value(26) - 0.75) < eps);
		ASSERT_TRUE(abs(var_dist->value(27) - 0.625) < eps);
		ASSERT_TRUE(abs(var_dist->value(28) - 0.5) < eps);
		ASSERT_TRUE(abs(var_dist->value(29) - 0.375) < eps);
		ASSERT_TRUE(abs(var_dist->value(30) - 0.25) < eps);
		ASSERT_TRUE(abs(var_dist->value(31) - 0.125) < eps);
		ASSERT_TRUE(abs(var_dist->value(32) - 0.904442) < eps);
		ASSERT_TRUE(abs(var_dist->value(33) - 0.43117) < eps);
		ASSERT_TRUE(abs(var_dist->value(34) - 0.125) < eps);
		ASSERT_TRUE(abs(var_dist->value(35) - 0.566217) < eps);
		ASSERT_TRUE(abs(var_dist->value(36) - 0.687878) < eps);
		ASSERT_TRUE(abs(var_dist->value(37) - 0.904728) < eps);
		ASSERT_TRUE(abs(var_dist->value(38) - 0.121599) < eps);
		ASSERT_TRUE(abs(var_dist->value(39) - 0.332775) < eps);
		ASSERT_TRUE(abs(var_dist->value(40) - 0.798784) < eps);
		ASSERT_TRUE(abs(var_dist->value(41) - 0.898221) < eps);
		ASSERT_TRUE(abs(var_dist->value(42) - 0.125) < eps);
		ASSERT_TRUE(abs(var_dist->value(43) - 0.186957) < eps);
		ASSERT_TRUE(abs(var_dist->value(44) - 0.904085) < eps);
		ASSERT_TRUE(abs(var_dist->value(45) - 0.785524) < eps);
		ASSERT_TRUE(abs(var_dist->value(46) - 0.783773) < eps);
		ASSERT_TRUE(abs(var_dist->value(47) - 0.689818) < eps);
		ASSERT_TRUE(abs(var_dist->value(48) - 0.693098) < eps);
		ASSERT_TRUE(abs(var_dist->value(49) - 0.567282) < eps);
		ASSERT_TRUE(abs(var_dist->value(50) - 0.567484) < eps);
		ASSERT_TRUE(abs(var_dist->value(51) - 0.47163) < eps);
		ASSERT_TRUE(abs(var_dist->value(52) - 0.474811) < eps);
		ASSERT_TRUE(abs(var_dist->value(53) - 0.350126) < eps);
		ASSERT_TRUE(abs(var_dist->value(54) - 0.348181) < eps);
		ASSERT_TRUE(abs(var_dist->value(55) - 0.350443) < eps);
		ASSERT_TRUE(abs(var_dist->value(56) - 0.475435) < eps);
		ASSERT_TRUE(abs(var_dist->value(57) - 0.258975) < eps);
		ASSERT_TRUE(abs(var_dist->value(58) - 0.35107) < eps);
		ASSERT_TRUE(abs(var_dist->value(59) - 0.264378) < eps);
		ASSERT_TRUE(abs(var_dist->value(60) - 0.266767) < eps);
		ASSERT_TRUE(abs(var_dist->value(61) - 0.357867) < eps);
		ASSERT_TRUE(abs(var_dist->value(62) - 0.460186) < eps);
		ASSERT_TRUE(abs(var_dist->value(63) - 0.567576) < eps);
		ASSERT_TRUE(abs(var_dist->value(64) - 0.682182) < eps);
		ASSERT_TRUE(abs(var_dist->value(65) - 0.693693) < eps);
		ASSERT_TRUE(abs(var_dist->value(66) - 0.787817) < eps);
		ASSERT_TRUE(abs(var_dist->value(67) - 0.774394) < eps);
		ASSERT_TRUE(abs(var_dist->value(68) - 0.31317) < eps);
		ASSERT_TRUE(abs(var_dist->value(69) - 0.475579) < eps);
		ASSERT_TRUE(abs(var_dist->value(70) - 0.562806) < eps);
		ASSERT_TRUE(abs(var_dist->value(71) - 0.656996) < eps);
		ASSERT_TRUE(abs(var_dist->value(72) - 0.78509) < eps);
		ASSERT_TRUE(abs(var_dist->value(73) - 0.125) < eps);
		ASSERT_TRUE(abs(var_dist->value(74) - 0.125) < eps);
		ASSERT_TRUE(abs(var_dist->value(75) - 0.908469) < eps);
		ASSERT_TRUE(abs(var_dist->value(76) - 0.568394) < eps);
		ASSERT_TRUE(abs(var_dist->value(77) - 0.258975) < eps);
		ASSERT_TRUE(abs(var_dist->value(78) - 0.812515) < eps);
		ASSERT_TRUE(abs(var_dist->value(79) - 0.109618) < eps);
		ASSERT_TRUE(abs(var_dist->value(80) - 0.194563) < eps);
		ASSERT_TRUE(abs(var_dist->value(81) - 0.902967) < eps);
		ASSERT_TRUE(abs(var_dist->value(82) - 0.267093) < eps);
		ASSERT_TRUE(abs(var_dist->value(83) - 0.254812) < eps);
		ASSERT_TRUE(abs(var_dist->value(84) - 0.689115) < eps);
		ASSERT_TRUE(abs(var_dist->value(85) - 0.43943) < eps);
		ASSERT_TRUE(abs(var_dist->value(86) - 0.469229) < eps);
		ASSERT_TRUE(abs(var_dist->value(87) - 0.689533) < eps);
		ASSERT_TRUE(abs(var_dist->value(88) - 0.909505) < eps);
		ASSERT_TRUE(abs(var_dist->value(89) - 0.910981) < eps);
		ASSERT_TRUE(abs(var_dist->value(90) - 0.0974435) < eps);
		ASSERT_TRUE(abs(var_dist->value(91) - 0.0974435) < eps);
		ASSERT_TRUE(abs(var_dist->value(92) - 0.567154) < eps);
		ASSERT_TRUE(abs(var_dist->value(93) - 0.17524) < eps);
		ASSERT_TRUE(abs(var_dist->value(94) - 0.17524) < eps);
		ASSERT_TRUE(abs(var_dist->value(95) - 0.17524) < eps);
		ASSERT_TRUE(abs(var_dist->value(96) - 0.175601) < eps);
		ASSERT_TRUE(abs(var_dist->value(97) - 0.174648) < eps);
	}

	/*
	for(auto id:m.nodes()){
	   std::cout << "ASSERT_TRUE(abs(var_dist->value( " << id << ") -  " << var_dist->value(id) << ") < eps);" << std::endl;
	}
	*/


	LevelSetCombined lsCombined(&m, markFrontNodesInt, markFrontNodesOut,
	                            m.getVariable<double,GMDS_NODE>("GMDS_Distance"),
	                            m.getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
	                            m.getVariable<double,GMDS_NODE>("GMDS_Distance_Out"));
	LevelSetCombined::STATUS resultComb = lsCombined.execute();
	ASSERT_EQ(LevelSetCombined::SUCCESS, resultComb);

	{
		ASSERT_TRUE(abs(var_dist->value( 0) -  0) < eps);
		ASSERT_TRUE(abs(var_dist->value( 1) -  0) < eps);
		ASSERT_TRUE(abs(var_dist->value( 2) -  1) < eps);
		ASSERT_TRUE(abs(var_dist->value( 3) -  1) < eps);
		ASSERT_TRUE(abs(var_dist->value( 4) -  0) < eps);
		ASSERT_TRUE(abs(var_dist->value( 5) -  0) < eps);
		ASSERT_TRUE(abs(var_dist->value( 6) -  0) < eps);
		ASSERT_TRUE(abs(var_dist->value( 7) -  0) < eps);
		ASSERT_TRUE(abs(var_dist->value( 8) -  0) < eps);
		ASSERT_TRUE(abs(var_dist->value( 9) -  0) < eps);
		ASSERT_TRUE(abs(var_dist->value( 10) -  0) < eps);
		ASSERT_TRUE(abs(var_dist->value( 11) -  0.125) < eps);
		ASSERT_TRUE(abs(var_dist->value( 12) -  0.25) < eps);
		ASSERT_TRUE(abs(var_dist->value( 13) -  0.375) < eps);
		ASSERT_TRUE(abs(var_dist->value( 14) -  0.5) < eps);
		ASSERT_TRUE(abs(var_dist->value( 15) -  0.625) < eps);
		ASSERT_TRUE(abs(var_dist->value( 16) -  0.75) < eps);
		ASSERT_TRUE(abs(var_dist->value( 17) -  0.875) < eps);
		ASSERT_TRUE(abs(var_dist->value( 18) -  1) < eps);
		ASSERT_TRUE(abs(var_dist->value( 19) -  1) < eps);
		ASSERT_TRUE(abs(var_dist->value( 20) -  1) < eps);
		ASSERT_TRUE(abs(var_dist->value( 21) -  1) < eps);
		ASSERT_TRUE(abs(var_dist->value( 22) -  1) < eps);
		ASSERT_TRUE(abs(var_dist->value( 23) -  1) < eps);
		ASSERT_TRUE(abs(var_dist->value( 24) -  1) < eps);
		ASSERT_TRUE(abs(var_dist->value( 25) -  0.875) < eps);
		ASSERT_TRUE(abs(var_dist->value( 26) -  0.75) < eps);
		ASSERT_TRUE(abs(var_dist->value( 27) -  0.625) < eps);
		ASSERT_TRUE(abs(var_dist->value( 28) -  0.5) < eps);
		ASSERT_TRUE(abs(var_dist->value( 29) -  0.375) < eps);
		ASSERT_TRUE(abs(var_dist->value( 30) -  0.25) < eps);
		ASSERT_TRUE(abs(var_dist->value( 31) -  0.125) < eps);
		ASSERT_TRUE(abs(var_dist->value( 32) -  0.878575) < eps);
		ASSERT_TRUE(abs(var_dist->value( 33) -  0.430478) < eps);
		ASSERT_TRUE(abs(var_dist->value( 34) -  0.121006) < eps);
		ASSERT_TRUE(abs(var_dist->value( 35) -  0.564683) < eps);
		ASSERT_TRUE(abs(var_dist->value( 36) -  0.686775) < eps);
		ASSERT_TRUE(abs(var_dist->value( 37) -  0.880632) < eps);
		ASSERT_TRUE(abs(var_dist->value( 38) -  0.118127) < eps);
		ASSERT_TRUE(abs(var_dist->value( 39) -  0.331873) < eps);
		ASSERT_TRUE(abs(var_dist->value( 40) -  0.796619) < eps);
		ASSERT_TRUE(abs(var_dist->value( 41) -  0.879855) < eps);
		ASSERT_TRUE(abs(var_dist->value( 42) -  0.121904) < eps);
		ASSERT_TRUE(abs(var_dist->value( 43) -  0.186657) < eps);
		ASSERT_TRUE(abs(var_dist->value( 44) -  0.878867) < eps);
		ASSERT_TRUE(abs(var_dist->value( 45) -  0.785214) < eps);
		ASSERT_TRUE(abs(var_dist->value( 46) -  0.78377) < eps);
		ASSERT_TRUE(abs(var_dist->value( 47) -  0.670876) < eps);
		ASSERT_TRUE(abs(var_dist->value( 48) -  0.67111) < eps);
		ASSERT_TRUE(abs(var_dist->value( 49) -  0.567279) < eps);
		ASSERT_TRUE(abs(var_dist->value( 50) -  0.567259) < eps);
		ASSERT_TRUE(abs(var_dist->value( 51) -  0.458602) < eps);
		ASSERT_TRUE(abs(var_dist->value( 52) -  0.459748) < eps);
		ASSERT_TRUE(abs(var_dist->value( 53) -  0.350124) < eps);
		ASSERT_TRUE(abs(var_dist->value( 54) -  0.348044) < eps);
		ASSERT_TRUE(abs(var_dist->value( 55) -  0.350428) < eps);
		ASSERT_TRUE(abs(var_dist->value( 56) -  0.460605) < eps);
		ASSERT_TRUE(abs(var_dist->value( 57) -  0.251323) < eps);
		ASSERT_TRUE(abs(var_dist->value( 58) -  0.350613) < eps);
		ASSERT_TRUE(abs(var_dist->value( 59) -  0.256132) < eps);
		ASSERT_TRUE(abs(var_dist->value( 60) -  0.259151) < eps);
		ASSERT_TRUE(abs(var_dist->value( 61) -  0.356467) < eps);
		ASSERT_TRUE(abs(var_dist->value( 62) -  0.450765) < eps);
		ASSERT_TRUE(abs(var_dist->value( 63) -  0.566838) < eps);
		ASSERT_TRUE(abs(var_dist->value( 64) -  0.668215) < eps);
		ASSERT_TRUE(abs(var_dist->value( 65) -  0.675217) < eps);
		ASSERT_TRUE(abs(var_dist->value( 66) -  0.786792) < eps);
		ASSERT_TRUE(abs(var_dist->value( 67) -  0.771363) < eps);
		ASSERT_TRUE(abs(var_dist->value( 68) -  0.312252) < eps);
		ASSERT_TRUE(abs(var_dist->value( 69) -  0.473718) < eps);
		ASSERT_TRUE(abs(var_dist->value( 70) -  0.560376) < eps);
		ASSERT_TRUE(abs(var_dist->value( 71) -  0.654425) < eps);
		ASSERT_TRUE(abs(var_dist->value( 72) -  0.785057) < eps);
		ASSERT_TRUE(abs(var_dist->value( 73) -  0.121101) < eps);
		ASSERT_TRUE(abs(var_dist->value( 74) -  0.121307) < eps);
		ASSERT_TRUE(abs(var_dist->value( 75) -  0.879649) < eps);
		ASSERT_TRUE(abs(var_dist->value( 76) -  0.56837) < eps);
		ASSERT_TRUE(abs(var_dist->value( 77) -  0.250701) < eps);
		ASSERT_TRUE(abs(var_dist->value( 78) -  0.811212) < eps);
		ASSERT_TRUE(abs(var_dist->value( 79) -  0.107494) < eps);
		ASSERT_TRUE(abs(var_dist->value( 80) -  0.194036) < eps);
		ASSERT_TRUE(abs(var_dist->value( 81) -  0.884004) < eps);
		ASSERT_TRUE(abs(var_dist->value( 82) -  0.263636) < eps);
		ASSERT_TRUE(abs(var_dist->value( 83) -  0.250619) < eps);
		ASSERT_TRUE(abs(var_dist->value( 84) -  0.673419) < eps);
		ASSERT_TRUE(abs(var_dist->value( 85) -  0.437155) < eps);
		ASSERT_TRUE(abs(var_dist->value( 86) -  0.458541) < eps);
		ASSERT_TRUE(abs(var_dist->value( 87) -  0.684564) < eps);
		ASSERT_TRUE(abs(var_dist->value( 88) -  0.903229) < eps);
		ASSERT_TRUE(abs(var_dist->value( 89) -  0.903371) < eps);
		ASSERT_TRUE(abs(var_dist->value( 90) -  0.0965897) < eps);
		ASSERT_TRUE(abs(var_dist->value( 91) -  0.0966563) < eps);
		ASSERT_TRUE(abs(var_dist->value( 92) -  0.564303) < eps);
		ASSERT_TRUE(abs(var_dist->value( 93) -  0.175233) < eps);
		ASSERT_TRUE(abs(var_dist->value( 94) -  0.17524) < eps);
		ASSERT_TRUE(abs(var_dist->value( 95) -  0.175171) < eps);
		ASSERT_TRUE(abs(var_dist->value( 96) -  0.175373) < eps);
		ASSERT_TRUE(abs(var_dist->value( 97) -  0.173965) < eps);
	}

	/*
	for(auto id:m.nodes()){
	   std::cout << "ASSERT_TRUE(abs(var_dist->value( " << id << ") -  " << var_dist->value(id) << ") < eps);" << std::endl;
	}
	 */

	m.unmarkAll<Node>(markFrontNodesInt);
	m.freeMark<Node>(markFrontNodesInt);
	m.unmarkAll<Node>(markFrontNodesOut);
	m.freeMark<Node>(markFrontNodesOut);

}

TEST(LevelSetTestClass, LevelSet_Test_Unit_C1_3D)
{
	Mesh m(MeshModel(DIM3 | R | F | E | N |
	                 R2N | F2N | E2N | R2F | F2R |
	                 F2E | E2F | R2E | N2R | N2F | N2E));
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/Aero/3D/C1_3D_0.3.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::R);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doctor(&m);
	doctor.buildFacesAndR2F();
	doctor.buildEdgesAndX2E();
	doctor.updateUpwardConnectivity();

	double eps = pow(10, -5);
	double eps_2 = pow(10, -3);

	int markFrontNodesInt = m.newMark<gmds::Node>();
	int markFrontNodesOut = m.newMark<gmds::Node>();

	for(auto id:m.nodes()){
		Node n = m.get<Node>(id);
		double coord_y = n.Y() ;
		double coord_x = n.X() ;
		double coord_z = n.Z() ;
		double rayon;
		rayon = sqrt( (pow(coord_x, 2) + pow(coord_y, 2) + pow(coord_z,2) ) ) ;
		if ( abs(rayon - 0.5) < pow(10,-3)) {
			// For this test case, the front to advance is the boundary where x²+y²=0.5
			m.mark<Node>(id,markFrontNodesInt);
		}
		else if (abs(rayon - 2.0) < pow(10,-3)) {
			m.mark<Node>(id,markFrontNodesOut);
		}
	}

	std::cout << "Fin de l'initialisation des marques" << std::endl ;

	Variable<double> *var_dist = m.newVariable<double, GMDS_NODE>("GMDS_Distance");
	Variable<double> *var_dist_int = m.newVariable<double,GMDS_NODE>("GMDS_Distance_Int");
	Variable<double> *var_dist_out = m.newVariable<double,GMDS_NODE>("GMDS_Distance_Out");

	LevelSetEloi ls(&m, markFrontNodesInt, var_dist);
	LevelSetEloi::STATUS result = ls.execute();
	ASSERT_EQ(LevelSetEloi::SUCCESS, result);

	/*
	for(auto id:m.nodes()){
	 if( id%10 == 0 ){
	   	std::cout << "ASSERT_NEAR(var_dist->value( " << id << "),  " << var_dist->value(id) << ", eps);" << std::endl;
		}
	}
	 */

	{
		ASSERT_NEAR(var_dist->value( 0),  0, eps);
		ASSERT_NEAR(var_dist->value( 10),  1.61241, eps);
		ASSERT_NEAR(var_dist->value( 20),  1.63179, eps);
		ASSERT_NEAR(var_dist->value( 30),  0, eps);
		ASSERT_NEAR(var_dist->value( 40),  0, eps);
		ASSERT_NEAR(var_dist->value( 50),  0, eps);
		ASSERT_NEAR(var_dist->value( 60),  0, eps);
		ASSERT_NEAR(var_dist->value( 70),  0, eps);
		ASSERT_NEAR(var_dist->value( 80),  0, eps);
		ASSERT_NEAR(var_dist->value( 90),  1.5962, eps);
		ASSERT_NEAR(var_dist->value( 100),  1.65789, eps);
		ASSERT_NEAR(var_dist->value( 110),  1.59394, eps);
		ASSERT_NEAR(var_dist->value( 120),  1.60694, eps);
		ASSERT_NEAR(var_dist->value( 130),  1.59233, eps);
		ASSERT_NEAR(var_dist->value( 140),  1.5918, eps);
		ASSERT_NEAR(var_dist->value( 150),  1.71121, eps);
		ASSERT_NEAR(var_dist->value( 160),  1.63226, eps);
		ASSERT_NEAR(var_dist->value( 170),  1.69452, eps);
		ASSERT_NEAR(var_dist->value( 180),  1.62043, eps);
		ASSERT_NEAR(var_dist->value( 190),  1.5768, eps);
		ASSERT_NEAR(var_dist->value( 200),  1.63219, eps);
		ASSERT_NEAR(var_dist->value( 210),  1.71106, eps);
		ASSERT_NEAR(var_dist->value( 220),  1.67479, eps);
		ASSERT_NEAR(var_dist->value( 230),  1.58522, eps);
		ASSERT_NEAR(var_dist->value( 240),  1.5993, eps);
		ASSERT_NEAR(var_dist->value( 250),  1.63547, eps);
		ASSERT_NEAR(var_dist->value( 260),  1.64956, eps);
		ASSERT_NEAR(var_dist->value( 270),  1.62404, eps);
		ASSERT_NEAR(var_dist->value( 280),  1.63969, eps);
		ASSERT_NEAR(var_dist->value( 290),  1.66174, eps);
		ASSERT_NEAR(var_dist->value( 300),  1.64382, eps);
		ASSERT_NEAR(var_dist->value( 310),  1.60292, eps);
		ASSERT_NEAR(var_dist->value( 320),  1.68418, eps);
		ASSERT_NEAR(var_dist->value( 330),  1.63033, eps);
		ASSERT_NEAR(var_dist->value( 340),  1.63582, eps);
		ASSERT_NEAR(var_dist->value( 350),  1.6229, eps);
		ASSERT_NEAR(var_dist->value( 360),  1.55687, eps);
		ASSERT_NEAR(var_dist->value( 370),  1.61603, eps);
		ASSERT_NEAR(var_dist->value( 380),  1.58987, eps);
		ASSERT_NEAR(var_dist->value( 390),  1.62534, eps);
		ASSERT_NEAR(var_dist->value( 400),  1.5847, eps);
		ASSERT_NEAR(var_dist->value( 410),  1.60256, eps);
		ASSERT_NEAR(var_dist->value( 420),  1.63073, eps);
		ASSERT_NEAR(var_dist->value( 430),  1.59585, eps);
		ASSERT_NEAR(var_dist->value( 440),  1.58763, eps);
		ASSERT_NEAR(var_dist->value( 450),  1.59926, eps);
		ASSERT_NEAR(var_dist->value( 460),  1.66821, eps);
		ASSERT_NEAR(var_dist->value( 470),  1.65529, eps);
		ASSERT_NEAR(var_dist->value( 480),  1.62263, eps);
		ASSERT_NEAR(var_dist->value( 490),  1.65274, eps);
		ASSERT_NEAR(var_dist->value( 500),  1.56724, eps);
		ASSERT_NEAR(var_dist->value( 510),  1.60224, eps);
		ASSERT_NEAR(var_dist->value( 520),  1.65441, eps);
		ASSERT_NEAR(var_dist->value( 530),  1.67438, eps);
		ASSERT_NEAR(var_dist->value( 540),  1.61979, eps);
		ASSERT_NEAR(var_dist->value( 550),  1.62692, eps);
		ASSERT_NEAR(var_dist->value( 560),  1.6159, eps);
		ASSERT_NEAR(var_dist->value( 570),  1.62201, eps);
		ASSERT_NEAR(var_dist->value( 580),  1.60073, eps);
		ASSERT_NEAR(var_dist->value( 590),  1.72696, eps);
		ASSERT_NEAR(var_dist->value( 600),  1.56032, eps);
		ASSERT_NEAR(var_dist->value( 610),  1.62631, eps);
		ASSERT_NEAR(var_dist->value( 620),  1.70401, eps);
		ASSERT_NEAR(var_dist->value( 630),  1.63168, eps);
		ASSERT_NEAR(var_dist->value( 640),  1.68273, eps);
		ASSERT_NEAR(var_dist->value( 650),  1.67028, eps);
		ASSERT_NEAR(var_dist->value( 660),  1.56639, eps);
		ASSERT_NEAR(var_dist->value( 670),  1.6016, eps);
		ASSERT_NEAR(var_dist->value( 680),  1.56633, eps);
		ASSERT_NEAR(var_dist->value( 690),  1.73682, eps);
		ASSERT_NEAR(var_dist->value( 700),  1.74781, eps);
		ASSERT_NEAR(var_dist->value( 710),  1.67324, eps);
		ASSERT_NEAR(var_dist->value( 720),  1.62059, eps);
		ASSERT_NEAR(var_dist->value( 730),  1.6386, eps);
		ASSERT_NEAR(var_dist->value( 740),  1.56522, eps);
		ASSERT_NEAR(var_dist->value( 750),  0.801593, eps);
		ASSERT_NEAR(var_dist->value( 760),  0.795046, eps);
		ASSERT_NEAR(var_dist->value( 770),  0.835696, eps);
		ASSERT_NEAR(var_dist->value( 780),  0.992935, eps);
		ASSERT_NEAR(var_dist->value( 790),  1.06457, eps);
		ASSERT_NEAR(var_dist->value( 800),  1.08728, eps);
		ASSERT_NEAR(var_dist->value( 810),  1.0638, eps);
		ASSERT_NEAR(var_dist->value( 820),  0.504526, eps);
		ASSERT_NEAR(var_dist->value( 830),  0.487168, eps);
		ASSERT_NEAR(var_dist->value( 840),  1.12722, eps);
		ASSERT_NEAR(var_dist->value( 850),  0.586665, eps);
		ASSERT_NEAR(var_dist->value( 860),  0.446255, eps);
		ASSERT_NEAR(var_dist->value( 870),  0.459978, eps);
		ASSERT_NEAR(var_dist->value( 880),  0.483766, eps);
		ASSERT_NEAR(var_dist->value( 890),  1.12859, eps);
		ASSERT_NEAR(var_dist->value( 900),  1.13422, eps);
		ASSERT_NEAR(var_dist->value( 910),  1.18768, eps);
		ASSERT_NEAR(var_dist->value( 920),  1.19167, eps);
		ASSERT_NEAR(var_dist->value( 930),  0.399557, eps);
		ASSERT_NEAR(var_dist->value( 940),  1.1755, eps);
		ASSERT_NEAR(var_dist->value( 950),  1.2125, eps);
		ASSERT_NEAR(var_dist->value( 960),  1.18804, eps);
		ASSERT_NEAR(var_dist->value( 970),  0.39235, eps);
		ASSERT_NEAR(var_dist->value( 980),  1.16178, eps);
		ASSERT_NEAR(var_dist->value( 990),  0.668155, eps);
		ASSERT_NEAR(var_dist->value( 1000),  0.354938, eps);
		ASSERT_NEAR(var_dist->value( 1010),  1.19716, eps);
		ASSERT_NEAR(var_dist->value( 1020),  1.2251, eps);
		ASSERT_NEAR(var_dist->value( 1030),  0.838927, eps);
		ASSERT_NEAR(var_dist->value( 1040),  0.67041, eps);
		ASSERT_NEAR(var_dist->value( 1050),  1.17535, eps);
		ASSERT_NEAR(var_dist->value( 1060),  1.21023, eps);
		ASSERT_NEAR(var_dist->value( 1070),  1.22418, eps);
		ASSERT_NEAR(var_dist->value( 1080),  0.736147, eps);
		ASSERT_NEAR(var_dist->value( 1090),  1.2817, eps);
		ASSERT_NEAR(var_dist->value( 1100),  0.849743, eps);
		ASSERT_NEAR(var_dist->value( 1110),  1.30965, eps);
		ASSERT_NEAR(var_dist->value( 1120),  1.24458, eps);
		ASSERT_NEAR(var_dist->value( 1130),  0.821591, eps);
		ASSERT_NEAR(var_dist->value( 1140),  0.326307, eps);
		ASSERT_NEAR(var_dist->value( 1150),  0.873712, eps);
		ASSERT_NEAR(var_dist->value( 1160),  0.893291, eps);
		ASSERT_NEAR(var_dist->value( 1170),  0.321761, eps);
		ASSERT_NEAR(var_dist->value( 1180),  0.742257, eps);
		ASSERT_NEAR(var_dist->value( 1190),  0.706967, eps);
		ASSERT_NEAR(var_dist->value( 1200),  1.21535, eps);
		ASSERT_NEAR(var_dist->value( 1210),  1.27691, eps);
		ASSERT_NEAR(var_dist->value( 1220),  1.30376, eps);
		ASSERT_NEAR(var_dist->value( 1230),  1.29467, eps);
		ASSERT_NEAR(var_dist->value( 1240),  1.36071, eps);
		ASSERT_NEAR(var_dist->value( 1250),  1.34309, eps);
		ASSERT_NEAR(var_dist->value( 1260),  1.27825, eps);
		ASSERT_NEAR(var_dist->value( 1270),  1.2892, eps);
		ASSERT_NEAR(var_dist->value( 1280),  1.40077, eps);
		ASSERT_NEAR(var_dist->value( 1290),  0.879879, eps);
		ASSERT_NEAR(var_dist->value( 1300),  1.00754, eps);
		ASSERT_NEAR(var_dist->value( 1310),  1.32511, eps);
		ASSERT_NEAR(var_dist->value( 1320),  1.38266, eps);
		ASSERT_NEAR(var_dist->value( 1330),  0.28421, eps);
		ASSERT_NEAR(var_dist->value( 1340),  0.923614, eps);
		ASSERT_NEAR(var_dist->value( 1350),  1.28072, eps);
		ASSERT_NEAR(var_dist->value( 1360),  0.778665, eps);
		ASSERT_NEAR(var_dist->value( 1370),  0.685197, eps);
		ASSERT_NEAR(var_dist->value( 1380),  0.967138, eps);
		ASSERT_NEAR(var_dist->value( 1390),  1.42285, eps);
		ASSERT_NEAR(var_dist->value( 1400),  0.256751, eps);
		ASSERT_NEAR(var_dist->value( 1410),  1.05374, eps);
	}

	LevelSetExtended lsExtended(&m, markFrontNodesInt, var_dist);
	LevelSetExtended::STATUS resultExtended = lsExtended.execute();
	ASSERT_EQ(LevelSetExtended::SUCCESS, resultExtended);

	/*
	for(auto id:m.nodes()){
		if( id%10 == 0 ){
			std::cout << "ASSERT_NEAR(var_dist->value( " << id << "),  " << var_dist->value(id) << ", eps);" << std::endl;
		}
	}
	 */

	{
		ASSERT_NEAR(var_dist->value( 0),  0, eps);
		ASSERT_NEAR(var_dist->value( 10),  1.53548, eps);
		ASSERT_NEAR(var_dist->value( 20),  1.53032, eps);
		ASSERT_NEAR(var_dist->value( 30),  0, eps);
		ASSERT_NEAR(var_dist->value( 40),  0, eps);
		ASSERT_NEAR(var_dist->value( 50),  0, eps);
		ASSERT_NEAR(var_dist->value( 60),  0, eps);
		ASSERT_NEAR(var_dist->value( 70),  0, eps);
		ASSERT_NEAR(var_dist->value( 80),  0, eps);
		ASSERT_NEAR(var_dist->value( 90),  1.54097, eps);
		ASSERT_NEAR(var_dist->value( 100),  1.536, eps);
		ASSERT_NEAR(var_dist->value( 110),  1.51927, eps);
		ASSERT_NEAR(var_dist->value( 120),  1.52062, eps);
		ASSERT_NEAR(var_dist->value( 130),  1.51842, eps);
		ASSERT_NEAR(var_dist->value( 140),  1.50843, eps);
		ASSERT_NEAR(var_dist->value( 150),  1.56375, eps);
		ASSERT_NEAR(var_dist->value( 160),  1.55341, eps);
		ASSERT_NEAR(var_dist->value( 170),  1.54642, eps);
		ASSERT_NEAR(var_dist->value( 180),  1.54171, eps);
		ASSERT_NEAR(var_dist->value( 190),  1.53341, eps);
		ASSERT_NEAR(var_dist->value( 200),  1.52609, eps);
		ASSERT_NEAR(var_dist->value( 210),  1.56359, eps);
		ASSERT_NEAR(var_dist->value( 220),  1.53177, eps);
		ASSERT_NEAR(var_dist->value( 230),  1.51291, eps);
		ASSERT_NEAR(var_dist->value( 240),  1.56595, eps);
		ASSERT_NEAR(var_dist->value( 250),  1.54097, eps);
		ASSERT_NEAR(var_dist->value( 260),  1.54466, eps);
		ASSERT_NEAR(var_dist->value( 270),  1.53272, eps);
		ASSERT_NEAR(var_dist->value( 280),  1.54491, eps);
		ASSERT_NEAR(var_dist->value( 290),  1.54771, eps);
		ASSERT_NEAR(var_dist->value( 300),  1.50661, eps);
		ASSERT_NEAR(var_dist->value( 310),  1.53549, eps);
		ASSERT_NEAR(var_dist->value( 320),  1.52007, eps);
		ASSERT_NEAR(var_dist->value( 330),  1.54585, eps);
		ASSERT_NEAR(var_dist->value( 340),  1.51734, eps);
		ASSERT_NEAR(var_dist->value( 350),  1.5389, eps);
		ASSERT_NEAR(var_dist->value( 360),  1.53579, eps);
		ASSERT_NEAR(var_dist->value( 370),  1.54084, eps);
		ASSERT_NEAR(var_dist->value( 380),  1.54097, eps);
		ASSERT_NEAR(var_dist->value( 390),  1.52327, eps);
		ASSERT_NEAR(var_dist->value( 400),  1.52623, eps);
		ASSERT_NEAR(var_dist->value( 410),  1.53683, eps);
		ASSERT_NEAR(var_dist->value( 420),  1.5443, eps);
		ASSERT_NEAR(var_dist->value( 430),  1.53576, eps);
		ASSERT_NEAR(var_dist->value( 440),  1.54332, eps);
		ASSERT_NEAR(var_dist->value( 450),  1.52772, eps);
		ASSERT_NEAR(var_dist->value( 460),  1.5361, eps);
		ASSERT_NEAR(var_dist->value( 470),  1.56596, eps);
		ASSERT_NEAR(var_dist->value( 480),  1.51088, eps);
		ASSERT_NEAR(var_dist->value( 490),  1.54447, eps);
		ASSERT_NEAR(var_dist->value( 500),  1.54668, eps);
		ASSERT_NEAR(var_dist->value( 510),  1.52736, eps);
		ASSERT_NEAR(var_dist->value( 520),  1.53566, eps);
		ASSERT_NEAR(var_dist->value( 530),  1.54158, eps);
		ASSERT_NEAR(var_dist->value( 540),  1.52833, eps);
		ASSERT_NEAR(var_dist->value( 550),  1.52716, eps);
		ASSERT_NEAR(var_dist->value( 560),  1.52266, eps);
		ASSERT_NEAR(var_dist->value( 570),  1.51141, eps);
		ASSERT_NEAR(var_dist->value( 580),  1.53216, eps);
		ASSERT_NEAR(var_dist->value( 590),  1.54135, eps);
		ASSERT_NEAR(var_dist->value( 600),  1.51866, eps);
		ASSERT_NEAR(var_dist->value( 610),  1.51926, eps);
		ASSERT_NEAR(var_dist->value( 620),  1.52582, eps);
		ASSERT_NEAR(var_dist->value( 630),  1.56913, eps);
		ASSERT_NEAR(var_dist->value( 640),  1.54712, eps);
		ASSERT_NEAR(var_dist->value( 650),  1.55872, eps);
		ASSERT_NEAR(var_dist->value( 660),  1.54097, eps);
		ASSERT_NEAR(var_dist->value( 670),  1.5403, eps);
		ASSERT_NEAR(var_dist->value( 680),  1.52689, eps);
		ASSERT_NEAR(var_dist->value( 690),  1.56117, eps);
		ASSERT_NEAR(var_dist->value( 700),  1.56339, eps);
		ASSERT_NEAR(var_dist->value( 710),  1.52901, eps);
		ASSERT_NEAR(var_dist->value( 720),  1.53793, eps);
		ASSERT_NEAR(var_dist->value( 730),  1.5576, eps);
		ASSERT_NEAR(var_dist->value( 740),  1.51283, eps);
		ASSERT_NEAR(var_dist->value( 750),  0.765228, eps);
		ASSERT_NEAR(var_dist->value( 760),  0.77472, eps);
		ASSERT_NEAR(var_dist->value( 770),  0.752711, eps);
		ASSERT_NEAR(var_dist->value( 780),  0.925648, eps);
		ASSERT_NEAR(var_dist->value( 790),  0.97754, eps);
		ASSERT_NEAR(var_dist->value( 800),  0.96888, eps);
		ASSERT_NEAR(var_dist->value( 810),  0.979349, eps);
		ASSERT_NEAR(var_dist->value( 820),  0.478643, eps);
		ASSERT_NEAR(var_dist->value( 830),  0.487168, eps);
		ASSERT_NEAR(var_dist->value( 840),  1.06428, eps);
		ASSERT_NEAR(var_dist->value( 850),  0.488837, eps);
		ASSERT_NEAR(var_dist->value( 860),  0.446255, eps);
		ASSERT_NEAR(var_dist->value( 870),  0.459978, eps);
		ASSERT_NEAR(var_dist->value( 880),  0.483766, eps);
		ASSERT_NEAR(var_dist->value( 890),  1.10488, eps);
		ASSERT_NEAR(var_dist->value( 900),  1.10751, eps);
		ASSERT_NEAR(var_dist->value( 910),  1.12721, eps);
		ASSERT_NEAR(var_dist->value( 920),  1.14436, eps);
		ASSERT_NEAR(var_dist->value( 930),  0.399557, eps);
		ASSERT_NEAR(var_dist->value( 940),  1.12639, eps);
		ASSERT_NEAR(var_dist->value( 950),  1.12337, eps);
		ASSERT_NEAR(var_dist->value( 960),  1.14909, eps);
		ASSERT_NEAR(var_dist->value( 970),  0.39235, eps);
		ASSERT_NEAR(var_dist->value( 980),  1.12843, eps);
		ASSERT_NEAR(var_dist->value( 990),  0.648654, eps);
		ASSERT_NEAR(var_dist->value( 1000),  0.354938, eps);
		ASSERT_NEAR(var_dist->value( 1010),  1.11989, eps);
		ASSERT_NEAR(var_dist->value( 1020),  1.19771, eps);
		ASSERT_NEAR(var_dist->value( 1030),  0.773162, eps);
		ASSERT_NEAR(var_dist->value( 1040),  0.64068, eps);
		ASSERT_NEAR(var_dist->value( 1050),  1.1583, eps);
		ASSERT_NEAR(var_dist->value( 1060),  1.18579, eps);
		ASSERT_NEAR(var_dist->value( 1070),  1.20493, eps);
		ASSERT_NEAR(var_dist->value( 1080),  0.701456, eps);
		ASSERT_NEAR(var_dist->value( 1090),  1.16043, eps);
		ASSERT_NEAR(var_dist->value( 1100),  0.791094, eps);
		ASSERT_NEAR(var_dist->value( 1110),  1.21254, eps);
		ASSERT_NEAR(var_dist->value( 1120),  1.2, eps);
		ASSERT_NEAR(var_dist->value( 1130),  0.816124, eps);
		ASSERT_NEAR(var_dist->value( 1140),  0.326307, eps);
		ASSERT_NEAR(var_dist->value( 1150),  0.811776, eps);
		ASSERT_NEAR(var_dist->value( 1160),  0.826004, eps);
		ASSERT_NEAR(var_dist->value( 1170),  0.321761, eps);
		ASSERT_NEAR(var_dist->value( 1180),  0.708221, eps);
		ASSERT_NEAR(var_dist->value( 1190),  0.586593, eps);
		ASSERT_NEAR(var_dist->value( 1200),  1.18588, eps);
		ASSERT_NEAR(var_dist->value( 1210),  1.19538, eps);
		ASSERT_NEAR(var_dist->value( 1220),  1.23188, eps);
		ASSERT_NEAR(var_dist->value( 1230),  1.21293, eps);
		ASSERT_NEAR(var_dist->value( 1240),  1.25122, eps);
		ASSERT_NEAR(var_dist->value( 1250),  1.24298, eps);
		ASSERT_NEAR(var_dist->value( 1260),  1.23037, eps);
		ASSERT_NEAR(var_dist->value( 1270),  1.21764, eps);
		ASSERT_NEAR(var_dist->value( 1280),  1.25552, eps);
		ASSERT_NEAR(var_dist->value( 1290),  0.804332, eps);
		ASSERT_NEAR(var_dist->value( 1300),  0.943324, eps);
		ASSERT_NEAR(var_dist->value( 1310),  1.25665, eps);
		ASSERT_NEAR(var_dist->value( 1320),  1.25354, eps);
		ASSERT_NEAR(var_dist->value( 1330),  0.28421, eps);
		ASSERT_NEAR(var_dist->value( 1340),  0.815433, eps);
		ASSERT_NEAR(var_dist->value( 1350),  1.24015, eps);
		ASSERT_NEAR(var_dist->value( 1360),  0.67796, eps);
		ASSERT_NEAR(var_dist->value( 1370),  0.546801, eps);
		ASSERT_NEAR(var_dist->value( 1380),  0.902189, eps);
		ASSERT_NEAR(var_dist->value( 1390),  1.24407, eps);
		ASSERT_NEAR(var_dist->value( 1400),  0.256751, eps);
		ASSERT_NEAR(var_dist->value( 1410),  0.945214, eps);
	}

	LevelSetCombined lsCombined(&m, markFrontNodesInt, markFrontNodesOut,
	                            var_dist,
	                            var_dist_int,
	                            var_dist_out);

	LevelSetCombined::STATUS resultComb = lsCombined.execute();
	ASSERT_EQ(LevelSetCombined::SUCCESS, resultComb);

	/*
	for(auto id:m.nodes()){
	 if( id%10 == 0 ){
	      std::cout << "ASSERT_NEAR(var_dist->value( " << id << "),  " << var_dist->value(id) << ", eps);" << std::endl;
	   }
	}
	 */
	{
		ASSERT_NEAR(var_dist->value( 0),  0, eps);
		ASSERT_NEAR(var_dist->value( 10),  1, eps);
		ASSERT_NEAR(var_dist->value( 20),  1, eps);
		ASSERT_NEAR(var_dist->value( 30),  0, eps);
		ASSERT_NEAR(var_dist->value( 40),  0, eps);
		ASSERT_NEAR(var_dist->value( 50),  0, eps);
		ASSERT_NEAR(var_dist->value( 60),  0, eps);
		ASSERT_NEAR(var_dist->value( 70),  0, eps);
		ASSERT_NEAR(var_dist->value( 80),  0, eps);
		ASSERT_NEAR(var_dist->value( 90),  1, eps);
		ASSERT_NEAR(var_dist->value( 100),  1, eps);
		ASSERT_NEAR(var_dist->value( 110),  1, eps);
		ASSERT_NEAR(var_dist->value( 120),  1, eps);
		ASSERT_NEAR(var_dist->value( 130),  1, eps);
		ASSERT_NEAR(var_dist->value( 140),  1, eps);
		ASSERT_NEAR(var_dist->value( 150),  1, eps);
		ASSERT_NEAR(var_dist->value( 160),  1, eps);
		ASSERT_NEAR(var_dist->value( 170),  1, eps);
		ASSERT_NEAR(var_dist->value( 180),  1, eps);
		ASSERT_NEAR(var_dist->value( 190),  1, eps);
		ASSERT_NEAR(var_dist->value( 200),  1, eps);
		ASSERT_NEAR(var_dist->value( 210),  1, eps);
		ASSERT_NEAR(var_dist->value( 220),  1, eps);
		ASSERT_NEAR(var_dist->value( 230),  1, eps);
		ASSERT_NEAR(var_dist->value( 240),  1, eps);
		ASSERT_NEAR(var_dist->value( 250),  1, eps);
		ASSERT_NEAR(var_dist->value( 260),  1, eps);
		ASSERT_NEAR(var_dist->value( 270),  1, eps);
		ASSERT_NEAR(var_dist->value( 280),  1, eps);
		ASSERT_NEAR(var_dist->value( 290),  1, eps);
		ASSERT_NEAR(var_dist->value( 300),  1, eps);
		ASSERT_NEAR(var_dist->value( 310),  1, eps);
		ASSERT_NEAR(var_dist->value( 320),  1, eps);
		ASSERT_NEAR(var_dist->value( 330),  1, eps);
		ASSERT_NEAR(var_dist->value( 340),  1, eps);
		ASSERT_NEAR(var_dist->value( 350),  1, eps);
		ASSERT_NEAR(var_dist->value( 360),  1, eps);
		ASSERT_NEAR(var_dist->value( 370),  1, eps);
		ASSERT_NEAR(var_dist->value( 380),  1, eps);
		ASSERT_NEAR(var_dist->value( 390),  1, eps);
		ASSERT_NEAR(var_dist->value( 400),  1, eps);
		ASSERT_NEAR(var_dist->value( 410),  1, eps);
		ASSERT_NEAR(var_dist->value( 420),  1, eps);
		ASSERT_NEAR(var_dist->value( 430),  1, eps);
		ASSERT_NEAR(var_dist->value( 440),  1, eps);
		ASSERT_NEAR(var_dist->value( 450),  1, eps);
		ASSERT_NEAR(var_dist->value( 460),  1, eps);
		ASSERT_NEAR(var_dist->value( 470),  1, eps);
		ASSERT_NEAR(var_dist->value( 480),  1, eps);
		ASSERT_NEAR(var_dist->value( 490),  1, eps);
		ASSERT_NEAR(var_dist->value( 500),  1, eps);
		ASSERT_NEAR(var_dist->value( 510),  1, eps);
		ASSERT_NEAR(var_dist->value( 520),  1, eps);
		ASSERT_NEAR(var_dist->value( 530),  1, eps);
		ASSERT_NEAR(var_dist->value( 540),  1, eps);
		ASSERT_NEAR(var_dist->value( 550),  1, eps);
		ASSERT_NEAR(var_dist->value( 560),  1, eps);
		ASSERT_NEAR(var_dist->value( 570),  1, eps);
		ASSERT_NEAR(var_dist->value( 580),  1, eps);
		ASSERT_NEAR(var_dist->value( 590),  1, eps);
		ASSERT_NEAR(var_dist->value( 600),  1, eps);
		ASSERT_NEAR(var_dist->value( 610),  1, eps);
		ASSERT_NEAR(var_dist->value( 620),  1, eps);
		ASSERT_NEAR(var_dist->value( 630),  1, eps);
		ASSERT_NEAR(var_dist->value( 640),  1, eps);
		ASSERT_NEAR(var_dist->value( 650),  1, eps);
		ASSERT_NEAR(var_dist->value( 660),  1, eps);
		ASSERT_NEAR(var_dist->value( 670),  1, eps);
		ASSERT_NEAR(var_dist->value( 680),  1, eps);
		ASSERT_NEAR(var_dist->value( 690),  1, eps);
		ASSERT_NEAR(var_dist->value( 700),  1, eps);
		ASSERT_NEAR(var_dist->value( 710),  1, eps);
		ASSERT_NEAR(var_dist->value( 720),  1, eps);
		ASSERT_NEAR(var_dist->value( 730),  1, eps);
		ASSERT_NEAR(var_dist->value( 740),  1, eps);
		ASSERT_NEAR(var_dist->value( 750),  0.5, eps);
		ASSERT_NEAR(var_dist->value( 760),  0.5, eps);
		ASSERT_NEAR(var_dist->value( 770),  0.5, eps);
		ASSERT_NEAR(var_dist->value( 780),  0.587174, eps);
		ASSERT_NEAR(var_dist->value( 790),  0.632301, eps);
		ASSERT_NEAR(var_dist->value( 800),  0.640251, eps);
		ASSERT_NEAR(var_dist->value( 810),  0.639604, eps);
		ASSERT_NEAR(var_dist->value( 820),  0.311419, eps);
		ASSERT_NEAR(var_dist->value( 830),  0.314836, eps);
		ASSERT_NEAR(var_dist->value( 840),  0.672465, eps);
		ASSERT_NEAR(var_dist->value( 850),  0.318264, eps);
		ASSERT_NEAR(var_dist->value( 860),  0.289654, eps);
		ASSERT_NEAR(var_dist->value( 870),  0.298527, eps);
		ASSERT_NEAR(var_dist->value( 880),  0.313214, eps);
		ASSERT_NEAR(var_dist->value( 890),  0.707843, eps);
		ASSERT_NEAR(var_dist->value( 900),  0.711937, eps);
		ASSERT_NEAR(var_dist->value( 910),  0.719597, eps);
		ASSERT_NEAR(var_dist->value( 920),  0.738993, eps);
		ASSERT_NEAR(var_dist->value( 930),  0.258639, eps);
		ASSERT_NEAR(var_dist->value( 940),  0.728993, eps);
		ASSERT_NEAR(var_dist->value( 950),  0.729708, eps);
		ASSERT_NEAR(var_dist->value( 960),  0.737001, eps);
		ASSERT_NEAR(var_dist->value( 970),  0.252264, eps);
		ASSERT_NEAR(var_dist->value( 980),  0.726585, eps);
		ASSERT_NEAR(var_dist->value( 990),  0.416543, eps);
		ASSERT_NEAR(var_dist->value( 1000),  0.229624, eps);
		ASSERT_NEAR(var_dist->value( 1010),  0.724718, eps);
		ASSERT_NEAR(var_dist->value( 1020),  0.758021, eps);
		ASSERT_NEAR(var_dist->value( 1030),  0.503182, eps);
		ASSERT_NEAR(var_dist->value( 1040),  0.415836, eps);
		ASSERT_NEAR(var_dist->value( 1050),  0.755683, eps);
		ASSERT_NEAR(var_dist->value( 1060),  0.761631, eps);
		ASSERT_NEAR(var_dist->value( 1070),  0.764088, eps);
		ASSERT_NEAR(var_dist->value( 1080),  0.46427, eps);
		ASSERT_NEAR(var_dist->value( 1090),  0.763127, eps);
		ASSERT_NEAR(var_dist->value( 1100),  0.521368, eps);
		ASSERT_NEAR(var_dist->value( 1110),  0.775383, eps);
		ASSERT_NEAR(var_dist->value( 1120),  0.775761, eps);
		ASSERT_NEAR(var_dist->value( 1130),  0.529295, eps);
		ASSERT_NEAR(var_dist->value( 1140),  0.209451, eps);
		ASSERT_NEAR(var_dist->value( 1150),  0.521958, eps);
		ASSERT_NEAR(var_dist->value( 1160),  0.537088, eps);
		ASSERT_NEAR(var_dist->value( 1170),  0.211398, eps);
		ASSERT_NEAR(var_dist->value( 1180),  0.467366, eps);
		ASSERT_NEAR(var_dist->value( 1190),  0.386448, eps);
		ASSERT_NEAR(var_dist->value( 1200),  0.769882, eps);
		ASSERT_NEAR(var_dist->value( 1210),  0.782715, eps);
		ASSERT_NEAR(var_dist->value( 1220),  0.791621, eps);
		ASSERT_NEAR(var_dist->value( 1230),  0.78874, eps);
		ASSERT_NEAR(var_dist->value( 1240),  0.79582, eps);
		ASSERT_NEAR(var_dist->value( 1250),  0.793539, eps);
		ASSERT_NEAR(var_dist->value( 1260),  0.791617, eps);
		ASSERT_NEAR(var_dist->value( 1270),  0.785033, eps);
		ASSERT_NEAR(var_dist->value( 1280),  0.803839, eps);
		ASSERT_NEAR(var_dist->value( 1290),  0.528676, eps);
		ASSERT_NEAR(var_dist->value( 1300),  0.597419, eps);
		ASSERT_NEAR(var_dist->value( 1310),  0.801816, eps);
		ASSERT_NEAR(var_dist->value( 1320),  0.798865, eps);
		ASSERT_NEAR(var_dist->value( 1330),  0.186454, eps);
		ASSERT_NEAR(var_dist->value( 1340),  0.536905, eps);
		ASSERT_NEAR(var_dist->value( 1350),  0.803924, eps);
		ASSERT_NEAR(var_dist->value( 1360),  0.448458, eps);
		ASSERT_NEAR(var_dist->value( 1370),  0.361844, eps);
		ASSERT_NEAR(var_dist->value( 1380),  0.587323, eps);
		ASSERT_NEAR(var_dist->value( 1390),  0.807606, eps);
		ASSERT_NEAR(var_dist->value( 1400),  0.16894, eps);
		ASSERT_NEAR(var_dist->value( 1410),  0.621509, eps);
	}


	m.unmarkAll<Node>(markFrontNodesInt);
	m.freeMark<Node>(markFrontNodesInt);
	m.unmarkAll<Node>(markFrontNodesOut);
	m.freeMark<Node>(markFrontNodesOut);

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("LevelSet_C1_3D.vtk");

}

/*----------------------------------------------------------------------------*/




























/*----------------------------------------------------------------------------*/
/*                   CAS TEST 2D CLASSE LevelSetEloi                          */
/*----------------------------------------------------------------------------*/
/*
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
 */

/*----------------------------------------------------------------------------*/




/*----------------------------------------------------------------------------*/
/*                   CAS TEST 3D CLASSE LevelSetEloi                          */
/*----------------------------------------------------------------------------*/
/*
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
*/
/*----------------------------------------------------------------------------*/




/*----------------------------------------------------------------------------*/
/*                   CAS TEST 2D CLASSE LevelSetExtended                      */
/*----------------------------------------------------------------------------*/
/*
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
*/
/*----------------------------------------------------------------------------*/




/*----------------------------------------------------------------------------*/
/*                   CAS TEST 3D CLASSE LevelSetExtended                      */
/*----------------------------------------------------------------------------*/
/*
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
*/
/*----------------------------------------------------------------------------*/




/*----------------------------------------------------------------------------*/
/*                    CAS TEST 2D CLASSE LevelSetCombined                     */
/*----------------------------------------------------------------------------*/
/*
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
*/
/*----------------------------------------------------------------------------*/




/*----------------------------------------------------------------------------*/
/*                    CAS TEST 3D CLASSE LevelSetCombined                     */
/*----------------------------------------------------------------------------*/
/*
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
*/
/*----------------------------------------------------------------------------*/





/*----------------------------------------------------------------------------*/
/*        Etude convergence en maillage LevelSetEloi/Extended                 */
/*----------------------------------------------------------------------------*/

/*
TEST(LevelSetTestClass, LevelSet_Cvg_2D_Test1)
{
	std::string dir(TEST_SAMPLES_DIR);

	std::vector<std::string> liste_fichiers_vtk;
	liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.5.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.25.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.15.vtk");
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
	//liste_fichiers_vtk.push_back(dir+"/Aero/Poubelle/Carre_0.075.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/Poubelle/Carre_0.06.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/Poubelle/Carre_0.05.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/Poubelle/Carre_0.04.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/Poubelle/Carre_0.02.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/Poubelle/Carre_0.015.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/Poubelle/Carre_0.01.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/Poubelle/Carre_0.009.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/Poubelle/Carre_0.008.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/Poubelle/Carre_0.005.vtk");

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

TEST(LevelSetTestClass, LevelSet_Cvg_2D_Test3)
{
	std::string dir(TEST_SAMPLES_DIR);

	std::vector<std::string> liste_fichiers_vtk;
	//liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.5.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.25.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.15.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.1.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.075.vtk");
	//liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.05.vtk");
	liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.025.vtk");	// Trop fin
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
			if (abs(sqrt(pow(coord_x, 2) + pow(coord_y, 2)) - 2.0) < err) {
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
		double rayon_ext(2.0);
		int nbr_noeuds(0);
		for (auto n_id : m.nodes()) {
			Node n = m.get<Node>(n_id);
			math::Point P = n.point();
			double dist_exacte = rayon_ext - sqrt(pow(P.X(), 2) + pow(P.Y(), 2));
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
		vtkWriter.write("LevelSet_Cvg_2D_Test2_Result.vtk");

		ASSERT_EQ(LevelSetEloi::SUCCESS, result_Eloi);
		ASSERT_EQ(LevelSetExtended::SUCCESS, result_Extended);
	}

	stream.close();

}

TEST(LevelSetTestClass, LevelSet_Cvg_2D_Test4)
{
	std::string dir(TEST_SAMPLES_DIR);

	std::vector<std::string> liste_fichiers_vtk;
	liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.5.vtk");
	liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.25.vtk");
	liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.15.vtk");
	liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.1.vtk");
	liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.075.vtk");
	liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.05.vtk");
	liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.025.vtk");	// Trop fin
	//liste_fichiers_vtk.push_back(dir+"/Aero/C1_2D_0.01.vtk"); // Trop fin

	// First, we create the file where we are going to store the info
	std::ofstream stream= std::ofstream("LSCombinedCvC12D.table", std::ios::out);
	//set the numerical precision (number of digits)
	stream.precision(15);
	//Header indicating which type of file it is
	stream << "SqrtNbrNoeuds NbrNoeuds ExtendedL2 \n";

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
		Variable<double> *var_dist_Comb = m.newVariable<double, GMDS_NODE>("GMDS_Distance");
		Variable<double> *var_dist_Int = m.newVariable<double, GMDS_NODE>("GMDS_Distance_Int");
		Variable<double> *var_dist_Out = m.newVariable<double, GMDS_NODE>("GMDS_Distance_Out");
		Variable<double> *var_dist_Exacte = m.newVariable<double, GMDS_NODE>("GMDS_Distance_Exacte");

		int markFrontNodesInt = m.newMark<gmds::Node>();
		int markFrontNodesOut = m.newMark<gmds::Node>();
		double rayon_int(0.5);
		double rayon_ext(2.0);
		for (auto id : bnd_node_ids) {
			Node n = m.get<Node>(id);
			double coord_y = n.Y();
			double coord_x = n.X();
			if (abs(sqrt(pow(coord_x, 2) + pow(coord_y, 2)) - rayon_ext) < err) {
				m.mark<Node>(id, markFrontNodesOut);
			}
			if (abs(sqrt(pow(coord_x, 2) + pow(coord_y, 2)) - rayon_int) < err) {
				m.mark<Node>(id, markFrontNodesInt);
			}
		}
		std::cout << "-> Fin initialisation bords" << std::endl;

		LevelSetCombined lsCombined(&m, markFrontNodesInt, markFrontNodesOut,
		                            m.getVariable<double,GMDS_NODE>("GMDS_Distance"),
		                            m.getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
		                            m.getVariable<double,GMDS_NODE>("GMDS_Distance_Out"));
		LevelSetCombined::STATUS result = lsCombined.execute();

		double err_nL1(0);
		double err_nL2(0);
		double err_max(0);

		double sum_dist_exacte(0);
		int nbr_noeuds(0);
		for (auto n_id : m.nodes()) {
			Node n = m.get<Node>(n_id);
			math::Point P = n.point();
			double dist_exacte_int =	sqrt(pow(P.X(), 2) + pow(P.Y(), 2)) - rayon_int ;
			double dist_exacte_ext = rayon_ext - sqrt(pow(P.X(), 2) + pow(P.Y(), 2)) ;
			double dist_exacte = dist_exacte_int/(dist_exacte_int + dist_exacte_ext);
			var_dist_Exacte->set(n_id, dist_exacte);
			sum_dist_exacte += pow(dist_exacte, 2);
			// Erreur en norme L2
			err_nL2 += pow(dist_exacte - var_dist_Comb->value(n_id), 2);
			// Erreur max
			double err_loc = abs(dist_exacte - var_dist_Comb->value(n_id)) / abs(dist_exacte);
			if (err_loc > err_max && (var_dist_Comb->value(n_id) > err)) {
				err_max = err_loc;
			}
			nbr_noeuds++;
		}
		err_nL2 = sqrt(err_nL2 / sum_dist_exacte);

		// Ecriture des erreurs dans le fichier
		stream << sqrt(nbr_noeuds) << " " << nbr_noeuds << " " <<  err_nL2 << "\n";

		m.unmarkAll<Node>(markFrontNodesInt);
		m.freeMark<Node>(markFrontNodesInt);
		m.unmarkAll<Node>(markFrontNodesOut);
		m.freeMark<Node>(markFrontNodesOut);

		gmds::VTKWriter vtkWriter(&ioService);
		vtkWriter.setCellOptions(gmds::N | gmds::F);
		vtkWriter.setDataOptions(gmds::N | gmds::F);
		vtkWriter.write("LevelSet_Cvg_2D_Test2_Result.vtk");

		ASSERT_EQ(LevelSetCombined::SUCCESS, result);

	}

	stream.close();

}
 */


/*----------------------------------------------------------------------------*/








