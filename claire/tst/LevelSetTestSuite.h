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

	double eps = pow(10, -6);
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
	   std::cout << "ASSERT_NEAR(var_dist->value( " << id << "),  " << var_dist->value(id) << ", eps);" << std::endl;
	}
	 */

	LevelSetExtended lsExtended(&m, markFrontNodesInt, var_dist);
	LevelSetExtended::STATUS resultExtended = lsExtended.execute();
	ASSERT_EQ(LevelSetExtended::SUCCESS, resultExtended);


	LevelSetCombined lsCombined(&m, markFrontNodesInt, markFrontNodesOut,
	                            var_dist,
	                            var_dist_int,
	                            var_dist_out);

	LevelSetCombined::STATUS resultComb = lsCombined.execute();
	ASSERT_EQ(LevelSetCombined::SUCCESS, resultComb);

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








