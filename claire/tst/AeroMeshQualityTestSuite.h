//
// Created by rochec on 20/05/2022.
//

/*----------------------------------------------------------------------------*/
#include <gmds/claire/AeroMeshQuality.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gtest/gtest.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

#ifndef GMDS_AEROMESHQUALITYTESTSUITE_H
#define GMDS_AEROMESHQUALITYTESTSUITE_H

TEST(ClaireTestClass, Test_ConditionQUAD)
{
	// Test
	Mesh m(MeshModel(DIM3 | R | F | N | F2N | N2F));

	GridBuilder gb(&m,2);
	gb.execute(3,1.0, 4, 1.0);

	ASSERT_EQ(m.getNbNodes(),12);
	ASSERT_EQ(m.getNbFaces(),6);

	//==================================================================
	// MESH PREPARATION
	//==================================================================
	MeshDoctor doctor(&m);
	doctor.updateUpwardConnectivity();

	double eps(pow(10,-6));

	for (auto f_id:m.faces()) {
		Face f = m.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double cond = math::AeroMeshQuality::ConditionQUAD(&m, nodes[0].id(), nodes[1].id(), nodes[2].id(), nodes[3].id());
		ASSERT_TRUE(abs(cond - 1.0) < eps);
	}
	//std::cout << "cond : " << cond << std::endl;


	// Test
	Mesh m2(MeshModel(DIM3 | R | F | N | F2N | N2F));

	GridBuilder gb2(&m2,2);
	gb2.execute(3,1.0, 4, 1.0);

	ASSERT_EQ(m2.getNbNodes(),12);
	ASSERT_EQ(m2.getNbFaces(),6);

	//==================================================================
	// MESH PREPARATION
	//==================================================================
	MeshDoctor doctor2(&m2);
	doctor2.updateUpwardConnectivity();

	Node n = m2.get<Node>(5);
	n.setX(0.5);

	n = m2.get<Node>(6);
	n.setY(1.6);

	std::map<TCellID, double> map_mq;
	map_mq[0] = 1.5;
	map_mq[1] = 1.43333;
	map_mq[2] = 1.11429;
	map_mq[3] = 1.16667;
	map_mq[4] = 2.2125;
	map_mq[5] = 1.11429;

	double eps2(pow(10,-3));

	for (auto f_id:m2.faces()) {
		Face f = m.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double cond = math::AeroMeshQuality::ConditionQUAD(&m2, nodes[0].id(), nodes[1].id(), nodes[2].id(), nodes[3].id());
		ASSERT_TRUE(abs(cond - map_mq[f_id]) < eps2);
	}

	//std::cout << "ASSERT_TRUE(abs(cond - " << cond << ") < eps); " << std::endl;

	IGMeshIOService ioService_geom(&m2);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("AeroMeshQuality_Condition.vtk");

}

TEST(ClaireTestClass, Test_EdgeRatioQUAD)
{
	// Test
	Mesh m(MeshModel(DIM3 | R | F | N | F2N | N2F));

	GridBuilder gb(&m,2);
	gb.execute(3,1.0, 4, 1.0);

	ASSERT_EQ(m.getNbNodes(),12);
	ASSERT_EQ(m.getNbFaces(),6);

	//==================================================================
	// MESH PREPARATION
	//==================================================================
	MeshDoctor doctor(&m);
	doctor.updateUpwardConnectivity();

	double eps(pow(10,-6));

	for (auto f_id:m.faces()) {
		Face f = m.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double cond = math::AeroMeshQuality::EdgeRatioQUAD(&m, nodes[0].id(), nodes[1].id(), nodes[2].id(), nodes[3].id());
		ASSERT_TRUE(abs(cond - 1.0) < eps);
	}
	//std::cout << "cond : " << cond << std::endl;


	// Test
	Mesh m2(MeshModel(DIM3 | R | F | N | F2N | N2F));

	GridBuilder gb2(&m2,2);
	gb2.execute(3,1.0, 4, 1.0);

	ASSERT_EQ(m2.getNbNodes(),12);
	ASSERT_EQ(m2.getNbFaces(),6);

	//==================================================================
	// MESH PREPARATION
	//==================================================================
	MeshDoctor doctor2(&m2);
	doctor2.updateUpwardConnectivity();

	Node n = m2.get<Node>(5);
	n.setX(0.5);

	n = m2.get<Node>(6);
	n.setY(1.6);

	std::map<TCellID, double> map_mq;
	map_mq[0] = 2.23607;
	map_mq[1] = 2.15407;
	map_mq[2] = 1.4;
	map_mq[3] = 1.5;
	map_mq[4] = 1.92055;
	map_mq[5] = 1.4;

	double eps2(pow(10,-3));

	for (auto f_id:m2.faces()) {
		Face f = m.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double mq = math::AeroMeshQuality::EdgeRatioQUAD(&m2, nodes[0].id(), nodes[1].id(), nodes[2].id(), nodes[3].id());
		//std::cout << f_id << " " << mq << std::endl;
		ASSERT_TRUE(abs(mq - map_mq[f_id]) < eps2);
	}

	IGMeshIOService ioService_geom(&m2);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("AeroMeshQuality_EdgeRatio.vtk");

}

TEST(ClaireTestClass, Test_JacobianQUAD)
{
	// Test
	Mesh m(MeshModel(DIM3 | R | F | N | F2N | N2F));

	GridBuilder gb(&m,2);
	gb.execute(3,1.0, 4, 1.0);

	ASSERT_EQ(m.getNbNodes(),12);
	ASSERT_EQ(m.getNbFaces(),6);

	//==================================================================
	// MESH PREPARATION
	//==================================================================
	MeshDoctor doctor(&m);
	doctor.updateUpwardConnectivity();

	double eps(pow(10,-6));

	for (auto f_id:m.faces()) {
		Face f = m.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double mq = math::AeroMeshQuality::JacobianQUAD(&m, nodes[0].id(), nodes[1].id(), nodes[2].id(), nodes[3].id());
		//std::cout << f_id <<  " " << mq << std::endl;
		ASSERT_TRUE(abs(mq - 1.0) < eps);
	}
	std::cout << "----------------" << std::endl;


	// Test
	Mesh m2(MeshModel(DIM3 | R | F | N | F2N | N2F));

	GridBuilder gb2(&m2,2);
	gb2.execute(3,1.0, 4, 1.0);

	ASSERT_EQ(m2.getNbNodes(),12);
	ASSERT_EQ(m2.getNbFaces(),6);

	//==================================================================
	// MESH PREPARATION
	//==================================================================
	MeshDoctor doctor2(&m2);
	doctor2.updateUpwardConnectivity();

	Node n = m2.get<Node>(5);
	n.setX(0.5);

	n = m2.get<Node>(6);
	n.setY(1.6);

	std::map<TCellID, double> map_mq;
	map_mq[0] = 0.5;
	map_mq[1] = 0.3;
	map_mq[2] = 1.0;
	map_mq[3] = 1.0;
	map_mq[4] = 0.4;
	map_mq[5] = 1.0;

	double eps2(pow(10,-3));

	for (auto f_id:m2.faces()) {
		Face f = m.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double mq = math::AeroMeshQuality::JacobianQUAD(&m2, nodes[0].id(), nodes[1].id(), nodes[2].id(), nodes[3].id());
		std::cout << f_id << " " << mq << std::endl;
		ASSERT_TRUE(abs(mq - map_mq[f_id]) < eps2);
	}

	IGMeshIOService ioService_geom(&m2);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("AeroMeshQuality_Jacobian.vtk");

}

TEST(ClaireTestClass, Test_ScaledJacobianQUAD)
{
	// Test
	Mesh m(MeshModel(DIM3 | R | F | N | F2N | N2F));

	GridBuilder gb(&m,2);
	gb.execute(3,1.0, 4, 1.0);

	ASSERT_EQ(m.getNbNodes(),12);
	ASSERT_EQ(m.getNbFaces(),6);

	//==================================================================
	// MESH PREPARATION
	//==================================================================
	MeshDoctor doctor(&m);
	doctor.updateUpwardConnectivity();

	double eps(pow(10,-6));

	for (auto f_id:m.faces()) {
		Face f = m.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double mq = math::AeroMeshQuality::ScaledJacobianQUAD(&m, nodes[0].id(), nodes[1].id(), nodes[2].id(), nodes[3].id());
		//std::cout << f_id <<  " " << mq << std::endl;
		ASSERT_TRUE(abs(mq - 1.0) < eps);
	}
	std::cout << "----------------" << std::endl;


	// Test
	Mesh m2(MeshModel(DIM3 | R | F | N | F2N | N2F));

	GridBuilder gb2(&m2,2);
	gb2.execute(3,1.0, 4, 1.0);

	ASSERT_EQ(m2.getNbNodes(),12);
	ASSERT_EQ(m2.getNbFaces(),6);

	//==================================================================
	// MESH PREPARATION
	//==================================================================
	MeshDoctor doctor2(&m2);
	doctor2.updateUpwardConnectivity();

	Node n = m2.get<Node>(5);
	n.setX(0.5);

	n = m2.get<Node>(6);
	n.setY(1.6);

	std::map<TCellID, double> map_mq;
	map_mq[0] = 0.894427;
	map_mq[1] = 0.768221;
	map_mq[2] = 0.928477;
	map_mq[3] = 0.894427;
	map_mq[4] = 0.475517;
	map_mq[5] = 0.928477;

	double eps2(pow(10,-3));

	for (auto f_id:m2.faces()) {
		Face f = m.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double mq = math::AeroMeshQuality::ScaledJacobianQUAD(&m2, nodes[0].id(), nodes[1].id(), nodes[2].id(), nodes[3].id());
		//std::cout << f_id << " " << mq << std::endl;
		ASSERT_TRUE(abs(mq - map_mq[f_id]) < eps2);
	}

	IGMeshIOService ioService_geom(&m2);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("AeroMeshQuality_Jacobian.vtk");

}

TEST(ClaireTestClass, Test_ShapeQUAD)
{
	// Test
	Mesh m(MeshModel(DIM3 | R | F | N | F2N | N2F));

	GridBuilder gb(&m,2);
	gb.execute(3,1.0, 4, 1.0);

	ASSERT_EQ(m.getNbNodes(),12);
	ASSERT_EQ(m.getNbFaces(),6);

	//==================================================================
	// MESH PREPARATION
	//==================================================================
	MeshDoctor doctor(&m);
	doctor.updateUpwardConnectivity();

	double eps(pow(10,-6));

	for (auto f_id:m.faces()) {
		Face f = m.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double mq = math::AeroMeshQuality::ShapeQUAD(&m, nodes[0].id(), nodes[1].id(), nodes[2].id(), nodes[3].id());
		//std::cout << f_id <<  " " << mq << std::endl;
		ASSERT_TRUE(abs(mq - 1.0) < eps);
	}
	std::cout << "----------------" << std::endl;


	// Test
	Mesh m2(MeshModel(DIM3 | R | F | N | F2N | N2F));

	GridBuilder gb2(&m2,2);
	gb2.execute(3,1.0, 4, 1.0);

	ASSERT_EQ(m2.getNbNodes(),12);
	ASSERT_EQ(m2.getNbFaces(),6);

	//==================================================================
	// MESH PREPARATION
	//==================================================================
	MeshDoctor doctor2(&m2);
	doctor2.updateUpwardConnectivity();

	Node n = m2.get<Node>(5);
	n.setX(0.5);

	n = m2.get<Node>(6);
	n.setY(1.6);

	std::map<TCellID, double> map_mq;
	map_mq[0] = 0.666667;
	map_mq[1] = 0.697674;
	map_mq[2] = 0.897436;
	map_mq[3] = 0.857143;
	map_mq[4] = 0.451977;
	map_mq[5] = 0.897436;

	double eps2(pow(10,-3));

	for (auto f_id:m2.faces()) {
		Face f = m.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double mq = math::AeroMeshQuality::ShapeQUAD(&m2, nodes[0].id(), nodes[1].id(), nodes[2].id(), nodes[3].id());
		//std::cout << f_id << " " << mq << std::endl;
		ASSERT_TRUE(abs(mq - map_mq[f_id]) < eps2);
	}

	IGMeshIOService ioService_geom(&m2);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("AeroMeshQuality_Jacobian.vtk");

}

TEST(ClaireTestClass, Test_SkewQUAD)
{
	// Test
	Mesh m(MeshModel(DIM3 | R | F | N | F2N | N2F));

	GridBuilder gb(&m,2);
	gb.execute(3,1.0, 4, 1.0);

	ASSERT_EQ(m.getNbNodes(),12);
	ASSERT_EQ(m.getNbFaces(),6);

	//==================================================================
	// MESH PREPARATION
	//==================================================================
	MeshDoctor doctor(&m);
	doctor.updateUpwardConnectivity();

	double eps(pow(10,-6));

	for (auto f_id:m.faces()) {
		Face f = m.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double mq = math::AeroMeshQuality::SkewQUAD(&m, nodes[0].id(), nodes[1].id(), nodes[2].id(), nodes[3].id());
		//std::cout << f_id <<  " " << mq << std::endl;
		ASSERT_TRUE(abs(mq - 0.0) < eps);
	}


	// Test
	Mesh m2(MeshModel(DIM3 | R | F | N | F2N | N2F));

	GridBuilder gb2(&m2,2);
	gb2.execute(3,1.0, 4, 1.0);

	ASSERT_EQ(m2.getNbNodes(),12);
	ASSERT_EQ(m2.getNbFaces(),6);

	//==================================================================
	// MESH PREPARATION
	//==================================================================
	MeshDoctor doctor2(&m2);
	doctor2.updateUpwardConnectivity();

	Node n = m2.get<Node>(5);
	n.setX(0.5);

	n = m2.get<Node>(6);
	n.setY(1.6);

	std::map<TCellID, double> map_mq;
	map_mq[0] = 0.242536;
	map_mq[1] = 0.0422699;
	map_mq[2] = 0.196116;
	map_mq[3] = 0.242536;
	map_mq[4] = 0.445328;
	map_mq[5] = 0.196116;

	double eps2(pow(10,-3));

	for (auto f_id:m2.faces()) {
		Face f = m.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double mq = math::AeroMeshQuality::SkewQUAD(&m2, nodes[0].id(), nodes[1].id(), nodes[2].id(), nodes[3].id());
		//std::cout << f_id << " " << mq << std::endl;
		ASSERT_TRUE(abs(mq - map_mq[f_id]) < eps2);
	}

	IGMeshIOService ioService_geom(&m2);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("AeroMeshQuality_Jacobian.vtk");

}

TEST(ClaireTestClass, Test_StretchQUAD)
{
	// Test
	Mesh m(MeshModel(DIM3 | R | F | N | F2N | N2F));

	GridBuilder gb(&m,2);
	gb.execute(3,1.0, 4, 1.0);

	ASSERT_EQ(m.getNbNodes(),12);
	ASSERT_EQ(m.getNbFaces(),6);

	//==================================================================
	// MESH PREPARATION
	//==================================================================
	MeshDoctor doctor(&m);
	doctor.updateUpwardConnectivity();

	double eps(pow(10,-6));

	for (auto f_id:m.faces()) {
		Face f = m.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double mq = math::AeroMeshQuality::StretchQUAD(&m, nodes[0].id(), nodes[1].id(), nodes[2].id(), nodes[3].id());
		//std::cout << f_id <<  " " << mq << std::endl;
		ASSERT_TRUE(abs(mq - 1.0) < eps);
	}


	// Test
	Mesh m2(MeshModel(DIM3 | R | F | N | F2N | N2F));

	GridBuilder gb2(&m2,2);
	gb2.execute(3,1.0, 4, 1.0);

	ASSERT_EQ(m2.getNbNodes(),12);
	ASSERT_EQ(m2.getNbFaces(),6);

	//==================================================================
	// MESH PREPARATION
	//==================================================================
	MeshDoctor doctor2(&m2);
	doctor2.updateUpwardConnectivity();

	Node n = m2.get<Node>(5);
	n.setX(0.5);

	n = m2.get<Node>(6);
	n.setY(1.6);

	std::map<TCellID, double> map_mq;
	map_mq[0] = 0.5;
	map_mq[1] = 0.606339;
	map_mq[2] = 0.821995;
	map_mq[3] = 0.784465;
	map_mq[4] = 0.612686;
	map_mq[5] = 0.821995;

	double eps2(pow(10,-3));

	for (auto f_id:m2.faces()) {
		Face f = m.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double mq = math::AeroMeshQuality::StretchQUAD(&m2, nodes[0].id(), nodes[1].id(), nodes[2].id(), nodes[3].id());
		//std::cout << f_id << " " << mq << std::endl;
		ASSERT_TRUE(abs(mq - map_mq[f_id]) < eps2);
	}

	IGMeshIOService ioService_geom(&m2);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("AeroMeshQuality_Jacobian.vtk");

}

#endif     // GMDS_AEROMESHQUALITYTESTSUITE_H
