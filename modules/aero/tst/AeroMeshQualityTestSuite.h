//
// Created by rochec on 20/05/2022.
//

/*----------------------------------------------------------------------------*/
#include <gmds/aero/AeroMeshQuality.h>
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

TEST(AeroTestClass, Test_minedgelenght)
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
		double min_edge = math::AeroMeshQuality::minlenghtedge(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point());
		ASSERT_TRUE(abs(min_edge - 1.0) < eps);
	}

}

TEST(AeroTestClass, Test_ConditionQUAD)
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
		double cond = math::AeroMeshQuality::ConditionQUAD(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point());
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
		Face f = m2.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double cond = math::AeroMeshQuality::ConditionQUAD(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point());
		ASSERT_TRUE(abs(cond - map_mq[f_id]) < eps2);
	}

	//std::cout << "ASSERT_TRUE(abs(cond - " << cond << ") < eps); " << std::endl;

	IGMeshIOService ioService_geom(&m2);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("AeroMeshQuality_Condition.vtk");

}

TEST(AeroTestClass, Test_EdgeRatioQUAD)
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
		double cond = math::AeroMeshQuality::EdgeRatioQUAD(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point());
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
		Face f = m2.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double mq = math::AeroMeshQuality::EdgeRatioQUAD(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point());
		ASSERT_TRUE(abs(mq - map_mq[f_id]) < eps2);
	}

	IGMeshIOService ioService_geom(&m2);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("AeroMeshQuality_EdgeRatio.vtk");

}

TEST(AeroTestClass, Test_JacobianQUAD)
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
		double mq = math::AeroMeshQuality::JacobianQUAD(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point());
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
		Face f = m2.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double mq = math::AeroMeshQuality::JacobianQUAD(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point());
		//std::cout << f_id << " " << mq << std::endl;
		ASSERT_TRUE(abs(mq - map_mq[f_id]) < eps2);
	}

	IGMeshIOService ioService_geom(&m2);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("AeroMeshQuality_Jacobian.vtk");

}

TEST(AeroTestClass, Test_ScaledJacobianQUAD)
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
		double mq = math::AeroMeshQuality::ScaledJacobianQUAD(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point());
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
		Face f = m2.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double mq = math::AeroMeshQuality::ScaledJacobianQUAD(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point());
		//std::cout << f_id << " " << mq << std::endl;
		ASSERT_TRUE(abs(mq - map_mq[f_id]) < eps2);
	}

	IGMeshIOService ioService_geom(&m2);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("AeroMeshQuality_Jacobian.vtk");

}

TEST(AeroTestClass, Test_ShapeQUAD)
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
		double mq = math::AeroMeshQuality::ShapeQUAD(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point());
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
		Face f = m2.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double mq = math::AeroMeshQuality::ShapeQUAD(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point());
		//std::cout << f_id << " " << mq << std::endl;
		ASSERT_TRUE(abs(mq - map_mq[f_id]) < eps2);
	}

	IGMeshIOService ioService_geom(&m2);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("AeroMeshQuality_Jacobian.vtk");

}

TEST(AeroTestClass, Test_SkewQUAD)
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
		double mq = math::AeroMeshQuality::SkewQUAD(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point());
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
		Face f = m2.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double mq = math::AeroMeshQuality::SkewQUAD(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point());
		//std::cout << f_id << " " << mq << std::endl;
		ASSERT_TRUE(abs(mq - map_mq[f_id]) < eps2);
	}

	IGMeshIOService ioService_geom(&m2);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("AeroMeshQuality_Jacobian.vtk");

}

TEST(AeroTestClass, Test_StretchQUAD)
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
		double mq = math::AeroMeshQuality::StretchQUAD(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point());
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
		Face f = m2.get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		double mq = math::AeroMeshQuality::StretchQUAD(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point());
		//std::cout << f_id << " " << mq << std::endl;
		ASSERT_TRUE(abs(mq - map_mq[f_id]) < eps2);
	}

	IGMeshIOService ioService_geom(&m2);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("AeroMeshQuality_Jacobian.vtk");

}

TEST(AeroTestClass, Test_ScaledJacobianHEX)
{
	math::Point p0(0,0,0);
	math::Point p1(1,0,0);
	math::Point p2(1,1,0);
	math::Point p3(0,1,0);
	math::Point p4(0,0,1);
	math::Point p5(1,0,1);
	math::Point p6(1,1,1);
	math::Point p7(0,1,1);

	//std::cout << "Scaled Jacobian Hex 1: " << math::AeroMeshQuality::ScaledJacobianHEX(p0,p1,p2,p3,p4,p5,p6,p7) << std::endl;
	//std::cout << "Scaled Jacobian Hex 2: " << math::AeroMeshQuality::ScaledJacobianHEX(p0,p3,p2,p1,p4,p7,p6,p5) << std::endl;

	std::cout << "SJ HEX: " << math::AeroMeshQuality::ScaledJacobianHEX(p0,p1,p2,p3,p4,p5,p6,p7) << std::endl;
	ASSERT_FLOAT_EQ( math::AeroMeshQuality::ScaledJacobianHEX(p0,p1,p2,p3,p4,p5,p6,p7) ,1.0);


	math::Point p10(0,0,0);
	math::Point p11(1,0,0);
	math::Point p12(1,1,0);
	math::Point p13(0,1,0);
	std::cout << "SJ QUAD: " << math::AeroMeshQuality::ScaledJacobianQUAD(p10,p11,p12,p13);

}

#endif     // GMDS_AEROMESHQUALITYTESTSUITE_H
