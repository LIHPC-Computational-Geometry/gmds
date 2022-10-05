/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include "gmds/blocking/CGNSWriter.h"
#include "gmds/ig/Blocking2D.h"
#include "gmds/ig/MeshDoctor.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKReader.h"
#include "gmds/io/VTKWriter.h"
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

TEST(CGNSTestSuite, test_monobloc){

	Blocking2D m;
	Node n1 = m.newBlockCorner(0,0);
	Node n2 = m.newBlockCorner(1,0);
	Node n3 = m.newBlockCorner(1,1);
	Node n4=  m.newBlockCorner(0,1);

	Blocking2D::Block b0 = m.newBlock(n1,n2,n3,n4);

	b0.setNbDiscretizationI(11);
	b0.setNbDiscretizationJ(11);

	m.initializeGridPoints();

	CGNSWriter writer(&m);
	writer.write("test_monobloc");

}
TEST(BlockingTestSuite, test_monobloc_2){

	Blocking2D m;
	Node n1 = m.newBlockCorner(0,0);
	Node n2 = m.newBlockCorner(1,0);
	Node n3 = m.newBlockCorner(1,1);
	Node n4=  m.newBlockCorner(0,1);

	Blocking2D::Block b0 = m.newBlock(n1,n2,n3,n4);

	b0.setNbDiscretizationI(4);
	b0.setNbDiscretizationJ(4);

	Node n5 = m.newBlockCorner(2,0);
	Node n6 = m.newBlockCorner(2,1.5);

	Blocking2D::Block b1 = m.newBlock(n6,n3,n2,n5);

	b1.setNbDiscretizationI(4);
	b1.setNbDiscretizationJ(4);

	b0.computeT();
	b1.computeT();

	std::cout<<"b0.T = ("<<b0.getT()[1][0]<<","<<b0.getT()[1][1]<<")"<<std::endl;
	std::cout<<"b1.T = ("<<b1.getT()[0][0]<<","<<b1.getT()[0][1]<<")"<<std::endl;

	m.initializeGridPoints();

	CGNSWriter writer(&m);
	writer.write("test_monobloc");

}
/*----------------------------------------------------------------------------*/
TEST(BlockingTestSuite, test_triangles){

	Mesh m(DIM2|N|F|N2F|F2N);

	std::string vtk_file = "/home/calderans/dev/test_cgns_tri.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.updateUpwardConnectivity();

	gmds::IGMeshIOService ioService1(&m);
	gmds::VTKWriter writer1(&ioService1);
	writer1.setCellOptions(N|F);
	writer1.setDataOptions(N|F);
	writer1.write("test_monobloc_tri.vtk");


	CGNSWriter writer(&m);
	writer.write("test_monobloc_tri");

}