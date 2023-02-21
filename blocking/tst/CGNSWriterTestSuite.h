/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include "gmds/ig/Blocking2D.h"
#include "gmds/ig/MeshDoctor.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKReader.h"
#include "gmds/io/VTKWriter.h"
#include <gmds/blocking/CGNSWriter.h>


#include <cgnslib.h>

/*----------------------------------------------------------------------------*/
//using namespace gmds;
/*----------------------------------------------------------------------------*/

TEST(CGNSWriterTestSuite, test_monobloc){

	gmds::Blocking2D m;
	gmds::Node n1 = m.newBlockCorner(0,0);
	gmds::Node n2 = m.newBlockCorner(1,0);
	gmds::Node n3 = m.newBlockCorner(1,1);
	gmds::Node n4=  m.newBlockCorner(0,1);

	gmds::Blocking2D::Block b0 = m.newBlock(n1,n2,n3,n4);

	b0.setNbDiscretizationI(11);
	b0.setNbDiscretizationJ(11);

	m.initializeGridPoints();

	gmds::blocking::CGNSWriter writer(&m);
	writer.write("test_monobloc");

}
TEST(CGNSWriterTestSuite, test_monobloc_2){

	gmds::Blocking2D m;
	gmds::Node n1 = m.newBlockCorner(0,0);
	gmds::Node n2 = m.newBlockCorner(1,0);
	gmds::Node n3 = m.newBlockCorner(1,1);
	gmds::Node n4=  m.newBlockCorner(0,1);

	gmds::Blocking2D::Block b0 = m.newBlock(n1,n2,n3,n4);

	b0.setNbDiscretizationI(4);
	b0.setNbDiscretizationJ(4);

	gmds::Node n5 = m.newBlockCorner(2,0);
	gmds::Node n6 = m.newBlockCorner(2,1.5);

	gmds::Blocking2D::Block b1 = m.newBlock(n6,n3,n2,n5);

	b1.setNbDiscretizationI(4);
	b1.setNbDiscretizationJ(4);

	b0.computeT();
	b1.computeT();

	std::cout<<"b0.T = ("<<b0.getT()[1][0]<<","<<b0.getT()[1][1]<<")"<<std::endl;
	std::cout<<"b1.T = ("<<b1.getT()[0][0]<<","<<b1.getT()[0][1]<<")"<<std::endl;

	m.initializeGridPoints();

	gmds::blocking::CGNSWriter writer(&m);
	writer.write("test_monobloc");

}
/*----------------------------------------------------------------------------*/
TEST(CGNSWriterTestSuite, DISABLED_test_triangles){

	gmds::Mesh m(gmds::DIM2|gmds::N|gmds::F|gmds::N2F|gmds::F2N);

	std::string vtk_file = "/home/calderans/dev/test_cgns_tri.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.updateUpwardConnectivity();

	gmds::IGMeshIOService ioService1(&m);
	gmds::VTKWriter writer1(&ioService1);
	writer1.setCellOptions(gmds::N|gmds::F);
	writer1.setDataOptions(gmds::N|gmds::F);
	writer1.write("test_monobloc_tri.vtk");

	gmds::blocking::CGNSWriter writer(&m);
	writer.write("test_monobloc_tri");

}
/*----------------------------------------------------------------------------*/