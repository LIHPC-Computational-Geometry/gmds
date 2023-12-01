/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/claire/Blocking3D.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(Blocking3DTestSuite, test_blocking3D_1)
{
	Blocking3D b;
	Node n1 = b.newBlockCorner(0,0, 0);
	Node n2 = b.newBlockCorner(1,0, 0);
	Node n3 = b.newBlockCorner(1,1, 0);
	Node n4=  b.newBlockCorner(0,1, 0);
	Node n5=  b.newBlockCorner(0,0, 1);
	Node n6=  b.newBlockCorner(1,0, 1);
	Node n7=  b.newBlockCorner(1,1, 1);
	Node n8=  b.newBlockCorner(0,1, 1);

	Blocking3D::Block b0 = b.newBlock(n1,n2,n3,n4,n5,n6,n7,n8);

	Node n9 = b.newBlockCorner(2,0,0);
	Node n10 = b.newBlockCorner(2,0,1);
	Node n11 = b.newBlockCorner(2,1,0);
	Node n12 = b.newBlockCorner(2,1,1);
	Blocking3D::Block b1 = b.newBlock(n2,n9,n11,n3,n6,n10,n12,n7);

	ASSERT_EQ(b0.id(), b.block(0).id());
	int nb_I(5);
	int nb_J(5);
	int nb_K(11);
	b0.setNbDiscretizationI(nb_I);
	b0.setNbDiscretizationJ(nb_J);
	b0.setNbDiscretizationK(nb_K);
	b1.setNbDiscretizationI(nb_I);
	b1.setNbDiscretizationJ(nb_J);
	b1.setNbDiscretizationK(nb_K);

	ASSERT_EQ(n1.id(),b0.origin());
	ASSERT_FLOAT_EQ(b0.getUnitVectorI().X(),1.0);
	ASSERT_FLOAT_EQ(b0.getUnitVectorI().Y(),0.0);
	ASSERT_FLOAT_EQ(b0.getUnitVectorI().Z(),0.0);
	ASSERT_FLOAT_EQ(b0.getUnitVectorJ().X(),0.0);
	ASSERT_FLOAT_EQ(b0.getUnitVectorJ().Y(),1.0);
	ASSERT_FLOAT_EQ(b0.getUnitVectorJ().Z(),0.0);
	ASSERT_FLOAT_EQ(b0.getUnitVectorK().X(),0.0);
	ASSERT_FLOAT_EQ(b0.getUnitVectorK().Y(),0.0);
	ASSERT_FLOAT_EQ(b0.getUnitVectorK().Z(),1.0);

	ASSERT_EQ(b0.getNbDiscretizationI(),nb_I);
	ASSERT_EQ(b0.getNbDiscretizationJ(),nb_J);
	ASSERT_EQ(b0.getNbDiscretizationK(),nb_K);

	b.initializeGridPoints();

	b0 = b.block(0);
	ASSERT_FLOAT_EQ(b0(1,1,1).X(), 0.25);
	ASSERT_FLOAT_EQ(b0(1,1,1).Y(), 0.25);
	ASSERT_FLOAT_EQ(b0(1,1,1).Z(), 0.1);

	ASSERT_FLOAT_EQ(b0(3,2,8).X(), 0.75);
	ASSERT_FLOAT_EQ(b0(3,2,8).Y(), 0.5);
	ASSERT_FLOAT_EQ(b0(3,2,8).Z(), 0.8);


	gmds::IGMeshIOService ioService(&b);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("Blocking3DTestSuite.vtk");

}
/*----------------------------------------------------------------------------*/