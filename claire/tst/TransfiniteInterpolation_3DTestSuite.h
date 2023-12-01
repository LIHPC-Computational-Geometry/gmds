//
// Created by rochec on 01/12/23.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/claire/TransfiniteInterpolation_3D.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(TransfiniteInterpolation_3DTest, testUnitHex)
{
	int nb_i(5);
	int nb_j(7);
	int nb_k(9);
	Array3D<math::Point> pnts(nb_i,nb_j,nb_k);

	double di = 1.0/(double(nb_i)-1.0) ;
	double dj = 1.0/(double(nb_j)-1.0) ;
	double dk = 1.0/(double(nb_k)-1.0) ;

	for (auto i=0;i<nb_i;i++)
	{
		for (auto j=0;j<nb_j;j++)
		{
			pnts(i,j,0).setXYZ(i*di, j*dj, 0.0) ;
			pnts(i,j,nb_k-1).setXYZ(i*di, j*dj, 1.0) ;
		}
	}

	for (auto i=0;i<nb_i;i++)
	{
		for (auto k=0;k<nb_k;k++)
		{
			pnts(i,0,k).setXYZ(i*di, 0.0, k*dk) ;
			pnts(i,nb_j-1,k).setXYZ(i*di, 1.0, k*dk) ;
		}
	}

	for (auto j=0;j<nb_j;j++)
	{
		for (auto k=0;k<nb_k;k++)
		{
			pnts(0,j,k).setXYZ(0.0, j*dj, k*dk) ;
			pnts(nb_i-1,j,k).setXYZ(1.0, j*dj, k*dk) ;
		}
	}

	EXPECT_TRUE(TransfiniteInterpolation_3D::computeHex(pnts));

	ASSERT_FLOAT_EQ(pnts(0,0,0).X(), 0.0);
	ASSERT_FLOAT_EQ(pnts(0,0,0).Y(), 0.0);
	ASSERT_FLOAT_EQ(pnts(0,0,0).Z(), 0.0);

	ASSERT_FLOAT_EQ(pnts(nb_i-1,nb_j-1,nb_k-1).X(), 1.0);
	ASSERT_FLOAT_EQ(pnts(nb_i-1,nb_j-1,nb_k-1).Y(), 1.0);
	ASSERT_FLOAT_EQ(pnts(nb_i-1,nb_j-1,nb_k-1).Z(), 1.0);

	ASSERT_FLOAT_EQ(pnts(1,1,1).X(), di);
	ASSERT_FLOAT_EQ(pnts(1,1,1).Y(), dj);
	ASSERT_FLOAT_EQ(pnts(1,1,1).Z(), dk);

	Mesh m = Mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	F2E | E2F | R2E | E2R | N2R | N2F | N2E));

	for (auto i=0;i<nb_i;i++)
	{
		for (auto j=0;j<nb_j;j++)
		{
			for (auto k=0;k<nb_k;k++)
			{
				m.newNode(pnts(i,j,k));
			}
		}
	}

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("TransfiniteInterpolation_3DTestSuite.vtk");

	/*

	EXPECT_TRUE(math::TransfiniteInterpolation::compute(pnts));

	EXPECT_DOUBLE_EQ(pnts[1][1].X(),0.25);
	EXPECT_DOUBLE_EQ(pnts[1][1].Y(),0.25);
	EXPECT_DOUBLE_EQ(pnts[1][1].Z(),0.00);

	EXPECT_DOUBLE_EQ(pnts[2][2].X(),0.5);
	EXPECT_DOUBLE_EQ(pnts[2][2].Y(),0.5);
	EXPECT_DOUBLE_EQ(pnts[2][2].Z(),0.00);
	 */
}
/*----------------------------------------------------------------------------*/