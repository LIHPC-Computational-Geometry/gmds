//
// Created by rochec on 01/12/23.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/aero/TransfiniteInterpolation_3D.h>
#include <gmds/math/TransfiniteInterpolation.h>
#include <iostream>
#include <fstream>
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
TEST(TransfiniteInterpolation_3DTest, testUnit_Transfinite2D_imgManuscrit)
{
	int nb_i(10);
	int nb_j(20);
	Array2D<math::Point> pnts(nb_i,nb_j);

	double di = 1.0/(double(nb_i)-1.0) ;
	double dj = 1.0/(double(nb_j)-1.0) ;

	for (auto i=0;i<nb_i;i++)
	{
		double y0 = sin(M_PI*(i*di))/10.0 ;
		double y1 = sin(2.0*M_PI*(i*di))/12.0 ;
		pnts(i,0).setXYZ(i*di, y0, 0.0) ;
		pnts(i,nb_j-1).setXYZ(i*di, 1+y1, 0.0) ;
	}

	for (auto j=0;j<nb_j;j++)
	{
		double x0 = sin(M_PI*(j*dj))/6.0 ;
		double x1 = sin(2.0*M_PI*(j*dj))/8.0 ;
		pnts(0,j).setXYZ(x0, j*dj, 0.0) ;
		pnts(nb_i-1,j).setXYZ(1+x1, j*dj, 0.0) ;
	}

	// First, we create the file where we are going to store the info
	//std::ofstream stream= std::ofstream("test.txt", std::ios::out);
	//set the numerical precision (number of digits)
	//stream.precision(15);

	Mesh m_bnd = Mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                        F2E | E2F | R2E | E2R | N2R | N2F | N2E));

	std::ofstream stream_bnd= std::ofstream("TFI2D_Bnd.table", std::ios::out);
	stream_bnd.precision(15);
	for (auto i=0;i<nb_i;i++)
	{
		for (auto j=0;j<nb_j;j++)
		{
			Node n = m_bnd.newNode(pnts(i,j));
			//stream << pnts(i,j).X() << " " << pnts(i,j).Y() << "\n";
			stream_bnd << n.X() << " " << n.Y() << "\n";
		}
	}
	stream_bnd.close();

	std::ofstream stream_param= std::ofstream("TFI2D_ParamSpace.table", std::ios::out);
	stream_param.precision(15);
	for (auto i=0;i<nb_i-1;i++)
	{
		for (auto j=0;j<nb_j-1;j++)
		{
			stream_param << i*di << " " << j*dj << "\n";
			stream_param << (i+1)*di << " " << j*dj << "\n";
			stream_param << (i+1)*di << " " << (j+1)*dj << "\n";
			stream_param << i*di << " " << (j+1)*dj << "\n";
			stream_param << i*di << " " << j*dj << "\n";
			stream_param << "\n";
		}
	}
	stream_param.close();

	gmds::IGMeshIOService ioService_bnd(&m_bnd);
	gmds::VTKWriter vtkWriter_bnd(&ioService_bnd);
	vtkWriter_bnd.setCellOptions(gmds::N|gmds::F);
	vtkWriter_bnd.setDataOptions(gmds::N|gmds::F);
	vtkWriter_bnd.write("TransfiniteInterpolation_2DTestSuite_Bnd.vtk");


	EXPECT_TRUE(math::TransfiniteInterpolation::computeQuad(pnts));

	std::ofstream stream_physic= std::ofstream("TFI2D_PhysicalSpace.table", std::ios::out);
	stream_physic.precision(15);
	for (auto i=0;i<nb_i-1;i++)
	{
		for (auto j=0;j<nb_j-1;j++)
		{
			stream_physic << pnts(i,j).X() << " " << pnts(i,j).Y() << "\n";
			stream_physic << pnts(i+1,j).X() << " " << pnts(i+1,j).Y() << "\n";
			stream_physic << pnts(i+1,j+1).X() << " " << pnts(i+1,j+1).Y() << "\n";
			stream_physic << pnts(i,j+1).X() << " " << pnts(i,j+1).Y() << "\n";
			stream_physic << pnts(i,j).X() << " " << pnts(i,j).Y() << "\n";
			stream_physic << "\n";
		}
	}
	stream_physic.close();

	/*
	ASSERT_FLOAT_EQ(pnts(0,0).X(), 0.0);
	ASSERT_FLOAT_EQ(pnts(0,0).Y(), 0.0);
	ASSERT_FLOAT_EQ(pnts(0,0).Z(), 0.0);

	ASSERT_FLOAT_EQ(pnts(nb_i-1,nb_j-1).X(), 1.0);
	ASSERT_FLOAT_EQ(pnts(nb_i-1,nb_j-1).Y(), 1.0);
	ASSERT_FLOAT_EQ(pnts(nb_i-1,nb_j-1).Z(), 0.0);

	ASSERT_FLOAT_EQ(pnts(1,1).X(), di);
	ASSERT_FLOAT_EQ(pnts(1,1).Y(), dj);
	 */

	Mesh m = Mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                        F2E | E2F | R2E | E2R | N2R | N2F | N2E));

	Array2D<TCellID> n_ids(nb_i,nb_j);
	for (auto i=0;i<nb_i;i++)
	{
		for (auto j=0;j<nb_j;j++)
		{
			Node n = m.newNode(pnts(i,j));
			n_ids(i,j) = n.id() ;
		}
	}

	for (auto i=0;i<nb_i-1;i++)
	{
		for (auto j=0;j<nb_j-1;j++)
		{
			m.newQuad(n_ids(i,j), n_ids(i+1,j), n_ids(i+1,j+1), n_ids(i,j+1));
		}
	}

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("TransfiniteInterpolation_2DTestSuite.vtk");

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
TEST(TransfiniteInterpolation_3DTest, testUnit_Transfinite3D_imgManuscrit)
{
	int nb_i(10);
	int nb_j(15);
	int nb_k(20);
	Array3D<math::Point> pnts(nb_i,nb_j,nb_k);

	double di = 1.0/(double(nb_i)-1.0) ;
	double dj = 1.0/(double(nb_j)-1.0) ;
	double dk = 1.0/(double(nb_k)-1.0) ;

	for (auto i=0;i<nb_i;i++)
	{
		for (auto j=0;j<nb_j;j++)
		{
			double z0 = sin(M_PI*(i*di))/10.0 ;
			double z1 = sin(M_PI*(i*di))/12.0 ;
			pnts(i,j,0).setXYZ(i*di, j*dj, 0.0+z0) ;
			pnts(i,j,nb_k-1).setXYZ(i*di, j*dj, 1.0+z1) ;
		}
	}

	for (auto i=0;i<nb_i;i++)
	{
		for (auto k=0;k<nb_k;k++)
		{
			double y0 = sin(M_PI*(k*dk*i*di))/10.0 ;
			double y1 = sin(M_PI*(i*di*k*dk))/25.0 ;
			pnts(i,0,k).setXYZ(i*di, 0.0+y0, k*dk) ;
			pnts(i,nb_j-1,k).setXYZ(i*di, 1.0+y1, k*dk) ;
		}
	}

	for (auto j=0;j<nb_j;j++)
	{
		for (auto k=0;k<nb_k;k++)
		{
			double x0 = sin(M_PI*(j*dj))/10.0 ;
			double x1 = sin(M_PI*(j*dj))/12.0 ;
			pnts(0,j,k).setXYZ(0.0+x0, j*dj, k*dk) ;
			pnts(nb_i-1,j,k).setXYZ(1.0+x1, j*dj, k*dk) ;
		}
	}


	EXPECT_TRUE(TransfiniteInterpolation_3D::computeHex(pnts));

	/*
	ASSERT_FLOAT_EQ(pnts(0,0,0).X(), 0.0);
	ASSERT_FLOAT_EQ(pnts(0,0,0).Y(), 0.0);
	ASSERT_FLOAT_EQ(pnts(0,0,0).Z(), 0.0);

	ASSERT_FLOAT_EQ(pnts(nb_i-1,nb_j-1,nb_k-1).X(), 1.0);
	ASSERT_FLOAT_EQ(pnts(nb_i-1,nb_j-1,nb_k-1).Y(), 1.0);
	ASSERT_FLOAT_EQ(pnts(nb_i-1,nb_j-1,nb_k-1).Z(), 1.0);

	ASSERT_FLOAT_EQ(pnts(1,1,1).X(), di);
	ASSERT_FLOAT_EQ(pnts(1,1,1).Y(), dj);
	ASSERT_FLOAT_EQ(pnts(1,1,1).Z(), dk);
	 */

	Mesh m_bnd = Mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                            F2E | E2F | R2E | E2R | N2R | N2F | N2E));

	for (auto i=0;i<nb_i;i++)
	{
		for (auto j=0;j<nb_j;j++)
		{
			for (auto k=0;k<nb_k;k++)
			{
				m_bnd.newNode(pnts(0,j,k));
				m_bnd.newNode(pnts(nb_i-1,j,k));
				m_bnd.newNode(pnts(i,0,k));
				m_bnd.newNode(pnts(i,nb_j-1,k));
				m_bnd.newNode(pnts(i,j,0));
				m_bnd.newNode(pnts(i,j,nb_k-1));
			}
		}
	}

	gmds::IGMeshIOService ioService_bnd(&m_bnd);
	gmds::VTKWriter vtkWriter_bnd(&ioService_bnd);
	vtkWriter_bnd.setCellOptions(gmds::N|gmds::R);
	vtkWriter_bnd.setDataOptions(gmds::N|gmds::R);
	vtkWriter_bnd.write("TransfiniteInterpolation_3DTestSuite_Bnd.vtk");

	Mesh m = Mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                        F2E | E2F | R2E | E2R | N2R | N2F | N2E));

	Array3D<TCellID> n_ids(nb_i,nb_j,nb_k);
	for (auto i=0;i<nb_i;i++)
	{
		for (auto j=0;j<nb_j;j++)
		{
			for (auto k=0;k<nb_k;k++)
			{
				Node n = m.newNode(pnts(i,j,k));
				n_ids(i,j,k) = n.id();
			}
		}
	}

	for (auto i=0;i<nb_i-1;i++)
	{
		for (auto j=0;j<nb_j-1;j++)
		{
			for (auto k=0;k<nb_k-1;k++)
			{
				m.newHex(n_ids(i,j,k), n_ids(i+1,j,k), n_ids(i+1,j+1,k), n_ids(i,j+1,k), n_ids(i,j,k+1), n_ids(i+1,j,k+1), n_ids(i+1,j+1,k+1), n_ids(i,j+1,k+1));
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