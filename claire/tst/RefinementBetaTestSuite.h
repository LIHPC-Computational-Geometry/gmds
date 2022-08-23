//
// Created by rochec on 27/07/2022.
//

/*----------------------------------------------------------------------------*/
#include <gmds/cadfac/FACManager.h>
#include "gmds/io/VTKReader.h"
#include <gmds/ig/Blocking2D.h>
#include <gmds/claire/RefinementBeta.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gtest/gtest.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

TEST(RefinementBetaTestClass, test_1)
{
	Blocking2D m;
	Node n1 = m.newBlockCorner(0,0);
	Node n2 = m.newBlockCorner(1,0);
	Node n3 = m.newBlockCorner(1,1);
	Node n4=  m.newBlockCorner(0,1);

	Blocking2D::Block b1 = m.newBlock(n1,n2,n3,n4);

	b1.setNbDiscretizationI(100);
	b1.setNbDiscretizationJ(100);

	m.initializeGridPoints();

	IGMeshIOService ios(&m);
	VTKWriter writer(&ios);
	writer.setCellOptions(N|F);
	writer.setDataOptions(N|F);
	writer.write("Test_RefinementBeta_Init.vtk");

	b1 = m.block(0);
	int Nx = b1.getNbDiscretizationI();
	int Ny = b1.getNbDiscretizationJ();

	for (int i=0; i < Nx; i++)
	{
		std::vector<math::Point> Points;
		for (int j = 0; j < Ny; j++)
		{
			Points.push_back(b1(i, j).point());
		}

		RefinementBeta ref(Points, pow(10, -8));
		ref.execute();
		Points = ref.GetNewPositions();

		for (int j = 1; j < Ny - 1; j++)
		{
			b1(i, j).setPoint({Points[j]});
		}

	}


	for (int i=0; i < Nx; i++)
	{
		double size_first_edge = (b1(i,1).point() - b1(i,0).point()).norm();
		ASSERT_NEAR(size_first_edge, pow(10,-8), pow(10,-8));
	}


	{
		ASSERT_NEAR(b1(0,0).Y(), 0, pow(10,-8));
		ASSERT_NEAR(b1(0,1).Y(), 9.99919e-09, pow(10,-8));
		ASSERT_NEAR(b1(0,2).Y(), 2.1928e-08, pow(10,-8));
		ASSERT_NEAR(b1(0,3).Y(), 3.61588e-08, pow(10,-8));
		ASSERT_NEAR(b1(0,4).Y(), 5.31358e-08, pow(10,-8));
		ASSERT_NEAR(b1(0,5).Y(), 7.3389e-08, pow(10,-8));
		ASSERT_NEAR(b1(0,6).Y(), 9.75505e-08, pow(10,-8));
		ASSERT_NEAR(b1(0,7).Y(), 1.26375e-07, pow(10,-8));
		ASSERT_NEAR(b1(0,8).Y(), 1.60761e-07, pow(10,-8));
		ASSERT_NEAR(b1(0,9).Y(), 2.01784e-07, pow(10,-8));
		ASSERT_NEAR(b1(0,10).Y(), 2.50722e-07, pow(10,-8));
		ASSERT_NEAR(b1(0,11).Y(), 3.09105e-07, pow(10,-8));
		ASSERT_NEAR(b1(0,12).Y(), 3.78755e-07, pow(10,-8));
		ASSERT_NEAR(b1(0,13).Y(), 4.61845e-07, pow(10,-8));
		ASSERT_NEAR(b1(0,14).Y(), 5.60969e-07, pow(10,-8));
		ASSERT_NEAR(b1(0,15).Y(), 6.79222e-07, pow(10,-8));
		ASSERT_NEAR(b1(0,16).Y(), 8.20296e-07, pow(10,-8));
		ASSERT_NEAR(b1(0,17).Y(), 9.88593e-07, pow(10,-8));
		ASSERT_NEAR(b1(0,18).Y(), 1.18937e-06, pow(10,-8));
		ASSERT_NEAR(b1(0,19).Y(), 1.42889e-06, pow(10,-8));
		ASSERT_NEAR(b1(0,20).Y(), 1.71463e-06, pow(10,-8));
		ASSERT_NEAR(b1(0,21).Y(), 2.05551e-06, pow(10,-8));
		ASSERT_NEAR(b1(0,22).Y(), 2.46218e-06, pow(10,-8));
		ASSERT_NEAR(b1(0,23).Y(), 2.94732e-06, pow(10,-8));
		ASSERT_NEAR(b1(0,24).Y(), 3.52608e-06, pow(10,-8));
		ASSERT_NEAR(b1(0,25).Y(), 4.21653e-06, pow(10,-8));
		ASSERT_NEAR(b1(0,26).Y(), 5.04022e-06, pow(10,-8));
		ASSERT_NEAR(b1(0,27).Y(), 6.02286e-06, pow(10,-8));
		ASSERT_NEAR(b1(0,28).Y(), 7.19513e-06, pow(10,-8));
		ASSERT_NEAR(b1(0,29).Y(), 8.59361e-06, pow(10,-8));
		ASSERT_NEAR(b1(0,30).Y(), 1.0262e-05, pow(10,-8));
		ASSERT_NEAR(b1(0,31).Y(), 1.22523e-05, pow(10,-8));
		ASSERT_NEAR(b1(0,32).Y(), 1.46267e-05, pow(10,-8));
		ASSERT_NEAR(b1(0,33).Y(), 1.74593e-05, pow(10,-8));
		ASSERT_NEAR(b1(0,34).Y(), 2.08384e-05, pow(10,-8));
		ASSERT_NEAR(b1(0,35).Y(), 2.48697e-05, pow(10,-8));
		ASSERT_NEAR(b1(0,36).Y(), 2.96789e-05, pow(10,-8));
		ASSERT_NEAR(b1(0,37).Y(), 3.54162e-05, pow(10,-8));
		ASSERT_NEAR(b1(0,38).Y(), 4.22606e-05, pow(10,-8));
		ASSERT_NEAR(b1(0,39).Y(), 5.04257e-05, pow(10,-8));
		ASSERT_NEAR(b1(0,40).Y(), 6.01664e-05, pow(10,-8));
		ASSERT_NEAR(b1(0,41).Y(), 7.17867e-05, pow(10,-8));
		ASSERT_NEAR(b1(0,42).Y(), 8.56492e-05, pow(10,-8));
		ASSERT_NEAR(b1(0,43).Y(), 0.000102187, pow(10,-8));
		ASSERT_NEAR(b1(0,44).Y(), 0.000121915, pow(10,-8));
		ASSERT_NEAR(b1(0,45).Y(), 0.00014545, pow(10,-8));
		ASSERT_NEAR(b1(0,46).Y(), 0.000173526, pow(10,-8));
		ASSERT_NEAR(b1(0,47).Y(), 0.000207019, pow(10,-8));
		ASSERT_NEAR(b1(0,48).Y(), 0.000246974, pow(10,-8));
		ASSERT_NEAR(b1(0,49).Y(), 0.000294638, pow(10,-8));
		ASSERT_NEAR(b1(0,50).Y(), 0.000351496, pow(10,-8));
		ASSERT_NEAR(b1(0,51).Y(), 0.000419322, pow(10,-8));
		ASSERT_NEAR(b1(0,52).Y(), 0.000500232, pow(10,-8));
		ASSERT_NEAR(b1(0,53).Y(), 0.000596746, pow(10,-8));
		ASSERT_NEAR(b1(0,54).Y(), 0.000711873, pow(10,-8));
		ASSERT_NEAR(b1(0,55).Y(), 0.0008492, pow(10,-8));
		ASSERT_NEAR(b1(0,56).Y(), 0.001013, pow(10,-8));
		ASSERT_NEAR(b1(0,57).Y(), 0.00120838, pow(10,-8));
		ASSERT_NEAR(b1(0,58).Y(), 0.00144141, pow(10,-8));
		ASSERT_NEAR(b1(0,59).Y(), 0.00171934, pow(10,-8));
		ASSERT_NEAR(b1(0,60).Y(), 0.00205081, pow(10,-8));
		ASSERT_NEAR(b1(0,61).Y(), 0.00244609, pow(10,-8));
		ASSERT_NEAR(b1(0,62).Y(), 0.00291745, pow(10,-8));
		ASSERT_NEAR(b1(0,63).Y(), 0.00347948, pow(10,-8));
		ASSERT_NEAR(b1(0,64).Y(), 0.00414956, pow(10,-8));
		ASSERT_NEAR(b1(0,65).Y(), 0.00494835, pow(10,-8));
		ASSERT_NEAR(b1(0,66).Y(), 0.00590047, pow(10,-8));
		ASSERT_NEAR(b1(0,67).Y(), 0.00703512, pow(10,-8));
		ASSERT_NEAR(b1(0,68).Y(), 0.00838706, pow(10,-8));
		ASSERT_NEAR(b1(0,69).Y(), 0.00999748, pow(10,-8));
	}

	{
		ASSERT_NEAR(b1(0,0).X(), 0, pow(10,-6));
		ASSERT_NEAR(b1(1,1).X(), 0.010101, pow(10,-6));
		ASSERT_NEAR(b1(2,2).X(), 0.020202, pow(10,-6));
		ASSERT_NEAR(b1(3,3).X(), 0.030303, pow(10,-6));
		ASSERT_NEAR(b1(4,4).X(), 0.040404, pow(10,-6));
		ASSERT_NEAR(b1(5,5).X(), 0.0505051, pow(10,-6));
		ASSERT_NEAR(b1(6,6).X(), 0.0606061, pow(10,-6));
		ASSERT_NEAR(b1(7,7).X(), 0.0707071, pow(10,-6));
		ASSERT_NEAR(b1(8,8).X(), 0.0808081, pow(10,-6));
		ASSERT_NEAR(b1(9,9).X(), 0.0909091, pow(10,-6));
		ASSERT_NEAR(b1(10,10).X(), 0.10101, pow(10,-6));
		ASSERT_NEAR(b1(11,11).X(), 0.111111, pow(10,-6));
		ASSERT_NEAR(b1(12,12).X(), 0.121212, pow(10,-6));
		ASSERT_NEAR(b1(13,13).X(), 0.131313, pow(10,-6));
		ASSERT_NEAR(b1(14,14).X(), 0.141414, pow(10,-6));
		ASSERT_NEAR(b1(15,15).X(), 0.151515, pow(10,-6));
		ASSERT_NEAR(b1(16,16).X(), 0.161616, pow(10,-6));
		ASSERT_NEAR(b1(17,17).X(), 0.171717, pow(10,-6));
		ASSERT_NEAR(b1(18,18).X(), 0.181818, pow(10,-6));
		ASSERT_NEAR(b1(19,19).X(), 0.191919, pow(10,-6));
		ASSERT_NEAR(b1(20,20).X(), 0.20202, pow(10,-6));
		ASSERT_NEAR(b1(21,21).X(), 0.212121, pow(10,-6));
		ASSERT_NEAR(b1(22,22).X(), 0.222222, pow(10,-6));
		ASSERT_NEAR(b1(23,23).X(), 0.232323, pow(10,-6));
		ASSERT_NEAR(b1(24,24).X(), 0.242424, pow(10,-6));
		ASSERT_NEAR(b1(25,25).X(), 0.252525, pow(10,-6));
		ASSERT_NEAR(b1(26,26).X(), 0.262626, pow(10,-6));
		ASSERT_NEAR(b1(27,27).X(), 0.272727, pow(10,-6));
		ASSERT_NEAR(b1(28,28).X(), 0.282828, pow(10,-6));
		ASSERT_NEAR(b1(29,29).X(), 0.292929, pow(10,-6));
		ASSERT_NEAR(b1(30,30).X(), 0.30303, pow(10,-6));
		ASSERT_NEAR(b1(31,31).X(), 0.313131, pow(10,-6));
		ASSERT_NEAR(b1(32,32).X(), 0.323232, pow(10,-6));
		ASSERT_NEAR(b1(33,33).X(), 0.333333, pow(10,-6));
		ASSERT_NEAR(b1(34,34).X(), 0.343434, pow(10,-6));
		ASSERT_NEAR(b1(35,35).X(), 0.353535, pow(10,-6));
		ASSERT_NEAR(b1(36,36).X(), 0.363636, pow(10,-6));
		ASSERT_NEAR(b1(37,37).X(), 0.373737, pow(10,-6));
		ASSERT_NEAR(b1(38,38).X(), 0.383838, pow(10,-6));
		ASSERT_NEAR(b1(39,39).X(), 0.393939, pow(10,-6));
		ASSERT_NEAR(b1(40,40).X(), 0.40404, pow(10,-6));
		ASSERT_NEAR(b1(41,41).X(), 0.414141, pow(10,-6));
		ASSERT_NEAR(b1(42,42).X(), 0.424242, pow(10,-6));
		ASSERT_NEAR(b1(43,43).X(), 0.434343, pow(10,-6));
		ASSERT_NEAR(b1(44,44).X(), 0.444444, pow(10,-6));
		ASSERT_NEAR(b1(45,45).X(), 0.454545, pow(10,-6));
		ASSERT_NEAR(b1(46,46).X(), 0.464646, pow(10,-6));
		ASSERT_NEAR(b1(47,47).X(), 0.474747, pow(10,-6));
		ASSERT_NEAR(b1(48,48).X(), 0.484848, pow(10,-6));
		ASSERT_NEAR(b1(49,49).X(), 0.494949, pow(10,-6));
		ASSERT_NEAR(b1(50,50).X(), 0.505051, pow(10,-6));
		ASSERT_NEAR(b1(51,51).X(), 0.515152, pow(10,-6));
		ASSERT_NEAR(b1(52,52).X(), 0.525253, pow(10,-6));
		ASSERT_NEAR(b1(53,53).X(), 0.535354, pow(10,-6));
		ASSERT_NEAR(b1(54,54).X(), 0.545455, pow(10,-6));
		ASSERT_NEAR(b1(55,55).X(), 0.555556, pow(10,-6));
		ASSERT_NEAR(b1(56,56).X(), 0.565657, pow(10,-6));
		ASSERT_NEAR(b1(57,57).X(), 0.575758, pow(10,-6));
		ASSERT_NEAR(b1(58,58).X(), 0.585859, pow(10,-6));
		ASSERT_NEAR(b1(59,59).X(), 0.59596, pow(10,-6));
		ASSERT_NEAR(b1(60,60).X(), 0.606061, pow(10,-6));
		ASSERT_NEAR(b1(61,61).X(), 0.616162, pow(10,-6));
		ASSERT_NEAR(b1(62,62).X(), 0.626263, pow(10,-6));
		ASSERT_NEAR(b1(63,63).X(), 0.636364, pow(10,-6));
		ASSERT_NEAR(b1(64,64).X(), 0.646465, pow(10,-6));
		ASSERT_NEAR(b1(65,65).X(), 0.656566, pow(10,-6));
		ASSERT_NEAR(b1(66,66).X(), 0.666667, pow(10,-6));
		ASSERT_NEAR(b1(67,67).X(), 0.676768, pow(10,-6));
		ASSERT_NEAR(b1(68,68).X(), 0.686869, pow(10,-6));
		ASSERT_NEAR(b1(69,69).X(), 0.69697, pow(10,-6));
		ASSERT_NEAR(b1(70,70).X(), 0.707071, pow(10,-6));
		ASSERT_NEAR(b1(71,71).X(), 0.717172, pow(10,-6));
		ASSERT_NEAR(b1(72,72).X(), 0.727273, pow(10,-6));
		ASSERT_NEAR(b1(73,73).X(), 0.737374, pow(10,-6));
		ASSERT_NEAR(b1(74,74).X(), 0.747475, pow(10,-6));
		ASSERT_NEAR(b1(75,75).X(), 0.757576, pow(10,-6));
		ASSERT_NEAR(b1(76,76).X(), 0.767677, pow(10,-6));
		ASSERT_NEAR(b1(77,77).X(), 0.777778, pow(10,-6));
		ASSERT_NEAR(b1(78,78).X(), 0.787879, pow(10,-6));
		ASSERT_NEAR(b1(79,79).X(), 0.79798, pow(10,-6));
		ASSERT_NEAR(b1(80,80).X(), 0.808081, pow(10,-6));
		ASSERT_NEAR(b1(81,81).X(), 0.818182, pow(10,-6));
		ASSERT_NEAR(b1(82,82).X(), 0.828283, pow(10,-6));
		ASSERT_NEAR(b1(83,83).X(), 0.838384, pow(10,-6));
		ASSERT_NEAR(b1(84,84).X(), 0.848485, pow(10,-6));
		ASSERT_NEAR(b1(85,85).X(), 0.858586, pow(10,-6));
		ASSERT_NEAR(b1(86,86).X(), 0.868687, pow(10,-6));
		ASSERT_NEAR(b1(87,87).X(), 0.878788, pow(10,-6));
		ASSERT_NEAR(b1(88,88).X(), 0.888889, pow(10,-6));
		ASSERT_NEAR(b1(89,89).X(), 0.89899, pow(10,-6));
		ASSERT_NEAR(b1(90,90).X(), 0.909091, pow(10,-6));
		ASSERT_NEAR(b1(91,91).X(), 0.919192, pow(10,-6));
		ASSERT_NEAR(b1(92,92).X(), 0.929293, pow(10,-6));
		ASSERT_NEAR(b1(93,93).X(), 0.939394, pow(10,-6));
		ASSERT_NEAR(b1(94,94).X(), 0.949495, pow(10,-6));
		ASSERT_NEAR(b1(95,95).X(), 0.959596, pow(10,-6));
		ASSERT_NEAR(b1(96,96).X(), 0.969697, pow(10,-6));
		ASSERT_NEAR(b1(97,97).X(), 0.979798, pow(10,-6));
		ASSERT_NEAR(b1(98,98).X(), 0.989899, pow(10,-6));
		ASSERT_NEAR(b1(99,99).X(), 1, pow(10,-6));
	}




	/*
	for (int i=0; i < Nx;i++)
	{
		for (int j=0; j < Ny; j++)
		{
			if (i == j) {
				std::cout << "ASSERT_NEAR(b1(" << i << "," << j << ").X(), " << b1(i, j).X() << ", pow(10,-6));" << std::endl;
			}
		}
	}
	*/

	/*
	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("Test_RefinementBeta_Result.vtk");
	*/

}
