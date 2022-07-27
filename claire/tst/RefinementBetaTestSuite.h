//
// Created by rochec on 27/07/2022.
//

/*----------------------------------------------------------------------------*/
#include <gmds/cad/FACManager.h>
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

	std::cout << "Size 1 : " << Points.size() << std::endl;

	RefinementBeta ref(Points, pow(10, -8));
	ref.execute();
	Points = ref.GetNewPositions();
	std::cout << "Size 2 : " << Points.size() << std::endl;

	for (int j = 1; j < Ny - 1; j++)
	{
		b1(i, j).setPoint({Points[j]});
	}

	}

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("Test_RefinementBeta_Result.vtk");

}
