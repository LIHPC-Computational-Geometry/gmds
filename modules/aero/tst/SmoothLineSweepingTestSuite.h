//
// Created by rochec on 24/08/2022.
//

#ifndef GMDS_SMOOTHLINESWEEPINGTESTSUITE_H
#define GMDS_SMOOTHLINESWEEPINGTESTSUITE_H

/*----------------------------------------------------------------------------*/
#include "gmds/io/VTKReader.h"
#include <gmds/aero/SmoothLineSweepingYao.h>
#include <gmds/aero/SmoothLineSweepingOrtho.h>
#include <gmds/ig/Mesh.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gtest/gtest.h>
#include <iostream>
#include <unit_test_config.h>
#include <random>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

TEST(AeroTestClass, testGrid_SmoothLineSweepingYao)
{
	Blocking2D m;
	Node n1 = m.newBlockCorner(0,0);
	Node n2 = m.newBlockCorner(1,0);
	Node n3 = m.newBlockCorner(1,1);
	Node n4=  m.newBlockCorner(0,1);

	Blocking2D::Block b1 = m.newBlock(n1,n2,n3,n4);

	b1.setNbDiscretizationI(30);
	b1.setNbDiscretizationJ(30);

	int Nx = b1.getNbDiscretizationI();
	int Ny = b1.getNbDiscretizationJ();

	m.initializeGridPoints();

	//==================================================================
	// PERFORM THE PERTURBATION
	//==================================================================
	constexpr int FLOAT_MIN = 0;
	constexpr int FLOAT_MAX = 1;
	std::random_device rd;
	std::default_random_engine eng(rd());
	std::uniform_real_distribution<float> distr(FLOAT_MIN, FLOAT_MAX);

	b1 = m.block(0);
	for(int i=1; i < Nx-1; i++){
		for (int j=1; j < Ny-1; j++)
		{
			b1(i,j).setXYZ(distr(eng), distr(eng), 0);
		}
	}


	IGMeshIOService ios(&m);
	VTKWriter writer(&ios);
	writer.setCellOptions(N|F);
	writer.setDataOptions(N|F);
	writer.write("testGrid_SmoothLineSweepingYao_Init.vtk");

	SmoothLineSweepingYao smoother( &b1, 1000, 0.0);
	smoother.execute();

	for (int i=0; i<Nx-1; i++)
	{
		for (int j=0; j<Ny-1; j++)
		{
			ASSERT_NEAR(b1(i,j).X(), i*(1.0/(Nx-1.0)), pow(10,-2));
			ASSERT_NEAR(b1(i,j).Y(), j*(1.0/(Ny-1.0)), pow(10,-2));
		}
	}

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("testGrid_SmoothLineSweepingYao_Result.vtk");

}


TEST(AeroTestClass, testGrid_SmoothLineSweepingOrtho)
{
	Blocking2D m;
	Node n1 = m.newBlockCorner(0,0);
	Node n2 = m.newBlockCorner(1,0);
	Node n3 = m.newBlockCorner(1.5,1);
	Node n4=  m.newBlockCorner(0.5,1);

	Blocking2D::Block b1 = m.newBlock(n1,n2,n3,n4);

	b1.setNbDiscretizationI(100);
	b1.setNbDiscretizationJ(100);

	m.initializeGridPoints();

	b1 = m.block(0);

	IGMeshIOService ios(&m);
	VTKWriter writer(&ios);
	writer.setCellOptions(N|F);
	writer.setDataOptions(N|F);
	writer.write("testGrid_SmoothLineSweepingOrtho.vtk");

	SmoothLineSweepingOrtho smoother( &b1, 300, 0.3);
	smoother.execute();

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("testGrid_SmoothLineSweepingOrtho_result.vtk");

}

#endif     // GMDS_SMOOTHLINESWEEPINGTESTSUITE_H
