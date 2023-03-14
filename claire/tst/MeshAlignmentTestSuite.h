//
// Created by rochec on 13/03/2023.
//

#include <gmds/claire/MeshAlignment_2D.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gtest/gtest.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

TEST(ClaireTestClass, MeshAlignment_2D_Test)
{
	// Test
	gmds::Mesh m_meshTri(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));

	// Create the triangles of the triangular mesh
	{
		Node n_0 = m_meshTri.newNode({0, 0, 0});
		Node n_1 = m_meshTri.newNode({1, 0, 0});
		Node n_2 = m_meshTri.newNode({2, 0, 0});
		Node n_3 = m_meshTri.newNode({0.5, 0.5, 0});
		Node n_4 = m_meshTri.newNode({1.5, 0.5, 0});
		Node n_5 = m_meshTri.newNode({0, 1, 0});
		Node n_6 = m_meshTri.newNode({1, 1, 0});
		Node n_7 = m_meshTri.newNode({2, 1, 0});
		Node n_8 = m_meshTri.newNode({0.5, 1.5, 0});
		Node n_9 = m_meshTri.newNode({1.5, 1.5, 0});
		Node n_10 = m_meshTri.newNode({0, 2, 0});
		Node n_11 = m_meshTri.newNode({1, 2, 0});
		Node n_12 = m_meshTri.newNode({2, 2, 0});

		m_meshTri.newTriangle(n_0, n_1, n_3);
		m_meshTri.newTriangle(n_1, n_4, n_3);
		m_meshTri.newTriangle(n_1, n_2, n_4);
		m_meshTri.newTriangle(n_0, n_3, n_5);
		m_meshTri.newTriangle(n_3, n_4, n_6);
		m_meshTri.newTriangle(n_4, n_2, n_7);
		m_meshTri.newTriangle(n_5, n_3, n_8);
		m_meshTri.newTriangle(n_3, n_6, n_8);
		m_meshTri.newTriangle(n_6, n_9, n_8);
		m_meshTri.newTriangle(n_6, n_4, n_9);
		m_meshTri.newTriangle(n_4, n_7, n_9);
		m_meshTri.newTriangle(n_5, n_8, n_10);
		m_meshTri.newTriangle(n_10, n_8, n_11);
		m_meshTri.newTriangle(n_8, n_9, n_11);
		m_meshTri.newTriangle(n_11, n_9, n_12);
		m_meshTri.newTriangle(n_9, n_7, n_12);
	}

	Variable<math::Vector3d>* var_vectors = m_meshTri.newVariable<math::Vector3d, GMDS_NODE>("VectorField");
	for (auto const& n_id:m_meshTri.nodes())
	{
		var_vectors->set(n_id, {0.5,0.5,0});
	}

	gmds::MeshDoctor doc(&m_meshTri);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	// Second mesh
	Mesh m_meshQuad(MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));
	GridBuilder gb(&m_meshQuad,2);
	gb.execute(4,0.5, 4, 0.5);
	gmds::MeshDoctor doc_quad(&m_meshQuad);
	doc_quad.buildEdgesAndX2E();
	doc_quad.updateUpwardConnectivity();

	// Check the Mesh Alignment
	MeshAlignment_2D align(&m_meshTri, var_vectors, &m_meshQuad);
	MeshAlignment_2D::STATUS align_result = align.execute();

	ASSERT_EQ(align_result, MeshAlignment_2D::SUCCESS);

	Variable<double>* var_angle_deviation = align.getVariableDeviation();

	ASSERT_NEAR(var_angle_deviation->value(5), 45.0, pow(10,-6));
	ASSERT_NEAR(var_angle_deviation->value(6), 45.0, pow(10,-6));
	ASSERT_NEAR(var_angle_deviation->value(9), 45.0, pow(10,-6));
	ASSERT_NEAR(var_angle_deviation->value(10), 45.0, pow(10,-6));

}