//
// Created by Paul Bourmaud on 18/03/2022.
//

#ifndef GMDS_PAULTESTSUITE_H
#define GMDS_PAULTESTSUITE_H
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gtest/gtest.h>
#include <iostream>
#include <gmds/paul/Grid.h>
#include <gmds/paul/Tools.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(PaulTestClass, testGrid_1)
{
	// Test
	Mesh m(MeshModel(DIM3 | R | F | N | F2N ));
	Mesh m2(MeshModel(DIM3 | R | F | N | F2N));
	Node n4=m2.newNode(0,0,0);
	m2.newTriangle(n4,n4,n4);
	Node n0 = m.newNode(0,0,0);
	Node n1 = m.newNode(0,1,0);
	Node n2 = m.newNode(1,1,0);
	Node n3 = m.newNode(1,0,0);
	Face f = m.newQuad(n0,n1,n2,n3);

	GridBuilderAround ga(&m,&m2,2);
	Tools tk(&ga);
	ASSERT_TRUE(tk.isValid(f));
}
#endif     // GMDS_PAULTESTSUITE_H
