#include <iostream>
#include "gmds/baptiste/RLBlockSet.h"
#include "gmds/ig/Mesh.h"
#include "gmds/io/IGMeshIOService.h"


//using namespace std;
using namespace gmds;

int main()
{
	//cout<<"Hello World";


	Mesh m(MeshModel(DIM3|F|N|F2N|N2F));
	//Node n0 = m.newNode(math::Point(0,0,0));
	//Node n1 = m.newNode(math::Point(0,3,0));
	//Node n2 = m.newNode(math::Point(3,3,0));


	RLBlockSet blockSet(&m);
	// RLBlockSet blockSet(m);

	/*
	double xMin = 0.2;
	double yMin = 0.5;
	double xMax = 5.3;
	double yMax = 4.7;

	blockSet.setFrame(xMin, yMin, xMax, yMax);


	int nbFaces = blockSet.m_mesh.getNbFaces();

	cout<<nbFaces;
	 */
	return 0;
}