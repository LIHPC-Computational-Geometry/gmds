//
// Created by chenyt on 26/09/24.
//
#include "gmds/ig/MeshDoctor.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKReader.h"
#include "gmds/io/VTKWriter.h"
#include "gmds/medialaxis/CrossField.h"
#include "gmds/medialaxis/MedialAxis2D.h"
#include "gmds/medialaxis/MedialAxis2DBuilder.h"
#include "gmds/medialaxis/MedialAxis3D.h"
#include "gmds/medialaxis/MedialAxis3DBuilder.h"
#include "gmds/medialaxis/MedialAxisMath.h"
#include "gmds/medialaxis/NonConformalHalfEdge.h"
#include "gmds/medialaxis/QuantizationGraph.h"
#include "gmds/medialaxis/QuantizationSolver.h"
#include <gmds/math/Point.h>
#include <gmds/math/Tetrahedron.h>
#include <iostream>
using namespace gmds;
using namespace math;


int main(int argc, char* argv[])
{
	std::cout << "============== Quantization ================" << std::endl;

	//==================================================================
	// PARAMETERS' PARSING
	//==================================================================
	std::string file_mesh, file_out;
	if (argc != 2) {
		std::cout << "Requires one parameters : \n";
		std::cout << "  - [IN ] a non conformal quad mesh (.vtk) \n"<<std::endl;
		throw gmds::GMDSException("Wrong number of parameters");
	}

	file_mesh = std::string (argv[1]);
	file_out = "block_decomp_out.vtk";
	std::cout << "Parameters " << std::endl;
	std::cout << "  - Input block decomposition file    : " << file_mesh << std::endl;
	std::cout << "  - Output block decomposition file  : " << file_out << std::endl;
	std::cout << "=======================================" << std::endl;


	//==================================================================
	// MESH READING
	//==================================================================

	std::cout<<"> Start mesh reading"<<std::endl;
	Mesh m(MeshModel(DIM3 | F | E | N | R |
	                 F2N | F2E |
	                 E2F | E2N | N2E | N2F));

	IGMeshIOService ioService(&m);
	VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(N| E| F);
	vtkReader.read(file_mesh);
	MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	// Test non conformal edges
//	Face f = m.get<Face>(6);
//	Edge e1 = m.get<Edge>(9);
//	Edge e2 = m.get<Edge>(10);
//	std::vector<Edge> edges;
//	edges.push_back(e1);
//	edges.push_back(e2);
//	NonConformalHalfEdge he(0,f,edges);
//	std::vector<NonConformalHalfEdge> half_edges;
//	half_edges.push_back(he);
//	std::vector<std::vector<Edge>> G = groupsOfAlignedEdges(f);
//	for (auto g:G)
//	{
//		std::cout<<"Test "<<std::endl;
//		for (auto e:g)
//			std::cout<<" "<<e.id()<<" ";
//		std::cout<<std::endl;
//	}

	// Build a quantization graph
	QuantizationSolver qs(m);
	qs.buildQuantizationGraph();
	QuantizationGraph* g = qs.getQuantizationGraph();
	g->updateConnectivity();
	g->buildQuantizationSolution();
	// g->display(); // Display the graph
	// g->displaySolution(); // Display the quantization solution
	// std::cout<<"Test "<<qgb.oppositeInQuad(6)<<std::endl;
	qs.setHalfEdgesLength();
	qs.buildQuantizedMeshNodesOnEdges();
	qs.buildQuantizedMeshInternalNodes();
	qs.buildQuantizedMeshFaces();
	qs.writeQuantizedMesh("quantized_mesh.vtk");

	// Write the output file
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(file_out);

//	// Create a test non-conformal mesh
//	Mesh m2(MeshModel(DIM3 | F | E | N |
//	                 F2N | F2E |
//	                 E2F | E2N | N2E | N2F));
//
//	IGMeshIOService ioService2(&m2);
//	m2.newNode(0.,0.);
//	m2.newNode(2.,0.);
//	m2.newNode(4.,0.);
//	m2.newNode(0.,1.);
//	m2.newNode(1.,1.);
//	m2.newNode(2.,1.);
//	m2.newNode(3.,1.);
//	m2.newNode(4.,1.);
//	m2.newNode(0.,2.);
//	m2.newNode(1.,2.);
//	m2.newNode(3.,2.);
//	m2.newNode(4.,2.);
//	m2.newNode(0.,3.);
//	m2.newNode(4.,3.);
//	m2.newFace({0,1,5,4,3});
//	m2.newFace({1,2,7,6,5});
//	m2.newFace({3,4,9,8});
//	m2.newFace({4,5,6,10,9});
//	m2.newFace({6,7,11,10});
//	m2.newFace({8,9,10,11,13,12});
//	VTKWriter vtkWriter2(&ioService2);
//	vtkWriter2.setCellOptions(N| E| F);
//	vtkWriter2.setDataOptions(N| E| F);
//	vtkWriter2.write("blocks_decomp_test.vtk");

}