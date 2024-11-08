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
#include "gmds/medialaxis/ConformalMeshBuilder.h"
#include "gmds/medialaxis/BlockStructureSimplifier.h"
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


	// First we quantize the non conformal mesh and build the corresponding conformal mesh
	QuantizationSolver qs(m);
	std::vector<std::vector<Node>> problematicCouples = qs.buildCompleteSolution();
	bool valid = problematicCouples.empty();
	if (valid)
	{
		// We can build the quantized mesh

		qs.setHalfEdgesLength();

		// Build the quantized mesh
		ConformalMeshBuilder cmb(m,qs.halfEdges(),qs.halfEdgesLengths());
		cmb.execute();
		cmb.writeQuantizedMesh("quantized_mesh.vtk");

		Mesh quantized_mesh = cmb.getQuantizedMesh();

	
		// Now we simplify the quantized mesh, ie we try to minimize the number of blocks
		BlockStructureSimplifier s(quantized_mesh);
		s.execute();
		s.writeSimplifiedMesh("simplified_mesh.vtk");
	}
	else
	{
		//We need to add a singularity dipole

		qs.buildFixedMesh(problematicCouples[0][0].id(),problematicCouples[0][1].id());
		qs.writeFixedMesh("fixed_mesh.vtk");
		qs.setFixedMeshConnectivity();
		Mesh fixed_mesh = qs.getFixedMesh();

		// Now we quantized the fixed mesh

		QuantizationSolver qs3(fixed_mesh);
		qs3.buildCompleteSolution();

		qs3.setHalfEdgesLength();

		// Build the quantized mesh
		ConformalMeshBuilder cmb(fixed_mesh,qs3.halfEdges(),qs3.halfEdgesLengths());
		cmb.execute();
		cmb.writeQuantizedMesh("quantized_mesh.vtk");
		
		Mesh quantized_mesh = cmb.getQuantizedMesh();

	
		// Now we simplify the quantized mesh, ie we try to minimize the number of blocks
		BlockStructureSimplifier s(quantized_mesh);
		s.execute();
		s.writeSimplifiedMesh("simplified_mesh.vtk");
	}

	

	// Write the output file
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(file_out);

	// // Create a test non-conformal mesh
	// Mesh m2(MeshModel(DIM3 | F | E | N |
	//                  F2N | F2E |
	//                  E2F | E2N | N2E | N2F));

	// IGMeshIOService ioService2(&m2);

	// m2.newNode(0.,0.);
	// m2.newNode(12.,0.);
	// m2.newNode(2.,1.);
	// m2.newNode(4.,2.);
	// m2.newNode(10.,2.);
	// m2.newNode(0.,5.);
	// m2.newNode(2.,5.);
	// m2.newNode(4.,5.);
	// m2.newNode(10.,7.);
	// m2.newNode(12.,7.);
	// m2.newNode(14.,7.);
	// m2.newNode(4.,10.);
	// m2.newNode(10.,10.);
	// m2.newNode(12.,11.);
	// m2.newNode(2.,12.);
	// m2.newNode(14.,12.);

	// m2.newFace({0,1,4,3,2});
	// m2.newFace({0,2,6,5});
	// m2.newFace({2,3,7,6});
	// m2.newFace({4,1,9,8});
	// m2.newFace({6,7,11,14});
	// m2.newFace({8,9,13,12});
	// m2.newFace({9,10,15,13});
	// m2.newFace({11,12,13,15,14});
	
	// VTKWriter vtkWriter2(&ioService2);
	// vtkWriter2.setCellOptions(N| E| F);
	// vtkWriter2.setDataOptions(N| E| F);
	// vtkWriter2.write("blocks_decomp_test.vtk");

}