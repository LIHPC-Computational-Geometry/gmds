#include "gmds/ig/MeshDoctor.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKReader.h"
#include "gmds/io/VTKWriter.h"
#include "gmds/io/MeshBWriter.h"
#include "gmds/medialaxis/MedialAxis2D.h"
#include "gmds/medialaxis/MedialAxis2DBuilder.h"
#include "gmds/medialaxis/MedialAxis3D.h"
#include "gmds/medialaxis/MedialAxis3DBuilder.h"
#include "gmds/medialaxis/CrossField.h"
#include "gmds/medialaxis/MedialAxisMath.h"
#include "gmds/medialaxis/MinDelaunayCleaner.h"
#include "gmds/medialaxis/MedaxBasedTMeshBuilder.h"
#include "gmds/medialaxis/QuantizationSolver.h"
#include "gmds/medialaxis/ConformalMeshBuilder.h"
#include "gmds/medialaxis/BlockStructureSimplifier.h"
#include "gmds/medialaxis/Conformalizer.h"
#include <gmds/math/Point.h>
#include <gmds/math/Tetrahedron.h>
#include <time.h> 
#include <iostream>
using namespace gmds;
using namespace math;


int main(int argc, char* argv[])
{
	std::cout << "============== Block decomposition using the medial axis ================" << std::endl;

	//==================================================================
	// PARAMETERS' PARSING
	//==================================================================
	std::string file_mesh, file_out, file_ma_out;
	int nb_iterations=0;
	if (argc != 2) {
		std::cout << "Requires one parameters : \n";
		std::cout << "  - [IN ] minimal Delaunay mesh (.vtk) (to build the medial axis) \n"<<std::endl;
		throw gmds::GMDSException("Wrong number of parameters");
	}

	file_mesh = std::string (argv[1]);
	file_out = "out.vtk";
	file_ma_out = "medax.vtk";
	std::cout << "Parameters " << std::endl;
	std::cout << "  - Mesh file    : " << file_mesh << std::endl;
	std::cout << "  - Output minimal Delaunay file  : " << file_out << std::endl;
	std::cout << "  - Output medial axis file  : " << file_ma_out << std::endl;
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

	// MinDelaunayCleaner cleaner(m);
	// cleaner.setFacesTypes();
	// cleaner.markSmallEdges();
	// cleaner.markFacesToDelete();
	// cleaner.buildCleanedMesh();
	// cleaner.setCleanedMeshConnectivity();

	// Mesh cleaned_min_del = cleaner.getCleanedMesh();



	// Create a 2D medial axis
	//medialaxis::MedialAxis2DBuilder mb(cleaned_min_del);
	medialaxis::MedialAxis2DBuilder mb(m);
	auto st = mb.execute();
	if (st == gmds::medialaxis::MedialAxis2DBuilder::SUCCESS)
	{
		// Get the medial axis
		medialaxis::MedialAxis2D* ma = mb.getMedialObject();
		// Write the medial axis
		ma->write(file_ma_out);
		// Get and write the smoothed medial axis
		medialaxis::MedialAxis2D* smoothed_ma = mb.getSmoothedMedialObject();
		smoothed_ma->write("smoothed_medax.vtk");

		Mesh medax = smoothed_ma->getMeshRepresentation();
		double mesh_size = 50.*smoothed_ma->meanMedEdgeLength();
		medialaxis::MedaxBasedTMeshBuilder tmb(medax,m);
		// Build a quad block decomposition
		tmb.setSectionID();
		tmb.computeSectionType();
		//tmb.transformGraySectionsIntoRed();
		//tmb.setTopMakerColoring();
		tmb.markRefinableSections();
		tmb.markForbiddenSingularPointsAndModifyBoundaryPoints();
		tmb.refineByAddingSingularNodes(mesh_size);
		tmb.buildTopoRepNodes();
		tmb.buildTopoRepEdges();
		tmb.setTopoRepConnectivity();
		tmb.optimizeMedaxColoring();
		tmb.markEdgesGeneratingTriangles();
		tmb.markCornersOnMinDel();
		tmb.buildTMeshNodesFromMinDelNodes();
		//tmb.buildBlockDecompMedialAndBoundaryNodes();
		tmb.buildSection2MedialAndBoundaryNodesAdjacency();
		tmb.buildMiddleNodes();
		tmb.buildBlocks();

		// tmb.writeTopoRep("topo_rep_test.vtk");
		// tmb.writeBlockDecomp("t_mesh_test.vtk");
		// return 0;

		tmb.transformDegenerateQuadsIntoTriangles();
		tmb.setTopoRepEdgesColor();
		tmb.writeTopoRep("topo_rep.vtk");
		// // Only to have nice pictures
		tmb.setBlockDecompConnectivity();
		// tmb.markBlocksSeparatingEdges();
		// tmb.addBigTJunctions();
		tmb.writeBlockDecomp("block_decomp.vtk");
		tmb.transformTrianglesIntoQuads();
		tmb.markTJunctions();

		tmb.writeBlockDecomp("block_decomp1.vtk");

		tmb.buildFinalTMesh();
		tmb.setFinalTMeshConnectivity();
		tmb.markInternalConstraintsOnFinalTMesh();
		tmb.writeFinalTMesh("t_mesh.vtk");
		
		Mesh t_mesh = tmb.getFinalTMesh();

		// Quantize the T-mesh
		QuantizationSolver qs(t_mesh, mesh_size);
		qs.buildCompleteSolution();
		// We can build the quantized mesh

		qs.setHalfEdgesLength();
		qs.setEdgesLength();

		// Build the conformal mesh
		Conformalizer conf(t_mesh,qs.halfEdges(),qs.halfEdgesLengths());
		conf.execute();
		Mesh conformal_mesh = conf.getConformalMesh();

		// Now we simplify the quantized mesh, ie we try to minimize the number of blocks
		BlockStructureSimplifier s(conformal_mesh);
		// s.execute();
		// s.writeSimplifiedMesh("simplified_mesh.vtk");
		s.markSeparatrices();
		s.setBlocksIDs();

		conf.projectOnBoundary(m);

		conf.writeConformalMesh("conformal_mesh_without_smoothing.vtk");

		clock_t t = clock();
		for (int i = 0; i < 100; i++)
			conf.smooth(m);
		t = clock()-t;
		std::cout<<"Smoothing time (s) : "<<double(t)/CLOCKS_PER_SEC<<std::endl;

		tmb.writeFinalTMesh("t_mesh.vtk");
		conf.writeConformalMesh("conformal_mesh.vtk");
	}

	// Write the output file
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(file_out);

	// cleaner.writeCleanedMesh("cleaned_min_del.vtk");
}

