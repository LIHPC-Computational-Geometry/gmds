#include "gmds/ig/MeshDoctor.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKReader.h"
#include "gmds/io/VTKWriter.h"
#include "gmds/medialaxis/MedialAxis2D.h"
#include "gmds/medialaxis/MedialAxis2DBuilder.h"
#include "gmds/medialaxis/MedialAxis3D.h"
#include "gmds/medialaxis/MedialAxis3DBuilder.h"
#include "gmds/medialaxis/CrossField.h"
#include "gmds/medialaxis/MedialAxisMath.h"
#include "gmds/medialaxis/MinDelaunayCleaner.h"
#include <gmds/math/Point.h>
#include <gmds/math/Tetrahedron.h>
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

	MinDelaunayCleaner cleaner(m);
	cleaner.setFacesTypes();
	cleaner.markSmallEdges();
	cleaner.markFacesToDelete();
	cleaner.buildCleanedMesh();
	cleaner.setCleanedMeshConnectivity();

	Mesh cleaned_min_del = cleaner.getCleanedMesh();



	// Create a 2D medial axis
	medialaxis::MedialAxis2DBuilder mb(cleaned_min_del);
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

		// Build a quad block decomposition
		smoothed_ma->setSectionID();
		smoothed_ma->computeSectionType();
		smoothed_ma->refineByAddingSingularNodes();
		smoothed_ma->buildTopoRepNodes();
		smoothed_ma->buildTopoRepEdges();
		smoothed_ma->setTopoRepConnectivity();
		smoothed_ma->buildBlockDecompMedialAndBoundaryNodes();
		smoothed_ma->buildSection2MedialAndBoundaryNodesAdjacency();
		smoothed_ma->buildMiddleNodes();
		smoothed_ma->buildBlocks();
		smoothed_ma->writeTopoRep("topo_rep.vtk");
		smoothed_ma->writeBlockDecomp("block_decomp.vtk");

	}

	// // Write the output file
	// VTKWriter vtkWriter(&ioService);
	// vtkWriter.setCellOptions(N| E| F);
	// vtkWriter.setDataOptions(N| E| F);
	// vtkWriter.write(file_out);
	cleaner.writeCleanedMesh("cleaned_min_del.vtk");
}

