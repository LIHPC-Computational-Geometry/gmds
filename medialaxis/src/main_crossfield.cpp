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
#include <gmds/math/Point.h>
#include <gmds/math/Tetrahedron.h>
#include <iostream>
using namespace gmds;
using namespace math;


int main(int argc, char* argv[])
{
	std::cout << "============== Build a cross field using the medial axis ================" << std::endl;

	//==================================================================
	// PARAMETERS' PARSING
	//==================================================================
	std::string file_mesh, file_out, file_ma_out;
	file_out = "out.vtk";
	file_ma_out = "medax.vtk";
	int nb_iterations=0;
	if (argc != 3) {
		std::cout << "Requires two parameters : \n";
		std::cout << "  - [IN ] minimal Delaunay mesh mesh (.vtk) (to build the medial axis) \n";
		std::cout << "  - [IN ] triangular mesh (on which the cross field is built) \n"<<std::endl;
		throw gmds::GMDSException("Wrong number of parameters");
	}

	file_mesh = std::string (argv[1]);
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

	std::cout << "NB nodes : " << m.getNbNodes() << std::endl;

	// Create a 2D medial axis
	medialaxis::MedialAxis2DBuilder mb(m);
	auto st = mb.execute(); // Warning : to build a cross field with the medial axis, activate the boundary connected components connexion in the execute function
	if (st == gmds::medialaxis::MedialAxis2DBuilder::SUCCESS)
	{
		// Get the medial axis
		medialaxis::MedialAxis2D* ma = mb.getMedialObject();
		// Write the medial axis
		ma->write(file_ma_out);
		// Get and write the smoothed medial axis
		medialaxis::MedialAxis2D* smoothed_ma = mb.getSmoothedMedialObject();
		smoothed_ma->write("smoothed_medax.vtk");
		// Get and display the singularities
		std::vector<Node> singularPoints = smoothed_ma->getSingularNodes();
		std::vector<double> singularityIndexes = smoothed_ma->getSingularityIndexes();
		for (int i = 0; i < singularPoints.size(); i++)
		{
			std::cout<<"Singularity at location "<<singularPoints[i].point()<<" of index "<<singularityIndexes[i]<<std::endl;
		}
		// Connect the boundary connected components using the medial axis
		std::vector<std::vector<math::Point>> connectedPoints = mb.m_connected_boundary_points;


		// Test topological representation
		//==================================================================
		// MESH READING
		//==================================================================
		std::string tri_mesh = std::string (argv[2]);
		Mesh m2(MeshModel(DIM3 | F | E | N | R |
		                  F2N | F2E |
		                  E2F | E2N | N2E | N2F));

		IGMeshIOService ioService2(&m2);
		VTKReader vtkReader2(&ioService2);
		vtkReader2.setCellOptions(N| E| F);
		vtkReader2.read(tri_mesh);
		MeshDoctor doc2(&m2);
		doc2.buildEdgesAndX2E();
		doc2.updateUpwardConnectivity();

		// Cross field
		medialaxis::CrossField cf_test(m2);
		cf_test.deleteIsolatedPoints();
		// Provide the triangular mesh with singularities located with the medial axis
		cf_test.setSingularitiesIndexes(singularPoints,singularityIndexes);
		// Connect the boundaries using the medial axis
		std::vector<std::vector<Node>> connexions = cf_test.connectedBoundaryNodes(connectedPoints);
		cf_test.setBoundariesConnexions(connexions);

		// Build the cross field
		cf_test.buildCrossField();

		// Write the cross field
		VTKWriter vtkWriter2(&ioService2);
		vtkWriter2.setCellOptions(N| E| F);
		vtkWriter2.setDataOptions(N| E| F);
		vtkWriter2.write("tri_out.vtk");

	}

	// Write the output file
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(file_out);
}