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
	std::cout << "============== Medial Axis ================" << std::endl;

	//==================================================================
	// PARAMETERS' PARSING
	//==================================================================
	std::string file_mesh, file_out, file_ma_out, dim;
	int nb_iterations=0;
	if (argc != 6) {
		std::cout << "Requires three parameters : \n";
		std::cout << "  - [IN ] tetrahedral mesh (.vtk) that describes the geometry, \n";
		std::cout << "  - [IN ] the dimension : 2D or 3D, \n";
		std::cout << "  - [OUT] the annotated minimal Delaunay mesh (.vtk). \n";
		std::cout << "  - [OUT] the medial axis mesh (.vtk). \n"<< std::endl;
		throw gmds::GMDSException("Wrong number of parameters");
	}

	file_mesh = std::string (argv[1]);
	dim = std::string(argv[2]);
	file_out = std::string(argv[3]);
	file_ma_out = std::string(argv[4]);
	std::cout << "Parameters " << std::endl;
	std::cout << "  - Mesh file    : " << file_mesh << std::endl;
	std::cout << "  - Output minimal Delaunay file  : " << file_out << std::endl;
	std::cout << "  - Output medial axis file  : " << file_ma_out << std::endl;
	std::cout << "=======================================" << std::endl;

	if (dim == "2D")
	{
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
			// Get and display the singularities
			std::vector<Node> singularPoints = smoothed_ma->getSingularNodes();
			std::vector<double> singularityIndexes = smoothed_ma->getSingularityIndexes();
			for (int i = 0; i < singularPoints.size(); i++)
			{
				std::cout<<"Singularity at location "<<singularPoints[i].point()<<" of index "<<singularityIndexes[i]<<std::endl;
			}
			// Connect the boundary connected components using the medial axis
			std::vector<std::vector<math::Point>> connectedPoints = mb.m_connected_boundary_points;


			///////// Test CrossField
			std::string tri_mesh = std::string (argv[5]);
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
			medialaxis::CrossField cf_test(m2);
			cf_test.setSingularitiesIndexes(singularPoints,singularityIndexes);
			std::vector<std::vector<Node>> connexions = cf_test.connectedBoundaryNodes(connectedPoints);
			cf_test.setBoundariesConnexions(connexions);
			std::vector<std::vector<Node>> connectedNodes = cf_test.connectedBoundaryNodes(connectedPoints);

			cf_test.buildCrossField();

			VTKWriter vtkWriter2(&ioService2);
			vtkWriter2.setCellOptions(N| E| F);
			vtkWriter2.setDataOptions(N| E| F);
			vtkWriter2.write("tri_out.vtk");

			/////////

			///////// Test topological representation
//			smoothed_ma->buildTopoRepNodes();
//			smoothed_ma->setSectionID();
//			smoothed_ma->buildTopoRepEdges();
//			smoothed_ma->setTopoRepConnectivity();
//
//			std::cout<<"Test Matrix : "<<std::endl;
//			Eigen::MatrixXd A = smoothed_ma->constraintMatrix();
//			int I = A.rows();
//			int J = A.cols();
//			for (int i = 0; i < I; i++)
//			{
//				for (int j = 0; j < J; j++)
//					std::cout<<A.coeff(i,j)<<" ";
//				std::cout<<std::endl;
//			}
//			std::cout<<"Test Kernel : "<<std::endl;
//			Eigen::MatrixXd ker = A.fullPivLu().kernel();
//			I = ker.rows();
//			J = ker.cols();
//			for (int i = 0; i < I; i++)
//			{
//				for (int j = 0; j < J; j++)
//					std::cout<<int(ker.coeff(i,j))<<" ";
//				std::cout<<std::endl;
//			}
//
//			Eigen::MatrixXd Res = A*ker;
//			double res = Res.norm();
//			std::cout<<"||A*ker||="<<res<<std::endl;
//
//
//
//
//			smoothed_ma->writeTopoRep("topo_rep.vtk");

			////////
		}

		// Write the output file
		VTKWriter vtkWriter(&ioService);
		vtkWriter.setCellOptions(N| E| F);
		vtkWriter.setDataOptions(N| E| F);
		vtkWriter.write(file_out);
	}


	if (dim == "3D")
	{
		//==================================================================
		// MESH READING
		//==================================================================

		std::cout<<"> Start mesh reading"<<std::endl;
		Mesh m(MeshModel(DIM3 | F | E | N | R |
		                 F2N | F2E | R2F | F2R | R2E | E2R | R2N | N2R |
		                 E2F | E2N | N2E | N2F));

		IGMeshIOService ioService(&m);
		VTKReader vtkReader(&ioService);
		vtkReader.setCellOptions(N| E| F| R);
		vtkReader.read(file_mesh);
		MeshDoctor doc(&m);
		doc.buildEdgesAndX2E();
		doc.buildFacesAndR2F();
		doc.updateUpwardConnectivity();

		std::cout << "NB nodes : " << m.getNbNodes() << std::endl;
		std::cout << "NB tetras : " << m.getNbRegions() << std::endl;

		// Create a 3D medial axis
			medialaxis::MedialAxis3DBuilder mb(m);
			auto st = mb.execute();
			if (st == gmds::medialaxis::MedialAxis3DBuilder::SUCCESS)
			{
				// Get the medial axis
				medialaxis::MedialAxis3D* ma = mb.getMedialObject();
				// Write the medial axis
				ma->write(file_ma_out);
			}



		// Write the output file
		VTKWriter vtkWriter(&ioService);
		vtkWriter.setCellOptions(N| E| F| R);
		vtkWriter.setDataOptions(N| E| F| R);
		vtkWriter.write(file_out);
	}

}