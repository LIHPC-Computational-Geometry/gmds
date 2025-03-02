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
#include "gmds/medialaxis/TrianglesRemover.h"
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
	
	// Get the ids of the boundary edges
	std::vector<TCellID> boundary_edges_ids;
	for (auto e_id:m.edges())
		boundary_edges_ids.push_back(e_id);

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
	// m = cleaned_min_del;

	// Remove isolated points from m
	for (auto n_id:m.nodes())
	{
		Node n = m.get<Node>(n_id);
		if (n.get<Face>().empty())
			m.deleteNode(n);
	}


	// Create a 2D medial axis
	//medialaxis::MedialAxis2DBuilder mb(cleaned_min_del,boundary_edges_ids);
	medialaxis::MedialAxis2DBuilder mb(m,boundary_edges_ids);
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
		double mesh_size = 80.*smoothed_ma->meanMedEdgeLength();
		medialaxis::MedaxBasedTMeshBuilder tmb(medax,m);
		// Build a quad block decomposition
		tmb.setSectionID();
		tmb.computeSectionType();
		//tmb.transformGraySectionsIntoRed();
		bool topMaker = false;
		if (topMaker)
			tmb.setTopMakerColoring();
		tmb.avoidIPsBecomingEPs();
		tmb.markRefinableSections();
		tmb.markForbiddenSingularPointsAndModifyBoundaryPoints();
		tmb.refineByAddingSingularNodes(mesh_size);
		tmb.buildTopoRepNodes();
		tmb.buildTopoRepEdges();
		tmb.setTopoRepConnectivity();
		if (topMaker)
		{
			tmb.computeDistanceToEndPoints();
			tmb.ensureType1ForTopoRepEndSections();
		}
		tmb.optimizeMedaxColoring();
		tmb.markEdgesGeneratingTriangles();
		tmb.markCornersOnMinDel();
		tmb.buildTMeshNodesFromMinDelNodes();
		//tmb.buildBlockDecompMedialAndBoundaryNodes();
		tmb.buildSection2MedialAndBoundaryNodesAdjacency();
		tmb.writeTopoRep("topo_rep.vtk");
		tmb.buildMiddleNodes();
		tmb.buildBlocks();

		

		tmb.transformDegenerateQuadsIntoTriangles();

		tmb.setTopoRepEdgesColor();
		tmb.writeTopoRep("topo_rep.vtk");
		tmb.setBlockDecompConnectivity();
		tmb.writeBlockDecomp("block_decomp2.vtk");

		// // Only to have nice pictures
		//std::cout<<"VIENS EN PAUSE MEC"<<std::endl;
		// tmb.markBlocksSeparatingEdges();
		// tmb.addBigTJunctions();


		//return 0;

		// If we want to proceed only with quads

		tmb.transformTrianglesIntoQuads();
		tmb.markTJunctions();

		tmb.writeBlockDecomp("block_decomp1.vtk");

		tmb.buildFinalTMesh();
		tmb.setFinalTMeshConnectivity();
		tmb.markInternalConstraintsOnFinalTMesh();
		tmb.writeFinalTMesh("t_mesh.vtk");
		
		Mesh t_mesh = tmb.getFinalTMesh();

		// // Quantize the T-mesh
		QuantizationSolver qs(t_mesh, mesh_size);
		//////////////// test
		// qs.buildQuantizationGraph();
		// tmb.writeFinalTMesh("t_mesh.vtk");
		////////////////
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
		if (topMaker)
		{
			s.traceTopMakerBlocksOutlines();
			s.setTopMakerBlocksIDs();
		}
			

		conf.projectOnBoundary(m);

		conf.writeConformalMesh("conformal_mesh_without_smoothing.vtk");


		bool remove_triangles = true;
		if (remove_triangles)
		{
			// for (int i = 0; i < 10; i++)
			// 	conf.smooth(m);
			// POST-PROCESSING TRIANGLES
			std::cout<<"========== Post-processing triangles =========="<<std::endl;
			// Building a T-mesh by placing quads templates in triangles fans
			TrianglesRemover triRem(conformal_mesh);
			triRem.buildDegeneratedNodesGroups();

			// std::vector<Face> fan0 = triRem.fan(4);
			// for (auto f:fan0)
			// 	std::cout<<"TEST "<<f.id()<<std::endl;
			// std::vector<Node> outline = triRem.outline(fan0);

			triRem.buildTMesh();
			triRem.writeTMesh("post_process_t_mesh.vtk");
			triRem.setTMeshConnectivity();
			triRem.markInternalConstraintsOnEdges();
			triRem.transformUnwantedTrianglesIntoQuads();
			triRem.writeTMesh("post_process_t_mesh.vtk");
			triRem.buildFinalTMesh();
			triRem.setFinalTMeshConnectivity();
			triRem.markInternalConstraintsOnFinalEdges();
			triRem.writeFinalTMesh("post_process_final_t_mesh.vtk");
			Mesh t_mesh2 = triRem.getFinalTMesh();
			// Quantize the T-mesh
			QuantizationSolver qs2(t_mesh2, mesh_size);
			qs2.buildCompleteSolution();
			// We can build the quantized mesh

			qs2.setHalfEdgesLength();
			qs2.setEdgesLength();

			triRem.writeTMesh("post_process_t_mesh.vtk");

			// Build the conformal mesh
			Conformalizer conf2(t_mesh2,qs2.halfEdges(),qs2.halfEdgesLengths());
			conf2.execute();
			Mesh conformal_mesh2 = conf2.getConformalMesh();
			
			// Now we extract the block structure
			BlockStructureSimplifier s2(conformal_mesh2);
			s2.markSeparatrices();
			s2.setBlocksIDs();

			conf2.projectOnBoundary(m);
			conf2.writeConformalMesh("post_processed_quad_mesh_without_smoothing.vtk");
			for (int i = 0; i < 10; i++)
				conf2.smooth(m);

			// Mesh quality
			auto scaledJacobian_with_post_process = conformal_mesh2.newVariable<double,GMDS_FACE>("scaled_jacobian");
			for (auto f_id:conformal_mesh2.faces())
			{
				Face f = conformal_mesh2.get<Face>(f_id);
				double sj = f.computeScaledJacobian2D();
				scaledJacobian_with_post_process->set(f_id,sj);
			}

			conf2.writeConformalMesh("post_processed_quad_mesh.vtk");
		}



		conf.writeConformalMesh("conformal_mesh_without_smoothing.vtk");
		conf.writeConformalMesh("contraintes_nous_not_smoothed.vtk");

		clock_t t = clock();
		for (int i = 0; i < 10; i++)
			conf.smooth(m);
		t = clock()-t;
		std::cout<<"Smoothing time (s) : "<<double(t)/CLOCKS_PER_SEC<<std::endl;

		// Mesh quality
		auto scaledJacobian_without_post_process = conformal_mesh.newVariable<double,GMDS_FACE>("scaled_jacobian");
		for (auto f_id:conformal_mesh.faces())
		{
			Face f = conformal_mesh.get<Face>(f_id);
			double sj = f.computeScaledJacobian2D();
			scaledJacobian_without_post_process->set(f_id,sj);
		}

		tmb.writeFinalTMesh("t_mesh.vtk");
		conf.writeConformalMesh("conformal_mesh.vtk");
		conf.writeConformalMesh("contraintes_nous_smoothed.vtk");

		// //Test if it is 2-manifold
		// std::cout<<"TEST "<<std::endl;
		// for (auto e_id:conformal_mesh.edges())
		// {
		// 	Edge e = conformal_mesh.get<Edge>(e_id);
		// 	if (e.get<Face>().size() != 2 && e.get<Face>().size() != 1)
		// 		std::cout<<"PROBLEM "<<e_id<<std::endl;
		// }
		
	}


	// Write the output file
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(file_out);


	//cleaner.writeCleanedMesh("cleaned_min_del.vtk");
}

