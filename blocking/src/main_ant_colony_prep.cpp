/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/cad/GeomPoint.h>
#include <gmds/cad/GeomCurve.h>
#include <gmds/cad/GeomSurface.h>
#include <gmds/blocking/CurvedBlocking.h>
#include <gmds/blocking/CurvedBlockingClassifier.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <gmds/igalgo/BoundaryOperator.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
void setUpGeom(cad::FACManager &AGeomManager, const std::string& AFileName)
{
	Mesh m_vol(MeshModel(DIM3 | R | F | E | N |
	                                 R2N | R2F | R2E | F2N |
	                                 F2R | F2E
	                                 | E2F | E2N | N2E));
	std::string vtk_file = AFileName;
	IGMeshIOService ioService(&m_vol);
	VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(N | R);
	vtkReader.read(vtk_file);
	MeshDoctor doc(&m_vol);
	doc.buildFacesAndR2F();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
	
	AGeomManager.initFrom3DMesh(&m_vol);
}


/*----------------------------------------------------------------------------*/
void setUpBlock(Mesh& AMesh, const std::string& AFileName)
{
	std::string vtk_file = AFileName;
	IGMeshIOService ioService(&AMesh);
	VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(N | R);
	vtkReader.read(vtk_file);
}
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
	std::cout << "============== Ant Colony Prep ================" << std::endl;

	//==================================================================
	// PARAMETERS' PARSING
	//==================================================================
	std::string file_geom, file_block, file_block_out;
	if (argc != 4) {
		std::cout << "Require four paramaters : \n";
		std::cout << "  - [IN ] tetrahedral mesh (.vtk) that describes the geometry, \n";
		std::cout << "  - [IN ] mesh (.vtk) that describes the block structure\n";
		std::cout << "  - [OUT] mesh (.vtk) that provides the block structure with ant colony info"<< std::endl;
		throw GMDSException("Wrong number of parameters");
	}

	file_geom = std::string(argv[1]);
	file_block = std::string(argv[2]);
	file_block_out = std::string(argv[3]);
	std::cout << "Parameters " << std::endl;
	std::cout << "  - Geometry file: " << file_geom << std::endl;
	std::cout << "  - Mesh file    : " << file_block << std::endl;
	std::cout << "  - Output block file  : " << file_block_out << std::endl;
	std::cout << "=======================================" << std::endl;

	cad::FACManager geom_model;
	setUpGeom(geom_model,file_geom);
	
	Mesh m(MeshModel(DIM3 | N |  R | R2N));
	setUpBlock(m,file_block);
	blocking::CurvedBlocking bl(&geom_model, false);
	bl.init_from_mesh(m);
	blocking::CurvedBlockingClassifier cl(&bl);
	cl.classify();
	
	//=== create the blcok out structure =========
	Mesh m_out(MeshModel(DIM3 | N | E | F | R | E2N | F2N | R2N | R2E| F2E|R2F |E2R|F2R));
	bl.convert_to_mesh(m_out);
	MeshDoctor doc(&m_out);
	doc.buildFacesAndR2F();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
	std::cout<<"OUT EDGES: "<<m_out.getNbEdges()<<std::endl;
	std::cout<<"OUT FACES: "<<m_out.getNbFaces()<<std::endl;

	auto var_edge_block = m_out.newVariable<int,GMDS_EDGE>("hex_number_edge");
	auto var_face_block = m_out.newVariable<int,GMDS_FACE>("hex_number_face");
	auto var_face_solution = m_out.newVariable<int,GMDS_FACE>("face_in_solution");
	var_face_solution->set(19,1);
	var_face_solution->set(27,1);/*
	var_face_solution->set(281,1);
	var_face_solution->set(304,1);
	var_face_solution->set(277,1);
	var_face_solution->set(300,1);
	var_face_solution->set(323,1);
	var_face_solution->set(327,1);
	var_face_solution->set(381,1);
	var_face_solution->set(378,1);
	var_face_solution->set(389,1);
	var_face_solution->set(385,1);*/
	for(auto f_id:m_out.faces()){
		Face fi = m_out.get<Face>(f_id);
		var_face_block->set(f_id,fi.get<Region>().size());
	}
	for(auto e_id:m_out.edges()){
		Edge ei = m_out.get<Edge>(e_id);
		var_edge_block->set(e_id,ei.get<Region>().size());
	}
	IGMeshIOService ioService(&m_out);
	VTKWriter writer(&ioService);
	writer.setCellOptions(N |E|F|R);
	writer.setDataOptions(N |E|F|R);
	writer.write(file_block_out);

	std::cout << "======== Task done =========" << std::endl;
}