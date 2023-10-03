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
	std::cout << "============== S3 ================" << std::endl;

	//==================================================================
	// PARAMETERS' PARSING
	//==================================================================
	std::string file_geom, file_block, file_block_out;
	if (argc != 2) {
		std::cout << "Require four paramaters : \n";
		std::cout << "  - [OUT] mesh (.vtk) that provides the output block structure"<< std::endl;
		throw GMDSException("Wrong number of parameters");
	}

	file_geom = std::string("S3_geom.vtk");
	file_block = std::string("S3_block.vtk");
	file_block_out = std::string(argv[1]);
	std::cout << "Parameters " << std::endl;
	std::cout << "  - Geometry file: " << file_geom << std::endl;
	std::cout << "  - Mesh file    : " << file_block << std::endl;
	std::cout << "  - Output block file  : " << file_block_out << std::endl;
	std::cout << "=======================================" << std::endl;

	//READ AND PREPARE GEOMETRY
	cad::FACManager geom_model;
	setUpGeom(geom_model,file_geom);

	//READ AND PREPARE BLOCK
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

	IGMeshIOService ioService(&m_out);
	VTKWriter writer(&ioService);
	writer.setCellOptions(N |E|F|R);
	writer.setDataOptions(N |E|F|R);
	writer.write(file_block_out);

	std::cout << "======== Task done =========" << std::endl;
}