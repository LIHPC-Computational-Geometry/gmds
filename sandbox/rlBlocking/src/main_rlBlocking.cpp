//
// Created by bourmaudp on 02/12/22.
//
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <string>


#include <gmds/ig/Mesh.h>
#include <gmds/ig/Node.h>
#include <gmds/quality/QuadQuality.h>
#include <gmds/quality/HexQuality.h>
#include <gmds/rlBlocking/BlockingQuality.h>
#include <gmds/igalgo/GridBuilder.h>

#include <gmds/blockMesher/BlockMesher.h>
#include <gmds/cad/GeomCurve.h>
#include <gmds/cad/GeomPoint.h>
#include <gmds/cad/GeomSurface.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/rlBlocking/LinkerBlockingGeom.h>
#include <gmds/rlBlocking/ValidBlocking.h>
#include <gmds/smoothy/LaplacianSmoother3C.h>

#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gmds/ig/MeshDoctor.h>

// Quality of classification include
#include <gmds/blocking/CurvedBlockingClassifier.h>


//================================================================================
using namespace gmds;
int main(int argc, char* argv[])
{
	std::cout << "============== TEST CurvedBlockingClassifier ================" << std::endl;

	cad::FACManager geom_model;
	gmds::Mesh vol_mesh(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::F | gmds::E | gmds::N | gmds::R2N | gmds::R2F | gmds::R2E | gmds::F2N | gmds::F2R | gmds::F2E
	                                    | gmds::E2F | gmds::E2N | gmds::N2E));

	std::string vtk_file = "/home/bourmaudp/Documents/mambo-master/Basic/vtk/cb3.vtk";
	gmds::IGMeshIOService ioServiceA(&vol_mesh);
	gmds::VTKReader vtkReaderA(&ioServiceA);
	vtkReaderA.setCellOptions(gmds::N | gmds::R);
	vtkReaderA.read(vtk_file);
	gmds::MeshDoctor docA(&vol_mesh);
	docA.buildFacesAndR2F();
	docA.buildEdgesAndX2E();
	docA.updateUpwardConnectivity();
	geom_model.initFrom3DMesh(&vol_mesh);

	gmds::blocking::CurvedBlocking bl(&geom_model,true);

	gmds::blocking::CurvedBlockingClassifier classifier(&bl);
	classifier.clear_classification();
	auto errors = classifier.classify();

	std::cout<<"Nb points Geom : "<<geom_model.getNbPoints()<<std::endl;
	std::cout<<"Nb points Blocking : "<<bl.get_all_nodes().size()<<std::endl;

	std::cout<<"liste points non capturés : "<<std::endl;
	for(auto pointsNoClass : errors.non_captured_points){
		std::cout<<"Point :"<<pointsNoClass<<std::endl;
	}
	std::cout<<"liste courbes non capturées : "<<std::endl;
	for(auto curvesNoClass : errors.non_captured_curves){
		std::cout<<curvesNoClass<<std::endl;
	}
	std::cout<<"liste surfaces non capturées : "<<std::endl;
	for(auto surfacesNoClass :errors.non_captured_surfaces){
		std::cout<<surfacesNoClass<<std::endl;
	}
	std::cout<<"Class do"<<std::endl;

	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::E|gmds::F|gmds::R|gmds::E2N|gmds::F2N|gmds::R2N));
	bl.convert_to_mesh(m);


	gmds::IGMeshIOService ios(&m);
	gmds::VTKWriter vtk_writer(&ios);
	vtk_writer.setCellOptions(gmds::N|gmds::R);
	vtk_writer.setDataOptions(gmds::N|gmds::R);
	vtk_writer.write("debug_blocking.vtk");
	gmds::VTKWriter vtk_writer_edges(&ios);
	vtk_writer_edges.setCellOptions(gmds::N|gmds::E);
	vtk_writer_edges.setDataOptions(gmds::N|gmds::E);
	vtk_writer_edges.write("debug_blocking_edges.vtk");





	return 0;

	std::cout << "============== Valid Blocks ================" << std::endl;

	//==================================================================
	// PARAMETERS' PARSING
	//==================================================================
	std::string file_geom, file_block;
	int nb_iterations=0;
	if (argc != 3) {
		std::cout << "Require two parameters : \n";
		std::cout << "  - [IN ] tetrahedral mesh (.vtk) that describes the geometry, \n";
		std::cout << "  - [IN ] mesh (.vtk) that describes the block structure \n";
		return 0;
	}

	file_geom = std::string(argv[1]);
	file_block = std::string(argv[2]);

	std::cout << "Parameters " << std::endl;
	std::cout << "  - Geometry file: " << file_geom << std::endl;
	std::cout << "  - Block file   : " << file_block << std::endl;
	std::cout << "=======================================" << std::endl;

	//==================================================================
	// GEOMETRY READING
	//==================================================================
	std::cout<<"> Start geometry reading"<<std::endl;
	Mesh geometry(MeshModel(DIM3 | R | F | E | N |
	                        R2N | R2F | R2E |
	                        F2N | F2R | F2E |
	                        E2F | E2N | N2E | N2R | N2F));

	IGMeshIOService ioService(&geometry);
	VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(N| R);
	vtkReader.read(file_geom);
	MeshDoctor doc(&geometry);
	doc.buildFacesAndR2F();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	cad::FACManager geom_manager;
	geom_manager.initFrom3DMesh(&geometry);

	std::cout<<"Geom ("<<geom_manager.getNbVolumes()<<", "
	          <<geom_manager.getNbSurfaces()<<", "
	          <<geom_manager.getNbCurves()<<", "
	          <<geom_manager.getNbPoints()<<")"<<std::endl;


	//==================================================================
	// MESH READING
	//==================================================================
	std::cout<<"> Start block reading"<<std::endl;
	//the used model is specified according to the geom smoother requirements.
	Mesh blocking(MeshModel(DIM3 | R | F | E | N |
	                        R2N | R2F | R2E | F2N | F2R | F2E | E2F | E2N | N2E | N2R | N2F));

	IGMeshIOService ioService2(&blocking);
	VTKReader vtkReader2(&ioService2);
	vtkReader2.setCellOptions(N|R);
	vtkReader2.read(file_block);
	MeshDoctor doc2(&blocking);
	doc2.buildFacesAndR2F();
	doc2.buildEdgesAndX2E();
	doc2.updateUpwardConnectivity();

	//==================================================================
	// TRY TO CLASSIFY MESH CELLS
	//==================================================================

	std::cout<<"> Start mesh->geometry classification"<<std::endl;
	LinkerBlockingGeom lkbg(&blocking,&geom_manager);
	cad::GeomMeshLinker linker(&blocking,&geom_manager);


	lkbg.execute(&linker);

	linker.writeVTKDebugMesh("linker_debug.vtk");

	ValidBlocking vb(&blocking,&geom_manager,&linker);
	vb.execute();
	std::cout<<"Return validity : "<<vb.execute()<<std::endl;




}
