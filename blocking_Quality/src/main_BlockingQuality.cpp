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
#include <gmds/blocking_Quality/BlockingQuality.h>
#include <gmds/igalgo/GridBuilder.h>

#include <gmds/blockMesher/BlockMesher.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/cad/GeomPoint.h>
#include <gmds/cad/GeomCurve.h>
#include <gmds/cad/GeomSurface.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/smoothy/LaplacianSmoother.h>
#include <gmds/blocking_Quality/LinkerBlockingGeom.h>
#include <gmds/blocking_Quality/ValidBlocking.h>


#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gmds/ig/MeshDoctor.h>

#include <gmds/blocking_Quality/mainBis1.h>


//================================================================================
using namespace gmds;
int main(int argc, char* argv[])
{
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

	std::vector<cad::GeomVolume*> vols;
	geom_manager.getVolumes(vols);
	TCoord bb_min[3], bb_max[3];
	vols[0]->computeBoundingBox(bb_min,bb_max);
	math::Point pmin(bb_min[0], bb_min[1],bb_min[2]);
	math::Point pmax(bb_max[0], bb_max[1],bb_max[2]);
	std::cout<<"Min "<<pmin<<", max "<<pmax<<std::endl;




	lkbg.execute(&linker);

	linker.writeVTKDebugMesh("linker_debug_bridge.vtk");

/*
	for(auto n_id : blocking.nodes()){

		Node n = blocking.get<Node>(n_id);
		std::cout<<n<<std::endl;
		std::cout<<"La dim : "<<linker.getGeomDim(n)<<" "<<"l'id : "<<linker.getGeomId(n)<<std::endl;
		if(linker.getGeomDim(n) == 1 ){
			auto p_pointeur = geom_manager.getPoint(linker.getGeomId(n));
			math::Point p = p_pointeur->point();
			auto p_id = p_pointeur->id();

			std::cout<<"L'id point : "<<p_id<<" coords : "<<p<<std::endl;

		}

	}*/
	std::map<std::vector<TCellID>,int> elements_No_Classified;
	ValidBlocking vb(&blocking,&geom_manager,&linker,&elements_No_Classified);
	//vb.execute();
	std::cout<<"Return validity : "<<vb.execute()<<std::endl;
	std::cout<<"NB Points Geometry: "<<geometry.getNbNodes()<<std::endl;
	std::cout<<"NB Points GeomManager "<<geom_manager.getNbPoints()<<std::endl;


	return 0;
}
