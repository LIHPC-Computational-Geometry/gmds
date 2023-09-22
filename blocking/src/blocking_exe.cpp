/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/blocking/Blocking.h>
#include <gmds/blocking/CurvedBlocking.h>
#include <gmds/blocking/CurvedBlockingClassifier.h>

#include "fstream"
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
void set_up_geom_model(gmds::cad::FACManager* AGeomModel, const std::string AFileName)
{
	gmds::Mesh vol_mesh(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::F | gmds::E | gmds::N | gmds::R2N | gmds::R2F | gmds::R2E | gmds::F2N | gmds::F2R | gmds::F2E
	                                    | gmds::E2F | gmds::E2N | gmds::N2E));
	//std::string dir(TEST_SAMPLES_DIR);
	//std::string vtk_file = dir +"/"+ AFileName;
	std::string vtk_file = AFileName;
	gmds::IGMeshIOService ioService(&vol_mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::R);
	vtkReader.read(vtk_file);
	gmds::MeshDoctor doc(&vol_mesh);
	doc.buildFacesAndR2F();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
	AGeomModel->initFrom3DMesh(&vol_mesh);

}

int main(int argc, char* argv[])
{
	std::cout << "============== CLASSIFICATION ================" << std::endl;

	/*
	 * The folder with the shapes to test needs to be defined like that :
	 * DataFolder/
	 * 	||
	 *   	\/
	 * -folderShape1 / folderShape1.vtk (the geometry file) folderShape1_blocking.vtk (the blocking file)
	 * -folderShape2 / folderShape2.vtk (the geometry file) folderShape2_blocking.vtk (the blocking file)
	 * ....
	 * -folderShapeX / folderShapeX.vtk (the geometry file) folderShapeX_blocking.vtk (the blocking file)
	 *
	 * The file on entry is just a list of the folder's name:
	 * folderShape1
	 * folderShape2
	 * ...
	 * folderShapeX
	*/
	if (argc != 2) {
		std::cout << "Require one paramater : \n";
		std::cout << "  - [IN ] the file with the name of the folder to go, \n";
		throw gmds::GMDSException("Wrong number of parameters");
	}
	auto file = std::string(argv[1]);
	//We browse the file until the end
	std::ifstream m_stream(file);
	std::string line;
	while(std::getline(m_stream,line)){

		//==================================================================
		// PARAMETERS PARSING
		//==================================================================
		std::string file_geom, file_mesh, file_out;

		//The path of the data folder
		std::string path_data_folder = "/home/bourmaudp/Documents/DATA/Easy_Shapes_Class/";
		//Path of the current folder with the geometry and the blocking
		std::string path_folder_shape = path_data_folder+line+"/";
		//Geometry file
		file_geom = path_folder_shape+line+".vtk";
		//Blocking file
		file_mesh = path_folder_shape+line+"_blocking.vtk";
		//The path to save the classification and the name of the file
		file_out = path_folder_shape+line+"_blocking_class_save.vtk";
		std::cout << "Parameters " << std::endl;
		std::cout << "  - Geometry file: " << file_geom << std::endl;
		std::cout << "  - Mesh file    : " << file_mesh << std::endl;
		std::cout << "  - Output file  : " << file_out << std::endl;
		std::cout << "=======================================" << std::endl;

		//==================================================================
		// GEOMETRY READING
		//==================================================================
		gmds::cad::FACManager geom_model;
		set_up_geom_model(&geom_model,file_geom);

		//==================================================================
		// MESH READING
		//==================================================================
		std::cout<<"> Start mesh reading"<<std::endl;
		//the used model is specified according to the class requirements.
		Mesh m(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::R|gmds::R2N));

		IGMeshIOService ioService2(&m);
		VTKReader vtkReader2(&ioService2);
		vtkReader2.setCellOptions(N|R);
		vtkReader2.read(file_mesh);
		MeshDoctor doc2(&m);
		doc2.updateUpwardConnectivity();

		std::cout<<"MESH Blocking INFO : N, "<<m.getNbNodes() <<" , R, "<<m.getNbRegions()<<std::endl;

		//==================================================================
		// CREATE THE BLOCKING FROM THE MESH
		//==================================================================
		gmds::blocking::CurvedBlocking bl(&geom_model, false);
		bl.init_from_mesh(m);
		std::cout<<"BL INFO : N, "<<bl.get_all_nodes().size() <<" , B, "<<bl.get_all_blocks().size()<<std::endl;


		//==================================================================
		// CLASSIFICATION BETWEEN THE BLOCKING AND THE GEOMETRY
		//==================================================================
		gmds::blocking::CurvedBlockingClassifier classifier(&bl);
		classifier.clear_classification();
		auto errors = classifier.classify();

		//==================================================================
		// SAVE CLASSIFICATION
		//==================================================================
		gmds::Mesh blockingMesh(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::E|gmds::F|gmds::R|gmds::E2N|gmds::F2N|gmds::R2N));
		bl.convert_to_mesh(blockingMesh);
		std::cout<<"BL CONVERT MESH INFO : N, "<<blockingMesh.getNbNodes() <<" , R, "<<blockingMesh.getNbRegions()<<std::endl;

		gmds::IGMeshIOService ios(&blockingMesh);
		gmds::VTKWriter vtk_writer(&ios);
		vtk_writer.setCellOptions(gmds::N|gmds::R);
		vtk_writer.setDataOptions(gmds::N|gmds::R);
		vtk_writer.write(file_out);
	}
}