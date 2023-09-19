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

	//==================================================================
	// PARAMETERS' PARSING
	//==================================================================
	std::string file_geom, file_mesh, file_out;
	if (argc != 4) {
		std::cout << "Require three paramaters : \n";
		std::cout << "  - [IN ] tetrahedral mesh (.vtk) that describes the geometry, \n";
		std::cout << "  - [IN ] mesh (.vtk) that describes the blocking to be classified, \n";
		std::cout << "  - [OUT] the name of the classified blocking (.vtk). \n" << std::endl;
		throw gmds::GMDSException("Wrong number of parameters");
	}

	file_geom = std::string(argv[1]);
	file_mesh = std::string(argv[2]);
	file_out = std::string(argv[3]);
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
	//the used model is specified according to the geom smoother requirements.
	Mesh m(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::R|gmds::R2N));

	IGMeshIOService ioService2(&m);
	VTKReader vtkReader2(&ioService2);
	vtkReader2.setCellOptions(N|R);
	vtkReader2.read(file_mesh);
	MeshDoctor doc2(&m);
	doc2.updateUpwardConnectivity();

	//==================================================================
	// CREATE THE BLOCKING FROM THE MESH
	//==================================================================
	gmds::blocking::CurvedBlocking bl(&geom_model, false);
	bl.init_from_mesh(m);

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

	gmds::IGMeshIOService ios(&blockingMesh);
	gmds::VTKWriter vtk_writer(&ios);
	vtk_writer.setCellOptions(gmds::N|gmds::R);
	vtk_writer.setDataOptions(gmds::N|gmds::R);
	vtk_writer.write(file_out);


}