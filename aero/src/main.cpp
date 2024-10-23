/*------------------------------------------------------------------------*/
// Created by Claire Roche on 21/10/2021.
/*------------------------------------------------------------------------*/
#ifdef USE_CGNS
   #include <gmds/blocking/CGNSWriter.h>
	#include <gmds/aero/CGNSWriter3D.h>
#endif
#include <gmds/aero/AbstractAeroPipeline.h>
#include <gmds/aero/AeroPipeline_2D.h>
#include <gmds/aero/AeroPipeline_3D.h>
#include <gmds/aero/Params.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace gmds;

/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
	std::cout << "=== AERO ALGO ====" << std::endl;

	if (argc != 5) {
		std::cout << "Merci de préciser <repertoire/de/travail> <fichier_param.ini> <fichier_sorti.cgns> <2D/3D>" << std::endl;
		exit(0);
	}

	std::string dir(argv[1]);
	std::string param_file(argv[2]);
	std::string output_file(argv[3]);
	std::string dimension(argv[4]);

	std::string input_file = dir +"/"+ param_file;

	if (dimension == "2D") {
		// Mesh Generation
		// Mesh Generation
		AeroPipeline_2D algo_aero2D(input_file, dir);
		AbstractAeroPipeline::STATUS aero2D_result = algo_aero2D.execute();

		if (aero2D_result == AbstractAeroPipeline::SUCCESS) {
#ifdef USE_CGNS
			blocking::CGNSWriter writer(algo_aero2D.getBlocking());
			writer.write(output_file, dir);
#else
			std::cout << "CGNS export is desactivated" << std::endl;
#endif
		}
	}
	else if (dimension == "3D") {
		// Mesh Generation
		dir = dir + "/";
		AeroPipeline_3D algo_aero3D(input_file, dir);
		AbstractAeroPipeline::STATUS aero_result = algo_aero3D.execute();

		if (aero_result == AbstractAeroPipeline::SUCCESS) {
#ifdef USE_CGNS
			gmds::aero::CGNSWriter3D writer(algo_aero3D.getBlocking());
			writer.write("", output_file, dir);
#else
			std::cout << "CGNS export is desactivated" << std::endl;
#endif
		}
	}
}
