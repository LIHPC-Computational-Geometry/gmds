/*------------------------------------------------------------------------*/
// Created by Claire Roche on 21/10/2021.
/*------------------------------------------------------------------------*/
#ifdef USE_CGNS
   #include <gmds/blocking/CGNSWriter.h>
#endif
#include <gmds/claire/AbstractAeroPipeline.h>
#include <gmds/claire/AeroPipeline_2D.h>
#include <gmds/claire/AeroPipeline_3D.h>
#include <gmds/claire/Params.h>
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

	if (argc != 4 )
	{
		std::cout << "Merci de préciser <repertoire/de/travail> <fichier_param.ini> <fichier_sorti.cgns>" << std::endl;
		exit(0);
	}

	std::string dir(argv[1]);
	std::string param_file(argv[2]);
	std::string output_file(argv[3]);

	std::string input_file=dir+param_file;

	// Mesh Generation
	AeroPipeline_2D algo_aero2D(input_file, dir);
	AbstractAeroPipeline::STATUS aero2D_result = algo_aero2D.execute();

	if(aero2D_result == AbstractAeroPipeline::SUCCESS) {
#ifdef USE_CGNS
		blocking::CGNSWriter writer(algo_aero2D.getBlocking());
		writer.write(output_file, dir);
#else
		std::cout<<"CGNS export is desactivated"<<std::endl;
#endif
	}else{
		std::cout<<"Erreur dans le pipeline Aero"<<std::endl;
	}
}