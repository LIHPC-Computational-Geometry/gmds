/*------------------------------------------------------------------------*/
// Created by Claire Roche on 21/10/2021.
/*------------------------------------------------------------------------*/
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
	std::cout << "=== CLAIRE ALGO ====" << std::endl;

	ParamsAero params;
	params.dim=ParamsAero::DIM_2D;
	params.input_file="...";

	AbstractAeroPipeline* algo = NULL;
	if(params.dim==ParamsAero::DIM_2D){
		algo = new AeroPipeline_2D(params);
	}
	else if(params.dim==ParamsAero::DIM_3D){
		algo = new AeroPipeline_3D(params);
	}
	else{
		std::cout<<" Wrong dimension "<<std::endl;
		exit(0);
	}
	algo->execute();

	delete algo;

}