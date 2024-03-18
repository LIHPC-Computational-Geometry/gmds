/*------------------------------------------------------------------------*/
#include "gmds/ig/MeshDoctor.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKReader.h"
#include "gmds/io/VTKWriter.h"
#include <gmds/morphMesh/MorphMesh.h>
#include <gmds/morphMesh/EllipticMorph.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/
int main(int argc, char* argv[]){

	std::string param_file(argv[1]);

	/*
	 * Mesh model needs to be changed to reduce mesh cells and memory usage
	 */
	gmds::Mesh m_mesh(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::F | gmds::E | gmds::N | gmds::N2E | gmds::N2F | gmds::N2R | gmds::E2N  | gmds::E2R | gmds::R2N | gmds::R2F | gmds::F2N | gmds::F2R));

	morphmesh::EllipticMorph emorph(param_file, &m_mesh);

	emorph.execute();
}