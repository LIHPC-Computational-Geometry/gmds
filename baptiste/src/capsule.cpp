#include <gmds/baptiste/capsule.h>
#include "gmds/io/IGMeshIOService.h"
#include <gmds/io/VTKReader.h>

using namespace gmds;

Capsule::Capsule() : m_mesh(MeshModel(DIM2|F|N|F2N|N2F)) {}

void Capsule::readMesh(std::string filename)
{
	gmds::IGMeshIOService ioService(&m_mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(filename);
}