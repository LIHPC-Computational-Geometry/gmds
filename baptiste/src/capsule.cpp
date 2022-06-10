#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKWriter.h"
#include <gmds/baptiste/capsule.h>
#include <gmds/io/VTKReader.h>

using namespace gmds;

Capsule::Capsule() : m_mesh(MeshModel(DIM2|F|N|F2N|N2F)) {}

Capsule::~Capsule() {}

void Capsule::readMesh(std::string filename)
{
	gmds::IGMeshIOService ioService(&m_mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(filename);
}

void Capsule::saveMesh(std::string filename)
{
	IGMeshIOService ioService(&m_mesh);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write(filename);
}