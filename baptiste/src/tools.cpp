#include <gmds/baptiste/tools.h>
#include "gmds/io/IGMeshIOService.h"
#include <gmds/io/VTKReader.h>
#include "gmds/igalgo/VolFracComputation.h"

using namespace gmds;

std::vector<double> gmds::LinearSpacedArray(double a, double b, std::size_t N)
{
	double h = (b - a) / static_cast<double>(N);
	std::vector<double> v(N+1);
	std::vector<double>::iterator x;
	double val;
	for (x = v.begin(), val = a; x != v.end(); ++x, val += h)
	{
		*x = val;
	}
	return v;
}

/*
Mesh gmds::readMesh(std::string filename)
{
	Mesh mesh = Mesh(MeshModel(DIM2|F|N|F2N));
	gmds::IGMeshIOService ioService(&mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(filename);
	return mesh;
}
 */

/*
void gmds::computeVolFrac(Mesh AMesh, Mesh targetMesh)
{
	AMesh.newVariable<double, GMDS_FACE>("volFrac");
	volfraccomputation_2d(&AMesh, &targetMesh, AMesh.getVariable<double, GMDS_FACE>("volFrac"));
}
 */

Variable<double>* gmds::getVolFrac(Mesh mesh)
{
	Variable<double>* res = mesh.getVariable<double, GMDS_FACE>("volFrac");
	return res;
}

Tools::Tools() : m_mesh(MeshModel(DIM2|F|N|F2N)) {}

void Tools::readMesh(std::string filename)
{
	gmds::IGMeshIOService ioService(&m_mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(filename);
}

void Tools::computeVolFrac(Mesh AMesh, Mesh targetMesh)
{
	AMesh.newVariable<double, GMDS_FACE>("volFrac");
	volfraccomputation_2d(&AMesh, &targetMesh, AMesh.getVariable<double, GMDS_FACE>("volFrac"));
}