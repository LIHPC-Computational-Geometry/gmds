#include <gmds/baptiste/tools.h>
#include "gmds/io/IGMeshIOService.h"
#include <gmds/io/VTKReader.h>
#include "gmds/igalgo/VolFracComputation.h"
#include "gmds/baptiste/RLBlockSet.h"

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

void gmds::cloneMesh(const Mesh &originalMesh, Mesh &newMesh)
{
	if (originalMesh.getDim() == 3 or newMesh.getDim() == 3)
	{
		throw GMDSException("Dimension must be 2");
	}

	for (int nodeID:originalMesh.nodes())
	{
		newMesh.newNode(originalMesh.get<Node>(nodeID).point());
	}
	std::set<TCellID> invalidIndexes;
	for (int faceID = 0; faceID <= originalMesh.getMaxLocalID(2); faceID++)
	{
		if (originalMesh.has<Face>(faceID))
		{
			newMesh.newFace(originalMesh.get<Face>(faceID).get<Node>());
		}
		else
		{
			newMesh.newFace({0, 0, 0, 0});
			invalidIndexes.insert(faceID);
		}
	}
	for (int faceID:invalidIndexes)
	{
		newMesh.deleteFace(faceID);
	}
}

void gmds::cloneBlockSet(const RLBlockSet &originalBlockSet, RLBlockSet &newBlockSet)
{
	cloneMesh(originalBlockSet.m_mesh, newBlockSet.m_mesh);
	newBlockSet.xSize = originalBlockSet.xSize;
	newBlockSet.ySize = originalBlockSet.ySize;
}


Mesh gmds::readMesh(std::string filename)
{
	Mesh mesh = Mesh(MeshModel(DIM2|F|N|F2N));
	gmds::IGMeshIOService ioService(&mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(filename);
	return mesh;
}
std::vector<Action *> gmds::getActionsVector()
{
	std::vector<Action *> actions;
	ActionDelete d = ActionDelete();
	actions.push_back(&d);
	for (int iV = 0; iV <= 1; iV++)
	{
		for (int iAxis = 0; iAxis <= 1; iAxis++) {
			std::vector<int> ranges = {-2, -1, 1, 2};
			for (int iRange : ranges)
			{
				ActionEdit e = ActionEdit(iV, iAxis, iRange);
				actions.push_back(&e);
			}
		}
	}
	return actions;
}