//
// Created by Baptiste Bony on 22/04/22.
//

#include <gmds/baptiste/RLBlockSet.h>
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKWriter.h"
#include <gmds/io/VTKReader.h>
#include "gmds/baptiste/tools.h"
#include "gmds/igalgo/VolFracComputation.h"

using namespace gmds;


RLBlockSet::RLBlockSet(MeshModel model)
	:m_mesh(model) {}

RLBlockSet::RLBlockSet()
  :m_mesh(MeshModel(DIM2|F|N|F2N|N2F)) {}

RLBlockSet::~RLBlockSet() {}


void RLBlockSet::setFrame(double xMin, double yMin, double xMax, double yMax, int nX, int nY)
{
	xSize = (xMax - xMin) / double(nX);
	ySize = (yMax - yMin) / double(nY);
	std::vector<double> xVector = LinearSpacedArray(xMin, xMax, nX);
	std::vector<double> yVector = LinearSpacedArray(yMin, yMax, nY);
	for (int i = 0; i < xVector.size()-1; ++i)
	{
		for (int j = 0; j < yVector.size()-1; ++j)
		{
			Node n0 = m_mesh.newNode(xVector[i], yVector[j]);
			Node n1 = m_mesh.newNode(xVector[i+1], yVector[j]);
			Node n2 = m_mesh.newNode(xVector[i+1], yVector[j+1]);
			Node n3 = m_mesh.newNode(xVector[i], yVector[j+1]);
			m_mesh.newQuad(n0, n1, n2, n3);
		}
	}
}

int RLBlockSet::countBlocks()
{
	return m_mesh.getNbFaces();
}

void RLBlockSet::deleteBlock(const int faceID)
{
	m_mesh.deleteFace(faceID);
}

void RLBlockSet::editCorner(const int faceID, bool v, std::string axis, int range)
{
	Face face = m_mesh.get<Face>(faceID);
	std::vector<Node> nodes = face.get<Node>();
	std::vector<Node>::iterator iter;
	if (v)
	{
		iter = std::min_element(nodes.begin(), nodes.end(),[](Node a, Node b){return a.X() <= b.X() and a.Y() <= b.Y();});
		int nodeID = iter->id();
		Node corner = m_mesh.get<Node>(nodeID);
		std::vector<Node>::iterator iterBis;
		if (axis == "x")
		{
			iterBis = std::min_element(nodes.begin(), nodes.end(),[](Node a, Node b){return a.X() <= b.X() and a.Y() >= b.Y();});
			int nodeIDBis = iterBis->id();
			Node cornerBis = m_mesh.get<Node>(nodeIDBis);
			corner.setX(corner.X()+range*xSize/10);
			cornerBis.setX(cornerBis.X()+range*xSize/10);
		}
		else if (axis == "y")
		{
			iterBis = std::min_element(nodes.begin(), nodes.end(),[](Node a, Node b){return a.X() >= b.X() and a.Y() <= b.Y();});
			int nodeIDBis = iterBis->id();
			Node cornerBis = m_mesh.get<Node>(nodeIDBis);
			corner.setY(corner.Y()+range*ySize/10);
			cornerBis.setY(cornerBis.Y()+range*ySize/10);
		}
	}
	else
	{
		iter = std::min_element(nodes.begin(), nodes.end(),[](Node a, Node b){return a.X() >= b.X() and a.Y() >= b.Y();});
		int nodeID = iter->id();
		Node corner = m_mesh.get<Node>(nodeID);
		std::vector<Node>::iterator iterBis;
		if (axis == "x")
		{
			iterBis = std::min_element(nodes.begin(), nodes.end(),[](Node a, Node b){return a.X() >= b.X() and a.Y() <= b.Y();});
			int nodeIDBis = iterBis->id();
			Node cornerBis = m_mesh.get<Node>(nodeIDBis);
			corner.setX(corner.X()+range*xSize/10);
			cornerBis.setX(cornerBis.X()+range*xSize/10);
		}
		else if (axis == "y")
		{
			iterBis = std::min_element(nodes.begin(), nodes.end(),[](Node a, Node b){return a.X() <= b.X() and a.Y() >= b.Y();});
			int nodeIDBis = iterBis->id();
			Node cornerBis = m_mesh.get<Node>(nodeIDBis);
			corner.setY(corner.Y()+range*ySize/10);
			cornerBis.setY(cornerBis.Y()+range*ySize/10);
		}
	}
}

void RLBlockSet::saveMesh(std::string title)
{
	IGMeshIOService ioService(&m_mesh);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write(title + ".vtk");
}

void RLBlockSet::setFromFile(std::string filename, int nX, int nY)
{
	Mesh targetMesh = Mesh(MeshModel(DIM3|F|N));
	gmds::IGMeshIOService ioService(&targetMesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(filename);

	double xMin, yMin, xMax, yMax;
	for (int nodeID : targetMesh.nodes())
	{
		Node node = targetMesh.get<Node>(nodeID);
		if (node.X() > xMax)
		{
			xMax = node.X();
		}
		if (node.X() < xMin)
		{
			xMin = node.X();
		}
		if (node.Y() > yMax)
		{
			yMax = node.Y();
		}
		if (node.Y() < yMin)
		{
			yMin = node.Y();
		}
	}
	setFrame(xMin, yMin, xMax, yMax, nX, nY);
}

std::vector<int> RLBlockSet::getAllFaces()
{
	std::vector<int> v;
	for (int faceID : m_mesh.faces())
	{
		v.push_back(faceID);
	}
	return v;
}

double RLBlockSet::getReward(Mesh &targetMesh)
{
	if (not m_mesh.hasVariable(GMDS_FACE, "volFrac"))
	{
		m_mesh.newVariable<double, GMDS_FACE>("volFrac");
	}

	volfraccomputation_2d(&m_mesh, &targetMesh, m_mesh.getVariable<double, GMDS_FACE>("volFrac"));
	Variable<double>* res = m_mesh.getVariable<double, GMDS_FACE>("volFrac");
	double sum = 0;
	for (int i =0; i < res->getNbValues(); i++)
	{
		sum += res->value(i);
	}
	return sum;
}