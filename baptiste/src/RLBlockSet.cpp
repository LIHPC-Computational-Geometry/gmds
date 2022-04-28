//
// Created by Baptiste Bony on 22/04/22.
//

#include <gmds/baptiste/RLBlockSet.h>
using namespace gmds;


RLBlockSet::RLBlockSet(MeshModel model)
	:m_mesh(model) {}

RLBlockSet::~RLBlockSet()
{}


std::vector<double> RLBlockSet::LinearSpacedArray(double a, double b, std::size_t N)
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

void RLBlockSet::setFrame(double xMin, double yMin, double xMax, double yMax, int nX, int nY)
{
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

int RLBlockSet::getNumberOfBlocks()
{
	int numberOfBlocks = m_mesh.getNbFaces();
	return numberOfBlocks;
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
	}
	else
	{
		iter = std::max_element(nodes.begin(), nodes.end(),[](Node a, Node b){return a.X() <= b.X() and a.Y() <= b.Y();});
	}
	int nodeID = iter->id();
	Node corner = m_mesh.get<Node>(nodeID);
	if (axis == "x")
	{
		corner.setX(corner.X()+range);
	}
	else if (axis == "y")
	{
		corner.setY(corner.Y()+range);
	}
}
