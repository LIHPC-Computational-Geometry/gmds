//
// Created by Baptiste Bony on 22/04/22.
//

#include <gmds/baptiste/RLBlockSet.h>
#include <gmds/utils/Variable.h>


using namespace gmds;


RLBlockSet::RLBlockSet(Mesh *AMesh)
	:m_mesh(*AMesh){;}

std::vector<double> LinearSpacedArray(double a, double b, std::size_t N)
{
	double h = (b - a) / static_cast<double>(N-1);
	std::vector<double> v(N);
	std::vector<double>::iterator x;
	double val;
	for (x = v.begin(), val = a; x != v.end(); ++x, val += h)
	{
		*x = val;
	}
	return v;
}

void RLBlockSet::createFourNodes(double xA, double yA, double xB, double yB)
{
	std::vector<double> AX = {xA, xB, xA, xB};
   std::vector<double> AY = {yA, yA, yB, yB};
   std::vector<double> AZ = {0, 0, 0, 0};
	Node n = m_mesh.newNode(xA, yA);
}

void RLBlockSet::setFrame(double xMin, double yMin, double xMax, double yMax, int nX, int nY)
{

}

int RLBlockSet::getNumberOfBlocks()
{

}

void RLBlockSet::deleteBlock(const int faceID)
{

}