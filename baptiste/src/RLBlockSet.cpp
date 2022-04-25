//
// Created by Baptiste Bony on 22/04/22.
//

#include </home/bonyb/Documents/GitHub/gmds/baptiste/inc/gmds/baptiste/RLBlockSet.h>
#include <gmds/utils/Variable.h>

using namespace gmds;

RLBlockSet::RLBlockSet(gmds::Mesh mesh)
{
    m_mesh = mesh;
}

void createFourNodes(double xA, double yA, double xB, double yB);
{
    std::vector<double> AX = {xA, xB, xA, xB};
    std::vector<double> AY = {yA, yA, yB, yB};
    std::vector<double> AZ = {0, 0, 0, 0};
    m_mesh.createNodes(AX, AY, AZ);
}

void setFrame(double xMin, double yMin, double xMax, double yMax)
{
    double xPitch = abs(xMax - xMin)/3;
    double yPitch = abs(yMax - yMin)/3;
    double x1 = xMin + xPitch;
    double x2 = xMin + 2*xPitch;
    double y1 = yMin + yPitch;
    double y2 = yMin + 2*yPitch;
    createFourNodes(xMin, yMin, x1, y1);
    createFourNodes(x1, yMin, x2, y1);
    createFourNodes(x2, yMin, xMax, y1);
    createFourNodes(xMin, y1, x1, y2);
    createFourNodes(x1, y1, x2, y2);
    createFourNodes(x2, y1, xMax, y2);
    createFourNodes(xMin, y2, x1, yMax);
    createFourNodes(x1, y2, x2, yMax);
    createFourNodes(x2, y2, xMax, yMax);
}

int getNumberOfBlocks()
{

}

void deleteBlock(const int faceID)
{

}