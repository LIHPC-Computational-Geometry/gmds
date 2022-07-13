//
// Created by Baptiste Bony on 22/04/22.
//

#include "gmds/baptiste/tools.h"
#include "gmds/igalgo/VolFracComputation.h"
#include "gmds/igalgo/r2d.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKWriter.h"
#include "gmds/math/Quadrilateral.h"
#include <gmds/baptiste/RLBlockSet.h>
#include <gmds/io/VTKReader.h>

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


Node RLBlockSet::findCorner(std::vector<Node> nodes, bool v)
{
	/*
	int nodeID = 0;
	Node corner = m_mesh.get<Node>(nodeID);
	std::cout << "X : " << corner.X() << ", Y : " << corner.Y() << "\n";
	std::cout << nodes.size();
	 */
	Node corner = nodes[0];
	for (auto node : nodes)
	{
		if (v)
		{
			if (node.X() <= corner.X() and  node.Y() <= corner.Y())
			{
				corner = node;
			}
		}
		else
		{
			if (node.X() >= corner.X() and  node.Y() >= corner.Y())
			{
				corner = node;
			}
		}
	}
	std::cout << "X : " << corner.X() << ", Y : " << corner.Y() << "\n";
	return corner;
}

Node RLBlockSet::findSecondCorner(Node corner, std::vector<Node> nodes, std::string axis)
{
	Node secondCorner = nodes[0];
	std::cout << "X : " << secondCorner.X() << ", Y : " << secondCorner.Y() << "\n";
	for (Node node : nodes)
	{
		secondCorner = node;
	}
	for (Node node : nodes)
	{
		if (axis == "x")
		{
			if (node.X() != corner.X() and  node.Y() == corner.Y())
			{
				secondCorner = node;
			}
		}
		else if (axis == "y")
		{
			if (node.X() == corner.X() and  node.Y() != corner.Y())
			{
				secondCorner = node;
			}
		}
	}
	std::cout << "X : " << secondCorner.X() << ", Y : " << secondCorner.Y() << "\n";
	return secondCorner;
}

void RLBlockSet::moveCorners(Node corner, Node secondCorner, std::string axis, int range)
{
	if (axis == "x")
	{
		corner.setX(corner.X()+range*xSize/10);
		secondCorner.setX(secondCorner.X()+range*xSize/10);
	}
	else if (axis == "y")
	{
		corner.setY(corner.Y()+range*ySize/10);
		secondCorner.setY(secondCorner.Y()+range*ySize/10);
	}
}

void RLBlockSet::editCornerBis(const int faceID, bool v, std::string axis, int range)
{
	Face face = m_mesh.get<Face>(faceID);
	std::cout << "xxxxxx " << face << "\n";
	std::vector<Node> nodes = face.get<Node>();
	std::cout << "$$$$$$$$$$$$$$$$$$$$ " << nodes.empty() << "\n";
	Node corner = findCorner(nodes, v);
	Node secondCorner = findSecondCorner(corner, nodes, axis);
	moveCorners(corner, secondCorner, axis,  range);
}

void RLBlockSet::saveMesh(std::string filename)
{
	IGMeshIOService ioService(&m_mesh);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write(filename);
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
	if (not targetMesh.hasVariable(GMDS_FACE, "volFrac"))
	{
		targetMesh.newVariable<double, GMDS_FACE>("volFrac");
	}
	computeVolFrac(&m_mesh, &targetMesh, m_mesh.getVariable<double, GMDS_FACE>("volFrac"), targetMesh.getVariable<double, GMDS_FACE>("volFrac"));

	Variable<double>* volFracVar = m_mesh.getVariable<double, GMDS_FACE>("volFrac");
	double a = 0;
    /*
	for (int i =0; i < volFracVar->getNbValues(); i++)
	{
		a += volFracVar->value(i);
	}
	a = a/countBlocks();
     */


	for (int faceID : m_mesh.faces())
	{
		Face face = m_mesh.get<Face>(faceID);
		a += face.area() * volFracVar->value(faceID);
	}
	a = a/getMeshArea(m_mesh);

	Variable<double>* volFracVarReverse = targetMesh.getVariable<double, GMDS_FACE>("volFrac");
	double b = 0;
    /*
	for (int i =0; i < volFracVarReverse->getNbValues(); i++)
	{
		b += volFracVarReverse->value(i);
	}
	b = b/targetMesh.getNbFaces();
     */

	for (int faceID : targetMesh.faces())
	{
		Face face = targetMesh.get<Face>(faceID);
		b += face.area() * volFracVarReverse->value(faceID);
	}
	b = b/getMeshArea(targetMesh);


	/*
	if (not m_mesh.hasVariable(GMDS_FACE, "volFrac"))
	{
		m_mesh.newVariable<double, GMDS_FACE>("volFrac");
	}
	volfraccomputation_2d(&m_mesh, &targetMesh, m_mesh.getVariable<double, GMDS_FACE>("volFrac"));
	Variable<double>* volFracVar = m_mesh.getVariable<double, GMDS_FACE>("volFrac");
	double a = 0;
	for (int i =0; i < volFracVar->getNbValues(); i++)
	{
		a += volFracVar->value(i);
	}
	a = a/countBlocks();

	if (not targetMesh.hasVariable(GMDS_FACE, "volFrac2"))
	{
		targetMesh.newVariable<double, GMDS_FACE>("volFrac2");
	}
	anotherVolFrac(&m_mesh, &targetMesh, targetMesh.getVariable<double, GMDS_FACE>("volFrac2"));
	Variable<double>* volFracVarReverse = targetMesh.getVariable<double, GMDS_FACE>("volFrac2");
	double b = 0;
	for (int i =0; i < volFracVarReverse->getNbValues(); i++)
	{
		b += volFracVarReverse->value(i);
	}
	b = b/targetMesh.getNbFaces();
	 */

	double alpha = 0.25;
	double beta = 4;
	double c = getOverlap();
    /*
    if (a < 0)
    {
        a = 0;
    }
    if (b < 0)
    {
        b = 0;
    }
    if (c < 0)
    {
        c = 0;
    }
    if (alpha * a + (1-alpha) * b - beta * c < 0)
    {
        return 0;
    }
     */
    std::cout << "Quad mesh volFrac : " << a << "\n";
    std::cout << "Target mesh volFrac : " << b << "\n";
    std::cout << "Overlap : " << c << "\n";
    double res = alpha * a + (1-alpha) * b - beta * c;
    res = std::fmax(0, res);
    return res;
}

bool RLBlockSet::isValid(Mesh &targetMesh)
{
    std::vector<Node> nodes;
    for (int faceID : targetMesh.faces())
    {
        Face face = targetMesh.get<Face>(faceID);
        std::vector<Node> nodesF = face.get<Node>();
        nodes.insert( nodes.end(), nodesF.begin(), nodesF.end() );
    }
    std::vector<Node>::iterator iter;
    iter = std::min_element(nodes.begin(), nodes.end(),[](Node a, Node b){return a.X() <= b.X();});
    double xMin = targetMesh.get<Node>(iter->id()).X();
    iter = std::min_element(nodes.begin(), nodes.end(),[](Node a, Node b){return a.X() >= b.X();});
    double xMax = targetMesh.get<Node>(iter->id()).X();
    iter = std::min_element(nodes.begin(), nodes.end(),[](Node a, Node b){return a.Y() <= b.Y();});
    double yMin = targetMesh.get<Node>(iter->id()).Y();
    iter = std::min_element(nodes.begin(), nodes.end(),[](Node a, Node b){return a.Y() >= b.Y();});
    double yMax = targetMesh.get<Node>(iter->id()).Y();
	bool res = true;
	for(auto faceID: m_mesh.faces())
	{
		gmds::Face f = m_mesh.get<Face>(faceID);
		std::vector<gmds::Node> n = f.get<gmds::Node>();

		gmds::math::Quadrilateral quad(n[0].point(), n[1].point(), n[2].point(), n[3].point());
		double sj = quad.computeScaledJacobian2D();
		if(sj < 0)
		{
			res = false;
			break;
		}
        if (f.area() <= 0)
        {
            res = false;
            break;
        }
        std::vector<Node> nodes = f.get<Node>();
        for (Node node : nodes)
        {
            if (node.X() > xMax or node.X() < xMin)
            {
                res = false;
                break;
            }
            if (node.Y() > yMax or node.Y() < yMin)
            {
                res = false;
                break;
            }
        }
	}
	return res;
}

std::string RLBlockSet::getStateID()
{
	/*
	std::ostringstream out;
	m_mesh.serialize(out);
	return out.str();
	 */
	std::string res;
	for (int nodeID : m_mesh.nodes())
	{
		Node node = m_mesh.get<Node>(nodeID);
		res += std::to_string(node.X());
		res += std::to_string(node.Y());
	}
	return res;
}

double RLBlockSet::overlap()
{
	double res = 0;
	for (int faceIDA : m_mesh.faces())
	{
		Face faceA = m_mesh.get<Face>(faceIDA);
		for (int faceIDB : m_mesh.faces())
		{
			Face faceB = m_mesh.get<Face>(faceIDB);
			if (faceIDA != faceIDB)
			{
				std::vector<Node> nodesA = faceA.get<Node>();
				std::vector<Node> nodesB = faceB.get<Node>();
				r2d_rvec2 vertices[4];
				vertices[0].x = nodesB[0].X();
				vertices[0].y = nodesB[0].Y();
				vertices[1].x = nodesB[1].X();
				vertices[1].y = nodesB[1].Y();
				vertices[2].x = nodesB[2].X();
				vertices[2].y = nodesB[2].Y();
				vertices[3].x = nodesB[3].X();
				vertices[3].y = nodesB[3].Y();

				r2d_plane planes[4];
				r2d_poly_faces_from_verts(planes,  vertices, 4);

				r2d_poly poly;
				r2d_rvec2 verts[4];

				verts[0].x = nodesA[0].X();
				verts[0].y = nodesA[0].Y();
				verts[1].x = nodesA[1].X();
				verts[1].y = nodesA[1].Y();
				verts[2].x = nodesA[2].X();
				verts[2].y = nodesA[2].Y();
				verts[3].x = nodesA[3].X();
				verts[3].y = nodesA[3].Y();

				r2d_init_poly(&poly, verts, 4);

				r2d_clip(&poly, planes, 4);
				r2d_int POLY_ORDER = 2;
				r2d_real om[R2D_NUM_MOMENTS(POLY_ORDER)];
				r2d_reduce(&poly, om, POLY_ORDER);

				res += om[0];
			}
		}
		res = res * faceA.area();
	}
	// return res / (xSize * ySize * countBlocks());
    double result;
    if (res / getMeshArea(m_mesh) > 1)
    {
        result = 0;
    }
    else
    {
        result = res / getMeshArea(m_mesh);
    }
    return res / getMeshArea(m_mesh);
}

double RLBlockSet::getOverlap()
{
    double res = 0;
    for (int faceID1 : m_mesh.faces())
    {
        for (int faceID2: m_mesh.faces())
        {
            if (faceID1 != faceID2)
            {
                Face face1 = m_mesh.get<Face>(faceID1);


                std::vector<Node> nodes1 = face1.get<Node>();
                std::vector<Node>::iterator iter1;
                iter1 = std::min_element(nodes1.begin(), nodes1.end(),[](Node a, Node b){return a.X() <= b.X() and a.Y() <= b.Y();});
                int nodeID1 = iter1->id();
                Node node1 = m_mesh.get<Node>(nodeID1);


                Face face2 = m_mesh.get<Face>(faceID2);
                std::vector<Node> nodes2 = face2.get<Node>();
                std::vector<Node>::iterator iter2;
                iter2 = std::min_element(nodes2.begin(), nodes2.end(),[](Node a, Node b){return a.X() >= b.X() and a.Y() >= b.Y();});
                int nodeID2 = iter2->id();
                Node node2 = m_mesh.get<Node>(nodeID2);

                double intersection = 0;
                if (node2.Y() > node1.Y() and node2.X() > node1.X())
                {
                    intersection = (node2.Y() - node1.Y()) * (node2.X() - node1.X());
                }
                if (intersection > face1.area() or intersection > face2.area())
                {
                    intersection = 0;
                }
                res += intersection;
            }
        }
    }
    return res/getMeshArea(m_mesh);
}

double RLBlockSet::getLocalIou(int faceID)
{
    return m_mesh.getVariable<double,GMDS_FACE>("volFrac")->value(faceID);
}

std::vector<double> RLBlockSet::getMinMaxCoordinates()
{
    std::vector<Node> nodes;
    for (int faceID : m_mesh.faces())
    {
        Face face = m_mesh.get<Face>(faceID);
        std::vector<Node> nodesF = face.get<Node>();
        nodes.insert( nodes.end(), nodesF.begin(), nodesF.end() );
    }
    std::vector<double> res;
    std::vector<Node>::iterator iter;
    iter = std::min_element(nodes.begin(), nodes.end(),[](Node a, Node b){return a.X() <= b.X();});
    res.push_back(m_mesh.get<Node>(iter->id()).X());
    iter = std::min_element(nodes.begin(), nodes.end(),[](Node a, Node b){return a.X() >= b.X();});
    res.push_back(m_mesh.get<Node>(iter->id()).X());
    iter = std::min_element(nodes.begin(), nodes.end(),[](Node a, Node b){return a.Y() <= b.Y();});
    res.push_back(m_mesh.get<Node>(iter->id()).Y());
    iter = std::min_element(nodes.begin(), nodes.end(),[](Node a, Node b){return a.Y() >= b.Y();});
    res.push_back(m_mesh.get<Node>(iter->id()).Y());
    return res;
}