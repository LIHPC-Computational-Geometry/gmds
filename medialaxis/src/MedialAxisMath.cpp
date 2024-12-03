//
// Created by chenyt on 19/04/24.
//

#include "gmds/medialaxis/MedialAxisMath.h"
#include <algorithm>
   using namespace gmds;

/*----------------------------------------------------------------------------*/
// Computes the maximum distance to a sphere 1 of a point of a sphere 2
double maxRange(const math::Point& ACenter1, const double& ARadius1, const math::Point& ACenter2, const double& ARadius2)
{
	math::Point A = ACenter2 + (-1.)*ACenter1;
	double max_range = vec(A).norm() + ARadius2;
	if (max_range < ARadius1)
		return 0.;
	else
		return max_range - ARadius1;
}

/*----------------------------------------------------------------------------*/
double oriented_angle(const math::Vector &AV1, const math::Vector &AV2)
{
	double cos_theta = AV1.dot(AV2)/(AV1.norm()*AV2.norm());
	if (fabs(cos_theta-1.) < 10e-6)
		return 0.;
	if (fabs(cos_theta+1.) < 10e-6)
		return M_PI;
	math::Vector AV3;
	AV3.setX(-AV1.Y());
	AV3.setY(AV1.X());
	AV3.setZ(0.);
	double sin_theta = AV3.dot(AV2)/(AV3.norm()*AV2.norm());
	double theta = acos(cos_theta);
	if (sin_theta < 0.)
		theta = - theta;
	return theta;
}

/*----------------------------------------------------------------------------*/
math::Vector edge2vec(const gmds::Edge &AE, const gmds::Node &AN)
{
	math::Point E = AE.get<Node>()[1].point() + (-1.)*AE.get<Node>()[0].point();
	if (AE.get<Node>()[0].id() == AN.id())
		return vec(E);
	if (AE.get<Node>()[1].id() == AN.id())
		return vec((-1.)*E);
	std::cout<<"edge2vec : the given node doesn't belong to the given edge"<<std::endl;
	return (vec(E));
}

/*----------------------------------------------------------------------------*/
Node getCommonNode(const Edge &AE1, const Edge &AE2)
{
	if (AE1.get<Node>()[0].id() == AE2.get<Node>()[0].id() || AE1.get<Node>()[0].id() == AE2.get<Node>()[1].id())
		return AE1.get<Node>()[0];
	if (AE1.get<Node>()[1].id() == AE2.get<Node>()[0].id() || AE1.get<Node>()[1].id() == AE2.get<Node>()[1].id())
		return AE1.get<Node>()[1];
	std::cout<<"getCommonNode : the given edges have no node in common"<<std::endl;
	Node n;
	return n;
}

/*----------------------------------------------------------------------------*/
std::vector<Edge> sortEdges(Node &AN, std::vector<Edge> &AV)
{
	// Angles
	std::vector<double> angles;
	math::Point X(1.,0.,0.);
	for (auto e:AV)
	{
		double angle = oriented_angle(vec(X), edge2vec(e,AN));
		angles.push_back(angle);
	}

	std::vector<Edge> sorted_edges;
	double prev_min = -10.;
	for (int i = 1; i <= AV.size(); i++)
	{
		// Find the ith min
		double min_angle = 10.;
		int min_pos;
		for (int j = 0; j < angles.size(); j++)
		{
			if (angles[j] < min_angle && angles[j] > prev_min)
			{
				min_angle = angles[j];
				min_pos = j;
			}
		}
		prev_min = min_angle;
		sorted_edges.push_back(AV[min_pos]);
	}
	return sorted_edges;
}

/*----------------------------------------------------------------------------*/
std::vector<Node> getNextNodes(Node &AN)
{
	std::vector<Node> nxt_nodes;
	for (auto &e:AN.get<Edge>())
	{
		if (e.get<Node>()[0].id() == AN.id())
			nxt_nodes.push_back(e.get<Node>()[1]);
	}
	return nxt_nodes;
}

/*----------------------------------------------------------------------------*/
void propagateValue(Node &AN, Variable<double> &AVar)
{
	std::vector<Node> nxt_nodes = getNextNodes(AN);
	double val1 = AVar.value(AN.id());
	double val2;
	for (auto &n:nxt_nodes)
	{
		val2 = AVar.value(n.id());
		AVar.set(n.id(),val2 + val1/double(nxt_nodes.size()));
	}
}

/*----------------------------------------------------------------------------*/
bool isASource(Node &AN)
{
	for (auto &e:AN.get<Edge>())
	{
		if (e.get<Node>()[1].id() == AN.id())
			return false;
	}
	return true;
}

/*----------------------------------------------------------------------------*/
bool isAWell(Node &AN)
{
	return getNextNodes(AN).empty();
}

/*----------------------------------------------------------------------------*/
int NbWells(Mesh* AMesh)
{
	int N = 0;
	for (auto n_id:AMesh->nodes())
	{
		Node n = AMesh->get<Node>(n_id);
		if (isAWell(n))
			N += 1;
	}
	return N;
}

/*----------------------------------------------------------------------------*/
std::vector<Edge> orientateEdges(Face &AF)
{
	std::vector<Edge> adj_edges = AF.get<Edge>();
	std::vector<Edge> oriented_edges;

	Node n0;
	Edge e1,e2,e_prev,e_nxt;
	double angle;

	// Find a node whose corner isn't flat
	for (auto node:AF.get<Node>())
	{
		for (auto edge:AF.get<Edge>())
		{
			if (edge.get<Node>()[0].id() == node.id() || edge.get<Node>()[1].id() == node.id())
			{
				e1 = edge;
				break;
			}
		}
		for (auto edge:AF.get<Edge>())
		{
			if (edge.get<Node>()[0].id() == node.id() || edge.get<Node>()[1].id() == node.id())
			{
				if (edge.id() != e1.id())
				{
					e2 = edge;
					break;
				}
			}
		}
		angle = oriented_angle(edge2vec(e1,node), edge2vec(e2,node));
		if (fabs(fabs(angle)-M_PI/2.) < M_PI/4)
		{
			n0 = node;
			break;
		}
	}

	if (angle > 0.)
	{
		oriented_edges.push_back(e2);
		oriented_edges.push_back(e1);
		e_prev = e1;
	}
	if (angle < 0.)
	{
		oriented_edges.push_back(e1);
		oriented_edges.push_back(e2);
		e_prev = e2;
	}
	// Erase e1 and e2
	int i1;
	for (int i = 0; i < adj_edges.size(); i++)
	{
		if (e1.id() == adj_edges[i].id())
			i1 = i;
	}
	adj_edges.erase(adj_edges.begin()+i1);
	for (int i = 0; i < adj_edges.size(); i++)
	{
		if (e2.id() == adj_edges[i].id())
			i1 = i;
	}
	adj_edges.erase(adj_edges.begin()+i1);

	// Add the following edges
	while(!adj_edges.empty())
	{
		// Find the edge following e_prev
		for (int i = 0; i < adj_edges.size(); i++)
		{
			e_nxt = adj_edges[i];
			if (e_prev.get<Node>()[0].id() == e_nxt.get<Node>()[0].id())
			{
				e_prev = e_nxt;
				oriented_edges.push_back(e_prev);
				adj_edges.erase(adj_edges.begin()+i);
				break;
			}
			if (e_prev.get<Node>()[0].id() == e_nxt.get<Node>()[1].id())
			{
				e_prev = e_nxt;
				oriented_edges.push_back(e_prev);
				adj_edges.erase(adj_edges.begin()+i);
				break;
			}
			if (e_prev.get<Node>()[1].id() == e_nxt.get<Node>()[0].id())
			{
				e_prev = e_nxt;
				oriented_edges.push_back(e_prev);
				adj_edges.erase(adj_edges.begin()+i);
				break;
			}
			if (e_prev.get<Node>()[1].id() == e_nxt.get<Node>()[1].id())
			{
				e_prev = e_nxt;
				oriented_edges.push_back(e_prev);
				adj_edges.erase(adj_edges.begin()+i);
				break;
			}
		}
	}
	return oriented_edges;
}

/*----------------------------------------------------------------------------*/
std::vector<std::vector<Edge>> groupsOfAlignedEdges(Face &AF)
{
	std::vector<std::vector<Edge>> groups;
	std::vector<Edge> oriented_edges = orientateEdges(AF);
	// Group of the first element
	std::vector<Edge> group0;
	Edge current_edge = oriented_edges[0];
	Edge nxt_edge;
	group0.push_back(current_edge);
	oriented_edges.erase(oriented_edges.begin());
	bool finished = false;
	while (!finished)
	{
		nxt_edge = oriented_edges[oriented_edges.size()-1];
		Node n = getCommonNode(current_edge,nxt_edge);
		double angle = oriented_angle(edge2vec(current_edge,n),edge2vec(nxt_edge,n));
		if (fabs(angle-M_PI) < 0.01 || fabs(angle+M_PI) < 0.01)
		{
			group0.push_back(nxt_edge);
			oriented_edges.erase(oriented_edges.begin()+oriented_edges.size()-1);
			current_edge = nxt_edge;
		}
		else
			finished = true;
	}
	std::reverse(group0.begin(),group0.end());
	groups.push_back(group0);
	// Build the other groups
	while(!oriented_edges.empty())
	{
		std::vector<Edge> newGroup;
		current_edge = oriented_edges[0];
		newGroup.push_back(current_edge);
		oriented_edges.erase(oriented_edges.begin());
		finished = false;
		while (!finished)
		{
			nxt_edge = oriented_edges[0];
			Node n = getCommonNode(current_edge,nxt_edge);
			double angle = oriented_angle(edge2vec(current_edge,n),edge2vec(nxt_edge,n));
			if (fabs(angle-M_PI) < 0.01 || fabs(angle+M_PI) < 0.01)
			{
				newGroup.push_back(nxt_edge);
				oriented_edges.erase(oriented_edges.begin());
				current_edge = nxt_edge;
			}
			else
				finished = true;
		}
		groups.push_back(newGroup);
	}
	return groups;
}

/*----------------------------------------------------------------------------*/
bool isOnSegment(math::Point AP0, math::Point AP1, math::Point AP2)
{
	math::Point U1 = AP1+(-1.)*AP0;
	math::Point U2 = AP2+(-1.)*AP0;
	if (vec(U1).norm() < 10e-6)
		return true;
	if (vec(U2).norm() < 10e-6)
		return true;
	double alpha = oriented_angle(vec(U1),vec(U2));
	return (fabs(alpha-M_PI)<10e-6 || fabs(alpha+M_PI)<10e-6);
}

/*----------------------------------------------------------------------------*/
std::vector<TCellID> insertPoint(Node AN, std::vector<Node> AV)
{
	std::vector<TCellID> new_list;
	math::Point P1,P2;
	bool added = false;
	for (int i = 0; i < AV.size(); i++)
	{
		if (added)
			new_list.push_back(AV[i].id());
		else
		{
			P2 = AV[i].point();
			if (i > 0)
				P1 = AV[i-1].point();
			else 
				P1 = AV[AV.size()-1].point();
			if (isOnSegment(AN.point(),P1,P2))
			{
				new_list.push_back(AN.id());
				added = true;
				new_list.push_back(AV[i].id());
			}
			else
				new_list.push_back(AV[i].id());				 
		}
	}
	return new_list;
}

/*----------------------------------------------------------------------------*/
std::vector<math::Point> merge(std::vector<math::Point> AV1, std::vector<math::Point> AV2)
{
	std::vector<math::Point> merged;
	for (auto point:AV1)
		merged.push_back(point);
	for (auto point:AV2)
	{
		bool already_seen = false;
		for (auto point2:AV1)
		{
			if (point.distance(point2) < 10e-6)
			{
				already_seen = true;
				break;
			}
		}
		if (!already_seen)
			merged.push_back(point);
	}
	return merged;
}

/*----------------------------------------------------------------------------*/
std::vector<Edge> order(std::vector<Edge> AVE, std::vector<double> AVX)
{
	std::vector<Edge> sorted_edges;
	while (!AVX.empty())
	{
		int min_pos;
		double min_x = 10e6;
		double x;
		for (int i = 0 ; i < AVX.size() ; i++)
		{
			x = AVX[i];
			if (x < min_x)
			{
				min_x =x;
				min_pos = i;
			}
		}
		sorted_edges.push_back(AVE[min_pos]);
		AVX.erase(AVX.begin() + min_pos);
		AVE.erase(AVE.begin() + min_pos);
	}
	return sorted_edges;
}

/*----------------------------------------------------------------------------*/
Edge getEdge(Node &AN1, Node &AN2)
{
	Edge e;
	for (auto e1:AN1.get<Edge>())
	{
		for (auto e2:AN2.get<Edge>())
		{
			if (e1.id() == e2.id())
			{
				e = e1;
				return e;
			}
		}
	}
	throw GMDSException("getEdge() : the two given nodes have no edge in common");
}

/*----------------------------------------------------------------------------*/
bool isInterior(Node &AN)
{
	bool isInt = true;
	for (auto e:AN.get<Edge>())
	{
		if (e.get<Face>().size() == 1)
		{
			isInt = false;
			break;
		}
	}
	return isInt;
}

/*----------------------------------------------------------------------------*/
double delimitedArea(std::vector<Edge> AV)
{
	double x = 0.;
	double y = 0.;
	// Compute the barycenter of the points
	for (auto e:AV)
	{
		for (auto n:e.get<Node>())
		{
			x += 1./2.*n.X();
			y += 1./2.*n.Y();
		}
	}
	x = x/double(AV.size());
	y = y/double(AV.size());
	math::Point Bar(x,y,0.);
	// Compute the area
	double area = 0.;
	for (auto e:AV)
	{
		math::Triangle T(e.get<Node>()[0].point(),e.get<Node>()[1].point(),Bar);
		area += T.area();
	}
	return area;
}

/*----------------------------------------------------------------------------*/
bool touchesBoundary(Edge &AE)
{
	return (!isInterior(AE.get<Node>()[0]) || !isInterior(AE.get<Node>()[1]));
}

/*----------------------------------------------------------------------------*/
Face getCommonFace(Edge &AE1, Edge &AE2)
{
	for (auto f1:AE1.get<Face>())
	{
		for (auto f2:AE2.get<Face>())
		{
			if (f1.id() == f2.id())
			{
				return f1;
			}
		}
	}
	throw GMDSException("getCommonFace() : the two input edges have no face in common");
}