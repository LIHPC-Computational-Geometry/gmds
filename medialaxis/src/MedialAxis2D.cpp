/*----------------------------------------------------------------------------*/
#include "gmds/medialaxis/MedialAxis2D.h"
#include "gmds/ig/MeshDoctor.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKWriter.h"
#include "gmds/medialaxis/MedialAxisMath.h"
#include <stack>
#include <queue>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace medialaxis {
/*----------------------------------------------------------------------------*/
MedialAxis2D::MedialAxis2D(){
	// Geometrical representation
	m_mesh_representation = new Mesh(MeshModel(DIM3 | F | E | N |
	                             F2N | E2N | N2E | N2F));
	// Dual medial point/primal triangle correspondence
	m_mesh_representation->newVariable<int,GMDS_NODE>("MP2triangle");
	// Medial points type: tangency of their associated triangle (1: end point, 2: regular point, 3: intersection point)
	m_mesh_representation->newVariable<int,GMDS_NODE>("medial_point_type");
	// Medial edges type: (1: end edge, 2: regular edge, 3: intersection edge)
	m_mesh_representation->newVariable<int,GMDS_EDGE>("medial_edge_type");
	// Medial points'branch type: 0 if on an interior branch, 1 if on an adjacent to the boundary branch, 2 if on several branches
	m_mesh_representation->newVariable<int,GMDS_NODE>("point_branch_type");
	// Medial edges'branch type: 0 if in an interior branch, 1 if in an adjacent to the boundary branch
	m_mesh_representation->newVariable<int,GMDS_EDGE>("edge_branch_type");
	// Medial edges'branch id
	m_mesh_representation->newVariable<int,GMDS_EDGE>("edge_branch_id");
	// Medial points'branch id
	m_mesh_representation->newVariable<int,GMDS_NODE>("point_branch_id");
	// Mark with 1 the edges whose branch corresponds to a detail, 0 the other edges
	m_mesh_representation->newVariable<int,GMDS_EDGE>("corresponds_to_details");
	// Mark with 1 the points whose branch corresponds to a detail, 0 the other points
	m_mesh_representation->newVariable<int,GMDS_NODE>("point_corresponds_to_details");
	// Cosine of the medial angles
	m_mesh_representation->newVariable<double,GMDS_NODE>("cos_medial_angle");
	// Medial radius
	m_mesh_representation->newVariable<double,GMDS_NODE>("medial_radius");
	// Medial radius orthogonality default
	m_mesh_representation->newVariable<double,GMDS_NODE>("medial_radius_orthogonality_default");
	// Locate the singularities (1 if there is a singularity at the point, 0 if not)
	m_mesh_representation->newVariable<int,GMDS_NODE>("singularity");
	// Optimum change in angle of a cross following a curve formed by two adjacent medial radii
	m_mesh_representation->newVariable<double,GMDS_NODE>("flux_through_medial_radii");
	// Flux residual across edges
	m_mesh_representation->newVariable<double,GMDS_EDGE>("flux_residual");
	// Flux residual at intersection points
	m_mesh_representation->newVariable<double,GMDS_NODE>("flux_residual_at_intersection_points");
	// Boundary points at which each medial point touches the boundary
	m_mesh_representation->newVariable<std::vector<math::Point>,GMDS_NODE>("touching_points");
	// List of triangles dual to the medial point
	m_mesh_representation->newVariable<std::vector<Face>,GMDS_NODE>("dual_triangles");
	// List of tengancy points of the medial points as nodes of the minimal triangulation 
	m_mesh_representation->newVariable<std::vector<Node>,GMDS_NODE>("tangent_boundary_nodes");
	// Medial point/group correspondence
	m_mesh_representation->newVariable<int,GMDS_NODE>("medialPoint2Group");
	// List of boundary connected components connected by this medial point
	m_mesh_representation->newVariable<std::vector<int>,GMDS_NODE>("boundary_connected_components_ids");
	// Mark with 1 dangles (end points corresponding to boundary circle arcs)
	m_mesh_representation->newVariable<int,GMDS_NODE>("is_a_dangle");
	// Mark with 1 medial points belonging to extensions (corresponding to boundary circle arcs)
	m_mesh_representation->newVariable<int,GMDS_NODE>("is_an_extension");
}
MedialAxis2D::~MedialAxis2D()
{
	if(m_mesh_representation!= nullptr)
		delete m_mesh_representation;
}

/*----------------------------------------------------------------------------*/
Mesh MedialAxis2D::getMeshRepresentation()
{
	return *m_mesh_representation;
}
/*----------------------------------------------------------------------------*/
Node
MedialAxis2D::newMedPoint(const math::Point &APnt)
{
	Node n = m_mesh_representation->newNode(APnt);
	std::vector<Node> nv = {n, n, n,};
	m_mesh_representation->newFace(nv);
	return n;
}

/*----------------------------------------------------------------------------*/
Node MedialAxis2D::getMedPoint(const gmds::TCellID APointID)
{
	return m_mesh_representation->get<Node>(APointID);
}

/*----------------------------------------------------------------------------*/
Edge
MedialAxis2D::newMedEdge(const TCellID &AN1, const TCellID &AN2)
{
	Edge e = m_mesh_representation->newEdge(AN1, AN2);
/*	auto var = m_mesh_representation->getVariable<int,GMDS_NODE>("toto");
	var->value(2);
	var->set(2,1234);
	*/
	return e;
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::setSingularities()
{
	auto sing = m_mesh_representation-> getVariable<int,GMDS_NODE>("singularity");
	for (auto n_id:m_mesh_representation->nodes())
	{
		if (sing->value(n_id) !=0)
		{
			Node n = m_mesh_representation->get<Node>(n_id);
			double index = (4.-double(4+sing->value(n_id)))/4.;
			m_singular_nodes.push_back(n);
			m_singularity_indexes.push_back(index);
		}
	}
}

/*----------------------------------------------------------------------------*/
std::vector<Node> MedialAxis2D::getSingularNodes()
{
	return m_singular_nodes;
}

/*----------------------------------------------------------------------------*/
std::vector<double> MedialAxis2D::getSingularityIndexes()
{
	return m_singularity_indexes;
}

/*----------------------------------------------------------------------------*/
double MedialAxis2D::minMedEdgeLength()
{
	double min_len = 10e5;
	for (auto e_id:m_mesh_representation->edges())
	{
		Edge e = m_mesh_representation->get<Edge>(e_id);
		if (e.length()<min_len)
			min_len = e.length();
	}
	return min_len;
}

/*----------------------------------------------------------------------------*/
double MedialAxis2D::maxMedEdgeLength()
{
	double max_len = 0.;
	for (auto e_id:m_mesh_representation->edges())
	{
		Edge e = m_mesh_representation->get<Edge>(e_id);
		if (e.length()>max_len)
			max_len = e.length();
	}
	return max_len;
}

/*----------------------------------------------------------------------------*/
double MedialAxis2D::meanMedEdgeLength()
{
	double mean_len = 0.;
	for (auto e_id:m_mesh_representation->edges())
	{
		Edge e = m_mesh_representation->get<Edge>(e_id);
		mean_len += e.length();
	}
	return mean_len/double(m_mesh_representation->getNbEdges());
}

/*----------------------------------------------------------------------------*/
int MedialAxis2D::getNbMedPoints()
{
	return m_mesh_representation->getNbNodes();
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::updateConnectivity()
{
	std::cout<<"> Updating medial axis connectivity..."<<std::endl;
	MeshDoctor doc(m_mesh_representation);
	doc.updateUpwardConnectivity();
	std::cout<<"... done. Medial axis topology may be used."<<std::endl;
}

/*----------------------------------------------------------------------------*/
std::vector<std::vector<Node>> MedialAxis2D::medialPointsGroups(const double &ATol)
{
	auto medPointType = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	// Medial point/group correspondence
	auto medialPoint2Group = m_mesh_representation->getVariable<int,GMDS_NODE>("medialPoint2Group");
	// Associate to each node the value 1 if already visited
	auto alreadyVisited = m_mesh_representation->newVariable<int,GMDS_NODE>("already_visited");
	for (auto n_id:m_mesh_representation->nodes())
		alreadyVisited->set(n_id,0);
	// Set of medial points groups
	std::vector<std::vector<Node>> groups;
	int groupID = 0;
	// Build the groups
	for (auto n_id:m_mesh_representation->nodes())
	{
		if (alreadyVisited->value(n_id) == 0)
		{
			std::vector<Node> newGroup;
			std::stack<Node> newGroupMembers;
			Node n = m_mesh_representation->get<Node>(n_id);
			newGroupMembers.push(n);
			alreadyVisited->set(n_id,1);
			while (!newGroupMembers.empty())
			{
				Node m = newGroupMembers.top();
				newGroupMembers.pop();
				newGroup.push_back(m);
				medialPoint2Group->set(m.id(),groupID);
				for (auto e:m.get<Edge>())
				{
					Node n1 = e.get<Node>()[0];
					Node n2 = e.get<Node>()[1];
					if (e.length() < ATol)
					{
						for (auto m1:e.get<Node>())
						{
							if (alreadyVisited->value(m1.id()) == 0)
							{
								newGroupMembers.push(m1);
								alreadyVisited->set(m1.id(),1);
							}
						}
					}
				}
			}
			groups.push_back(newGroup);
			groupID += 1;
		}
	}
	m_mesh_representation->deleteVariable(GMDS_NODE,alreadyVisited);
	return groups;
}

/*----------------------------------------------------------------------------*/
int MedialAxis2D::medialPoint2Group(const gmds::Node &ANode)
{
	auto medialPoint2Group = m_mesh_representation->getVariable<int,GMDS_NODE>("medialPoint2Group");
	return medialPoint2Group->value(ANode.id());
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::refine(const double& ATol)
{
	auto tangentNodes = m_mesh_representation->getVariable<std::vector<Node>,GMDS_NODE>("tangent_boundary_nodes");
	auto touchingPoints = m_mesh_representation->getVariable<std::vector<math::Point>,GMDS_NODE>("touching_points");
	auto medRadius = m_mesh_representation->getVariable<double,GMDS_NODE>("medial_radius");
	auto dualTriangles = m_mesh_representation->getVariable<std::vector<Face>,GMDS_NODE>("dual_triangles");
	auto iae = m_mesh_representation->getVariable<int,GMDS_NODE>("is_an_extension");
	std::cout<<"> Refining medial axis"<<std::endl;
	double min_len = 10e5;
	double max_len = 0.;
	double mean_len = 0.;
	for (auto e_id:m_mesh_representation->edges())
	{
		double len = m_mesh_representation->get<Edge>(e_id).length();
		if (len < min_len)
			min_len = len;
		if (len > max_len)
			max_len = len;
		mean_len += len;
	}
	mean_len = mean_len/m_mesh_representation->getNbEdges();
	std::cout<<"Min edge length : "<<min_len<<"; Max edge length : "<<max_len<<"; Mean edge length : "<<mean_len<<std::endl;
	// Adding new medial points where they lack
	int NbPointsAdded = 0;
	int NbEdgesAdded = 0;
	int NbEdgesDeleted = 0;
	for (auto e_id:m_mesh_representation->edges())
	{
		Edge e = m_mesh_representation->get<Edge>(e_id);
		Node n1 = e.get<Node>()[0];
		Node n2 = e.get<Node>()[1];
		math::Point direction = n2.point()+(-1.)*n1.point();
		direction = direction*(1./vec(direction).norm());
		if ((touchingPoints->value(n1.id()).size() != 2 && touchingPoints->value(n2.id()).size() != 2))
		{
			if (touchingPoints->value(n1.id()).size() == 1 && touchingPoints->value(n2.id()).size() == 1 && e.length() > ATol)
				std::cout<<"Could not refine an edge because it is incident to two EPs."<<std::endl;
			if (touchingPoints->value(n1.id()).size() == 1 && touchingPoints->value(n2.id()).size() > 2 || touchingPoints->value(n1.id()).size() > 2 && touchingPoints->value(n2.id()).size() == 1 && e.length() > ATol)
				std::cout<<"Could not refine an edge because it is incident to an EP and an IP."<<std::endl;
			if (touchingPoints->value(n1.id()).size() > 2 && touchingPoints->value(n2.id()).size() > 2)
			{
				math::Point M = (1./2.)*(n1.point()+n2.point());
				Node newPoint = newMedPoint(M);
				// Update tangency points
				double max_theta = -10.;
				double min_theta = 10.;
				double theta;
				Node left_node;
				Node right_node;
				for (auto n:tangentNodes->value(n1.id()))
				{
					theta = oriented_angle(edge2vec(e,n1),vec(n.point()+(-1.)*n1.point()));
					if (theta > 0. && theta < min_theta)
					{
						left_node = n;
						min_theta = theta;
					}
					if (theta < 0. && theta > max_theta)
					{
						right_node = n;
						max_theta = theta;
					}
				}
				std::vector<Node> t;
				std::vector<math::Point> p;
				t.push_back(left_node);
				t.push_back(right_node);
				p.push_back(left_node.point());
				p.push_back(right_node.point());
				tangentNodes->set(newPoint.id(),t);
				touchingPoints->set(newPoint.id(),p);
				if (iae->value(n1.id()) == 1 && iae->value(n2.id()) == 1)
					iae->set(newPoint.id(),1);
				Edge newEdge1 = newMedEdge(n1.id(),newPoint.id());
				Edge newEdge2 = newMedEdge(newPoint.id(),n2.id());
				NbEdgesAdded += 2;
				NbEdgesDeleted += 1;
				NbPointsAdded += 1;
				m_mesh_representation->deleteEdge(e);
			}
			continue;
		}
		if (e.length() > ATol)
		{
			int NbPointsToAdd = int (e.length()/mean_len);
			if (fabs(e.length()/mean_len - double (NbPointsToAdd))<0.5)
				NbPointsToAdd -= 1;
			Node previousPoint = n1;
			Node newPoint;
			for (int i = 1; i<=NbPointsToAdd; i++)
			{
				math::Point P = n1.point() + double (i)*mean_len*direction;
				newPoint = newMedPoint(P);
				NbPointsAdded += 1;
				// Updating the touching points, the medial radii and the dual triangles of the new medial points
				if (i < NbPointsToAdd / 2)
				{
					if (touchingPoints->value(n1.id()).size() == 2)
					{
						touchingPoints->set(newPoint.id(),touchingPoints->value(n1.id()));
						dualTriangles->set(newPoint.id(),dualTriangles->value(n1.id()));
						tangentNodes->set(newPoint.id(),tangentNodes->value(n1.id()));
					}
					else
					{
						touchingPoints->set(newPoint.id(),touchingPoints->value(n2.id()));
						dualTriangles->set(newPoint.id(),dualTriangles->value(n2.id()));
						tangentNodes->set(newPoint.id(),tangentNodes->value(n2.id()));
					}
					medRadius->set(newPoint.id(),newPoint.point().distance(touchingPoints->value(n1.id())[0]));
				}
				else
				{
					if (touchingPoints->value(n2.id()).size() == 2)
					{
						touchingPoints->set(newPoint.id(),touchingPoints->value(n2.id()));
						dualTriangles->set(newPoint.id(),dualTriangles->value(n2.id()));
						tangentNodes->set(newPoint.id(),tangentNodes->value(n2.id()));
					}
					else
					{
						touchingPoints->set(newPoint.id(),touchingPoints->value(n1.id()));
						dualTriangles->set(newPoint.id(),dualTriangles->value(n1.id()));
						tangentNodes->set(newPoint.id(),tangentNodes->value(n1.id()));
					}
					medRadius->set(newPoint.id(),newPoint.point().distance(touchingPoints->value(n2.id())[0]));
				}
				if (iae->value(n1.id()) == 1 && iae->value(n2.id()) == 1)
					iae->set(newPoint.id(),1);
				newMedEdge(previousPoint.id(),newPoint.id());
				NbEdgesAdded += 1;
				previousPoint = newPoint;
			}
			newMedEdge(previousPoint.id(),n2.id());
			NbEdgesAdded += 1;
			m_mesh_representation->deleteEdge(e);
			NbEdgesDeleted += 1;
		}
	}
	//updateConnectivity();
	std::cout<<NbPointsAdded<<" medial points added, "<<NbEdgesAdded<<" medial edges added and "<<NbEdgesDeleted<<" medial edges deleted."<<std::endl;
	std::cout<<"Done refining."<<std::endl;
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::setPrimalTriangleID(gmds::TCellID AMedPointID, gmds::TCellID ATriID)
{
	auto IDs = m_mesh_representation->getVariable<int,GMDS_NODE>("MP2triangle");
	IDs->set(AMedPointID,ATriID);
}

/*----------------------------------------------------------------------------*/
TCellID MedialAxis2D::primalTriangleID(gmds::TCellID AMedPointID)
{
	auto IDs = m_mesh_representation->getVariable<int,GMDS_NODE>("MP2triangle");
	return IDs->value(AMedPointID);
}

/*----------------------------------------------------------------------------*/
void
MedialAxis2D::setMedialPointType()
{
	// Fill up N2E connectivity
	//MeshDoctor doc(m_mesh_representation);
	//doc.updateUpwardConnectivity();
	auto var = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	for (auto n_id:m_mesh_representation->nodes())
	{
		Node n = m_mesh_representation->get<Node>(n_id);
		std::vector<Edge> adj_edges = n.get<Edge>();
		int size = adj_edges.size();
		var->set(n.id(),size);
	}
}

/*----------------------------------------------------------------------------*/
int MedialAxis2D::getMedialPointType(gmds::TCellID AId)
{
	auto var = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	return var->value(AId);
}

/*----------------------------------------------------------------------------*/
void
MedialAxis2D::setMedialEdgeType()
{
	// Requires medial point type already set
	auto medPointType = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	auto medEdgeType = m_mesh_representation->getVariable<int,GMDS_EDGE>("medial_edge_type");
	for (auto n_id:m_mesh_representation->edges())
	{
		Edge e = m_mesh_representation->get<Edge>(n_id);
		int type = 2;
		for (auto n:e.get<Node>())
		{
			if (medPointType->value(n.id()) == 1)
			{
				type = 1;
				break;
			}
			if (medPointType->value(n.id()) == 3)
				type = 3;
		}
		medEdgeType->set(e.id(),type);
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::setBoundaryConnectedComponentsIDs(gmds::TCellID AId, std::vector<int> &AV)
{
	auto IDs = m_mesh_representation->getVariable<std::vector<int>,GMDS_NODE>("boundary_connected_components_ids");
	IDs->set(AId,AV);
}

/*----------------------------------------------------------------------------*/
Eigen::MatrixXi MedialAxis2D::findOptimalBoundaryConnectedComponentsConnexion(int &AN)
{
	std::cout<<"> Finding optimal connexion through the medial axis between components"<<std::endl;

	auto cosMedAngle = m_mesh_representation->getVariable<double,GMDS_NODE>("cos_medial_angle");
	auto medRadius = m_mesh_representation->getVariable<double,GMDS_NODE>("medial_radius");
	auto componentsIDs = m_mesh_representation->getVariable<std::vector<int>,GMDS_NODE>("boundary_connected_components_ids");
	auto medialPointType = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");

	Eigen::MatrixXd optimalMedAngle, optimalMedRadius;
	Eigen::MatrixXi optimalMedPoint;
	optimalMedAngle.resize(AN,AN);
	optimalMedRadius.resize(AN,AN);
	optimalMedPoint.resize(AN,AN);
	for (int I = 1; I < AN; I ++)
	{
		for (int J = 0; J < I; J ++)
		{
			optimalMedPoint(I,J) = -1;
			optimalMedRadius(I,J) = 10e5;
			optimalMedAngle(I,J) = 2.;
		}
	}

	// Find the optimal medial angles between the components
	for (int I = 1; I < AN; I ++)
	{
		for (int J = 0; J < I; J ++)
		{
			for (auto n_id:m_mesh_representation->nodes())
			{
				if (medialPointType->value(n_id) == 2)
				{
					std::vector<int> ids = componentsIDs->value(n_id);
					if (ids.size() == 2)
					{
						if (((ids[0] == I) && (ids[1] == J)) || ((ids[1] == I) && (ids[0] == J)))
						{
							if (cosMedAngle->value(n_id) < optimalMedAngle(I,J))
								optimalMedAngle(I,J) = cosMedAngle->value(n_id);
						}

					}
				}
			}
		}
	}

	// Find the optimal medial radius and point between the components
	for (int I = 1; I < AN; I ++)
	{
		for (int J = 0; J < I; J ++)
		{
			for (auto n_id:m_mesh_representation->nodes())
			{
				if (medialPointType->value(n_id) == 2)
				{
					std::vector<int> ids = componentsIDs->value(n_id);
					if (ids.size() == 2)
					{
						if (((ids[0] == I) && (ids[1] == J)) || ((ids[1] == I) && (ids[0] == J)))
						{
							if (cosMedAngle->value(n_id) < optimalMedAngle(I,J) + 0.05)
							{
								if (medRadius->value(n_id) < optimalMedRadius(I,J))
								{
									optimalMedRadius(I,J) = medRadius->value(n_id);
									optimalMedPoint(I,J) = int(n_id);
								}
							}
						}

					}
				}
			}
		}
	}

	return optimalMedPoint;
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::setCosMedialAngle()
{
	auto cosMedAngle = m_mesh_representation->getVariable<double,GMDS_NODE>("cos_medial_angle");
	auto touchingPoints = m_mesh_representation->getVariable<std::vector<math::Point>,GMDS_NODE>("touching_points");
	auto medPointType = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	for (auto n_id:m_mesh_representation->nodes())
	{
		if (medPointType->value(n_id) != 2)
			cosMedAngle->set(n_id,-2.);
		else
		{
			Node MP = m_mesh_representation->get<Node>(n_id);
			math::Point P1 = touchingPoints->value(n_id)[0];
			math::Point P2 = touchingPoints->value(n_id)[1];
			P1 = P1 + (-1.)*MP.point();
			P2 = P2 + (-1.)*MP.point();
			double cosMA = vec(P1).dot(vec(P2))/(vec(P1).norm()*vec(P2).norm());
			cosMedAngle->set(n_id,cosMA);
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::smoothCosMedialAngle()
{
	// Test
	//auto isMarked = m_mesh_representation->newVariable<int,GMDS_NODE>("is_to_smooth");

	// Mark the nodes to smooth (we don't smooth the nodes with too high orthogonality default, IP, EP and local extremum)
	auto isToSmooth = m_mesh_representation->newMark<Node>();
	auto orthogonalityDefault = m_mesh_representation->getVariable<double,GMDS_NODE>("medial_radius_orthogonality_default");
	auto medPointType = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	auto cosMedAngle = m_mesh_representation->getVariable<double,GMDS_NODE>("cos_medial_angle");
	for (auto n_id:m_mesh_representation->nodes())
	{
		if (medPointType->value(n_id) == 2)
		{
			if (orthogonalityDefault->value(n_id) < 0.1)
			{
				// Check if it is a local extremum
				Node n = m_mesh_representation->get<Node>(n_id);
				std::vector<Node> closeNeighbours;
				Edge e = n.get<Edge>()[0];
				Node nxt = n;
				int N = 10; // Nb of neighbours on each side to check
				for (int i = 0; i < N; i++)
				{
					nxt = getNextPoint(nxt.id(),e.id());
					if (medPointType->value(nxt.id()) != 2)
						break;
					if (orthogonalityDefault->value(n_id) >= 0.1)
						break;
					closeNeighbours.push_back(nxt);
					e = getNextEdge(e.id(),nxt.id());
				}
				e = n.get<Edge>()[1];
				nxt = n;
				for (int i = 0; i < N; i++)
				{
					nxt = getNextPoint(nxt.id(),e.id());
					if (medPointType->value(nxt.id()) != 2)
						break;
					if (orthogonalityDefault->value(n_id) >= 0.1)
						break;
					closeNeighbours.push_back(nxt);
					e = getNextEdge(e.id(),nxt.id());
				}
				bool isExtremum = true;
				if (cosMedAngle->value(n_id) < cosMedAngle->value(closeNeighbours[0].id()))
				{
					for (auto n1:closeNeighbours)
					{
						if (cosMedAngle->value(n_id) > cosMedAngle->value(n1.id()))
						{
							isExtremum = false;
							break;
						}
					}
				}
				if (cosMedAngle->value(n_id) > cosMedAngle->value(closeNeighbours[0].id()))
				{
					for (auto n1:closeNeighbours)
					{
						if (cosMedAngle->value(n_id) < cosMedAngle->value(n1.id()))
						{
							isExtremum = false;
							break;
						}
					}
				}
				if (!isExtremum)
				{
					m_mesh_representation->mark<Node>(n_id,isToSmooth);
					//isMarked->set(n_id,1);
				}
			}
		}
	}

	// Smooth
	for (auto n_id:m_mesh_representation->nodes())
	{
		Node n = m_mesh_representation->get<Node>(n_id);
		if (m_mesh_representation->isMarked<Node>(n_id,isToSmooth))
		{
			std::vector<Node> neighbours;
			for (auto e:n.get<Edge>())
			{
				for (auto n2:e.get<Node>())
				{
					if ((n2.id() != n.id()) && (medPointType->value(n2.id()) == 2))
						neighbours.push_back(n2);
				}
			}
			if (m_mesh_representation->isMarked<Node>(neighbours[0].id(),isToSmooth))
			{
				if (m_mesh_representation->isMarked<Node>(neighbours[1].id(),isToSmooth))
				{
					double cosTheta1 = cosMedAngle->value(neighbours[0].id());
					double cosTheta2 = cosMedAngle->value(neighbours[1].id());
					cosMedAngle->set(n_id, (cosTheta1 + cosTheta2) / 2.);
				}
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::setTouchingPoints(const gmds::TCellID &APointID, std::vector<math::Point> APoints)
{
	auto touchingPoints = m_mesh_representation->getVariable<std::vector<math::Point>,GMDS_NODE>("touching_points");
	touchingPoints->set(APointID,APoints);
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::setTangentNodes(const gmds::TCellID &APointID, std::vector<Node> ANodes)
{
	auto tbn = m_mesh_representation->getVariable<std::vector<Node>,GMDS_NODE>("tangent_boundary_nodes");
	tbn->set(APointID,ANodes);
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::setDualTriangles(const gmds::TCellID &APointID, std::vector<Face> ATriangles)
{
	auto triangles = m_mesh_representation->getVariable<std::vector<Face>,GMDS_NODE>("dual_triangles");
	triangles->set(APointID,ATriangles);
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::setDangle(const TCellID &APointID, bool ADangleValue)
{
	auto iad = m_mesh_representation->getVariable<int,GMDS_NODE>("is_a_dangle");
	if (ADangleValue)
		iad->set(APointID,1);
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::markAsExtension(const TCellID &APointID)
{
	auto iae = m_mesh_representation->getVariable<int,GMDS_NODE>("is_an_extension");
	iae->set(APointID,1);
}

/*----------------------------------------------------------------------------*/
std::vector<math::Point> MedialAxis2D::getTouchingPoints(const gmds::TCellID &APointID)
{
	auto touchingPoints = m_mesh_representation->getVariable<std::vector<math::Point>,GMDS_NODE>("touching_points");
	return touchingPoints->value(APointID);
}

/*----------------------------------------------------------------------------*/
std::vector<Node> MedialAxis2D::getTangentNodes(const gmds::TCellID &APointID)
{
	auto tbn = m_mesh_representation->getVariable<std::vector<Node>,GMDS_NODE>("tangent_boundary_nodes");
	return tbn->value(APointID);
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::setFluxThroughMedialRadii()
{
	std::cout<<"> Computing the flux through medial radii"<<std::endl;
	auto cosMedAngle = m_mesh_representation->getVariable<double,GMDS_NODE>("cos_medial_angle");
	auto fluxThroughRadii = m_mesh_representation->getVariable<double,GMDS_NODE>("flux_through_medial_radii");
	auto medPointType = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	double epsilon = 0.09; // a small perturbation to avoid proliferation of singularities on branches where theta_M=pi/4 or 3pi/4
	epsilon = 0.;
	for (auto n_id:m_mesh_representation->nodes())
	{
		int type = medPointType->value(n_id);
		if (type == 2)
		{
			double flux;
			double medAngle;
			double cosMedialAngle = cosMedAngle->value(n_id);
			// Compute the medial angle
			if (fabs(cosMedialAngle+1) < 10e-6)
				medAngle = M_PI;
			else
			{
				if (fabs(cosMedialAngle-1) < 10e-6)
					medAngle = 0.;
				else
					medAngle = acos(cosMedAngle->value(n_id));
			}
			// Compute the flux
			if (medAngle < M_PI/4.+epsilon)
				flux = -medAngle;
			else
			{
				if (medAngle < 3.*M_PI/4.+epsilon)
					flux = M_PI/2.-medAngle;
				else
					flux = M_PI-medAngle;
			}
			fluxThroughRadii->set(n_id,flux);
		}
	}
}

/*----------------------------------------------------------------------------*/
double MedialAxis2D::orientation(const math::Vector &AU, const math::Vector &AV, const math::Vector &AE)
{
	if (AE.isZero())
		return 0.;
	math::Vector E = AE/AE.norm();
	math::Vector U = AU/AU.norm();
	math::Vector V = AV/AV.norm();
	if (E.dot(U) >= V.dot(U) && E.dot(V) >= U.dot(V) && (E.dot(U) >= 0. || E.dot(V) >= 0.))
		return 1.;
	else
		return -1.;
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::setFluxResidual()
{
	std::cout<<"> Computing the flux residual across medial edges and at intersection points"<<std::endl;
	auto fluxResidual = m_mesh_representation->getVariable<double,GMDS_EDGE>("flux_residual");
	auto fluxResAtIP=m_mesh_representation->getVariable<double,GMDS_NODE>("flux_residual_at_intersection_points");
	auto touchingPoints = m_mesh_representation->getVariable<std::vector<math::Point>,GMDS_NODE>("touching_points");
	auto flux = m_mesh_representation->getVariable<double,GMDS_NODE>("flux_through_medial_radii");
	auto medPointType = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	// Flux residual across edges
	for (auto e_id:m_mesh_representation->edges())
	{
		Edge medEdge = m_mesh_representation->get<Edge>(e_id);
		Node medPoint1 = medEdge.get<Node>()[0];
		Node medPoint2 = medEdge.get<Node>()[1];
		// We only deal with edges whose medial points both have type 2
		if (medPointType->value(medPoint1.id()) == 2 && medPointType->value(medPoint2.id()) == 2)
		{
			std::vector<math::Point> touchingPoints1 = touchingPoints->value(medPoint1.id());
			std::vector<math::Point> touchingPoints2 = touchingPoints->value(medPoint2.id());
			// Flux through medPoint1
			math::Vector E1 = vec(medPoint2.point() + (-1.)*medPoint1.point());
			math::Vector U1 = vec(touchingPoints1[0] + (-1.)*medPoint1.point());
			math::Vector V1 = vec(touchingPoints1[1] + (-1.)*medPoint1.point());
			double orientation1 = orientation(U1,V1,E1);
			double Phi1 = -orientation1*flux->value(medPoint1.id());
			// Flux through medPoint2
			double orientation2 = -orientation1;
			double Phi2 = -orientation2*flux->value(medPoint2.id());
			// Flux residual
			double deltaPhi = Phi1 + Phi2;
			fluxResidual->set(medEdge.id(),deltaPhi);
		}
	}
	// FLux residual at intersection points
	for (auto n_id:m_mesh_representation->nodes())
	{
		if (medPointType->value(n_id) > 2)
		{
			Node IP = m_mesh_representation->get<Node>(n_id);
			// Adjacent nodes
			std::vector<Node> adj_nodes;
			for (auto e:IP.get<Edge>())
			{
				for (auto n:e.get<Node>())
				{
					if (n.id() != IP.id())
						adj_nodes.push_back(n);
				}
			}
			// FLux residual
			double deltaPhi = 0.;
			for (auto n:adj_nodes)
			{
				std::vector<math::Point> touchingPoints1 = touchingPoints->value(n.id());
				if (touchingPoints1.size() == 2)
				{
					// Flux through n
					math::Vector E = vec(IP.point() + (-1.)*n.point());
					math::Vector U = vec(touchingPoints1[0] + (-1.)*n.point());
					math::Vector V = vec(touchingPoints1[1] + (-1.)*n.point());
					double orient = orientation(U,V,E);
					double Phi = -orient*flux->value(n.id());
					deltaPhi += Phi;
				}
			}
			fluxResAtIP->set(IP.id(),deltaPhi);
		}
	}
	// Test
//	for (auto e_id:m_mesh_representation->edges())
//	{
//		Edge e = m_mesh_representation->get<Edge>(e_id);
//		Node n1 = e.get<Node>()[0];
//		Node n2 = e.get<Node>()[1];
//		double flux1 = flux->value(n1.id());
//		double flux2 = flux->value(n2.id());
//		double residual = fluxResidual->value(e.id());
//		std::cout<<"Res : "<<residual<<" flux1 : "<<flux1<<" flux2 : "<<flux2<<std::endl;
//	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::setMedialRadius(const TCellID &APointID, double AValue)
{
	auto medRadius = m_mesh_representation->getVariable<double,GMDS_NODE>("medial_radius");
	Node n = m_mesh_representation->get<Node>(APointID);
	medRadius->set(n.id(),AValue);
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::setMedialRadiusOrthogonalityDefault(const gmds::TCellID &APointID, double AValue)
{
	auto orthoDef = m_mesh_representation->getVariable<double,GMDS_NODE>("medial_radius_orthogonality_default");
	orthoDef->set(APointID,AValue);
}

/*----------------------------------------------------------------------------*/
double MedialAxis2D::getMedialRadius(const gmds::TCellID &APointID)
{
	auto medRadius = m_mesh_representation->getVariable<double,GMDS_NODE>("medial_radius");
	return medRadius->value(APointID);
}

/*----------------------------------------------------------------------------*/
double MedialAxis2D::getMedialRadiusOrthogonalityDefault(const gmds::TCellID &APointID)
{
	auto orthoDef = m_mesh_representation->getVariable<double,GMDS_NODE>("medial_radius_orthogonality_default");
	return orthoDef->value(APointID);
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::placeSingularities(const double& AMedRadiusFraction)
{
	std::cout<<"> Placing singularities"<<std::endl;
	auto fluxResAtIP=m_mesh_representation->getVariable<double,GMDS_NODE>("flux_residual_at_intersection_points");
	auto fluxRes = m_mesh_representation->getVariable<double,GMDS_EDGE>("flux_residual");
	auto sing = m_mesh_representation-> getVariable<int,GMDS_NODE>("singularity");
	auto medRadius=m_mesh_representation->getVariable<double,GMDS_NODE>("medial_radius");
	int NbSingularities = 0;
	// Singularities at intersection points
	for (auto n_id:m_mesh_representation->nodes())
	{
		double res = fluxResAtIP->value(n_id);
		int singType = int (round(2.*res/M_PI));
		sing->set(n_id,singType);
		if (singType != 0)
			NbSingularities += 1;
	}
	std::cout<<"NB singularities at intersection points : "<<NbSingularities<<std::endl;
	// Singularities on medial branches
	for (auto e_id:m_mesh_representation->edges())
	{
		Edge e = m_mesh_representation->get<Edge>(e_id);
		double res = fluxRes->value(e_id);
		int singType = int (round(2.*res/M_PI));
		// If singType != 0, there isn't already a singularity here and we are not too close to the boundary, a singularity is placed
		if (singType != 0 && sing->value(e.get<Node>()[0].id()) == 0 && e.length() < medRadius->value(e.get<Node>()[0].id())/AMedRadiusFraction)
		{
			sing->set(e.get<Node>()[0].id(),singType);
			NbSingularities += 1;
		}
	}
	std::cout<<"Total NB singularities : "<<NbSingularities<<std::endl;
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::removeSingularityDipoles()
{
	std::cout<<"> Removing singularity dipoles"<<std::endl;
	auto sing = m_mesh_representation-> getVariable<int,GMDS_NODE>("singularity");
	// Remove dipoles
	for (auto e_id:m_mesh_representation->edges())
	{
		Edge e = m_mesh_representation->get<Edge>(e_id);
		Node n1 = e.get<Node>()[0];
		Node n2 = e.get<Node>()[1];
		if (sing->value(n1.id()) != 0 && sing->value(n1.id()) == -sing->value(n2.id()))
		{
			sing->set(n1.id(),0);
			sing->set(n2.id(),0);
		}
	}
	// Remove double singularities
	for (auto e_id:m_mesh_representation->edges())
	{
		Edge e = m_mesh_representation->get<Edge>(e_id);
		Node n1 = e.get<Node>()[0];
		Node n2 = e.get<Node>()[1];
		if (sing->value(n1.id()) != 0 && sing->value(n1.id()) == sing->value(n2.id()))
		{
			sing->set(n2.id(),0);
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::checkSingularities(const double &AOrthogonalityDefaultTol)
{
	std::cout<<"> Checking singularities"<<std::endl;
	auto sing = m_mesh_representation-> getVariable<int,GMDS_NODE>("singularity");
	auto orthoDef = m_mesh_representation->getVariable<double,GMDS_NODE>("medial_radius_orthogonality_default");
	int NbSing = 0;
	for (auto n_id:m_mesh_representation->nodes())
	{
		if (sing->value(n_id) != 0)
		{
			NbSing += 1;
			if (orthoDef->value(n_id) > AOrthogonalityDefaultTol)
			{
				sing->set(n_id,0);
				NbSing -= 1;
			}
		}
	}
	std::cout<<"NB singularities after check : "<<NbSing<<std::endl;
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::moveSingularitiesToIPs(const double &ANbPoints)
{
	std::cout<<"> Moving singularities to intersection points"<<std::endl;
	auto sing = m_mesh_representation-> getVariable<int,GMDS_NODE>("singularity");
	auto medPointType = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	int N = 0;
	for (auto n_id:m_mesh_representation->nodes())
	{
		if (sing->value(n_id) != 0)
		{
			if (medPointType->value(n_id) == 2)
			{
				Node n = m_mesh_representation->get<Node>(n_id);
				Edge e;
				Node nxt = n;
				bool found = false;
				// Look for an IP in the first direction
				e = n.get<Edge>()[0];
				for (int i = 1; i <= ANbPoints; i++)
				{
					nxt = getNextPoint(nxt.id(),e.id());
					if (medPointType->value(nxt.id()) == 1)
						break;
					if (medPointType->value(nxt.id()) > 2)
					{
						sing->set(nxt.id(),sing->value(n.id()));
						sing->set(n.id(),0);
						found = true;
						N += 1;
						break;
					}
					if (medPointType->value(nxt.id()) == 2)
						e = getNextEdge(e.id(),nxt.id());
				}
				// If not found in the first direction, try the second
				if (!found)
				{
					e = n.get<Edge>()[1];
					nxt = n;
					for (int i = 1; i <= ANbPoints; i++)
					{
						nxt = getNextPoint(nxt.id(),e.id());
						if (medPointType->value(nxt.id()) == 1)
							break;
						if (medPointType->value(nxt.id()) > 2)
						{
							sing->set(nxt.id(),sing->value(n.id()));
							sing->set(n.id(),0);
							N += 1;
							break;
						}
						if (medPointType->value(nxt.id()) == 2)
							e = getNextEdge(e.id(),nxt.id());
					}
				}
			}
		}
	}
	std::cout<<"NB singularities moved : "<<N<<std::endl;
}

/*----------------------------------------------------------------------------*/
Node
MedialAxis2D::getNextPoint(const TCellID &APointID, const TCellID &AEdgeID)
{
	// Requires E2N connections
	// To use only if the point APointID belongs to the edge AEdgeID
   Edge e = m_mesh_representation-> get<Edge>(AEdgeID);
	std::vector<Node> adj_nodes = e.get<Node>();
	Node n1 = adj_nodes[0];
	Node n2 = adj_nodes[1];
	if (n1.id() == APointID)
		return n2;
	else
		return n1;
}

/*----------------------------------------------------------------------------*/
Edge
MedialAxis2D::getNextEdge(const TCellID &AEdgeID, const TCellID &APointID)
{
	// Requires N2E connections
	// To use only if the point APointID belongs to the edge AEdgeID and has type 2
	Node n = m_mesh_representation-> get<Node>(APointID);
	std::vector<Edge> adj_edges = n.get<Edge>();
	Edge e1 = adj_edges[0];
	Edge e2 = adj_edges[1];
	if (e1.id() == AEdgeID)
		return e2;
	else
		return e1;
}

/*----------------------------------------------------------------------------*/
Node
MedialAxis2D::getExtremPoint(const TCellID &APointID, const TCellID &AEdgeID)
{
	// Requires medial points type set
	// Requires E2E connections
	// To use only if the point APointID belongs to the edge AEdgeID
	auto var = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	Node nxtPoint = getNextPoint(APointID, AEdgeID);
	int type = var->value(nxtPoint.id());
	Edge nxtEdge = m_mesh_representation->get<Edge>(AEdgeID);
	while (type == 2)
	{
		nxtEdge = getNextEdge(nxtEdge.id(), nxtPoint.id());
		nxtPoint = getNextPoint(nxtPoint.id(), nxtEdge.id());
		type = var->value(nxtPoint.id());
	}
	return nxtPoint;
}

/*----------------------------------------------------------------------------*/
void
MedialAxis2D::setBranchTypeOnPoints()
{
	// WARNING: this function is not optimized, its average complexity is in NbrPoints²/NbrBranches
	// Requires medial points type set
	auto branchType = m_mesh_representation->getVariable<int,GMDS_NODE>("point_branch_type");
	auto edgeBranchType = m_mesh_representation->getVariable<int,GMDS_EDGE>("edge_branch_type");
	for (auto n_id:m_mesh_representation->nodes())
	{
		Node n = m_mesh_representation->get<Node>(n_id);
		std::vector<Edge> adj_edges = n.get<Edge>();
		if (adj_edges.size() > 0)
			branchType->set(n_id,edgeBranchType->value(adj_edges[0].id()));
	}
}

/*----------------------------------------------------------------------------*/
void
MedialAxis2D::setBranchTypeOnEdges()
{
	// Requires branch IDs and medial edges type already set
	auto medEdgeType = m_mesh_representation->getVariable<int,GMDS_EDGE>("medial_edge_type");
	auto BranchType = m_mesh_representation->getVariable<int,GMDS_EDGE>("edge_branch_type");
	auto BranchID = m_mesh_representation->getVariable<int,GMDS_EDGE>("edge_branch_id");
	for (auto n_id:m_mesh_representation->edges())
	{
		BranchType->set(n_id,0);
	}
	for (auto n_id:m_mesh_representation->edges())
	{
		if (medEdgeType->value(n_id) == 1)
		{
			for (auto m_id:m_mesh_representation->edges())
			{
				if (BranchID->value(m_id) == BranchID->value(n_id))
					BranchType->set(m_id,1);
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
int
MedialAxis2D::setBranchIdOnEdges()
{
	// Requires medial points type set
	auto mpType = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	auto branchID = m_mesh_representation->getVariable<int,GMDS_EDGE>("edge_branch_id");
	// Mark if an edge has already been visited
	auto alreadyVisited = m_mesh_representation->newVariable<int,GMDS_EDGE>("already_visited");
	for (auto n_id:m_mesh_representation->edges())
	{
		alreadyVisited->set(n_id,0);
	}
	// Successive IDs of the branches
	int id = 0;
	for (auto n_id:m_mesh_representation->edges())
	{
		// Testing if this edge has already been visited
		if (alreadyVisited->value(n_id) == 0)
		{
			std::stack<TCellID> branch;
			branch.push(n_id);
			alreadyVisited->set(n_id,1);
			while (!branch.empty())
			{
				TCellID i = branch.top();
				branch.pop();
				branchID->set(i,id);
				Edge e = m_mesh_representation->get<Edge>(i);
				for (auto n:e.get<Node>())
				{
					if (mpType->value(n.id()) == 2)
					{
						for (auto edge:n.get<Edge>())
						{
							if (alreadyVisited->value(edge.id()) == 0)
							{
								branch.push(edge.id());
								alreadyVisited->set(edge.id(),1);
							}
						}
					}
				}
			}
			id += 1;
		}
	}
	std::cout<<"NB medial branches: "<<id<<std::endl;
	m_mesh_representation->deleteVariable(GMDS_EDGE,alreadyVisited);
	return id;
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::setBranchIdOnPoints()
{
	// Requires branch IDs on edges
	auto edgeBranchID = m_mesh_representation->getVariable<int,GMDS_EDGE>("edge_branch_id");
	auto pointBranchID = m_mesh_representation->getVariable<int,GMDS_NODE>("point_branch_id");
	for (auto n_id:m_mesh_representation->nodes())
	{
		Node n = m_mesh_representation->get<Node>(n_id);
		std::vector<Edge> adj_edges = n.get<Edge>();
		if (adj_edges.size() == 1 || adj_edges.size() == 2)
			pointBranchID->set(n.id(),edgeBranchID->value(adj_edges[0].id()));
		else
			pointBranchID->set(n.id(),-1);
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::identifyDetailsBranches(const double &ATol)
{
	// Requires medial branches IDs, medial points type, medial radii
	auto pointBranchID = m_mesh_representation->getVariable<int,GMDS_NODE>("point_branch_id");
	auto edgeBranchID = m_mesh_representation->getVariable<int,GMDS_EDGE>("edge_branch_id");
	auto medPointType = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	auto medRadius = m_mesh_representation->getVariable<double,GMDS_NODE>("medial_radius");
	auto correspondsToDetails = m_mesh_representation->getVariable<int,GMDS_EDGE>("corresponds_to_details");
	auto pointCorrespondsToDetails = m_mesh_representation->getVariable<int,GMDS_NODE>("point_corresponds_to_details");
	// Number of branches
	int NbBranches = setBranchIdOnEdges();
	// Vector stocking the medial points of type >= 3 surrounding the branch for each branch
	std::vector<std::vector<Node>> intersectionPoints(NbBranches);
	for (auto n_id:m_mesh_representation->nodes())
	{
		if (medPointType->value(n_id) >= 3)
		{
			Node n = m_mesh_representation->get<Node>(n_id);
			std::vector<Edge> adj_edges = n.get<Edge>();
			std::vector<int> medialBranches;
			for (auto edge:adj_edges)
			{
				int branchID = edgeBranchID->value(edge.id());
				bool alreadySeen = false;
				for (auto ID:medialBranches)
				{
					if (ID == branchID)
						alreadySeen = true;
				}
				if (!alreadySeen)
					medialBranches.push_back(branchID);
			}
			for (auto ID:medialBranches)
				intersectionPoints[ID].push_back(n);
		}
	}

	// Vector stocking the max ranges for each branch
	std::vector<double> maxRanges(NbBranches);
	// Initializing max ranges to 0
	for (int i = 0; i<NbBranches; i++)
		maxRanges[i] = 0.;
	// Computing the max ranges
	for (auto n_id:m_mesh_representation->nodes())
	{
		if (medPointType->value(n_id) == 2 || medPointType->value(n_id) == 1)
		{
			Node n = m_mesh_representation->get<Node>(n_id);
			int branchID = pointBranchID->value(n_id);
			std::vector<Node> intPoints = intersectionPoints[branchID];
			if (intPoints.size() > 0)
			{
				double range = maxRange(intPoints[0].point(),medRadius->value(intPoints[0].id()),n.point(),medRadius->value(n.id()));
				for (auto intPoint:intPoints)
				{
					double range1 = maxRange(intPoint.point(),medRadius->value(intPoint.id()),n.point(),medRadius->value(n.id()));
					if (range1 < range)
						range = range1;
				}
				if (range > maxRanges[branchID])
					maxRanges[branchID] = range;
			}
			else
				// This means that this branch has no intersection point, so it doesn't correspond to a detail
				maxRanges[branchID] = ATol + 1.;

		}
	}

	// Marking the edges whose branch corresponds to details
	for (auto e_id:m_mesh_representation->edges())
	{
		if (maxRanges[edgeBranchID->value(e_id)] < ATol)
			correspondsToDetails->set(e_id,1);
		else
			correspondsToDetails->set(e_id,0);
	}
	// Marking the points whose branch corresponds to details
	for (auto n_id:m_mesh_representation->nodes())
	{
		Node n = m_mesh_representation->get<Node>(n_id);
		if (n.get<Edge>().size() != 1 && n.get<Edge>().size() != 2)
			pointCorrespondsToDetails->set(n.id(),0);
		else
		{
			pointCorrespondsToDetails->set(n.id(),correspondsToDetails->value(n.get<Edge>()[0].id()));
		}
	}
}

/*----------------------------------------------------------------------------*/
int MedialAxis2D::correspondsToDetail(const gmds::TCellID &APointID)
{
	// Requires points_correspond_to_detail
	auto pointCorrespondsToDetails = m_mesh_representation->getVariable<int,GMDS_NODE>("point_corresponds_to_details");
	return (pointCorrespondsToDetails->value(APointID));
}

/*----------------------------------------------------------------------------*/
void
MedialAxis2D::write(std::basic_string<char> AFileName)
{
	IGMeshIOService ioService(m_mesh_representation);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(AFileName);
}
/*----------------------------------------------------------------------------*/
}  // end namespace medialaxis
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/