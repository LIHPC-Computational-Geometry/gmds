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
	// Medial point/group correspondence
	m_mesh_representation->newVariable<int,GMDS_NODE>("medialPoint2Group");
	// List of boundary connected components connected by this medial point
	m_mesh_representation->newVariable<std::vector<int>,GMDS_NODE>("boundary_connected_components_ids");

	// Topological representation
	m_topological_representation = new Mesh(MeshModel(DIM3 | F | E | N |
	                                           F2N | E2N | N2E | N2F));
	// Type of the medial axis section: 0 if it is a quad edge // to the boundary, 1 if it is a quad diagonal, 2 if it is an edge orthogonal to the boundary
	m_topological_representation->newVariable<int,GMDS_EDGE>("section_type");
	// Type of the node: 0 if it is internal (there is an equation at this point), 1 if not (there is a border condition at this point)
	m_topological_representation->newVariable<int,GMDS_NODE>("node_type");
	// ID of the section to which the point belongs
	m_mesh_representation->newVariable<int,GMDS_NODE>("section_id");
	// ID of the subsection to which the point belongs
	m_mesh_representation->newVariable<int,GMDS_NODE>("subsection_id");
	// Type of the subsection to which the point belongs
	m_mesh_representation->newVariable<int,GMDS_NODE>("medial_section_type");
	// Section/section ID correspondence
	m_topological_representation->newVariable<int,GMDS_EDGE>("section_to_section_id");
	// Correspondence medial point/node of the topological representation
	m_mesh_representation->newVariable<int,GMDS_NODE>("med_point_to_sing_node");
	// Correspondence medial point/node of the topological representation
	m_topological_representation->newVariable<int,GMDS_NODE>("sing_node_to_med_point");
	// Position of the wings on a section (1 if the wings point in the same direction as the oriented section, -1 if in the opposite direction, 0 if no wing)
	m_topological_representation->newVariable<int,GMDS_EDGE>("wings_position");
	// Position of each section in the vector of degrees of freedom
	m_topological_representation->newVariable<int,GMDS_EDGE>("section2matrix");
	// Quantization degrees of freedom graph representation
	// Dof graph
	m_dof_graph= new Mesh(MeshModel(DIM3 | E | N |
	                                                    E2N | N2E));
	// Correspondence degree of freedom/vertex in the dof graph
	m_topological_representation->newVariable<int,GMDS_EDGE>("section2LeftQuadLength");
	m_topological_representation->newVariable<int,GMDS_EDGE>("section2RightQuadLength");
	m_topological_representation->newVariable<int,GMDS_EDGE>("section2QuadHeight");
	m_topological_representation->newVariable<int,GMDS_EDGE>("section2LeftDiagoQuadLength");
	m_topological_representation->newVariable<int,GMDS_EDGE>("section2RightDiagoQuadLength");
	// Solution to the quantization problem
	m_dof_graph->newVariable<double,GMDS_NODE>("quantization_solution");
	// Topology of a medial axis based block decomposition
	m_ma_block_decomposition = new Mesh(MeshModel(DIM3 | F | E | N |
	                                                  F2N | E2N | N2E | N2F));
	// Going from topo rep singu node to its corresponding block decomp medial node
	m_topological_representation->newVariable<int,GMDS_NODE>("singuNode2MedNode");
	// Going from topo rep singu node to its corresponding block decomp boundary nodes
	m_topological_representation->newVariable<std::vector<int>,GMDS_NODE>("singuNode2BoundNodes");
	// Going from topo rep singu node to its corresponding block decomp middle nodes
	m_topological_representation->newVariable<std::vector<int>,GMDS_NODE>("singuNode2MiddleNode");
	// Directions of the medial section
	m_topological_representation->newVariable<math::Vector,GMDS_EDGE>("tailDirection");
	m_topological_representation->newVariable<math::Vector,GMDS_EDGE>("headDirection");
	// Correspondence section/block decomposition nodes
	m_topological_representation->newVariable<int,GMDS_EDGE>("tailNode");
	m_topological_representation->newVariable<int,GMDS_EDGE>("headNode");
	m_topological_representation->newVariable<int,GMDS_EDGE>("leftTailBoundaryNode");
	m_topological_representation->newVariable<int,GMDS_EDGE>("rightTailBoundaryNode");
	m_topological_representation->newVariable<int,GMDS_EDGE>("leftHeadBoundaryNode");
	m_topological_representation->newVariable<int,GMDS_EDGE>("rightHeadBoundaryNode");
	m_topological_representation->newVariable<int,GMDS_EDGE>("leftHeadMiddleNode");
	m_topological_representation->newVariable<int,GMDS_EDGE>("rightHeadMiddleNode");
	m_topological_representation->newVariable<int,GMDS_EDGE>("leftTailMiddleNode");
	m_topological_representation->newVariable<int,GMDS_EDGE>("rightTailMiddleNode");
}
MedialAxis2D::~MedialAxis2D()
{
	if(m_mesh_representation!= nullptr)
		delete m_mesh_representation;
	if(m_topological_representation!= nullptr)
		delete m_topological_representation;
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
					if (e.length() < ATol || (medPointType->value(n1.id()) == 3 && medPointType->value(n2.id()) == 3))
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
	auto touchingPoints = m_mesh_representation->getVariable<std::vector<math::Point>,GMDS_NODE>("touching_points");
	auto medRadius = m_mesh_representation->getVariable<double,GMDS_NODE>("medial_radius");
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
		if ((e.length() > ATol) && ((touchingPoints->value(n1.id()).size() != 2) && (touchingPoints->value(n2.id()).size() != 2)))
		{
			std::cout<<"Could not refine an edge because it is incident to no regular medial point (medial point type 2)."<<std::endl;
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
				// Updating the touching points and the medial radii of the new medial points
				if (i < NbPointsToAdd / 2)
				{
					if (touchingPoints->value(n1.id()).size() == 2)
						touchingPoints->set(newPoint.id(),touchingPoints->value(n1.id()));
					else
						touchingPoints->set(newPoint.id(),touchingPoints->value(n2.id()));
					medRadius->set(newPoint.id(),newPoint.point().distance(touchingPoints->value(n1.id())[0]));
				}
				else
				{
					if (touchingPoints->value(n2.id()).size() == 2)
						touchingPoints->set(newPoint.id(),touchingPoints->value(n2.id()));
					else
						touchingPoints->set(newPoint.id(),touchingPoints->value(n1.id()));
					medRadius->set(newPoint.id(),newPoint.point().distance(touchingPoints->value(n2.id())[0]));
				}
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
std::vector<math::Point> MedialAxis2D::getTouchingPoints(const gmds::TCellID &APointID)
{
	auto touchingPoints = m_mesh_representation->getVariable<std::vector<math::Point>,GMDS_NODE>("touching_points");
	return touchingPoints->value(APointID);
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::setFluxThroughMedialRadii()
{
	std::cout<<"> Computing the flux through medial radii"<<std::endl;
	auto cosMedAngle = m_mesh_representation->getVariable<double,GMDS_NODE>("cos_medial_angle");
	auto fluxThroughRadii = m_mesh_representation->getVariable<double,GMDS_NODE>("flux_through_medial_radii");
	auto medPointType = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
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
			if (medAngle < M_PI/4.)
				flux = -medAngle;
			else
			{
				if (medAngle < 3.*M_PI/4.)
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
	auto medPointType = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	for (auto n_id:m_mesh_representation->nodes())
	{
		Node n = m_mesh_representation->get<Node>(n_id);
		std::vector<Edge> adj_edges = n.get<Edge>();
		int size = adj_edges.size();
		if (size == 3)
			branchType->set(n.id(),2);
		if (size == 2)
		{
			Edge e1 = adj_edges[0];
			Node extremPoint1 = getExtremPoint(n.id(),e1.id());
			Edge e2 = adj_edges[1];
			Node extremPoint2 = getExtremPoint(n.id(),e2.id());
			int mpt1 = medPointType->value(extremPoint1.id());
			int mpt2 = medPointType->value(extremPoint2.id());
			if ((mpt1 == 1)||(mpt2 == 1))
				branchType->set(n.id(),1);
			else
				branchType->set(n.id(),0);
		}
		if (size == 1)
		{
			branchType->set(n.id(),1);
		}
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
void MedialAxis2D::setSectionID()
{
	std::cout<<"> Setting sections IDs"<<std::endl;
	auto medPoint2singNode = m_mesh_representation->getVariable<int,GMDS_NODE>("med_point_to_sing_node");
	auto sectionID = m_mesh_representation->getVariable<int,GMDS_NODE>("section_id");
	auto subSectionID = m_mesh_representation->getVariable<int,GMDS_NODE>("subsection_id");
	auto alreadyVisited = m_mesh_representation->newVariable<int,GMDS_NODE>("already_visited");
	auto medPointType = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	auto singNode2medPoint = m_topological_representation->getVariable<int,GMDS_NODE>("sing_node_to_med_point");
	auto nodeType = m_topological_representation->getVariable<int,GMDS_NODE>("node_type");
	auto singu = m_mesh_representation->getVariable<int,GMDS_NODE>("singularity");
	int ID = 0;
	std::vector<Node> section_start; // Vector storing a starting point for each section
	std::vector<Edge> section_direction; // Vector storing a starting edge containing the starting point for each section
	for (auto n_id:m_mesh_representation->nodes())
	{
		if ((medPointType->value(n_id) == 2) && (singu->value(n_id) == 0) && (alreadyVisited->value(n_id) == 0))
		{
			// We build a new section (= new edge of the topo rep)

			// Find the extreme nodes of the section
			Node n = m_mesh_representation->get<Node>(n_id);
			alreadyVisited->set(n_id,1);
			sectionID->set(n.id(),ID);
			Node prev, nxt;
			Edge e;
			// We go to the left
			prev = n;
			e = n.get<Edge>()[0];
			nxt = getNextPoint(prev.id(),e.id());
			while ((medPointType->value(nxt.id()) == 2) && (singu->value(nxt.id()) == 0))
			{
				alreadyVisited->set(nxt.id(),1);
				sectionID->set(nxt.id(),ID);
				e = getNextEdge(e.id(),nxt.id());
				prev = nxt;
				nxt = getNextPoint(prev.id(),e.id());
			}
			section_start.push_back(nxt);
			section_direction.push_back(e);
			//Node n1 = nxt;
			// We go to the right
			prev = n;
			e = n.get<Edge>()[1];
			nxt = getNextPoint(prev.id(),e.id());
			while ((medPointType->value(nxt.id()) == 2) && (singu->value(nxt.id()) == 0))
			{
				alreadyVisited->set(nxt.id(),1);
				sectionID->set(nxt.id(),ID);
				e = getNextEdge(e.id(),nxt.id());
				prev = nxt;
				nxt = getNextPoint(prev.id(),e.id());
			}
			//Node n2 = nxt;
			ID += 1;
		}
	}

	// // Now we build the subsections, in order to refine the mesh
	// int subSecID = 0;
	// for (int id = 0; id < ID; id++)
	// {
	// 	// Chose one medial point every 10 points on the section to refine it
	// 	Node prev = section_start[id];
	// 	Edge e = section_direction[id];
	// 	Node nxt = getNextPoint(prev.id(),e.id());
	// 	int steps = 1;
	// 	bool Continue = true;
	// 	while (sectionID->value(nxt.id()) == id && Continue)
	// 	{
	// 		e = getNextEdge(e.id(),nxt.id());
	// 		prev = nxt;
	// 		if (medPointType->value(prev.id()) == 1)
	// 			Continue = false;
	// 		nxt = getNextPoint(prev.id(),e.id());
	// 		if (steps < 15)
	// 		{
	// 			subSectionID->set(prev.id(),subSecID);
	// 			steps += 1;
	// 		}
	// 		else
	// 		{
	// 			if (sectionID->value(nxt.id()) == id)
	// 			{
	// 				Node newNode = m_topological_representation->newNode(prev.point());
	// 				medPoint2singNode->set(prev.id(),newNode.id());
	// 				singNode2medPoint->set(newNode.id(),prev.id());
	// 				nodeType->set(newNode.id(),2);
	// 				subSecID += 1;
	// 				steps = 0;
	// 			}
	// 		}


	// 	}
	// }
	

	std::cout<<"NB sections : "<<ID<<std::endl;
	m_mesh_representation->deleteVariable(GMDS_NODE,alreadyVisited);
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::computeSectionType()
{
	auto sectionID = m_mesh_representation->getVariable<int,GMDS_NODE>("section_id");
	auto cosMedAngle = m_mesh_representation->getVariable<double,GMDS_NODE>("cos_medial_angle");
	auto type = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_section_type");
	// Compute the highest section ID
	int maxID = 0;
	for (auto n_id:m_mesh_representation->nodes())
	{
		if (sectionID->value(n_id) > maxID)
			maxID = sectionID->value(n_id);
	}
	// Vectors storing the number of medial point of each section for each type
	std::vector<int> type0(maxID+1);
	std::vector<int> type1(maxID+1);
	for (int i = 0; i <= maxID; i++)
	{
		type0[i] = 0;
		type1[i] = 0;
	}
	for (auto n_id:m_mesh_representation->nodes())
	{
		int id = sectionID->value(n_id);
		if (cosMedAngle->value(n_id) < -sqrt(2.)/2.)
			type0[id] = type0[id] + 1;
		else
		{
			if (fabs(cosMedAngle->value(n_id)) < sqrt(2.)/2.)
				type1[id] = type1[id] + 1;
			else
				type0[id] = type0[id] + 1;
		}
	}
	for (auto n_id:m_mesh_representation->nodes())
	{
		int id = sectionID->value(n_id);
		if (type0[id] > type1[id])
			type->set(n_id,0);
		else
			type->set(n_id,1);
	}
}

/*----------------------------------------------------------------------------*/
std::vector<Node> MedialAxis2D::sectionNodes(int AID)
{
	auto sectionID = m_mesh_representation->getVariable<int,GMDS_NODE>("section_id");
	auto medPointType = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	auto singu = m_mesh_representation->getVariable<int,GMDS_NODE>("singularity");
	// Find a node of the section
	Node n0;
	for (auto n_id:m_mesh_representation->nodes())
	{
		if ((medPointType->value(n_id) == 2) && (singu->value(n_id) == 0) && (sectionID->value(n_id) == AID))
		{
			n0 = m_mesh_representation->get<Node>(n_id);
			break;
		}
	}
	// Build the two halfs of the section
	std::vector<Node> sec;
	Node prev, nxt;
	Edge e;
	// We go to the left
	prev = n0;
	e = n0.get<Edge>()[0];
	nxt = getNextPoint(prev.id(),e.id());
	while ((medPointType->value(nxt.id()) == 2) && (singu->value(nxt.id()) == 0))
	{
		sec.push_back(prev);
		e = getNextEdge(e.id(),nxt.id());
		prev = nxt;
		nxt = getNextPoint(prev.id(),e.id());
	}
	if (!sec.empty())
		sec.erase(sec.begin());
	std::reverse(sec.begin(),sec.end());
	// We go to the right
	prev = n0;
	e = n0.get<Edge>()[1];
	nxt = getNextPoint(prev.id(),e.id());
	while ((medPointType->value(nxt.id()) == 2) && (singu->value(nxt.id()) == 0))
	{
		sec.push_back(prev);
		e = getNextEdge(e.id(),nxt.id());
		prev = nxt;
		nxt = getNextPoint(prev.id(),e.id());
	}
	return sec;
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::refineByAddingSingularNodes()
{
	std::cout<<"> Refining the future block decomposition by adding singular nodes"<<std::endl;
	auto sectionID = m_mesh_representation->getVariable<int,GMDS_NODE>("section_id");
	auto singu = m_mesh_representation->getVariable<int,GMDS_NODE>("singularity");
	// Compute the highest section ID
	int maxID = 0;
	for (auto n_id:m_mesh_representation->nodes())
	{
		if (sectionID->value(n_id) > maxID)
			maxID = sectionID->value(n_id);
	}
	std::vector<Node> section;
	for (int id = 0; id <= maxID; id++)
	{
		section = sectionNodes(id);
		int Nb = section.size();
		int q = Nb/10;
		for (int i = 1; i < q; i++)
			singu->set(section[10*i].id(),10);
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::buildTopoRepNodes()
{
	std::cout<<"> Building nodes of the topological representation"<<std::endl;
	std::cout<<"> WARNING: if the medial axis has two neighbouring singular points (for example a singu and an IP), the topo rep is not valid"<<std::endl;
	auto singu = m_mesh_representation->getVariable<int,GMDS_NODE>("singularity");
	auto medPointType = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	auto medPoint2singNode = m_mesh_representation->getVariable<int,GMDS_NODE>("med_point_to_sing_node");
	auto singNode2medPoint = m_topological_representation->getVariable<int,GMDS_NODE>("sing_node_to_med_point");
	auto nodeType = m_topological_representation->getVariable<int,GMDS_NODE>("node_type");

	for (auto n_id:m_mesh_representation->nodes())
	{
		if (medPointType->value(n_id) == 2)
		{
			if (singu->value(n_id) != 0)
			{
				Node n = m_mesh_representation->get<Node>(n_id);
				Node newNode = m_topological_representation->newNode(n.point());
				medPoint2singNode->set(n_id,newNode.id());
				singNode2medPoint->set(newNode.id(),n.id());
				nodeType->set(newNode.id(),2);
			}
			else
				medPoint2singNode->set(n_id,-1);
		}
		else
		{
			Node n = m_mesh_representation->get<Node>(n_id);
			Node newNode = m_topological_representation->newNode(n.point());
			medPoint2singNode->set(n.id(),newNode.id());
			singNode2medPoint->set(newNode.id(),n.id());
			nodeType->set(newNode.id(),medPointType->value(n_id));
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::buildTopoRepEdges()
{
	std::cout<<"> Building edges of the topological representation"<<std::endl;
	std::cout<<"> WARNING: if the medial axis has only one section which is a cycle, the topo rep is not valid"<<std::endl;
	auto medPointType = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	auto cosMedAngle = m_mesh_representation->getVariable<double,GMDS_NODE>("cos_medial_angle");
	auto sectionType = m_topological_representation->getVariable<int,GMDS_EDGE>("section_type");
	auto secTypeOnNodes = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_section_type");
	auto medPoint2singNode = m_mesh_representation->getVariable<int,GMDS_NODE>("med_point_to_sing_node");
	auto wings = m_topological_representation->getVariable<int,GMDS_EDGE>("wings_position");
	auto touchingPoints = m_mesh_representation->getVariable<std::vector<math::Point>,GMDS_NODE>("touching_points");
	auto section2sectionID = m_topological_representation->getVariable<int,GMDS_EDGE>("section_to_section_id");
	auto sectionID = m_mesh_representation->getVariable<int,GMDS_NODE>("section_id");
	//auto subSectionID = m_mesh_representation->getVariable<int,GMDS_NODE>("subsection_id");
	auto tailDir = m_topological_representation->getVariable<math::Vector,GMDS_EDGE>("tailDirection");
	auto headDir = m_topological_representation->getVariable<math::Vector,GMDS_EDGE>("headDirection");
	for (auto n_id:m_mesh_representation->nodes())
	{
		if (medPoint2singNode->value(n_id) == -1 && medPointType->value(n_id) == 2)
		{
			// We build a new section (= new edge of the topo rep)
			// Type of the section
			int type = secTypeOnNodes->value(n_id);
			// if (cosMedAngle->value(n_id) < -sqrt(2.)/2.)
			// 	type = 0;
			// else
			// {
			// 	if (fabs(cosMedAngle->value(n_id)) < sqrt(2.)/2.)
			// 		type = 1;
			// 	else
			// 		type = 0;
			// }

			//int sectID = subSectionID->value(n_id);
			int sectID = sectionID->value(n_id);

			// Find the extreme nodes of the section
			Node n = m_mesh_representation->get<Node>(n_id);
			medPoint2singNode->set(n_id,-2);
			Node prev, nxt;
			Edge e;
			// We go to the left
			prev = n;
			e = n.get<Edge>()[0];
			nxt = getNextPoint(prev.id(),e.id());
			while (medPoint2singNode->value(nxt.id()) < 0)
			{
				medPoint2singNode->set(nxt.id(),-2);
				e = getNextEdge(e.id(),nxt.id());
				prev = nxt;
				nxt = getNextPoint(prev.id(),e.id());
			}
			// Tail direction of the section
			math::Vector tail_dir = edge2vec(e,nxt);
			Node n1 = nxt;
			// We go to the right
			prev = n;
			e = n.get<Edge>()[1];
			nxt = getNextPoint(prev.id(),e.id());
			while (medPoint2singNode->value(nxt.id()) < 0)
			{
				medPoint2singNode->set(nxt.id(),-2);
				e = getNextEdge(e.id(),nxt.id());
				prev = nxt;
				nxt = getNextPoint(prev.id(),e.id());
			}
			// Head direction of the section
			math::Vector head_dir = edge2vec(e,prev);
			Node n2 = nxt;
			TCellID N1 = medPoint2singNode->value(n1.id());
			TCellID N2 = medPoint2singNode->value(n2.id());
			Edge newSection = m_topological_representation->newEdge(N1,N2);
			// Directions of the section
			tailDir->set(newSection.id(),tail_dir);
			headDir->set(newSection.id(),head_dir);
			// All end sections must be of type 1
			if (medPointType->value(n1.id()) == 1 || medPointType->value(n2.id()) == 1)
				type = 1;
			sectionType->set(newSection.id(),type);
			section2sectionID->set(newSection.id(),sectID);
			// Wings position
			Node P1 = newSection.get<Node>()[0];
			math::Point S = getNextPoint(n_id,n.get<Edge>()[1].id()).point() + (-1.)*n.point();
			math::Point R = touchingPoints->value(n_id)[0] + (-1.)*n.point();
			double orient = vec(S).dot(vec(R));
			if (type == 1)
			{
				if (orient >= 0)
					wings->set(newSection.id(),-1);
				else
					wings->set(newSection.id(),1);
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
void
MedialAxis2D::writeTopoRep(std::basic_string<char> AFileName)
{
	std::cout<<"> Writing the topological representation"<<std::endl;
	IGMeshIOService ioService(m_topological_representation);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(AFileName);
}

/*----------------------------------------------------------------------------*/
int MedialAxis2D::wings(gmds::Node &AN, gmds::Edge &AE)
{
	auto wings = m_topological_representation->getVariable<int,GMDS_EDGE>("wings_position");
	if (orientation(AN,AE) == wings->value(AE.id()))
		return fabs(orientation(AN,AE));
	else
		return 0;
}

/*----------------------------------------------------------------------------*/
int MedialAxis2D::orientation(gmds::Node &AN, gmds::Edge &AE)
{
	if (AN.id() == AE.get<Node>()[0].id())
		return 1;
	if (AN.id() == AE.get<Node>()[1].id())
		return -1;
	return 0;
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::setTopoRepConnectivity()
{
	std::cout<<"> Setting topological representation connectivity"<<std::endl;
	MeshDoctor doc(m_topological_representation);
	doc.updateUpwardConnectivity();
}

/*----------------------------------------------------------------------------*/
Node MedialAxis2D::getNextSingularNode(gmds::Node &AN, gmds::Edge &AE)
{
	if (AE.get<Node>()[0].id() == AN.id())
		return AE.get<Node>()[1];
	if (AE.get<Node>()[1].id() == AN.id())
		return AE.get<Node>()[0];
	std::cout<<"getNextSingularNode : error, the given Node does not belong to the given edge"<<std::endl;
	return AN;
}

/*----------------------------------------------------------------------------*/
Edge MedialAxis2D::getNextMedialSection(gmds::Edge &AE, gmds::Node &AN)
{
	if (AN.get<Edge>().size() != 2)
	{
		std::cout<<"getNextMedialSection : error, the given node doesn't have valence 2"<<std::endl;
		return AE;
	}
	if (AN.get<Edge>()[0].id() == AE.id())
		return AN.get<Edge>()[1];
	if (AN.get<Edge>()[1].id() == AE.id())
		return AN.get<Edge>()[0];
	std::cout<<"getNextMedialSection : error, the given Node does not belong to the given edge"<<std::endl;
	return AE;
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::browseTopoRep()
{
	Node current_node = m_topological_representation->get<Node>(0);
	Node prev_node;
	Edge current_edge = current_node.get<Edge>()[0];
	TCellID  initial_edge_id = current_edge.id();
	while (true)
	{
		std::cout<<"We are currently at node "<<current_node.id()<<" and at edge "<<current_edge.id()<<std::endl;
		prev_node = current_node;
		current_node = getNextSingularNode(current_node,current_edge);
		if (current_node.get<Edge>().size() == 2)
			current_edge = getNextMedialSection(current_edge,current_node);
		if (current_node.get<Edge>().size() >= 3)
		{
			Edge nxt_edge;
			double theta = 10.;
			for (auto e:current_node.get<Edge>())
			{
				if (e.id() != current_edge.id())
				{
					if (oriented_angle(edge2vec(current_edge,prev_node), edge2vec(e,current_node)) < theta)
					{
						nxt_edge = e;
						theta = oriented_angle(edge2vec(current_edge,prev_node), edge2vec(e,current_node));
					}
				}
			}
			current_edge = nxt_edge;
		}
		if (current_node.id() == 0 && current_edge.id() == initial_edge_id)
			break;
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::buildBlockDecompMedialAndBoundaryNodes()
{
	std::cout<<"> Building medial and boundary nodes of the medial axis based block decomposition"<<std::endl;
	auto nodeType = m_topological_representation->getVariable<int,GMDS_NODE>("node_type");
	auto singNode2medPoint = m_topological_representation->getVariable<int,GMDS_NODE>("sing_node_to_med_point");
	auto touchingPoints = m_mesh_representation->getVariable<std::vector<math::Point>,GMDS_NODE>("touching_points");
	auto sn2mn = m_topological_representation->getVariable<int,GMDS_NODE>("singuNode2MedNode");
	auto sn2bn = m_topological_representation->getVariable<std::vector<int>,GMDS_NODE>("singuNode2BoundNodes");
	for (auto n_id:m_topological_representation->nodes())
	{
		Node singu_node = m_topological_representation->get<Node>(n_id);
		Node newNode;
		newNode = m_ma_block_decomposition->newNode(singu_node.point());
		sn2mn->set(n_id,newNode.id());
		if (nodeType->value(n_id) == 1)
			continue;
		Node med_point = m_mesh_representation->get<Node>(singNode2medPoint->value(n_id));
		std::vector<int> newNodes;
		std::vector<math::Point> boundaryTangencyPoints = touchingPoints->value(med_point.id());

		if (nodeType->value(n_id) >= 2)
		{
			for (auto p:boundaryTangencyPoints)
			{
				newNode = m_ma_block_decomposition->newNode(p);
				newNodes.push_back(newNode.id());
			}
		}
		sn2bn->set(n_id,newNodes);
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::writeBlockDecomp(std::basic_string<char> AFileName)
{
	std::cout<<"> Writing the medial axis based block decomposition"<<std::endl;
	IGMeshIOService ioService(m_ma_block_decomposition);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(AFileName);
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::buildSection2MedialAndBoundaryNodesAdjacency()
{
	auto TN = m_topological_representation->getVariable<int,GMDS_EDGE>("tailNode");
	auto HN = m_topological_representation->getVariable<int,GMDS_EDGE>("headNode");
	auto LTBN = m_topological_representation->getVariable<int,GMDS_EDGE>("leftTailBoundaryNode");
	auto RTBN = m_topological_representation->getVariable<int,GMDS_EDGE>("rightTailBoundaryNode");
	auto LHBN = m_topological_representation->getVariable<int,GMDS_EDGE>("leftHeadBoundaryNode");
	auto RHBN = m_topological_representation->getVariable<int,GMDS_EDGE>("rightHeadBoundaryNode");
	auto sn2mn = m_topological_representation->getVariable<int,GMDS_NODE>("singuNode2MedNode");
	auto sn2bn = m_topological_representation->getVariable<std::vector<int>,GMDS_NODE>("singuNode2BoundNodes");
	auto nodeType = m_topological_representation->getVariable<int,GMDS_NODE>("node_type");
	auto tailDir = m_topological_representation->getVariable<math::Vector,GMDS_EDGE>("tailDirection");
	auto headDir = m_topological_representation->getVariable<math::Vector,GMDS_EDGE>("headDirection");

	for (auto n_id:m_topological_representation->nodes())
	{
		Node n = m_topological_representation->get<Node>(n_id);
		int medial_node = sn2mn->value(n_id);
		std::vector<int> boundary_nodes = sn2bn->value(n_id);
		for (auto section:n.get<Edge>())
		{
			if (orientation(n,section) == 1)
			{
				TN->set(section.id(),medial_node);
				if (nodeType->value(n_id) == 1)
				{
					LTBN->set(section.id(),-1);
					RTBN->set(section.id(),-1);
				}
				if (nodeType->value(n_id) > 1)
				{
					// Find the left and right boundary points
					double max_theta = -10.;
					double min_theta = 10.;
					int left_node;
					int right_node;
					for (auto ind:boundary_nodes)
					{
						Node boundary_node = m_ma_block_decomposition->get<Node>(ind);
						math::Point R = boundary_node.point()+(-1.)*n.point();
						//double theta = oriented_angle(edge2vec(section,n),vec(R));
						double theta = oriented_angle(tailDir->value(section.id()),vec(R));
						if (theta < 0. && theta > max_theta)
						{
							max_theta = theta;
							right_node = ind;
						}
						if (theta > 0. && theta < min_theta)
						{
							min_theta = theta;
							left_node = ind;
						}
					}
					LTBN->set(section.id(),left_node);
					RTBN->set(section.id(),right_node);
				}
			}

			else
			{
				HN->set(section.id(),medial_node);
				if (nodeType->value(n_id) == 1)
				{
					LHBN->set(section.id(),-1);
					RHBN->set(section.id(),-1);
				}
				if (nodeType->value(n_id) > 1)
				{
					// Find the left and right boundary points
					double max_theta = -10.;
					double min_theta = 10.;
					int left_node;
					int right_node;
					for (auto ind:boundary_nodes)
					{
						Node boundary_node = m_ma_block_decomposition->get<Node>(ind);
						math::Point R = boundary_node.point()+(-1.)*n.point();
						//double theta = oriented_angle(edge2vec(section,n),vec(R));
						double theta = oriented_angle(-headDir->value(section.id()),vec(R));
						if (theta < 0. && theta > max_theta)
						{
							max_theta = theta;
							left_node = ind;
						}
						if (theta > 0. && theta < min_theta)
						{
							min_theta = theta;
							right_node = ind;
						}
					}
					LHBN->set(section.id(),left_node);
					RHBN->set(section.id(),right_node);
				}
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::buildMiddleNodes()
{
	std::cout<<"> Building middle nodes of the block decomposition"<<std::endl;
	auto LHMN = m_topological_representation->getVariable<int,GMDS_EDGE>("leftHeadMiddleNode");
	auto RHMN = m_topological_representation->getVariable<int,GMDS_EDGE>("rightHeadMiddleNode");
	auto LTMN = m_topological_representation->getVariable<int,GMDS_EDGE>("leftTailMiddleNode");
	auto RTMN = m_topological_representation->getVariable<int,GMDS_EDGE>("rightTailMiddleNode");
	//auto TN = m_topological_representation->getVariable<int,GMDS_EDGE>("tailNode");
	//auto HN = m_topological_representation->getVariable<int,GMDS_EDGE>("headNode");
	auto LTBN = m_topological_representation->getVariable<int,GMDS_EDGE>("leftTailBoundaryNode");
	auto RTBN = m_topological_representation->getVariable<int,GMDS_EDGE>("rightTailBoundaryNode");
	auto LHBN = m_topological_representation->getVariable<int,GMDS_EDGE>("leftHeadBoundaryNode");
	auto RHBN = m_topological_representation->getVariable<int,GMDS_EDGE>("rightHeadBoundaryNode");
	auto sectionType = m_topological_representation->getVariable<int,GMDS_EDGE>("section_type");
	auto nodeType = m_topological_representation->getVariable<int,GMDS_NODE>("node_type");
	for (auto s_id:m_topological_representation->edges())
	{
		LHMN->set(s_id,-1);
		RHMN->set(s_id,-1);
		LTMN->set(s_id,-1);
		RTMN->set(s_id,-1);
	}
	for (auto s_id:m_topological_representation->edges())
	{
		if (sectionType->value(s_id) == 1)
		{
			Edge section = m_topological_representation->get<Edge>(s_id);
			if (nodeType->value(section.get<Node>()[0].id()) > 1 && nodeType->value(section.get<Node>()[1].id()) > 1)
			{
				// Get the node where the middle points appear
				Node n;
				if (wings(section.get<Node>()[0],section) == 0)
					n = section.get<Node>()[0];
				else
					n = section.get<Node>()[1];
				// Get the neighbouring sections
				//std::vector<Edge> neighbours = neighbouringEdges(section,n);
				std::vector<Edge> neighbours = orderedNeigbourSections(section,n);
				// Build the medial nodes
				if (orientation(n,section) == 1)
				{
					Node ln = m_ma_block_decomposition->get<Node>(LTBN->value(s_id));
					Node rn = m_ma_block_decomposition->get<Node>(RTBN->value(s_id));
					Node newNode;
					if (LTMN->value(s_id) < 0)
					{
						// Placement of the new node
						math::Point dir = ln.point()+(-1.)*n.point();
						dir = (1./vec(dir).norm())*dir;
						double l = (1./sqrt(2.))*section.length();
						if (l > vec((ln.point()+n.point())*(1./2.)+(-1.)*n.point()).norm())
							l = vec((ln.point()+n.point())*(1./2.)+(-1.)*n.point()).norm();
						//newNode = m_ma_block_decomposition->newNode((ln.point()+n.point())*(1./2.));
						newNode = m_ma_block_decomposition->newNode(n.point()+l*dir);
						LTMN->set(s_id,newNode.id());
						// Update the neighbour middle node
						if (neighbours.size() == 1)
						{
							Edge s = neighbours[0];
							if (orientation(n,s) == 1)
								RTMN->set(s.id(),newNode.id());
							if (orientation(n,s) == -1)
								LHMN->set(s.id(),newNode.id());
						}
						if (neighbours.size() > 1)
						{
							Edge s = neighbours[1];
							if (orientation(n,s) == 1)
								RTMN->set(s.id(),newNode.id());
							if (orientation(n,s) == -1)
								LHMN->set(s.id(),newNode.id());
						}
					}
					if (RTMN->value(s_id) < 0)
					{
						// Placement of the new node
						math::Point dir = rn.point()+(-1.)*n.point();
						dir = (1./vec(dir).norm())*dir;
						double l = (1./sqrt(2.))*section.length();
						if (l > vec((rn.point()+n.point())*(1./2.)+(-1.)*n.point()).norm())
							l = vec((rn.point()+n.point())*(1./2.)+(-1.)*n.point()).norm();
						//newNode = m_ma_block_decomposition->newNode((rn.point()+n.point())*(1./2.));
						newNode = m_ma_block_decomposition->newNode(n.point()+l*dir);
						RTMN->set(s_id,newNode.id());
						// Update the neighbour middle node
						if (neighbours.size() == 1)
						{
							Edge s = neighbours[0];
							if (orientation(n,s) == 1)
								LTMN->set(s.id(),newNode.id());
							if (orientation(n,s) == -1)
								RHMN->set(s.id(),newNode.id());
						}
						if (neighbours.size() > 1)
						{
							Edge s = neighbours[0];
							if (orientation(n,s) == 1)
								LTMN->set(s.id(),newNode.id());
							if (orientation(n,s) == -1)
								RHMN->set(s.id(),newNode.id());
						}
					}

				}
				else // if (orientation(n,section) == -1)
				{
					Node ln = m_ma_block_decomposition->get<Node>(LHBN->value(s_id));
					Node rn = m_ma_block_decomposition->get<Node>(RHBN->value(s_id));
					Node newNode;
					if (LHMN->value(s_id) < 0)
					{
						// Placement of the new node
						math::Point dir = ln.point()+(-1.)*n.point();
						dir = (1./vec(dir).norm())*dir;
						double l = (1./sqrt(2.))*section.length();
						if (l > vec((ln.point()+n.point())*(1./2.)+(-1.)*n.point()).norm())
							l = vec((ln.point()+n.point())*(1./2.)+(-1.)*n.point()).norm();
						//newNode = m_ma_block_decomposition->newNode((ln.point()+n.point())*(1./2.));
						newNode = m_ma_block_decomposition->newNode(n.point()+l*dir);
						LHMN->set(s_id,newNode.id());
						// Update the neighbour middle node
						if (neighbours.size() == 1)
						{
							Edge s = neighbours[0];
							if (orientation(n,s) == 1)
								LTMN->set(s.id(),newNode.id());
							if (orientation(n,s) == -1)
								RHMN->set(s.id(),newNode.id());
						}
						if (neighbours.size() > 1)
						{
							Edge s = neighbours[0];
							if (orientation(n,s) == 1)
								LTMN->set(s.id(),newNode.id());
							if (orientation(n,s) == -1)
								RHMN->set(s.id(),newNode.id());
						}
					}
					if (RHMN->value(s_id) < 0)
					{
						// Placement of the new node
						math::Point dir = rn.point()+(-1.)*n.point();
						dir = (1./vec(dir).norm())*dir;
						double l = (1./sqrt(2.))*section.length();
						if (l > vec((rn.point()+n.point())*(1./2.)+(-1.)*n.point()).norm())
							l = vec((rn.point()+n.point())*(1./2.)+(-1.)*n.point()).norm();
						//newNode = m_ma_block_decomposition->newNode((rn.point()+n.point())*(1./2.));
						newNode = m_ma_block_decomposition->newNode(n.point()+l*dir);
						RHMN->set(s_id,newNode.id());
						// Update the neighbour middle node
						if (neighbours.size() == 1)
						{
							Edge s = neighbours[0];
							if (orientation(n,s) == 1)
								RTMN->set(s.id(),newNode.id());
							if (orientation(n,s) == -1)
								LHMN->set(s.id(),newNode.id());
						}
						if (neighbours.size() > 1)
						{
							Edge s = neighbours[1];
							if (orientation(n,s) == 1)
								RTMN->set(s.id(),newNode.id());
							if (orientation(n,s) == -1)
								LHMN->set(s.id(),newNode.id());
						}
					}
				}
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2D::buildBlocks()
{
	auto TN = m_topological_representation->getVariable<int,GMDS_EDGE>("tailNode");
	auto HN = m_topological_representation->getVariable<int,GMDS_EDGE>("headNode");
	auto LTBN = m_topological_representation->getVariable<int,GMDS_EDGE>("leftTailBoundaryNode");
	auto RTBN = m_topological_representation->getVariable<int,GMDS_EDGE>("rightTailBoundaryNode");
	auto LHBN = m_topological_representation->getVariable<int,GMDS_EDGE>("leftHeadBoundaryNode");
	auto RHBN = m_topological_representation->getVariable<int,GMDS_EDGE>("rightHeadBoundaryNode");
	auto LHMN = m_topological_representation->getVariable<int,GMDS_EDGE>("leftHeadMiddleNode");
	auto RHMN = m_topological_representation->getVariable<int,GMDS_EDGE>("rightHeadMiddleNode");
	auto LTMN = m_topological_representation->getVariable<int,GMDS_EDGE>("leftTailMiddleNode");
	auto RTMN = m_topological_representation->getVariable<int,GMDS_EDGE>("rightTailMiddleNode");
	auto sectionType = m_topological_representation->getVariable<int,GMDS_EDGE>("section_type");
	auto nodeType = m_topological_representation->getVariable<int,GMDS_NODE>("node_type");

	int tn, hn, ltbn, rtbn, lhbn, rhbn, lhmn, rhmn, ltmn, rtmn;
	for (auto s_id:m_topological_representation->edges())
	{
		if (sectionType->value(s_id) == 0)
		{
			tn = TN->value(s_id);
			hn = HN->value(s_id);
			ltbn = LTBN->value(s_id);
			rtbn = RTBN->value(s_id);
			lhbn = LHBN->value(s_id);
			rhbn = RHBN->value(s_id);
			lhmn = LHMN->value(s_id);
			rhmn = RHMN->value(s_id);
			ltmn = LTMN->value(s_id);
			rtmn = RTMN->value(s_id);

			std::vector<TCellID> block1;
			block1.push_back(hn);
			if (lhmn >= 0)
				block1.push_back(lhmn);
			block1.push_back(lhbn);
			block1.push_back(ltbn);
			if (ltmn >= 0)
				block1.push_back(ltmn);
			block1.push_back(tn);
			m_ma_block_decomposition->newFace(block1);

			std::vector<TCellID> block2;
			block2.push_back(tn);
			if (rtmn >= 0)
				block2.push_back(rtmn);
			block2.push_back(rtbn);
			block2.push_back(rhbn);
			if (rhmn >= 0)
				block2.push_back(rhmn);
			block2.push_back(hn);
			m_ma_block_decomposition->newFace(block2);
		}

		if (sectionType->value(s_id) == 1)
		{
			tn = TN->value(s_id);
			hn = HN->value(s_id);
			ltbn = LTBN->value(s_id);
			rtbn = RTBN->value(s_id);
			lhbn = LHBN->value(s_id);
			rhbn = RHBN->value(s_id);
			lhmn = LHMN->value(s_id);
			rhmn = RHMN->value(s_id);
			ltmn = LTMN->value(s_id);
			rtmn = RTMN->value(s_id);
			Edge section = m_topological_representation->get<Edge>(s_id);

			if (nodeType->value(section.get<Node>()[0].id()) == 1 || nodeType->value(section.get<Node>()[1].id()) == 1)
			{
				std::vector<TCellID> block;
				block.push_back(hn);
				if (lhmn >= 0)
					block.push_back(lhmn);
				block.push_back(lhbn+ltbn+1);
				if (ltmn >= 0)
					block.push_back(ltmn);
				block.push_back(tn);
				if (rtmn >= 0)
					block.push_back(rtmn);
				block.push_back(rhbn+rtbn+1);
				if (rhmn >= 0)
					block.push_back(rhmn);
				m_ma_block_decomposition->newFace(block);
			}

			else
			{
				if (wings(section.get<Node>()[0],section) == 1)
				{
					m_ma_block_decomposition->newQuad(hn,lhmn,tn,rhmn);
					std::vector<TCellID> block;
					block.push_back(tn);
					block.push_back(lhmn);
					block.push_back(lhbn);
					block.push_back(ltbn);
					if (ltmn >= 0)
						block.push_back(ltmn);
					m_ma_block_decomposition->newFace(block);
					block.clear();
					block.push_back(tn);
					if (rtmn >= 0)
						block.push_back(rtmn);
					block.push_back(rtbn);
					block.push_back(rhbn);
					block.push_back(rhmn);
					m_ma_block_decomposition->newFace(block);
				}
				else
				{
					m_ma_block_decomposition->newQuad(hn,ltmn,tn,rtmn);
					std::vector<TCellID> block;
					block.push_back(hn);
					if (lhmn >= 0)
						block.push_back(lhmn);
					block.push_back(lhbn);
					block.push_back(ltbn);
					block.push_back(ltmn);
					m_ma_block_decomposition->newFace(block);
					block.clear();
					block.push_back(hn);
					block.push_back(rtmn);
					block.push_back(rtbn);
					block.push_back(rhbn);
					if (rhmn >= 0)
						block.push_back(rhmn);
					m_ma_block_decomposition->newFace(block);
					// m_ma_block_decomposition->newQuad(hn,lhbn,ltbn,ltmn);
					// m_ma_block_decomposition->newQuad(hn,rtmn,rtbn,rhbn);
				}
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
std::vector<Edge> MedialAxis2D::sortedAdjacentEdges(gmds::Node &AN)
{
	std::vector<Edge> adj_edges = AN.get<Edge>();
	std::vector<Edge> sorted_edges;
	math::Vector X;
	X.setX(1.);
	X.setY(0.);
	X.setZ(0.);
	while (!adj_edges.empty())
	{
		int min_pos;
		double min_angle = 10.;
		double angle;
		for (int i = 0 ; i < adj_edges.size() ; i++)
		{
			Edge e = adj_edges[i];
			angle = oriented_angle(X,edge2vec(e,AN));
			if (angle < min_angle)
			{
				min_angle = angle;
				min_pos = i;
			}
		}
		sorted_edges.push_back(adj_edges[min_pos]);
		adj_edges.erase(adj_edges.begin() + min_pos);
	}
	return sorted_edges;
}
/*----------------------------------------------------------------------------*/
std::vector<Edge> MedialAxis2D::neighbouringEdges(gmds::Edge &AE, gmds::Node &AN)
{
	std::vector<Edge> neighbours;
	std::vector<Edge> sorted_edges = sortedAdjacentEdges(AN);

	int pos = -1;
	for (int i = 0; i < sorted_edges.size(); i++)
	{
		if (sorted_edges[i].id() == AE.id())
		{
			pos = i;
			break;
		}
	}
	if (pos == -1)
	{
		std::cout<<"neighbouringEdges() : Warning, the given node doesn't belong to the given edge"<<std::endl;
		return neighbours;
	}
	if (sorted_edges.size() == 2)
	{
		neighbours.push_back(sorted_edges[1-pos]);
		return neighbours;
	}
	if (pos > 0 && pos < sorted_edges.size() - 1)
	{
		neighbours.push_back(sorted_edges[pos-1]);
		neighbours.push_back(sorted_edges[pos+1]);
		return neighbours;
	}
	if (pos == 0)
	{
		neighbours.push_back(sorted_edges[sorted_edges.size() - 1]);
		neighbours.push_back(sorted_edges[1]);
		return neighbours;
	}
	if (pos == sorted_edges.size() - 1)
	{
		neighbours.push_back(sorted_edges[pos-1]);
		neighbours.push_back(sorted_edges[0]);
		return neighbours;
	}
	return neighbours;

}

/*----------------------------------------------------------------------------*/
std::vector<Edge> MedialAxis2D::orderedNeigbourSections(Edge &ASection, Node &AN)
{
	auto tailDir = m_topological_representation->getVariable<math::Vector,GMDS_EDGE>("tailDirection");
	auto headDir = m_topological_representation->getVariable<math::Vector,GMDS_EDGE>("headDirection");
	
	std::vector<Edge> adj_edges = AN.get<Edge>();
	std::vector<double> angles;
	math::Vector X;
	X.setX(1.);
	X.setY(0.);
	X.setZ(0.);
	for (auto s:adj_edges)
	{
		double angle;
		if (orientation(AN,s) == 1)
			angle = oriented_angle(X,tailDir->value(s.id()));
		else
			angle = oriented_angle(X,-headDir->value(s.id()));
		angles.push_back(angle);
	}
	std::vector<Edge> neighbours;
	std::vector<Edge> sorted_edges = order(adj_edges,angles);

	int pos = -1;
	for (int i = 0; i < sorted_edges.size(); i++)
	{
		if (sorted_edges[i].id() == ASection.id())
		{
			pos = i;
			break;
		}
	}
	if (pos == -1)
	{
		std::cout<<"neighbouringEdges() : Warning, the given node doesn't belong to the given edge"<<std::endl;
		return neighbours;
	}
	if (sorted_edges.size() == 2)
	{
		neighbours.push_back(sorted_edges[1-pos]);
		return neighbours;
	}
	if (pos > 0 && pos < sorted_edges.size() - 1)
	{
		neighbours.push_back(sorted_edges[pos-1]);
		neighbours.push_back(sorted_edges[pos+1]);
		return neighbours;
	}
	if (pos == 0)
	{
		neighbours.push_back(sorted_edges[sorted_edges.size() - 1]);
		neighbours.push_back(sorted_edges[1]);
		return neighbours;
	}
	if (pos == sorted_edges.size() - 1)
	{
		neighbours.push_back(sorted_edges[pos-1]);
		neighbours.push_back(sorted_edges[0]);
		return neighbours;
	}

	return neighbours;

}


/*----------------------------------------------------------------------------*/
}  // end namespace medialaxis
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/