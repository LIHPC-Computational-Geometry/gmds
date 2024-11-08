//
// Created by chenyt on 22/04/24.
//

/*----------------------------------------------------------------------------*/
#include "gmds/medialaxis/MedialAxis2DBuilder.h"
#include "gmds/medialaxis/MedialAxisMath.h"
#include "gmds/ig/MeshDoctor.h"
#include <gmds/math/Triangle.h>
#include <gmds/math/Point.h>
#include <stack>
#include <queue>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace medialaxis {
/*----------------------------------------------------------------------------*/
MedialAxis2DBuilder::MedialAxis2DBuilder(Mesh &AMesh){
	m_mesh = &AMesh;
	// Delaunay triangle / medial point correspondence
	m_mesh->newVariable<int,GMDS_FACE>("triangle2MP");
	// Delaunay internal edge / medial edge correspondence
	m_mesh->newVariable<int,GMDS_EDGE>("intEdge2MedEdge");
	// Edges classification between boundary (1) and internal (0) edges
	m_mesh->newVariable<int,GMDS_EDGE>("edge_class");
	// Edges classification between corner (1) and non corner (0) edges
	m_mesh->newVariable<int,GMDS_EDGE>("corner");
	// Faces classification (tangency = discrete number of tangency points of the circumcenter)
	m_mesh->newVariable<int,GMDS_FACE>("tangency");
	// Points classification (1 if the point is considered a detail, 0 if not)
	m_mesh->newVariable<int,GMDS_NODE>("is_detail");
	// Edge classification (1 if the edge is considered a detail, 0 if not)
	m_mesh->newVariable<int,GMDS_EDGE>("edge_is_detail");
	// Geometry sides Ids
	m_mesh->newVariable<int,GMDS_EDGE>("side_id");
	// Boundary connected component ID
	m_mesh->newVariable<int,GMDS_NODE>("boundary_connected_component_id");
	// Boundary connected components IDs
	m_mesh->newVariable<std::vector<int>,GMDS_FACE>("boundary_connected_components_ids");

	// Voronoï medial axis
	m_voronoi_medax = new MedialAxis2D();
	// Smoothed and refined medial axis
	m_smoothed_medax = new MedialAxis2D();

	// Graph of the boundary connected components
	m_boundary_connected_components_graph = new Mesh(MeshModel(DIM3 | F | E | N |
	                                                  F2N | E2N | N2E | N2F));
	// Primal triangle corresponding to the optimal medial point connecting two components
	m_boundary_connected_components_graph->newVariable<int,GMDS_EDGE>("optimal_primal_triangle");
}
/*----------------------------------------------------------------------------*/
MedialAxis2DBuilder::~MedialAxis2DBuilder()
{
	if(m_voronoi_medax != nullptr)
		delete m_voronoi_medax;
}
/*----------------------------------------------------------------------------*/
MedialAxis2D * MedialAxis2DBuilder::getMedialObject()
{
	return m_voronoi_medax;
}

/*----------------------------------------------------------------------------*/
MedialAxis2D * MedialAxis2DBuilder::getSmoothedMedialObject()
{
	return m_smoothed_medax;
}

/*-------------------------------------------------------------------------*/
Mesh* MedialAxis2DBuilder::getMesh()
{
	return m_mesh;
}

/*----------------------------------------------------------------------------*/
void MedialAxis2DBuilder::classEdges()
{
	auto var = m_mesh->getVariable<int,GMDS_EDGE>("edge_class");
	for (auto n_id:m_mesh->edges()){
		Edge e = m_mesh->get<Edge>(n_id);
		std::vector<Node> ne = e.get<Node>();
		std::vector<TCellID> adj_faces = m_mesh->getCommonFaces(ne[0],ne[1]);
		// Boundary edges are marked with 1
		if (adj_faces.size() == 1){
			var->set(e.id(),1);
		}
		// Interior edges are marked with 0
		else {
			var->set(e.id(), 0);
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2DBuilder::classFaces()
{
	this->classEdges();
	auto var_edge = m_mesh->getVariable<int,GMDS_EDGE>("edge_class");
	auto var_face = m_mesh->getVariable<int,GMDS_FACE>("tangency");
	for (auto n_id:m_mesh->faces()){
		Face f = m_mesh->get<Face>(n_id);
		std::vector<Edge> ef = f.get<Edge>();
		int tangency = 3 - var_edge->value(ef[0].id()) - var_edge->value(ef[1].id()) - var_edge->value(ef[2].id());
		var_face->set(f.id(), tangency);
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2DBuilder::markCorners()
{
	auto corners = m_mesh->getVariable<int,GMDS_EDGE>("corner");
	for (auto e_id:m_mesh->edges())
	{
		Edge e = m_mesh->get<Edge>(e_id);
		// If it is a boundary edge
		if (e.get<Face>().size() == 1)
		{
			Face f = e.get<Face>()[0];
			int NbNeighbours = 0;
			for (auto edge:f.get<Edge>())
				NbNeighbours += (edge.get<Face>().size() - 1);
			if (NbNeighbours == 1)
				corners->set(e_id,1);
			else
				corners->set(e_id,0);
		}
		else
			corners->set(e_id,0);
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2DBuilder::setSideId()
{
	std::cout<<"> Setting geometry sides IDs"<<std::endl;
	auto corners = m_mesh->getVariable<int,GMDS_EDGE>("corner");
	auto sideId = m_mesh->getVariable<int,GMDS_EDGE>("side_id");
	// Initialize side IDs at -1 (so that internal edges belong to side -1)
	for (auto e_id:m_mesh->edges())
		sideId->set(e_id,-1);
	auto alreadyVisited = m_mesh->newVariable<int,GMDS_EDGE>("alreadyVisited");
	// Mark internal edges as already visited
	for (auto e_id:m_mesh->edges())
	{
		Edge e = m_mesh->get<Edge>(e_id);
		if (e.get<Face>().size() == 2)
			alreadyVisited->set(e_id,1);
		else
			alreadyVisited->set(e_id,0);
	}
	// Sides IDs
	int Id = 0;
	for (auto e_id:m_mesh->edges())
	{
		if ((alreadyVisited->value(e_id) == 0) && (corners->value(e_id) == 0))
		{
			std::stack<TCellID> side;
			side.push(e_id);
			alreadyVisited->set(e_id,1);
			while (!side.empty())
			{
				Edge e = m_mesh->get<Edge>(side.top());
				side.pop();
				sideId->set(e.id(),Id);
				if (corners->value(e.id()) == 0)
				{
					for (auto n:e.get<Node>())
					{
						for (auto edge:n.get<Edge>())
						{
							if (alreadyVisited->value(edge.id()) == 0)
							{
								side.push(edge.id());
								alreadyVisited->set(edge.id(),1);
							}
						}
					}
				}
			}
			Id += 1;
		}
	}
	// Find the forgotten sides (those made of only one edge)
	for (auto e_id:m_mesh->edges())
	{
		if (alreadyVisited->value(e_id) == 0)
		{
			sideId->set(e_id,Id);
			Id += 1;
		}
	}
	std::cout<<"NB sides : "<<Id<<std::endl;
	m_mesh->deleteVariable(GMDS_EDGE, alreadyVisited);
}

/*----------------------------------------------------------------------------*/
void MedialAxis2DBuilder::setMedialRadiiOrthogonalityDefault(const double &ABoundaryCurvatureTol)
{
	auto triangle2MP = m_mesh->getVariable<int,GMDS_FACE>("triangle2MP");
	for (auto t_id:m_mesh->faces())
	{
		Face tri = m_mesh->get<Face>(t_id);
		TCellID n_id = triangle2MP->value(t_id);
		Node n = m_voronoi_medax->getMedPoint(n_id);
		int type = 0;
		for (auto e:tri.get<Edge>())
			type += e.get<Face>().size() - 1;
		if (type == 2)
		{
			Node A, B;
			math::Point P1, P2, E1, E2, E3;
			// Get the boundary edge of tri
			for (auto e:tri.get<Edge>())
			{
				if (e.get<Face>().size() == 1)
				{
					A = e.get<Node>()[0];
					B = e.get<Node>()[1];
					E1 = B.point() + (-1.)*A.point();
					P1 = (1./2.)*A.point() + (1./2.)*B.point();
					break;
				}
			}
			// Get the alone point of tri and its adjacent boundary edges
			std::vector<Edge> adj_bound_edges;
			for (auto n1:tri.get<Node>())
			{
				if (n1.id() != A.id()  && n1.id() != B.id())
				{
					P2 = n1.point();
					for (auto e:n1.get<Edge>())
					{
						if (e.get<Face>().size() == 1)
						{
							adj_bound_edges.push_back(e);
						}
					}
					break;
				}
			}
			if (adj_bound_edges.size() != 2)
				break;
			E2 = adj_bound_edges[0].get<Node>()[1].point() + (-1.)*adj_bound_edges[0].get<Node>()[0].point();
			E3 = adj_bound_edges[1].get<Node>()[1].point() + (-1.)*adj_bound_edges[1].get<Node>()[0].point();
			// Normalized medial radii
			math::Point R1 = P1 + (-1.)*n.point();
			R1 = R1*(1./vec(R1).norm());
			math::Point R2 = P2 + (-1.)*n.point();
			R2 = R2*(1./vec(R2).norm());
			// Normalized boundary edges
			E1 = E1*(1./vec(E1).norm());
			E2 = E2*(1./vec(E2).norm());
			E3 = E3*(1./vec(E3).norm());
			// Compute the inner products between the medial radii and the boundary
			double ip1 = fabs(vec(R1).dot(vec(E1)));
			double ip2 = fabs(vec(R2).dot(vec(E2)));
			double ip3 = fabs(vec(R2).dot(vec(E3)));
			double max_ip = ip1;
			if (ip2 > max_ip)
				max_ip = ip2;
			if (ip3 > max_ip)
				max_ip = ip3;
			// If the medial point correspond to a concave corner of the geometry (ie E2 and E3 are not aligned)
			// then the tangency of the medial radius doesn't have sense
			if (fabs(vec(E2).dot(vec(E3))-1) > ABoundaryCurvatureTol && fabs(vec(E2).dot(vec(E3))+1) > ABoundaryCurvatureTol)
				max_ip = 0.;
			m_voronoi_medax->setMedialRadiusOrthogonalityDefault(n_id,max_ip);
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2DBuilder::setTouchingPoints()
{
	auto triangle2MP = m_mesh->getVariable<int,GMDS_FACE>("triangle2MP");
	for (auto f_id:m_mesh->faces())
	{
		Face face = m_mesh->get<Face>(f_id);
		TCellID medPoint = triangle2MP->value(face.id());
		std::vector<math::Point> touchingPoints;
		std::vector<Edge> adj_edges = face.get<Edge>();
		// Medial point type computation
		int type = 0;
		for (auto e:adj_edges)
			type += e.get<Face>().size() -1;
		// If it is an end point, it has one touching point, approximated by its triangle's barycenter
		if (type == 1)
		{
			math::Point P1 = face.get<Node>()[0].point();
			math::Point P2 = face.get<Node>()[1].point();
			math::Point P3 = face.get<Node>()[2].point();
			math::Point P = (P1+P2+P3)*(1./3.);
			touchingPoints.push_back(P);
		}
		// If it is an intersection point, all its triangle points are touching points
		if (type >= 3)
		{
			for (auto P:face.get<Node>())
				touchingPoints.push_back(P.point());
		}
		// If it is a regular point, it has two touching points, approximated by the barycenter of the two
		// neighbouring points and the remaining point of its triangle
		if (type == 2)
		{
			// Find the only boundary edge of the triangle
			Edge e;
			for (auto edge:face.get<Edge>())
			{
				if (edge.get<Face>().size() == 1)
				{
					e = edge;
					break;
				}
			}
			math::Point P1 = (e.get<Node>()[0].point()+e.get<Node>()[1].point())*(1./2.);
			touchingPoints.push_back(P1);
			math::Point P2;
			for (auto n:face.get<Node>())
			{
				if (n.id() != e.get<Node>()[0].id() && n.id() != e.get<Node>()[1].id())
					P2 = n.point();
			}
			touchingPoints.push_back(P2);
		}
		m_voronoi_medax->setTouchingPoints(medPoint,touchingPoints);
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2DBuilder::setMedialRadius()
{
	// Requires the medial axis already constructed
	auto medial_point = m_mesh->getVariable<int,GMDS_FACE>("triangle2MP");
	for (auto t_id:m_mesh->faces())
	{
		Face Triangle = m_mesh->get<Face>(t_id);
		std::vector<Node> nf = Triangle.get<Node>();
		math::Triangle T(nf[0].point(), nf[1].point(), nf[2].point());
		math::Point MP = T.getCircumcenter();
		double r = MP.distance(nf[0].point());
		m_voronoi_medax->setMedialRadius(medial_point->value(Triangle.id()),r);
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2DBuilder::markGeometryDetailsOnPoints()
{
	auto medial_point = m_mesh->getVariable<int,GMDS_FACE>("triangle2MP");
	auto isDetail = m_mesh->getVariable<int,GMDS_NODE>("is_detail");
	for (auto n_id:m_mesh->nodes())
		isDetail->set(n_id,0);
	for (auto t_id:m_mesh->faces())
	{
		TCellID medPoint = medial_point->value(t_id);
		if (m_voronoi_medax->correspondsToDetail(medPoint) == 1)
		{
			Face tri = m_mesh->get<Face>(t_id);
			for (auto n:tri.get<Node>())
			{
				isDetail->set(n.id(),1);
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2DBuilder::markGeometryDetailsOnEdges()
{
	// Requires details marked on points
	auto edgeIsDetail = m_mesh->getVariable<int,GMDS_EDGE>("edge_is_detail");
	auto isDetail = m_mesh->getVariable<int,GMDS_NODE>("is_detail");
	for (auto e_id:m_mesh->edges())
		edgeIsDetail->set(e_id,0);
	for (auto e_id:m_mesh->edges())
	{
		Edge e = m_mesh->get<Edge>(e_id);
		Node n1 = e.get<Node>()[0];
		Node n2 = e.get<Node>()[1];
		// An edge is considered as detail if its two nodes are, and if it is a boundary edge
		if (isDetail->value(n1.id()) == 1 && isDetail->value(n2.id()) == 1 && e.get<Face>().size() == 1)
			edgeIsDetail->set(e.id(),1);
	}
}

/*----------------------------------------------------------------------------*/
double MedialAxis2DBuilder::maxDelEdgeLength()
{
	double maxLength = 0.;
	for (auto e_id:m_mesh->edges())
	{
		Edge e = m_mesh->get<Edge>(e_id);
		if (e.length() > maxLength)
			maxLength = e.length();
	}
	return maxLength;
}

/*----------------------------------------------------------------------------*/
void MedialAxis2DBuilder::buildVoronoiMedialPointsAndEdges()
{
	// Delaunay triangle / medial point correspondence
	auto medial_point = m_mesh->getVariable<int,GMDS_FACE>("triangle2MP");
	auto edge_correspondance = m_mesh->getVariable<int,GMDS_EDGE>("intEdge2MedEdge");
	// Adding the medial points
	for (auto n_id:m_mesh->faces()){
		Face f = m_mesh->get<Face>(n_id);
		std::vector<Node> nf = f.get<Node>();
		math::Triangle T(nf[0].point(), nf[1].point(), nf[2].point());
		math::Point MP = T.getCircumcenter();
		Node n = m_voronoi_medax->newMedPoint(MP);
		medial_point->set(f.id(), n.id());
		m_voronoi_medax->setPrimalTriangleID(n.id(),f.id());
	}
	// Create the medial edges
	for (auto e_id:m_mesh->edges())
	{
		Edge e = m_mesh->get<Edge>(e_id);
		std::vector<Face> adj_faces = e.get<Face>();
		if (adj_faces.size() == 2)
		{
			Face f1 = adj_faces[0];
			Face f2 = adj_faces[1];
			TCellID mp_id1 = medial_point->value(f1.id());
			TCellID mp_id2 = medial_point->value(f2.id());
			Edge me = m_voronoi_medax->newMedEdge(mp_id1, mp_id2);
			edge_correspondance->set(e.id(),me.id());
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2DBuilder::buildVoronoiMedialAxis()
{
	std::cout<<" "<<std::endl;
	std::cout<<"===== Building the Voronoï medial axis ====="<<std::endl;
	// Build points and edges
	buildVoronoiMedialPointsAndEdges();

	// Attach to each medial points its corresponding boundary points
	setTouchingPoints();

	// Set medial radii
	setMedialRadius();

	// Set medial radii orthogonality default
	setMedialRadiiOrthogonalityDefault(0.05);

	// Update the medial axis connectivity
	m_voronoi_medax->updateConnectivity();

	// Set the medial points type values
	m_voronoi_medax->setMedialPointType();

	// Set the medial edges type values
	m_voronoi_medax->setMedialEdgeType();

	// Set the branches IDs on the edges
	m_voronoi_medax->setBranchIdOnEdges();

	// Set the branches IDs on the points
	m_voronoi_medax->setBranchIdOnPoints();

	// Set the branches type on the edges
	m_voronoi_medax->setBranchTypeOnEdges();

	// Set the branches type on the points
	m_voronoi_medax->setBranchTypeOnPoints();

	std::cout<<"NB medial points : "<<m_voronoi_medax->getNbMedPoints()<<std::endl;
	std::cout<<"Min medial edge length : "<<m_voronoi_medax->minMedEdgeLength()<<std::endl;
	std::cout<<"Max medial edge length : "<<m_voronoi_medax->maxMedEdgeLength()<<std::endl;
	std::cout<<"Mean medial edge length : "<<m_voronoi_medax->meanMedEdgeLength()<<std::endl;

	std::cout<<"============================================"<<std::endl;
	std::cout<<" "<<std::endl;
}

/*----------------------------------------------------------------------------*/
void MedialAxis2DBuilder::buildSmoothedMedaxFromVoronoi(const std::vector<std::vector<Node>> AGroups)
{
	std::cout<<" "<<std::endl;
	std::cout<<"===== Building the smoothed medial axis ====="<<std::endl;
	// Adjacency between the groups
	std::vector<std::vector<int>> setsOfNeighbours;
	for (int groupID=0; groupID<AGroups.size(); groupID++)
	{
		std::vector<int> neighbours;
		for (auto n:AGroups[groupID])
		{
			for (auto e:n.get<Edge>())
			{
				for (auto m:e.get<Node>())
				{
					int groupID2 = m_voronoi_medax->medialPoint2Group(m);
					if (groupID2 < groupID)
					{
						// Check if we already saw this group
						bool alreadySeen = false;
						for (auto groupID3:neighbours)
						{
							if (groupID3 == groupID2)
							{
								alreadySeen = true;
								break;
							}
						}
						// If not add it to the neighbours list
						if (!alreadySeen)
							neighbours.push_back(groupID2);
					}
				}
			}
		}
		setsOfNeighbours.push_back(neighbours);
	}
	// Groups/new medial points correspondence
	std::vector<TCellID> group2NewMedPoint(AGroups.size());
	// Build the medial points
	for (int groupID=0; groupID<AGroups.size(); groupID++)
	{
		Node newPoint = m_smoothed_medax->newMedPoint(AGroups[groupID][0].point());
		group2NewMedPoint[groupID] = newPoint.id();
	}
	// Build the medial edges
	for (int groupID=0; groupID<AGroups.size(); groupID++)
	{
		for (auto groupID2:setsOfNeighbours[groupID])
		{
			m_smoothed_medax->newMedEdge(group2NewMedPoint[groupID2],group2NewMedPoint[groupID]);
		}
	}

	// Set the touching points, medial radii and medial radii orthogonality default
	for (int groupID=0; groupID<AGroups.size(); groupID++)
	{
		TCellID newPoint = group2NewMedPoint[groupID];
		TCellID oldPoint = AGroups[groupID][0].id();
		for (auto n:AGroups[groupID])
		{
			if (m_voronoi_medax->getMedialPointType(n.id()) > 2)
				oldPoint = n.id();
		}
		// Computing touching points
		// Nb of intersection points in the group
		int NbIPs = 0;
		for (auto node:AGroups[groupID])
		{
			if (m_voronoi_medax->getTouchingPoints(node.id()).size() > 2)
				NbIPs += 1;
		}
		if (NbIPs == 0)
			m_smoothed_medax->setTouchingPoints(newPoint, m_voronoi_medax->getTouchingPoints(oldPoint));
		if (NbIPs == 1)
		{
			for (auto node:AGroups[groupID])
			{
				if (m_voronoi_medax->getTouchingPoints(node.id()).size() > 2)
					{
						oldPoint = node.id();
						break;
					}
			}
			m_smoothed_medax->setTouchingPoints(newPoint, m_voronoi_medax->getTouchingPoints(oldPoint));
		}
		if (NbIPs >= 2)
		{
			TCellID ip1,ip2;
			for (auto node:AGroups[groupID])
			{
				if (m_voronoi_medax->getTouchingPoints(node.id()).size() > 2)
					{
						ip1 = node.id();
						break;
					}
			}
			for (auto node:AGroups[groupID])
			{
				if (m_voronoi_medax->getTouchingPoints(node.id()).size() > 2 && node.id() != ip1)
					{
						ip2 = node.id();
						break;
					}
			}
			std::vector<math::Point> merged_tangency_points = merge(m_voronoi_medax->getTouchingPoints(ip1),m_voronoi_medax->getTouchingPoints(ip2));
			m_smoothed_medax->setTouchingPoints(newPoint,merged_tangency_points);
		}
		m_smoothed_medax->setMedialRadius(newPoint, m_voronoi_medax->getMedialRadius(oldPoint));
		m_smoothed_medax->setMedialRadiusOrthogonalityDefault(newPoint, m_voronoi_medax->getMedialRadiusOrthogonalityDefault(oldPoint));
	}

	// Refining the medial axis
	double mean_edge_length = m_smoothed_medax->meanMedEdgeLength();
	m_smoothed_medax->refine(2*mean_edge_length);

	// Update the medial axis connectivity
	m_smoothed_medax->updateConnectivity();

	// Set the medial points type values
	m_smoothed_medax->setMedialPointType();

	// Set the medial edges type values
	m_smoothed_medax->setMedialEdgeType();

	// Set the branches IDs on the edges
	m_smoothed_medax->setBranchIdOnEdges();

	// Set the branches IDs on the points
	m_smoothed_medax->setBranchIdOnPoints();

	// Set the branches type on the edges
	m_smoothed_medax->setBranchTypeOnEdges();

	// Set the branches type on the points
	m_smoothed_medax->setBranchTypeOnPoints();

	std::cout<<"NB groups : "<<AGroups.size()<<std::endl;
	std::cout<<"NB medial points : "<<m_smoothed_medax->getNbMedPoints()<<std::endl;
	std::cout<<"Min medial edge length : "<<m_smoothed_medax->minMedEdgeLength()<<std::endl;
	std::cout<<"Max medial edge length : "<<m_smoothed_medax->maxMedEdgeLength()<<std::endl;
	std::cout<<"Mean medial edge length : "<<m_smoothed_medax->meanMedEdgeLength()<<std::endl;

	std::cout<<"============================================"<<std::endl;
	std::cout<<" "<<std::endl;
}

/*----------------------------------------------------------------------------*/
void MedialAxis2DBuilder::identifyDetails()
{
	std::cout<<" "<<std::endl;
	std::cout<<"== Identifying details on the geometry using the medial axis =="<<std::endl;

	// Identify branches corresponding to details
	m_voronoi_medax->identifyDetailsBranches(maxDelEdgeLength()/20.);

	// Mark the details on the points of the input geometry
	markGeometryDetailsOnPoints();

	// Mark the details on the edges of the input geometry
	markGeometryDetailsOnEdges();

	std::cout<<"==============================================================="<<std::endl;
	std::cout<<" "<<std::endl;
}

/*----------------------------------------------------------------------------*/
void MedialAxis2DBuilder::placeSingularities()
{
	std::cout<<" "<<std::endl;
	std::cout<<"== Placing singularities using the medial axis =="<<std::endl;

	// Set the cosine of medial angle
	m_smoothed_medax->setCosMedialAngle();

	// Smooth cos(medial angle) (works approximately)
	std::cout<<"> Smoothing cosine of the medial angle"<<std::endl;
	for (int i = 0; i < 0; i++)
		m_smoothed_medax->smoothCosMedialAngle();

	// Set the optimal cross flux through curves formed by two adjacent medial radii
	m_smoothed_medax->setFluxThroughMedialRadii();

	// Set the flux residual across edges
	m_smoothed_medax->setFluxResidual();

	// Place singularities
	m_smoothed_medax->placeSingularities(5.0);

	// Check singularities
	m_smoothed_medax->checkSingularities(0.1);

	// Move singularities to intersection points
	m_smoothed_medax->moveSingularitiesToIPs(5);

	// Write singularities as attributes
	m_smoothed_medax->setSingularities();

	std::cout<<"==================================================="<<std::endl;
	std::cout<<" "<<std::endl;
}

/*----------------------------------------------------------------------------*/
void MedialAxis2DBuilder::deleteIsolatedPoints()
{
	std::cout<<"> Deleting isolated nodes from the primal Delaunay mesh"<<std::endl;
	int N = 0;
	for (auto n_id:m_mesh->nodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		if (n.get<Edge>().empty())
		{
			m_mesh->deleteNode(n);
			N += 1;
		}
	}
	std::cout<<N<<" isolated nodes deleted"<<std::endl;
}

/*----------------------------------------------------------------------------*/
void MedialAxis2DBuilder::setBoundaryConnectedComponents(int &AN)
{
	std::cout<<"> Identifying the connected components of the boundary"<<std::endl;
	auto IDs = m_mesh->getVariable<int,GMDS_NODE>("boundary_connected_component_id");
	for (auto n_id:m_mesh->nodes())
		IDs->set(n_id,-1);
	int ID = 0;
	for (auto n_id:m_mesh->nodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		if (IDs->value(n_id) == -1)
		{
			std::stack<Node> nodesToAdd;
			IDs->set(n_id,ID);
			nodesToAdd.push(n);
			while (!nodesToAdd.empty())
			{
				Node n2 = nodesToAdd.top();
				nodesToAdd.pop();
				for (auto e:n2.get<Edge>())
				{
					if (e.get<Face>().size() == 1)
					{
						for (auto n3:e.get<Node>())
						{
							if (IDs->value(n3.id()) == -1)
							{
								IDs->set(n3.id(),ID);
								nodesToAdd.push(n3);
							}
						}
					}
				}
			}
			ID += 1;
		}
	}
	AN = ID;
	std::cout<<"NB boundary connected components : "<<ID<<std::endl;
}

/*----------------------------------------------------------------------------*/
void MedialAxis2DBuilder::setBoundaryConnectedComponentsOnFaces()
{
	std::cout<<"> Attaching to each Delaunay triangle its corresponding connected components"<<std::endl;
	auto IDs = m_mesh->getVariable<int,GMDS_NODE>("boundary_connected_component_id");
	auto IDsOnFaces = m_mesh->getVariable<std::vector<int>,GMDS_FACE>("boundary_connected_components_ids");
	for (auto f_id:m_mesh->faces())
	{
		Face f = m_mesh->get<Face>(f_id);
		std::vector<int> ids;
		int id1 = IDs->value(f.get<Node>()[0].id());
		int id2 = IDs->value(f.get<Node>()[1].id());
		int id3 = IDs->value(f.get<Node>()[2].id());
		ids.push_back(id1);
		if (id2 != id1)
			ids.push_back(id2);
		if ((id3 != id1) && (id3 != id2))
			ids.push_back(id3);
		IDsOnFaces->set(f_id,ids);
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2DBuilder::setBoundaryConnectedComponentsOnMedPoints()
{
	std::cout<<"> Attaching to each medial point its corresponding connected components"<<std::endl;
	auto IDsOnFaces = m_mesh->getVariable<std::vector<int>,GMDS_FACE>("boundary_connected_components_ids");
	auto tri2MP = m_mesh->getVariable<int,GMDS_FACE>("triangle2MP");
	for (auto f_id:m_mesh->faces())
	{
		m_voronoi_medax->setBoundaryConnectedComponentsIDs(tri2MP->value(f_id),IDsOnFaces->value(f_id));
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis2DBuilder::buildBoundaryComponentsGraph(Eigen::MatrixXi &AM)
{
	std::cout<<"> Building boundary connected components graph"<<std::endl;
	auto tri = m_boundary_connected_components_graph->getVariable<int,GMDS_EDGE>("optimal_primal_triangle");
	int N = AM.rows();
	for (int I = 0; I < N; I++)
	{
		math::Point P;
		m_boundary_connected_components_graph->newNode(P);
	}
	for (int I = 1; I < N; I++)
	{
		for (int J = 0; J < I; J++)
		{
			if (AM(I,J) >= 0)
			{
				Edge e = m_boundary_connected_components_graph->newEdge(I,J);
				tri->set(e.id(),m_voronoi_medax->primalTriangleID(AM(I,J)));
			}
		}
	}

	// Connectivity of the graph
	MeshDoctor doc(m_boundary_connected_components_graph);
	doc.updateUpwardConnectivity();
}

/*----------------------------------------------------------------------------*/
std::vector<Edge> MedialAxis2DBuilder::connectedComponentsSpanningTree()
{
	auto previous = m_boundary_connected_components_graph->newVariable<int,GMDS_NODE>("previous");
	for (auto  n_id:m_boundary_connected_components_graph->nodes())
		previous->set(n_id,-2);

	std::queue<Node> front;
	std::vector<Edge> tree;

	// Root of the tree
	Node n = m_boundary_connected_components_graph->get<Node>(0);
	front.push(n);
	previous->set(0,-1);
	// Build the tree
	while (!front.empty())
	{
		n = front.front();
		front.pop();
		for (auto e:n.get<Edge>())
		{
			for (auto n2:e.get<Node>())
			{
				if (previous->value(n2.id()) == -2)
				{
					front.push(n2);
					previous->set(n2.id(),n.id());
					tree.push_back(e);
				}
			}
		}
	}
	m_boundary_connected_components_graph->deleteVariable(GMDS_NODE,previous);
	return tree;
}

/*----------------------------------------------------------------------------*/
std::vector<std::vector<math::Point>> MedialAxis2DBuilder::connectedBoundaryPoints(std::vector<Edge> &AST)
{
	std::cout<<"> Computing the connected boundary points"<<std::endl;
	auto tri = m_boundary_connected_components_graph->getVariable<int,GMDS_EDGE>("optimal_primal_triangle");
	std::vector<std::vector<math::Point>> connectedPoints;
	for (auto e:AST)
	{
		int f_id = tri->value(e.id());
		Face f = m_mesh->get<Face>(f_id);
		Node n1, n2, n3;
		for (auto e:f.get<Edge>())
		{
			if (e.get<Face>().size() == 1)
			{
				n2 = e.get<Node>()[0];
				n3 = e.get<Node>()[1];
				break;
			}
		}
		for (auto n:f.get<Node>())
		{
			if ((n.id() != n2.id()) && (n.id() != n3.id()))
			{
				n1 = n;
				break;
			}
		}
		std::vector<math::Point> points;
		points.push_back(n1.point());
		points.push_back((n2.point()+n3.point())*(1./2.));
		connectedPoints.push_back(points);
	}

	return connectedPoints;
}

/*----------------------------------------------------------------------------*/
void MedialAxis2DBuilder::connectBoundaryConnectedComponents()
{
	std::cout<<" "<<std::endl;
	std::cout<<"== Connecting connected components of the boundary =="<<std::endl;

	// Delete isolated points
	deleteIsolatedPoints();

	// Set the boundary connected components IDs
	int N; // NB of connected components
	setBoundaryConnectedComponents(N);
	// Attach to each primal triangle and its dual medial point its adjacent connected components IDs
	setBoundaryConnectedComponentsOnFaces();
	setBoundaryConnectedComponentsOnMedPoints();

	// Set the cosine of medial angle
	m_voronoi_medax->setCosMedialAngle();

	// Find optimal boundary connected components
	Eigen::MatrixXi optimalMedPoint = m_voronoi_medax->findOptimalBoundaryConnectedComponentsConnexion(N);

	// Build boundary components graph
	buildBoundaryComponentsGraph(optimalMedPoint);

	// Spanning tree of the graph
	std::vector<Edge> tree = connectedComponentsSpanningTree();

	// Connected points
	std::vector<std::vector<math::Point>> connectedPoints = connectedBoundaryPoints(tree);
	m_connected_boundary_points = connectedPoints;

	std::cout<<"==================================================="<<std::endl;
	std::cout<<" "<<std::endl;
}


/*----------------------------------------------------------------------------*/
MedialAxis2DBuilder::STATUS

MedialAxis2DBuilder::execute()
{
	// Mark corners
	markCorners();

	// Set geometry sides IDs
	setSideId();

	// Voronoi medial axis
	buildVoronoiMedialAxis();

	// Smoothed medial axis
	std::vector<std::vector<Node>> groups = m_voronoi_medax->medialPointsGroups(m_voronoi_medax->meanMedEdgeLength()/2.);
	buildSmoothedMedaxFromVoronoi(groups);

	// Identify details on geometry using the (Voronoï) medial axis
	//identifyDetails();

	// Place singularities using the (smoothed) medial axis
	placeSingularities();

	// Connect the connected components of the boundary using the medial axis
	connectBoundaryConnectedComponents();


	return MedialAxis2DBuilder::SUCCESS;
}
/*----------------------------------------------------------------------------*/
}  // end namespace medialaxis
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
