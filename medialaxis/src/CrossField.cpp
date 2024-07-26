//
// Created by chenyt on 10/07/24.
//

#include "gmds/medialaxis/CrossField.h"
#include <stack>
#include <queue>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace medialaxis {
/*----------------------------------------------------------------------------*/
CrossField::CrossField(Mesh &AMesh)
{
	m_mesh = &AMesh;
	// Boundary connected component ID
	m_mesh->newVariable<int,GMDS_NODE>("boundary_connected_component_id");
	// Edges connecting boundaries
	m_mesh->newVariable<int,GMDS_EDGE>("connects_boundaries");
	// Singularities indexes
	m_mesh->newVariable<double,GMDS_FACE>("singularity_index");
	// Crosses reference angle
	m_mesh->newVariable<double,GMDS_NODE>("cross_reference_angle");
	// Adjustment angles
	m_mesh->newVariable<double,GMDS_EDGE>("adjustment_angle");
	// Crosses components
	m_mesh->newVariable<math::Vector3d,GMDS_NODE>("U1");
	m_mesh->newVariable<math::Vector3d,GMDS_NODE>("U2");
	m_mesh->newVariable<math::Vector3d,GMDS_NODE>("U3");
	m_mesh->newVariable<math::Vector3d,GMDS_NODE>("U4");
	// Crosses
	m_cross_field_2D = m_mesh->newVariable<math::Cross2D,GMDS_NODE>("cross");
}

/*----------------------------------------------------------------------------*/
CrossField::~CrossField()
{
	// When we define a destructor as below, it doesn't work. Why?????
	//if (m_mesh != nullptr)
		//delete m_mesh;
}

/*----------------------------------------------------------------------------*/
void CrossField::setBoundaryConnectedComponents(int &AN)
{
	std::cout<<"> Identifying the connected components of the boundary"<<std::endl;
	auto IDs = m_mesh->getVariable<int,GMDS_NODE>("boundary_connected_component_id");
	for (auto n_id:m_mesh->nodes())
		IDs->set(n_id,-1);
	int ID = 0;
	for (auto n_id:m_mesh->nodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		if (!isInterior(n) && (IDs->value(n_id) == -1))
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
					for (auto n3:e.get<Node>())
					{
						if (!isInterior(n3) && (IDs->value(n3.id()) == -1))
						{
							IDs->set(n3.id(),ID);
							nodesToAdd.push(n3);
						}
					}
				}
			}
			ID += 1;
		}
	}
	AN = ID;
	std::cout<<"NB connected components : "<<ID<<std::endl;
}

/*----------------------------------------------------------------------------*/
void CrossField::setBoundariesConnexions(std::vector<std::vector<Node>> AV)
{
	m_boundaries_connexions = AV;
}

/*----------------------------------------------------------------------------*/
std::vector<std::vector<Node>> CrossField::connectedBoundaryNodes(std::vector<std::vector<math::Point>> &AV)
{
	std::cout<<"> Connecting the boundaries"<<std::endl;
	std::vector<std::vector<Node>> connectedNodes;
	for (auto points:AV)
	{
		std::vector<Node> nodes;
		TCellID f_id1 = locateInMesh(points[0]);
		Face f1 = m_mesh->get<Face>(f_id1);
		for (auto n:f1.get<Node>())
		{
			if (!isInterior(n))
			{
				nodes.push_back(n);
				break;
			}
		}
		TCellID f_id2 = locateInMesh(points[1]);
		Face f2 = m_mesh->get<Face>(f_id2);
		for (auto n:f2.get<Node>())
		{
			if (!isInterior(n))
			{
				nodes.push_back(n);
				break;
			}
		}
		connectedNodes.push_back(nodes);
	}
	return connectedNodes;
}

/*----------------------------------------------------------------------------*/
std::vector<Edge> CrossField::findPath(Node &AN1, Node &AN2)
{
	// Mark the edges belonging to the path
	auto connectsBoundaries = m_mesh->getVariable<int,GMDS_EDGE>("connects_boundaries");
	// Previous vertex
	auto prev = m_mesh->newVariable<Node,GMDS_NODE>("prev");
	// Edge connecting to the previous vertex
	auto prev_edge = m_mesh->newVariable<Edge,GMDS_NODE>("prev_edge");
	// Mark if the shortest path to a vertex is already found
	auto alreadyFound = m_mesh->newVariable<int,GMDS_NODE>("already_found");
	std::queue<Node> front;
	front.push(AN1);
	alreadyFound->set(AN1.id(),1);
	while(alreadyFound->value(AN2.id()) == 0)
	{
		Node n = front.front();
		front.pop();
		for (auto e:n.get<Edge>())
		{
			if (isInterior(e))
			{
				for (auto n2:e.get<Node>())
				{
					if (alreadyFound->value(n2.id()) == 0)
					{
						alreadyFound->set(n2.id(),1);
						prev->set(n2.id(),n);
						prev_edge->set(n2.id(),e);
						front.push(n2);
					}
				}
			}
		}
	}
	// Find the edges connecting the two nodes
	std::vector<Edge> path;
	Node n = AN2;
	Node nxt;
	while (n.id() != AN1.id())
	{
		nxt = prev->value(n.id());
		path.push_back(prev_edge->value(n.id()));
		connectsBoundaries->set(prev_edge->value(n.id()).id(),1);
		n = nxt;
	}


	m_mesh->deleteVariable(GMDS_NODE,prev);
	m_mesh->deleteVariable(GMDS_NODE,prev_edge);
	m_mesh->deleteVariable(GMDS_NODE,alreadyFound);


	return path;
}

/*----------------------------------------------------------------------------*/
double CrossField::signedArea(const math::Point &AP1, const math::Point &AP2, const math::Point &AP3)
{
	double a = (AP2.X()-AP1.X())*(AP3.Y()-AP1.Y()) - (AP3.X()-AP1.X())*(AP2.Y()-AP1.Y());
	return (1./2.)*a;
}

/*----------------------------------------------------------------------------*/
std::vector<Node> CrossField::orientateEdge(const Edge AE, const gmds::Face AF)
{
	std::vector<Node> orientedEdge;
	Node N1 = AE.get<Node>()[0];
	Node N2 = AE.get<Node>()[1];
	Node N0;
	for (auto n:AF.get<Node>())
	{
		if (n.id() != N1.id() && n.id() != N2.id())
		{
			N0 = n;
			break;
		}
	}
	if (signedArea(N0.point(),N1.point(),N2.point()) >= 0.)
	{
		orientedEdge.push_back(N1);
		orientedEdge.push_back(N2);
	}
	else
	{
		orientedEdge.push_back(N2);
		orientedEdge.push_back(N1);
	}
	return orientedEdge;
}

/*----------------------------------------------------------------------------*/
int CrossField::orientation(const gmds::Edge AE, const gmds::Face AF)
{
	Node N1 = AE.get<Node>()[0];
	Node N2 = AE.get<Node>()[1];
	Node N0;
	for (auto n:AF.get<Node>())
	{
		if (n.id() != N1.id() && n.id() != N2.id())
		{
			N0 = n;
			break;
		}
	}
	if (signedArea(N0.point(),N1.point(),N2.point()) >= 0.)
	{
		return 1;
	}
	else
	{
		return  -1;
	}
}

/*----------------------------------------------------------------------------*/
std::vector<double> CrossField::barycentricCoordinates(const math::Point AP, const gmds::Face AF)
{
	std::vector<double> barCoord;
	for (auto e:AF.get<Edge>())
	{
		std::vector<Node> orientedEdge = orientateEdge(e,AF);
		double bar = signedArea(orientedEdge[0].point(),orientedEdge[1].point(),AP);
		barCoord.push_back(bar);
	}
	return barCoord;
}

/*----------------------------------------------------------------------------*/
TCellID  CrossField::locateInMesh(const math::Point AP)
{
	TCellID Id;
	for (auto f_id:m_mesh->faces())
	{
		Face f = m_mesh->get<Face>(f_id);
		if (barycentricCoordinates(AP,f)[0] >= -0. && barycentricCoordinates(AP,f)[1] >= -0. && barycentricCoordinates(AP,f)[2] >= -0.)
		{
			Id = f_id;
			break;
		}
	}
	return Id;
}

/*----------------------------------------------------------------------------*/
void CrossField::setSingularitiesIndexes(const std::vector<Node> ANodes, const std::vector<double> AIndexes)
{
	auto indexes = m_mesh->getVariable<double,GMDS_FACE>("singularity_index");
	for (auto n_id:m_mesh->nodes())
		indexes->set(n_id,0.);
	for (int i = 0; i<ANodes.size(); i++)
	{
		Node n = ANodes[i];
		// Find the closest vertex to n
		TCellID f_id = locateInMesh(n.point());
		indexes->set(f_id,AIndexes[i]);
	}
}

/*----------------------------------------------------------------------------*/
bool CrossField::isInterior(const gmds::Edge &AE)
{
	return (int(AE.get<Face>().size()) == 2);
}

/*----------------------------------------------------------------------------*/
bool CrossField::isInterior(const gmds::Node &AN)
{
	bool isInt = true;
	for (auto e:AN.get<Edge>())
	{
		if (!isInterior(e))
		{
			isInt = false;
			break;
		}
	}
	return isInt;
}

/*----------------------------------------------------------------------------*/
int CrossField::NbInteriorNodes()
{
	int N = 0;
	for (auto n_id:m_mesh->nodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		if (isInterior(n))
			N += 1;
	}
	std::cout<<"Interior nodes : "<<N<<" total : "<<m_mesh->getNbNodes()<<std::endl;
	return N;
}

/*----------------------------------------------------------------------------*/
int CrossField::NbInteriorEdges()
{
	int N = 0;
	for (auto e_id:m_mesh->edges())
	{
		Edge e = m_mesh->get<Edge>(e_id);
		if (isInterior(e))
			N += 1;
	}
	std::cout<<"Interior edges : "<<N<<" total : "<<m_mesh->getNbEdges()<<std::endl;
	return N;
}

/*----------------------------------------------------------------------------*/
std::complex<double> CrossField::point2complex(math::Point &AP)
{
	std::complex<double> c;
	c.real(AP.X());
	c.imag(AP.Y());
	return c;
}

/*----------------------------------------------------------------------------*/
void CrossField::initializeCrossReferenceAngle()
{
	auto refAngle = m_mesh->getVariable<double,GMDS_NODE>("cross_reference_angle");
	for (auto n_id:m_mesh->nodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		if (!isInterior(n))
		{
			// Get the adjacent boundary edges
			Edge e1, e2;
			for (auto e:n.get<Edge>())
			{
				if (!isInterior(e))
				{
					e1 = e;
					break;
				}
			}
			for (auto e:n.get<Edge>())
			{
				if (!isInterior(e) && (e.id() != e1.id()))
				{
					e2 = e;
					break;
				}
			}
			// Complex representation
			math::Point E1 = e1.get<Node>()[1].point() + (-1.)*e1.get<Node>()[0].point();
			math::Point E2 = e2.get<Node>()[1].point() + (-1.)*e2.get<Node>()[0].point();
			std::complex<double> c1 = point2complex(E1);
			std::complex<double> c2 = point2complex(E2);
			c1 = pow(c1,4);
			c2 = pow(c2,4);
			std::complex<double> c = (c1+c2)/2.;
			c = c/abs(c);
			double alpha = arg(c);
			refAngle->set(n_id,alpha);
		}
	}
}

/*----------------------------------------------------------------------------*/
void CrossField::initializeAdjustmentAngles()
{
	std::cout<<"> Initializing adjustment angles on the boundary"<<std::endl;
	auto adjAngles = m_mesh->getVariable<double,GMDS_EDGE>("adjustment_angle");
	auto refAngle = m_mesh->getVariable<double,GMDS_NODE>("cross_reference_angle");
	double tot = 0.;
	for (auto e_id:m_mesh->edges())
	{
		Edge e = m_mesh->get<Edge>(e_id);
		if (isInterior(e))
			adjAngles->set(e_id,0.);
		else
		{
			Node n1 = e.get<Node>()[0];
			Node n2 = e.get<Node>()[1];
			double omega = refAngle->value(n2.id()) - refAngle->value(n1.id());
			if (omega <= -M_PI)
				omega += 2.*M_PI;
			if (omega > M_PI)
				omega -= 2.*M_PI;
			omega = omega/4.;
			adjAngles->set(e_id,omega);
			Face f = e.get<Face>()[0];
			tot += double(orientation(e,f))*omega;
		}
	}
	// Check if the angles are well initialized
	std::cout<<"Sum of the adjustment angles on the boundary : "<<tot<<std::endl;
}

/*----------------------------------------------------------------------------*/
void CrossField::computeAdjustmentAngles()
{
	// Correspondence cycles/faces and internal edges/edges
	int NbCycles = m_mesh->getNbFaces();
	int NbIntEdges = NbInteriorEdges();
	auto face2cycle = m_mesh->newVariable<int,GMDS_FACE>("face2cycle");
	auto edge2intEdge = m_mesh->newVariable<int,GMDS_EDGE>("edge2intEdge");
	auto singuIndex = m_mesh->getVariable<double,GMDS_FACE>("singularity_index");
	auto adjAngles = m_mesh->getVariable<double,GMDS_EDGE>("adjustment_angle");
	std::vector<TCellID> cycle2face(NbCycles);
	std::vector<TCellID> intEdge2edge(NbIntEdges);
	int id = 0;
	for (auto f_id:m_mesh->faces())
	{
		face2cycle->set(f_id,id);
		cycle2face[id] = f_id;
		id += 1;
	}
	id = 0;
	for (auto e_id:m_mesh->edges())
	{
		Edge e = m_mesh->get<Edge>(e_id);
		if (isInterior(e))
		{
			edge2intEdge->set(e_id,id);
			intEdge2edge[id] = e_id;
			id += 1;
		}
		else
			edge2intEdge->set(e_id,-1);
	}

	// Identify connected components of the boundary
	int NbBoundConnectedComponents;
	setBoundaryConnectedComponents(NbBoundConnectedComponents);
	//NbBoundConnectedComponents = 1;

	// Assemble the system
	// Singularities
	Eigen::VectorXd b(NbCycles + NbBoundConnectedComponents - 1);
	// Set rows corresponding to cycles
	for (int i = 0; i < NbCycles; i++)
	{
		TCellID f_id = cycle2face[i];
		Face f = m_mesh->get<Face>(f_id);
		double sing = singuIndex->value(f_id);
		double constraint = 0.;
		for (auto e:f.get<Edge>())
			constraint += double(orientation(e,f))*adjAngles->value(e.id());
		b[i] = 2.*M_PI*sing - constraint;
	}

	// Cycle basis
	Eigen::SparseMatrix<double> d0(NbCycles + NbBoundConnectedComponents - 1, NbIntEdges);
	// Set the rows corresponding to cycles
	for (auto f_id:m_mesh->faces())
	{
		Face f = m_mesh->get<Face>(f_id);
		int cycleId = face2cycle->value(f_id);
		for (auto e:f.get<Edge>())
		{
			int intEdgeId = edge2intEdge->value(e.id());
			if (intEdgeId >= 0)
				d0.insert(cycleId,intEdgeId) = double(orientation(e,f));
		}
	}

	// Set rows corresponding to paths between connected components of the boundary

	//auto componentsIDs = m_mesh->getVariable<int,GMDS_NODE>("boundary_connected_component_id");
	auto refAngles = m_mesh->getVariable<double,GMDS_NODE>("cross_reference_angle");
	int i = 0;
	for (auto nodes:m_boundaries_connexions)
	{
		Node n0 = nodes[0];
		Node n1 = nodes[1];
		// Update b
		double alpha = refAngles->value(n1.id())/4. - refAngles->value(n0.id())/4.;
		if (alpha <= -M_PI/4.)
			alpha += M_PI/2.;
		if (alpha > M_PI/4.)
			alpha -= M_PI/2.;
		b[NbCycles + i] = alpha;
		// Get the path between the two nodes and update d0
		std::vector<Edge> path = findPath(n0,n1);
		Node a = n1;
		for (auto e:path)
		{
			double orient;
			if (a.id() == e.get<Node>()[0].id())
			{
				orient = 1.;
				a = e.get<Node>()[1];
			}
			else
			{
				orient = -1.;
				a = e.get<Node>()[0];
			}
			d0.insert(NbCycles + i,edge2intEdge->value(e.id())) = orient;
		}
		i += 1;
	}


//	Node n0;
//	for (auto n_id:m_mesh->nodes())
//	{
//		if (componentsIDs->value(n_id) == 0)
//		{
//			n0 = m_mesh->get<Node>(n_id);
//			break;
//		}
//	}
//	// Construct the rows
//	for (int ID = 1; ID < NbBoundConnectedComponents; ID ++)
//	{
//		// Find a node on the current component
//		Node n;
//		for (auto n_id:m_mesh->nodes())
//		{
//			if (componentsIDs->value(n_id) == ID)
//			{
//				n = m_mesh->get<Node>(n_id);
//				break;
//			}
//		}
//		// Update b
//		 double alpha = refAngles->value(n.id())/4. - refAngles->value(n0.id())/4.;
//		 if (alpha <= -M_PI/4.)
//			 alpha += M_PI/2.;
//		 if (alpha > M_PI/4.)
//			 alpha -= M_PI/2.;
//		 // Number of rotations to add on the path: we have to decide!
//		 double k = 0.;
//		 alpha += k*M_PI/2.;
//		 b[NbCycles + ID - 1] = alpha;
//		// Get the path between the two nodes and update d0
//		std::vector<Edge> path = findPath(n0,n);
//		Node a = n;
//		for (auto e:path)
//		{
//			double orient;
//			if (a.id() == e.get<Node>()[0].id())
//			{
//				orient = 1.;
//				a = e.get<Node>()[1];
//			}
//			else
//			{
//				orient = -1.;
//				a = e.get<Node>()[0];
//			}
//			d0.insert(NbCycles + ID - 1,edge2intEdge->value(e.id())) = orient;
//		}
//	}


	// Solving the system
	std::cout<<"> Solving for adjustment angles"<<std::endl;
	Eigen::SparseMatrix<double> A(NbCycles + NbBoundConnectedComponents - 1,NbCycles + NbBoundConnectedComponents - 1);
	A = d0*d0.transpose();
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(A);
	Eigen::VectorXd lambda = chol.solve(b);
	Eigen::VectorXd omega = d0.transpose()*lambda;
	double prob_residual = (A*lambda - b).norm();
	double constraint_residual = (d0*omega-b).norm();
	std::cout<<"> Check if the system is well solved. ||d0*omega-b||="<<constraint_residual<<" and ||A*lambda-b||="<<prob_residual<<std::endl;
	if ((prob_residual < 10e-6) && (constraint_residual < 10e-6))
		std::cout<<"The problem is well solved"<<std::endl;
	else
		std::cout<<"The problem in not well solved"<<std::endl;
	for (int i = 0; i < NbIntEdges; i++)
	{
		TCellID e_id = intEdge2edge[i];
		adjAngles->set(e_id,omega[i]);
	}


	m_mesh->deleteVariable(GMDS_FACE,face2cycle);
	m_mesh->deleteVariable(GMDS_EDGE,edge2intEdge);
}

/*----------------------------------------------------------------------------*/
void CrossField::propagateRefAngleFromBoundary()
{
	// ==========================================================
	// Build a spanning tree of the nodes whith root on the boundary
	// ==========================================================
	// Find the root
	Node root;
	for (auto n_id:m_mesh->nodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		if (!isInterior(n))
		{
			root = n;
			break;
		}
	}
	// Build the tree
	auto sons = m_mesh->newVariable<std::vector<TCellID>,GMDS_NODE>("sons");
	auto alreadyVisited = m_mesh->newVariable<int,GMDS_NODE>("already_visited");
	stack<Node> currentLeaves;
	currentLeaves.push(root);
	alreadyVisited->set(root.id(),1);
	while(!currentLeaves.empty())
	{
		Node n = currentLeaves.top();
		std::vector<TCellID> sonsList;
		currentLeaves.pop();
		for (auto e:n.get<Edge>())
		{
			for (auto n2:e.get<Node>())
			{
				if (alreadyVisited->value(n2.id()) != 1)
				{
					sonsList.push_back(n2.id());
					currentLeaves.push(n2);
					alreadyVisited->set(n2.id(),1);
				}
			}
		}
		sons->set(n.id(),sonsList);
	}


	// ==========================================================
	// Propagate the crosses reference angle from the root
	// ==========================================================
	auto refAngles = m_mesh->getVariable<double,GMDS_NODE>("cross_reference_angle");
	auto adjAngles = m_mesh->getVariable<double,GMDS_EDGE>("adjustment_angle");
	// Put the real ref angle on the root
	refAngles->set(root.id(),refAngles->value(root.id())/4.);
	// Propagate from the root
	stack<Node> positionOnTheGraph;
	positionOnTheGraph.push(root);
	while (!positionOnTheGraph.empty())
	{
		Node n = positionOnTheGraph.top();
		positionOnTheGraph.pop();
		for (auto son:sons->value(n.id()))
		{
			// Find the edge connecting n to its son
			Edge e;
			double orient;
			for (auto e2:n.get<Edge>())
			{
				if ((e2.get<Node>()[0].id() == n.id()) && (e2.get<Node>()[1].id() == son))
				{
					e = e2;
					orient = 1.;
					break;
				}
				if ((e2.get<Node>()[1].id() == n.id()) && (e2.get<Node>()[0].id() == son))
				{
					e = e2;
					orient = -1.;
					break;
				}
			}
			double angle = refAngles->value(n.id()) + orient*adjAngles->value(e.id());
			if (angle <= -M_PI)
				angle += 2.*M_PI;
			if (angle > M_PI)
				angle -= 2.*M_PI;
			refAngles->set(son,angle);
			Node s = m_mesh->get<Node>(son);
			positionOnTheGraph.push(s);
		}
	}



	// Delete the local variables
	m_mesh->deleteVariable(GMDS_NODE,alreadyVisited);
	m_mesh->deleteVariable(GMDS_NODE,sons);
}

/*----------------------------------------------------------------------------*/
void CrossField::buildCrossesFromRefAngles()
{
	auto refAngles = m_mesh->getVariable<double,GMDS_NODE>("cross_reference_angle");
	auto U1 = m_mesh->getVariable<math::Vector3d,GMDS_NODE>("U1");
	auto U2 = m_mesh->getVariable<math::Vector3d,GMDS_NODE>("U2");
	auto U3 = m_mesh->getVariable<math::Vector3d,GMDS_NODE>("U3");
	auto U4 = m_mesh->getVariable<math::Vector3d,GMDS_NODE>("U4");
	for (auto n_id:m_mesh->nodes())
	{
		double alpha = refAngles->value(n_id);
		math::Vector3d u;
		u.setX(cos(alpha));
		u.setY(sin(alpha));
		math::Cross2D c(4.*alpha);
		m_cross_field_2D->set(n_id,c);
		U1->set(n_id,c.componentVectors()[0]);
		U2->set(n_id,c.componentVectors()[1]);
		U3->set(n_id,c.componentVectors()[2]);
		U4->set(n_id,c.componentVectors()[3]);
	}
}

/*----------------------------------------------------------------------------*/
void CrossField::buildCrossField()
{
	std::cout<<" "<<std::endl;
	std::cout<<"==== Building cross field ===="<<std::endl;
	initializeCrossReferenceAngle();
	initializeAdjustmentAngles();
	computeAdjustmentAngles();
	propagateRefAngleFromBoundary();
	buildCrossesFromRefAngles();
	std::cout<<"=============================="<<std::endl;
	std::cout<<" "<<std::endl;
}


/*----------------------------------------------------------------------------*/
}  // end namespace medialaxis
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/