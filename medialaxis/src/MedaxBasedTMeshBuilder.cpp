/*----------------------------------------------------------------------------*/
#include "gmds/medialaxis/MedaxBasedTMeshBuilder.h"
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
MedaxBasedTMeshBuilder::MedaxBasedTMeshBuilder(Mesh &AMedax, Mesh &AMinDel){
	// Medial axis
	m_medax = &AMedax;

	// Corresponding minimal Delaunay triangulation
	m_min_delaunay = &AMinDel;
	// min tri nodes/T-mesh nodes correspondence
	m_min_delaunay->newVariable<int,GMDS_NODE>("minTriNode2TMeshNode");

	// Topological representation
	m_topological_representation = new Mesh(MeshModel(DIM3 | F | E | N |
	                                           F2N | E2N | N2E | N2F));
	// Type of the medial axis section: 0 if it is a quad edge // to the boundary, 1 if it is a quad diagonal, 2 if it is an edge orthogonal to the boundary
	m_topological_representation->newVariable<int,GMDS_EDGE>("section_type");
	// Type of the node: 0 if it is internal (there is an equation at this point), 1 if not (there is a border condition at this point)
	m_topological_representation->newVariable<int,GMDS_NODE>("node_type");
	// ID of the section to which the point belongs
	m_medax->newVariable<int,GMDS_NODE>("section_id");
	// ID of the subsection to which the point belongs
	m_medax->newVariable<int,GMDS_NODE>("subsection_id");
	// Type of the subsection to which the point belongs
	m_medax->newVariable<int,GMDS_NODE>("medial_section_type");
	// Section/section ID correspondence
	m_topological_representation->newVariable<int,GMDS_EDGE>("section_to_section_id");
	// Correspondence medial point/node of the topological representation
	m_medax->newVariable<int,GMDS_NODE>("med_point_to_sing_node");
	// Correspondence medial point/node of the topological representation
	m_topological_representation->newVariable<int,GMDS_NODE>("sing_node_to_med_point");
	// Position of the wings on a section (1 if the wings point in the same direction as the oriented section, -1 if in the opposite direction, 0 if no wing)
	m_topological_representation->newVariable<int,GMDS_EDGE>("wings_position");
	// Position of each section in the vector of degrees of freedom
	m_topological_representation->newVariable<int,GMDS_EDGE>("section2matrix");
	
	// Correspondence degree of freedom/vertex in the dof graph
	m_topological_representation->newVariable<int,GMDS_EDGE>("section2LeftQuadLength");
	m_topological_representation->newVariable<int,GMDS_EDGE>("section2RightQuadLength");
	m_topological_representation->newVariable<int,GMDS_EDGE>("section2QuadHeight");
	m_topological_representation->newVariable<int,GMDS_EDGE>("section2LeftDiagoQuadLength");
	m_topological_representation->newVariable<int,GMDS_EDGE>("section2RightDiagoQuadLength");
	// Topology of a medial axis based block decomposition
	m_t_mesh = new Mesh(MeshModel(DIM3 | F | E | N | F2E |
	                                                  F2N | E2F | E2N | N2E | N2F));
	// Mark with 1 edges of the topological representation that will generate a triangle
	m_topological_representation->newVariable<int,GMDS_EDGE>("generates_a_triangle");
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

	// Correspondance face/id of the medial section
	m_t_mesh->newVariable<int,GMDS_FACE>("face2sectionID");
	// Mark with 1 edges separating different blocks
	m_t_mesh->newVariable<int,GMDS_EDGE>("separates_blocks");
	// Correspondance face/type of its corresponding medial section
	m_t_mesh->newVariable<int,GMDS_FACE>("face2sectionType");
	// ID of the connected component of the section
	m_t_mesh->newVariable<int,GMDS_EDGE>("boundary_connected_component_id");
	// T-mesh node/min tri node correspondance
	m_t_mesh->newVariable<int,GMDS_NODE>("TMeshNode2MinTriNode"); 
	// Mark with 1 T-junctions
	m_t_mesh->newVariable<int,GMDS_NODE>("is_a_T_junction"); 
	// Mark with 1 triangle verticies
	m_t_mesh->newVariable<int,GMDS_NODE>("is_a_triangle_vertex"); 
	// Mark with 1 triangles
	m_t_mesh->newVariable<int,GMDS_FACE>("is_a_triangle");
	// Mark with 1 nodes belonging to an internal constraint
	m_t_mesh->newVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");

	// Final T-mesh
	m_final_t_mesh = new Mesh(MeshModel(DIM3 | F | E | N | F2E |
	                                                  F2N | E2F | E2N | N2E | N2F));
	// Mark with 1 edges corresponding to internal constraints
	m_final_t_mesh->newVariable<int,GMDS_EDGE>("internal_constraint");
	// Mark with 1 nodes belonging to an internal constraint
	m_final_t_mesh->newVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
}

/*----------------------------------------------------------------------------*/
Node
MedaxBasedTMeshBuilder::getNextPoint(const TCellID &APointID, const TCellID &AEdgeID)
{
	// Requires E2N connections
	// To use only if the point APointID belongs to the edge AEdgeID
   Edge e = m_medax-> get<Edge>(AEdgeID);
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
MedaxBasedTMeshBuilder::getNextEdge(const TCellID &AEdgeID, const TCellID &APointID)
{
	// Requires N2E connections
	// To use only if the point APointID belongs to the edge AEdgeID and has type 2
	Node n = m_medax-> get<Node>(APointID);
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
MedaxBasedTMeshBuilder::getExtremPoint(const TCellID &APointID, const TCellID &AEdgeID)
{
	// Requires medial points type set
	// Requires E2E connections
	// To use only if the point APointID belongs to the edge AEdgeID
	auto var = m_medax->getVariable<int,GMDS_NODE>("medial_point_type");
	Node nxtPoint = getNextPoint(APointID, AEdgeID);
	int type = var->value(nxtPoint.id());
	Edge nxtEdge = m_medax->get<Edge>(AEdgeID);
	while (type == 2)
	{
		nxtEdge = getNextEdge(nxtEdge.id(), nxtPoint.id());
		nxtPoint = getNextPoint(nxtPoint.id(), nxtEdge.id());
		type = var->value(nxtPoint.id());
	}
	return nxtPoint;
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::setSectionID()
{
	std::cout<<"> Setting sections IDs"<<std::endl;
	auto medPoint2singNode = m_medax->getVariable<int,GMDS_NODE>("med_point_to_sing_node");
	auto sectionID = m_medax->getVariable<int,GMDS_NODE>("section_id");
	auto subSectionID = m_medax->getVariable<int,GMDS_NODE>("subsection_id");
	auto alreadyVisited = m_medax->newVariable<int,GMDS_NODE>("already_visited");
	auto medPointType = m_medax->getVariable<int,GMDS_NODE>("medial_point_type");
	auto singNode2medPoint = m_topological_representation->getVariable<int,GMDS_NODE>("sing_node_to_med_point");
	auto nodeType = m_topological_representation->getVariable<int,GMDS_NODE>("node_type");
	auto singu = m_medax->getVariable<int,GMDS_NODE>("singularity");
	int ID = 0;
	std::vector<Node> section_start; // Vector storing a starting point for each section
	std::vector<Edge> section_direction; // Vector storing a starting edge containing the starting point for each section
	for (auto n_id:m_medax->nodes())
	{
		if ((medPointType->value(n_id) == 2) && (singu->value(n_id) == 0) && (alreadyVisited->value(n_id) == 0))
		{
			// We build a new section (= new edge of the topo rep)

			// Find the extreme nodes of the section
			Node n = m_medax->get<Node>(n_id);
			alreadyVisited->set(n_id,1);
			sectionID->set(n.id(),ID);
			Node prev, nxt;
			Edge e;
			// We go to the left
			prev = n;
			e = n.get<Edge>()[0];
			nxt = getNextPoint(prev.id(),e.id());
			while ((medPointType->value(nxt.id()) == 2) && (singu->value(nxt.id()) == 0) && alreadyVisited->value(nxt.id()) == 0)
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
			while ((medPointType->value(nxt.id()) == 2) && (singu->value(nxt.id()) == 0) && alreadyVisited->value(nxt.id()) == 0)
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
	
	m_nb_medial_sections = ID;
	std::cout<<"NB sections : "<<ID<<std::endl;
	m_medax->deleteVariable(GMDS_NODE,alreadyVisited);
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::computeSectionType()
{
	auto sectionID = m_medax->getVariable<int,GMDS_NODE>("section_id");
	auto cosMedAngle = m_medax->getVariable<double,GMDS_NODE>("cos_medial_angle");
	auto type = m_medax->getVariable<int,GMDS_NODE>("medial_section_type");
	// Compute the highest section ID
	int maxID = 0;
	for (auto n_id:m_medax->nodes())
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
	for (auto n_id:m_medax->nodes())
	{
		int id = sectionID->value(n_id);
		if (cosMedAngle->value(n_id) < -sqrt(2.)/2.)
			type0[id] = type0[id] + 1;
		else
		{
			type1[id] = type1[id] + 1;
			// if (fabs(cosMedAngle->value(n_id)) < sqrt(2.)/2.)
				
			// else
			// 	type0[id] = type0[id] + 1;
		}
	}
	for (auto n_id:m_medax->nodes())
	{
		int id = sectionID->value(n_id);
		if (type0[id] > type1[id])
			type->set(n_id,0);
		else
			type->set(n_id,1);
	}
}

/*----------------------------------------------------------------------------*/
std::vector<Node> MedaxBasedTMeshBuilder::sectionNodes(int AID)
{
	auto sectionID = m_medax->getVariable<int,GMDS_NODE>("section_id");
	auto medPointType = m_medax->getVariable<int,GMDS_NODE>("medial_point_type");
	auto singu = m_medax->getVariable<int,GMDS_NODE>("singularity");
	auto visited = m_medax->newVariable<int,GMDS_NODE>("visited");
	// Find a node of the section
	Node n0;
	for (auto n_id:m_medax->nodes())
	{
		if ((medPointType->value(n_id) == 2) && (singu->value(n_id) == 0) && (sectionID->value(n_id) == AID))
		{
			n0 = m_medax->get<Node>(n_id);
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
	while ((medPointType->value(nxt.id()) == 2) && (singu->value(nxt.id()) == 0) && visited->value(nxt.id()) == 0)
	{
		if (visited->value(prev.id()) == 0)
			sec.push_back(prev);
		visited->set(prev.id(),1);
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
	while ((medPointType->value(nxt.id()) == 2) && (singu->value(nxt.id()) == 0) && visited->value(nxt.id()) == 0)
	{
		if (visited->value(prev.id()) == 0)
			sec.push_back(prev);
		visited->set(prev.id(),1);
		e = getNextEdge(e.id(),nxt.id());
		prev = nxt;
		nxt = getNextPoint(prev.id(),e.id());
	}
	m_medax->deleteVariable(GMDS_NODE,visited);
	return sec;
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::refineByAddingSingularNodes()
{
	std::cout<<"> Refining the future block decomposition by adding singular nodes"<<std::endl;
	auto sectionID = m_medax->getVariable<int,GMDS_NODE>("section_id");
	auto singu = m_medax->getVariable<int,GMDS_NODE>("singularity");
	// Compute the highest section ID
	int maxID = 0;
	for (auto n_id:m_medax->nodes())
	{
		if (sectionID->value(n_id) > maxID)
			maxID = sectionID->value(n_id);
	}
	std::vector<Node> section;
	for (int id = 0; id <= maxID; id++)
	{
		section = sectionNodes(id);
		int Nb = section.size();
		int q = Nb/30;
		for (int i = 1; i < q; i++)
			singu->set(section[30*i].id(),10);
	}
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::buildTopoRepNodes()
{
	std::cout<<"> Building nodes of the topological representation"<<std::endl;
	std::cout<<"> WARNING: if the medial axis has two neighbouring singular points (for example a singu and an IP), the topo rep is not valid"<<std::endl;
	auto singu = m_medax->getVariable<int,GMDS_NODE>("singularity");
	auto medPointType = m_medax->getVariable<int,GMDS_NODE>("medial_point_type");
	auto medPoint2singNode = m_medax->getVariable<int,GMDS_NODE>("med_point_to_sing_node");
	auto singNode2medPoint = m_topological_representation->getVariable<int,GMDS_NODE>("sing_node_to_med_point");
	auto nodeType = m_topological_representation->getVariable<int,GMDS_NODE>("node_type");
	auto touchingPoints = m_medax->getVariable<std::vector<math::Point>,GMDS_NODE>("touching_points");

	for (auto n_id:m_medax->nodes())
	{
		if (medPointType->value(n_id) == 2)
		{
			if (singu->value(n_id) != 0)
			{
				Node n = m_medax->get<Node>(n_id);
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
			if (medPointType->value(n_id) == 1)
			{
				Node n = m_medax->get<Node>(n_id);
				Node newNode = m_topological_representation->newNode(touchingPoints->value(n_id)[0]);
				medPoint2singNode->set(n_id,newNode.id());
				singNode2medPoint->set(newNode.id(),n_id);
				nodeType->set(newNode.id(),medPointType->value(n_id));
			}
			else
			{
				Node n = m_medax->get<Node>(n_id);
				Node newNode = m_topological_representation->newNode(n.point());
				medPoint2singNode->set(n.id(),newNode.id());
				singNode2medPoint->set(newNode.id(),n.id());
				nodeType->set(newNode.id(),medPointType->value(n_id));
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::buildTopoRepEdges()
{
	std::cout<<"> Building edges of the topological representation"<<std::endl;
	std::cout<<"> WARNING: if the medial axis has only one section which is a cycle, the topo rep is not valid"<<std::endl;
	auto medPointType = m_medax->getVariable<int,GMDS_NODE>("medial_point_type");
	auto cosMedAngle = m_medax->getVariable<double,GMDS_NODE>("cos_medial_angle");
	auto sectionType = m_topological_representation->getVariable<int,GMDS_EDGE>("section_type");
	auto secTypeOnNodes = m_medax->getVariable<int,GMDS_NODE>("medial_section_type");
	auto medPoint2singNode = m_medax->getVariable<int,GMDS_NODE>("med_point_to_sing_node");
	auto wings = m_topological_representation->getVariable<int,GMDS_EDGE>("wings_position");
	auto touchingPoints = m_medax->getVariable<std::vector<math::Point>,GMDS_NODE>("touching_points");
	auto section2sectionID = m_topological_representation->getVariable<int,GMDS_EDGE>("section_to_section_id");
	auto sectionID = m_medax->getVariable<int,GMDS_NODE>("section_id");
	//auto subSectionID = m_medax->getVariable<int,GMDS_NODE>("subsection_id");
	auto tailDir = m_topological_representation->getVariable<math::Vector,GMDS_EDGE>("tailDirection");
	auto headDir = m_topological_representation->getVariable<math::Vector,GMDS_EDGE>("headDirection");
	for (auto n_id:m_medax->nodes())
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
			Node n = m_medax->get<Node>(n_id);
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
			{
				if (type != 1)
				{
					type = 1;
					sectID = m_nb_medial_sections;
					m_nb_medial_sections = m_nb_medial_sections + 1;
				}
			}
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
MedaxBasedTMeshBuilder::writeTopoRep(std::basic_string<char> AFileName)
{
	std::cout<<"> Writing the topological representation"<<std::endl;
	IGMeshIOService ioService(m_topological_representation);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(AFileName);
}

/*----------------------------------------------------------------------------*/
int MedaxBasedTMeshBuilder::wings(gmds::Node &AN, gmds::Edge &AE)
{
	auto wings = m_topological_representation->getVariable<int,GMDS_EDGE>("wings_position");
	if (orientation(AN,AE) == wings->value(AE.id()))
		return fabs(orientation(AN,AE));
	else
		return 0;
}

/*----------------------------------------------------------------------------*/
int MedaxBasedTMeshBuilder::orientation(gmds::Node &AN, gmds::Edge &AE)
{
	if (AN.id() == AE.get<Node>()[0].id())
		return 1;
	if (AN.id() == AE.get<Node>()[1].id())
		return -1;
	return 0;
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::setTopoRepConnectivity()
{
	std::cout<<"> Setting topological representation connectivity"<<std::endl;
	MeshDoctor doc(m_topological_representation);
	doc.updateUpwardConnectivity();
}

/*----------------------------------------------------------------------------*/
Node MedaxBasedTMeshBuilder::getNextSingularNode(gmds::Node &AN, gmds::Edge &AE)
{
	if (AE.get<Node>()[0].id() == AN.id())
		return AE.get<Node>()[1];
	if (AE.get<Node>()[1].id() == AN.id())
		return AE.get<Node>()[0];
	std::cout<<"getNextSingularNode : error, the given Node does not belong to the given edge"<<std::endl;
	return AN;
}

/*----------------------------------------------------------------------------*/
Edge MedaxBasedTMeshBuilder::getNextMedialSection(gmds::Edge &AE, gmds::Node &AN)
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
void MedaxBasedTMeshBuilder::browseTopoRep()
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
std::vector<Node> MedaxBasedTMeshBuilder::shortestPathAlongBoundaryOrConstraints(Node &AN1, Node &AN2)
{
	auto constr = m_min_delaunay->getVariable<int,GMDS_EDGE>("internal_constraint");
	auto visited = m_min_delaunay->newVariable<int,GMDS_NODE>("visited");
	std::vector<Node> path;
	std::vector<Node> new_path;
	std::vector<Node> shortestPath;
	path.push_back(AN1);
	visited->set(AN1.id(),1);
	std::queue<std::vector<Node>> toBeContinued;
	toBeContinued.push(path);
	bool finished = false;
	while (!toBeContinued.empty() && !finished)
	{
		path = toBeContinued.front();
		toBeContinued.pop();
		Node last_node = path[path.size()-1];
		if (last_node.id() == AN2.id())
		{
			shortestPath = path;
			finished = true;
		}
		else
		{
			for (auto e:last_node.get<Edge>())
			{
				if (e.get<Face>().size() == 1 || constr->value(e.id()) == 1)
				{
					for (auto n:e.get<Node>())
					{
						if (visited->value(n.id()) == 0)
						{
							new_path = path;
							new_path.push_back(n);
							visited->set(n.id(),1);
							toBeContinued.push(new_path);
						}
					}
				}
			}
		}
	}
	m_min_delaunay->deleteVariable(GMDS_NODE,visited);
	return shortestPath;
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::markEdgesGeneratingTriangles()
{
	auto gt = m_topological_representation->getVariable<int,GMDS_EDGE>("generates_a_triangle");
	auto boundaryTangentNodes = m_medax->getVariable<std::vector<Node>,GMDS_NODE>("tangent_boundary_nodes");
	auto singNode2medPoint = m_topological_representation->getVariable<int,GMDS_NODE>("sing_node_to_med_point");
	for (auto e_id:m_topological_representation->edges())
	{
		Edge e = m_topological_representation->get<Edge>(e_id);
		Node n1 = e.get<Node>()[0];
		Node n2 = e.get<Node>()[1];
		std::vector<Node> tbn1 = boundaryTangentNodes->value(singNode2medPoint->value(n1.id()));
		std::vector<Node> tbn2 = boundaryTangentNodes->value(singNode2medPoint->value(n2.id()));
		bool tri = false;
		for (auto n:tbn1)
		{
			for (auto m:tbn2)
			{
				if (n.id() == m.id())
					tri = true;
			}
		}
		if (tri)
			gt->set(e_id,1);
	}
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::buildTMeshNodesFromMinDelNodes()
{
	std::cout<<"> Building medial and boundary nodes of the medial axis based T-mesh"<<std::endl;
	auto singNode2medPoint = m_topological_representation->getVariable<int,GMDS_NODE>("sing_node_to_med_point");
	auto boundaryTangentNodes = m_medax->getVariable<std::vector<Node>,GMDS_NODE>("tangent_boundary_nodes");
	auto minTriNode2TMeshNode = m_min_delaunay->getVariable<int,GMDS_NODE>("minTriNode2TMeshNode");
	auto TMeshNode2MinTriNode = m_t_mesh->getVariable<int,GMDS_NODE>("TMeshNode2MinTriNode");
	auto sn2mn = m_topological_representation->getVariable<int,GMDS_NODE>("singuNode2MedNode");
	auto nodeType = m_topological_representation->getVariable<int,GMDS_NODE>("node_type");
	auto sn2bn = m_topological_representation->getVariable<std::vector<int>,GMDS_NODE>("singuNode2BoundNodes");
	auto del_nodes_constr = m_min_delaunay->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
	auto tmesh_nodes_constr = m_t_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
	// Mark nodes of the minimal triangulation appearing in the T-mesh
	auto appears = m_min_delaunay->newMark<Node>();
	for (auto n_id:m_topological_representation->nodes())
	{
		
		int med_point = singNode2medPoint->value(n_id);
		std::vector<Node> tangent_nodes = boundaryTangentNodes->value(med_point);
		for (auto n:tangent_nodes)
			m_min_delaunay->mark<Node>(n.id(),appears);
		
	}
	// Add these nodes to the T-mesh
	for (auto n_id:m_min_delaunay->nodes())
	{
		if (m_min_delaunay->isMarked<Node>(n_id,appears))
		{
			Node n = m_min_delaunay->get<Node>(n_id);
			Node new_node = m_t_mesh->newNode(n.point());
			minTriNode2TMeshNode->set(n_id,new_node.id());
			TMeshNode2MinTriNode->set(new_node.id(),n_id);
			if (del_nodes_constr->value(n.id()) == 1)
				tmesh_nodes_constr->set(new_node.id(),1);
		}
		else
		{
			minTriNode2TMeshNode->set(n_id,-1);
		}
	}
	// Add medial nodes to the T-mesh
	for (auto n_id:m_topological_representation->nodes())
	{
		if (nodeType->value(n_id) >= 2)
		{
			// Build medial node
			Node singu_node = m_topological_representation->get<Node>(n_id);
			Node newNode;
			newNode = m_t_mesh->newNode(singu_node.point());
			sn2mn->set(n_id,newNode.id());
			// Update sn2bn
			std::vector<int> bound_nodes;
			int med_point = singNode2medPoint->value(n_id);
			std::vector<Node> tangent_nodes = boundaryTangentNodes->value(med_point);
			for (auto n:tangent_nodes)
			{
				bound_nodes.push_back(minTriNode2TMeshNode->value(n.id()));
			}
			sn2bn->set(n_id,bound_nodes);
		}
		if (nodeType->value(n_id) == 1)
		{
			sn2mn->set(n_id,minTriNode2TMeshNode->value(boundaryTangentNodes->value(singNode2medPoint->value(n_id))[0].id()));
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::buildBlockDecompMedialAndBoundaryNodes()
{
	std::cout<<"> Building medial and boundary nodes of the medial axis based block decomposition"<<std::endl;
	auto nodeType = m_topological_representation->getVariable<int,GMDS_NODE>("node_type");
	auto singNode2medPoint = m_topological_representation->getVariable<int,GMDS_NODE>("sing_node_to_med_point");
	auto touchingPoints = m_medax->getVariable<std::vector<math::Point>,GMDS_NODE>("touching_points");
	auto boundaryTangentNodes = m_medax->getVariable<std::vector<Node>,GMDS_NODE>("tangent_boundary_nodes");
	auto sn2mn = m_topological_representation->getVariable<int,GMDS_NODE>("singuNode2MedNode");
	auto sn2bn = m_topological_representation->getVariable<std::vector<int>,GMDS_NODE>("singuNode2BoundNodes");
	auto dualTriangles = m_medax->getVariable<std::vector<Face>,GMDS_NODE>("dual_triangles");
	for (auto n_id:m_topological_representation->nodes())
	{
		// // Build medial node
		// Node singu_node = m_topological_representation->get<Node>(n_id);
		// Node newNode;
		// newNode = m_t_mesh->newNode(singu_node.point());
		// sn2mn->set(n_id,newNode.id());
		// if (nodeType->value(n_id) == 1)
		// 	continue;
		// // Build boundary nodes
		// Node med_point = m_medax->get<Node>(singNode2medPoint->value(n_id));
		// std::vector<int> newNodes;
		// std::vector<math::Point> boundaryTangencyPoints = touchingPoints->value(med_point.id());

		// if (nodeType->value(n_id) >= 2)
		// {
		// 	for (auto p:boundaryTangencyPoints)
		// 	{
		// 		newNode = m_t_mesh->newNode(p);
		// 		newNodes.push_back(newNode.id());
		// 	}
		// }
		// sn2bn->set(n_id,newNodes);

		// Build medial node
		Node singu_node = m_topological_representation->get<Node>(n_id);
		Node newNode;
		newNode = m_t_mesh->newNode(singu_node.point());
		sn2mn->set(n_id,newNode.id());
		if (nodeType->value(n_id) == 1)
			continue;
		// Build boundary nodes
		Node med_point = m_medax->get<Node>(singNode2medPoint->value(n_id));
		std::vector<int> newNodes;
		std::vector<Node> tangentNodes = boundaryTangentNodes->value(med_point.id());

		if (nodeType->value(n_id) >= 2)
		{
			for (auto p:tangentNodes)
			{
				newNode = m_t_mesh->newNode(p.point());
				newNodes.push_back(newNode.id());
			}
		}
		sn2bn->set(n_id,newNodes);
	}
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::writeBlockDecomp(std::basic_string<char> AFileName)
{
	std::cout<<"> Writing the medial axis based block decomposition"<<std::endl;
	IGMeshIOService ioService(m_t_mesh);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(AFileName);
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::buildSection2MedialAndBoundaryNodesAdjacency()
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
						Node boundary_node = m_t_mesh->get<Node>(ind);
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
						Node boundary_node = m_t_mesh->get<Node>(ind);
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
void MedaxBasedTMeshBuilder::buildMiddleNodes()
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
					Node ln = m_t_mesh->get<Node>(LTBN->value(s_id));
					Node rn = m_t_mesh->get<Node>(RTBN->value(s_id));
					Node newNode;
					if (LTMN->value(s_id) < 0)
					{
						// Placement of the new node
						math::Point dir = ln.point()+(-1.)*n.point();
						dir = (1./vec(dir).norm())*dir;
						double l = (1./sqrt(2.))*section.length();
						if (l > vec((ln.point()+n.point())*(1./2.)+(-1.)*n.point()).norm())
							l = vec((ln.point()+n.point())*(1./2.)+(-1.)*n.point()).norm();
						//newNode = m_t_mesh->newNode((ln.point()+n.point())*(1./2.));
						newNode = m_t_mesh->newNode(n.point()+l*dir);
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
						//newNode = m_t_mesh->newNode((rn.point()+n.point())*(1./2.));
						newNode = m_t_mesh->newNode(n.point()+l*dir);
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
					Node ln = m_t_mesh->get<Node>(LHBN->value(s_id));
					Node rn = m_t_mesh->get<Node>(RHBN->value(s_id));
					Node newNode;
					if (LHMN->value(s_id) < 0)
					{
						// Placement of the new node
						math::Point dir = ln.point()+(-1.)*n.point();
						dir = (1./vec(dir).norm())*dir;
						double l = (1./sqrt(2.))*section.length();
						if (l > vec((ln.point()+n.point())*(1./2.)+(-1.)*n.point()).norm())
							l = vec((ln.point()+n.point())*(1./2.)+(-1.)*n.point()).norm();
						//newNode = m_t_mesh->newNode((ln.point()+n.point())*(1./2.));
						newNode = m_t_mesh->newNode(n.point()+l*dir);
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
						//newNode = m_t_mesh->newNode((rn.point()+n.point())*(1./2.));
						newNode = m_t_mesh->newNode(n.point()+l*dir);
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
void MedaxBasedTMeshBuilder::buildBlocks()
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
	auto f2sid = m_t_mesh->getVariable<int,GMDS_FACE>("face2sectionID");
	auto f2st = m_t_mesh->getVariable<int,GMDS_FACE>("face2sectionType");
	auto section2sectionID = m_topological_representation->getVariable<int,GMDS_EDGE>("section_to_section_id");
	auto t_junct = m_t_mesh->getVariable<int,GMDS_NODE>("is_a_T_junction"); 
	auto tri_ver = m_t_mesh->getVariable<int,GMDS_NODE>("is_a_triangle_vertex"); 
	auto tri = m_t_mesh->getVariable<int,GMDS_FACE>("is_a_triangle"); 
	auto minTriNode2TMeshNode = m_min_delaunay->getVariable<int,GMDS_NODE>("minTriNode2TMeshNode");
	auto TMeshNode2MinTriNode = m_t_mesh->getVariable<int,GMDS_NODE>("TMeshNode2MinTriNode");

	int tn, hn, ltbn, rtbn, lhbn, rhbn, lhmn, rhmn, ltmn, rtmn;
	Face new_face;
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
			{
				block1.push_back(lhmn);
				t_junct->set(lhmn,1);
			}
			block1.push_back(lhbn);
			Node n1 = m_min_delaunay->get<Node>(TMeshNode2MinTriNode->value(lhbn));
			Node n2 = m_min_delaunay->get<Node>(TMeshNode2MinTriNode->value(ltbn));
			std::vector<Node> inbetweenBoundNodes = shortestPathAlongBoundaryOrConstraints(n1,n2);
			if (inbetweenBoundNodes.size() >= 3)
			{
				for (int i = 1; i < inbetweenBoundNodes.size()-1; i++)
				{
					if (minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()) >= 0)
					{
						block1.push_back(minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()));
						t_junct->set(minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()),1);
					}
				}
			}
			block1.push_back(ltbn);
			if (ltmn >= 0)
			{
				block1.push_back(ltmn);
				t_junct->set(ltmn,1);
			}
			block1.push_back(tn);
			new_face = m_t_mesh->newFace(block1);
			f2sid->set(new_face.id(),section2sectionID->value(s_id));
			f2st->set(new_face.id(),0);
			bool isATri = (lhbn == ltbn);
			if (isATri)
			{
				tri->set(new_face.id(),1);
				tri_ver->set(lhbn,1);
			}

			std::vector<TCellID> block2;
			block2.push_back(tn);
			if (rtmn >= 0)
			{
				block2.push_back(rtmn);
				t_junct->set(rtmn,1);
			}
			block2.push_back(rtbn);
			n1 = m_min_delaunay->get<Node>(TMeshNode2MinTriNode->value(rtbn));
			n2 = m_min_delaunay->get<Node>(TMeshNode2MinTriNode->value(rhbn));
			inbetweenBoundNodes = shortestPathAlongBoundaryOrConstraints(n1,n2);
			if (inbetweenBoundNodes.size() >= 3)
			{
				for (int i = 1; i < inbetweenBoundNodes.size()-1; i++)
				{
					if (minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()) >= 0)
					{
						block2.push_back(minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()));
						t_junct->set(minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()),1);
					}
				}
			}
			block2.push_back(rhbn);
			if (rhmn >= 0)
			{
				block2.push_back(rhmn);
				t_junct->set(rhmn,1);
			}
			block2.push_back(hn);
			new_face = m_t_mesh->newFace(block2);
			f2sid->set(new_face.id(),section2sectionID->value(s_id));
			f2st->set(new_face.id(),0);
			isATri = (rtbn == rhbn);
			if (isATri)
			{
				tri->set(new_face.id(),1);
				tri_ver->set(rtbn,1);
			}
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
				{
					block.push_back(lhmn);
					t_junct->set(lhmn,1);
				}
				if (lhbn == -1)
				{
					Node n1 = m_min_delaunay->get<Node>(TMeshNode2MinTriNode->value(hn));
					Node n2 = m_min_delaunay->get<Node>(TMeshNode2MinTriNode->value(ltbn));
					std::vector<Node> inbetweenBoundNodes = shortestPathAlongBoundaryOrConstraints(n1,n2);
					if (inbetweenBoundNodes.size() >= 3)
					{
						for (int i = 1; i < inbetweenBoundNodes.size()-1; i++)
						{
							if (minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()) >= 0)
							{
								block.push_back(minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()));
								t_junct->set(minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()),1);
							}
						}
					}
				}
				block.push_back(lhbn+ltbn+1);
				if (ltmn >= 0)
				{
					block.push_back(ltmn);
					t_junct->set(ltmn,1);
				}
				if (ltbn == -1)
				{
					Node n1 = m_min_delaunay->get<Node>(TMeshNode2MinTriNode->value(lhbn));
					Node n2 = m_min_delaunay->get<Node>(TMeshNode2MinTriNode->value(tn));
					std::vector<Node> inbetweenBoundNodes = shortestPathAlongBoundaryOrConstraints(n1,n2);
					if (inbetweenBoundNodes.size() >= 3)
					{
						for (int i = 1; i < inbetweenBoundNodes.size()-1; i++)
						{
							if (minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()) >= 0)
							{
								block.push_back(minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()));
								t_junct->set(minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()),1);
							}
						}
					}
				}
				block.push_back(tn);
				if (rtmn >= 0)
				{
					block.push_back(rtmn);
					t_junct->set(rtmn,1);
				}
				if (rtbn == -1)
				{
					Node n1 = m_min_delaunay->get<Node>(TMeshNode2MinTriNode->value(tn));
					Node n2 = m_min_delaunay->get<Node>(TMeshNode2MinTriNode->value(rhbn));
					std::vector<Node> inbetweenBoundNodes = shortestPathAlongBoundaryOrConstraints(n1,n2);
					if (inbetweenBoundNodes.size() >= 3)
					{
						for (int i = 1; i < inbetweenBoundNodes.size()-1; i++)
						{
							if (minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()) >= 0)
							{
								block.push_back(minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()));
								t_junct->set(minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()),1);
							}
						}
					}
				}
				block.push_back(rhbn+rtbn+1);
				if (rhmn >= 0)
				{
					block.push_back(rhmn);
					t_junct->set(rhmn,1);
				}
				if (rhbn == -1)
				{
					Node n1 = m_min_delaunay->get<Node>(TMeshNode2MinTriNode->value(rtbn));
					Node n2 = m_min_delaunay->get<Node>(TMeshNode2MinTriNode->value(hn));
					std::vector<Node> inbetweenBoundNodes = shortestPathAlongBoundaryOrConstraints(n1,n2);
					if (inbetweenBoundNodes.size() >= 3)
					{
						for (int i = 1; i < inbetweenBoundNodes.size()-1; i++)
						{
							if (minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()) >= 0)
							{
								block.push_back(minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()));
								t_junct->set(minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()),1);
							}
						}
					}
				}
				new_face = m_t_mesh->newFace(block);
				f2sid->set(new_face.id(),section2sectionID->value(s_id));
				f2st->set(new_face.id(),1);
			}

			else
			{
				if (wings(section.get<Node>()[0],section) == 1)
				{
					new_face = m_t_mesh->newQuad(hn,lhmn,tn,rhmn);
					f2sid->set(new_face.id(),section2sectionID->value(s_id));
					f2st->set(new_face.id(),1);
					std::vector<TCellID> block;
					block.push_back(tn);
					block.push_back(lhmn);
					block.push_back(lhbn);
					Node n1 = m_min_delaunay->get<Node>(TMeshNode2MinTriNode->value(lhbn));
					Node n2 = m_min_delaunay->get<Node>(TMeshNode2MinTriNode->value(ltbn));
					std::vector<Node> inbetweenBoundNodes = shortestPathAlongBoundaryOrConstraints(n1,n2);
					if (inbetweenBoundNodes.size() >= 3)
					{
						for (int i = 1; i < inbetweenBoundNodes.size()-1; i++)
						{
							if (minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()) >= 0)
							{
								block.push_back(minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()));
								t_junct->set(minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()),1);
							}
						}
					}
					block.push_back(ltbn);
					if (ltmn >= 0)
					{
						block.push_back(ltmn);
						t_junct->set(ltmn,1);
					}
					new_face = m_t_mesh->newFace(block);
					f2sid->set(new_face.id(),section2sectionID->value(s_id));
					f2st->set(new_face.id(),1);
					bool isATri = (lhbn == ltbn);
					if (isATri)
					{
						tri->set(new_face.id(),1);
						tri_ver->set(lhbn,1);
					}
					block.clear();
					block.push_back(tn);
					if (rtmn >= 0)
					{
						block.push_back(rtmn);
						t_junct->set(rtmn,1);
					}
					block.push_back(rtbn);
					n1 = m_min_delaunay->get<Node>(TMeshNode2MinTriNode->value(rtbn));
					n2 = m_min_delaunay->get<Node>(TMeshNode2MinTriNode->value(rhbn));
					inbetweenBoundNodes = shortestPathAlongBoundaryOrConstraints(n1,n2);
					if (inbetweenBoundNodes.size() >= 3)
					{
						for (int i = 1; i < inbetweenBoundNodes.size()-1; i++)
						{
							if (minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()) >= 0)
							{
								block.push_back(minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()));
								t_junct->set(minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()),1);
							}
						}
					}
					block.push_back(rhbn);
					block.push_back(rhmn);
					new_face = m_t_mesh->newFace(block);
					f2sid->set(new_face.id(),section2sectionID->value(s_id));
					f2st->set(new_face.id(),1);
					isATri = (rtbn == rhbn);
					if (isATri)
					{
						tri->set(new_face.id(),1);
						tri_ver->set(rtbn,1);
					}
				}
				else
				{
					new_face = m_t_mesh->newQuad(hn,ltmn,tn,rtmn);
					f2sid->set(new_face.id(),section2sectionID->value(s_id));
					f2st->set(new_face.id(),1);
					std::vector<TCellID> block;
					block.push_back(hn);
					if (lhmn >= 0)
					{
						block.push_back(lhmn);
						t_junct->set(lhmn,1);
					}
					block.push_back(lhbn);
					Node n1 = m_min_delaunay->get<Node>(TMeshNode2MinTriNode->value(lhbn));
					Node n2 = m_min_delaunay->get<Node>(TMeshNode2MinTriNode->value(ltbn));
					std::vector<Node> inbetweenBoundNodes = shortestPathAlongBoundaryOrConstraints(n1,n2);
					if (inbetweenBoundNodes.size() >= 3)
					{
						for (int i = 1; i < inbetweenBoundNodes.size()-1; i++)
						{
							if (minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()) >= 0)
							{
								block.push_back(minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()));
								t_junct->set(minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()),1);
							}
						}
					}
					block.push_back(ltbn);
					block.push_back(ltmn);
					new_face = m_t_mesh->newFace(block);
					f2sid->set(new_face.id(),section2sectionID->value(s_id));
					f2st->set(new_face.id(),1);
					bool isATri = (lhbn == ltbn);
					if (isATri)
					{
						tri->set(new_face.id(),1);
						tri_ver->set(lhbn,1);
					}
					block.clear();
					block.push_back(hn);
					block.push_back(rtmn);
					block.push_back(rtbn);
					n1 = m_min_delaunay->get<Node>(TMeshNode2MinTriNode->value(rtbn));
					n2 = m_min_delaunay->get<Node>(TMeshNode2MinTriNode->value(rhbn));
					inbetweenBoundNodes = shortestPathAlongBoundaryOrConstraints(n1,n2);
					if (inbetweenBoundNodes.size() >= 3)
					{
						for (int i = 1; i < inbetweenBoundNodes.size()-1; i++)
						{
							if (minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()) >= 0)
							{
								block.push_back(minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()));
								t_junct->set(minTriNode2TMeshNode->value(inbetweenBoundNodes[i].id()),1);
							}
						}
					}
					block.push_back(rhbn);
					if (rhmn >= 0)
					{
						block.push_back(rhmn);
						t_junct->set(rhmn,1);
					}
					new_face = m_t_mesh->newFace(block);
					f2sid->set(new_face.id(),section2sectionID->value(s_id));
					f2st->set(new_face.id(),1);
					isATri = (rtbn == rhbn);
					if (isATri)
					{
						tri->set(new_face.id(),1);
						tri_ver->set(rtbn,1);
					}
					// m_t_mesh->newQuad(hn,lhbn,ltbn,ltmn);
					// m_t_mesh->newQuad(hn,rtmn,rtbn,rhbn);
				}
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::transformDegenerateQuadsIntoTriangles()
{
	std::cout<<"> Transforming degenerate quads into triangles"<<std::endl;
	auto tri = m_t_mesh->getVariable<int,GMDS_FACE>("is_a_triangle"); 
	auto tri_ver = m_t_mesh->getVariable<int,GMDS_NODE>("is_a_triangle_vertex"); 
	for (auto f_id:m_t_mesh->faces())
	{
		if (tri->value(f_id) == 1)
		{
			Face f = m_t_mesh->get<Face>(f_id);
			std::vector<TCellID> nodes;
			bool double_vertex_added = false;
			for (auto n:f.get<Node>())
			{
				if (tri_ver->value(n.id()) == 0)
				{
					nodes.push_back(n.id());
				}
				else
				{
					if (!double_vertex_added)
					{
						nodes.push_back(n.id());
						double_vertex_added = true;
					}
				}
			}
			Face new_face = m_t_mesh->newFace(nodes);
			tri->set(new_face.id(),2);
			m_t_mesh->deleteFace(f_id);
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::transformTrianglesIntoQuads()
{
	std::cout<<"> Transforming triangles into non degenerate quads"<<std::endl;
	auto tri = m_t_mesh->getVariable<int,GMDS_FACE>("is_a_triangle"); 
	auto tri_ver = m_t_mesh->getVariable<int,GMDS_NODE>("is_a_triangle_vertex");
	auto tri_edge = m_t_mesh->newVariable<int,GMDS_EDGE>("belongs_to_triangle");
	auto tmesh_nodes_constr = m_t_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
	// Mark with 1 edges belonging to a quad
	for (auto e_id:m_t_mesh->edges())
	{
		Edge e = m_t_mesh->get<Edge>(e_id);
		for (auto f:e.get<Face>())
		{
			if (tri->value(f.id()) == 2)
			{
				tri_edge->set(e_id,1);
				break;
			}
		}
	}
	// For each triangle vertex, transform the triangles of its fan into quads
	for (auto n_id:m_t_mesh->nodes())
	{
		if (tri_ver->value(n_id) == 1)
		{
			Node n = m_t_mesh->get<Node>(n_id);
			int NbNonTriEdges = 0;
			for (auto e:n.get<Edge>())
			{
				if (tri_edge->value(e.id()) == 0)
				{
					NbNonTriEdges += 1;
				}
			}
			if (NbNonTriEdges == 1)
			{
				std::vector<Edge> adj_edges = n.get<Edge>();
				adj_edges = sortEdges(n,adj_edges);
				std::vector<Edge> tri_edges;
				int regular_edge_pos;
				Edge regular_edge;
				for (int i = 0; i < adj_edges.size(); i++)
				{
					if (tri_edge->value(adj_edges[i].id()) == 0)
					{
						regular_edge_pos = i;
						regular_edge = adj_edges[i];
						break;
					}
				}
				for (int i = 1; i < adj_edges.size(); i++)
				{
					tri_edges.push_back(adj_edges[(regular_edge_pos+i)%(adj_edges.size())]);
				}
				if (tri_edges.size() == 2)
					throw GMDSException("transformTrianglesIntoQuads() : a constraint triangle fan has only one triangle");
				// Create the new nodes and quads
				Node prev, nxt;
				Face new_face;
				for (int i = 0; i < tri_edges.size()-1; i++)
				{
					if (i == 0)
					{
						prev = n;
					}
					else
					{
						prev = nxt;
					}
					if (i < tri_edges.size()-2)
					{
						nxt = m_t_mesh->newNode(n.point());
						if (tmesh_nodes_constr->value(n.id()) == 1)
							tmesh_nodes_constr->set(nxt.id(),1);
					}
					else
					{
						nxt = n;
					}
					Face f = getCommonFace(tri_edges[i],tri_edges[i+1]);
					std::vector<TCellID> nodes;
					for (auto n1:f.get<Node>())
					{
						if (n1.id() != n_id)
							nodes.push_back(n1.id());
						else
						{
							nodes.push_back(nxt.id());
							nodes.push_back(prev.id());
						}
					}
					new_face = m_t_mesh->newFace(nodes);
					m_t_mesh->deleteFace(f.id());
				}
			}
			if (NbNonTriEdges == 2)
			{
				std::vector<Edge> adj_edges = n.get<Edge>();
				adj_edges = sortEdges(n,adj_edges);
				std::vector<Edge> tri_edges1;
				std::vector<Edge> tri_edges2;
				int regular_edge_pos1,regular_edge_pos2;
				Edge regular_edge1,regular_edge2;
				for (int i = 0; i < adj_edges.size(); i++)
				{
					if (tri_edge->value(adj_edges[i].id()) == 0)
					{
						regular_edge_pos1 = i;
						regular_edge1 = adj_edges[i];
						break;
					}
				}
				for (int i = 0; i < adj_edges.size(); i++)
				{
					if (tri_edge->value(adj_edges[i].id()) == 0 && adj_edges[i].id() != regular_edge1.id())
					{
						regular_edge_pos2 = i;
						regular_edge2 = adj_edges[i];
						break;
					}
				}
				// Build the two fans of edges
				for (int i = 1; i < regular_edge_pos2 - regular_edge_pos1; i++)
				{
					tri_edges1.push_back(adj_edges[(regular_edge_pos1+i)%adj_edges.size()]);
				}
				for (int i = 1; i < regular_edge_pos1 - regular_edge_pos2 + adj_edges.size(); i++)
				{
					tri_edges2.push_back(adj_edges[(regular_edge_pos2+i)%adj_edges.size()]);
				}
				if (tri_edges1.empty() || tri_edges2.empty())
				{
					std::vector<Edge> tri_edges;
					Edge e1,e2;
					if (tri_edges1.empty())
					{
						tri_edges = tri_edges2;
						e1 = regular_edge2;
						e2 = regular_edge1;
					}
					else
					{
						tri_edges = tri_edges1;
						e1 = regular_edge1;
						e2 = regular_edge2;
					}
					int NbNewPoints = tri_edges.size()-1;
					math::Vector dir = edge2vec(e2,n);
					std::vector<Node> new_points;
					for (int i = 0; i < NbNewPoints; i++)
					{
						Node new_node = m_t_mesh->newNode(n.point()+(double(i+1)/double(NbNewPoints+1))*dir);
						new_points.push_back(new_node);
						if (tmesh_nodes_constr->value(n.id()) == 1)
							tmesh_nodes_constr->set(new_node.id(),1);
					}
					for (int i = 0; i < tri_edges.size()-1; i++)
					{
						Face f = getCommonFace(tri_edges[i],tri_edges[i+1]);
						std::vector<TCellID> nodes;
						Node prev,nxt;
						if (i == 0)
							prev = n;
						else
							prev = new_points[i-1];
						nxt = new_points[i];
						for (auto n1:f.get<Node>())
						{
							if (n1.id() != n_id)
								nodes.push_back(n1.id());
							else
							{
								nodes.push_back(nxt.id());
								nodes.push_back(prev.id());
							}
						}
						Face new_face = m_t_mesh->newFace(nodes);
						m_t_mesh->deleteFace(f.id());
					}
					// Change the neigbouring face
					Face f = getCommonFace(tri_edges[tri_edges.size()-1],e2);
					std::vector<TCellID> nodes;
					for (auto n1:f.get<Node>())
					{
						if (n1.id() != n_id)
							nodes.push_back(n1.id());
						else
						{
							nodes.push_back(new_points[new_points.size()-1].id());
						}
					}
					Face new_face = m_t_mesh->newFace(nodes);
					m_t_mesh->deleteFace(f.id());
					// Change the face behind
					if (isInterior(n))
					{
						f = getCommonFace(e1,e2);
						nodes.clear();
						for (auto n1:f.get<Node>())
						{
							if (n1.id() != n_id)
								nodes.push_back(n1.id());
							else
							{
								nodes.push_back(n.id());
								for (int i = 0; i < new_points.size(); i++)
									nodes.push_back(new_points[i].id());
							}
						}
						new_face = m_t_mesh->newFace(nodes);
						m_t_mesh->deleteFace(f.id());
					}
				}
				else
				{
					// If both fans are not empty
					// Deal with the first fan
					int NbNewPoints = tri_edges1.size()-1;
					math::Vector dir = edge2vec(regular_edge2,n);
					std::vector<Node> new_points;
					for (int i = 0; i < NbNewPoints; i++)
					{
						Node new_node = m_t_mesh->newNode(n.point()+(double(i+1)/double(NbNewPoints+1))*dir);
						new_points.push_back(new_node);
						if (tmesh_nodes_constr->value(n.id()) == 1)
							tmesh_nodes_constr->set(new_node.id(),1);
					}
					for (int i = 0; i < tri_edges1.size()-1; i++)
					{
						Face f = getCommonFace(tri_edges1[i],tri_edges1[i+1]);
						std::vector<TCellID> nodes;
						Node prev,nxt;
						if (i == 0)
							prev = n;
						else
							prev = new_points[i-1];
						nxt = new_points[i];
						for (auto n1:f.get<Node>())
						{
							if (n1.id() != n_id)
								nodes.push_back(n1.id());
							else
							{
								nodes.push_back(nxt.id());
								nodes.push_back(prev.id());
							}
						}
						Face new_face = m_t_mesh->newFace(nodes);
						m_t_mesh->deleteFace(f.id());
					}
					// Change the neigbouring face
					Face f = getCommonFace(tri_edges1[tri_edges1.size()-1],regular_edge2);
					std::vector<TCellID> nodes;
					for (auto n1:f.get<Node>())
					{
						if (n1.id() != n_id)
							nodes.push_back(n1.id());
						else
						{
							nodes.push_back(new_points[new_points.size()-1].id());
						}
					}
					Face new_face = m_t_mesh->newFace(nodes);
					m_t_mesh->deleteFace(f.id());
					// Change the face behind
					if (isInterior(n))
					{
						f = getCommonFace(tri_edges2[0],regular_edge2);
						nodes.clear();
						for (auto n1:f.get<Node>())
						{
							if (n1.id() != n_id)
								nodes.push_back(n1.id());
							else
							{
								nodes.push_back(n.id());
								for (int i = 0; i < new_points.size(); i++)
									nodes.push_back(new_points[i].id());
							}
						}
						new_face = m_t_mesh->newFace(nodes);
						m_t_mesh->deleteFace(f.id());
					}
					// Deal with the second fan
					NbNewPoints = tri_edges2.size()-1;
					dir = edge2vec(regular_edge1,n);
					new_points.clear();
					for (int i = 0; i < NbNewPoints; i++)
					{
						Node new_node = m_t_mesh->newNode(n.point()+(double(i+1)/double(NbNewPoints+1))*dir);
						new_points.push_back(new_node);
						if (tmesh_nodes_constr->value(n.id()) == 1)
							tmesh_nodes_constr->set(new_node.id(),1);
					}
					for (int i = 0; i < tri_edges2.size()-1; i++)
					{
						Face f = getCommonFace(tri_edges2[i],tri_edges2[i+1]);
						std::vector<TCellID> nodes;
						Node prev,nxt;
						if (i == 0)
							prev = n;
						else
							prev = new_points[i-1];
						nxt = new_points[i];
						for (auto n1:f.get<Node>())
						{
							if (n1.id() != n_id)
								nodes.push_back(n1.id());
							else
							{
								nodes.push_back(nxt.id());
								nodes.push_back(prev.id());
							}
						}
						Face new_face = m_t_mesh->newFace(nodes);
						m_t_mesh->deleteFace(f.id());
					}
					// Change the neigbouring face
					f = getCommonFace(tri_edges2[tri_edges2.size()-1],regular_edge1);
					nodes.clear();
					for (auto n1:f.get<Node>())
					{
						if (n1.id() != n_id)
							nodes.push_back(n1.id());
						else
						{
							nodes.push_back(new_points[new_points.size()-1].id());
						}
					}
					new_face = m_t_mesh->newFace(nodes);
					m_t_mesh->deleteFace(f.id());
					// Change the face behind
					if (isInterior(n))
					{
						f = getCommonFace(tri_edges1[0],regular_edge1);
						nodes.clear();
						for (auto n1:f.get<Node>())
						{
							if (n1.id() != n_id)
								nodes.push_back(n1.id());
							else
							{
								nodes.push_back(n.id());
								for (int i = 0; i < new_points.size(); i++)
									nodes.push_back(new_points[i].id());
							}
						}
						new_face = m_t_mesh->newFace(nodes);
						m_t_mesh->deleteFace(f.id());
					}
				}
			}
			if (NbNonTriEdges == 3)
			{
				std::vector<Edge> adj_edges = n.get<Edge>();
				adj_edges = sortEdges(n,adj_edges);
				std::vector<Edge> tri_edges;
				int regular_edge_pos1,regular_edge_pos2,regular_edge_pos3;
				Edge regular_edge1,regular_edge2,regular_edge3;
				if (tri_edge->value(adj_edges[0].id()) == 0 && tri_edge->value(adj_edges[1].id()) == 1)
				{
					regular_edge1 = adj_edges[adj_edges.size()-2];
					regular_edge_pos1 = adj_edges.size()-2;
					regular_edge2 = adj_edges[adj_edges.size()-1];
					regular_edge_pos2 = adj_edges.size()-1;
					regular_edge3 = adj_edges[0];
					regular_edge_pos3 = 0;
				}
				else if (tri_edge->value(adj_edges[adj_edges.size()-1].id()) == 0 && tri_edge->value(adj_edges[adj_edges.size()-2].id()) == 1)
				{
					regular_edge1 = adj_edges[adj_edges.size()-1];
					regular_edge_pos1 = adj_edges.size()-1;
					regular_edge2 = adj_edges[0];
					regular_edge_pos2 = 0;
					regular_edge3 = adj_edges[1];
					regular_edge_pos3 = 1;
				}
				else
				{
					for (int i = 0; i < adj_edges.size(); i++)
					{
						if (tri_edge->value(adj_edges[i].id()) == 0)
						{
							regular_edge_pos1 = i;
							regular_edge1 = adj_edges[i];
							break;
						}
					}
					for (int i = 0; i < adj_edges.size(); i++)
					{
						if (tri_edge->value(adj_edges[i].id()) == 0 && adj_edges[i].id() != regular_edge1.id())
						{
							regular_edge_pos2 = i;
							regular_edge2 = adj_edges[i];
							break;
						}
					}
					for (int i = 0; i < adj_edges.size(); i++)
					{
						if (tri_edge->value(adj_edges[i].id()) == 0 && adj_edges[i].id() != regular_edge1.id() && adj_edges[i].id() != regular_edge2.id())
						{
							regular_edge_pos3 = i;
							regular_edge3 = adj_edges[i];
							break;
						}
					}
				}
				// Build the fan of edges
				if (regular_edge_pos3 > regular_edge_pos1)
				{
					for (int i = 1; i < adj_edges.size() + regular_edge_pos1 - regular_edge_pos3; i++)
					{
						tri_edges.push_back(adj_edges[(regular_edge_pos3+i)%adj_edges.size()]);
					}
				}
				else
				{
					for (int i = 1; i < regular_edge_pos1 - regular_edge_pos3; i++)
					{
						tri_edges.push_back(adj_edges[(regular_edge_pos3+i)%adj_edges.size()]);
					}
				}
				// Add the new verticies and transform the triangles into quads
				int NbNewPoints = tri_edges.size()-1;
				math::Vector dir = edge2vec(regular_edge1,n);
				std::vector<Node> new_points;
				for (int i = 0; i < NbNewPoints; i++)
				{
					Node new_node = m_t_mesh->newNode(n.point()+(double(i+1)/double(NbNewPoints+1))*dir);
					new_points.push_back(new_node);
					if (tmesh_nodes_constr->value(n.id()) == 1)
						tmesh_nodes_constr->set(new_node.id(),1);
				}
				for (int i = 0; i < tri_edges.size()-1; i++)
				{
					Face f = getCommonFace(tri_edges[i],tri_edges[i+1]);
					std::vector<TCellID> nodes;
					Node prev,nxt;
					if (i == 0)
						prev = n;
					else
						prev = new_points[i-1];
					nxt = new_points[i];
					for (auto n1:f.get<Node>())
					{
						if (n1.id() != n_id)
							nodes.push_back(n1.id());
						else
						{
							nodes.push_back(nxt.id());
							nodes.push_back(prev.id());
						}
					}
					Face new_face = m_t_mesh->newFace(nodes);
					m_t_mesh->deleteFace(f.id());
				}
				// Change the neigbouring face
				Face f = getCommonFace(tri_edges[tri_edges.size()-1],regular_edge1);
				std::vector<TCellID> nodes;
				for (auto n1:f.get<Node>())
				{
					if (n1.id() != n_id)
						nodes.push_back(n1.id());
					else
					{
						nodes.push_back(new_points[new_points.size()-1].id());
					}
				}
				Face new_face = m_t_mesh->newFace(nodes);
				m_t_mesh->deleteFace(f.id());
				// Change the face behind
				if (isInterior(n))
				{
					f = getCommonFace(regular_edge1,regular_edge2);
					nodes.clear();
					for (auto n1:f.get<Node>())
					{
						if (n1.id() != n_id)
							nodes.push_back(n1.id());
						else
						{
							nodes.push_back(n.id());
							for (int i = 0; i < new_points.size(); i++)
								nodes.push_back(new_points[i].id());
						}
					}
					new_face = m_t_mesh->newFace(nodes);
					m_t_mesh->deleteFace(f.id());
				}
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
std::vector<Edge> MedaxBasedTMeshBuilder::sortedAdjacentEdges(gmds::Node &AN)
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
std::vector<Edge> MedaxBasedTMeshBuilder::neighbouringEdges(gmds::Edge &AE, gmds::Node &AN)
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
std::vector<Edge> MedaxBasedTMeshBuilder::orderedNeigbourSections(Edge &ASection, Node &AN)
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
void MedaxBasedTMeshBuilder::setBlockDecompConnectivity()
{
	std::cout<<"> Setting block decomposition connectivity"<<std::endl;
	MeshDoctor doc(m_t_mesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::markBlocksSeparatingEdges()
{
	std::cout<<"> Marking separating edges"<<std::endl;
	auto f2sID = m_t_mesh->getVariable<int,GMDS_FACE>("face2sectionID");
	auto f2st = m_t_mesh->getVariable<int,GMDS_FACE>("face2sectionType");
	auto sb = m_t_mesh->getVariable<int,GMDS_EDGE>("separates_blocks");
	for (auto e_id:m_t_mesh->edges())
	{
		Edge e = m_t_mesh->get<Edge>(e_id);
		if (e.get<Face>().size() == 2)
		{
			Face f1 = e.get<Face>()[0];
			Face f2 = e.get<Face>()[1];
			if (f2sID->value(f1.id()) != f2sID->value(f2.id()))
				sb->set(e_id,1);
			else
			{
				if (f2st->value(f1.id()) == 0 && !touchesBoundary(e))
					sb->set(e_id,1);
			}
		}
	}
	// Unmark the superfluous marks
	for (auto e_id:m_t_mesh->edges())
	{
		Edge e = m_t_mesh->get<Edge>(e_id);
		if (e.get<Face>().size() == 2 && sb->value(e_id) == 1)
		{
			Face f1 = e.get<Face>()[0];
			Face f2 = e.get<Face>()[1];
			if (f2sID->value(f1.id()) != f2sID->value(f2.id()) && f2st->value(f1.id()) == 1 && f2st->value(f2.id()) == 1)
			{
				bool isACorner1, isACorner2;
				int NbBoundaryEdges = 0;
				for (auto e1:f1.get<Edge>())
				{
					if (e1.get<Face>().size() == 1)
						NbBoundaryEdges += 1;
				}
				isACorner1 = (NbBoundaryEdges == 2);
				NbBoundaryEdges = 0;
				for (auto e1:f2.get<Edge>())
				{
					if (e1.get<Face>().size() == 1)
						NbBoundaryEdges += 1;
				}
				isACorner2 = (NbBoundaryEdges == 2);
				if (isACorner1 || isACorner2)
					sb->set(e_id,0);
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::addBigTJunctions()
{
	std::cout<<"> Adding big T-junctions"<<std::endl;
	auto sb = m_t_mesh->getVariable<int,GMDS_EDGE>("separates_blocks");
	auto f2sID = m_t_mesh->getVariable<int,GMDS_FACE>("face2sectionID");
	for (auto n_id:m_t_mesh->nodes())
	{
		Node n = m_t_mesh->get<Node>(n_id);
		if (isInterior(n) && n.get<Face>().size() == 5)
		{
			for (auto e:n.get<Edge>())
			{
				if (sb->value(e.id()) == 0)
				{
					int sectID = f2sID->value(e.get<Face>()[0].id());
					Edge target;
					double min_angle = 15.;
					for (auto e_id2:m_t_mesh->edges())
					{
						Edge e2 = m_t_mesh->get<Edge>(e_id2);
						if (e2.get<Face>().size() == 2)
						{
							if ((f2sID->value(e2.get<Face>()[0].id()) == sectID || f2sID->value(e2.get<Face>()[1].id()) == sectID) && f2sID->value(e2.get<Face>()[0].id()) != f2sID->value(e2.get<Face>()[1].id()))
							{
								math::Point C = (1./2.)*(e2.get<Node>()[0].point()+e2.get<Node>()[1].point());
								double alpha = oriented_angle(edge2vec(e,n),vec(C+(-1.)*n.point()));
								if (fabs(alpha) < min_angle)
								{
									min_angle = fabs(alpha);
									target = e2;
								}
							}
						}
					}
					math::Point bigTJunction;
					if ((n.point().distance(target.get<Node>()[0].point()) < n.point().distance(target.get<Node>()[1].point()) && isInterior(target.get<Node>()[0])) || !isInterior(target.get<Node>()[1]))
						bigTJunction = target.get<Node>()[0].point();
					else
						bigTJunction = target.get<Node>()[1].point();
					if (n.point().distance((1./2.)*(target.get<Node>()[0].point()+target.get<Node>()[1].point())) < n.point().distance(bigTJunction))	
						bigTJunction = (1./2.)*(target.get<Node>()[0].point()+target.get<Node>()[1].point());
					Node new_node = m_t_mesh->newNode(bigTJunction);
					Edge new_edge = m_t_mesh->newEdge(n.id(),new_node.id());
					sb->set(new_edge.id(),1);
				}
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::ensureConnectivityThroughConstraints()
{
	std::cout<<"> Ensuring connectivity through constraints"<<std::endl;
	// create set of non zero boundary edges an of boundary nodes
	auto visited = m_t_mesh->newMark<Node>();
	std::vector<Edge> boundary_edges;
	std::vector<Node> boundary_nodes;
	for (auto e_id:m_t_mesh->edges())
	{
		Edge e = m_t_mesh->get<Edge>(e_id);
		if (e.get<Face>().size() == 1 && e.length() > 1e-6)
		{
			boundary_edges.push_back(e);
			for (auto n:e.get<Node>())
			{
				if (!m_t_mesh->isMarked<Node>(n.id(),visited))
				{
					boundary_nodes.push_back(n);
					m_t_mesh->mark<Node>(n.id(),visited);
				}
			}
		}
	}
	// Mofify faces on constraints
	for (auto e:boundary_edges)
	{
		Face f = e.get<Face>()[0];
		std::vector<Node> adj_nodes = f.get<Node>();
		std::vector<TCellID> ids;
		for (auto n:adj_nodes)
			ids.push_back(n.id());
		for (auto n:boundary_nodes)
		{
			if (n.id()!=e.get<Node>()[0].id() && n.id()!=e.get<Node>()[1].id() && isOnSegment(n.point(),e.get<Node>()[0].point(),e.get<Node>()[1].point()))
			{
				ids = insertPoint(n,adj_nodes);
				adj_nodes.clear();
				for (auto id:ids)
					adj_nodes.push_back(m_t_mesh->get<Node>(id));
			}
		}
		m_t_mesh->newFace(ids);
		m_t_mesh->deleteFace(f.id());
		m_t_mesh->deleteEdge(e.id());
	}
	Face f = m_t_mesh->get<Face>(265);
	std::cout<<"hfiehrt "<<f.get<Node>().size()<<std::endl;
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::setBoundaryConnectedComponents()
{
	std::cout<<"> Identifying the connected components of the boundary of the T-mesh"<<std::endl;
	std::vector<std::vector<Edge>> components;
	auto IDs = m_t_mesh->getVariable<int,GMDS_EDGE>("boundary_connected_component_id");
	for (auto e_id:m_t_mesh->edges())
		IDs->set(e_id,-1);
	int ID = 0;
	for (auto e_id:m_t_mesh->edges())
	{
		Edge e = m_t_mesh->get<Edge>(e_id);
		if (e.get<Face>().size() == 1 && IDs->value(e_id) == -1)
		{
			// Then we build a new connected component
			std::vector<Edge> component;
			std::stack<Edge> toAdd;
			toAdd.push(e);
			IDs->set(e_id,ID);
			while (!toAdd.empty())
			{
				e = toAdd.top();
				toAdd.pop();
				component.push_back(e);
				for (auto n:e.get<Node>())
				{
					for (auto e1:n.get<Edge>())
					{
						if (e1.get<Face>().size() == 1 && IDs->value(e1.id()) == -1)
						{
							toAdd.push(e1);
							IDs->set(e1.id(),ID);
						}
					}
				}
			}
			components.push_back(component);
			ID += 1;
		}
	}
	std::cout<<"NB boundary connected components : "<<ID<<std::endl;
	m_boundary_connected_components = components;
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::buildFinalTMesh()
{
	std::cout<<"> Building the final T-mesh"<<std::endl;
	auto tmesh_nodes_constr = m_t_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
	auto final_tmesh_nodes_constr = m_final_t_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
	for (auto n_id:m_t_mesh->nodes())
	{
		Node n = m_t_mesh->get<Node>(n_id);
		Node new_node = m_final_t_mesh->newNode(n.point());
		if (tmesh_nodes_constr->value(n_id) == 1)
			final_tmesh_nodes_constr->set(new_node.id(),1);
	}
	for (auto f_id:m_t_mesh->faces())
	{
		Face f = m_t_mesh->get<Face>(f_id);
		std::vector<TCellID> nodes;
		for (auto n:f.get<Node>())
		{
			nodes.push_back(n.id());
		}
		Face new_face = m_final_t_mesh->newFace(nodes);
	}
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::writeFinalTMesh(std::basic_string<char> AFileName)
{
	std::cout<<"> Writing the final T-mesh"<<std::endl;
	IGMeshIOService ioService(m_final_t_mesh);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(AFileName);
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::setFinalTMeshConnectivity()
{
	std::cout<<"> Setting connectivity of the final T-mesh"<<std::endl;
	MeshDoctor doc(m_final_t_mesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
}

/*----------------------------------------------------------------------------*/
Mesh MedaxBasedTMeshBuilder::getFinalTMesh()
{
	return *m_final_t_mesh;
}

/*----------------------------------------------------------------------------*/
void MedaxBasedTMeshBuilder::markInternalConstraintsOnFinalTMesh()
{
	std::cout<<"> Marking internal constraints on the final T-mesh"<<std::endl;
	auto nodes_constr = m_final_t_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
	auto int_constr = m_final_t_mesh->getVariable<int,GMDS_EDGE>("internal_constraint");
	//auto tmn2mtn = m_t_mesh->getVariable<int,GMDS_NODE>("TMeshNode2MinTriNode"); 
	for (auto e_id:m_final_t_mesh->edges())
	{
		Edge e = m_final_t_mesh->get<Edge>(e_id);
		Node n1 = e.get<Node>()[0];
		Node n2 = e.get<Node>()[1];
		if (nodes_constr->value(n1.id()) == 1 && nodes_constr->value(n2.id()) == 1)
			int_constr->set(e.id(),1);
		else if (nodes_constr->value(n1.id()) == 1 && isInterior(n1) && !isInterior(n2))
			int_constr->set(e.id(),1);	
		else if (nodes_constr->value(n2.id()) == 1 && isInterior(n2) && !isInterior(n1))
			int_constr->set(e.id(),1);	
	}
}
/*----------------------------------------------------------------------------*/
}  // end namespace medialaxis
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/