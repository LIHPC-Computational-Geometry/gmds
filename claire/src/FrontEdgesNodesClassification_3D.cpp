//
// Created by rochec on 20/12/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/claire/FrontEdgesNodesClassification_3D.h>
#include <gmds/claire/NodeNeighbourhoodOnFront_3D.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <iostream>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/
FrontEdgesNodesClassification_3D::FrontEdgesNodesClassification_3D(Mesh *AMesh, Front_3D *AFront, Variable<int>* A_EdgesClassification, Variable<int>* A_NodesClassification,
                                                                   Mesh *AMesh_T, FastLocalize *Afl, Variable<math::Vector3d>* A_VectorField) :
  m_mesh(AMesh),
  m_mesh_T(AMesh_T),
  m_fl(Afl),
  m_Front(AFront)
{
	m_EdgesClassification = A_EdgesClassification;
	m_NodesClassification = A_NodesClassification;
	m_mark_EdgesForTemplates = m_mesh->newMark<Edge>();
	m_mark_NodesForTemplates = m_mesh->newMark<Node>();
	m_VectorField = A_VectorField;
}
/*------------------------------------------------------------------------*/
FrontEdgesNodesClassification_3D::~FrontEdgesNodesClassification_3D()
{
	m_mesh->unmarkAll<Edge>(m_mark_EdgesForTemplates);
	m_mesh->freeMark<Edge>(m_mark_EdgesForTemplates);
	m_mesh->unmarkAll<Node>(m_mark_NodesForTemplates);
	m_mesh->freeMark<Node>(m_mark_NodesForTemplates);
}
/*------------------------------------------------------------------------*/
FrontEdgesNodesClassification_3D::STATUS
FrontEdgesNodesClassification_3D::execute()
{
	FrontEdgesClassification();	// Fill the variable m_EdgesClassification

	gmds::IGMeshIOService ioService(m_mesh);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::E);
	vtkWriter.setDataOptions(gmds::N|gmds::E);
	vtkWriter.write("AeroEdgesClassification_3D_"+std::to_string(m_Front->getFrontID())+".vtk");

	FrontNodesClassification();	// Fill the variable m_NodesClassification

	ComputeValid_GFE();

	// Reclassify all the edges and nodes
	//FrontEdgesClassification();
	//FrontNodesClassification();
	ComputeValidLoop_GFE();

	return FrontEdgesNodesClassification_3D::SUCCESS;
}
/*------------------------------------------------------------------------*/
TInt
FrontEdgesNodesClassification_3D::getMarkEdgesTemplates()
{
	return m_mark_EdgesForTemplates;
}
/*------------------------------------------------------------------------*/
TInt
FrontEdgesNodesClassification_3D::getMarkNodesTemplates()
{
	return m_mark_NodesForTemplates;
}
/*------------------------------------------------------------------------*/
int
FrontEdgesNodesClassification_3D::SingleEdgeClassification(TCellID e_id)
{
	int edge_classification(0);

	Variable<int>* var_node_couche_id = m_mesh->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_face_couche_id = m_mesh->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");

	Edge e = m_mesh->get<Edge>(e_id);

	std::vector<Face> e_faces = e.get<Face>();
	std::vector<Face> e_faces_on_front;
	for (auto const& f:e_faces)
	{
		if (var_face_couche_id->value(f.id()) == m_Front->getFrontID())
		{
			e_faces_on_front.push_back(f);
		}
	}
	// Each edge of the front is connected to two faces of the front.
	// So the size of the vector e_faces_on_front is 2.

	math::Vector3d n1 = e_faces_on_front[0].normal();
	math::Vector3d n2 = e_faces_on_front[1].normal();

	math::Point f1_center = e_faces_on_front[0].center();
	math::Point f2_center = e_faces_on_front[1].center();

	math::Vector3d v1;
	math::Vector3d v2;

	if (m_Front->getFrontID() > 0) {
		// Each face of the front is connected to a unique region
		// So the size of the vector e_faces_on_front[0].get<Region>() is 1.
		Region r1 = (e_faces_on_front[0].get<Region>())[0];
		Region r2 = (e_faces_on_front[1].get<Region>())[0];

		math::Point r1_center = r1.center();
		math::Point r2_center = r2.center();

		v1 = (f1_center - r1.center()).normalize();
		v2 = (f2_center - r2.center()).normalize();
	}
	else
	{	// If we are on the first front, there is no region in the mesh.
		gmds::Cell::Data data_1 = m_fl->find(e_faces_on_front[0].center());
		gmds::Cell::Data data_2 = m_fl->find(e_faces_on_front[1].center());

		v1 = m_VectorField->value(data_1.id);
		v2 = m_VectorField->value(data_2.id);
	}

	// Compute if the vectors n1 and n2 are well oriented.
	if ( n1.dot(v1) < 0 )
	{
		n1=-n1;
	}
	if ( n2.dot(v2) < 0)
	{
		n2=-n2;
	}

	math::Vector3d v_f1_tan = (f1_center-e.center()).normalize() ;
	math::Vector3d v_f2_tan = (f2_center-e.center()).normalize() ;

	/*
	double angle(0);
	std::vector<Region> e_regions = e.get<Region>() ;
	for (auto r:e_regions)
	{
	   std::vector<Face> two_adj_faces = math::Utils::getFacesAdjToEdgeInHexa(m_meshH, e.id(), r.id()) ;
	   Edge e_opp_0 = math::Utils::oppositeEdgeInFace(m_meshH, e.id(), two_adj_faces[0].id());
	   Edge e_opp_1 = math::Utils::oppositeEdgeInFace(m_meshH, e.id(), two_adj_faces[1].id());

	   math::Point center_e = e.center();
	   math::Vector3d v0 = (e_opp_0.center()-center_e).normalize() ;
	   math::Vector3d v1 = (e_opp_1.center()-center_e).normalize() ;

	   angle += acos(v0.dot(v1));
	}
	 */

	// Edge classification.
	// 0 : Side
	// 1 : Corner
	// 2 : End
	// 3 : Reversal

	double angle = acos(n1.dot(n2));

	if (M_PI/4.0 <= angle && angle < 3.0*M_PI/4.0)
	{

		//edge_classification = 1;
		if ( v_f1_tan.dot(n2) <= 0 && v_f2_tan.dot(n1) <= 0)
		{
			edge_classification = 1;
		}
		else
		{
			/*
			std::cout << "--------------------------" << std::endl;
			std::cout << "Edge: " << e.get<Node>()[0].id() << ", " << e.get<Node>()[1].id() << std::endl;
			std::cout << v_f1_tan.dot(n2) << std::endl;
			std::cout << v_f2_tan.dot(n1) << std::endl;
			 */
			edge_classification = 2;
		}

	}
	else if (5.0*M_PI/4.0 <= angle && angle < 7.0*M_PI/4.0)
	{
		edge_classification = 2;
	}
	else if (3.0*M_PI/4.0 <= angle && angle < 5.0*M_PI/4.0)
	{
		edge_classification = 3;
	}

	return edge_classification;
}
/*------------------------------------------------------------------------*/
void
FrontEdgesNodesClassification_3D::FrontEdgesClassification()
{
	Variable<int>* var_node_couche_id = m_mesh->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_face_couche_id = m_mesh->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");

	for (auto e_id:m_mesh->edges())
	{
		Edge e = m_mesh->get<Edge>(e_id);
		std::vector<Node> e_nodes = e.get<Node>();
		if (var_node_couche_id->value(e_nodes[0].id()) == m_Front->getFrontID()
		    && var_node_couche_id->value(e_nodes[1].id()) == m_Front->getFrontID())
		{
			// So, the edge is on the front.
			int edge_classification = SingleEdgeClassification(e_id);
			m_EdgesClassification->set(e_id, edge_classification);
		}
	}

}
/*------------------------------------------------------------------------*/
void
FrontEdgesNodesClassification_3D::FrontNodesClassification()
{
	Variable<int>* var_node_couche_id = m_mesh->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");

	for (auto n_id:m_Front->getNodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		std::vector<Edge> n_edges = n.get<Edge>();
		int compteur_local_feature_edges(0);
		for (auto const& e:n_edges)
		{
			std::vector<Node> e_nodes = e.get<Node>();
			if (var_node_couche_id->value(e_nodes[0].id()) == m_Front->getFrontID()
			    && var_node_couche_id->value(e_nodes[1].id()) == m_Front->getFrontID()
			    && m_EdgesClassification->value(e.id()) >= 1)
			{
				compteur_local_feature_edges++;
			}
		}

		m_NodesClassification->set(n_id, compteur_local_feature_edges);
	}

}
/*------------------------------------------------------------------------*/
int
FrontEdgesNodesClassification_3D::singleNodeClassification(TCellID n_id){
	int classification(0);

	NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n_id);
	n_neighbourhood.execute();

	// Compute the number or corner, end and reversal edges around the node n
	// on the front
	int compteur_corner(0);
	int compteur_end(0);
	int compteur_reversal(0);
	for (auto e_id: n_neighbourhood.getOrderedEdges())
	{
		Edge e = m_mesh->get<Edge>(e_id);
		if (m_EdgesClassification->value(e_id) == 1)
		{
			compteur_corner +=1;
		}
		else if (m_EdgesClassification->value(e_id) == 2)
		{
			compteur_end +=1;
		}
		else if (m_EdgesClassification->value(e_id) == 3)
		{
			compteur_reversal +=1;
		}
	}

	// Check the ordering of the singular edges around the node n
	if (compteur_corner == 3 && compteur_end == 0 && compteur_reversal == 0)
	{
		classification = 1;
	}
	else if (compteur_corner==2 && compteur_end==1)
	{
		classification = 2;
	}
	else if (compteur_corner==1 && compteur_end==2)
	{
		classification = 3;
	}
	else if (compteur_end == 3 && n_neighbourhood.getOrderedEdges().size()==3)
	{
		classification = 4;
	}
	else if (n_neighbourhood.getOrderedEdges().size() == 6
	         && ( ( m_EdgesClassification->value(n_neighbourhood.getOrderedEdges()[0]) == 1
	              && m_EdgesClassification->value(n_neighbourhood.getOrderedEdges()[2]) == 1
	              && m_EdgesClassification->value(n_neighbourhood.getOrderedEdges()[4]) == 1
	              && m_EdgesClassification->value(n_neighbourhood.getOrderedEdges()[1]) == 2
	              && m_EdgesClassification->value(n_neighbourhood.getOrderedEdges()[3]) == 2
	              && m_EdgesClassification->value(n_neighbourhood.getOrderedEdges()[5]) == 2 )
	             || ( m_EdgesClassification->value(n_neighbourhood.getOrderedEdges()[1]) == 1
	                 && m_EdgesClassification->value(n_neighbourhood.getOrderedEdges()[3]) == 1
	                 && m_EdgesClassification->value(n_neighbourhood.getOrderedEdges()[5]) == 1
	                 && m_EdgesClassification->value(n_neighbourhood.getOrderedEdges()[0]) == 2
	                 && m_EdgesClassification->value(n_neighbourhood.getOrderedEdges()[2]) == 2
	                 && m_EdgesClassification->value(n_neighbourhood.getOrderedEdges()[4]) == 2 ) ) )
	{
		classification = 5;
	}
	else if (n_neighbourhood.getOrderedEdges().size() == 6
	         && compteur_corner == 2 && compteur_end==2 && compteur_reversal==0)
	{
		classification = 6;
	}
	else if (n_neighbourhood.getOrderedEdges().size() == 3
	         && compteur_corner == 2 && compteur_end==0 && compteur_reversal==1)
	{
		classification = 7;
	}
	else if (n_neighbourhood.getOrderedEdges().size() == 6
	         && compteur_corner == 0 && compteur_end==2 && compteur_reversal==1)
	{
		classification = 8;
	}

	return classification;
}
/*------------------------------------------------------------------------*/
bool
FrontEdgesNodesClassification_3D::isValidNodeForTemplate(TCellID n_id)
{
	bool isValid(false);

	std::vector<TCellID> n_ordered_edges = m_Front->orderedFrontEdgesAroundNode(m_mesh, n_id);
	int compteur_corner(0);
	int compteur_end(0);
	int compteur_reversal(0);
	// Compute the number of each singular edges around the node n
	for (auto e_id:n_ordered_edges) {
		Edge e = m_mesh->get<Edge>(e_id);
		if (m_EdgesClassification->value(e_id) == 1) {
			compteur_corner += 1;
		}
		else if (m_EdgesClassification->value(e_id) == 2) {
			compteur_end += 1;
		}
		else if (m_EdgesClassification->value(e_id) == 3) {
			compteur_reversal += 1;
		}
	}

	// Check the valid configurations
	if (compteur_corner == 3 && compteur_end == 0 && compteur_reversal == 0)
	{
		isValid = true;	// Template 3 CORNER
	}
	else if (compteur_corner==2 && compteur_end==1 && compteur_reversal==0)
	{
		isValid = true;	// Template 2 CORNER, 1 END
	}
	else if (compteur_corner==1 && compteur_end==2 && compteur_reversal==0)
	{
		std::vector<TCellID> end_edges;
		TCellID corner_edge;
		for (auto e_loc_id:n_ordered_edges)
		{
			if (m_EdgesClassification->value(e_loc_id)==1)
			{
				corner_edge = e_loc_id;
			}
			else if (m_EdgesClassification->value(e_loc_id)==2)
			{
				end_edges.push_back(e_loc_id);
			}
		}
		NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n_id);
		n_neighbourhood.execute();
		std::vector<TCellID> faces_id = n_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(end_edges[0], end_edges[1], corner_edge);
		if (faces_id.size()==3) {
			isValid = true;     // Template 1 CORNER, 2 END
		}
	}
	else if (compteur_end == 3 && n_ordered_edges.size()==3)
	{
		isValid = true;	// Template 3 END
	}
	else if (n_ordered_edges.size() == 6
	         && ( ( m_EdgesClassification->value(n_ordered_edges[0]) == 1
	              && m_EdgesClassification->value(n_ordered_edges[2]) == 1
	              && m_EdgesClassification->value(n_ordered_edges[4]) == 1
	              && m_EdgesClassification->value(n_ordered_edges[1]) == 2
	              && m_EdgesClassification->value(n_ordered_edges[3]) == 2
	              && m_EdgesClassification->value(n_ordered_edges[5]) == 2 )
	             || ( m_EdgesClassification->value(n_ordered_edges[1]) == 1
	                 && m_EdgesClassification->value(n_ordered_edges[3]) == 1
	                 && m_EdgesClassification->value(n_ordered_edges[5]) == 1
	                 && m_EdgesClassification->value(n_ordered_edges[0]) == 2
	                 && m_EdgesClassification->value(n_ordered_edges[2]) == 2
	                 && m_EdgesClassification->value(n_ordered_edges[4]) == 2 ) ) )
	{
		isValid = true;	// Template 3 CORNER, 3 END
	}
	else if (n_ordered_edges.size() == 6
	         && compteur_corner == 2 && compteur_end==2 && compteur_reversal==0)
	{
		isValid = true;	// Template 2 CORNER, 2 END
	}
	else if (n_ordered_edges.size() == 3
	         && compteur_corner == 2 && compteur_end==0 && compteur_reversal==1)
	{
		isValid = true;	// Template 2 CORNER, 1 REVERSAL
	}
	else if (n_ordered_edges.size() == 6
	         && compteur_corner == 0 && compteur_end==2 && compteur_reversal==1)
	{
		isValid = true;	// Template 2 END, 1 REVERSAL
	}

	return isValid;
}
/*------------------------------------------------------------------------*/
FrontEdgesNodesClassification_3D::Global_Feature_Edge
FrontEdgesNodesClassification_3D::ComputeOneGFE(TCellID n_id, TCellID e_id)
{
	Global_Feature_Edge new_GFE;
	new_GFE.Start_n_id = n_id;
	new_GFE.edges_id.push_back(e_id);

	TInt mark_EdgeUsed = m_mesh->newMark<Edge>();

	Edge e_loc = m_mesh->get<Edge>(e_id);
	Node n_loc = e_loc.getOppositeNode(m_mesh->get<Node>(n_id));
	while (m_NodesClassification->value(n_loc.id()) == 2
	       && !m_mesh->isMarked(e_loc, mark_EdgeUsed)
	       && n_loc.id() != new_GFE.Start_n_id
	       && m_EdgesClassification->value(e_loc.id()) == m_EdgesClassification->value(new_GFE.edges_id[0]))
	{
	   m_mesh->mark(e_loc, mark_EdgeUsed);
	   std::vector<Edge> n_loc_edges = n_loc.get<Edge>();
	   bool nextEdgeFound(false);
	   for (auto const& e_adj:n_loc_edges)
	   {
	      if ( e_adj.id() != e_loc.id()
	          && m_EdgesClassification->value(e_adj.id()) > 0
	          && !nextEdgeFound)
	      {
	         e_loc = e_adj ;
	         new_GFE.edges_id.push_back(e_loc.id());
	         n_loc = e_loc.getOppositeNode(n_loc);
	         nextEdgeFound = true;
	      }
	   }
	}

	new_GFE.End_n_id = n_loc.id() ;

	m_mesh->unmarkAll<Edge>(mark_EdgeUsed);
	m_mesh->freeMark<Edge>(mark_EdgeUsed);

	return new_GFE;

}
/*------------------------------------------------------------------------*/
std::vector<FrontEdgesNodesClassification_3D::Global_Feature_Edge>
FrontEdgesNodesClassification_3D::ComputeAllGFE()
{
	std::vector<Global_Feature_Edge> All_GFE;
	TInt mark_UsedEdges = m_mesh->newMark<Edge>();

	for (auto n_id:m_Front->getNodes())
	{
		if (m_NodesClassification->value(n_id) >= 3
		    && isValidNodeForTemplate(n_id))
		{
			Node n = m_mesh->get<Node>(n_id);
			std::vector<Edge> n_edges = n.get<Edge>();
			for (auto const& e:n_edges)
			{
				if (m_EdgesClassification->value(e.id()) > 0
				    && !m_mesh->isMarked(e, mark_UsedEdges))
				{
					Global_Feature_Edge new_GFE = ComputeOneGFE(n_id, e.id()) ;
					All_GFE.push_back(new_GFE);
					for (auto e_loc_id:new_GFE.edges_id)
					{
						m_mesh->mark(m_mesh->get<Edge>(e_loc_id), mark_UsedEdges);
					}
				}
			}
		}
	}

	m_mesh->unmarkAll<Edge>(mark_UsedEdges);
	m_mesh->freeMark<Edge>(mark_UsedEdges);

	return All_GFE;
}
/*------------------------------------------------------------------------*/
bool
FrontEdgesNodesClassification_3D::isThisPathValidForTemplates(Global_Feature_Edge& GFE)
{
	bool isValid(true);

	if (m_NodesClassification->value(GFE.Start_n_id) < 3
	    || m_NodesClassification->value(GFE.End_n_id) < 3)
	{
		isValid = false;
	}
	if (!isValidNodeForTemplate(GFE.Start_n_id)
	    || !isValidNodeForTemplate(GFE.End_n_id))
	{
		isValid = false;
	}
	if (GFE.edges_id.size() < 3)	// To improve: if the number of edges is under 3, in some case, the algorithm can perform. (ex: if the 2 opposite nodes are corners). Need to list the valid cases.
	{
		isValid = false;
	}
	if (GFE.edges_id.size() < 3
	    && singleNodeClassification(GFE.Start_n_id) == 1
	    && singleNodeClassification(GFE.End_n_id) == 1)
	{
		isValid = true;
	}
	if (GFE.Start_n_id == GFE.End_n_id
	    && GFE.edges_id.size() > 4)		// Loop paths
	{
		isValid = true;
	}

	return isValid;
}
/*------------------------------------------------------------------------*/
void
FrontEdgesNodesClassification_3D::ComputeValid_GFE()
{
	// Init
	m_All_global_feature_edges = ComputeAllGFE();

	bool isAllTreated(false);

	while (!isAllTreated && !m_All_global_feature_edges.empty())
	{
		isAllTreated = true;
		for (auto GFE:m_All_global_feature_edges)
		{
			if (!isThisPathValidForTemplates(GFE))
			{
				m_NodesClassification->set(GFE.Start_n_id, m_NodesClassification->value(GFE.Start_n_id)-1);
				m_NodesClassification->set(GFE.End_n_id, m_NodesClassification->value(GFE.End_n_id)-1);
				for (auto edge_id:GFE.edges_id)
				{
					m_EdgesClassification->set(edge_id, 0);
				}
				isAllTreated = false;
			}
		}
		// Re-compute the list of Global Feature Edges
		m_All_global_feature_edges = ComputeAllGFE();
	}

	// Mark the final selected nodes and edges for the templates
	for (auto const &GFE:m_All_global_feature_edges)
	{
		m_mesh->mark(m_mesh->get<Node>(GFE.Start_n_id), m_mark_NodesForTemplates);
		m_mesh->mark(m_mesh->get<Node>(GFE.End_n_id), m_mark_NodesForTemplates);
		for (auto edge_id:GFE.edges_id)
		{
			m_mesh->mark(m_mesh->get<Edge>(edge_id), m_mark_EdgesForTemplates);
		}
	}

	/*
	for (auto GFE:m_All_global_feature_edges)
	{
		std::cout << "---------" << std::endl;
		std::cout << "Starting point: " << GFE.Start_n_id << std::endl;
		std::cout << "Ending point: " << GFE.End_n_id << std::endl;
		std::cout << "Nbr edges: " << GFE.edges_id.size() << std::endl;
		std::cout << "Node " << GFE.Start_n_id << ", adj to " << m_NodesClassification->value(GFE.Start_n_id) << std::endl;
		std::cout << "Node " << GFE.End_n_id << ", adj to " << m_NodesClassification->value(GFE.End_n_id) << std::endl;
		std::cout << "Valid ? " << isThisPathValidForTemplates(GFE) << std::endl;
		std::cout << "---------" << std::endl;
	}
	*/

}
/*------------------------------------------------------------------------*/
void
FrontEdgesNodesClassification_3D::ComputeValidLoop_GFE()
{
	Variable<int>* var_node_couche_id = m_mesh->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	TInt mark_EdgesUsed = m_mesh->newMark<Edge>();
	for (auto e_id:m_mesh->edges())
	{
		Edge e = m_mesh->get<Edge>(e_id);
		std::vector<Node> e_nodes = e.get<Node>();
		if (var_node_couche_id->value(e_nodes[0].id()) == m_Front->getFrontID()
		    && var_node_couche_id->value(e_nodes[1].id()) == m_Front->getFrontID()
		    && !m_mesh->isMarked(e, mark_EdgesUsed)
		    && !m_mesh->isMarked(e, m_mark_EdgesForTemplates))
		{
			Global_Feature_Edge GFE = ComputeOneGFE(e_nodes[0].id(), e_id);
			if (isThisPathValidForTemplates(GFE))
			{
				m_All_global_feature_edges.push_back(GFE);
				for (auto e_loc_id:GFE.edges_id)
				{
					m_mesh->mark(m_mesh->get<Edge>(e_loc_id), m_mark_EdgesForTemplates);
					m_mesh->mark(m_mesh->get<Edge>(e_loc_id), mark_EdgesUsed);
				}
			}
			else
			{
				for (auto e_loc_id:GFE.edges_id)
				{
					//std::cout << e_loc_id << std::endl;
				}
			}
		}
	}
	m_mesh->unmarkAll<Edge>(mark_EdgesUsed);
	m_mesh->freeMark<Edge>(mark_EdgesUsed);

	/*
	std::cout << "Feature edges loop..." << std::endl;
	for (auto GFE:m_All_global_feature_edges)
	{
	   std::cout << "---------" << std::endl;
	   std::cout << "Starting point: " << GFE.Start_n_id << std::endl;
	   std::cout << "Ending point: " << GFE.End_n_id << std::endl;
	}
	*/
}
/*------------------------------------------------------------------------*/