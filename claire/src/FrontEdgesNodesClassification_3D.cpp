//
// Created by rochec on 20/12/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/claire/FrontEdgesNodesClassification_3D.h>
#include <gmds/claire/NodeNeighbourhoodOnFront_3D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/
FrontEdgesNodesClassification_3D::FrontEdgesNodesClassification_3D(Mesh *AMesh, Front_3D *AFront, Variable<int>* A_EdgesClassification, Variable<int>* A_NodesClassification) :
  m_mesh(AMesh),
  m_Front(AFront)
{
	m_EdgesClassification = A_EdgesClassification;
	m_NodesClassification = A_NodesClassification;
	m_mark_semiEdges = m_mesh->newMark<Edge>();
	m_mark_semiNodes = m_mesh->newMark<Node>();
	m_mark_EdgesForTemplates = m_mesh->newMark<Edge>();
	m_mark_NodesForTemplates = m_mesh->newMark<Node>();
}
/*------------------------------------------------------------------------*/
FrontEdgesNodesClassification_3D::~FrontEdgesNodesClassification_3D()
{
	m_mesh->unmarkAll<Edge>(m_mark_semiEdges);
	m_mesh->freeMark<Edge>(m_mark_semiEdges);
	m_mesh->unmarkAll<Node>(m_mark_semiNodes);
	m_mesh->freeMark<Node>(m_mark_semiNodes);
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
	FrontNodesClassification();	// Fill the variable m_NodesClassification
	MarkSemiEdgesandNodes();		// Mark the semi nodes, and the semi edges
	/*
	m_global_feature_edges = ComputeAllGlobalFeatureEdge();
	ComputeNodesEdgesForTemplates();
	 */
	m_All_global_feature_edges = ComputeAllGFE();
	ComputeValid_GFE();

	return FrontEdgesNodesClassification_3D::SUCCESS;
}
/*------------------------------------------------------------------------*/
TInt
FrontEdgesNodesClassification_3D::getMarkSemiEdges()
{
	return m_mark_semiEdges;
}
/*------------------------------------------------------------------------*/
TInt
FrontEdgesNodesClassification_3D::getMarkSemiNodes()
{
	return m_mark_semiNodes;
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
std::vector<std::vector<TCellID>>
FrontEdgesNodesClassification_3D::getGlobalFeatureEdge()
{
	return m_global_feature_edges;
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

	// Each face of the front is connected to a unique region
	// So the size of the vector e_faces_on_front[0].get<Region>() is 1.
	Region r1 = (e_faces_on_front[0].get<Region>())[0];
	Region r2 = (e_faces_on_front[1].get<Region>())[0];

	math::Point r1_center = r1.center();
	math::Point r2_center = r2.center();

	math::Point f1_center = e_faces_on_front[0].center();
	math::Point f2_center = e_faces_on_front[1].center();

	math::Vector3d v1 = (f1_center-r1.center()).normalize() ;
	math::Vector3d v2 = (f2_center-r2.center()).normalize() ;

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
void
FrontEdgesNodesClassification_3D::MarkSemiEdgesandNodes()
{
	for (auto n_id:m_Front->getNodes())
	{
		if (m_NodesClassification->value(n_id) == 1)
		{
			Node n = m_mesh->get<Node>(n_id);
			m_mesh->mark(n, m_mark_semiNodes);	// Mark the semi node
			std::vector<Edge> n_edges = n.get<Edge>();
			for (auto const& e:n_edges)
			{
				if (m_EdgesClassification->value(e.id()) > 0)
				{
					m_mesh->mark(e, m_mark_semiEdges);	// Mark the semi edge
				}
			}
		}
	}
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
	       && !m_mesh->isMarked(e_loc, mark_EdgeUsed))
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
		if (m_NodesClassification->value(n_id) >= 3)
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
void
FrontEdgesNodesClassification_3D::ComputeValid_GFE()
{
	// Init
	for (auto const& GFE:m_All_global_feature_edges)
	{
		if (m_NodesClassification->value(GFE.End_n_id) >= 3
		    && isValidNodeForTemplate(GFE.End_n_id))
		{
			m_mesh->mark(m_mesh->get<Node>(GFE.End_n_id), m_mark_NodesForTemplates);
		}
		if (m_NodesClassification->value(GFE.Start_n_id) >= 3
		    && isValidNodeForTemplate(GFE.Start_n_id))
		{
			m_mesh->mark(m_mesh->get<Node>(GFE.Start_n_id), m_mark_NodesForTemplates);
		}
		for (auto e_loc_id:GFE.edges_id)
		{
			m_mesh->mark(m_mesh->get<Edge>(e_loc_id), m_mark_EdgesForTemplates);
		}
	}

	bool isAllTreated(false);

	while (!isAllTreated && !m_All_global_feature_edges.empty())
	{
		isAllTreated = true;
		auto it=m_All_global_feature_edges.begin();
		while (it != m_All_global_feature_edges.end())
		{
			Global_Feature_Edge GFE = *it;
			if (!m_mesh->isMarked(m_mesh->get<Node>(GFE.Start_n_id), m_mark_NodesForTemplates)
			    || !m_mesh->isMarked(m_mesh->get<Node>(GFE.End_n_id), m_mark_NodesForTemplates)
			    || (GFE.edges_id.size() <= 1 && m_Front->getFrontID() > 2))		// In this work, we can't ensure the topology validity of the patterns in case there is only one edge between two singular nodes
			{
					it = m_All_global_feature_edges.erase(it);
				   m_mesh->unmark(m_mesh->get<Node>(GFE.Start_n_id), m_mark_NodesForTemplates);
				   m_mesh->unmark(m_mesh->get<Node>(GFE.End_n_id), m_mark_NodesForTemplates);
				   for (auto e_loc_id:GFE.edges_id)
				   {
					   m_mesh->unmark(m_mesh->get<Edge>(e_loc_id), m_mark_EdgesForTemplates);
				   }
				   isAllTreated = false;
			}
			else
			{
				it++;
			}
		}
	}
	/*
	for (auto GFE:m_All_global_feature_edges)
	{
		std::cout << "---------" << std::endl;
		std::cout << "Starting point: " << GFE.Start_n_id << std::endl;
		std::cout << "Ending point: " << GFE.End_n_id << std::endl;
	}
	 */
}
/*------------------------------------------------------------------------*/