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

	m_global_feature_edges = ComputeAllGlobalFeatureEdge();

	ComputeNodesEdgesForTemplates();

	return FrontEdgesNodesClassification_3D::SUCCESS;
}
/*------------------------------------------------------------------------*/
int
FrontEdgesNodesClassification_3D::getMarkSemiEdges()
{
	return m_mark_semiEdges;
}
/*------------------------------------------------------------------------*/
int
FrontEdgesNodesClassification_3D::getMarkSemiNodes()
{
	return m_mark_semiNodes;
}
/*------------------------------------------------------------------------*/
int
FrontEdgesNodesClassification_3D::getMarkEdgesTemplates()
{
	return m_mark_EdgesForTemplates;
}
/*------------------------------------------------------------------------*/
int
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
	for (auto f:e_faces)
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
		std::vector<Node> e_nodes = e.get<Node>() ;
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
		for (auto e:n_edges)
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
			for (auto e:n_edges)
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
std::vector<TCellID>
FrontEdgesNodesClassification_3D::ComputeOneGlobalFeatureEdge(TCellID n_id, TCellID e_id)
{
	std::vector<TCellID> global_feature_edge;
	global_feature_edge.push_back(e_id);

	int mark_UsedEdges = m_mesh->newMark<Edge>();

	Edge e_feature_loc = m_mesh->get<Edge>(e_id);
	Node n_loc = e_feature_loc.getOppositeNode(m_mesh->get<Node>(n_id));
	while (!m_mesh->isMarked(e_feature_loc, m_mark_semiEdges)
	       && m_NodesClassification->value(n_loc.id()) == 2
	       && !m_mesh->isMarked(e_feature_loc, mark_UsedEdges))
	{
		m_mesh->mark(e_feature_loc, mark_UsedEdges);
		std::vector<Edge> n_loc_edges = n_loc.get<Edge>();
		bool nextEdgeFound(false);
		for (auto e_adj:n_loc_edges)
		{
			if ( e_adj.id() != e_feature_loc.id()
			    && m_EdgesClassification->value(e_adj.id()) > 0
			    && !nextEdgeFound)
			{
				e_feature_loc = e_adj ;
				global_feature_edge.push_back(e_feature_loc.id());
				n_loc = e_feature_loc.getOppositeNode(n_loc);
				nextEdgeFound = true;
			}
		}
	}

	m_mesh->unmarkAll<Edge>(mark_UsedEdges);
	m_mesh->freeMark<Edge>(mark_UsedEdges);

	return global_feature_edge;
}
/*------------------------------------------------------------------------*/
std::vector<std::vector<TCellID>>
FrontEdgesNodesClassification_3D::ComputeAllGlobalFeatureEdge()
{
	std::vector<std::vector<TCellID>> all_glob_feature_edge;
	int mark_UsedEdges = m_mesh->newMark<Edge>();

	for (auto n_id:m_Front->getNodes())
	{
		if (m_NodesClassification->value(n_id) >= 3)
		{
			Node n = m_mesh->get<Node>(n_id);
			std::vector<Edge> n_edges = n.get<Edge>();
			for (auto e:n_edges)
			{
				if (m_EdgesClassification->value(e.id()) > 0
				    && !m_mesh->isMarked(e, mark_UsedEdges))
				{
					std::vector<TCellID> global_feature_edge = ComputeOneGlobalFeatureEdge(n_id, e.id()) ;
					bool isSemiGlobalEdge(false);
					for (auto e_loc_id:global_feature_edge)
					{
						m_mesh->mark( m_mesh->get<Edge>(e_loc_id), mark_UsedEdges);
						if (m_mesh->isMarked(m_mesh->get<Edge>(e_loc_id), m_mark_semiEdges))
						{
							isSemiGlobalEdge = true;
						}
					}
					if (!isSemiGlobalEdge)
					{
						all_glob_feature_edge.push_back(global_feature_edge);
					}
				}
			}
		}
	}

	m_mesh->unmarkAll<Edge>(mark_UsedEdges);
	m_mesh->freeMark<Edge>(mark_UsedEdges);

	return all_glob_feature_edge;
}
/*------------------------------------------------------------------------*/
void
FrontEdgesNodesClassification_3D::ComputeNodesEdgesForTemplates()
{
	// Initialization
	for (auto global_feature_edge:m_global_feature_edges)
	{
		for (auto edge_id:global_feature_edge)
		{
			m_mesh->mark(m_mesh->get<Edge>(edge_id), m_mark_EdgesForTemplates);
			std::vector<Node> nodes = (m_mesh->get<Edge>(edge_id)).get<Node>();
			if (m_NodesClassification->value( nodes[0].id() ) >= 3)
			{
				m_mesh->mark(m_mesh->get<Node>(nodes[0].id()), m_mark_NodesForTemplates);
			}
			if (m_NodesClassification->value( nodes[1].id() ) >= 3)
			{
				m_mesh->mark(m_mesh->get<Node>(nodes[1].id()), m_mark_NodesForTemplates);
			}
		}
	}

	bool isAllMarked(false);

	while (!isAllMarked)
	{
		isAllMarked=true;
		for (auto n_id:m_Front->getNodes())
		{
			if (m_mesh->isMarked(m_mesh->get<Node>(n_id), m_mark_NodesForTemplates))
			{
				NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n_id);
				n_neighbourhood.execute();
				int compteur_marked_edges(0);
				for (auto edge_id:n_neighbourhood.getOrderedEdges())
				{
					if (m_mesh->isMarked(m_mesh->get<Edge>(edge_id), m_mark_EdgesForTemplates))
					{
						compteur_marked_edges++;
					}
				}
				if (compteur_marked_edges < 3)
				{
					m_mesh->unmark(m_mesh->get<Node>(n_id), m_mark_NodesForTemplates);
					isAllMarked = false;
				}
			}
		}

		for (auto global_feature_edge:m_global_feature_edges)
		{
			Edge e_0 = m_mesh->get<Edge>(global_feature_edge[0]);
			Edge e_last = m_mesh->get<Edge>(global_feature_edge[global_feature_edge.size()-1]);

			std::vector<Node> e_0_nodes = e_0.get<Node>();
			std::vector<Node> e_last_nodes = e_last.get<Node>();

			bool isGlobFeatEdgeAlreadyMarked = (m_mesh->isMarked(e_0, m_mark_EdgesForTemplates)
			    || m_mesh->isMarked(e_last, m_mark_EdgesForTemplates) ) ;

			bool isStartGFEMarked = ( m_mesh->isMarked(e_0_nodes[0], m_mark_NodesForTemplates)
			    || m_mesh->isMarked(e_0_nodes[1], m_mark_NodesForTemplates) ) ;

			bool isEndGBFMarked = ( m_mesh->isMarked(e_last_nodes[0], m_mark_NodesForTemplates)
			        || m_mesh->isMarked(e_last_nodes[1], m_mark_NodesForTemplates) ) ;

			if ( isGlobFeatEdgeAlreadyMarked
			    && ( !isStartGFEMarked || !isEndGBFMarked) )
			{
				isAllMarked=false;
				for (auto edge_loc_id:global_feature_edge)
				{
					m_mesh->unmark(m_mesh->get<Edge>(edge_loc_id), m_mark_EdgesForTemplates);
					m_mesh->unmark(e_0_nodes[0], m_mark_NodesForTemplates);
					m_mesh->unmark(e_0_nodes[1], m_mark_NodesForTemplates);
					m_mesh->unmark(e_last_nodes[0], m_mark_NodesForTemplates);
					m_mesh->unmark(e_last_nodes[1], m_mark_NodesForTemplates);
				}
			}
		}

	}

}