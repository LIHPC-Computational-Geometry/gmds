//
// Created by rochec on 20/12/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/claire/FrontEdgesNodesClassification_3D.h>
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
}
/*------------------------------------------------------------------------*/
FrontEdgesNodesClassification_3D::~FrontEdgesNodesClassification_3D()
{
	m_mesh->unmarkAll<Edge>(m_mark_semiEdges);
	m_mesh->freeMark<Edge>(m_mark_semiEdges);
	m_mesh->unmarkAll<Node>(m_mark_semiNodes);
	m_mesh->freeMark<Node>(m_mark_semiNodes);
}
/*------------------------------------------------------------------------*/
FrontEdgesNodesClassification_3D::STATUS
FrontEdgesNodesClassification_3D::execute()
{
	FrontEdgesClassification();	// Fill the variable m_EdgesClassification
	FrontNodesClassification();	// Fill the variable m_NodesClassification

	return FrontEdgesNodesClassification_3D::SUCCESS;
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
		edge_classification = 1;
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
