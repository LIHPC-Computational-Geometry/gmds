//
// Created by rochec on 21/12/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/claire/NodeNeighbourhoodOnFront_3D.h>
#include <gmds/claire/Utils.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/
NodeNeighbourhoodOnFront_3D::NodeNeighbourhoodOnFront_3D(Mesh *AMesh, Front_3D *AFront, TCellID An_id) :
  m_mesh(AMesh),
  m_Front(AFront)
{
	m_n_id = An_id;
}
/*------------------------------------------------------------------------*/
NodeNeighbourhoodOnFront_3D::~NodeNeighbourhoodOnFront_3D()
{

}
/*------------------------------------------------------------------------*/
NodeNeighbourhoodOnFront_3D::STATUS
NodeNeighbourhoodOnFront_3D::execute()
{
	orderedFrontEdgesAroundNode();	// Fill the vector m_orderedEdges
	orderedFrontFacesAroundNode();	// Fill the vector m_orderedFaces

	return NodeNeighbourhoodOnFront_3D::SUCCESS;
}
/*------------------------------------------------------------------------*/
std::vector<TCellID>
NodeNeighbourhoodOnFront_3D::getOrderedEdges()
{
	return m_orderedEdges;
}
/*------------------------------------------------------------------------*/
std::vector<TCellID>
NodeNeighbourhoodOnFront_3D::getOrderedFaces()
{
	return m_orderedFaces;
}
/*------------------------------------------------------------------------*/
std::vector<TCellID>
NodeNeighbourhoodOnFront_3D::adjFacesToEdge(TCellID e_id)
{
	std::vector<TCellID> adj_faces;
	int position(0);
	for (int i=0; i<= m_orderedEdges.size()-1;i++)
	{
		if (m_orderedEdges[i]==e_id)
		{
			position = i;
		}
	}

	if (position==0)
	{
		adj_faces.push_back(m_orderedFaces[m_orderedFaces.size()-1]);
		adj_faces.push_back(m_orderedFaces[0]);
	}
	else
	{
		adj_faces.push_back(m_orderedFaces[position-1]);
		adj_faces.push_back(m_orderedFaces[position]);
	}

	return adj_faces;
}
/*------------------------------------------------------------------------*/
TCellID
NodeNeighbourhoodOnFront_3D::nextEdgeOfFace(TCellID f_id, TCellID e_id)
{
	TCellID next_edge_id;
	std::vector<Edge> f_edges = (m_mesh->get<Face>(f_id)).get<Edge>();
	for (auto n_edges_id:m_orderedEdges)
	{
		if ( n_edges_id != e_id
		    && (n_edges_id == f_edges[0].id() || n_edges_id == f_edges[1].id()
		        || n_edges_id == f_edges[2].id() || n_edges_id == f_edges[3].id()))
		{
			next_edge_id = n_edges_id ;
		}
	}
	return next_edge_id;
}
/*------------------------------------------------------------------------*/
TCellID
NodeNeighbourhoodOnFront_3D::adjFaceToEdge1InEdge2SideAvoidingEdge3(TCellID e1_id, TCellID e2_id, TCellID e3_id)
{
	TCellID face_id;

	std::vector<TCellID> adj_faces = (*this).adjFacesToEdge(e1_id);

	TCellID f_id = adj_faces[0];
	TCellID e_id = e1_id ;
	bool faceFound(false);
	while (!faceFound)
	{
		TCellID next_edge_id = (*this).nextEdgeOfFace(f_id, e_id);
		if (next_edge_id==e2_id)
		{
			face_id=adj_faces[0];
			faceFound = true;
		}
		else if (next_edge_id==e3_id)
		{
			face_id=adj_faces[1];
			faceFound = true;
		}
		else
		{
			e_id = next_edge_id;
			std::vector<TCellID> e_faces = (*this).adjFacesToEdge(e_id) ;
			if (e_faces[0] == f_id)
			{
				f_id = e_faces[1];
			}
			else
			{
				f_id = e_faces[0];
			}
		}
	}

	return face_id;
}
/*------------------------------------------------------------------------*/
std::vector<TCellID>
NodeNeighbourhoodOnFront_3D::facesBtwEdge1nEdge2AvoidingEdge3(TCellID e1_id, TCellID e2_id, TCellID e3_id)
{
	std::vector<TCellID> faces;

	std::vector<TCellID> adj_faces = (*this).adjFacesToEdge(e1_id);
	std::vector<TCellID> list_faces;

	TCellID f_id = adj_faces[0];
	TCellID e_id = e1_id ;
	bool faceFound(false);
	int max_iter(0);
	while (!faceFound && max_iter < 100)
	{
		TCellID next_edge_id = (*this).nextEdgeOfFace(f_id, e_id);
		list_faces.push_back(f_id);
		if (next_edge_id==e2_id)
		{
			faces = list_faces;
			faceFound = true;
		}
		else if (next_edge_id==e3_id)
		{
			f_id=adj_faces[1];
			faceFound = true;
		}
		else
		{
			e_id = next_edge_id;
			std::vector<TCellID> e_faces = (*this).adjFacesToEdge(e_id) ;
			if (e_faces[0] == f_id)
			{
				f_id = e_faces[1];
			}
			else
			{
				f_id = e_faces[0];
			}
		}
		max_iter++;
	}

	if (max_iter == 100)
	{
		std::cout << "ATTENTION NodeNeighbourhoodOnFront: max iteration" << std::endl;
	}

	if (f_id==adj_faces[1] && faceFound)
	{
		list_faces.clear();
		f_id = adj_faces[1];
		e_id = e1_id ;
		faceFound = false;
		while (!faceFound)
		{
			TCellID next_edge_id = (*this).nextEdgeOfFace(f_id, e_id);
			list_faces.push_back(f_id);
			if (next_edge_id==e2_id)
			{
				faces = list_faces;
				faceFound = true;
			}
			else
			{
				e_id = next_edge_id;
				std::vector<TCellID> e_faces = (*this).adjFacesToEdge(e_id) ;
				if (e_faces[0] == f_id)
				{
					f_id = e_faces[1];
				}
				else
				{
					f_id = e_faces[0];
				}
			}
		}
	}

	return faces;
}
/*------------------------------------------------------------------------*/
void
NodeNeighbourhoodOnFront_3D::orderedFrontEdgesAroundNode()
{
	Node n = m_mesh->get<Node>(m_n_id);
	std::vector<Edge> n_edges = n.get<Edge>();
	Variable<int>* var_node_couche_id = m_mesh->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");

	// Get the front edges connected to n
	std::vector<Edge> n_edges_on_Front;
	for (auto const& e:n_edges)
	{
		Node n_opp = e.getOppositeNode(n);
		if (var_node_couche_id->value(n_opp.id()) == m_Front->getFrontID())
		{
			n_edges_on_Front.push_back(e);
		}
	}

	TInt mark_isTreated = m_mesh->newMark<Face>();
	m_orderedEdges.push_back(n_edges_on_Front[0].id());	// Choose a first edge, and a first face

	for (int i=1;i<=n_edges_on_Front.size()-1;i++)
	{
		TCellID e_id = m_orderedEdges[m_orderedEdges.size()-1];
		Edge e = m_mesh->get<Edge>(e_id);
		Node n_opp = e.getOppositeNode(m_n_id);

		Face f;
		std::vector<Face> e_faces = e.get<Face>() ;
		for (auto const& f_loc:e_faces)
		{
			std::vector<Node> f_loc_nodes = f_loc.get<Node>();

			if (!m_mesh->isMarked(f_loc, mark_isTreated)
			    && var_node_couche_id->value(f_loc_nodes[0].id()) == m_Front->getFrontID()
			    && var_node_couche_id->value(f_loc_nodes[1].id()) == m_Front->getFrontID()
			    && var_node_couche_id->value(f_loc_nodes[2].id()) == m_Front->getFrontID()
			    && var_node_couche_id->value(f_loc_nodes[3].id()) == m_Front->getFrontID())
			{
				f = f_loc;
			}
		}
		m_mesh->mark(f, mark_isTreated);

		std::vector<Edge> f_edges = f.get<Edge>();
		for (auto const& e_loc:f_edges)
		{
			std::vector<Node> e_loc_nodes = e_loc.get<Node>();
			if ( (e_loc_nodes[0].id() == m_n_id && e_loc_nodes[1].id() != n_opp.id() )
			    || (e_loc_nodes[1].id() == m_n_id && e_loc_nodes[0].id() != n_opp.id() ) )
			{
				m_orderedEdges.push_back(e_loc.id());
			}
		}

	}

	m_mesh->unmarkAll<Face>(mark_isTreated);
	m_mesh->freeMark<Face>(mark_isTreated);

}
/*------------------------------------------------------------------------*/
void
NodeNeighbourhoodOnFront_3D::orderedFrontFacesAroundNode()
{
	for (int i=0; i<m_orderedEdges.size()-1; i++)
	{
		Edge e1 = m_mesh->get<Edge>(m_orderedEdges[i]) ;
		Edge e2 = m_mesh->get<Edge>(m_orderedEdges[i+1]) ;

		Node n_opp_1 = e1.getOppositeNode(m_mesh->get<Node>(m_n_id));
		Node n_opp_2 = e2.getOppositeNode(m_mesh->get<Node>(m_n_id));

		TCellID f_id = math::Utils::CommonFace3Nodes(m_mesh, n_opp_1.id(), m_n_id, n_opp_2.id()) ;
		m_orderedFaces.push_back(f_id);
	}

	Edge e1 = m_mesh->get<Edge>(m_orderedEdges[m_orderedEdges.size()-1]) ;
	Edge e2 = m_mesh->get<Edge>(m_orderedEdges[0]) ;

	Node n_opp_1 = e1.getOppositeNode(m_mesh->get<Node>(m_n_id));
	Node n_opp_2 = e2.getOppositeNode(m_mesh->get<Node>(m_n_id));

	TCellID f_id = math::Utils::CommonFace3Nodes(m_mesh, n_opp_1.id(), m_n_id, n_opp_2.id()) ;
	m_orderedFaces.push_back(f_id);

}
/*------------------------------------------------------------------------*/
std::vector<TCellID>
NodeNeighbourhoodOnFront_3D::facesBtwEdge1nEdge2inFaceSide(TCellID e1_id, TCellID e2_id, TCellID f_id)
{
	// Init the vector with the first face
	std::vector<TCellID> faces_id;
	faces_id.push_back(f_id);

	TCellID current_face_id = f_id;
	TCellID current_edge_id = e1_id;
	TCellID e_next_id = (*this).nextEdgeOfFace(f_id, e1_id);
	while (e_next_id != e2_id)
	{
		std::vector<TCellID> adj_faces = (*this).adjFacesToEdge(e_next_id);
		if (adj_faces[0] == current_face_id )
		{
			current_face_id = adj_faces[1];
		}
		else
		{
			current_face_id = adj_faces[0];
		}
		faces_id.push_back(current_face_id);
		current_edge_id = e_next_id;
		e_next_id = (*this).nextEdgeOfFace(current_face_id, current_edge_id);
	}

	return faces_id;
}