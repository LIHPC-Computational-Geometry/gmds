//
// Created by rochec on 09/05/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/claire/Utils.h>
#include <gmds/claire/SmoothingPaving_2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

SmoothingPaving_2D::SmoothingPaving_2D(Mesh *AMesh, Front AFront) {
	m_mesh = AMesh;
	m_Front = AFront;
}


/*------------------------------------------------------------------------*/
SmoothingPaving_2D::STATUS
SmoothingPaving_2D::execute()
{
	std::vector<TCellID> nodes_id = m_Front.getNodes();
	for (auto n_id:nodes_id)
	{
		m_mesh->mark<Node>(n_id, m_markFrontNodes);
	}

	FrontNodeSmoothing();
	FrontNodeSmoothing();
	FrontNodeSmoothing();

	return SmoothingPaving_2D::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void SmoothingPaving_2D::UpdateOldCoords()
{
	for (auto n_id:m_mesh->nodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		m_map_old_coords[n_id] = n.point() ;
	}
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void
SmoothingPaving_2D::FrontNodeSmoothing()
{
	UpdateOldCoords();

	std::vector<TCellID> front_nodes = m_Front.getNodes();
	for (auto n_id:front_nodes)
	{
		Node n = m_mesh->get<Node>(n_id);
		std::vector<Face> adj_faces = n.get<Face>();
		if (adj_faces.size()!=2) {
			math::Vector3d Da = ComputeDa(n_id);
			Node n = m_mesh->get<Node>(n_id);
			n.setPoint(n.point() + Da);
		}
		else
		{
			math::Vector3d Db = ComputeDb(n_id);
			Node n = m_mesh->get<Node>(n_id);
			n.setPoint(n.point() + Db);
		}
	}
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
math::Vector3d
SmoothingPaving_2D::ComputeDa(TCellID n_id)
{
	math::Point n_new_pos;
	Node n = m_mesh->get<Node>(n_id);
	std::vector<Face> faces = n.get<Face>();

	for (auto f:faces)
	{
		std::vector<Node> nodes = f.get<Node>();
		std::vector<Node> ord_nodes;
		for (auto n:nodes)
		{
			if ( (n.id() != n_id)
			    && (math::Utils::CommonEdge(m_mesh, n.id(), n_id) != NullID) )
			{
				n_new_pos = n_new_pos + m_map_old_coords[n.id()];
			}
			else if (n.id() != n_id)
			{
				n_new_pos = n_new_pos - m_map_old_coords[n.id()];
			}
		}


	}

	n_new_pos.setX(n_new_pos.X()/faces.size() );
	n_new_pos.setY(n_new_pos.Y()/faces.size() );
	n_new_pos.setZ(n_new_pos.Z()/faces.size() );

	math::Vector3d Da;
	Da = n_new_pos - n.point();

	return Da;

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
math::Vector3d
SmoothingPaving_2D::ComputeDb(TCellID n_id)
{
	Node n = m_mesh->get<Node>(n_id);
	std::vector<Node> adj_nodes = math::Utils::AdjacentNodes(m_mesh, n);
	Node nj;
	for (auto n_adj:adj_nodes)
	{
		if (!m_mesh->isMarked<Node>(n_adj.id(),m_markFrontNodes))
		{
			nj = n_adj;
		}
	}

	math::Vector3d Da = ComputeDa(n_id);
	math::Vector3d Vi = n.point();
	math::Vector3d Vj = nj.point();
	math::Vector3d Vip = Vi + Da;
	double ld = (Vi-Vj).norm() ;
	double la = (Vip-Vj).norm();

	math::Vector3d Db;
	Db = Vj-Vi + (Da+Vi-Vj)*(ld/la);
	return Db;

}
/*------------------------------------------------------------------------*/
