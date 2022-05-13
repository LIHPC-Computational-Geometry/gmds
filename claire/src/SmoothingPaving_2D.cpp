//
// Created by rochec on 09/05/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/math/Line.h>
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

	InitializeMapLd();

	FrontNodeSmoothing();
	//FrontNodeSmoothing();
	//FrontNodeSmoothing();

	m_mesh->unmarkAll<Node>(m_markFrontNodes);
	m_mesh->freeMark<Node>(m_markFrontNodes);

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
void SmoothingPaving_2D::InitializeMapLd()
{
	std::vector<TCellID> nodes_id = m_Front.getNodes();

	for (auto n_id:nodes_id)
	{
		Node n = m_mesh->get<Node>(n_id);
		std::vector<Node> adj_nodes = math::Utils::AdjacentNodes(m_mesh, n);
		Node nj;
		std::vector<Node> ni_adj;
		for (auto n_adj:adj_nodes)
		{
			if (!m_mesh->isMarked<Node>(n_adj.id(),m_markFrontNodes))
			{
				nj = n_adj;
				math::Vector3d Vd = n.point()-nj.point() ;
				m_map_ld[n_id] = Vd.norm();
			}
		}
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
			math::Vector3d Dc = ComputeDc(n_id);
			Node n = m_mesh->get<Node>(n_id);
			n.setPoint(n.point() + (Db+Dc)/2.0);
		}

	}
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
math::Vector3d
SmoothingPaving_2D::ComputeDa(TCellID n_id)
{
	math::Vector3d Vip(0,0,0);
	Node n = m_mesh->get<Node>(n_id);
	std::vector<Face> faces = n.get<Face>();

	for (auto f:faces)
	{
		std::vector<Node> nodes = f.get<Node>();

		for (auto n:nodes)
		{
			if ( (n.id() != n_id)
			    && (math::Utils::CommonEdge(m_mesh, n.id(), n_id) != NullID) )
			{
				Vip = Vip + m_map_old_coords[n.id()];
			}
			else if (n.id() != n_id)
			{
				Vip = Vip - m_map_old_coords[n.id()];
			}
		}


	}

	Vip = Vip / faces.size();

	math::Vector3d Da;
	Da = Vip - n.point();

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
	math::Vector3d Vi = m_map_old_coords[n.id()];
	math::Vector3d Vj = m_map_old_coords[nj.id()];
	double la = Da.norm();
	//double ld = (Vi-Vj).norm() ;
	double ld = m_map_ld[n_id] ;

	math::Vector3d Db;
	Db = Vj-Vi + (Da+Vi-Vj)*(ld/la);
	return Db;

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
math::Vector3d
SmoothingPaving_2D::ComputeDc(TCellID n_id)
{
	Node n = m_mesh->get<Node>(n_id);
	std::vector<Node> adj_nodes = math::Utils::AdjacentNodes(m_mesh, n);
	Node nj;
	std::vector<Node> ni_adj;
	for (auto n_adj:adj_nodes)
	{
		if (!m_mesh->isMarked<Node>(n_adj.id(),m_markFrontNodes))
		{
			nj = n_adj;
		}
		else if(n_adj.id() != n.id())
		{
			ni_adj.push_back(n_adj);
		}
	}

	math::Vector3d Pi = n.point()-nj.point() ;
	math::Vector3d Pim1 = m_map_old_coords[ni_adj[0].id()] - nj.point() ;
	math::Vector3d Pip1 = m_map_old_coords[ni_adj[1].id()] - nj.point() ;
	math::Vector3d Pb1 = Pim1 + Pip1 ;

	// Calcul de la direction
	math::Vector3d Pb2 = Pi + Pb1;
	Pb2.normalize();

	// Calcul de la norme
	math::Point Q;
	math::Line L1(nj.point(),nj.point() + Pb1);
	math::Line L2(m_map_old_coords[ni_adj[0].id()],m_map_old_coords[ni_adj[2].id()]);
	double param;
	L1.intersect2D(L2, Q, param);
	math::Vector3d QNj = Q - nj.point() ;
	math::Vector3d NiNj = nj.point() - n.point();
	//double ld = NiNj.norm() ;
	double ld = m_map_ld[n_id] ;
	double lq = QNj.norm() ;

	if(ld > lq)
	{
		Pb2 = Pb2*( lq+ld )/2.0;
	}
	else
	{
		Pb2 = ld*Pb2;
	}

	math::Vector3d Dc = Pb2-Pi;

	return Dc;

}
/*------------------------------------------------------------------------*/
