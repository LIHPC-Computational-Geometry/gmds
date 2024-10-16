//
// Created by rochec on 09/05/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/math/Line.h>
#include <gmds/aero/Utils.h>
#include <gmds/aero/SmoothingPaving_2D.h>
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

	// Lissage des noeuds du front
	//FrontNodeSmoothing();
	//FrontNodeSmoothing();
	//FrontNodeSmoothing();
	//FrontNodeSmoothing();
	//FrontNodeSmoothing();

	// Lissage des noeuds intérieurs au front
	//InteriorNodeSmoothing();
	//InteriorNodeSmoothing();

	m_mesh->unmarkAll<Node>(m_markFrontNodes);
	m_mesh->freeMark<Node>(m_markFrontNodes);

	return SmoothingPaving_2D::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void
SmoothingPaving_2D::UpdateOldCoords()
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
SmoothingPaving_2D::InitializeMapLd()
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

		if (adj_faces.size()!=2)
		{
			math::Vector3d Da = ComputeDa(n_id);
			n.setPoint(n.point() + Da);
		}

		if (adj_faces.size()==2)
		{
			math::Vector3d Db = ComputeDb(n_id);
			math::Vector3d Dc = ComputeDc(n_id);
			n.setPoint(n.point() + (Db+Dc)/2.0);
		}

	}

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
math::Vector3d
SmoothingPaving_2D::ComputeDa(TCellID n_id)
{
	math::Vector3d Vip({0,0,0});
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
				math::Point tmp_p= m_map_old_coords[n.id()];
				Vip += {tmp_p.X(),tmp_p.Y(),tmp_p.Z()};
			}
			else if ( (n.id() != n_id)
			         && (math::Utils::CommonEdge(m_mesh, n.id(), n_id) == NullID) )
			{
				math::Point tmp_p= m_map_old_coords[n.id()];
				Vip -= {tmp_p.X(),tmp_p.Y(),tmp_p.Z()};
			}
		}


	}

	Vip = Vip / faces.size();

	math::Vector3d Da;
	Da = Vip - math::Vector3d({n.X(),n.Y(),n.Z()});

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
	math::Vector3d Vi = math::vec(m_map_old_coords[n.id()]);
	math::Vector3d Vj = math::vec(m_map_old_coords[nj.id()]);
	math::Vector3d Vip = Vi+Da;
	double la = (Vip-Vj).norm();
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

	math::Vector3d Pi = n.point()-m_map_old_coords[nj.id()] ;
	math::Vector3d Pim1 = m_map_old_coords[ni_adj[0].id()] - m_map_old_coords[nj.id()] ;
	math::Vector3d Pip1 = m_map_old_coords[ni_adj[1].id()] - m_map_old_coords[nj.id()] ;
	math::Vector3d Pb1 = Pim1 + Pip1 ;

	// Calcul de la direction
	math::Vector3d Pb2 = Pi + Pb1;
	Pb2.normalize();

	// Calcul de la norme
	math::Point Q;
	math::Line L1(m_map_old_coords[nj.id()],m_map_old_coords[nj.id()] + Pb2);
	math::Line L2(m_map_old_coords[ni_adj[0].id()],m_map_old_coords[ni_adj[1].id()]);
	double param;
	bool intersection = L1.intersect2D(L2, Q, param);
	//std::cout << "Intersection trouvée ? " << intersection << std::endl;
	if (!intersection)
	{
		Q = n.point();
	}
	math::Vector3d QNj = Q - m_map_old_coords[nj.id()] ;
	/*
	std::cout << "------------" << std::endl;
	std::cout << "n id : " << n_id << std::endl;
	std::cout << "param : " << param << std::endl;
	std::cout << "Q : " << Q << std::endl;
	std::cout << "Pb1 : " << Pb1 << std::endl;
	std::cout << "Autre vecteur : " << m_map_old_coords[ni_adj[1].id()] - m_map_old_coords[ni_adj[0].id()] << std::endl;
	std::cout << "QNj : " << QNj << std::endl;
	 */

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


/*------------------------------------------------------------------------*/
void
SmoothingPaving_2D::InteriorNodeSmoothing()
{
	UpdateOldCoords();

	std::vector<TCellID> front_nodes = m_Front.getNodes();
	Variable<int> * var_couche = m_mesh->getVariable<int,GMDS_NODE>("GMDS_Couche_Id");

	Node n_front = m_mesh->get<Node>(front_nodes[0]);
	int id_front = var_couche->value(n_front.id() ) ;

	for (auto n_id:m_mesh->nodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		std::vector<Node> adj_nodes = math::Utils::AdjacentNodes(m_mesh, n);
		int couche_n = var_couche->value(n.id());

		if (couche_n > 1 && couche_n < id_front ) {

			math::Vector3d Di_num({0,0,0});
			double Di_denom(0);

			for (auto n_adj : adj_nodes) {

				math::Vector3d Vj = n_adj.point() - n.point();
				math::Vector3d Cj = Vj;
				int couche = var_couche->value(n_adj.id());

				if (couche == id_front) {
					math::Vector3d Dc = ComputeDc(n_id);
					Cj = Cj + Dc;
				}

				double Cj_norm = Cj.norm();
				Di_num = Di_num + Cj_norm*Cj;
				Di_denom = Di_denom + Cj_norm;

			}

			n.setPoint(n.point() + Di_num/Di_denom );

		}
	}

}
/*------------------------------------------------------------------------*/
