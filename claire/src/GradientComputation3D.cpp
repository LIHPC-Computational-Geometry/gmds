//
// Created by rochec on 21/01/2022.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/GradientComputation3D.h>
#include <limits>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/


GradientComputation3D::GradientComputation3D(Mesh *AMesh, Variable<double>* Adistance) {
	m_mesh = AMesh;
	m_distance = Adistance;
	m_gradient3D = m_mesh->newVariable<math::Vector3d,GMDS_REGION>("gradient_3D");
}


/*------------------------------------------------------------------------*/
GradientComputation3D::STATUS GradientComputation3D::execute()
{

	for(auto region_id:m_mesh->regions()) {
		math::Vector3d Gradient = computeGradientOnSimpleRegion(region_id);
		m_gradient3D->set(region_id, Gradient) ;
	}

	return GradientComputation3D::SUCCESS;
}
/*------------------------------------------------------------------------*/



/*------------------------------------------------------------------------*/
math::Vector3d GradientComputation3D::computeGradientOnSimpleRegion(TCellID region_id){
	math::Vector3d Gradient ;
	Region reg = m_mesh->get<Region>(region_id);
	std::vector<TCellID> region_nodes_ids = reg.getIDs<Node>();

	double dist_n0 = m_distance->value(region_nodes_ids[0]) ;
	double dist_n1 = m_distance->value(region_nodes_ids[1]) ;
	double dist_n2 = m_distance->value(region_nodes_ids[2]) ;
	double dist_n3 = m_distance->value(region_nodes_ids[3]) ;
	double Vt = reg.volume();

	Node n0 = m_mesh->get<Node>(region_nodes_ids[0]);
	Node n1 = m_mesh->get<Node>(region_nodes_ids[1]);
	Node n2 = m_mesh->get<Node>(region_nodes_ids[2]);
	Node n3 = m_mesh->get<Node>(region_nodes_ids[3]);

	math::Point p0 = n0.point();
	math::Point p1 = n1.point();
	math::Point p2 = n2.point();
	math::Point p3 = n3.point();

	math::Vector3d Vec_loc_1;
	math::Vector3d Vec_loc_2;

	Vec_loc_1 = p0-p2 ;
	Vec_loc_2 = p3-p2 ;
	math::Vector3d Vec_1 = Vec_loc_2.cross(Vec_loc_1);

	Vec_loc_1 = p0-p3 ;
	Vec_loc_2 = p1-p3 ;
	math::Vector3d Vec_2 = Vec_loc_2.cross(Vec_loc_1);

	Vec_loc_1 = p2-p0 ;
	Vec_loc_2 = p1-p0 ;
	math::Vector3d Vec_3 = Vec_loc_2.cross(Vec_loc_1);

	Gradient = Vec_1*(dist_n1-dist_n0)/(2.0*Vt) + Vec_2*(dist_n2-dist_n0)/(2.0*Vt) + Vec_3*(dist_n3-dist_n0)/(2.0*Vt) ;

	/*
	std::cout << "----------------------------------" << std::endl ;
	std::cout << "ID région traitée : " << region_id << std::endl;
	std::cout << "Produit vectoriel 1 :" << Vec_1 << std::endl;
	std::cout << "Produit vectoriel 2 :" << Vec_2 << std::endl;
	std::cout << "Produit vectoriel 3 :" << Vec_3 << std::endl;
	std::cout << "Gradient :" << Gradient << std::endl;
	*/

	return Gradient;
}
/*------------------------------------------------------------------------*/
