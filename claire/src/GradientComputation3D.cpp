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
	math::Vector3d Gradient(0,0,0) ;
	Region reg = m_mesh->get<Region>(region_id);
	std::vector<TCellID> region_nodes_ids = reg.getIDs<Node>();

	double di = m_distance->value(region_nodes_ids[0]) ;
	double dj = m_distance->value(region_nodes_ids[1]) ;
	double dk = m_distance->value(region_nodes_ids[2]) ;
	double dh = m_distance->value(region_nodes_ids[3]) ;
	double Vt = reg.volume();

	Node ni = m_mesh->get<Node>(region_nodes_ids[0]);
	Node nj = m_mesh->get<Node>(region_nodes_ids[1]);
	Node nk = m_mesh->get<Node>(region_nodes_ids[2]);
	Node nh = m_mesh->get<Node>(region_nodes_ids[3]);

	math::Point pi = ni.point();
	math::Point pj = nj.point();
	math::Point pk = nk.point();
	math::Point ph = nh.point();

	math::Vector3d vki = pi-pk ;
	math::Vector3d vkh = ph-pk ;
	math::Vector3d vhi = pi-ph ;
	math::Vector3d vhj = pj-ph ;
	math::Vector3d vik = pk-pi ;
	math::Vector3d vij = pj-pi ;

	math::Vector3d Vec_1 = vki.cross(vkh);
	math::Vector3d Vec_2 = vhi.cross(vhj);
	math::Vector3d Vec_3 = vik.cross(vij);

	Gradient = (Vec_1*(dj-di)/(2.0*Vt) + Vec_2*(dk-di)/(2.0*Vt) + Vec_3*(dh-di)/(2.0*Vt)) ;


	std::cout << "----------------------------------" << std::endl ;
	std::cout << "ID région traitée : " << region_id << std::endl;
	//std::cout << "Produit vectoriel 1 :" << Vec_1 << std::endl;
	//std::cout << "Produit vectoriel 2 :" << Vec_2 << std::endl;
	//std::cout << "Produit vectoriel 3 :" << Vec_3 << std::endl;
	std::cout << "Gradient :" << Gradient << std::endl;


	return Gradient;
}
/*------------------------------------------------------------------------*/
