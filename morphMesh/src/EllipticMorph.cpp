//
// Created by calderans on 30/08/23.
//

#include "gmds/morphMesh/EllipticMorph.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/morphMesh/FastLocalize.h"
/*----------------------------------------------------------------------------*/

using namespace gmds;
using  namespace morphmesh;
/*----------------------------------------------------------------------------*/
EllipticMorph::EllipticMorph(Mesh *AMesh): m_mesh(AMesh) {

		m_locked = m_mesh->newMark<Node>();
	   markLockedCells();

}
/*----------------------------------------------------------------------------*/
EllipticMorph::~EllipticMorph() {}
/*----------------------------------------------------------------------------*/
void EllipticMorph::execute(const std::vector<std::vector<double>>& AEllipses)
{

	std::vector<Node> nodes_int;
	std::vector<Node> nodes_ext;

	/*for(auto n : m_mesh->nodes()){
		Node node = m_mesh->get<Node>(n);
		if(m_mesh->isMarked<Node>(n, m_locked)){
			nodes_int.push_back(node);
		}
		if((node.getIDs<Region>().size() == 4 || node.getIDs<Region>().size() == 2) && !m_mesh->isMarked<Node>(n,m_locked)){
			nodes_ext.push_back(node);
		}
	}*/

	FastLocalize fl_int(m_lockedIn);
	FastLocalize fl_ext(m_lockedOut);
	FastLocalize fl_morphed(m_morphed);

	std::map<TCellID, math::Vector3d> vecs;
	std::set<TCellID> modified_nodes;

	for(int i = 0; i < AEllipses.size()-1; i++) {

		double Xmin = AEllipses[i][0];
		double Xmax = AEllipses[i+1][0];

		std::vector<Node> execution;
		std::vector<Node> interval_morphed;

		for(auto n : m_morphed){
			if (n.X() >= Xmin && n.X() < Xmax) {
				interval_morphed.push_back(n);
			}
		}


		for(auto n : m_mesh->nodes()){

			Node node = m_mesh->get<Node>(n);

			if (node.X() >= Xmin && node.X() < Xmax && !m_mesh->isMarked<Node>(n, m_locked)) {
				execution.push_back(node);
			}
		}

		double coefa_min = AEllipses[i][1];
		double coefb_min = AEllipses[i][2];
		double coefc_min = AEllipses[i][3];

		double coefa_max = AEllipses[i+1][1];
		double coefb_max = AEllipses[i+1][2];
		double coefc_max = AEllipses[i+1][3];


		for(auto n : interval_morphed){
			math::Point p = n.point();

			math::Point p_int = m_mesh->get<Node>(fl_int.find(p)).point();

			double distXmin = p.X()-Xmin;
			double distXmax = Xmax-Xmin;
			double w = distXmin/distXmax;

			double coefy_min = p.Y() >= 0 ? coefa_min : coefc_min;
			double coefy_max = p.Y() >= 0 ? coefa_max : coefc_max;

			double coefy_final = coefy_min * (1/(1 + exp(-5*( ( (1 - w ) *2)-1) ) )) + coefy_max * (1/(1 + exp(-5*( ( w *2)-1) ) ));
			double coefz_final = coefb_min * (1/(1 + exp(-5*( ( (1 - w ) *2)-1) ) )) + coefb_max * (1/(1 + exp(-5*( ( w *2)-1) ) ));

			double coefy_test = coefy_min * (3*pow(1-w,2) - 2*pow(1-w,3)) + coefy_max * (3*pow(w,2) - 2*pow(w,3));
			double coefz_test = coefb_min * (3*pow(1-w,2) - 2*pow(1-w,3)) + coefb_max * (3*pow(w,2) - 2*pow(w,3));

			//std::cout<<"coefy = "<<coefy_final<<"; coefz = "<<coefz_final<<std::endl;


			double theta = atan(n.Y() / n.Z());
			if (theta < 0) {
				theta *= -1;
			}

			double sinus = sin(theta);
			double cosinus = cos(theta);


			if (n.Y() < 0) {
				//sinus *= -1;
			}
			if (n.Z() < 0) {
				//cosinus *= -1;
			}

			std::cout<<"coefy = "<<coefy_test<<std::endl;
			//std::cout<<"coefy = "<<fabs(1 - coefy_final) * sinus <<"; coefz = "<<fabs(1 - coefz_final) * cosinus<<std::endl;

			double newY = (p.Y() * ((coefy_test * sinus) + coefz_test*(1-sinus)));
			double newZ = (p.Z() * ((coefz_test * cosinus) + coefy_test*(1-cosinus)));

			math::Vector3d vec;
			vec.setXYZ(0, newY - p.Y(), newZ - p.Z());
			vecs[n.id()] = vec;
			modified_nodes.insert(n.id());
		}

		for (auto n : execution) {
			modified_nodes.insert(n.id());

			math::Point p = n.point();

			Node nearest_morphed = m_mesh->get<Node>(fl_morphed.find(p));

			math::Point p_int = m_mesh->get<Node>(fl_int.find(p)).point();
			math::Point p_ext = m_mesh->get<Node>(fl_ext.find(p)).point();
			math::Point p_m = nearest_morphed.point();
			math::Point p_axe(p.X(), 0, 0);

			math::Point p_origin = p.distance(p_int) < p.distance(p_axe) ? p_int : p_axe;

			double dist_p = p.distance(p_axe);
			double dist_pm = p_m.distance(p_axe);

			double dist1;
			double dist2;
			double omega;


			if(dist_p < dist_pm){
				dist1 = p.distance(p_int);
				dist2 = p_m.distance(p_int);
				omega = dist1 / dist2;
			}else if(dist_p > dist_pm){
				dist1 = p.distance(p_ext);
				dist2 = p_m.distance(p_ext);
				omega = dist1 / dist2;
			}else{
				std::cout<<"error where is p ?"<<std::endl;
			}

			if(omega > 1){
				std::cout<<"error"<<std::endl;
			}

			double distXmin = p.X()-Xmin;
			double distXmax = Xmax-Xmin;
			double w = distXmin/distXmax;

			double coefy_min = p.Y() >= 0 ? coefa_min : coefc_min;
			double coefy_max = p.Y() >= 0 ? coefa_max : coefc_max;

			double coefy_test = coefy_min * (3*pow(1-w,2) - 2*pow(1-w,3)) + coefy_max * (3*pow(w,2) - 2*pow(w,3));
			double coefz_test = coefb_min * (3*pow(1-w,2) - 2*pow(1-w,3)) + coefb_max * (3*pow(w,2) - 2*pow(w,3));

			math::Vector3d vec;
			vec = vecs[nearest_morphed.id()] * omega;

			vec = omega * vec;

			vecs[n.id()] = vec;
		}

		/*for (auto n : execution) {
			modified_nodes.insert(n.id());

			math::Point p = n.point();

			Node nearest_morphed = m_mesh->get<Node>(fl_morphed.find(p));

			math::Point p_int = m_mesh->get<Node>(fl_int.find(p)).point();
			math::Point p_ext = m_mesh->get<Node>(fl_ext.find(p)).point();
			math::Point p_m = nearest_morphed.point();
			math::Point p_axe(p.X(), 0, 0);

			math::Point p_origin = p.distance(p_int) < p.distance(p_axe) ? p_int : p_axe;

			double dist_p = p.distance(p_axe);
			double dist_pm = p_m.distance(p_axe);

			double dist1;
			double dist2;
			double omega;


			if(dist_p < dist_pm){
				dist1 = p.distance(p_int);
				dist2 = p_m.distance(p_int);
				omega = dist1 / dist2;
			}else if(dist_p > dist_pm){
				dist1 = p.distance(p_ext);
				dist2 = p_m.distance(p_ext);
				omega = dist1 / dist2;
			}else{
				std::cout<<"error where is p ?"<<std::endl;
			}

			if(omega > 1){
				std::cout<<"error"<<std::endl;
			}


			double distXmin = p.X()-Xmin;
			double distXmax = Xmax-Xmin;
			double w = distXmin/distXmax;

			// math::Point p_origin = p.distance(p_int) < p.distance(p_axe) ? ;

			double coefy_min = p.Y() >= 0 ? coefa_min : coefc_min;
			double coefy_max = p.Y() >= 0 ? coefa_max : coefc_max;

			double coefy_final = coefy_min * (1/(1 + exp(-5*( ( (1 - w ) *2)-1) ) )) + coefy_max * (1/(1 + exp(-5*( ( w *2)-1) ) ));
			double coefz_final = coefb_min * (1/(1 + exp(-5*( ( (1 - w ) *2)-1) ) )) + coefb_max * (1/(1 + exp(-5*( ( w *2)-1) ) ));

			//double coefy_final = coefy_min * (1 - w) + coefy_max * (w);
			//double coefz_final = coefb_min * (1 - w) + coefb_max * (w);


			double theta = atan(n.Y() / n.Z());
			if (theta < 0) {
				theta *= -1;
			}

			double dist1 = p.distance(p_int);
			double dist2 = p_ext.distance(p_int);
			double omega = dist1 / dist2;

			// double z_new = omega*(coefa * cos(theta));
			// double y_new = omega*(coefb * sin(theta));
			double sinus = sin(theta);
			double cosinus = cos(theta);

			if (n.Y() < 0) {
				sinus *= -1;
			}
			if (n.Z() < 0) {
				cosinus *= -1;
			}

			math::Vector3d vec;
			vec = vecs[nearest_morphed.id()] * omega;

			vec = omega * vec;

			vecs[n.id()] = vec;
		}*/
	}

	for(auto n : modified_nodes){

		math::Vector3d vec(vecs[n]);
		Node node = m_mesh->get<Node>(n);
		double y = node.Y();
		double z = node.Z();

		node.setY(y+vec.Y());
		node.setZ(z+vec.Z());
	}

}
/*----------------------------------------------------------------------------*/
void EllipticMorph::markLockedCells()
{
	/*Variable<int>* locked_regions = m_mesh->newVariable<int,GMDS_REGION>("m_locked");


	std::set<Region> regions;
	std::set<TCellID> nodes;
	nodes.insert(27006);
	m_mesh->mark<Node>(27006,m_locked);
	for (int i = 0; i < 15; ++i) {
		for(auto n : nodes){
			Node node = m_mesh->get<Node>(n);

			for(const auto& r : node.get<Region>()){
				regions.emplace(r);
			}
		}
		nodes.clear();
		for(const auto& r : regions){
			locked_regions->set(r.id(),1);
			for(auto n : r.getIDs<Node>()){
				if(!m_mesh->isMarked<Node>(n,m_locked)){
					nodes.emplace(n);
					m_mesh->mark<Node>(n,m_locked);

				}
			}
		}
		regions.clear();
	}
	 for(auto n : m_mesh->nodes()){
  Node node = m_mesh->get<Node>(n);
  math::Point p(node.point().X(),0,0);
  if(node.point().distance(p) < 0.2){
	  m_mesh->mark<Node>(n,m_locked);
  }
}
	 */

	Variable<int>* locked_in = m_mesh->getVariable<int,GMDS_REGION>("Hors_Groupe_3D");
	Variable<int>* locked_out = m_mesh->getVariable<int,GMDS_REGION>("ext");
	Variable<int>* morphed = m_mesh->getVariable<int,GMDS_REGION>("couche");

	std::set<TCellID> set_li;
	std::set<TCellID> set_lo;
	std::set<TCellID> set_mo;

	for(auto r : m_mesh->regions()){
		Region region = m_mesh->get<Region>(r);
		if(locked_in->value(r) == 1){
			for(auto n : region.getIDs<Node>()){
				set_li.insert(n);
			}
		}else if(locked_out->value(r) == 1){
			for(auto n : region.getIDs<Node>()){
				set_lo.insert(n);
			}
		}else if(morphed->value(r) == 1){
			for(auto n : region.getIDs<Node>()){
				set_mo.insert(n);
			}
		}
	}

	for(auto n : set_li){
		Node node = m_mesh->get<Node>(n);
		m_lockedIn.push_back(node);
		m_mesh->mark<Node>(n,m_locked);
	}
	for(auto n : set_lo){
		Node node = m_mesh->get<Node>(n);
		m_lockedOut.push_back(node);
		m_mesh->mark<Node>(n,m_locked);
	}
	for(auto n : set_mo){
		Node node = m_mesh->get<Node>(n);
		m_morphed.push_back(node);
		m_mesh->mark<Node>(n,m_locked);
	}
}