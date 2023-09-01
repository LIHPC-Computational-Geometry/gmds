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


	for(auto n : m_mesh->nodes()){
		Node node = m_mesh->get<Node>(n);
		if(m_mesh->isMarked<Node>(n, m_locked)){
			nodes_int.push_back(node);
		}
		if((node.getIDs<Region>().size() == 4 || node.getIDs<Region>().size() == 2) && !m_mesh->isMarked<Node>(n,m_locked)){
			nodes_ext.push_back(node);
		}
	}
	std::map<TCellID, math::Vector3d> vecs;
	std::set<TCellID> modified_nodes;

	for(int i = 0; i < AEllipses.size()-1; i++) {

		double Xmin = AEllipses[i][0];
		double Xmax = AEllipses[i+1][0];

		std::vector<Node> execution;

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

		FastLocalize fl_int(nodes_int);
		FastLocalize fl_ext(nodes_ext);


		for (auto n : execution) {
			modified_nodes.insert(n.id());

			math::Point p = n.point();

			Node n_ext = m_mesh->get<Node>(fl_ext.find(p));
			math::Point p_int = m_mesh->get<Node>(fl_int.find(p)).point();
			math::Point p_ext = n_ext.point();
			math::Point p_axe(p.X(), 0, 0);

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
			vec.setXYZ(0, (coefy_final * sinus) - p_ext.Y(), (coefz_final * cosinus) - p_ext.Z());

			vec = omega * vec;

			vecs[n.id()] = vec;
		}
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
	Variable<int>* locked_regions = m_mesh->newVariable<int,GMDS_REGION>("m_locked");


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
}