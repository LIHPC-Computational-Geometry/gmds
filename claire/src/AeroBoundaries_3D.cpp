//
// Created by rochec on 24/03/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AeroBoundaries_3D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AeroBoundaries_3D::AeroBoundaries_3D(Mesh *AMesh) :
  AbstractAeroBoundaries(AMesh)
{

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroBoundaries_3D::MarkBoundariesNodes(){

	//Get the boundary node ids
	for (auto f_id:m_mesh->faces())
	{
		Face f = m_mesh->get<Face>(f_id);
		if (f.get<Region>().size() == 1)
		{
			std::vector<Node> face_nodes = f.get<Node>() ;
			for (auto n:face_nodes){
				m_bnd_nodes_ids.push_back(n.id());
			}
		}
	}

	// Initialise une marque sur les noeuds du bord
	for (auto n_id:m_bnd_nodes_ids){
		Node n = m_mesh->get<Node>(n_id);
		m_mesh->mark(n, m_markBoundaryNodes);
	}

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroBoundaries_3D::WhichColorIsAmont(){

	// Calcul des boites englobantes
	std::vector<double> x_min(m_nbrBords,0.0);
	std::vector<double> y_min(m_nbrBords,0.0);
	std::vector<double> z_min(m_nbrBords,0.0);
	std::vector<double> x_max(m_nbrBords,0.0);
	std::vector<double> y_max(m_nbrBords,0.0);
	std::vector<double> z_max(m_nbrBords,0.0);

	for (auto n_id:m_bnd_nodes_ids ) {
		Node n = m_mesh->get<Node>(n_id);
		math::Point p = n.point();
		// A quel bord appartient le noeud n_id
		int couleur = m_var_color_bords->value(n_id);
		if(p.X() < x_min[couleur-1]){
			x_min[couleur-1] = p.X();
		}
		if(p.X() > x_max[couleur-1]){
			x_max[couleur-1] = p.X();
		}
		if(p.Y() < y_min[couleur-1]){
			y_min[couleur-1] = p.Y();
		}
		if(p.Y() > y_max[couleur-1]){
			y_max[couleur-1] = p.Y();
		}
		if(p.Z() < z_min[couleur-1]){
			z_min[couleur-1] = p.Z();
		}
		if(p.Z() > z_max[couleur-1]){
			z_max[couleur-1] = p.Z();
		}
	}

	// On recherche la plus grosse boite englobante
	for(int i=2;i<=m_nbrBords;i++){
		// On compte le nombre de boites que la couleur i englobe
		int nbr_boite_englobe = 0;
		for (int j=1;j<=m_nbrBords;j++){
			if (i != j){
				if ( (x_min[i-1] < x_min[j-1]) &&
				    (y_min[i-1] < y_min[j-1]) &&
				    (z_min[i-1] < z_min[j-1]) &&
				    (x_max[j-1] < x_max[i-1]) &&
				    (y_max[j-1] < y_max[i-1]) &&
				    (z_max[j-1] < z_max[i-1]) ) {
					nbr_boite_englobe += 1;
				}
			}
		}
		if(nbr_boite_englobe == m_nbrBords-1){
			m_color_Amont = i;
		}
	}



}
/*------------------------------------------------------------------------*/
