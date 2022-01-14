//
// Created by rochec on 14/01/2022.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/DistanceMap.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

DistanceMap::DistanceMap() {

}



/*-------------------------------------------------------------------*/
void DistanceMap::add(double v0, TCellID n_id){
	m_map[v0].push_back(n_id);
};
/*-------------------------------------------------------------------*/




/*-------------------------------------------------------------------*/
void DistanceMap::remove(double v0, TCellID n_id){

	std::remove(m_map[v0].begin(),m_map[v0].end(),n_id);

	// Si le premier vecteur est nul, il faut alors retirer le couple
	// (clé, vecteur d'id) de la map.
	if(m_map[v0].empty()) {
		m_map.erase(v0);
	}
};
/*-------------------------------------------------------------------*/





/*-------------------------------------------------------------------*/
void DistanceMap::getAndRemoveFirst(double &AMinDist, TCellID &AMinId){
	AMinDist = m_map.begin()->first;
	if (!m_map[AMinDist].empty()) {
		AMinId = m_map[AMinDist][0];
	}

	remove(AMinDist, AMinId);

};
/*-------------------------------------------------------------------*/





/*-------------------------------------------------------------------*/
void DistanceMap::update(double v0_old, double v0_new, TCellID n_id){

	// L'ancien élèment de la map est enlevé
	remove(v0_old, n_id);
	// Le nouvel élèment dans la map est ajouté
	add(v0_new, n_id);

};
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
bool DistanceMap::check(){
	bool test(true);
	// Ce parcours n'est pas optimisé. Il vaudrait mieux ajouter un
	// while true, comme ça, à la détection du premier pb (si il
	// y en a un), on arrête d'itérer.
	for(auto val:m_map){
		for(auto id:val.second) {
			int compteur = 0;
			for(auto val_2:m_map){
				for(auto id_2:val.second) {
					if (id==id_2){
						compteur +=1;
					}
				}
			}
			if (compteur!=1){
				test = false;
			}
		}
	}

	return test;
};
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
void DistanceMap::getNbrIds(double v0, int &nbr){
	// Si la clé existe
	if (m_map.find(v0) != m_map.end()) {
		nbr = m_map[v0].size();
	}
	else{
		nbr = 0;
	}
};
/*-------------------------------------------------------------------*/









namespace  gmds{
	std::ostream& operator<<(std::ostream& AOS, const DistanceMap& ADM){
	   for(auto val:ADM.m_map){
		   AOS<<val.first<<": ";
		   for(auto id:val.second)
			   AOS<<id<<" ";
		   AOS<<"\n";
	   }
	   return AOS;}
}