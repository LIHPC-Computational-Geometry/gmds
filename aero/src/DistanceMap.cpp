//
// Created by rochec on 14/01/2022.
//
/*------------------------------------------------------------------------*/
#include <gmds/aero/DistanceMap.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

DistanceMap::DistanceMap() {

}



/*-------------------------------------------------------------------*/
void DistanceMap::add(double v0, TCellID n_id){
	m_map[v0].push_back(n_id);
}
/*-------------------------------------------------------------------*/
void DistanceMap::remove(double v0, TCellID n_id){

	// Cette manière n'a pas l'air de fonctionner
	//std::remove(m_map[v0].begin(),m_map[v0].end(),n_id);

	// Option non optimisée, il faudrait un while pour arrêter le parcours
	// lorsque l'indice est trouvé
	// ---
	int indice(-1);
	int increm(0);
	/*
	for (auto i:m_map[v0]){
		if (m_map[v0][i]==n_id){
			indice = i;
		}
	} */

	while (increm < m_map[v0].size() && indice == -1) {
		if (m_map[v0][increm]==n_id){
			indice = increm;
		}
		increm += 1;
	}
	// ---

	m_map[v0].erase(m_map[v0].begin()+indice);

	// Si le premier vecteur est nul, il faut alors retirer le couple
	// (clé, vecteur d'id) de la map.
	if(m_map[v0].empty()) {
		m_map.erase(v0);
	}
}
/*-------------------------------------------------------------------*/
void DistanceMap::getAndRemoveFirst(double &AMinDist, TCellID &AMinId){
	AMinDist = m_map.begin()->first;
	if (!m_map[AMinDist].empty()) {
		AMinId = m_map[AMinDist][0];
	}

	remove(AMinDist, AMinId);

}
/*-------------------------------------------------------------------*/
void DistanceMap::update(double v0_old, double v0_new, TCellID n_id){

	// L'ancien élèment de la map est enlevé
	remove(v0_old, n_id);
	// Le nouvel élèment dans la map est ajouté
	add(v0_new, n_id);

}
/*-------------------------------------------------------------------*/
bool DistanceMap::check(){
	bool test(true);
	// Ce parcours n'est pas optimisé. Il vaudrait mieux ajouter un
	// while true, comme ça, à la détection du premier pb (si il
	// y en a un), on arrête d'itérer.
	for(auto const& val:m_map){
		for(auto id:val.second) {
			int compteur = 0;
			//std::cout << "Noeud traité :" << id << std::endl;
			for(auto const& val_2:m_map){
				for(auto id_2:val_2.second) {
					//std::cout << id_2 << std::endl;
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
}
/*-------------------------------------------------------------------*/
void DistanceMap::getNbrIds(double v0, int &nbr){
	// Si la clé existe
	if (m_map.find(v0) != m_map.end()) {
		nbr = m_map[v0].size();
	}
	else{
		nbr = 0;
	}
}
/*-------------------------------------------------------------------*/
void DistanceMap::getNbrIdsTotal(int &nbr){
   nbr = 0;
	int nbr_local = 0;
	double v0;
	for (auto const& paire:m_map){
		v0 = paire.first ;
		getNbrIds(v0, nbr_local);
		nbr += nbr_local;
	}


}
/*-------------------------------------------------------------------*/
bool DistanceMap::isEmpty(){
	return m_map.empty() ;
}
/*-------------------------------------------------------------------*/



namespace  gmds{
	std::ostream& operator<<(std::ostream& AOS, const DistanceMap& ADM){
	   for(auto const& val:ADM.m_map){
		   AOS<<val.first<<": ";
		   for(auto id:val.second)
			   AOS<<id<<" ";
		   AOS<<"\n";
	   }
	   return AOS;}
}