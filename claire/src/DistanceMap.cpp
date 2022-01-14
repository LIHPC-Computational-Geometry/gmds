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


};
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
bool DistanceMap::check(){

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