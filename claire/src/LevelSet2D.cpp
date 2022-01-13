//
// Created by rochec on 13/01/2022.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/LevelSet2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

/*
// OPTION 1
LevelSet2D::LevelSet2D(Mesh *AMesh, std::vector<TCellID> Afront_nodes_Ids) {
	m_mesh = AMesh;
	m_front_nodes_Ids = Afront_nodes_Ids;
	m_distance = m_mesh->newVariable<double,GMDS_NODE>("distance");
}
*/
// OPTION 2
LevelSet2D::LevelSet2D(Mesh *AMesh, Variable<int> *Afront_nodes_Ids) {
	m_mesh = AMesh;
	m_front_nodes_Ids = Afront_nodes_Ids;
	m_distance = m_mesh->newVariable<double,GMDS_NODE>("distance");

}

/*------------------------------------------------------------------------*/
LevelSet2D::STATUS LevelSet2D::execute()
{
	// OPTION 1
	// Initialisation du champs de distances
	// InitialisationDistances();

	// Mise en place du tas
	InitialisationNodeMinHeap() ;
	std::cout << "------- TEST --------" << std::endl ;
	std::cout << "Nombre de noeuds total :" << m_mesh->getNbNodes() << std::endl;
	std::cout << m_NodeMinHeap[0].size() << std::endl ;
	std::cout << m_NodeMinHeap[10^6].size() << std::endl ;

	return LevelSet2D::SUCCESS;
}
/*------------------------------------------------------------------------*/



// OPTION 2
/*------------------------------------------------------------------------*/
/*
// Fonction InitialisationDistances : Permet d'itialiser les valeurs du champs de distances.
// Je décide ici de mettre à -1 les distances à l'initialisation, et à 0 les distances des noeuds sur le front regardé.
// En entrée :
// En sortie :
void LevelSet2D::InitialisationDistances() {

	for(auto id:m_mesh->nodes()){
		m_distance->value(id)=10^6;
	}

	for(auto id:m_front_nodes_Ids){
		m_distance->value(id)=0;
	}

	// REMARQUE : Si on avait choisi l'option 2, alors on n'aurait pas eu à réitérer
	// sur les id des noeuds dans le vecteur front, on aurait pu tout initialiser
	// en un coup.

}
*/
/*------------------------------------------------------------------------*/





/*------------------------------------------------------------------------*/
// Fonction InitialisationDistances : Permet d'itialiser le tas. En fait, ici, je préférerais qu'on ait choisi
// l'option 2.
// En entrée :
// En sortie :
void LevelSet2D::InitialisationNodeMinHeap() {

	std::vector<TCellID> front_Ids;
	std::vector<TCellID> other_Ids;

	for(auto id:m_mesh->nodes()){
		if(m_front_nodes_Ids->value(id)==1){
			front_Ids.push_back(id);
		}
		else{
			other_Ids.push_back(id);
		}
	}

	m_NodeMinHeap[0]=front_Ids;
	m_NodeMinHeap[10^6]=other_Ids; // Valeur arbitraire pour l'instant, juste pour tester en attendant de mettre l'infini

}
/*------------------------------------------------------------------------*/