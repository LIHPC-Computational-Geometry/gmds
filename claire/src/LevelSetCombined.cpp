//
// Created by rochec on 19/01/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/LevelSetCombined.h>
#include <gmds/claire/LevelSetExtended.h>
//#include <limits>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

LevelSetCombined::LevelSetCombined(Mesh *AMesh, int AmarkFrontNodesInt, int AmarkFrontNodesOut,
                                   Variable<double>* Adistance, Variable<double>* Adistance_Int, Variable<double>* Adistance_Out) {
	m_mesh = AMesh;
	m_markFrontNodesInt = AmarkFrontNodesInt;
	m_markFrontNodesOut = AmarkFrontNodesOut;
	m_distance = Adistance;
	m_distance_Int = Adistance_Int;
	m_distance_Out = Adistance_Out;

}





/*------------------------------------------------------------------------*/
LevelSetCombined::STATUS
LevelSetCombined::execute()
{
	initialisationDistancesInt();
	initialisationDistancesOut();

	combineDistanceFields();

	return LevelSetCombined::SUCCESS;
}
/*------------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
void
LevelSetCombined::setValue(TCellID n_id, double v0){
	m_distance->value(n_id) = v0 ;
};
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
void
LevelSetCombined::combineDistanceFields() {
	double distInt;
	double distOut;
	for (auto id:m_mesh->nodes()){
		distInt = m_distance_Int->value(id);
		distOut = m_distance_Out->value(id);
		if (distInt+distOut == 0) {
			setValue(id, -1 );
		}
		else {
			setValue(id, distInt / (distInt + distOut));     // Attention, pas de barrière si distInt+dOut =0 (c'est à dire si les deux fronts ont un noeud commun)
		}
	}
};
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
void
LevelSetCombined::initialisationDistancesInt() {
	LevelSetExtended lsInt(m_mesh, m_markFrontNodesInt, m_distance_Int);
	lsInt.execute();
	double distInt;
	for (auto id:m_mesh->nodes()){
		lsInt.getValue(id, distInt);
		m_distance_Int->set(id, distInt);
	}
	//m_mesh->deleteVariable(GMDS_NODE, "GMDS_Distance");
};
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
void
LevelSetCombined::initialisationDistancesOut() {
	LevelSetExtended lsOut(m_mesh, m_markFrontNodesOut, m_distance_Out);
	lsOut.execute();
	double distOut;
	for (auto id:m_mesh->nodes()){
		lsOut.getValue(id, distOut);
		m_distance_Out->set(id, distOut);
	}
	//m_mesh->deleteVariable(GMDS_NODE, "GMDS_Distance");
};
/*-------------------------------------------------------------------*/