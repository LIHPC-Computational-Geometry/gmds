//
// Created by rochec on 19/01/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/aero/LevelSetCombined.h>
#include <gmds/aero/LevelSetExtended.h>
//#include <limits>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

LevelSetCombined::LevelSetCombined(Mesh *AMesh, TInt AmarkFrontNodesInt, TInt AmarkFrontNodesOut,
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
	double distInt;
	double distOut;

	// Calcul du level set de l'intérieur vers l'extérieur
	LevelSetExtended lsInt(m_mesh, m_markFrontNodesInt, m_distance_Int);
	lsInt.execute();

	// Calcul du level set de l'extérieur vers l'intérieur
	LevelSetExtended lsOut(m_mesh, m_markFrontNodesOut, m_distance_Out);
	lsOut.execute();

	for (auto id:m_mesh->nodes()){

		lsInt.getValue(id, distInt);
		m_distance_Int->set(id, distInt);

		lsOut.getValue(id, distOut);
		m_distance_Out->set(id, distOut);

		if (distInt+distOut == 0) {
			setValue(id, -1 );
		}
		else {
			setValue(id, distInt / (distInt + distOut));     // Attention, pas de barrière si distInt+dOut =0 (c'est à dire si les deux fronts ont un noeud commun)
		}
	}

	return LevelSetCombined::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void
LevelSetCombined::setValue(TCellID n_id, double v0){
	m_distance->value(n_id) = v0 ;
}
/*-------------------------------------------------------------------*/
