//
// Created by rochec on 13/01/2022.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/LevelSet2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

// OPTION 1
LevelSet2D::LevelSet2D(Mesh *AMesh, std::vector<TCellID> Afront_nodes_Ids) {
	m_mesh = AMesh;
	m_front_nodes_Ids = Afront_nodes_Ids;
	m_distance = m_mesh->newVariable<double,GMDS_NODE>("distance");

}
/*
// OPTION 2
LevelSet2D::LevelSet2D(Mesh *AMesh, Variable<int> *Afront_nodes_Ids) {
	m_mesh = AMesh;
	m_front_nodes_Ids = Afront_nodes_Ids;
	m_distance = m_mesh->newVariable<double,GMDS_NODE>("distance");

}
 */

/*------------------------------------------------------------------------*/
LevelSet2D::STATUS LevelSet2D::execute()
{
	std::cout<<"HELLO"<<std::endl;

	return LevelSet2D::SUCCESS;
}
/*------------------------------------------------------------------------*/


