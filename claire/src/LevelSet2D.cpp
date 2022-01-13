//
// Created by rochec on 13/01/2022.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/LevelSet2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/
LevelSet2D::LevelSet2D(Mesh *AMesh) {
	m_mesh = AMesh;
	m_distance = m_mesh->newVariable<double,GMDS_NODE>("distance");

}

/*------------------------------------------------------------------------*/
LevelSet2D::STATUS LevelSet2D::execute()
{
	std::cout<<"HELLO"<<std::endl;

	return LevelSet2D::SUCCESS;
}
/*------------------------------------------------------------------------*/


