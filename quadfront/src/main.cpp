#include <iostream>
#include "gmds/quadfront/Quadfront.h"

/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
namespace quadfront {

int
main()
{
	std::cout << "============== Début Maillage par avancée de front ==============" << std::endl;

	gmds::quadfront::Quadfront quad = gmds::quadfront::Quadfront("/home/pagea/Documents/Travail/data/t1-2.vtk");

	for (auto e : quad.m_)

	// On parcourt les nodes à partir du m_edgeBoundary
	// Si l'angle < 3*M_PI/4
	//		-> on change à 1
	// Sinon
	//		-> on cherche le side edge
	// Si le side edge est inférieur à M_PI/6
	//		-> On séléctione le side edge comme coté du quad
	// Sinon

	// Changer la fonction angle pour les cas limites

	std::cout << "============== Fin Maillage par avancée de front ==============" << std::endl;
}
}
}