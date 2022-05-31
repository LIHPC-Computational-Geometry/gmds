//
// Created by bourmaudp on 31/05/22.
//

#include "gmds/paul/Politique.h"

using namespace gmds;
Politique::Politique(Environment *AEnv)
   :env(*AEnv){;}

void Politique::initQTable()
{
	for (int i = 0; i <5;i++ ){
		std::vector<double> ligne;
		for (int j = 0; j<5;j++){
			ligne.push_back(0);
		}
		Q_Table.push_back(ligne);
	}

	for (auto l : Q_Table){
		std::cout<<"une ligne de la qtable : "<<std::endl;
		for (auto d : l ){
			std::cout<<d<<std::endl;
		}
	}
}

std::vector<std::vector<double>> Politique::getQTable()
{
	return Q_Table;
}

