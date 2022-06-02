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
/*
	for (auto l : Q_Table){
		std::cout<<"une ligne de la qtable : "<<std::endl;
		for (auto d : l ){
			std::cout<<d<<std::endl;
		}
	}*/
}

std::vector<std::vector<double>> Politique::getQTable()
{
	return Q_Table;
}

int Politique::getNextAction(int intervalIoU)
{
	auto theQTable = getQTable();
	int actionSelect;
	double maxQvalue;
	bool firstIte = true;
	float epsilon = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);;
	std::cout<<"Epsilon : "<< epsilon<<std::endl;
	if(epsilon < 0.9)
	{
		std::cout<<"Dans le if"<<std::endl;
		for (auto i : theQTable[intervalIoU]) {
			if (firstIte) {
				actionSelect = i;
				maxQvalue = theQTable[intervalIoU][i];
				firstIte = false;
			}
			else {
				if (maxQvalue < theQTable[intervalIoU][i]) {
					actionSelect = i;
					maxQvalue = theQTable[intervalIoU][i];
				}
			}
		}
		std::cout<<"Action Select : "<<actionSelect<<std::endl;
		return actionSelect;
	}
	else{
		int randomAction = rand() % 5;
		std::cout<<"Action alea "<< randomAction<<std::endl;
		return randomAction;

	}

}

double Politique::maxQValue(int intervalIoU)
{
	auto theQTable = getQTable();
	double maxQValue;
	bool firstIte = true;
	for (auto i : theQTable[intervalIoU]) {
		if (firstIte) {
			maxQValue = theQTable[intervalIoU][i];
			firstIte = false;
		}
		else {
			if (maxQValue < theQTable[intervalIoU][i]) {
				maxQValue = theQTable[intervalIoU][i];
			}
		}
	}
	return maxQValue;

}

void Politique::updateQTable(int intervalIoU, int  actionSelect, double newQValue)
{
	Q_Table[intervalIoU][actionSelect]=newQValue;

}

int Politique::getInterval(double localIoU)
{
	int interval;
	if(localIoU==0){
		interval = 0;
		return interval;
	}
	else if(localIoU>0 && localIoU<=0.3){
		interval = 1;
		return interval;
	}
	else if(localIoU>0.3 && localIoU<=0.6){
		interval = 2;
		return interval;
	}
	else if(localIoU>0.6 && localIoU<1){
		interval = 3;
		return interval;
	}
	else if(localIoU==1){
		interval = 4;
		return interval;
	}
}