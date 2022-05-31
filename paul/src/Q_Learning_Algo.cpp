//
// Created by bourmaudp on 30/05/22.
//
/*---------------------------------------------------------*/
#include "gmds/paul/Q_Learning_Algo.h"
/*---------------------------------------------------------*/

using namespace gmds;

void gmds::executeAlgoQlearning(Environment environment){
	double epsilon = 0.9;
	double discount_factor= 0.9;
	double learning_rate= 0.9;


	for (int i=0; i<=10;i++){
		while(environment.globalIoU()<0.9){

			Face faceSelected = environment.faceSelect();
			auto allFaces = environment.g_grid.m_mesh.faces();
			std::vector<Face> oldFaces;
			for(auto f : allFaces){
				oldFaces.push_back(environment.g_grid.m_mesh.get<Face>(f));
			}
			auto allNodes = environment.g_grid.m_mesh.nodes();
			std::vector<Face> oldNodes;
			for(auto n : allNodes){
				oldNodes.push_back(environment.g_grid.m_mesh.get<Face>(n));
			}





		}
	}
	std::cout<<"Training complete !"<<std::endl;

}

int gmds::getNextAction(){
	int ok = 1;
return ok;
}

void gmds::initQTable(){
	std::vector<std::vector<double>> Q_Table;
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