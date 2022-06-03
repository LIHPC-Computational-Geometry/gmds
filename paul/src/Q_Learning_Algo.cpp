//
// Created by bourmaudp on 30/05/22.
//
/*---------------------------------------------------------*/
#include "gmds/paul/Q_Learning_Algo.h"
/*---------------------------------------------------------*/

using namespace gmds;

void gmds::executeTrainQlearning(Environment environmentInit){
	double epsilon = 0.9;
	double discount_factor= 0.9;
	double learning_rate= 0.9;

	Politique politique(&environmentInit);
	politique.initQTable();

	Tools toolInit(&environmentInit.g_grid);


	for (int i=0; i<=4;i++){
		std::cout <<"Dans le for : "<<i<<std::endl;



		Mesh copyGridMesh(environmentInit.g_grid.m_mesh.getModel());
		GridBuilderAround copyGrid(&copyGridMesh,&environmentInit.g_grid.meshTarget,2);
		//copyGrid.executeGrid2D(5);
		toolInit.cloneMesh(environmentInit.g_grid,copyGrid);
		//Tools toolsCopy(&copyGrid);
		Actions actionQLearning(&copyGrid);
		Environment environment(&copyGrid,&environmentInit.g_grid.meshTarget,&actionQLearning);
		/*
		auto facesAllCopy = copyGrid.m_mesh.faces();
		for (auto f : facesAllCopy){
			std::cout<<"La Face : "<<copyGrid.m_mesh.get<Face>(f)<<std::endl;
		}*/





		while(environment.globalIoU()<0.9){

			std::cout<<"Valeur env Global IoU : \n"<<environment.globalIoU()<<std::endl;

			Face faceSelected = environment.faceSelect();

			std::cout<<"Face Select : "<<faceSelected<< " avec valeur activate "<<copyGrid.getActivate(faceSelected)<<" et vol Frac "
			          <<copyGrid.m_mesh.getVariable<double,GMDS_FACE>("volFrac")->value(faceSelected.id())<<std::endl;
			double localIoU = environment.localIoU(faceSelected);
			int intervalIoU = politique.getInterval(localIoU);
			int actionIndex = politique.getNextAction(intervalIoU);
			double oldQValue = politique.Q_Table[intervalIoU][actionIndex];
			double maxQValue = politique.maxQValue(intervalIoU);

			environment.executeAction(faceSelected,actionIndex);
			double reward = environment.reward(faceSelected);

			double temporal_difference = reward + (0.9 * maxQValue) - oldQValue;
			double newQValue = oldQValue + (0.9 * temporal_difference);

			politique.updateQTable(intervalIoU,actionIndex,newQValue);
		}

	}
	std::cout<<"Training complete !"<<std::endl;
	std::cout<<"LA Q TABLE "<<std::endl;
	for (auto l : politique.Q_Table ){
		std::cout<<"une ligne de la qtable : "<<std::endl;
		for (auto d : l ){
			std::cout<<d<<std::endl;
		}
	}
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