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

	Politique *politique=new Politique(&environmentInit);
	politique->initQTable();

	Tools toolInit(&environmentInit.g_grid);


	for (int i=0; i<=100;i++){
		std::cout <<"Dans le for : "<<i<<std::endl;


		std::cout<<"Copy mesh"<<std::endl;
		Mesh *copyGridMesh= new Mesh(environmentInit.g_grid.m_mesh.getModel());
		std::cout<<"Copy grid"<<std::endl;

		GridBuilderAround *copyGrid=new GridBuilderAround(copyGridMesh,&environmentInit.g_grid.meshTarget,2);
		//copyGrid.executeGrid2D(5);
		std::cout<<"Clone"<<std::endl;

		toolInit.cloneMesh(environmentInit.g_grid,*copyGrid);
		//Tools toolsCopy(&copyGrid);
		std::cout<<"Action"<<std::endl;

		Actions *actionQLearning=new Actions(copyGrid);
		Tools *toolsQlearning=new Tools(copyGrid);
		std::cout<<"Env"<<std::endl;
		Environment *environment= new Environment(copyGrid,&environmentInit.g_grid.meshTarget,actionQLearning,toolsQlearning);
		/*
		auto facesAllCopy = copyGrid.m_mesh.faces();
		for (auto f : facesAllCopy){
			std::cout<<"La Face : "<<copyGrid.m_mesh.get<Face>(f)<<std::endl;
		}*/

		//Face selectFaceInit = environmentInit.faceSelect();
		//Face selectFace = environment.faceSelect();

		//environmentInit.executeAction(selectFaceInit,1);
		//environment.executeAction(selectFace,1);

		int countIte=0;
		while(environment->globalIoU()<0.96 && countIte<30){

			environment->calcVolFrac();
			//std::cout<<"Valeur env Global IoU : \n"<<environment->globalIoU()<<std::endl;
			//std::cout<<"Select Face"<<std::endl;
			Face faceSelected = environment->faceSelect();

			//std::cout<<"Face Select : "<<faceSelected<< " avec valeur activate "<<copyGrid->getActivate(faceSelected)<<" et vol Frac "
			          //<<copyGrid->m_mesh.getVariable<double,GMDS_FACE>("volFrac")->value(faceSelected.id())<<std::endl;
			//std::cout<<"Calcul IoU"<<std::endl;

			double localIoU = environment->localIoU(faceSelected);
			//std::cout<<"Calcul intervalIoU"<<std::endl;
			int intervalIoU = politique->getInterval(localIoU);
			//std::cout<<"Calcul actionIndex"<<std::endl;
			int actionIndex = politique->getNextAction(intervalIoU);
			//std::cout<<"oldQValue"<<std::endl;
			double oldQValue = politique->Q_Table[intervalIoU][actionIndex];
			//std::cout<<"maxQValue"<<std::endl;
			double maxQValue = politique->maxQValue(intervalIoU);

			//std::cout<<"executeAction"<<std::endl;
			if(localIoU==0){
				std::cout<<"Action delete"<<std::endl;
				std::cout<<"local IoU : "<<localIoU<<std::endl;
				environment->executeAction(faceSelected, 0);
				actionIndex=0;
			}
			else {
				environment->executeAction(faceSelected, actionIndex);
			}
			//std::cout<<"Calcul reward"<<std::endl;
			double reward = environment->reward(faceSelected);
			//std::cout<<"Calcul TD"<<std::endl;
			double temporal_difference = reward + (0.9 * maxQValue) - oldQValue;
			//std::cout<<"Calcul newQvalue"<<std::endl;
			double newQValue = oldQValue + (0.9 * temporal_difference);
			//std::cout << "#nbface " << faceSelected.nbNodes() << std::endl;
			//std::cout<<"Calcul updateQTable"<<std::endl;
			politique->updateQTable(intervalIoU,actionIndex,newQValue);
			//std::cout<<"Count Ite : "<<countIte<<std::endl;

			countIte+=1;
		}

		std::cout << "Fin du for, Bientot le nouveau"<<std::endl;

		gmds::IGMeshIOService ioService_write_Copy_grid(copyGridMesh);
		gmds::VTKWriter vtkWriterGba(&ioService_write_Copy_grid);
		vtkWriterGba.setCellOptions(gmds::N|gmds::F);
		vtkWriterGba.setDataOptions(gmds::N|gmds::F);
		vtkWriterGba.write("Q_learning_test.vtk");

/*
		environment.deleteEnv();
		copyGrid.deleteGridAround();
		actionQLearning.deleteActions();
		*/


	}

	std::cout<<"Training complete !"<<std::endl;
	std::cout<<"LA Q TABLE "<<std::endl;
	for (auto l : politique->Q_Table ){
		std::cout<<"une ligne de la qtable : "<<l.front()<<std::endl;
		for (auto d : l ){
			std::cout<<d<<std::endl;
		}
	}
	exit(0);
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