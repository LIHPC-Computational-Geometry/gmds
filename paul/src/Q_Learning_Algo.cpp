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


	for (int i=0; i<=5;i++){
		std::cout <<"Dans le for : "<<i<<std::endl;


		std::cout<<"Copy mesh"<<std::endl;
		Mesh *copyGridMesh= new Mesh(environmentInit.g_grid.m_mesh->getModel());
		std::cout<<"Copy grid"<<std::endl;

		GridBuilderAround *copyGrid=new GridBuilderAround(copyGridMesh,environmentInit.g_grid.meshTarget,2);
		//copyGrid.executeGrid2D(5);
		std::cout<<"Clone"<<std::endl;

		toolInit.cloneMesh(environmentInit.g_grid,*copyGrid);
		//Tools toolsCopy(&copyGrid);
		std::cout<<"Action"<<std::endl;

		Actions *actionQLearning=new Actions(copyGrid);
		Tools *toolsQlearning=new Tools(copyGrid);
		std::cout<<"Env"<<std::endl;
		Environment *environment= new Environment(copyGrid,environmentInit.g_grid.meshTarget,actionQLearning,toolsQlearning);
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
		while(environment->globalIoU()<0.95 && countIte<130){

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
				//std::cout<<"Action delete"<<std::endl;
				//std::cout<<"local IoU : "<<localIoU<<std::endl;
				environment->executeAction(faceSelected, 0);
				actionIndex=0;
			}
			else {
				if(actionIndex == 1){
					std::vector<Node> listNodesFaceSelected = toolsQlearning->getListNodesOfFace(faceSelected.id());
					int cutDirection = toolsQlearning->bestCutDirection(faceSelected);
					environment->executeAction(faceSelected, actionIndex);
					if(cutDirection==0){
						auto listFaceCutHorizontal0 = toolsQlearning->getFacesCommon(listNodesFaceSelected[0].id(),listNodesFaceSelected[1].id());
						auto listFaceCutHorizontal1 = toolsQlearning->getFacesCommon(listNodesFaceSelected[3].id(),listNodesFaceSelected[2].id());

						std::vector<Face> listNewFaces;
						for(auto f0:listFaceCutHorizontal0){
							for(auto f1:listFaceCutHorizontal1){
								if(toolsQlearning->commonNodesFaces(f0,f1)){
									listNewFaces.push_back(f0);
									listNewFaces.push_back(f1);
								}
							}
						}

						if (copyGrid->m_mesh->getVariable<double,GMDS_FACE>("volFrac")->value(listNewFaces[0].id())<
							 copyGrid->m_mesh->getVariable<double,GMDS_FACE>("volFrac")->value(listNewFaces[1].id())){
							faceSelected = listNewFaces[1];
						}
						else{
							faceSelected = listNewFaces[0];
						}





					}else{

						auto listFaceCutVertical0 = toolsQlearning->getFacesCommon(listNodesFaceSelected[0].id(),listNodesFaceSelected[3].id());
						auto listFaceCutVertical1 = toolsQlearning->getFacesCommon(listNodesFaceSelected[1].id(),listNodesFaceSelected[2].id());

						std::vector<Face> listNewFaces;
						for(auto f0:listFaceCutVertical0){
							for(auto f1:listFaceCutVertical1){
								if(toolsQlearning->commonNodesFaces(f0,f1)){
									listNewFaces.push_back(f0);
									listNewFaces.push_back(f1);
								}
							}
						}
						if (copyGrid->m_mesh->getVariable<double,GMDS_FACE>("volFrac")->value(listNewFaces[0].id())<
							 copyGrid->m_mesh->getVariable<double,GMDS_FACE>("volFrac")->value(listNewFaces[1].id())){
							faceSelected = listNewFaces[1];
						}
						else{
							faceSelected = listNewFaces[0];
						}

					}
				}
				else {
					environment->executeAction(faceSelected, actionIndex);
				}
			}


			std::cout<<"La face select : "<<faceSelected<<std::endl;
			std::cout<<"Noeuds la face  : "<<std::endl;
			for (auto n : toolsQlearning->getListNodesOfFace(faceSelected.id())){
				std::cout<<n<<std::endl;
			}
			std::cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
			std::cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
			std::cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
			//std::cout<<"Calcul reward"<<std::endl;
			if(environment->localIoU(faceSelected)<0){
				std::cout<<"Local IoU negatif : "<<environment->localIoU(faceSelected)<<std::endl;
				std::cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
				std::cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
				std::cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
				break ;

			}
			double reward = environment->reward(faceSelected);
			//std::cout<<"Calcul TD"<<std::endl;
			double temporal_difference = reward + (0.9 * maxQValue) - oldQValue;
			//std::cout<<"Calcul newQvalue"<<std::endl;
			double newQValue = oldQValue + (0.9 * temporal_difference);
			//std::cout << "#nbface " << faceSelected.nbNodes() << std::endl;
			//std::cout<<"Calcul updateQTable"<<std::endl;
			politique->updateQTable(intervalIoU,actionIndex,newQValue);

			std::string s = std::to_string(i);

			std::string name = "training_results/"+s+".vtk";

			toolInit.saveMesh(copyGridMesh,name);

			std::cout<<"Count Ite : "<<countIte<<std::endl;

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
	saveQTable(politique->Q_Table,"Q_Table_Save.txt");
	exit(0);
}

void gmds::executeQLearning(Environment environment)
{
	Politique *qLearningPolicy=new Politique(&environment);
	qLearningPolicy->initQTable();
	qLearningPolicy->Q_Table=readQTable(qLearningPolicy->getQTable(),"Q_Table_Save.txt");



	std::cout<<"La lecture de la QTable"<<std::endl;
	for (auto l : qLearningPolicy->Q_Table ){
		std::cout<<"une ligne de la qtable : "<<l.front()<<std::endl;
		for (auto d : l ){
			std::cout<<d<<std::endl;
		}
	}

	while(environment.globalIoU()<0.9){
		environment.calcVolFrac();

		Face faceSelected = environment.faceSelect();

		std::cout<<"Face select "<<faceSelected<<std::endl;

		double localIoU = environment.localIoU(faceSelected);

		int intervalIoU = qLearningPolicy->getInterval(localIoU);

		int actionIndex =  qLearningPolicy->getNextActionQLearning(intervalIoU);

		std::cout<<"interval IoU "<<intervalIoU<<std::endl;

		environment.executeAction(faceSelected,actionIndex);

		gmds::IGMeshIOService ioService_write_Copy_grid(environment.g_grid.m_mesh);
		gmds::VTKWriter vtkWriterGba(&ioService_write_Copy_grid);
		vtkWriterGba.setCellOptions(gmds::N|gmds::F);
		vtkWriterGba.setDataOptions(gmds::N|gmds::F);
		vtkWriterGba.write("Q_learning_result.vtk");

	}




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

void gmds::saveQTable(std::vector<std::vector<double>> theQTable,std::string NameFile)
{
	std::ofstream myfile;
	myfile.open(NameFile);
	for(auto i :theQTable){
		for(auto j : i){
			myfile<<j<<"\n";
		}
	}
	myfile.close();
}

std::vector<std::vector<double>> gmds::readQTable(std::vector<std::vector<double>> theQTable, std::string NameFile)
{
	std::fstream myfile;
	myfile.open(NameFile,std::ios::in);
	int x;
	if (!myfile) {
		std::cout << "Unable to open file";
		exit(1); // terminate with error
	}

	std::string myString;
	if(myfile.is_open()){
		std::string tp;
		int nbLigne = 0;
		int ligneQtable=0;
		while(std::getline(myfile,tp)) {
			std::cout<<"Numero ligne "<<ligneQtable<<std::endl;

			std::stringstream ss;

			if(nbLigne != 0 && nbLigne%3 == 0){
				ligneQtable++;
			}

			std::cout << "My string " << tp << std::endl;
			ss << tp; // send it to the string stream

			double x;
			if(ss >> x) // send it to a double, test for correctness
			{
				int colonne = nbLigne %3;

				std::cout<<"La colonne : "<<colonne<<std::endl;
				std::cout << "success, " << " x = " << x << std::endl;
				theQTable[ligneQtable][colonne]=x;
			}
			else
			{
				std::cout << "error converting " << tp << std::endl;
			}


			nbLigne++;
		}
	}

	for (auto l : theQTable){
		std::cout<<"une ligne de la qtable : "<<l.front()<<std::endl;
		for (auto d : l ){
			std::cout<<d<<std::endl;
		}
	}

	return theQTable;

}