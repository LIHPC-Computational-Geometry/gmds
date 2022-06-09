//
// Created by bourmaudp on 10/05/22.
//
#include "gmds/paul/Environment.h"
using namespace gmds;

Environment::Environment(GridBuilderAround *AGrid, Mesh *AMeshTarget,Actions *AAction)
	:g_grid(*AGrid),m_mesh(*AMeshTarget),action(*AAction){;}


//Environment::~Environment(){}

void Environment::deleteEnv()
{
	delete this;
}

double Environment::globalIoU()
{
	double globalIoU=0;
	int nbActivateFace=0;
	for (auto f : g_grid.m_mesh.faces()){
		if(g_grid.getActivate(g_grid.m_mesh.get<Face>(f))==1){
			nbActivateFace+=1;
			globalIoU+=Environment::localIoU(g_grid.m_mesh.get<Face>(f));
			//std::cout<<"Calcul global IoU, face : "<<f<<" avec local IoU de "<<Environment::localIoU(g_grid.m_mesh.get<Face>(f))<<std::endl;
		}
	}
	globalIoU=globalIoU/nbActivateFace;
	//std::cout<<"Global IoU : "<<globalIoU<<std::endl;
	return globalIoU;
}
double Environment::localIoU(const gmds::Face AFace)
{
	double localIoU = g_grid.m_mesh.getVariable<double,GMDS_FACE>("volFrac")->value(AFace.id());
	//std::cout<<"Local IoU : "<<localIoU<<std::endl;
	if (localIoU<0){
		std::cout<<"local IoU inf 0"<<std::endl;
		volfraccomputation_2d(&g_grid.m_mesh,&g_grid.meshTarget,g_grid.m_mesh.getVariable<double,GMDS_FACE>("volFrac"));
		localIoU = g_grid.m_mesh.getVariable<double,GMDS_FACE>("volFrac")->value(AFace.id());
		std::cout<<"Apres maj local IoU"<<localIoU<<std::endl;
	}
	return localIoU;
}

double Environment::reward(const Face AFace)
{
	double alpha = 0.1;
	double reward = globalIoU() + alpha * localIoU(AFace);

	return reward;
}

Face Environment::faceSelect()
{
	Face AFace;
	bool firstIte = true;
	double minLocalIoU;

	for (auto f : g_grid.m_mesh.faces()){
		Face actualFace = g_grid.m_mesh.get<Face>(f);
		double actualLocalIoU = Environment::localIoU(actualFace);
		if(g_grid.getActivate(g_grid.m_mesh.get<Face>(f))==1) {
			if (firstIte) {
				minLocalIoU = actualLocalIoU;
				AFace = actualFace;
				firstIte = false;
			}
			else if (actualLocalIoU < minLocalIoU) {
				minLocalIoU = actualLocalIoU;
				AFace = actualFace;
			}
		}
	}
	//std::cout << AFace << std::endl;
	return AFace;
}

void Environment::executeAction(Face AFace,int numberAction)
{
	if(numberAction==0){
	       action.executeDeleteFace(AFace.id());
	}
	else if(numberAction==1){
	    action.executeCutFace(AFace,1);
	}
	else if(numberAction==2){
		action.executeCutFace(AFace,0);
	}
	else if(numberAction==3){
	    action.executeGlideMinNodeFace(AFace);
	}
	else if(numberAction==4){
	    action.executeGlideMaxNodeFace(AFace);
	}

	volfraccomputation_2d(&g_grid.m_mesh,&g_grid.meshTarget,g_grid.m_mesh.getVariable<double,GMDS_FACE>("volFrac"));

}

void Environment::calcVolFrac()
{
	volfraccomputation_2d(&g_grid.m_mesh,&g_grid.meshTarget,g_grid.m_mesh.getVariable<double,GMDS_FACE>("volFrac"));
}
