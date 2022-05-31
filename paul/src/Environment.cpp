//
// Created by bourmaudp on 10/05/22.
//
#include "gmds/paul/Environment.h"
using namespace gmds;

Environment::Environment(GridBuilderAround *AGrid, Mesh *AMesh,Actions *AAction)
	:g_grid(*AGrid),m_mesh(*AMesh),action(*AAction){;}

double Environment::globalIoU()
{
	double globalIoU=0;
	int nbActivateFace=0;
	for (auto f : g_grid.m_mesh.faces()){
		if(g_grid.getActivate(g_grid.m_mesh.get<Face>(f))==1){
			nbActivateFace+=1;
			globalIoU+=Environment::localIoU(g_grid.m_mesh.get<Face>(f));
		}
	}
	globalIoU=globalIoU/nbActivateFace;
	//std::cout<<"Global IoU : "<<globalIoU<<std::endl;
	return globalIoU;
}
double Environment::localIoU(const gmds::Face AFace)
{
	volfraccomputation_2d(&g_grid.m_mesh,&g_grid.meshTarget,g_grid.m_mesh.getVariable<double,GMDS_FACE>("volFrac"));
	double localIoU = g_grid.m_mesh.getVariable<double,GMDS_FACE>("volFrac")->value(AFace.id());
	//std::cout<<"Local IoU : "<<localIoU<<std::endl;
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
	    action.executeGlideMaxNodeFace(AFace);
	}
	else if(numberAction==4){
	    action.executeGlideMinNodeFace(AFace);
	}

}
