//
// Created by bourmaudp on 10/05/22.
//
#include "gmds/paul/Environment.h"
using namespace gmds;

Environment::Environment(GridBuilderAround *AGrid, Mesh *AMesh)
	:g_grid(*AGrid),m_mesh(*AMesh){;}

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
	std::cout<<"Global IoU : "<<globalIoU<<std::endl;
	return globalIoU;
}
double Environment::localIoU(const gmds::Face AFaceID)
{
	double localIoU = g_grid.m_mesh.getVariable<double,GMDS_FACE>("volFrac")->value(AFaceID.id());
	std::cout<<"Local IoU : "<<localIoU<<std::endl;
	return localIoU;
}
