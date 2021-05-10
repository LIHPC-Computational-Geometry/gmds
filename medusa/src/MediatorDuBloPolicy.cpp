/*----------------------------------------------------------------------------*/
#include <medusa/control/MediatorDuBloPolicy.h>
/*----------------------------------------------------------------------------*/
#include <medusa/view/GraphicView.h>
#include <medusa/view/TextView.h>
#include <gmds/utils/Exception.h>
#include <medusa/model/MedusaBackEnd.h>
/*----------------------------------------------------------------------------*/
using namespace medusa;
/*----------------------------------------------------------------------------*/
MediatorDubloPolicy::MediatorDubloPolicy() {}
/*----------------------------------------------------------------------------*/
MediatorDubloPolicy::~MediatorDubloPolicy() {}
/*----------------------------------------------------------------------------*/
int MediatorDubloPolicy::
selectVTKCell(TextView *AView, const int ADim, const vtkIdType AID)
{   AView->selectVTKCell(ADim,AID);
    return AView->getInputAxis();}
/*----------------------------------------------------------------------------*/
int MediatorDubloPolicy::
selectVTKCell(GraphicView *AView, const int ADim, const vtkIdType AID, int AAxis)
{
    AView->removeAxis();
    std::vector<vtkIdType> surf;
    int sheet_ID = 0;
    if(AAxis == 3){
        surf = MedusaBackend::getInstance()->createBoundary(ADim, AID, sheet_ID);
        AView->updateMesh();
    } else{
        surf = MedusaBackend::getInstance()->createSurface(ADim, AID, sheet_ID, AAxis);
    }
    if(sheet_ID != -1) {
        AView->createSurface(surf, sheet_ID);
        AView->createPlanSurface(sheet_ID, MedusaBackend::getInstance()->grids().back());
    }

    return sheet_ID;
}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::
help(TextView *AView)
{AView->help();}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::
generateDual(GraphicView *AView)
{AView->updateDual(MedusaBackend::getInstance()->generateDual());}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::
generateDual(TextView *AView)
{AView->generateDual();}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::
remove(int ASheetID)
{MedusaBackend::getInstance()->removeSheet(ASheetID);}
/*----------------------------------------------------------------------------*/
/*void MediatorDubloPolicy::
remove(GraphicView *AView, int ASheetID)
{
  MedusaBackend::getInstance()->removeSheet(ASheetID);
    AView->remove(ASheetID);
}*/
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::
showAxis(GraphicView *AView, vtkIdType ACellID)
{
    //create a table to store the coordinates of the frame vectors
    //coords[0] is the center of the cell and coords[1,2,3] the three directions of the frame chart in the mesh cell
    double coords[3][3];
    MedusaBackend::getInstance()->getFrameAxis(coords, ACellID);
    AView->showAxis(coords);
}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::generateBlocks(GraphicView *AView) {
    MedusaBackend::getInstance()->generateBlocks();
    AView->viewBlocks();
}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::smoothBlocks(medusa::GraphicView *AView) {
    MedusaBackend::getInstance()->smoothBlocks();
    AView->viewBlocks();
}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::surfaceMode(GraphicView *AView) {
    AView->viewSurface();
}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::resetDual(GraphicView *AView) {
    MedusaBackend::getInstance()->resetDual();
    AView->setModeDomain();
}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::resetBlock(GraphicView *AView) {
    MedusaBackend::getInstance()->resetBlock();
    AView->toggleBlocks();
}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::resetBlockSmoothing(medusa::GraphicView *AView) {
    MedusaBackend::getInstance()->resetBlockSmoothing();
    AView->toggleBlocks();
}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::singularityGraphToggle(GraphicView *AView) {
    AView->toggleSingGraph();
}