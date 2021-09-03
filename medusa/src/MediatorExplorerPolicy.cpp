/*----------------------------------------------------------------------------*/
#include <medusa/control/MediatorExplorerPolicy.h>
/*----------------------------------------------------------------------------*/
#include <medusa/view/GraphicView.h>
#include <medusa/view/TextView.h>
/*----------------------------------------------------------------------------*/
using namespace medusa;
/*----------------------------------------------------------------------------*/
MediatorExplorerPolicy::MediatorExplorerPolicy() {}
/*----------------------------------------------------------------------------*/
MediatorExplorerPolicy::~MediatorExplorerPolicy() {}
/*----------------------------------------------------------------------------*/
int MediatorExplorerPolicy::
selectVTKCell(TextView *AView, const int ADim, const vtkIdType AID)
{
    AView->selectVTKCell(ADim,AID);}
/*----------------------------------------------------------------------------*/
int MediatorExplorerPolicy::
selectVTKCell(GraphicView *AView, const int ADim, const vtkIdType AID, int AAxis)
{AView->selectVTKCell(ADim,AID);
    return 0;}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::
help(TextView *AView)
{AView->help();}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::remove(int ASheetID) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::generateDual(GraphicView *AView) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::generateDual(TextView *AView) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::showAxis(GraphicView *AView, vtkIdType ACellID) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::generateBlocks(medusa::GraphicView *AView) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::smoothBlocks(medusa::GraphicView *AView) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::surfaceMode(GraphicView *AView){
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::resetDual(GraphicView *AView) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::resetBlock(GraphicView *AView) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::resetBlockSmoothing(medusa::GraphicView *AView) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::singularityGraphToggle(GraphicView *AView) {
    std::cout<<"Not implemented yet"<<std::endl;
}
