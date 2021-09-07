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
void MediatorExplorerPolicy::remove(GraphicView *AView) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::removeItem(int AID, GraphicView *AView) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::removeItem(int AID, TextView *AView) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::generate(GraphicView *AView) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::generate(TextView *AView) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::undoGenerate(GraphicView *AView) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::showAxis(GraphicView *AView, vtkIdType ACellID) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::surfaceMode(medusa::GraphicView *AView) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::singularityGraphToggle(GraphicView *AView) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::opacity(GraphicView *AView) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorExplorerPolicy::color(GraphicView *AView) {
    std::cout<<"Not implemented yet"<<std::endl;
}
