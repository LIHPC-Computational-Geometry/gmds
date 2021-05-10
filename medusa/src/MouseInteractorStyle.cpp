/*----------------------------------------------------------------------------*/
#include <medusa/control/MediatorControl.h>
#include "medusa/control/MouseInteractorStyle.h"
/*----------------------------------------------------------------------------*/
using namespace medusa;
/*----------------------------------------------------------------------------*/
vtkStandardNewMacro(MouseInteractorStyle);
/*----------------------------------------------------------------------------*/
MouseInteractorStyle::MouseInteractorStyle():
        vtkInteractorStyleTrackballCamera() {

}
/*----------------------------------------------------------------------------*/
void MouseInteractorStyle::OnChar() {

    vtkInteractorStyleTrackballCamera::OnChar();
    switch (this->Interactor->GetKeyCode())
    {
        case 'p':
        case 'P':
            pick();
            break;

        case 'w':
        case 'W':
            wireframe();
            break;

        case 's':
        case 'S':
            solid();
            break;

        case 'o':
        case 'O':
            opacity();
            break;

        case 'h':
        case 'H':
            help();
            break;

        case 'g':
        case 'G':
            generateDual();
            break;

        case 'd':
        case 'D':
            viewDual();
            break;

        case 'r':
        case 'R':
            remove();
            break;

        case 'b':
        case 'B':
            generateBlocks();
            break;

        case 'v':
        case 'V':
            surfaceMode();
            break;

        case 'u':
        case 'U':
            undo();
            break;

        case 'm':
        case 'M':
            finalMesh();
            break;
    }
}
//----------------------------------------------------------------------------
void MouseInteractorStyle::pick() {

    // Get the location of the click (in window coordinates)
    int *pos = this->GetInteractor()->GetEventPosition();

    vtkSmartPointer<vtkCellPicker> picker =
            vtkSmartPointer<vtkCellPicker>::New();
    picker->SetTolerance(0.0005);

    // Pick from this location.
    picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

    double *worldPosition = picker->GetPickPosition();

    if (picker->GetCellId() != -1){
        MediatorControl::getInstance()->select(3, picker->GetCellId());
    }
    // Forward events
    vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
}
/*----------------------------------------------------------------------------*/
void MouseInteractorStyle::setPolicy(){

    //MediatorControl::setPolicy();
}
/*----------------------------------------------------------------------------*/
void MouseInteractorStyle::wireframe(){
    MediatorControl::getInstance()->wireframe();
}
/*----------------------------------------------------------------------------*/
void MouseInteractorStyle::solid(){
    MediatorControl::getInstance()->solid();
}
/*----------------------------------------------------------------------------*/
void MouseInteractorStyle::opacity(){
    MediatorControl::getInstance()->opacity();
}
/*----------------------------------------------------------------------------*/
void MouseInteractorStyle::help(){
    MediatorControl::getInstance()->help();
}
/*----------------------------------------------------------------------------*/
void MouseInteractorStyle::generateDual(){
    MediatorControl::getInstance()->generateDual();
}
/*----------------------------------------------------------------------------*/
void MouseInteractorStyle::viewDual(){
    MediatorControl::getInstance()->toggleDual();
}
/*----------------------------------------------------------------------------*/
void MouseInteractorStyle::remove(){
    MediatorControl::getInstance()->remove();

    int *pos = this->GetInteractor()->GetEventPosition();

    vtkSmartPointer<vtkCellPicker> picker =
            vtkSmartPointer<vtkCellPicker>::New();
    picker->SetTolerance(0.0005);

    // Pick from this location.
    picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

    if (picker->GetCellId() != -1){

        std::cout<<"Picker cell id "<<picker->GetCellId()<<std::endl;

        MediatorControl::getInstance()->removeSheet(MediatorControl::getInstance()->removeActor(picker->GetActor()));
    }
    vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
}
/*----------------------------------------------------------------------------*/
void MouseInteractorStyle::generateBlocks() {
    MediatorControl::getInstance()->generateBlocks();
}
/*----------------------------------------------------------------------------*/
void MouseInteractorStyle::surfaceMode() {
    MediatorControl::getInstance()->surfaceMode();
}
/*----------------------------------------------------------------------------*/
void MouseInteractorStyle::finalMesh()  {
    MediatorControl::getInstance()->finalMesh();
}
/*----------------------------------------------------------------------------*/
void MouseInteractorStyle::undo() {
    MediatorControl::getInstance()->undo();
}
/*----------------------------------------------------------------------------*/
void MouseInteractorStyle::singularityGraph() {
    MediatorControl::getInstance()->singularityGraphToggle();
}
/*----------------------------------------------------------------------------*/