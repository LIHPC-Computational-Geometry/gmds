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
            generate();
            break;

        case 'd':
        case 'D':
            viewDual();
            break;

        case 's':
        case 'S':
            singularityGraph();
            break;

        case 'r':
        case 'R':
            remove();
            break;

        case 'v':
        case 'V':
            surfaceMode();
            break;

        case 'u':
        case 'U':
            undo();
            break;

        case 'c':
        case 'C':
            color();
            break;

        case 'm':
        case 'M':
            mode();
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
void MouseInteractorStyle::opacity(){
    MediatorControl::getInstance()->opacity();
}
/*----------------------------------------------------------------------------*/
void MouseInteractorStyle::help(){
    MediatorControl::getInstance()->help();
}
/*----------------------------------------------------------------------------*/
void MouseInteractorStyle::generate(){
    MediatorControl::getInstance()->generate();
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

    vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    if (picker->GetCellId() != -1){

        std::cout<<"Picker cell id "<<picker->GetCellId()<<std::endl;

        MediatorControl::getInstance()->removeItem(MediatorControl::getInstance()->removeActor(picker->GetActor()));
    }

}
/*----------------------------------------------------------------------------*/
void MouseInteractorStyle::surfaceMode() {
    MediatorControl::getInstance()->surfaceMode();
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
void MouseInteractorStyle::color() {
    MediatorControl::getInstance()->color();
}
/*----------------------------------------------------------------------------*/
void MouseInteractorStyle::mode() {
    MediatorControl::getInstance()->changePolicy();
}