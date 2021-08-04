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
    double coords[3][3];
    int sheet_ID = 0;
    if(MedusaBackend::getInstance()->getFrameAxis(coords,AID)) {

        AView->showAxis(coords);
        char axis = '0';
        int axis_index = 0;
        bool wrong_answer = true;
        while (wrong_answer) {
            std::cout
                    << "Select dual surface normal direction : '1' for red, '2' for green or '3' to create a boundary surface"
                    << std::endl;
            cin >> axis;
            switch (axis) {
                case '1':
                    wrong_answer = false;
                    axis_index = 1;
                    break;
                case '2':
                    wrong_answer = false;
                    axis_index = 2;
                    break;
                case '3':
                    wrong_answer = false;
                    axis_index = 3;
                    break;
                case 'q':
                    wrong_answer = false;
                    axis_index = 0;
                    break;
                default:
                    std::cout << "Wrong input, please retry." << std::endl;
                    cin.clear();
            }
        }
        std::cout << "Test" << std::endl;
        AView->removeAxis();
        std::vector<vtkIdType> surf;

        if (axis_index == 3) {
            surf = MedusaBackend::getInstance()->createBoundary(ADim, AID, sheet_ID);
            AView->updateMesh();
        } else {
            surf = MedusaBackend::getInstance()->createSurface(ADim, AID, sheet_ID, axis_index);
        }
        if (sheet_ID != -1) {
            AView->createSurface(surf, sheet_ID);
            AView->createPlanSurface(sheet_ID, MedusaBackend::getInstance()->grids().back());
        }
    }
    else{
        sheet_ID = -1;
    }

    return sheet_ID;
}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::
help(TextView *AView)
{AView->help();}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::
generate(GraphicView *AView)
{
    if(m_mode == 0) {
        AView->updateDual(MedusaBackend::getInstance()->generateDual());
        m_mode = 1;
    } else if(m_mode == 1){
        MedusaBackend::getInstance()->generateBlocks();
        AView->viewBlocks();
        m_mode = 2;
    }
}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::
generate(TextView *AView)
{
    if(m_mode == 0) {
        AView->generateDual();
        //m_mode = 1;
    } else if(m_mode == 1){
        //m_mode = 2;
    } else if(m_mode == 2){
        std::cout<<"Blocks already generated"<<std::endl;
    }
}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::remove(GraphicView *AView) {
    AView->unlockSheet();
}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::
removeItem(int AID, GraphicView *AView){
    MedusaBackend::getInstance()->removeSheet(AID);
    AView->removeActor(AID);
}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::
removeItem(int AID, TextView *AView)
{MedusaBackend::getInstance()->removeSheet(AID);}
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
void MediatorDubloPolicy::surfaceMode(GraphicView *AView) {
    AView->viewSurface();
}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::undoGenerate(GraphicView *AView) {
    if(m_mode == 1) {
        MedusaBackend::getInstance()->resetDual();
        AView->resetDual();
        AView->setModeDomain();
        m_mode = 0;
    } else if(m_mode == 2){
        MedusaBackend::getInstance()->resetBlock();
        AView->toggleBlocks();
        m_mode = 1;
    }
}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::singularityGraphToggle(GraphicView *AView) {
    AView->toggleSingGraph();
}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::opacity(GraphicView *AView) {
    AView->opacity();
}
/*----------------------------------------------------------------------------*/
void MediatorDubloPolicy::color(GraphicView *AView) {
    if(m_mode == 1) {
        int nbZones;
        std::cout<<"Press \"e\" to erase wrong dual regions or \"c\" to continue with this dual structure"<<std::endl;
        char input;
        cin>>input;
        if(input == 'e'){
            nbZones = MedusaBackend::getInstance()->correctDualRegions();
            AView->updateDual(nbZones);
        }else if (input == 'c') {
            nbZones = MedusaBackend::getInstance()->reGenerateDual();
        } else{

            throw std::string("Error! wrong input");
        }

        AView->updateDual(nbZones);
        cin.clear();
    }
}
