/*----------------------------------------------------------------------------*/
#include "medusa/control/MediatorGeomPolicy.h"
/*----------------------------------------------------------------------------*/
#include <medusa/view/GraphicView.h>
#include <medusa/view/TextView.h>
#include <medusa/model/MedusaBackEnd.h>
/*----------------------------------------------------------------------------*/
using namespace medusa;
/*----------------------------------------------------------------------------*/
MediatorGeomPolicy::MediatorGeomPolicy() {}
/*----------------------------------------------------------------------------*/
MediatorGeomPolicy::~MediatorGeomPolicy() {}
/*----------------------------------------------------------------------------*/
int MediatorGeomPolicy::
selectVTKCell(TextView *AView, const int ADim, const vtkIdType AID)
{AView->selectVTKCell(ADim,AID);}
/*----------------------------------------------------------------------------*/
int MediatorGeomPolicy::
selectVTKCell(GraphicView *AView, const int ADim, const vtkIdType AID, int AAxis)
{

    std::vector<vtkIdType> ids;
    int surf_id = MedusaBackend::getInstance()->pickSurface(AID, ids);
    if(surf_id == -1)
        return -1;
    if(m_mode == 1){
        MedusaBackend::getInstance()->removeSurfFromSelection();
        MedusaBackend::getInstance()->unconstrainFrameTets(MedusaBackend::getInstance()->vtkToGMDSIDs(ids));
    }
    AView->surfacePicked(ids);

    return surf_id;
}
/*----------------------------------------------------------------------------*/
void MediatorGeomPolicy::
help(TextView *AView)
{AView->help();}
/*----------------------------------------------------------------------------*/
void MediatorGeomPolicy::remove(GraphicView *AView) {
    AView->unlockLines();
}
/*----------------------------------------------------------------------------*/
void MediatorGeomPolicy::removeItem(int AID, GraphicView *AView) {
    if(m_mode == 1){

        std::vector<gmds::TCellID> offset =
                MedusaBackend::getInstance()->getOffset(1, MedusaBackend::getInstance()->pickSingLine(AID));

        AView->showSetOfTets(MedusaBackend::getInstance()->gmdsToVtkIDs(offset));

        char more = '0';
        bool again = true;
        while(again){
            std::cout<<"Press '+' to add increase the zone"<<std::endl;
            cin >> more;
            switch(more){
                case '+':
                    offset =MedusaBackend::getInstance()->getOffset(1, offset);
                    AView->showSetOfTets(MedusaBackend::getInstance()->gmdsToVtkIDs(offset));
                    break;
                default:
                    cin.clear();
                    again = false;

            }
        }

        MedusaBackend::getInstance()->unconstrainFrameTets(offset);
        AView->addTetSetActor();

    }else if(m_mode == 2) {
        MedusaBackend::getInstance()->removeSurfFromSelection();
        AView->removeLastSelection();
    }
}
/*----------------------------------------------------------------------------*/
void MediatorGeomPolicy::removeItem(int AID, TextView *AView) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorGeomPolicy::generate(GraphicView *AView) {


    //m_mode = 1;

    if(m_mode == 0){
        AView->showPoints(MedusaBackend::getInstance()->wrongGeomPoint());
        //MedusaBackend::getInstance()->generateFrameField();
        int nb_sing_lines = MedusaBackend::getInstance()->getNbSingLines();
        for(int i = 0; i<nb_sing_lines; i++){
            AView->buildSingLine(i, MedusaBackend::getInstance()->gmdsToVtkIDs(MedusaBackend::getInstance()->pickSingLine(i)));
        }
        m_mode = 1;
        AView->textMode(1);
    }else if(m_mode == 1){
        AView->update();
        //Etape pour gérer le frame field après sa génération
        AView->hideSetOfTets();
        AView->removeSelection();

        MedusaBackend::getInstance()->generateFrameField();
        AView->resetSingLine();

        int nb_sing_lines = MedusaBackend::getInstance()->getNbSingLines();
        std::cout<<"nb lines "<<nb_sing_lines<<std::endl;
        for(int i = 0; i<nb_sing_lines; i++){
            AView->buildSingLine(i, MedusaBackend::getInstance()->gmdsToVtkIDs(MedusaBackend::getInstance()->pickSingLine(i)));
        }
        MedusaBackend::getInstance()->noCut();

        AView->showPoints(MedusaBackend::getInstance()->wrongGeomPoint());
        /*char more = '0';
            std::cout<<"Do you want to go to geometry cutting session ? (y/n)"<<std::endl;
            cin >> more;
            switch(more){
                case 'n':
                    m_mode = 1;
                    break;
                case 'y':
                    m_mode = 2;
                    break;
                default:
                    cin.clear();

            }
            cin.clear();
        //m_mode = 2;
        AView->textMode(m_mode);*/
    }
    /*else if(m_mode == 2){
        //Generation de la coupe
        if(MedusaBackend::getInstance()->generateCut()==0) {
            AView->updateCut();
            m_mode = 3;
        }else{
            MedusaBackend::getInstance()->noCut();
            std::cout<<"ERROR: No surface selected"<<std::endl;
        }
    }*/else{

    }
}
/*----------------------------------------------------------------------------*/
void MediatorGeomPolicy::generate(TextView *AView) {
    if(m_mode < 2){
        std::cout<<"Frame field generated"<<std::endl;
    }else if(m_mode == 2){
        //Generation de la coupe
        std::cout<<"Cut generated"<<std::endl;
    }else{

    }
}
/*----------------------------------------------------------------------------*/
void MediatorGeomPolicy::undoGenerate(GraphicView *AView) {
    MedusaBackend::getInstance()->undoCut();
    AView->toggleBRep();
}
/*----------------------------------------------------------------------------*/
void MediatorGeomPolicy::showAxis(GraphicView *AView, vtkIdType ACellID) {
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorGeomPolicy::surfaceMode(GraphicView *AView){
    std::cout<<"Not implemented yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorGeomPolicy::singularityGraphToggle(GraphicView *AView) {
    AView->toggleSingGraph();
}
/*----------------------------------------------------------------------------*/
void MediatorGeomPolicy::opacity(GraphicView *AView) {
    if(m_mode == 0){

    } else if(m_mode == 1){
        AView->opacity();
    } else if(m_mode == 2){
        AView->toggleCutOpacity();
    }

}
/*----------------------------------------------------------------------------*/
void MediatorGeomPolicy::color(GraphicView *AView) {
    AView->toggleBRep();
}