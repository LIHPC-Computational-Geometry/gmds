/*----------------------------------------------------------------------------*/
#include <medusa/control/MediatorControl.h>

#include <medusa/model/MedusaBackEnd.h>
#include <medusa/control/MediatorExplorerPolicy.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
/*----------------------------------------------------------------------------*/
#include <gmds/utils/Exception.h>
#include <medusa/control/MedusaCommandGenerate.h>
#include <medusa/control/MediatorDuBloPolicy.h>
/*----------------------------------------------------------------------------*/
using namespace medusa;
using namespace gmds;
/*----------------------------------------------------------------------------*/
MediatorControl* MediatorControl::m_instance=NULL;
/*----------------------------------------------------------------------------*/
MediatorControl* MediatorControl::getInstance() {return m_instance;}
/*----------------------------------------------------------------------------*/
void MediatorControl::initialize() {
    if(m_instance==NULL)
        m_instance=new MediatorControl();
}/*----------------------------------------------------------------------------*/
MediatorControl::MediatorControl():m_policy(new MediatorExplorerPolicy()) {;}
/*----------------------------------------------------------------------------*/
MediatorControl::~MediatorControl() {;}
/*----------------------------------------------------------------------------*/
void MediatorControl::setPolicy(medusa::MediatorPolicyItf *APolicy) {
    if(APolicy==NULL)
        throw GMDSException("Error, null policy!!!");
    m_policy=APolicy;
}
/*----------------------------------------------------------------------------*/
void MediatorControl::update() {
    std::cout<<"GlobalEventListener: receive a notification from the back end"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MediatorControl::addCollaborator(GraphicView *AView) {
    m_graphic_collaborator.push_back(AView);
}
/*----------------------------------------------------------------------------*/
void MediatorControl::addCollaborator(TextView *AView) {
    m_text_collaborator.push_back(AView);
}

/*----------------------------------------------------------------------------*/
void MediatorControl::select(int ADIm, vtkIdType ACellID) {

    //if(MedusaBackend::getInstance()->grids().size()!=1)
    //throw GMDSException("Warning, picking is only available for a single mesh");

    MedusaCommandPick* c = new MedusaCommandPick(m_policy, m_graphic_collaborator[0],m_text_collaborator[0],ADIm, ACellID);
    MedusaCommand::State sc = c->execute();
    if(sc == MedusaCommand::SUCCES) {
        m_commands.push_back(c);
    }else{
        std::cout << "FAILED command : [pick]" << std::endl;
    }

}
/*----------------------------------------------------------------------------*/
void MediatorControl::help() {
    for(auto v:m_text_collaborator){
        m_policy->help(v);
    }
}
/*----------------------------------------------------------------------------*/
void MediatorControl::generate() {
    MedusaCommandGenerate *c = new MedusaCommandGenerate(m_policy, m_graphic_collaborator[0], m_text_collaborator[0]);
    MedusaCommand::State sc = c->execute();
    if (sc == MedusaCommand::SUCCES){
        m_commands.push_back(c);
    }else{
        std::cout << "FAILED command : [generate]" << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
void MediatorControl::toggleDual() {
    if(m_dual_state == -1){
        m_text_collaborator[0]->dualError();
    }
    else if(m_dual_state == 0) {
        m_graphic_collaborator[0]->setDualVisibilityON();
        m_dual_state = 1;
    }else{
        m_graphic_collaborator[0]->setDualVisibilityOFF();
        m_dual_state = 0;
    }
}
/*----------------------------------------------------------------------------*/
void MediatorControl::wireframe(){
    for(auto v:m_graphic_collaborator){
        v->wireframe();
    }
}
/*----------------------------------------------------------------------------*/
void MediatorControl::opacity(){
    m_policy->opacity(m_graphic_collaborator[0]);
}
/*----------------------------------------------------------------------------*/
void MediatorControl::remove(){
    m_policy->remove(m_graphic_collaborator[0]);
}
/*----------------------------------------------------------------------------*/
int MediatorControl::removeActor(vtkSmartPointer<vtkActor> ADeletedActor){


    //=========================================================================
    //Ici il faut descendre l'appel aux policy
    //=========================================================================

    return m_graphic_collaborator[0]->selectLine(ADeletedActor);
    return m_graphic_collaborator[0]->removeActor(ADeletedActor);
}
/*----------------------------------------------------------------------------*/
void MediatorControl::removeItem(int AItemID){
    m_text_collaborator[0]->remove(AItemID);
    m_policy->removeItem(AItemID, m_graphic_collaborator[0]);
}
/*----------------------------------------------------------------------------*/
void MediatorControl::surfaceMode() {
    m_policy->surfaceMode(m_graphic_collaborator[0]);
}
/*----------------------------------------------------------------------------*/
void MediatorControl::undo() {
    if(!m_commands.empty()) {
        MedusaCommand *c = m_commands.back();
        c->undo();
        m_commands.pop_back();
    }
}
/*----------------------------------------------------------------------------*/
void MediatorControl::singularityGraphToggle() {
    m_policy->singularityGraphToggle(m_graphic_collaborator[0]);
}
/*----------------------------------------------------------------------------*/
void MediatorControl::color() {
    m_policy->color(m_graphic_collaborator[0]);
}
/*----------------------------------------------------------------------------*/
void MediatorControl::changePolicy() {

    m_policy = new medusa::MediatorDubloPolicy;
    MedusaBackend::getInstance()->modeBlocks();
    m_graphic_collaborator[0]->textMode(3);
    //m_graphic_collaborator[0]->updateCut();
}
/*----------------------------------------------------------------------------*/
