/*----------------------------------------------------------------------------*/
#include <medusa/control/MediatorControl.h>

#include <medusa/model/MedusaBackEnd.h>
#include <medusa/control/MediatorExplorerPolicy.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
/*----------------------------------------------------------------------------*/
#include <gmds/utils/Exception.h>
#include <medusa/control/MedusaCommandDual.h>
#include <medusa/control/MedusaCommandBlock.h>
#include <medusa/control/MedusaCommandBlockSmoothing.h>
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

    MedusaCommandSheet* c = new MedusaCommandSheet(m_policy, m_graphic_collaborator[0],m_text_collaborator[0],ADIm, ACellID);
    MedusaCommand::State sc = c->execute();
    if(sc == MedusaCommand::SUCCES) {
        std::cout << "SUCCEED to create dual sheet" << std::endl;
        m_commands.push_back(c);
    }else{
        std::cout << "FAILED to create dual sheet" << std::endl;
    }

}
/*----------------------------------------------------------------------------*/
void MediatorControl::help() {
    for(auto v:m_text_collaborator){
        m_policy->help(v);
    }
}
/*----------------------------------------------------------------------------*/
void MediatorControl::generateDual() {
    MedusaCommandDual* c = new MedusaCommandDual(m_policy, m_graphic_collaborator[0],m_text_collaborator[0]);
    MedusaCommand::State sc = c->execute();
    if(sc == MedusaCommand::SUCCES)
        m_commands.push_back(c);
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
void MediatorControl::solid(){
    for(auto v:m_graphic_collaborator){
        v->solid();
    }
}
/*----------------------------------------------------------------------------*/
void MediatorControl::opacity(){
    for(auto v:m_graphic_collaborator){
        v->opacity();
    }
}
/*----------------------------------------------------------------------------*/
void MediatorControl::remove(){
    m_graphic_collaborator[0]->remove();
}
/*----------------------------------------------------------------------------*/
int MediatorControl::removeActor(vtkSmartPointer<vtkActor> ADeletedActor){
    return m_graphic_collaborator[0]->removeActor(ADeletedActor);
}
/*----------------------------------------------------------------------------*/
void MediatorControl::removeSheet(int ASheetID){
    m_text_collaborator[0]->remove(ASheetID);
    m_policy->remove(ASheetID);
}
/*----------------------------------------------------------------------------*/
void MediatorControl::generateBlocks() {
    MedusaCommandBlock* c = new MedusaCommandBlock(m_policy, m_graphic_collaborator[0],m_text_collaborator[0]);
    MedusaCommand::State sc = c->execute();
    if(sc == MedusaCommand::SUCCES) {
        m_commands.push_back(c);
    }
}
/*----------------------------------------------------------------------------*/
void MediatorControl::finalMesh() {
    MedusaCommandBlockSmoothing* c = new MedusaCommandBlockSmoothing(m_policy, m_graphic_collaborator[0], m_text_collaborator[0]);
    MedusaCommand::State sc = c->execute();
    if(sc == MedusaCommand::SUCCES) {
        m_commands.push_back(c);
    }
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
