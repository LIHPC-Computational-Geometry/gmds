//
// Created by calderans on 28/01/20.
//

#include "medusa/control/MedusaCommandPick.h"

using namespace medusa;

MedusaCommandPick::MedusaCommandPick(MediatorPolicyItf* APolicy, GraphicView* AGraphicView, TextView* ATextView,int ADIM, vtkIdType ACellID):MedusaCommand(APolicy,AGraphicView,ATextView) {

    m_DIM = ADIM;
    m_CellID = ACellID;
}
/*----------------------------------------------------------------------------*/

MedusaCommand::State MedusaCommandPick::execute() {
    //m_policy->showAxis(m_graphic_collaborator,m_CellID);

    int axis = 0;
    //axis = m_policy->selectVTKCell(m_text_collaborator,m_DIM,m_CellID);


    //if(axis != 0) {
        //We pick on the grid[0] which is the tet mesh
        m_surfID = m_policy->selectVTKCell(m_graphic_collaborator, m_DIM, m_CellID, axis);

        if(m_surfID != -1) {
            std::cout<<"SUCCES: [pick]"<<std::endl;
            return SUCCES;
        }
    //}
    return FAILURE;
}
/*----------------------------------------------------------------------------*/
MedusaCommand::State MedusaCommandPick::undo() {

    m_policy->removeItem(m_surfID, m_graphic_collaborator);
    m_policy->removeItem(m_surfID, m_text_collaborator);
}
/*----------------------------------------------------------------------------*/
