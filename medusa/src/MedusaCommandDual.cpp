//
// Created by calderans on 28/01/20.
//

#include "medusa/control/MedusaCommandDual.h"

using namespace medusa;

MedusaCommandDual::MedusaCommandDual(MediatorPolicyItf* APolicy, GraphicView* AGraphicView, TextView* ATextView):MedusaCommand(APolicy,AGraphicView,ATextView) {

}
/*----------------------------------------------------------------------------*/

MedusaCommand::State MedusaCommandDual::execute() {

    m_policy->generateDual(m_text_collaborator);
    m_policy->generateDual(m_graphic_collaborator);
    return SUCCES;
}
/*----------------------------------------------------------------------------*/
MedusaCommand::State MedusaCommandDual::undo() {
    std::cout<<"Undo dual"<<std::endl;
    m_policy->resetDual(m_graphic_collaborator);
}
/*----------------------------------------------------------------------------*/
