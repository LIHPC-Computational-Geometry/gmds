//
// Created by calderans on 28/01/20.
//

#include "medusa/control/MedusaCommandGenerate.h"

using namespace medusa;

MedusaCommandGenerate::MedusaCommandGenerate(MediatorPolicyItf* APolicy, GraphicView* AGraphicView, TextView* ATextView):MedusaCommand(APolicy,AGraphicView,ATextView) {

}
/*----------------------------------------------------------------------------*/

MedusaCommand::State MedusaCommandGenerate::execute() {

    m_policy->generate(m_text_collaborator);
    m_policy->generate(m_graphic_collaborator);
    return SUCCES;
}
/*----------------------------------------------------------------------------*/
MedusaCommand::State MedusaCommandGenerate::undo() {
    std::cout<<"Undo generate"<<std::endl;
    m_policy->undoGenerate(m_graphic_collaborator);
}
/*----------------------------------------------------------------------------*/
