//
// Created by calderans on 28/01/20.
//

#include "medusa/control/MedusaCommandBlock.h"

using namespace medusa;

MedusaCommandBlock::MedusaCommandBlock(MediatorPolicyItf* APolicy, GraphicView* AGraphicView, TextView* ATextView):MedusaCommand(APolicy,AGraphicView,ATextView) {

}
/*----------------------------------------------------------------------------*/

MedusaCommand::State MedusaCommandBlock::execute() {
    m_policy->generateBlocks(m_graphic_collaborator);
    return SUCCES;
}
/*----------------------------------------------------------------------------*/
MedusaCommand::State MedusaCommandBlock::undo() {
    m_policy->resetBlock(m_graphic_collaborator);
}
/*----------------------------------------------------------------------------*/
